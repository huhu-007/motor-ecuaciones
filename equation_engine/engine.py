"""EquationEngine: fachada que orquesta parser -> dominio -> búsqueda -> Bolzano -> bisección ->
reconocimiento -> verificación exacta.  Devuelve estructuras JSON-serializables."""
from __future__ import annotations
import json
import time
from dataclasses import dataclass, field, asdict
from fractions import Fraction
from typing import Any, Dict, List, Optional

from .parser import parse_equation, ParseError, ParsedEquation
from .nodes import to_text, variables, contains_var
from .domain import DomainAnalyzer
from .evaluator import FunctionEvaluator
from .backends import BackendPool
from .errors import DomainError, NotExact
from .exact import ExactBackend, Rad, ConstMono
from .bisection import BisectionDigits, make_sign_fn, PrecisionExhausted, UndefinedInside
from .bolzano import BOLZANO_NO_APLICABLE
from .rootfinder import RootFinder, SearchConfig, RootInterval
from .polynomial import (reduced_ratfunc, poly_text, sturm_chain, count_roots, pgcd, pderiv, pdivmod,
                         deg as pdeg, psign, to_int_coeffs, trim)
from .recognizers import ExactFormRecognizer, Candidate, rad_to_str, PolynomialAnalyzer


class ExpressionParserFacade:
    @staticmethod
    def parse(text):
        from .parser import ExpressionParser
        return ExpressionParser().parse(text)


def constmono_to_str(c, pretty=True) -> str:
    """c·π^a·e^q  ->  '2π', 'π/2', '-(3π)/4', 'e^3'..."""
    from .recognizers import _qtxt
    q = c.c
    parts = []
    if c.pi:
        parts.append(("π" if pretty else "pi") + ("" if c.pi == 1 else f"^{c.pi}"))
    if c.e:
        parts.append("e" if c.e == 1 else f"e^{_qtxt(c.e) if c.e.denominator == 1 else '(' + _qtxt(c.e) + ')'}")
    if not parts:
        return _qtxt(q)
    body = ("·" if pretty else "*").join(parts)
    n, d = abs(q.numerator), q.denominator
    num = body if n == 1 else (f"{n}{body}" if pretty else f"{n}*{body}")
    if d > 1:
        num = f"({num})/{d}" if n > 1 else f"{num}/{d}"
    return ("-" if q < 0 else "") + num


@dataclass
class EngineConfig:
    limit: Fraction = Fraction(100)
    step: Fraction = Fraction(1, 10)
    digits: int = 30                 # decimales mostrados por raíz
    analysis_digits: int = 60        # decimales usados para reconocer formas exactas
    max_analysis_digits: int = 120   # precisión máxima si hace falta subir dinámicamente
    analysis_budget: float = 8.0     # segundos por raíz en el reconocimiento
    max_radical_index: int = 100
    max_poly_degree: int = 12
    max_analyzed_roots: int = 10
    analyze: bool = True


def fraction_to_decimal(fr: Fraction, digits: int) -> str:
    sign = "-" if fr < 0 else ""
    fr = abs(fr)
    ip = fr.numerator // fr.denominator
    rem = fr - ip
    scaled = (rem * 10 ** digits).numerator // (rem * 10 ** digits).denominator
    return f"{sign}{ip}" + (f".{scaled:0{digits}d}" if digits > 0 else "")


class EquationEngine:
    def __init__(self, config: Optional[EngineConfig] = None, **overrides):
        self.cfg = config or EngineConfig()
        for k, v in overrides.items():
            setattr(self.cfg, k, v)
        self.pool = BackendPool()
        self.search_cfg = SearchConfig(limit=Fraction(self.cfg.limit), step=Fraction(self.cfg.step))
        self.finder = RootFinder(self.pool, self.search_cfg)
        self.recognizer = ExactFormRecognizer(self.cfg.max_radical_index, self.cfg.max_poly_degree,
                                              self.cfg.analysis_budget)
        self.exact = ExactBackend()

    # ============================================================== API pública
    def solve(self, text: str) -> Dict[str, Any]:
        out = self._skeleton(text)
        try:
            eq = parse_equation(text)
        except ParseError as e:
            out.update(success=False, status="invalid_input", messages=[str(e)], error=str(e))
            return out
        ev = FunctionEvaluator(eq.function, eq.variable)
        cons = DomainAnalyzer(eq.variable).analyze(eq.function)
        out.update(variable=eq.variable, relation=eq.relation,
                   normalized_function=to_text(eq.function),
                   domain={"text": DomainAnalyzer.text(cons), "constraints": [c.to_dict() for c in cons]})
        if eq.relation != "=":
            out.update(success=False, status="unsupported_relation",
                       messages=["La relación se ha interpretado correctamente, pero en esta versión solo se "
                                 "resuelven ecuaciones (=). Las inecuaciones no se resuelven todavía."])
            return out
        if not contains_var(eq.function):
            return self._constant_case(out, ev)
        rf = reduced_ratfunc(eq.function, eq.variable)
        if rf is not None:
            return self._solve_rational(out, eq, ev, rf)
        return self._solve_generic(out, eq, ev)

    def analyze_expression(self, text: str, digits: int = 60) -> Dict[str, Any]:
        """Analiza una expresión CONSTANTE (p. ej. '√3+√2', '√2/√3', '⁴√16', '¹⁰⁰√(7/3)').
        1) Si se puede evaluar con aritmética exacta, la forma exacta es una DEMOSTRACIÓN.
        2) Si no, se evalúa con alta precisión y se reconoce numéricamente (solo evidencia)."""
        out = {"success": False, "expression": text, "exact_form": None, "exact_form_text": None,
               "exactness": "none", "decimal_form": None, "digits": digits, "warnings": [], "recognition": None}
        try:
            node = ExpressionParserFacade.parse(text)
        except ParseError as e:
            out.update(status="invalid_input", error=str(e))
            return out
        if variables(node):
            out.update(status="invalid_input", error="La expresión contiene variables; usa solve() para ecuaciones.")
            return out
        ev = FunctionEvaluator(node, "x")
        be = self.pool.mp(digits + 30)
        try:
            val = ev.evaluate(be, be.from_fraction(Fraction(0)))
        except DomainError as e:
            out.update(status="undefined", error=f"fuera del dominio real: {e}")
            return out
        out["decimal_form"] = be.ctx.nstr(val, digits + 1)
        try:
            z = ev.evaluate(self.exact, Rad.rat(0))
            pretty = rad_to_str(z) if isinstance(z, Rad) else constmono_to_str(z, True)
            text_ = rad_to_str(z, False) if isinstance(z, Rad) else constmono_to_str(z, False)
            out.update(success=True, status="exact", exact_form=pretty, exact_form_text=text_,
                       exactness="exact_proven",
                       proof={"method": "evaluación simbólica exacta (Q-álgebra de radicales / π, e)"})
            return out
        except (NotExact, ZeroDivisionError):
            pass
        except DomainError as e:
            out.update(status="undefined", error=str(e))
            return out
        rep = self.recognizer.recognize(self._mp_to_fraction(val), digits)
        out.update(success=True, status="numerical_only" if rep.best is None else "numerical_evidence",
                   recognition=rep.to_dict())
        if rep.best is not None:
            out["candidate_form"] = rep.best.pretty
            out["exactness"] = rep.best.level
        out["warnings"].append("No se pudo evaluar de forma exacta: lo anterior es evidencia numérica, no una demostración.")
        return out

    def analyze_number(self, value, digits: int) -> Dict[str, Any]:
        """Analiza un número ya calculado (cadena decimal, Fraction o mpf) con `digits` decimales ciertos.
        Si es una cadena, los dígitos están TRUNCADOS: se centra el valor sumando media unidad del último
        decimal, y `digits` nunca puede superar los decimales que realmente trae la cadena."""
        note = None
        if isinstance(value, str):
            txt = value.strip()
            neg = txt.startswith("-")
            ip, _, fp = txt.lstrip("-+").partition(".")
            if digits > len(fp):
                note = (f"se pidieron {digits} dígitos pero la cadena solo trae {len(fp)} decimales: "
                        f"se usan {len(fp)}")
                digits = len(fp)
            fr = Fraction(int(ip or 0))
            if fp:
                fr += Fraction(int(fp), 10 ** len(fp)) + Fraction(1, 2 * 10 ** len(fp))
            value = -fr if neg else fr
        rep = self.recognizer.recognize(value, digits)
        d = rep.to_dict()
        if note:
            d["notes"].insert(0, note)
        return d

    @staticmethod
    def _mp_to_fraction(v):
        sign, man, exp, _ = v._mpf_
        f = Fraction(int(man)) * (Fraction(2) ** int(exp))
        return -f if sign else f

    def to_json(self, result: dict, **kw) -> str:
        return json.dumps(result, ensure_ascii=False, default=str, **kw)

    # ---------------- preparación para el CLI interactivo (intervalos + flujo de dígitos)
    def prepare(self, text: str):
        """Devuelve (equation, evaluator, intervals, result_skeleton). Para uso interactivo."""
        eq = parse_equation(text)
        ev = FunctionEvaluator(eq.function, eq.variable)
        rf = reduced_ratfunc(eq.function, eq.variable)
        if rf is not None and pdeg(rf[0]) >= 1:
            res = self.finder.scan_polynomial(rf[0])["res"]
        else:
            res = self.finder.scan_generic(ev, make_sign_fn(ev, self.pool))
        return eq, ev, res

    # ============================================================== casos
    def _skeleton(self, text):
        return {"success": False, "status": None, "equation": text, "variable": None, "relation": None,
                "normalized_function": None, "domain": None,
                "search_range": [float(-self.cfg.limit), float(self.cfg.limit)],
                "method": None, "bolzano_applicable": False, "bolzano_message": None,
                "roots": [], "unverified_candidates": [], "rejected_sign_changes": [],
                "proven_no_real_solution": False, "warnings": [], "messages": [], "diagnostics": {}}

    def _constant_case(self, out, ev):
        try:
            z = ev.evaluate(self.exact, Rad.rat(0))
            s = z.sign()
            if s == 0:
                out.update(success=True, status="identity", method="exact_constant",
                           messages=["La igualdad es una identidad: se cumple para todo x del dominio."])
            else:
                out.update(success=True, status="no_real_solution", proven_no_real_solution=True,
                           method="exact_constant",
                           messages=["La ecuación no contiene la incógnita y es falsa (demostrado de forma exacta): "
                                     "no tiene solución real."])
            return out
        except (NotExact, ZeroDivisionError):
            pass
        except DomainError as e:
            out.update(success=False, status="undefined_in_range", method="exact_constant",
                       messages=[f"La expresión no está definida en ℝ: {e}"])
            return out
        be = self.pool.mp(40)
        try:
            v = ev.evaluate(be, be.from_fraction(Fraction(0)))
            nz = abs(v) > be.ctx.mpf(10) ** -30
            out.update(success=True, status="no_real_solution" if nz else "identity", method="numerical_constant",
                       warnings=["Resultado basado en evaluación numérica (40 dígitos), no en una demostración exacta."],
                       messages=["Constante no nula: no hay solución." if nz else "Constante ≈ 0: posible identidad."])
        except DomainError as e:
            out.update(success=False, status="undefined_in_range", messages=[f"No definida: {e}"])
        return out

    # ---------------------------------------------------------------- racional / polinómico
    def _solve_rational(self, out, eq, ev, rf):
        num, den, common = rf
        if not trim(num):
            out.update(success=True, status="identity", method="exact_polynomial",
                       messages=["La ecuación es una identidad en su dominio."])
            return out
        if pdeg(num) == 0:
            out.update(success=True, status="no_real_solution", proven_no_real_solution=True,
                       method="exact_polynomial",
                       messages=["El numerador es una constante no nula: la ecuación no tiene solución real (demostrado)."])
            return out
        t0 = time.monotonic()
        scan = self.finder.scan_polynomial(num)
        res = scan["res"]
        sq = pdivmod(num, pgcd(num, pderiv(num)))[0]
        chain_all = sturm_chain(sq)
        out["method"] = "bolzano_bisection_exact_polynomial"
        out["diagnostics"] = {"route": "rational_function",
                              "numerator": poly_text(num), "denominator": poly_text(den),
                              "square_free_factors": scan["factors"],
                              "real_roots_total": scan["total_real"], "real_roots_outside_range": scan["outside"]}
        if pdeg(common) > 0:
            out["warnings"].append(
                f"Factor común {poly_text(common)} entre numerador y denominador: sus ceros NO están en el "
                "dominio (no son soluciones).")
        if pdeg(den) > 0 or pdeg(common) > 0:
            out["warnings"].append("Bolzano se aplica al numerador polinómico (continuo en ℝ); los polos de f "
                                   "(ceros del denominador) no son raíces y no pueden producir falsas raíces.")
        ctx = {"evaluator": ev, "chain_all": chain_all, "poly": True, "_sqprod": sq}
        roots = self._build_roots(res, ev, ctx)
        out["roots"] = roots
        out["rejected_sign_changes"] = res.rejected
        if scan["outside"]:
            out["warnings"].append(f"{scan['outside']} raíz/raíces real(es) quedan FUERA del rango de búsqueda "
                                   f"[{float(-self.cfg.limit):g}, {float(self.cfg.limit):g}].")
        if roots:
            out.update(success=True, status="solved", bolzano_applicable=True)
            if any(r["multiplicity"] and r["multiplicity"] > 1 for r in roots):
                out["messages"].append("Hay raíces de multiplicidad > 1: f no cambia de signo en las de multiplicidad par; "
                                       "se localizan aplicando Bolzano al factor libre de cuadrados (Yun) y se "
                                       "certifican con Sturm.")
        elif scan["total_real"] == 0:
            out.update(success=True, status="no_real_solution", proven_no_real_solution=True,
                       messages=["La ecuación no tiene soluciones reales (demostrado: secuencia de Sturm)."])
        else:
            out.update(success=True, status="no_root_found_in_range",
                       messages=["No se ha encontrado ninguna raíz en el rango de búsqueda "
                                 "(existen raíces reales fuera de él)."])
        out["diagnostics"]["seconds"] = round(time.monotonic() - t0, 2)
        return out

    # ---------------------------------------------------------------- general
    def _solve_generic(self, out, eq, ev):
        t0 = time.monotonic()
        sign_fn = make_sign_fn(ev, self.pool)
        res = self.finder.scan_generic(ev, sign_fn)
        out["method"] = "bolzano_bisection_interval_arithmetic"
        out["diagnostics"] = {"route": "generic", "grid_points": res.total_points,
                              "undefined_grid_points": res.undefined_points}
        out["rejected_sign_changes"] = res.rejected
        ctx = {"evaluator": ev, "chain_all": None, "poly": False}
        roots = self._build_roots(res, ev, ctx)
        # raíces tangentes (sin cambio de signo): solo se aceptan si se VERIFICAN de forma exacta
        tang, unver = self._touching_roots(res.touch_candidates, roots, ev)
        roots += tang
        roots.sort(key=lambda r: r["_sortkey"])
        out["roots"], out["unverified_candidates"] = roots, unver
        out["warnings"] += res.notes
        if res.rejected:
            out["warnings"].append(f"{len(res.rejected)} cambio(s) de signo descartados por discontinuidad "
                                   "(polos / puntos no definidos): no son raíces.")
        if unver:
            out["warnings"].append("Hay posibles raíces de multiplicidad par sin cambio de signo que Bolzano no "
                                   "puede certificar (ver 'unverified_candidates').")
        if any(r["method"] != "tangent_root_exact" for r in roots):
            out["bolzano_applicable"] = True
        if roots:
            out.update(success=True, status="solved")
            if not out["bolzano_applicable"]:
                out["bolzano_message"] = BOLZANO_NO_APLICABLE
            out["warnings"].append(f"Solo se busca en [{float(-self.cfg.limit):g}, {float(self.cfg.limit):g}]; "
                                   "pueden existir más raíces fuera de ese rango.")
        else:
            self._no_roots_generic(out, res, ev)
        out["diagnostics"]["seconds"] = round(time.monotonic() - t0, 2)
        return out

    def _no_roots_generic(self, out, res, ev):
        rng = f"[{float(-self.cfg.limit):g}, {float(self.cfg.limit):g}]"
        if res.undefined_points == res.total_points:
            out.update(success=False, status="undefined_in_range",
                       messages=[f"La función no está definida en ningún punto de {rng} (dominio: "
                                 f"{out['domain']['text']})."])
            return
        proven = (res.undefined_points == 0 and not res.rejected and not res.touch_candidates
                  and self.finder.prove_no_root_in_range(ev))
        if proven:
            out.update(success=True, status="no_root_found_in_range",
                       messages=[f"No existe ninguna raíz en {rng}: demostrado con aritmética de intervalos. "
                                 "Esto NO implica que no haya soluciones fuera de ese rango."])
            out["diagnostics"]["proven_no_root_in_range"] = True
            return
        reasons = []
        if res.undefined_points:
            reasons.append("la función no está definida en parte del rango")
        if res.rejected:
            reasons.append("hay cambios de signo debidos a discontinuidades (no válidos para Bolzano)")
        if res.touch_candidates:
            reasons.append("hay posibles raíces sin cambio de signo (multiplicidad par)")
        out.update(success=True, status="bolzano_not_applicable", bolzano_applicable=False,
                   bolzano_message=BOLZANO_NO_APLICABLE,
                   messages=[BOLZANO_NO_APLICABLE,
                             "Esto NO significa que la ecuación no tenga solución: solo que no se halló un "
                             f"intervalo de Bolzano válido en {rng}."
                             + (" Motivos: " + "; ".join(reasons) + "." if reasons else "")])

    # ============================================================== construcción de raíces
    def _build_roots(self, res, ev, ctx) -> List[dict]:
        roots = []
        order = sorted(range(len(res.intervals)),
                       key=lambda i: abs(float(res.intervals[i].lo + res.intervals[i].hi) / 2))
        allowed = set(order[: self.cfg.max_analyzed_roots]) if self.cfg.analyze else set()
        for i, iv in enumerate(res.intervals):
            roots.append(self._one_root(iv, ev, ctx, analyze=i in allowed))
        roots.sort(key=lambda r: r["_sortkey"])
        for k, r in enumerate(roots, 1):
            r["index"] = k
        return roots

    def _one_root(self, iv: RootInterval, ev, ctx, analyze: bool) -> dict:
        D = self.cfg.digits
        r = {"index": None, "method": None, "multiplicity": iv.multiplicity,
             "bolzano": iv.report.to_dict() if iv.report else None,
             "interval": [str(iv.lo), str(iv.hi)], "origin": iv.origin,
             "decimal_form": None, "digits": 0, "exact_form": None, "exact_form_text": None,
             "exactness": "none", "candidate_form": None, "proof": None, "periodicity": None,
             "recognition": None, "minimal_polynomial": None, "analysis_skipped": False}
        # ---- raíz exacta en un punto (f(x)=0 comprobado con aritmética exacta)
        if iv.exact_point:
            x = iv.lo
            r.update(method="exact_zero_point", decimal_form=fraction_to_decimal(x, D), digits=D,
                     exact_form=self._q(x), exact_form_text=self._q(x), exactness="exact_proven",
                     proof={"method": "f(x)=0 en aritmética exacta", "unique_in_interval": True},
                     _sortkey=float(x))
            return r
        # ---- bisección (Bolzano) con dígitos confirmados
        bd = BisectionDigits(iv.sign_fn, iv.lo, iv.hi)
        try:
            bd.run(D)
            want = self.cfg.analysis_digits if (analyze and not bd.exact and not bd.periodo) else 0
            if want > len(bd.dec):
                bd.run(want, stop_on_period=True)
        except (PrecisionExhausted, UndefinedInside) as e:
            r.update(method="bolzano_bisection_failed", proof={"error": str(e)}, _sortkey=float(iv.lo))
            return r
        lo, hi = bd.bracket()
        r.update(method="bolzano_bisection", digits=len(bd.dec), decimal_form=bd.value_text(),
                 interval=[str(lo), str(hi)], _sortkey=float((lo + hi) / 2))
        if bd.exact or bd.periodo:
            fr = (lo if bd.exact else bd.rational_from_period())
            cand = Candidate("rational", self._q(fr), self._q(fr), Rad.rat(fr))
            cand.evidence = {"source": "colapso del intervalo a un punto exacto" if bd.exact
                             else "periodicidad de los decimales (≥20 dígitos repetidos)"}
            if bd.periodo:
                from .bisection import texto_periodo
                r["periodicity"] = {"period_notation": texto_periodo(bd.signo, bd.entero, bd.periodo),
                                    "rational": self._q(fr),
                                    "status": "numerical_evidence (verificada abajo con aritmética exacta)"}
            cand = self._make_verifier(ev, (lo, hi), ctx)(cand)
            if cand.level != "refuted":
                r.update(exact_form=cand.pretty if cand.level != "numerical_evidence" else None,
                         exact_form_text=cand.text if cand.level != "numerical_evidence" else None,
                         candidate_form=cand.pretty, exactness=cand.level, proof=cand.proof,
                         decimal_form=fraction_to_decimal(fr, D), digits=D)
            return r
        if not analyze:
            r["analysis_skipped"] = True
            return r
        self._recognize_into(r, bd, ev, ctx)
        return r

    def _recognize_into(self, r, bd, ev, ctx):
        D = len(bd.dec)
        for attempt in range(2):
            lo, hi = bd.bracket()
            x = (lo + hi) / 2
            verifier = self._make_verifier(ev, (lo, hi), ctx)
            polyctx = self._make_poly_context((lo, hi), ctx, ev) if ctx["poly"] else None
            rep = self.recognizer.recognize(x, D, verifier=verifier, poly_context=polyctx)
            poly_known = bool(rep.polynomial and rep.polynomial.get("found"))
            low_precision = any("insuficientes" in n or "no se puede comprobar" in n for n in rep.notes) \
                or (rep.polynomial is not None and not rep.polynomial.get("found")
                    and rep.polynomial.get("max_degree_tried", 0) < self.cfg.max_poly_degree)
            if rep.best is not None or poly_known or not low_precision or D >= self.cfg.max_analysis_digits:
                break
            new = min(self.cfg.max_analysis_digits, 2 * D)
            try:
                bd.run(new, stop_on_period=False)
            except (PrecisionExhausted, UndefinedInside):
                break
            D = len(bd.dec)
        r["recognition"] = rep.to_dict()
        r["digits"] = min(len(bd.dec), self.cfg.digits)
        r["digits_computed"] = len(bd.dec)
        r["decimal_form"] = self._trunc(bd, self.cfg.digits)
        if rep.polynomial and rep.polynomial.get("found"):
            pl = rep.polynomial
            r["minimal_polynomial"] = {"polynomial": pl["polynomial"], "degree": pl["degree"],
                                       "status": (pl.get("proof") or {}).get("status", "numerical_evidence"),
                                       "irreducible_over_Q": pl.get("irreducible_over_Q"),
                                       "radical_remark": pl.get("radical_remark")}
        if rep.best is not None:
            b = rep.best
            r["candidate_form"] = b.pretty
            r["exactness"] = b.level
            r["proof"] = b.proof
            if b.level in ("exact_proven", "exact_verified"):
                r["exact_form"], r["exact_form_text"] = b.pretty, b.text

    # ---------------- raíces tangentes
    def _touching_roots(self, cands, found, ev):
        tang, unver = [], []
        for c in cands[: 40]:
            x = c["x"]
            xf = float(x)
            if any(abs(float(r["_sortkey"]) - xf) < 1e-6 for r in found):
                continue
            rec = None
            if self.cfg.analyze and c["digits"] >= 8:
                from mpmath.ctx_mp import MPContext
                cx = MPContext(); cx.dps = 50
                xv = x if hasattr(x, "_mpf_") else cx.mpf(x)
                xfr = Fraction(int(cx.floor(xv * 10 ** 30)), 10 ** 30)
                rec_eng = ExactFormRecognizer(self.cfg.max_radical_index, 6, min(3.0, self.cfg.analysis_budget))
                verifier = self._make_verifier(ev, None, {"chain_all": None, "poly": False, "evaluator": ev})
                rep = rec_eng.recognize(xfr, c["digits"], verifier=verifier)
                if rep.best is not None and rep.best.level in ("exact_verified", "exact_proven"):
                    b = rep.best
                    val = b.numeric(cx) if b.numeric else xv
                    tang.append({"index": None, "method": "tangent_root_exact", "multiplicity": None,
                                 "bolzano": {"applicable": False,
                                             "reason": "f no cambia de signo en esta raíz: Bolzano no la detecta; "
                                                       "se ha verificado con aritmética exacta f(x)=0"},
                                 "interval": None, "origin": "touching_minimum",
                                 "decimal_form": cx.nstr(val, self.cfg.digits + 1), "digits": self.cfg.digits,
                                 "exact_form": b.pretty, "exact_form_text": b.text, "exactness": "exact_verified",
                                 "candidate_form": b.pretty, "proof": b.proof, "periodicity": None,
                                 "recognition": {"digits_available": c["digits"], "statement": rep.statement()},
                                 "analysis_skipped": False, "_sortkey": float(val)})
                    continue
            unver.append({"x_approx": float(x), "abs_f_at_minimum": c["abs_f"],
                          "reason": c["reason"] + "; no se pudo verificar de forma exacta (puede ser una raíz o solo un mínimo muy pequeño)."})
        return tang, unver

    # ============================================================== verificación exacta
    @staticmethod
    def _trunc(bd, n):
        txt = f"{bd.signo}{bd.entero}"
        return txt + ("." + "".join(map(str, bd.dec[:n])) if bd.dec and n > 0 else "")

    @staticmethod
    def _q(q: Fraction) -> str:
        return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"

    def _make_verifier(self, ev, bracket, ctx):
        chain_all = ctx.get("chain_all")
        pool = self.pool

        def unique_in_bracket():
            if bracket is None:
                return None
            if ctx.get("poly") and chain_all is not None:
                lo, hi = bracket
                return count_roots(chain_all, lo, hi) == 1
            return None

        def verify(cand: Candidate) -> Candidate:
            value = cand.value
            if value is None:
                cand.proof = {"f_candidate_is_exactly_zero": None,
                              "reason": "el candidato no es representable con el backend exacto: solo evidencia numérica"}
                return cand
            try:
                z = ev.evaluate(self.exact, value)
                zs = z.sign()
            except (NotExact, ZeroDivisionError) as e:
                cand.proof = {"f_candidate_is_exactly_zero": None,
                              "reason": f"f no se puede evaluar de forma exacta en el candidato ({e})"}
                return cand
            except DomainError as e:
                cand.level, cand.proof = "refuted", {"f_candidate_is_exactly_zero": False,
                                                     "reason": f"el candidato no pertenece al dominio: {e}"}
                return cand
            if zs != 0:
                cand.level = "refuted"
                cand.proof = {"f_candidate_is_exactly_zero": False,
                              "reason": "f(candidato) ≠ 0 exactamente: el candidato NO es solución "
                                        "(coincidencia numérica descartada)"}
                return cand
            proof = {"f_candidate_is_exactly_zero": True, "method": "evaluación exacta f(candidato)=0"}
            if bracket is not None:
                lo, hi = bracket
                iv = pool.iv(max(60, len(str(hi.denominator)) + 40))
                try:
                    enc = value.to_iv(iv.ctx)
                    proof["candidate_inside_solver_interval"] = bool(enc.a >= iv.from_fraction(lo).a
                                                                     and enc.b <= iv.from_fraction(hi).b)
                except Exception:
                    proof["candidate_inside_solver_interval"] = None
            u = unique_in_bracket()
            proof["unique_in_interval"] = u
            if u is True and proof.get("candidate_inside_solver_interval"):
                cand.level = "exact_proven"
            elif proof.get("candidate_inside_solver_interval") is False:
                cand.level = "numerical_evidence"
                proof["note"] = "es raíz exacta pero no cae en el intervalo calculado (no se afirma que sea esta raíz)"
            else:
                cand.level = "exact_verified"
                proof["note"] = ("f(candidato)=0 exacto; la unicidad de la raíz en el intervalo no está "
                                 "demostrada (función no polinómica)" if u is None else "")
            cand.proof = proof
            return cand
        return verify

    def _make_poly_context(self, bracket, ctx, ev):
        """Demuestra que el polinomio mínimo q (hallado numéricamente) tiene la raíz buscada:
        q | P (división exacta) y Sturm cuenta 1 raíz de q en el intervalo."""
        chain_all = ctx.get("chain_all")
        lo, hi = bracket

        def prove(coeffs):
            q = [Fraction(c) for c in reversed(coeffs)]
            if chain_all is None:
                return None
            # P = producto libre de cuadrados de todos los factores: q | P  y  raíz de q en (lo,hi]
            from .polynomial import sturm_chain as sc
            sqP = ctx["_sqprod"]
            _, rem = pdivmod(sqP, q)
            if rem:
                return {"status": "q no divide al polinomio de la ecuación: no demostrado"}
            cq = sc(q)
            nq = count_roots(cq, lo, hi)
            uniq = count_roots(chain_all, lo, hi) == 1
            if nq >= 1 and uniq:
                return {"status": "exact_proven",
                        "statement": "q divide exactamente al polinomio libre de cuadrados de la ecuación y "
                                     "la raíz es la única de éste en el intervalo (Sturm): q(raíz)=0 es EXACTO",
                        "minimal_polynomial_irreducibility": "no demostrada (salvo grado ≤ 3)"}
            return {"status": "numerical_evidence"}
        return prove
