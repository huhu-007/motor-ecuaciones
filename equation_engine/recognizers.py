"""ExactFormRecognizer, RadicalRecognizer, ConstantRecognizer y PolynomialAnalyzer.

Política anti falsos positivos
------------------------------
* Toda identificación numérica se acepta solo si  (nº de dígitos disponibles) excede lo que
  «cuesta» la hipótesis: se acota la probabilidad de coincidencia fortuita
  (≈ 0.3·Q²·2ε para racionales de denominador ≤ Q; relaciones enteras con |coef| ≤ M sobre m números
  exigen D ≥ m·log10(M) + margen).
* Un candidato es solo «numerical_evidence». Pasa a «exact_verified/exact_proven» ÚNICAMENTE si el
  verificador exacto (aritmética exacta + Sturm) lo demuestra.
"""
from __future__ import annotations
import itertools
import math
import time
from dataclasses import dataclass, field
from fractions import Fraction
from math import gcd, lcm, floor
from typing import Callable, Dict, List, Optional

from .exact import Rad, ConstMono
from .errors import NotExact, DomainError
from .numtheory import factorize, primes_upto, split_perfect_power, is_prime
from .polynomial import to_int_coeffs, poly_text, pdivmod, peval, deg as pdeg

HYP_MARGIN = 1e-9       # probabilidad global tolerada de coincidencia fortuita por familia de hipótesis
_SUPD = "⁰¹²³⁴⁵⁶⁷⁸⁹"


# ------------------------------------------------------------------ candidatos
@dataclass
class Candidate:
    kind: str
    pretty: str
    text: str
    value: object = None                     # Rad | ConstMono (forma exacta evaluable) o None
    node: object = None                      # AST (p. ej. radicales anidados)
    annihilator: Optional[List[int]] = None  # polinomio entero (mayor grado primero) con candidato como raíz
    evidence: dict = field(default_factory=dict)
    equivalent_forms: List[str] = field(default_factory=list)
    level: str = "numerical_evidence"        # | exact_verified | exact_proven
    proof: Optional[dict] = None
    numeric: Optional[Callable] = None       # ctx -> mpf  (valor numérico del candidato, para validar)

    def to_dict(self):
        return {"kind": self.kind, "exact_form": self.pretty, "exact_form_text": self.text,
                "level": self.level, "evidence": self.evidence, "proof": self.proof,
                "equivalent_forms": self.equivalent_forms, "annihilating_polynomial": self.annihilator}


# ------------------------------------------------------------------ formato
def root_symbol(m: int) -> str:
    return {2: "√", 3: "∛", 4: "∜"}.get(m) or ("".join(_SUPD[int(c)] for c in str(m)) + "√")


def _qtxt(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def mono_parts(mono):
    """(índice m, radicando b) de un monomio: prod p^e = b^(1/m)."""
    m = 1
    for _, e in mono:
        m = lcm(m, e.denominator)
    b = 1
    for p, e in mono:
        b *= p ** int(e * m)
    return m, b


def _term(c: Fraction, mono, pretty: bool) -> str:
    """Texto de |c| * monomio (sin signo)."""
    c = abs(c)
    if not mono:
        return _qtxt(c)
    m, b = mono_parts(mono)
    n, d = c.numerator, c.denominator
    if pretty:
        rad = root_symbol(m) + str(b)
        numer = ("" if n == 1 else str(n)) + rad
        return (f"({numer})/{d}" if n > 1 else f"{numer}/{d}") if d > 1 else numer
    rad = f"sqrt({b})" if m == 2 else f"root_{m}({b})"
    numer = (f"{n}*" if n != 1 else "") + rad
    return f"{numer}/{d}" if d > 1 else numer


def _compact_single(r: Rad, pretty: bool) -> Optional[str]:
    """Para un monomio c·b^(1/m) con c>0: la forma m√(c^m·b) (p. ej. ¹⁰⁰√(7/3)) si es mucho más corta."""
    (mono, c), = r.t.items()
    if not mono:
        return None
    m, b = mono_parts(mono)
    rr = (abs(c) ** m) * b
    inner = _qtxt(rr)
    if pretty:
        return ("-" if c < 0 else "") + root_symbol(m) + (inner if rr.denominator == 1 and len(inner) <= 3 else f"({inner})")
    rad = "sqrt" if m == 2 else f"root_{m}"
    return ("-" if c < 0 else "") + f"{rad}({inner})"


def rad_to_str(r: Rad, pretty=True) -> str:
    if r.is_zero():
        return "0"
    if len(r.t) == 1:
        alt = _compact_single(r, pretty)
        std = _single_std(r, pretty)
        if alt and len(alt) + 3 < len(std):
            return alt
        return std
    return _multi_str(r, pretty)


def _single_std(r, pretty):
    (mono, c), = r.t.items()
    return ("-" if c < 0 else "") + _term(c, mono, pretty)


def _multi_str(r: Rad, pretty=True) -> str:
    items = sorted(r.t.items(), key=lambda kv: (len(kv[0]) > 0, mono_parts(kv[0]) if kv[0] else (0, 0)))
    out = ""
    for i, (mono, c) in enumerate(items):
        t = _term(c, mono, pretty)
        if i == 0:
            out += ("-" if c < 0 else "") + t
        else:
            out += (" - " if c < 0 else " + ") + t
    return out


# ------------------------------------------------------------------ utilidades numéricas
def to_fraction(v) -> Fraction:
    """mpf -> Fraction exacta (conserva el signo)."""
    sign, man, exp, _bc = v._mpf_
    f = Fraction(int(man)) * (Fraction(2) ** int(exp))
    return -f if sign else f


def rational_fit(ctx, v, err, D=None, margin=HYP_MARGIN):
    """Fraction p/q con |v - p/q| <= 2·err y denominador ≤ Q, o None.
    Q = min( sqrt(margin/2err) , 10^(D/4) ).  El segundo tope es imprescindible: irracionales
    cuadráticos (p. ej. (√3+√2)^20) tienen cocientes parciales enormes y admiten convergentes de
    altura ~10^30 con error ~10^-70 POR ESTRUCTURA; limitando la altura a 10^(D/4) una
    coincidencia así exigiría cocientes parciales > 10^(D/2), imposible en estos casos.
    Devuelve (fracción, Q, cota de falso positivo) o None."""
    err = abs(err)
    if err <= 0:
        err = ctx.mpf(10) ** (-ctx.dps)
    Q = int(ctx.sqrt(margin / (2 * err))) if err < margin / 2 else 0
    if D is not None:
        Q = min(Q, 10 ** max(1, D // 4))
    if Q < 1:
        return None
    fr = to_fraction(v).limit_denominator(Q)
    if abs(v - ctx.mpf(fr.numerator) / fr.denominator) <= 2 * err:
        return fr, Q, float(0.6 * Q * Q * err)
    return None


# ------------------------------------------------------------------ reconocedores
class RadicalRecognizer:
    def __init__(self, max_index=100, max_poly_dps_margin=15):
        self.max_index = max_index

    # --- potencias racionales: x^n ∈ Q  (n = 2..100) ---
    def powers(self, ctx, x, D, deadline):
        hits, limited = {}, False
        ax = abs(x)
        eps = ctx.mpf(10) ** (-D)
        for n in range(2, self.max_index + 1):
            if time.monotonic() > deadline:
                break
            v = ax ** n
            err = n * ax ** (n - 1) * eps * 1.5
            fit = rational_fit(ctx, v, err, D)
            if fit is None:
                if err >= HYP_MARGIN / 8:
                    limited = True
                continue
            if fit[0] > 0:
                hits[n] = fit
        return hits, limited

    def monomial(self, x, hits) -> Optional[Candidate]:
        """x = a·ⁿ√b / c, radicales fraccionarios, productos y cocientes: todo es un monomio radical."""
        neg = x < 0
        for n in sorted(hits):
            fr, Q, fp = hits[n]
            try:
                r = Rad._factor_fraction_root(fr, n)
            except NotExact:
                continue
            if r.is_rational():                    # x^n racional y x racional: lo gestiona el caso racional
                continue
            r = -r if neg else r
            cand = Candidate("single_radical", rad_to_str(r), rad_to_str(r, False), r)
            cand.numeric = lambda cx, e=r: e.to_mpf(cx)
            cand.evidence = {"x_pow_n_rational": f"{fr.numerator}/{fr.denominator}", "n": n,
                             "coincidence_probability_bound": fp}
            p, q = fr.numerator, fr.denominator
            eq = {cand.pretty}
            rs = root_symbol(n)
            eq.add(("-" if neg else "") + f"{rs}{p}" + (f"/{rs}{q}" if q > 1 else ""))
            a = round(p ** (1 / n)) if p < 10 ** 15 else None
            if q > 1 and a and a ** n == p:
                eq.add(("-" if neg else "") + f"{a}/{rs}{q}")
            b = round(q ** (1 / n)) if q < 10 ** 15 else None
            if q > 1 and b and b ** n == q:
                eq.add(("-" if neg else "") + f"{rs}{p}/{b}")
            if q > 1:
                eq.add(("-" if neg else "") + f"{rs}({p}/{q})")
            cand.equivalent_forms = sorted(eq, key=lambda t: (len(t), t))
            best = cand.equivalent_forms[0]
            if len(best) < len(cand.pretty) and best.replace("(", "").replace(")", "").replace("/", "").lstrip("-").lstrip(rs).isdigit() is False:
                pass
            if len(cand.equivalent_forms[0]) < len(cand.pretty) - 8:
                cand.pretty = cand.equivalent_forms[0]
            return cand
        return None

    # --- sumas/diferencias de radicales (relaciones enteras PSLQ) ---
    @staticmethod
    def _families(D):
        """Cada familia = tupla de (base, N): generadores base^(1/N); la base del caso
        (1, √s) puede ser compuesta libre de cuadrados.  Orden: dimensión creciente."""
        sqfree = [s for s in range(2, 201) if all(e == 1 for e in factorize(s).values())]
        P7, P13, P31 = primes_upto(7), primes_upto(13), primes_upto(31)
        fams = []
        for s_ in sqfree:
            fams.append(((s_, 2),))                                           # dim 2
        for N in (3, 4, 5, 6, 7, 8, 9, 10, 12):
            for p in P13:
                fams.append(((p, N),))                                        # dim N
        for p, q in itertools.combinations(P31, 2):
            fams.append(((p, 2), (q, 2)))                                     # dim 4
        for a, b in ((2, 3), (2, 4), (3, 4), (2, 5), (3, 5), (2, 6)):         # índices mixtos
            for p in P7:
                for q in P7:
                    if p != q:
                        fams.append(((p, a), (q, b)))                         # dim a·b
        for trio in itertools.combinations(P13, 3):
            fams.append(tuple((p, 2) for p in trio))                          # dim 8
        for N in (3, 4):
            for p, q in itertools.combinations(P7, 2):
                fams.append(((p, N), (q, N)))                                 # dim N²
        fams.sort(key=lambda f: math.prod(n for _, n in f))
        return fams

    @staticmethod
    def _basis(spec):
        gens = [Rad._factor_fraction_root(Fraction(b), N) for b, N in spec]
        out = []
        for js in itertools.product(*[range(N) for _, N in spec]):
            e = Rad.rat(1)
            for g, j in zip(gens, js):
                if j:
                    e = e * g.pow_int(j)
            out.append(e)
        return out

    def combination(self, ctx, x, D, deadline, notes, min_dim=0, max_dim=10 ** 9,
                    degree=None) -> Optional[Candidate]:
        """degree: grado del polinomio mínimo de x si se conoce. Si x ∈ K y [K:Q] = m entonces
        grado(x) | m: solo se prueban familias cuya dimensión es múltiplo de ese grado."""
        skipped = 0
        for spec in self._families(D):
            dim = math.prod(n for _, n in spec)
            if not (min_dim <= dim <= max_dim) or (degree and dim % degree):
                continue
            if time.monotonic() > deadline:
                notes.append("presupuesto de tiempo agotado durante la búsqueda de combinaciones de radicales")
                return None
            basis = self._basis(spec)
            m = len(basis) + 1
            k = min(8, (D - 15) // m)
            if k < 3:
                skipped += 1
                continue
            M = 10 ** k
            vals = [b.to_mpf(ctx) for b in basis]
            rel = ctx.pslq([x] + vals, tol=ctx.mpf(10) ** (-(D - 4)), maxcoeff=M, maxsteps=20000)
            if not rel or rel[0] == 0:
                continue
            expr = Rad()
            for c, b in zip(rel[1:], basis):
                expr = expr + b.scale(Fraction(-c, rel[0]))
            if expr.is_zero() or abs(expr.to_mpf(ctx) - x) > ctx.mpf(10) ** (-(D - 3)):
                continue
            kind = "single_radical" if (len(expr.t) < 2 and not expr.is_rational()) else "radical_combination"
            c = Candidate(kind, rad_to_str(expr), rad_to_str(expr, False), expr)
            c.numeric = lambda cx, e=expr: e.to_mpf(cx)
            c.evidence = {"pslq_dimension": m, "max_coefficient_searched": M, "digits_used": D,
                          "coincidence_probability_bound": float(M) ** m * 10.0 ** (-D)}
            return c
        if skipped:
            notes.append(f"{skipped} familias de radicales no comprobadas: dígitos insuficientes")
        return None

    # --- radicales anidados ⁿ√(α + β√d) ---
    def nested(self, ctx, x, D, deadline, notes, max_n=12, degree=None) -> Optional[Candidate]:
        ax = abs(x)
        eps = ctx.mpf(10) ** (-D)
        for n in range(2, max_n + 1):
            if time.monotonic() > deadline:
                return None
            if degree and (2 * n) % degree:
                continue
            v = ax ** n
            Dv = D - int(math.log10(n * float(ax) ** (n - 1) + 1)) - 1
            k = (Dv - 15) // 3
            if k < 2:
                continue
            q = ctx.findpoly(v, 2, maxcoeff=10 ** k, tol=ctx.mpf(10) ** (-(Dv - 4)), maxsteps=20000)
            if not q or len(q) != 3 or q[0] == 0:
                continue
            A, B, C = q
            disc = B * B - 4 * A * C
            if disc <= 0:
                continue
            f, d = split_perfect_power(disc, 2)
            if d == 1:
                continue                          # v racional: ya cubierto por x^n ∈ Q
            alpha = Fraction(-B, 2 * A)
            beta_abs = Fraction(f, 2 * abs(A))
            for s in (1, -1):
                val = ctx.mpf(alpha.numerator) / alpha.denominator + s * ctx.mpf(beta_abs.numerator) / beta_abs.denominator * ctx.sqrt(d)
                if abs(val - v) < ctx.mpf(10) ** (-(Dv - 3)):
                    beta = s * beta_abs
                    break
            else:
                continue
            ann = [1] + [0] * (2 * n - 1) + [0]
            # y^{2n} - 2·alpha·y^n + (alpha² - beta²·d) = 0   (con denominadores limpiados)
            c2, c1, c0 = Fraction(1), -2 * alpha, alpha * alpha - beta * beta * d
            L = lcm(c1.denominator, c0.denominator)
            coeffs = [0] * (2 * n + 1)
            coeffs[0], coeffs[n], coeffs[2 * n] = int(c2 * L), int(c1 * L), int(c0 * L)   # grado creciente
            ann = list(reversed(coeffs))
            inner = f"{_qtxt(alpha)} {'+' if beta > 0 else '-'} " + (
                "" if abs(beta) == 1 else _qtxt(abs(beta)) + "") + f"√{d}"
            inner_t = f"{_qtxt(alpha)}{'+' if beta > 0 else '-'}" + ("" if abs(beta) == 1 else _qtxt(abs(beta)) + "*") + f"sqrt({d})"
            sgn = "-" if x < 0 else ""
            pretty = f"{sgn}{root_symbol(n)}({inner})"
            text = f"{sgn}" + (f"sqrt({inner_t})" if n == 2 else f"root_{n}({inner_t})")
            cand = Candidate("nested_radical", pretty, text, None)
            cand.numeric = lambda cx, n_=n, a_=alpha, b_=beta, d_=d, sg=(-1 if x < 0 else 1): sg * cx.root(
                cx.mpf(a_.numerator) / a_.denominator + cx.mpf(b_.numerator) / b_.denominator * cx.sqrt(d_), n_)
            cand.annihilator = ann if x > 0 else [c if (i % 2 == 0) else -c for i, c in enumerate(ann)][::1]
            cand.evidence = {"n": n, "x_pow_n_minimal_polynomial": [A, B, C], "digits_used": Dv,
                             "coincidence_probability_bound": float(10 ** k) ** 3 * 10.0 ** (-Dv)}
            # denesting para n = 2:  √(α+β√d) = √((α+γ)/2) ± √((α-γ)/2) si α²-β²d = γ² racional
            if n == 2:
                g2 = alpha * alpha - beta * beta * d
                if g2 >= 0:
                    gn, gd = split_perfect_power(g2.numerator, 2)[0], split_perfect_power(g2.denominator, 2)[0]
                    if gn ** 2 == g2.numerator and gd ** 2 == g2.denominator:
                        gam = Fraction(gn, gd)
                        u, w = (alpha + gam) / 2, (alpha - gam) / 2
                        if u >= 0 and w >= 0:
                            try:
                                E = Rad._factor_fraction_root(u, 2) + (Rad._factor_fraction_root(w, 2) if w > 0 else Rad()).scale(1 if beta > 0 else -1)
                                if (E * E - (Rad.rat(alpha) + Rad._factor_fraction_root(Fraction(d), 2).scale(beta))).is_zero():
                                    E = -E if x < 0 else E
                                    cand.equivalent_forms = [rad_to_str(E)]
                                    cand.value = E
                                    cand.pretty, cand.text = rad_to_str(E) + "", rad_to_str(E, False)
                                    cand.kind = "radical_combination"
                                    cand.equivalent_forms = [pretty]
                                    cand.numeric = lambda cx, e=E: e.to_mpf(cx)
                            except (NotExact, DomainError):
                                pass
            return cand
        return None


class ConstantRecognizer:
    """q·π, q·e, ln(q).  (π y e se tratan como constantes simbólicas)."""

    def recognize(self, ctx, x, D) -> List[Candidate]:
        out = []
        eps = ctx.mpf(10) ** (-D)
        if x == 0:
            return out
        for name, val, tag in (("π", ctx.pi, "pi"), ("e", ctx.e, "e")):
            r = x / val
            fit = rational_fit(ctx, r, eps / val * 1.5, D)
            if fit and fit[0] != 0:
                q = fit[0]
                pretty = (("-" if q < 0 else "") + ("" if abs(q) == 1 else _qtxt(abs(q).numerator and Fraction(abs(q).numerator, 1)) if False else ""))
                n_, d_ = abs(q).numerator, abs(q).denominator
                num = name if n_ == 1 else f"{n_}{name}"
                pretty = ("-" if q < 0 else "") + (num if d_ == 1 else (f"({num})/{d_}" if n_ > 1 else f"{num}/{d_}"))
                numt = tag if n_ == 1 else f"{n_}*{tag}"
                text = ("-" if q < 0 else "") + (numt if d_ == 1 else f"{numt}/{d_}")
                cm = ConstMono(q, 1 if tag == "pi" else 0, 0 if tag == "pi" else 1)
                c = Candidate("constant", pretty, text, cm)
                c.numeric = lambda cx, m_=cm: m_.to_mpf(cx)
                c.evidence = {"digits_used": D, "coincidence_probability_bound": fit[2]}
                out.append(c)
        if x > 0:
            lx = ctx.ln(x)
            fit = rational_fit(ctx, lx, eps / x * 1.5, D)
            if fit and fit[0] != 0 and fit[0] != 1:
                q = fit[0]
                cm = ConstMono(1, 0, q)
                c = Candidate("constant", f"e^{_qtxt(q) if q.denominator == 1 else '(' + _qtxt(q) + ')'}",
                              f"e^({_qtxt(q)})", cm)
                c.numeric = lambda cx, q_=q: cx.exp(cx.mpf(q_.numerator) / q_.denominator)
                c.evidence = {"digits_used": D, "coincidence_probability_bound": fit[2],
                              "note": "ln(x) es racional"}
                out.append(c)
        ex = ctx.exp(x)
        fit = rational_fit(ctx, ex, ex * eps * 1.5, D)
        if fit and fit[0] > 0 and fit[0] != 1:
            q = fit[0]
            c = Candidate("constant", f"ln({_qtxt(q)})", f"ln({_qtxt(q)})", None)
            c.numeric = lambda cx, q_=q: cx.ln(cx.mpf(q_.numerator) / q_.denominator)
            c.evidence = {"digits_used": D, "coincidence_probability_bound": fit[2],
                          "note": "exp(x) es racional"}
            out.append(c)
        return out


# ------------------------------------------------------------------ polinomio mínimo
class PolynomialAnalyzer:
    def __init__(self, max_degree=12):
        self.max_degree = max_degree

    @staticmethod
    def reliable_degree(D, min_k=2):
        """Mayor grado d con coeficientes ≤ 10^min_k comprobable con D dígitos (D ≥ 12 + k(d+1))."""
        return max(0, (D - 12) // min_k - 1)

    def analyze(self, ctx, x, D, deadline, notes) -> dict:
        res = {"found": False, "max_reliable_degree": self.reliable_degree(D),
               "max_degree_tried": 0, "precision_sufficient_for_requested_degree": None}
        dmax = min(self.max_degree, res["max_reliable_degree"])
        res["precision_sufficient_for_requested_degree"] = self.max_degree <= res["max_reliable_degree"]
        for d in range(1, dmax + 1):
            if time.monotonic() > deadline:
                notes.append("presupuesto de tiempo agotado buscando el polinomio mínimo")
                break
            k = (D - 12) // (d + 1)
            q = ctx.findpoly(x, d, maxcoeff=10 ** k, tol=ctx.mpf(10) ** (-(D - 4)), maxsteps=50000)
            res["max_degree_tried"] = d
            if q:
                coeffs = [int(c) for c in q]
                while coeffs and coeffs[0] == 0:
                    coeffs.pop(0)
                # Criterio honesto: si x tiene error <= 10^-D entonces |P(x)| <= 10^-D·|P'(x)|.
                # (evita falsos positivos por cuasi-raíces múltiples, donde P y P' son minúsculos)
                val, dval = ctx.polyval(coeffs, x, derivative=True)
                if dval != 0 and abs(val) <= 2 * ctx.mpf(10) ** (-D) * abs(dval):
                    res.update(found=True, degree=len(coeffs) - 1, coefficients=coeffs,
                               polynomial=poly_text([Fraction(c) for c in reversed(coeffs)]) + " = 0",
                               coefficient_bound_searched=10 ** k)
                    res.update(self.radical_remarks(coeffs))
                    return res
        if res["max_degree_tried"] < self.max_degree:
            res["note"] = (f"con {D} dígitos solo se puede descartar hasta grado {res['max_reliable_degree']} "
                           f"con coeficientes ≤ 100; se pidió {self.max_degree}")
        return res

    @staticmethod
    def has_rational_root(coeffs) -> bool:
        """Test exacto (teorema de la raíz racional) para polinomio entero, mayor grado primero."""
        c = list(reversed(coeffs))            # grado creciente
        if c[0] == 0:
            return True
        divs = lambda n: [d for d in range(1, abs(n) + 1) if n % d == 0] if abs(n) < 10 ** 7 else None
        a0, an = divs(c[0]), divs(c[-1])
        if a0 is None or an is None:
            return None
        P = [Fraction(v) for v in c]
        for p in a0:
            for q in an:
                for s in (1, -1):
                    if peval(P, Fraction(s * p, q)) == 0:
                        return True
        return False

    @staticmethod
    def radical_remarks(coeffs) -> dict:
        d = len(coeffs) - 1
        out = {}
        if d == 1:
            out["radical_remark"] = "raíz racional"
        elif d == 2:
            out["radical_remark"] = "raíz de un polinomio cuadrático: se expresa con una raíz cuadrada (fórmula general)"
        else:
            rr = PolynomialAnalyzer.has_rational_root(coeffs) if d <= 3 else None
            if d == 3 and rr is False:
                out["irreducible_over_Q"] = True
                a, b, c_, dd = coeffs
                disc = 18 * a * b * c_ * dd - 4 * b ** 3 * dd + b * b * c_ * c_ - 4 * a * c_ ** 3 - 27 * a * a * dd * dd
                out["discriminant"] = disc
                out["radical_remark"] = (
                    "cúbica irreducible con 3 raíces reales: NO se puede escribir con radicales reales (casus irreducibilis)"
                    if disc > 0 else
                    "cúbica irreducible con una sola raíz real: se puede escribir con radicales (Cardano), no construido aquí")
            else:
                out["irreducible_over_Q"] = "no demostrado"
                if d & (d - 1):
                    out["radical_remark"] = ("Si este polinomio es el mínimo, el grado no es potencia de 2: es IMPOSIBLE "
                                             "escribirlo solo con raíces cuadradas")
                elif d == 4:
                    out["radical_remark"] = "grado 4: resoluble por radicales (Ferrari); no construido aquí"
                else:
                    out["radical_remark"] = "la resolubilidad por radicales depende del grupo de Galois (no calculado)"
        return out


# ------------------------------------------------------------------ fachada
@dataclass
class RecognitionReport:
    digits_available: int
    candidates: List[Candidate] = field(default_factory=list)
    best: Optional[Candidate] = None
    rational: Optional[dict] = None
    square_rational: Optional[bool] = None
    rational_powers: Dict[int, str] = field(default_factory=dict)
    polynomial: Optional[dict] = None
    precision: dict = field(default_factory=dict)
    notes: List[str] = field(default_factory=list)

    def to_dict(self):
        return {"digits_available": self.digits_available,
                "best": self.best.to_dict() if self.best else None,
                "candidates": [c.to_dict() for c in self.candidates],
                "rational": self.rational, "square_is_rational": self.square_rational,
                "rational_powers": self.rational_powers, "polynomial": self.polynomial,
                "precision": self.precision, "notes": self.notes,
                "statement": self.statement()}

    def statement(self):
        if self.best is None:
            return "No se ha identificado ninguna forma exacta: solo hay aproximación numérica."
        lv = {"exact_proven": "FORMA EXACTA DEMOSTRADA",
              "exact_verified": "FORMA EXACTA VERIFICADA (f(candidato)=0 exacto; unicidad en el intervalo no demostrada)",
              "numerical_evidence": "EVIDENCIA NUMÉRICA (no es una demostración)"}[self.best.level]
        return f"{lv}: {self.best.pretty}"


class ExactFormRecognizer:
    def __init__(self, max_radical_index=100, max_poly_degree=12, budget_seconds=20.0):
        self.radicals = RadicalRecognizer(max_radical_index)
        self.constants = ConstantRecognizer()
        self.poly = PolynomialAnalyzer(max_poly_degree)
        self.budget = budget_seconds

    def recognize(self, x, digits: int, verifier: Optional[Callable] = None,
                  poly_context: Optional[Callable] = None) -> RecognitionReport:
        """x: Fraction | mpf | str con error absoluto <= 10^-digits.
        verifier(candidate)->candidate actualiza level/proof (demostración exacta).
        poly_context(coeffs)->dict con la prueba exacta del polinomio mínimo.
        Trabaja en un contexto mpmath PROPIO (sin efectos laterales sobre el llamador)."""
        from mpmath.ctx_mp import MPContext
        t0 = time.monotonic()
        deadline = t0 + self.budget
        D = int(digits)
        ctx = MPContext()
        ctx.dps = D + 40
        if isinstance(x, Fraction):
            x = ctx.mpf(x.numerator) / x.denominator
        elif isinstance(x, str):
            x = ctx.mpf(x)
        else:
            x = ctx.mpf(to_fraction(x).numerator) / to_fraction(x).denominator if hasattr(x, "man_exp") else ctx.mpf(x)
        rep = RecognitionReport(D)
        H = max(5, D // 8)                      # dígitos de reserva para validar (holdout)
        Did = D - H
        eps = ctx.mpf(10) ** (-Did)
        rep.precision = {"digits_available": D, "digits_used_to_identify": Did,
                         "holdout_digits": H,
                         "reliable_polynomial_degree": PolynomialAnalyzer.reliable_degree(Did)}
        hold_tol = ctx.mpf(10) ** (-(D - 1))

        def accept(c: Candidate) -> bool:
            # validación con los dígitos de reserva: un candidato verdadero sigue coincidiendo
            if c.numeric is not None:
                try:
                    ok = abs(c.numeric(ctx) - x) <= hold_tol
                except Exception:
                    ok = False
                if not ok:
                    rep.notes.append(f"candidato {c.pretty} descartado: no coincide con los {H} dígitos de reserva")
                    return False
                c.evidence["holdout_digits_confirmed"] = H
            if verifier:
                c = verifier(c)
            rep.candidates.append(c)
            if rep.best is None or _rank(c) > _rank(rep.best):
                rep.best = c
            return c.level != "numerical_evidence"

        # 1) racional
        fit = rational_fit(ctx, x, eps * 1.5, Did)
        if fit:
            fr, Q, fp = fit
            c = Candidate("rational", _qtxt(fr), _qtxt(fr), Rad.rat(fr))
            c.numeric = lambda cx, f_=fr: cx.mpf(f_.numerator) / f_.denominator
            c.evidence = {"digits_used": Did, "coincidence_probability_bound": fp}
            if accept(c):
                pass
            if rep.best is not None:
                rep.rational = {"value": _qtxt(fr), "coincidence_probability_bound": fp}
                rep.square_rational = True
                rep.precision["seconds"] = round(time.monotonic() - t0, 3)
                return rep
        # 2) potencias racionales x^n (n = 2..100)
        hits, limited = self.radicals.powers(ctx, x, Did, deadline)
        rep.rational_powers = {n: _qtxt(h[0]) for n, h in hits.items()}
        rep.square_rational = 2 in hits
        if limited:
            rep.notes.append("con esta precisión no se puede comprobar si algunas potencias altas son racionales")
        done = False
        # 3) constantes π, e, ln
        for c in self.constants.recognize(ctx, x, Did):
            done = accept(c) or done
        # 4) polinomio mínimo (barato): da el GRADO, que poda las familias de radicales
        poly = None
        if not done and rep.best is None:
            poly = self.poly.analyze(ctx, x, Did, deadline, rep.notes)
            if poly.get("found") and poly_context:
                poly["proof"] = poly_context(poly["coefficients"])
            rep.polynomial = poly
        degree = poly.get("degree") if poly and poly.get("found") else None
        algebraic_known_absent = bool(poly) and not poly.get("found") and \
            poly.get("max_degree_tried", 0) >= self.poly.max_degree
        if not done and rep.best is None and not (poly and poly.get("found") is False and algebraic_known_absent
                                                   and degree is None and False):
            # 5) un solo radical (coef·radical, fracciones, productos, cocientes)
            c = self.radicals.monomial(x, hits)
            if c:
                done = accept(c)
            # 6) sumas de radicales pequeñas (dim <= 4), 7) anidados, 8) sumas mayores
            if not done and rep.best is None:
                c = self.radicals.combination(ctx, x, Did, deadline, rep.notes, max_dim=4, degree=degree)
                if c:
                    done = accept(c)
            if not done and rep.best is None:
                c = self.radicals.nested(ctx, x, Did, deadline, rep.notes, degree=degree)
                if c:
                    done = accept(c)
            if not done and rep.best is None:
                c = self.radicals.combination(ctx, x, Did, deadline, rep.notes, min_dim=5, degree=degree)
                if c:
                    done = accept(c)
        rep.precision["seconds"] = round(time.monotonic() - t0, 3)
        if rep.best is None:
            rep.notes.append("Sin forma exacta identificada con los dígitos disponibles.")
        return rep


def _rank(c: Candidate):
    return {"exact_proven": 3, "exact_verified": 2, "numerical_evidence": 1}[c.level]
