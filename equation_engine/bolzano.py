"""BolzanoAnalyzer: comprobación rigurosa de las hipótesis del teorema de Bolzano / TVI."""
from __future__ import annotations
from dataclasses import dataclass, field, asdict
from fractions import Fraction
from typing import Optional, Tuple

from .errors import DomainError, Uncertified, NotExact
from .nodes import to_text

BOLZANO_NO_APLICABLE = "No se puede aplicar el teorema de Bolzano en este tipo de ecuaciones."


@dataclass
class BolzanoReport:
    a: str
    b: str
    f_a_defined: bool
    f_b_defined: bool
    continuity: str            # 'certified' | 'polynomial' | 'not_certified' | 'undefined'
    sign_a: Optional[int]
    sign_b: Optional[int]
    product: Optional[str]     # '<0' | '=0' | '>0'
    applicable: bool
    strict: bool
    reason: str
    applied_to: str = "f"      # 'f' o la descripción del factor libre de cuadrados usado

    def to_dict(self):
        return asdict(self)


class BolzanoAnalyzer:
    def __init__(self, pool, max_depth=8, dps=30):
        self.pool, self.max_depth, self.dps = pool, max_depth, dps

    # ---------- continuidad en [a,b] ----------
    def certify_continuity(self, evaluator, a: Fraction, b: Fraction) -> Tuple[str, str]:
        """f es composición de funciones elementales (continuas donde están definidas);
        por tanto es continua en [a,b] si TODAS las restricciones de dominio se cumplen en
        todo [a,b]. Se certifica con aritmética de intervalos y subdivisión adaptativa."""
        be = self.pool.iv(self.dps)
        stack = [(a, b, 0)]
        while stack:
            lo, hi, d = stack.pop()
            try:
                evaluator(be, be.interval(lo, hi))
            except Uncertified as e:
                if d >= self.max_depth:
                    return "not_certified", f"no se pudo garantizar el dominio en [{float(lo):g}, {float(hi):g}]: {e}"
                m = (lo + hi) / 2
                stack += [(lo, m, d + 1), (m, hi, d + 1)]
            except DomainError as e:
                return "undefined", str(e)
            except (OverflowError, ZeroDivisionError, ValueError) as e:
                return "not_certified", f"desbordamiento/indeterminación: {e}"
        return "certified", "todas las restricciones de dominio se cumplen en [a,b] (intervalos)"

    def excludes_zero(self, evaluator, a: Fraction, b: Fraction, depth=4) -> Optional[bool]:
        """True si se PRUEBA que f no se anula en [a,b]; None si no se puede decidir."""
        be = self.pool.iv(self.dps)
        stack = [(a, b, 0)]
        while stack:
            lo, hi, d = stack.pop()
            try:
                v = evaluator(be, be.interval(lo, hi))
            except DomainError:
                return None
            if be.sign(v) is None:
                if d >= depth:
                    return None
                m = (lo + hi) / 2
                stack += [(lo, m, d + 1), (m, hi, d + 1)]
        return True

    # ---------- informe completo ----------
    def analyze(self, evaluator, a, b, sa, sb, guard_ok=True, applied_to="f"):
        defined = sa is not None and sb is not None
        if not defined:
            return BolzanoReport(str(a), str(b), sa is not None, sb is not None, "undefined",
                                 sa, sb, None, False, False, "f no está definida en un extremo")
        prod = "<0" if sa * sb < 0 else ("=0" if sa * sb == 0 else ">0")
        cont, why = self.certify_continuity(evaluator, a, b)
        ok = cont == "certified" and prod in ("<0", "=0") and guard_ok
        reason = why
        if cont == "certified" and prod == ">0":
            reason = "f(a)·f(b) > 0: no hay cambio de signo"
        elif cont == "certified" and not guard_ok:
            reason = "falló la comprobación de aproximación a cero (posible polo)"
        return BolzanoReport(str(a), str(b), True, True, cont, sa, sb, prod, ok, prod == "<0",
                             reason if not ok else
                             ("f(a)·f(b) < 0 y f continua en [a,b]: existe al menos una raíz (Bolzano)"
                              if prod == "<0" else "f(a)·f(b) = 0: extremo raíz"),
                             applied_to)

    @staticmethod
    def polynomial_report(a, b, sa, sb, factor_text):
        return BolzanoReport(str(a), str(b), True, True, "polynomial", sa, sb,
                             "<0" if sa * sb < 0 else "=0", True, sa * sb < 0,
                             "polinomio: continuo en ℝ; f(a)·f(b) < 0 ⇒ raíz (Bolzano)", factor_text)
