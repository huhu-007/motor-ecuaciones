"""RootFinder: búsqueda automática de intervalos (cambio de signo) + Bolzano.

Dos rutas:
  * scan_generic    : cualquier combinación de funciones (polinomios, racionales, exp, log, trig, raíces)
  * scan_polynomial : f racional -> numerador polinómico exacto; Yun (multiplicidades) + Sturm
                      (garantiza no perder raíces ni contar de más).
"""
from __future__ import annotations
from dataclasses import dataclass, field
from fractions import Fraction
from typing import Callable, List, Optional

from .errors import DomainError, NotExact
from .exact import ExactBackend, Rad
from .bolzano import BolzanoAnalyzer, BolzanoReport
from .polynomial import (Poly, yun, sturm_chain, count_roots, count_all_real, cauchy_bound,
                         isolate, peval, psign, poly_text, deg)


@dataclass
class SearchConfig:
    limit: Fraction = Fraction(100)        # se buscan raíces en [-limit, limit]
    step: Fraction = Fraction(1, 10)       # paso de la rejilla
    scan_dps: int = 30
    guard_steps: int = 40                  # bisecciones de la comprobación anti-polo
    probe_depth: int = 12                  # sondeo cerca de fronteras de dominio
    touch_threshold_exp: int = 20          # |f|min < 10^-20  => candidata de raíz tangente
    max_touch_refinements: int = 200


@dataclass
class RootInterval:
    lo: Fraction
    hi: Fraction
    sign_lo: int
    sign_hi: int
    origin: str                              # 'grid' | 'boundary_probe' | 'exact_grid_zero' | 'sturm' | ...
    report: Optional[BolzanoReport] = None
    exact_point: bool = False
    multiplicity: Optional[int] = None
    factor: Optional[Poly] = None
    sign_fn: Optional[Callable] = None       # signo exacto/riguroso para refinar
    touching: bool = False


@dataclass
class ScanResult:
    intervals: List[RootInterval] = field(default_factory=list)
    rejected: List[dict] = field(default_factory=list)
    touch_candidates: List[dict] = field(default_factory=list)
    undefined_points: int = 0
    total_points: int = 0
    notes: List[str] = field(default_factory=list)


def _sgn(v):
    return 1 if v > 0 else (-1 if v < 0 else 0)


class RootFinder:
    def __init__(self, pool, config: SearchConfig = None):
        self.pool, self.cfg = pool, config or SearchConfig()
        self.bolzano = BolzanoAnalyzer(pool, dps=self.cfg.scan_dps)

    # ------------------------------------------------------------------ utilidades
    def _grid(self):
        n = int(self.cfg.limit / self.cfg.step)
        return [i * self.cfg.step for i in range(-n, n + 1)]

    def _mp_value(self, evaluator, x, dps=None):
        be = self.pool.mp(dps or self.cfg.scan_dps)
        try:
            return evaluator(be, be.from_fraction(x))
        except (DomainError, ZeroDivisionError, OverflowError, ValueError):
            return None

    def exact_zero(self, evaluator, x: Fraction):
        """True / False / None (desconocido) : ¿f(x) = 0 EXACTAMENTE?"""
        try:
            return evaluator(ExactBackend(), Rad.rat(x)).sign() == 0
        except (NotExact, DomainError, ZeroDivisionError):
            return None

    def looks_like_root(self, evaluator, a, fa, b, fb):
        """Protección anti-polo (original parece_raiz): tras N bisecciones |f| debe decrecer."""
        inicial = min(abs(fa), abs(fb))
        for _ in range(self.cfg.guard_steps):
            m = (a + b) / 2
            fm = self._mp_value(evaluator, m, self.cfg.scan_dps + 20)
            if fm is None:
                return False
            if fm == 0:
                return True
            if fa * fm < 0:
                b, fb = m, fm
            else:
                a, fa = m, fm
        return max(abs(fa), abs(fb)) < inicial

    # ------------------------------------------------------------------ ruta general
    def scan_generic(self, evaluator, sign_fn) -> ScanResult:
        res = ScanResult()
        xs = self._grid()
        vals = [self._mp_value(evaluator, x) for x in xs]
        res.total_points = len(xs)
        res.undefined_points = sum(v is None for v in vals)
        eps = Fraction(1, 10 ** 12)
        zero_flags = [False] * len(xs)
        for i, (x, v) in enumerate(zip(xs, vals)):
            if v is not None and abs(v) < 1e-12:
                ez = self.exact_zero(evaluator, x)
                if ez:
                    zero_flags[i] = True
                    res.intervals.append(RootInterval(x, x, 0, 0, "exact_grid_zero",
                                                      exact_point=True, sign_fn=sign_fn))
                elif ez is None:
                    res.touch_candidates.append({"x": float(x), "reason": "f(x) ≈ 0 en un punto de la rejilla; sin prueba exacta",
                                                 "abs_f": float(abs(v)), "digits": 10})
        # --- cambios de signo entre puntos consecutivos ---
        for i in range(len(xs) - 1):
            v0, v1 = vals[i], vals[i + 1]
            if v0 is None or v1 is None or zero_flags[i] or zero_flags[i + 1]:
                continue
            if _sgn(v0) * _sgn(v1) < 0:
                self._consider(evaluator, sign_fn, xs[i], v0, xs[i + 1], v1, "grid", res)
        # --- sondeo junto a las fronteras del dominio ---
        for i in range(len(xs) - 1):
            if (vals[i] is None) != (vals[i + 1] is None):
                self._probe_boundary(evaluator, sign_fn, xs, vals, i, res)
        # --- raíces tangentes (multiplicidad par): sin cambio de signo ---
        self._touching(evaluator, xs, vals, zero_flags, res)
        res.intervals.sort(key=lambda r: (r.lo, r.hi))
        dedup = []
        for r in res.intervals:
            if dedup and dedup[-1].lo == r.lo and dedup[-1].hi == r.hi:
                continue
            dedup.append(r)
        res.intervals = dedup
        return res

    def _consider(self, evaluator, sign_fn, a, fa, b, fb, origin, res):
        guard = self.looks_like_root(evaluator, a, fa, b, fb)
        rep = self.bolzano.analyze(evaluator, a, b, _sgn(fa), _sgn(fb), guard_ok=guard)
        if rep.applicable:
            res.intervals.append(RootInterval(a, b, _sgn(fa), _sgn(fb), origin, rep, sign_fn=sign_fn))
        else:
            res.rejected.append({"interval": [float(a), float(b)],
                                 "reason": ("cambio de signo por discontinuidad (polo / función no definida): no es una raíz"
                                            if not guard or rep.continuity != "certified" else rep.reason),
                                 "continuity": rep.continuity, "detail": rep.reason})

    def _defined(self, evaluator, x):
        return self._mp_value(evaluator, x) is not None

    def _probe_boundary(self, evaluator, sign_fn, xs, vals, i, res):
        if vals[i] is not None:
            p0, u = xs[i], xs[i + 1]
        else:
            p0, u = xs[i + 1], xs[i]
        d = p0                               # extremo definido; u es el indefinido
        for _ in range(45):                  # localiza la frontera del dominio
            m = (d + u) / 2
            if self._defined(evaluator, m):
                d = m
            else:
                u = m
        beta = d
        pts = [p0]
        for k in range(1, self.cfg.probe_depth + 1):
            p = beta + (p0 - beta) / Fraction(10) ** k
            pts.append(p)
        pts.append(beta)
        prev = None
        for p in pts:
            v = self._mp_value(evaluator, p, self.cfg.scan_dps + 20)
            cur = (p, v) if v is not None else None
            if prev and cur and _sgn(prev[1]) * _sgn(cur[1]) < 0:
                a, fa, b, fb = (prev[0], prev[1], cur[0], cur[1]) if prev[0] < cur[0] else (cur[0], cur[1], prev[0], prev[1])
                self._consider(evaluator, sign_fn, a, fa, b, fb, "boundary_probe", res)
            if cur:
                prev = cur

    def _touching(self, evaluator, xs, vals, zero_flags, res):
        cand = 0
        for i in range(1, len(xs) - 1):
            v0, v1, v2 = vals[i - 1], vals[i], vals[i + 1]
            if None in (v0, v1, v2) or zero_flags[i]:
                continue
            if not (_sgn(v0) == _sgn(v1) == _sgn(v2) != 0):
                continue
            if not (abs(v1) <= abs(v0) and abs(v1) <= abs(v2)):
                continue
            if self.bolzano.excludes_zero(evaluator, xs[i - 1], xs[i + 1]):
                continue                     # PROBADO: f no se anula aquí
            cand += 1
            if cand > self.cfg.max_touch_refinements:
                res.notes.append("demasiados mínimos locales de |f|; no se refinaron todos")
                break
            x, fmin = self._minimize_abs(evaluator, xs[i - 1], xs[i + 1])
            if fmin is not None and fmin < Fraction(1, 10 ** self.cfg.touch_threshold_exp):
                import math
                digits = max(6, int(-math.log10(float(fmin) + 1e-300) / 2) - 2)
                res.touch_candidates.append({"x": x, "abs_f": float(fmin), "digits": digits,
                                             "reason": "mínimo de |f| ≈ 0 sin cambio de signo (posible raíz de multiplicidad par)"})

    def _minimize_abs(self, evaluator, lo, hi):
        be = self.pool.mp(50)
        a, b = be.ctx.mpf(lo.numerator) / lo.denominator, be.ctx.mpf(hi.numerator) / hi.denominator
        g = (be.ctx.sqrt(5) - 1) / 2

        def f(t):
            try:
                return abs(evaluator(be, t))
            except (DomainError, ZeroDivisionError, OverflowError, ValueError):
                return be.ctx.inf
        c, d = b - g * (b - a), a + g * (b - a)
        fc, fd = f(c), f(d)
        for _ in range(110):
            if fc < fd:
                b, d, fd = d, c, fc
                c = b - g * (b - a)
                fc = f(c)
            else:
                a, c, fc = c, d, fd
                d = a + g * (b - a)
                fd = f(d)
        x = (a + b) / 2
        fx = f(x)
        if fx == be.ctx.inf:
            return None, None
        return x, Fraction(str(be.ctx.nstr(fx, 5))) if fx != 0 else Fraction(0)

    def prove_no_root_in_range(self, evaluator):
        """True si, en TODO el rango, f está definida y se prueba (intervalos) que no se anula."""
        xs = self._grid()
        for a, b in zip(xs, xs[1:]):
            if not self.bolzano.excludes_zero(evaluator, a, b, depth=5):
                return False
        return True

    # ------------------------------------------------------------------ ruta polinómica
    def scan_polynomial(self, num: Poly) -> dict:
        """Devuelve dict(res=ScanResult, total_real=int, outside=int, factors=[...])."""
        res = ScanResult()
        L = self.cfg.limit
        eps = Fraction(1, 10 ** 6)
        lo0 = -L - eps
        while any(peval(num, lo0) == 0 for _ in (0,)):
            lo0 -= eps
        factors, total, inside = [], 0, 0
        for fac, mult in yun(num):
            chain = sturm_chain(fac)
            tot = count_all_real(chain)
            ins = count_roots(chain, lo0, L)
            factors.append({"factor": poly_text(fac), "multiplicity": mult, "real_roots": tot,
                            "real_roots_in_range": ins})
            total += tot
            inside += ins
            sign_fn = (lambda p: (lambda x: psign(p, x)))(fac)
            found = self._grid_polynomial(fac, mult, sign_fn)
            if len(found) != ins:              # Sturm detecta raíces que la rejilla no separó
                found = []
                for lo, hi in isolate(chain, lo0, L):
                    if peval(fac, hi) == 0:
                        found.append(RootInterval(hi, hi, 0, 0, "sturm", exact_point=True,
                                                  multiplicity=mult, factor=fac, sign_fn=sign_fn))
                    else:
                        found.append(RootInterval(lo, hi, psign(fac, lo), psign(fac, hi), "sturm",
                                                  multiplicity=mult, factor=fac, sign_fn=sign_fn))
                res.notes.append("la rejilla no separaba todas las raíces; se completó con aislamiento de Sturm")
            for r in found:
                if not r.exact_point:
                    r.report = BolzanoAnalyzer.polynomial_report(r.lo, r.hi, r.sign_lo, r.sign_hi, poly_text(fac))
                else:
                    r.report = BolzanoAnalyzer.polynomial_report(r.lo, r.hi, 0, 0, poly_text(fac))
                    r.report.product = "=0"
                    r.report.reason = "raíz exacta en un punto de la rejilla (f(x)=0 exacto)"
                res.intervals.append(r)
        res.intervals.sort(key=lambda r: (r.lo, r.hi))
        return {"res": res, "total_real": total, "outside": total - inside, "factors": factors}

    def _grid_polynomial(self, fac, mult, sign_fn):
        xs = self._grid()
        sg = [psign(fac, x) for x in xs]
        out = []
        for i, x in enumerate(xs):
            if sg[i] == 0:
                out.append(RootInterval(x, x, 0, 0, "exact_grid_zero", exact_point=True,
                                        multiplicity=mult, factor=fac, sign_fn=sign_fn))
            elif i + 1 < len(xs) and sg[i + 1] != 0 and sg[i] * sg[i + 1] < 0:
                out.append(RootInterval(x, xs[i + 1], sg[i], sg[i + 1], "grid", multiplicity=mult,
                                        factor=fac, sign_fn=sign_fn))
        return out
