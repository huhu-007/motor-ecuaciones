"""Backends numéricos: MPBackend (alta precisión) e IntervalBackend (rigor)."""
from __future__ import annotations
from fractions import Fraction
from mpmath.ctx_mp import MPContext
from mpmath.ctx_iv import MPIntervalContext

from .errors import DomainError, Uncertified


class MPBackend:
    kind = "mp"

    def __init__(self, dps=30):
        self.ctx = MPContext()
        self.ctx.dps = dps

    def from_fraction(self, q):
        return self.ctx.mpf(q.numerator) / q.denominator

    def const(self, name):
        return +(self.ctx.pi if name == "pi" else self.ctx.e)

    neg = lambda self, a: -a
    add = lambda self, a, b: a + b
    sub = lambda self, a, b: a - b
    mul = lambda self, a, b: a * b
    div = lambda self, a, b: a / b
    abs = lambda self, a: abs(a)

    def pow_int(self, a, n):
        return a ** n

    def pow_general(self, a, b):
        return self.ctx.exp(b * self.ctx.log(a))

    def root(self, a, n):
        if a == 0:
            return a
        if a < 0:
            if n % 2 == 0:
                raise DomainError("raíz par de número negativo")
            return -self.ctx.root(-a, n)
        return self.ctx.root(a, n)

    def sin(self, a): return self.ctx.sin(a)
    def cos(self, a): return self.ctx.cos(a)
    def exp(self, a): return self.ctx.exp(a)
    def log(self, a): return self.ctx.log(a)
    def log_base(self, v, b): return self.ctx.log(v) / self.ctx.log(b)

    def require_nonzero(self, v, msg="división por cero"):
        if v == 0:
            raise DomainError(msg)

    def require_positive(self, v, msg="se requiere valor > 0"):
        if not v > 0:
            raise DomainError(msg)

    def require_nonnegative(self, v, msg="se requiere valor >= 0"):
        if v < 0:
            raise DomainError(msg)

    def sign(self, v):
        return 1 if v > 0 else (-1 if v < 0 else 0)

    def is_zero(self, v):
        return v == 0

    def to_mpf(self, v, ctx=None):
        return v


class IntervalBackend:
    """Aritmética de intervalos de mpmath: los resultados son ENCERRAMIENTOS rigurosos."""
    kind = "iv"

    def __init__(self, dps=30):
        self.ctx = MPIntervalContext()
        self.ctx.dps = dps

    def from_fraction(self, q):
        return self.ctx.mpf(q.numerator) / self.ctx.mpf(q.denominator)

    def interval(self, lo: Fraction, hi: Fraction):
        a, b = self.from_fraction(lo), self.from_fraction(hi)
        return self.ctx.mpf([a.a, b.b])

    def const(self, name):
        return self.ctx.pi if name == "pi" else self.ctx.e

    neg = lambda self, a: -a
    add = lambda self, a, b: a + b
    sub = lambda self, a, b: a - b
    mul = lambda self, a, b: a * b
    abs = lambda self, a: abs(a)

    def div(self, a, b):
        if b.a <= 0 <= b.b:
            raise Uncertified("denominador puede anularse")
        return a / b

    def pow_int(self, a, n):
        return a ** n if n >= 0 else 1 / (a ** (-n))

    def pow_general(self, a, b):
        return self.ctx.exp(b * self.ctx.log(a))

    def _pos_root(self, lo, hi, n):          # 0 <= lo <= hi
        c = self.ctx
        up = c.exp(c.log(c.mpf(hi)) / n) if hi > 0 else c.mpf(0)
        dn = c.exp(c.log(c.mpf(lo)) / n) if lo > 0 else c.mpf(0)
        return c.mpf([dn.a, up.b])

    def root(self, a, n):
        lo, hi = a.a, a.b
        if lo >= 0:
            return self._pos_root(lo, hi, n)
        if n % 2 == 0:
            raise Uncertified("raíz par de un intervalo con negativos")
        if hi <= 0:
            return -self._pos_root(-hi, -lo, n)
        r1, r2 = self._pos_root(0, -lo, n), self._pos_root(0, hi, n)
        return self.ctx.mpf([-r1.b, r2.b])

    def sin(self, a): return self.ctx.sin(a)
    def cos(self, a): return self.ctx.cos(a)
    def exp(self, a): return self.ctx.exp(a)
    def log(self, a): return self.ctx.log(a)
    def log_base(self, v, b): return self.ctx.log(v) / self.ctx.log(b)

    def _viol(self, v, cond_ok, cond_bad, msg):
        if cond_ok:
            return
        raise (DomainError if cond_bad else Uncertified)(msg)

    def require_nonzero(self, v, msg="división por cero"):
        self._viol(v, v.a > 0 or v.b < 0, v.a == 0 and v.b == 0, msg)

    def require_positive(self, v, msg="se requiere valor > 0"):
        self._viol(v, v.a > 0, v.b <= 0, msg)

    def require_nonnegative(self, v, msg="se requiere valor >= 0"):
        self._viol(v, v.a >= 0, v.b < 0, msg)

    def sign(self, v):
        """+1/-1 si el encierro lo garantiza; None si contiene al 0."""
        if v.a > 0:
            return 1
        if v.b < 0:
            return -1
        return None

    def is_zero(self, v):
        return v.a == 0 and v.b == 0

    def to_mpf(self, v, ctx=None):
        return v.mid


class BackendPool:
    """Cachea contextos de mpmath por precisión (crear contextos es costoso).
    Cada instancia del motor tiene su propio pool: no hay estado global compartido."""

    def __init__(self):
        self._mp, self._iv = {}, {}

    def mp(self, dps=30) -> MPBackend:
        if dps not in self._mp:
            self._mp[dps] = MPBackend(dps)
        return self._mp[dps]

    def iv(self, dps=30) -> IntervalBackend:
        if dps not in self._iv:
            self._iv[dps] = IntervalBackend(dps)
        return self._iv[dps]
