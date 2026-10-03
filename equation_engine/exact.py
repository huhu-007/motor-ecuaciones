"""Aritmética EXACTA.

Rad:       elementos de Q[p^(1/N), ...] con p primo. Representación canónica por
           monomios  coef * prod p^(e_p),  0 < e_p < 1.  Por el teorema de Besicovitch
           (los p^(j/N) con primos distintos son linealmente independientes sobre Q)
           la representación es única: «vector nulo  <=>  número nulo».
ConstMono: coef * pi^a * e^q  (un solo monomio; suma solo entre monomios iguales).
ExactBackend: backend del evaluador para demostraciones exactas (f(candidato) = 0).
Lanza NotExact cuando el resultado no es representable: jamás se aproxima.
"""
from __future__ import annotations
from fractions import Fraction
from math import floor, gcd
from typing import Dict, Tuple

from .numtheory import factorize
from .errors import DomainError, NotExact

Mono = Tuple[Tuple[int, Fraction], ...]     # ((primo, exponente en (0,1)), ...) ordenado


def _mono_mul(a: Mono, b: Mono):
    """Devuelve (monomio, factor_racional)."""
    d: Dict[int, Fraction] = dict(a)
    for p, e in b:
        d[p] = d.get(p, Fraction(0)) + e
    factor = Fraction(1)
    out = []
    for p, e in sorted(d.items()):
        k = floor(e)
        if k:
            factor *= Fraction(p) ** k
            e -= k
        if e:
            out.append((p, e))
    return tuple(out), factor


class Rad:
    __slots__ = ("t",)

    def __init__(self, terms: Dict[Mono, Fraction] | None = None):
        self.t = {m: c for m, c in (terms or {}).items() if c != 0}

    # --- construcción ---
    @staticmethod
    def rat(q) -> "Rad":
        return Rad({(): Fraction(q)})

    @staticmethod
    def _factor_fraction_root(q: Fraction, n: int) -> "Rad":
        """q^(1/n) para q racional positivo."""
        if q <= 0:
            raise DomainError("radicando no positivo")
        coef, mono = Fraction(1), {}
        for num, sgn in ((q.numerator, 1), (q.denominator, -1)):
            for p, e in factorize(num).items():
                if not factorize.complete and p > 10 ** 12 and not _is_prime_safe(p):
                    raise NotExact("no se pudo factorizar el radicando")
                k, r = divmod(e, n)
                coef *= Fraction(p) ** (sgn * k)
                if r:
                    mono[p] = mono.get(p, Fraction(0)) + sgn * Fraction(r, n)
        m, f = (), Fraction(1)
        for p, e in mono.items():
            m, g = _mono_mul(m, ((p, e % 1),)) if e % 1 else (m, Fraction(1))
            f *= g * Fraction(p) ** floor(e)
        return Rad({m: coef * f})

    # --- consultas ---
    def is_zero(self):
        return not self.t

    def is_rational(self):
        return all(m == () for m in self.t)

    def rational(self) -> Fraction:
        return self.t.get((), Fraction(0))

    def is_monomial(self):
        return len(self.t) == 1

    # --- aritmética ---
    def __add__(self, o):
        d = dict(self.t)
        for m, c in o.t.items():
            d[m] = d.get(m, Fraction(0)) + c
        return Rad(d)

    def __neg__(self):
        return Rad({m: -c for m, c in self.t.items()})

    def __sub__(self, o):
        return self + (-o)

    def __mul__(self, o):
        d: Dict[Mono, Fraction] = {}
        for m1, c1 in self.t.items():
            for m2, c2 in o.t.items():
                m, f = _mono_mul(m1, m2)
                d[m] = d.get(m, Fraction(0)) + c1 * c2 * f
        return Rad(d)

    def scale(self, q):
        return Rad({m: c * q for m, c in self.t.items()})

    def pow_int(self, n: int) -> "Rad":
        if n < 0:
            return self.inverse().pow_int(-n)
        if self.is_rational():                       # evita desbordar la memoria con 7^(10^6)
            q = self.rational()
            if max(abs(q.numerator).bit_length(), q.denominator.bit_length()) * n > 2_000_000:
                raise NotExact("potencia demasiado grande para cálculo exacto")
        elif n > 2000:
            raise NotExact("potencia demasiado grande para cálculo exacto")
        r, b = Rad.rat(1), self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
            if len(r.t) > 2000 or len(b.t) > 2000:
                raise NotExact("expresión demasiado grande")
        return r

    def inverse(self) -> "Rad":
        if self.is_zero():
            raise ZeroDivisionError
        if self.is_monomial():
            (m, c), = self.t.items()
            inv, f = (), Fraction(1)
            for p, e in m:                       # p^-e = p^-1 * p^(1-e)
                inv, g = _mono_mul(inv, ((p, 1 - e),))
                f *= g / p
            return Rad({inv: f / c})
        # caso general: resolver alpha*y = 1 en la base de monomios
        primes = sorted({p for m in self.t for p, _ in m})
        N = 1
        for m in self.t:
            for _, e in m:
                N = N * e.denominator // gcd(N, e.denominator)
        size = N ** len(primes)
        if size > 400:
            raise NotExact("inversión demasiado grande")
        import itertools
        basis = []
        for js in itertools.product(range(N), repeat=len(primes)):
            basis.append(tuple((p, Fraction(j, N)) for p, j in zip(primes, js) if j))
        index = {m: i for i, m in enumerate(basis)}
        cols = []
        for m in basis:
            prod = self * Rad({m: Fraction(1)})
            col = [Fraction(0)] * size
            for mm, c in prod.t.items():
                col[index[mm]] = c
            cols.append(col)
        A = [[cols[j][i] for j in range(size)] + [Fraction(1 if i == index[()] else 0)]
             for i in range(size)]
        for c in range(size):                      # Gauss-Jordan exacto
            piv = next((r for r in range(c, size) if A[r][c] != 0), None)
            if piv is None:
                raise ZeroDivisionError
            A[c], A[piv] = A[piv], A[c]
            pv = A[c][c]
            A[c] = [v / pv for v in A[c]]
            for r in range(size):
                if r != c and A[r][c] != 0:
                    f = A[r][c]
                    A[r] = [a - f * b for a, b in zip(A[r], A[c])]
        return Rad({basis[i]: A[i][size] for i in range(size)})

    def root(self, n: int) -> "Rad":
        """Raíz n-ésima real; solo si el radicando es un monomio."""
        if self.is_zero():
            return self
        if not self.is_monomial():
            raise NotExact("raíz de una suma de radicales (no implementada de forma exacta)")
        (m, c), = self.t.items()
        if c < 0:
            if n % 2 == 0:
                raise DomainError("raíz par de un número negativo")
            return -((-self).root(n))
        res = Rad._factor_fraction_root(c, n)
        for p, e in m:
            res = res * Rad({((p, (e / n) % 1),) if (e / n) % 1 else (): Fraction(1)})
        return res

    # --- signo / valor numérico rigurosos ---
    def to_iv(self, iv):
        tot = iv.mpf(0)
        for m, c in self.t.items():
            term = iv.mpf(c.numerator) / c.denominator
            for p, e in m:
                term = term * iv.exp(iv.mpf(e.numerator) / e.denominator * iv.log(iv.mpf(p)))
            tot = tot + term
        return tot

    def sign(self) -> int:
        if self.is_zero():
            return 0
        if self.is_rational():
            return 1 if self.rational() > 0 else -1
        from mpmath.ctx_iv import MPIntervalContext
        iv = MPIntervalContext()
        for dps in (30, 60, 120, 240, 480, 960):
            iv.dps = dps
            v = self.to_iv(iv)
            if v.a > 0:
                return 1
            if v.b < 0:
                return -1
        raise NotExact("signo indeterminado")

    def to_mpf(self, ctx):
        tot = ctx.mpf(0)
        for m, c in self.t.items():
            term = ctx.mpf(c.numerator) / c.denominator
            for p, e in m:
                term *= ctx.exp(ctx.mpf(e.numerator) / e.denominator * ctx.log(p))
            tot += term
        return tot

    def __eq__(self, o):
        return isinstance(o, Rad) and self.t == o.t

    def __hash__(self):
        return hash(frozenset(self.t.items()))


def _is_prime_safe(p):
    from .numtheory import is_prime
    return is_prime(p)


class ConstMono:
    __slots__ = ("c", "pi", "e")

    def __init__(self, c, pi=0, e=0):
        self.c, self.pi, self.e = Fraction(c), int(pi), Fraction(e)

    def sign(self):
        return (self.c > 0) - (self.c < 0)

    def to_mpf(self, ctx):
        return ctx.mpf(self.c.numerator) / self.c.denominator * ctx.pi ** self.pi * ctx.exp(
            ctx.mpf(self.e.numerator) / self.e.denominator)

    def to_iv(self, iv):
        return (iv.mpf(self.c.numerator) / self.c.denominator * iv.pi ** self.pi
                * iv.exp(iv.mpf(self.e.numerator) / self.e.denominator))


def _cm(v) -> ConstMono:
    if isinstance(v, ConstMono):
        return v
    if isinstance(v, Rad) and v.is_rational():
        return ConstMono(v.rational())
    raise NotExact("mezcla de radicales con π/e")


def _simplify(c: ConstMono):
    if c.pi == 0 and c.e == 0:
        return Rad.rat(c.c)
    if c.c == 0:
        return Rad.rat(0)
    return c


class ExactBackend:
    """Backend exacto para FunctionEvaluator. Valores: Rad | ConstMono."""
    kind = "exact"

    def from_fraction(self, q):
        return Rad.rat(q)

    def const(self, name):
        return ConstMono(1, 1, 0) if name == "pi" else ConstMono(1, 0, 1)

    def neg(self, a):
        return -a if isinstance(a, Rad) else ConstMono(-a.c, a.pi, a.e)

    def _both_rad(self, a, b):
        return isinstance(a, Rad) and isinstance(b, Rad)

    def add(self, a, b):
        if self._both_rad(a, b):
            return a + b
        a, b = _cm(a), _cm(b)
        if a.c == 0:
            return b
        if b.c == 0:
            return a
        if (a.pi, a.e) != (b.pi, b.e):
            raise NotExact("suma de constantes distintas")
        return _simplify(ConstMono(a.c + b.c, a.pi, a.e))

    def sub(self, a, b):
        return self.add(a, self.neg(b))

    def mul(self, a, b):
        if self._both_rad(a, b):
            return a * b
        a, b = _cm(a), _cm(b)
        return _simplify(ConstMono(a.c * b.c, a.pi + b.pi, a.e + b.e))

    def div(self, a, b):
        if self._both_rad(a, b):
            return a * b.inverse()
        a, b = _cm(a), _cm(b)
        return _simplify(ConstMono(a.c / b.c, a.pi - b.pi, a.e - b.e))

    def pow_int(self, a, n):
        if isinstance(a, Rad):
            return a.pow_int(n)
        if abs(n) > 2000:
            raise NotExact("potencia demasiado grande para cálculo exacto")
        return _simplify(ConstMono(a.c ** n, a.pi * n, a.e * n))

    def root(self, a, n):
        if isinstance(a, Rad):
            return a.root(n)
        if a.pi == 0 and a.c == 1 and (a.e / n).denominator <= 10 ** 6:
            return ConstMono(1, 0, a.e / n)
        raise NotExact("raíz de constante trascendente")

    def pow_general(self, a, b):
        if isinstance(b, Rad) and b.is_rational():
            q = b.rational()
            if isinstance(a, ConstMono) and a.c == 1 and a.pi == 0:
                return _simplify(ConstMono(1, 0, a.e * q))
            if isinstance(a, Rad):
                return self.pow_int(self.root(a, q.denominator), q.numerator)
        raise NotExact("potencia no representable")

    def exp(self, a):
        if isinstance(a, Rad) and a.is_rational():
            return _simplify(ConstMono(1, 0, a.rational()))
        raise NotExact("exp no representable")

    def log(self, a):
        if isinstance(a, Rad) and a.is_rational() and a.rational() == 1:
            return Rad.rat(0)
        if isinstance(a, ConstMono) and a.c == 1 and a.pi == 0:
            return Rad.rat(a.e)
        raise NotExact("log no representable")

    def log_base(self, v, b):
        if isinstance(v, Rad) and isinstance(b, Rad) and v.is_rational() and b.is_rational():
            bv, vv = b.rational(), v.rational()
            for k in range(-64, 65):
                if bv ** k == vv:
                    return Rad.rat(k)
        return self.div(self.log(v), self.log(b))

    def _angle(self, a):
        if isinstance(a, Rad) and a.is_zero():
            return Fraction(0)
        if isinstance(a, ConstMono) and a.pi == 1 and a.e == 0:
            k = a.c * 2
            if k.denominator == 1:
                return a.c
        raise NotExact("ángulo no notable")

    def sin(self, a):
        k = int(self._angle(a) * 2)
        return Rad.rat([0, 1, 0, -1][k % 4])

    def cos(self, a):
        k = int(self._angle(a) * 2)
        return Rad.rat([1, 0, -1, 0][k % 4])

    def abs(self, a):
        return self.neg(a) if self._sign(a) < 0 else a

    def _sign(self, a):
        return a.sign()

    def _chk(self, a):
        s = self._sign(a)
        return s

    def require_nonzero(self, v, msg="división por cero"):
        if self._sign(v) == 0:
            raise DomainError(msg)

    def require_positive(self, v, msg="se requiere valor > 0"):
        if self._sign(v) <= 0:
            raise DomainError(msg)

    def require_nonnegative(self, v, msg="se requiere valor >= 0"):
        if self._sign(v) < 0:
            raise DomainError(msg)

    def is_zero(self, v):
        return self._sign(v) == 0

    def to_mpf(self, v, ctx):
        return v.to_mpf(ctx)
