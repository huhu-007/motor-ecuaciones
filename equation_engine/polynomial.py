"""Polinomios exactos sobre Q (listas de Fraction, grado creciente), Yun y Sturm."""
from __future__ import annotations
from fractions import Fraction
from math import gcd
from typing import List, Optional, Tuple

from .nodes import Node, Num, Const, Var, Neg, BinOp, Root, Func, const_fraction

Poly = List[Fraction]


def trim(p: Poly) -> Poly:
    p = list(p)
    while p and p[-1] == 0:
        p.pop()
    return p


def deg(p: Poly) -> int:
    return len(trim(p)) - 1


def padd(a, b):
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def pneg(a):
    return [-c for c in a]


def psub(a, b):
    return padd(a, pneg(b))


def pmul(a, b):
    if not a or not b:
        return []
    r = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] += x * y
    return trim(r)


def pscale(a, q):
    return trim([c * q for c in a])


def pderiv(a):
    return trim([i * a[i] for i in range(1, len(a))])


def pdivmod(a: Poly, b: Poly) -> Tuple[Poly, Poly]:
    a, b = trim(a), trim(b)
    if not b:
        raise ZeroDivisionError
    q = [Fraction(0)] * max(len(a) - len(b) + 1, 0)
    r = list(a)
    while len(r) >= len(b) and r:
        c = r[-1] / b[-1]
        k = len(r) - len(b)
        q[k] = c
        for i, bc in enumerate(b):
            r[k + i] -= c * bc
        r = trim(r)
    return trim(q), trim(r)


def pgcd(a: Poly, b: Poly) -> Poly:
    a, b = trim(a), trim(b)
    while b:
        a, b = b, pdivmod(a, b)[1]
    return pscale(a, 1 / a[-1]) if a else []


def peval(p: Poly, x: Fraction) -> Fraction:
    r = Fraction(0)
    for c in reversed(p):
        r = r * x + c
    return r


def psign(p: Poly, x: Fraction) -> int:
    v = peval(p, x)
    return (v > 0) - (v < 0)


def to_int_coeffs(p: Poly) -> List[int]:
    p = trim(p)
    if not p:
        return []
    from math import lcm
    L = 1
    for c in p:
        L = lcm(L, c.denominator)
    ints = [int(c * L) for c in p]
    g = 0
    for c in ints:
        g = gcd(g, abs(c))
    ints = [c // g for c in ints]
    if ints[-1] < 0:
        ints = [-c for c in ints]
    return ints


def poly_text(p: Poly, var="x") -> str:
    """Antes poli_txt. Coeficientes enteros (se escala)."""
    ints = to_int_coeffs(p)
    n, s = len(ints) - 1, ""
    for i in range(n, -1, -1):
        c = ints[i]
        if c == 0:
            continue
        t = "" if (abs(c) == 1 and i > 0) else str(abs(c))
        t += "" if i == 0 else (var if i == 1 else f"{var}^{i}")
        s += (("-" if c < 0 else "") + t) if not s else ((" - " if c < 0 else " + ") + t)
    return s or "0"


# ---------------- Yun: descomposición libre de cuadrados ----------------
def yun(p: Poly) -> List[Tuple[Poly, int]]:
    """p = c * prod a_i^i con a_i libres de cuadrados y coprimos. Devuelve [(a_i, i)]."""
    p = trim(p)
    if deg(p) < 1:
        return []
    dp = pderiv(p)
    a0 = pgcd(p, dp)
    b, c = pdivmod(p, a0)[0], pdivmod(dp, a0)[0]
    d = psub(c, pderiv(b))
    out, i = [], 1
    while deg(b) > 0:
        a = pgcd(b, d)
        if deg(a) > 0:
            out.append((pscale(a, 1 / a[-1]), i))
        b, c = pdivmod(b, a)[0], pdivmod(d, a)[0]
        d = psub(c, pderiv(b))
        i += 1
    return out


# ---------------- Sturm ----------------
def sturm_chain(p: Poly) -> List[Poly]:
    chain = [trim(p), pderiv(p)]
    while chain[-1]:
        r = pdivmod(chain[-2], chain[-1])[1]
        if not r:
            break
        chain.append(pneg(r))
    return [c for c in chain if c]


def _var_at(chain, x: Optional[Fraction], side=0):
    """Variaciones de signo en x (o en ±∞ si x es None y side=±1)."""
    signs = []
    for c in chain:
        if x is None:
            v = c[-1] * (1 if side > 0 else (-1) ** (len(c) - 1))
        else:
            v = peval(c, x)
        if v != 0:
            signs.append(v > 0)
    return sum(1 for i in range(len(signs) - 1) if signs[i] != signs[i + 1])


def count_roots(chain, a: Fraction, b: Fraction) -> int:
    """Raíces reales distintas en (a, b] (polinomio libre de cuadrados)."""
    return _var_at(chain, a) - _var_at(chain, b)


def count_all_real(chain) -> int:
    return _var_at(chain, None, -1) - _var_at(chain, None, +1)


def cauchy_bound(p: Poly) -> Fraction:
    p = trim(p)
    return 1 + max((abs(c / p[-1]) for c in p[:-1]), default=Fraction(0))


def isolate(chain, a: Fraction, b: Fraction, max_depth=200):
    """Intervalos (lo, hi] con exactamente una raíz cada uno (Sturm + bisección)."""
    out = []

    def rec(lo, hi, depth):
        n = count_roots(chain, lo, hi)
        if n == 0:
            return
        if n == 1 or depth >= max_depth:
            out.append((lo, hi))
            return
        m = (lo + hi) / 2
        rec(lo, m, depth + 1)
        rec(m, hi, depth + 1)

    rec(a, b, 0)
    return out


# ---------------- función racional desde el AST ----------------
def ratfunc(n: Node, var: str) -> Optional[Tuple[Poly, Poly]]:
    """(numerador, denominador) si n es una función racional de var; si no, None."""
    if isinstance(n, Num):
        return trim([n.value]), [Fraction(1)]
    if isinstance(n, Var):
        return ([Fraction(0), Fraction(1)], [Fraction(1)]) if n.name == var else None
    if isinstance(n, Neg):
        r = ratfunc(n.arg, var)
        return None if r is None else (pneg(r[0]), r[1])
    if isinstance(n, BinOp):
        if n.op == "^":
            q = const_fraction(n.right)
            if q is None or q.denominator != 1 or abs(q) > 200:
                return None
            l = ratfunc(n.left, var)
            if l is None:
                return None
            k = int(q)
            num, den = [Fraction(1)], [Fraction(1)]
            for _ in range(abs(k)):
                num, den = pmul(num, l[0]), pmul(den, l[1])
            if k < 0:
                if not trim(l[0]):
                    return None
                num, den = den, num
            return num, den
        l, r = ratfunc(n.left, var), ratfunc(n.right, var)
        if l is None or r is None:
            return None
        if n.op == "+":
            return padd(pmul(l[0], r[1]), pmul(r[0], l[1])), pmul(l[1], r[1])
        if n.op == "-":
            return psub(pmul(l[0], r[1]), pmul(r[0], l[1])), pmul(l[1], r[1])
        if n.op == "*":
            return pmul(l[0], r[0]), pmul(l[1], r[1])
        if n.op == "/":
            if not trim(r[0]):
                return None
            return pmul(l[0], r[1]), pmul(l[1], r[0])
    return None


def reduced_ratfunc(n: Node, var: str):
    """(num, den, común) con num/den coprimos; 'común' = gcd cancelado (puntos no definidos)."""
    r = ratfunc(n, var)
    if r is None:
        return None
    num, den = trim(r[0]), trim(r[1])
    if not num:
        return [], [Fraction(1)], []
    g = pgcd(num, den)
    if deg(g) > 0:
        num, den = pdivmod(num, g)[0], pdivmod(den, g)[0]
    else:
        g = []
    return num, den, g
