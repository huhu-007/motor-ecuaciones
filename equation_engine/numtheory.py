"""Teoría de números mínima: primalidad, factorización (Pollard rho), utilidades."""
import math
import random
from typing import Dict

_SMALL = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37]


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    for p in _SMALL:
        if n % p == 0:
            return n == p
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in _SMALL:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True


def _rho(n: int, limit: int):
    if n % 2 == 0:
        return 2
    rnd = random.Random(n)
    for _ in range(8):
        c, x, y, d, it = rnd.randrange(1, n), rnd.randrange(2, n), 0, 1, 0
        y = x
        while d == 1 and it < limit:
            x = (x * x + c) % n
            y = (y * y + c) % n
            y = (y * y + c) % n
            d = math.gcd(abs(x - y), n)
            it += 1
        if 1 < d < n:
            return d
    return None


def factorize(n: int, rho_limit: int = 200000) -> Dict[int, int]:
    """Factorización de n>0. Si un cofactor no se consigue factorizar se devuelve tal cual
    (clave compuesta) y `factorize.complete` queda en False."""
    factorize.complete = True
    out: Dict[int, int] = {}
    if n < 1:
        raise ValueError("n debe ser positivo")
    for p in range(2, 10000):
        if p * p > n:
            break
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
    stack = [n] if n > 1 else []
    while stack:
        m = stack.pop()
        if m == 1:
            continue
        if is_prime(m):
            out[m] = out.get(m, 0) + 1
            continue
        d = _rho(m, rho_limit)
        if d is None:
            factorize.complete = False
            out[m] = out.get(m, 0) + 1
            continue
        stack += [d, m // d]
    return out


factorize.complete = True


def primes_upto(n: int):
    return [p for p in range(2, n + 1) if is_prime(p)]


def split_perfect_power(n: int, k: int = 2):
    """n = a^k * m con m sin factores k-ésimos. Devuelve (a, m).  (antes: separar_cuadrado)"""
    a, m = 1, 1
    for p, e in factorize(n).items():
        a *= p ** (e // k)
        m *= p ** (e % k)
    return a, m
