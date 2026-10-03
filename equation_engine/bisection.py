"""Bisección (Bolzano) y generación de decimales con detección de periodicidad.

BisectionDigits conserva el algoritmo original: los extremos son Fraction exactos, un decimal se
confirma cuando floor(a·10^d) == floor(b·10^d), y si la frontera decimal contenida en [a,b] es raíz
exacta el intervalo colapsa (los decimales siguientes son ceros -> periodo 0 -> racional)."""
from __future__ import annotations
from fractions import Fraction
from math import floor
from typing import Callable, List, Optional, Tuple

from .errors import DomainError, NotExact
from .exact import Rad

MIN_DIGITOS = 20     # mínimo de dígitos repetidos para dar una periodicidad por buena
MIN_REPS = 3         # y al menos 3 periodos completos


class PrecisionExhausted(Exception):
    pass


class UndefinedInside(Exception):
    pass


def make_sign_fn(evaluator, pool, exact_backend=None):
    """Signo RIGUROSO de f(x), x Fraction: intervalos con precisión creciente; si el encierro
    contiene 0, se intenta decidir con aritmética exacta (f(x) = 0 exacto)."""
    from .exact import ExactBackend
    eb = exact_backend or ExactBackend()

    def sign_at(x: Fraction) -> int:
        dps = max(30, int(x.denominator.bit_length() * 0.31) + 30)
        for attempt in range(6):
            be = pool.iv(dps)
            try:
                v = evaluator(be, be.from_fraction(x))
            except DomainError as e:
                raise UndefinedInside(str(e))
            s = be.sign(v)
            if s is not None:
                return s
            if attempt == 0:
                try:
                    z = evaluator(eb, Rad.rat(x))
                    if z.sign() == 0:
                        return 0
                except (NotExact, DomainError, ZeroDivisionError):
                    pass
            dps *= 2
        raise PrecisionExhausted(f"no se pudo decidir el signo de f en {float(x):g}")
    return sign_at


def buscar_periodo(dec):
    """Devuelve (parte_no_periodica, periodo) o None.  (función original, sin cambios)"""
    n = len(dec)
    for k in range(1, n // MIN_REPS + 1):
        run, j = 0, n - 1
        while j - k >= 0 and dec[j] == dec[j - k]:
            run += 1
            j -= 1
        if run >= max(MIN_DIGITOS, (MIN_REPS - 1) * k):
            inicio = n - run - k
            return dec[:inicio], dec[inicio:inicio + k]
    return None


def fraccion_periodica(signo, entero, periodo):
    previo, per = periodo
    p, k = len(previo), len(per)
    num_previo = int("".join(map(str, previo)) or 0)
    num_per = int("".join(map(str, per)) or 0)
    fr = Fraction(entero) + Fraction(num_previo, 10 ** p) + Fraction(num_per, 10 ** p * (10 ** k - 1))
    return -fr if signo else fr


def texto_periodo(signo, entero, periodo):
    previo, per = periodo
    return f"{signo}{entero}.{''.join(map(str, previo))}({''.join(map(str, per))})"


class BisectionDigits:
    """Iterador de decimales confirmados de la raíz de f en [lo, hi] (f(lo)·f(hi) < 0 o raíz exacta)."""

    def __init__(self, sign_at: Callable[[Fraction], int], lo: Fraction, hi: Fraction):
        lo, hi = Fraction(lo), Fraction(hi)
        while lo < 0 < hi:                       # normaliza: el intervalo no puede contener al 0
            s = sign_at(Fraction(0))
            if s == 0:
                lo = hi = Fraction(0)
                break
            sl = sign_at(lo)
            if s * sl < 0:
                hi = Fraction(0)
            else:
                lo = Fraction(0)
        if lo < 0 or (lo == hi and lo < 0):      # raíz negativa: se calcula la de f(-x)
            self.signo, self.g, self.a, self.b = "-", (lambda x: sign_at(-x)), -hi, -lo
        else:
            self.signo, self.g, self.a, self.b = "", sign_at, lo, hi
        self.fa = self.g(self.a) if self.a != self.b else 0
        self.d = 0
        self.entero: Optional[int] = None
        self.dec: List[int] = []
        self.periodo = None
        self.exact = self.a == self.b

    def _step(self):
        m = (self.a + self.b) / 2
        sm = self.g(m)
        if sm == 0:
            self.a = self.b = m
            self.exact = True
        elif self.fa * sm < 0:
            self.b = m
        else:
            self.a, self.fa = m, sm
        esc = 10 ** self.d
        if floor(self.a * esc) != floor(self.b * esc):     # ¿frontera decimal exacta es la raíz?
            c = Fraction(floor(self.a * esc) + 1, esc)
            if self.g(c) == 0:
                self.a = self.b = c
                self.exact = True

    def next_digit(self) -> Optional[int]:
        """Siguiente dígito confirmado (el primero es la parte entera)."""
        while floor(self.a * 10 ** self.d) != floor(self.b * 10 ** self.d):
            if self.a == self.b:
                break
            self._step()
        v = floor(self.a * 10 ** self.d)
        if self.d == 0:
            self.entero = v
            out = v
        else:
            out = v % 10
            self.dec.append(out)
            self.periodo = buscar_periodo(self.dec)
        self.d += 1
        return out

    def run(self, digits: int, stop_on_period=True):
        """Calcula hasta `digits` decimales o hasta detectar periodicidad."""
        if self.entero is None:
            self.next_digit()
        while len(self.dec) < digits and not (stop_on_period and self.periodo):
            self.next_digit()
        return self

    def bracket(self) -> Tuple[Fraction, Fraction]:
        """Intervalo encerrante en la variable original (deshace la reflexión)."""
        a, b = self.a, self.b
        return (-b, -a) if self.signo == "-" else (a, b)

    def value_text(self):
        return f"{self.signo}{self.entero}" + ("." + "".join(map(str, self.dec)) if self.dec else "")

    def rational_from_period(self):
        return fraccion_periodica(self.signo, self.entero, self.periodo) if self.periodo else None
