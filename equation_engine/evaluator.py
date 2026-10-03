"""FunctionEvaluator: un único recorrido del AST, parametrizado por backend."""
from __future__ import annotations
from fractions import Fraction
from math import gcd

from .nodes import Node, Num, Const, Var, Neg, BinOp, Root, Func, const_fraction, contains_var
from .functions import REGISTRY
from .errors import DomainError, NotExact

MAX_EXACT_EXP = 10 ** 6


class FunctionEvaluator:
    def __init__(self, node: Node, variable: str = "x"):
        self.node, self.var = node, variable

    def __call__(self, backend, x):
        return self.ev(self.node, backend, x)

    def evaluate(self, backend, x):
        return self.ev(self.node, backend, x)

    def ev(self, n, be, x):
        if isinstance(n, Num):
            return be.from_fraction(n.value)
        if isinstance(n, Const):
            return be.const(n.name)
        if isinstance(n, Var):
            return x
        if isinstance(n, Neg):
            return be.neg(self.ev(n.arg, be, x))
        if isinstance(n, BinOp):
            if n.op == "^":
                return self._pow(n, be, x)
            a, b = self.ev(n.left, be, x), self.ev(n.right, be, x)
            if n.op == "+":
                return be.add(a, b)
            if n.op == "-":
                return be.sub(a, b)
            if n.op == "*":
                return be.mul(a, b)
            if n.op == "/":
                be.require_nonzero(b, "división por cero")
                return be.div(a, b)
        if isinstance(n, Root):
            v = self.ev(n.radicand, be, x)
            if n.index % 2 == 0:
                be.require_nonnegative(v, f"raíz de índice par de un radicando negativo")
            return be.root(v, n.index)
        if isinstance(n, Func):
            v = self.ev(n.arg, be, x)
            if n.name == "log":
                be.require_positive(v, "log(x) exige x > 0")
                if n.base is None:
                    b = be.from_fraction(Fraction(10))
                else:
                    b = self.ev(n.base, be, x)
                    be.require_positive(b, "la base del logaritmo debe ser > 0")
                    be.require_nonzero(be.sub(b, be.from_fraction(Fraction(1))),
                                       "la base del logaritmo no puede ser 1")
                return be.log_base(v, b)
            return REGISTRY[n.name].evaluate(be, v)
        raise TypeError(n)

    def _pow(self, n, be, x):
        q = const_fraction(n.right)
        a = self.ev(n.left, be, x)
        if q is not None:
            if q.denominator == 1:
                p = int(q)
                if be.kind == "exact" and abs(p) > MAX_EXACT_EXP:
                    raise NotExact("exponente demasiado grande para cálculo exacto")
                if p <= 0:
                    be.require_nonzero(a, "0 elevado a un exponente <= 0 no está definido")
                return be.pow_int(a, p) if p != 0 else be.from_fraction(Fraction(1))
            r = q.denominator
            if r % 2 == 0:
                (be.require_positive if q < 0 else be.require_nonnegative)(
                    a, "base negativa con exponente fraccionario de denominador par")
            elif q < 0:
                be.require_nonzero(a, "0 elevado a exponente negativo")
            return be.pow_int(be.root(a, r), q.numerator)
        b = self.ev(n.right, be, x)
        be.require_positive(a, "la base debe ser > 0 si el exponente no es un entero/racional constante")
        return be.pow_general(a, b)
