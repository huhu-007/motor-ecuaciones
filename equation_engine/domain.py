"""DomainAnalyzer: restricciones de dominio reales de una expresión."""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from typing import List, Optional, Tuple

from .nodes import (Node, Num, Const, Var, Neg, BinOp, Root, Func, walk, to_text,
                    const_fraction, contains_var)
from .functions import REGISTRY

_SYM = {">": ">", ">=": "≥", "!=": "≠"}


@dataclass
class Constraint:
    expr: Node
    relation: str            # '>', '>=', '!='  (expr REL 0)
    reason: str
    description: str = ""

    def to_dict(self):
        return {"expression": to_text(self.expr), "relation": self.relation,
                "reason": self.reason, "description": self.description}


def linear_coeffs(n: Node, var: str) -> Optional[Tuple[Fraction, Fraction]]:
    """(a, b) si n == a*var + b de forma estructural y exacta; si no, None."""
    c = const_fraction(n)
    if c is not None:
        return Fraction(0), c
    if isinstance(n, Var):
        return (Fraction(1), Fraction(0)) if n.name == var else None
    if isinstance(n, Neg):
        r = linear_coeffs(n.arg, var)
        return None if r is None else (-r[0], -r[1])
    if isinstance(n, BinOp):
        l, r = linear_coeffs(n.left, var), linear_coeffs(n.right, var)
        if n.op in "+-" and l and r:
            s = 1 if n.op == "+" else -1
            return l[0] + s * r[0], l[1] + s * r[1]
        if n.op == "*" and l and r:
            if l[0] == 0:
                return l[1] * r[0], l[1] * r[1]
            if r[0] == 0:
                return r[1] * l[0], r[1] * l[1]
        if n.op == "/" and l and r and r[0] == 0 and r[1] != 0:
            return l[0] / r[1], l[1] / r[1]
    return None


def _q(v: Fraction) -> str:
    return str(v.numerator) if v.denominator == 1 else f"{v.numerator}/{v.denominator}"


class DomainAnalyzer:
    def __init__(self, variable: str = "x"):
        self.var = variable

    def analyze(self, node: Node) -> List[Constraint]:
        out: List[Constraint] = []
        seen = set()

        def add(expr, rel, reason):
            if not contains_var(expr):
                return    # constante (2, e, π...): no restringe a la incógnita; si fuese falsa, el evaluador lo detecta
            key = (to_text(expr), rel)
            if key in seen:
                return
            seen.add(key)
            out.append(Constraint(expr, rel, reason, self._describe(expr, rel)))

        for n in walk(node):
            if isinstance(n, BinOp) and n.op == "/":
                add(n.right, "!=", "denominador ≠ 0")
            elif isinstance(n, BinOp) and n.op == "^":
                q = const_fraction(n.right)
                if q is not None:
                    if q.denominator == 1:
                        if q <= 0:
                            add(n.left, "!=", "base ≠ 0 (exponente ≤ 0)")
                    elif q.denominator % 2 == 0:
                        add(n.left, ">" if q < 0 else ">=", "exponente con denominador par")
                    elif q < 0:
                        add(n.left, "!=", "base ≠ 0 (exponente negativo)")
                else:
                    add(n.left, ">", "base > 0 con exponente variable/irracional")
            elif isinstance(n, Root):
                if n.index % 2 == 0:
                    add(n.radicand, ">=", f"raíz de índice par {n.index}")
            elif isinstance(n, Func):
                if n.name == "log":
                    add(n.arg, ">", "argumento del logaritmo > 0")
                    if n.base is not None:
                        add(n.base, ">", "base del logaritmo > 0")
                        add(BinOp("-", n.base, Num(Fraction(1))), "!=", "base del logaritmo ≠ 1")
                else:
                    for e, rel in REGISTRY[n.name].constraints(n.arg):
                        add(e, rel, f"dominio de {n.name}")
        return out

    def _describe(self, expr: Node, rel: str) -> str:
        if isinstance(expr, Func) and expr.name == "cos" and rel == "!=":
            lin = linear_coeffs(expr.arg, self.var)
            if lin and lin[0] != 0:
                a, b = lin
                if a == 1 and b == 0:
                    return f"{self.var} ≠ π/2 + kπ (k entero)"
            return f"cos({to_text(expr.arg)}) ≠ 0"
        if (rel == ">" and isinstance(expr, BinOp) and expr.op == "^"):
            q = const_fraction(expr.right)
            if q is not None and q.denominator == 1 and q > 0 and q % 2 == 0:
                return self._describe(expr.left, "!=")      # u^(2k) > 0  <=>  u ≠ 0
        lin = linear_coeffs(expr, self.var)
        if lin and lin[0] != 0:
            a, b = lin
            r = -b / a
            sym = _SYM[rel]
            if a < 0 and rel != "!=":
                sym = {">": "<", "≥": "≤"}[sym]
            return f"{self.var} {sym} {_q(r)}"
        return f"{to_text(expr)} {_SYM[rel]} 0"

    @staticmethod
    def text(constraints: List[Constraint]) -> str:
        return "; ".join(c.description for c in constraints) if constraints else "ℝ (sin restricciones)"
