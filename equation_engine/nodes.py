"""Nodos del AST y utilidades (recorrido, texto, constantes racionales)."""
from __future__ import annotations
from dataclasses import dataclass
from fractions import Fraction
from typing import Optional, Set


class Node:
    __slots__ = ()


@dataclass(frozen=True)
class Num(Node):
    value: Fraction


@dataclass(frozen=True)
class Const(Node):
    name: str            # 'pi' | 'e'  (constantes simbólicas, NO variables)


@dataclass(frozen=True)
class Var(Node):
    name: str


@dataclass(frozen=True)
class Neg(Node):
    arg: Node


@dataclass(frozen=True)
class BinOp(Node):
    op: str              # + - * / ^
    left: Node
    right: Node


@dataclass(frozen=True)
class Root(Node):
    index: int           # 2..100
    radicand: Node


@dataclass(frozen=True)
class Func(Node):
    name: str
    arg: Node
    base: Optional[Node] = None   # solo para log_b


def children(n: Node):
    if isinstance(n, Neg):
        return (n.arg,)
    if isinstance(n, BinOp):
        return (n.left, n.right)
    if isinstance(n, Root):
        return (n.radicand,)
    if isinstance(n, Func):
        return (n.arg,) + ((n.base,) if n.base is not None else ())
    return ()


def walk(n: Node):
    yield n
    for c in children(n):
        yield from walk(c)


def variables(n: Node) -> Set[str]:
    return {m.name for m in walk(n) if isinstance(m, Var)}


def contains_var(n: Node) -> bool:
    return any(isinstance(m, Var) for m in walk(n))


def node_count(n: Node) -> int:
    return sum(1 for _ in walk(n))


def const_fraction(n: Node) -> Optional[Fraction]:
    """Valor racional exacto si el subárbol es una constante racional pura."""
    if isinstance(n, Num):
        return n.value
    if isinstance(n, Neg):
        v = const_fraction(n.arg)
        return None if v is None else -v
    if isinstance(n, BinOp):
        a, b = const_fraction(n.left), const_fraction(n.right)
        if a is None or b is None:
            return None
        if n.op == "+":
            return a + b
        if n.op == "-":
            return a - b
        if n.op == "*":
            return a * b
        if n.op == "/":
            return None if b == 0 else a / b
        if n.op == "^":
            if b.denominator != 1 or abs(b) > 200 or (a == 0 and b <= 0):
                return None
            return a ** int(b)
    return None


def sub(a: Node, b: Node) -> Node:
    if isinstance(b, Num) and b.value == 0:
        return a
    return BinOp("-", a, b)


# ---------------- texto (re-analizable por el propio parser) ----------------
def _wrap(t, cond):
    return f"({t})" if cond else t


def _fmt(n: Node):
    if isinstance(n, Num):
        v = n.value
        if v.denominator == 1:
            return str(v.numerator), (5 if v >= 0 else 3)
        return f"{v.numerator}/{v.denominator}", (2 if v > 0 else 3)
    if isinstance(n, Const):
        return n.name, 5
    if isinstance(n, Var):
        return n.name, 5
    if isinstance(n, Neg):
        t, p = _fmt(n.arg)
        return "-" + _wrap(t, p <= 3), 3
    if isinstance(n, BinOp):
        lt, lp = _fmt(n.left)
        rt, rp = _fmt(n.right)
        neg_r = rt.startswith("-")
        if n.op == "+":
            return f"{_wrap(lt, lp < 1)}+{_wrap(rt, rp < 1 or neg_r)}", 1
        if n.op == "-":
            return f"{_wrap(lt, lp < 1)}-{_wrap(rt, rp < 2 or neg_r)}", 1
        if n.op == "*":
            return f"{_wrap(lt, lp < 2)}*{_wrap(rt, rp <= 2 or neg_r)}", 2
        if n.op == "/":
            return f"{_wrap(lt, lp < 2)}/{_wrap(rt, rp < 4 or neg_r)}", 2
        if n.op == "^":
            return f"{_wrap(lt, lp <= 4 or lt.startswith('-'))}^{_wrap(rt, rp < 4 or neg_r)}", 4
    if isinstance(n, Root):
        t, _ = _fmt(n.radicand)
        return (f"sqrt({t})" if n.index == 2 else f"root_{n.index}({t})"), 5
    if isinstance(n, Func):
        t, _ = _fmt(n.arg)
        if n.name == "log" and n.base is not None:
            bt, bp = _fmt(n.base)
            b = bt if (isinstance(n.base, Num) and n.base.value.denominator == 1
                       and n.base.value > 0) else f"({bt})"
            return f"log_{b}({t})", 5
        return f"{n.name}({t})", 5
    raise TypeError(n)


def to_text(n: Node) -> str:
    return _fmt(n)[0]
