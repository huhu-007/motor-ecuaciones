"""ExpressionParser: tokenizador + descenso recursivo. Sin eval() ni ast.parse().

Gramática (la multiplicación implícita tiene la misma precedencia que '*'):
  relación := expr [ (=|<|>|<=|>=|!=) expr ]
  expr     := término { (+|-) término }
  término  := unario { (*|/|implícita) unario }
  unario   := (+|-) unario | potencia
  potencia := primario [ ^ unario ]              (asociativa por la derecha)
  primario := número | x | π | e | (expr) | f(arg) | log_b(arg) | √arg | ⁿ√arg | root_n(arg)
Un argumento SIN paréntesis debe ser un átomo aislado (sin x, √2); si sería ambiguo
(sin 2x, √x^2, √2x) se rechaza y se pide usar paréntesis.
"""
from __future__ import annotations
import re
from dataclasses import dataclass
from fractions import Fraction
from typing import List, NamedTuple, Optional

from .nodes import (Node, Num, Const, Var, Neg, BinOp, Root, Func, variables, sub)
from .functions import FUNCTION_NAMES

MAX_LEN = 500
MAX_DEPTH = 80
MAX_INDEX = 100

_SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"
_SUPMAP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻", "0123456789-")
_SUB = "₀₁₂₃₄₅₆₇₈₉"
_SUBMAP = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
_REPL = [("−", "-"), ("–", "-"), ("—", "-"), ("×", "*"), ("·", "*"), ("⋅", "*"),
         ("∙", "*"), ("÷", "/"), ("∕", "/"), ("**", "^"), ("≤", "<="), ("≥", ">="),
         ("≠", "!="), ("=<", "<="), ("=>", ">="), ("[", "("), ("]", ")"),
         ("{", "("), ("}", ")")]
_NAMES = sorted(list(FUNCTION_NAMES) + ["sqrt", "cbrt", "root", "pi"], key=lambda s: -len(s))


class ParseError(ValueError):
    pass


class Token(NamedTuple):
    kind: str      # NUM VAR CONST FUNC ROOT ROOTN OP LP RP REL UND
    value: object
    pos: int


def tokenize(src: str) -> List[Token]:
    if len(src) > MAX_LEN:
        raise ParseError(f"La ecuación es demasiado larga (máximo {MAX_LEN} caracteres).")
    s = src.strip().lower()
    for a, b in _REPL:
        s = s.replace(a, b)
    if re.search(r"\d\s*e\s*[+-]?\d", s.replace(" ", "")) and re.search(r"\d[e]\d", s.replace(" ", "")):
        raise ParseError("No se admite notación científica (3e5): escribe 3*10^5.")
    toks: List[Token] = []
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c.isspace():
            i += 1
            continue
        m = re.compile(r"\d+([.,]\d+)?|[.,]\d+").match(s, i)
        if m and (c.isdigit() or c in ".,"):
            toks.append(Token("NUM", Fraction(m.group().replace(",", ".")), i))
            i = m.end()
            continue
        if c in _SUP:
            j = i
            while j < n and s[j] in _SUP:
                j += 1
            digits = s[i:j].translate(_SUPMAP)
            if j < n and s[j] in "√":                       # ⁿ√ : índice de raíz
                if not digits.isdigit():
                    raise ParseError("Índice de raíz no válido.")
                toks.append(Token("ROOT", int(digits), i))
                i = j + 1
            else:                                            # x² , x⁻¹
                neg = digits.startswith("-")
                d = digits.lstrip("-")
                if not d.isdigit():
                    raise ParseError("Exponente no válido.")
                toks += [Token("OP", "^", i), Token("LP", "(", i)]
                if neg:
                    toks.append(Token("OP", "-", i))
                toks += [Token("NUM", Fraction(int(d)), i), Token("RP", ")", i)]
                i = j
            continue
        if c in _SUB:
            j = i
            while j < n and s[j] in _SUB:
                j += 1
            toks += [Token("UND", "_", i), Token("NUM", Fraction(int(s[i:j].translate(_SUBMAP))), i)]
            i = j
            continue
        if c in "√∛∜":
            toks.append(Token("ROOT", {"√": 2, "∛": 3, "∜": 4}[c], i))
            i += 1
            continue
        if c == "π":
            toks.append(Token("CONST", "pi", i))
            i += 1
            continue
        if c.isalpha():
            for name in _NAMES:
                if s.startswith(name, i):
                    i += len(name)
                    if name == "pi":
                        toks.append(Token("CONST", "pi", i))
                    elif name == "sqrt":
                        toks.append(Token("ROOT", 2, i))
                    elif name == "cbrt":
                        toks.append(Token("ROOT", 3, i))
                    elif name == "root":
                        toks.append(Token("ROOTN", None, i))
                    else:
                        toks.append(Token("FUNC", name, i))
                    break
            else:
                toks.append(Token("CONST" if c == "e" else "VAR", "e" if c == "e" else c, i))
                i += 1
            continue
        if c == "_":
            toks.append(Token("UND", "_", i)); i += 1; continue
        if c in "+-*/^":
            toks.append(Token("OP", c, i)); i += 1; continue
        if c == "(":
            toks.append(Token("LP", c, i)); i += 1; continue
        if c == ")":
            toks.append(Token("RP", c, i)); i += 1; continue
        if s.startswith(("<=", ">=", "!="), i):
            toks.append(Token("REL", s[i:i + 2], i)); i += 2; continue
        if c in "=<>":
            toks.append(Token("REL", c, i)); i += 1; continue
        if c == ",":
            raise ParseError("Coma no válida (la coma solo se admite como separador decimal, 1,5).")
        if c == "|":
            raise ParseError("Usa abs(...) en lugar de barras |...|.")
        raise ParseError(f"Carácter no admitido: '{c}'.")
    return toks


_START = {"NUM", "VAR", "CONST", "FUNC", "ROOT", "ROOTN", "LP"}


class ExpressionParser:
    def parse_relation(self, text: str):
        """Devuelve (izquierda, operador, derecha) ; operador None si no hay relación."""
        self.toks = tokenize(text)
        self.i = 0
        self.depth = 0
        if not self.toks:
            raise ParseError("La ecuación está vacía.")
        left = self.expr()
        op, right = None, None
        if self._kind() == "REL":
            op = self.toks[self.i].value
            self.i += 1
            right = self.expr()
        if self.i < len(self.toks):
            t = self.toks[self.i]
            if t.kind == "REL":
                raise ParseError("Solo puede haber un signo de relación (=, <, >, <=, >=).")
            raise ParseError(f"Símbolo inesperado: '{t.value}'.")
        return left, op, right

    def parse(self, text: str) -> Node:
        left, op, _ = self.parse_relation(text)
        if op is not None:
            raise ParseError("Se esperaba una expresión, no una relación.")
        return left

    # ---- utilidades ----
    def _kind(self, k=0):
        j = self.i + k
        return self.toks[j].kind if j < len(self.toks) else None

    def _val(self):
        return self.toks[self.i].value

    def _expect(self, kind, what):
        if self._kind() != kind:
            raise ParseError(f"Se esperaba {what}.")
        t = self.toks[self.i]
        self.i += 1
        return t

    def _enter(self):
        self.depth += 1
        if self.depth > MAX_DEPTH:
            raise ParseError("Expresión demasiado anidada.")

    def _leave(self):
        self.depth -= 1

    # ---- gramática ----
    def expr(self) -> Node:
        self._enter()
        node = self.term()
        while self._kind() == "OP" and self._val() in "+-":
            op = self._val(); self.i += 1
            node = BinOp(op, node, self.term())
        self._leave()
        return node

    def term(self) -> Node:
        node = self.unary()
        while True:
            k = self._kind()
            if k == "OP" and self._val() in "*/":
                op = self._val(); self.i += 1
                node = BinOp(op, node, self.unary())
            elif k in _START:
                if k == "NUM" and self.toks[self.i - 1].kind == "NUM":
                    raise ParseError("Dos números seguidos: falta un operador.")
                node = BinOp("*", node, self.unary())
            else:
                return node

    def unary(self) -> Node:
        if self._kind() == "OP" and self._val() in "+-":
            op = self._val(); self.i += 1
            self._enter()
            inner = self.unary()
            self._leave()
            return Neg(inner) if op == "-" else inner
        return self.power()

    def power(self) -> Node:
        base = self.primary()
        if self._kind() == "OP" and self._val() == "^":
            self.i += 1
            self._enter()
            exp = self.unary()
            self._leave()
            return BinOp("^", base, exp)
        return base

    def _arg(self, what) -> Node:
        """Argumento de función/raíz: (expr) o átomo aislado."""
        if self._kind() == "LP":
            self.i += 1
            e = self.expr()
            self._expect("RP", "')'")
            return e
        k = self._kind()
        if k == "NUM":
            a = Num(self._val()); self.i += 1
        elif k == "VAR":
            a = Var(self._val()); self.i += 1
        elif k == "CONST":
            a = Const(self._val()); self.i += 1
        elif k in ("FUNC", "ROOT", "ROOTN"):
            a = self.primary()
        else:
            raise ParseError(f"Falta el argumento de {what}.")
        nxt = self._kind()
        if nxt in ("NUM", "VAR", "CONST", "LP") or (nxt == "OP" and self._val() == "^"):
            raise ParseError(f"Ambiguo: usa paréntesis en el argumento de {what} (p. ej. sin(2x)).")
        return a

    def _index(self) -> int:
        self._expect("UND", "'_' seguido del índice")
        if self._kind() == "LP":
            self.i += 1
            t = self._expect("NUM", "un entero")
            self._expect("RP", "')'")
        else:
            t = self._expect("NUM", "un entero")
        v = t.value
        if v.denominator != 1 or not (2 <= v <= MAX_INDEX):
            raise ParseError(f"El índice de la raíz debe ser un entero entre 2 y {MAX_INDEX}.")
        return int(v)

    def primary(self) -> Node:
        k = self._kind()
        if k is None:
            raise ParseError("La expresión termina de forma inesperada.")
        t = self.toks[self.i]
        if k == "NUM":
            self.i += 1
            return Num(t.value)
        if k == "VAR":
            self.i += 1
            return Var(t.value)
        if k == "CONST":
            self.i += 1
            return Const(t.value)
        if k == "LP":
            self.i += 1
            e = self.expr()
            self._expect("RP", "')'")
            return e
        if k == "ROOT":
            self.i += 1
            if not (2 <= t.value <= MAX_INDEX):
                raise ParseError(f"El índice de la raíz debe estar entre 2 y {MAX_INDEX}.")
            return Root(t.value, self._arg("la raíz"))
        if k == "ROOTN":
            self.i += 1
            idx = self._index()
            return Root(idx, self._arg("root_n"))
        if k == "FUNC":
            self.i += 1
            name = t.value
            base = None
            if name == "log" and self._kind() == "UND":
                self.i += 1
                if self._kind() == "LP":
                    self.i += 1
                    base = self.expr()
                    self._expect("RP", "')'")
                elif self._kind() == "NUM":
                    base = Num(self._val()); self.i += 1
                elif self._kind() == "CONST":
                    base = Const(self._val()); self.i += 1
                else:
                    raise ParseError("Base de logaritmo no válida: log_2(x), log_(x+1)(y).")
            power = None
            if self._kind() == "OP" and self._val() == "^" and self._kind(1) == "NUM":
                power = self.toks[self.i + 1].value          # sin^2(x) = (sin x)^2
                if power.denominator != 1:
                    raise ParseError("Exponente no válido en sin^n(x).")
                self.i += 2
                if self._kind() != "LP":
                    raise ParseError("Tras sin^2 se requiere un argumento entre paréntesis.")
            f = Func(name, self._arg(name), base)
            return BinOp("^", f, Num(power)) if power is not None else f
        raise ParseError(f"Símbolo inesperado: '{t.value}'.")


@dataclass
class ParsedEquation:
    source: str
    relation: str              # '=', '<', '>', '<=', '>=', '!='
    left: Node
    right: Node
    function: Node             # left - right
    variable: str
    implicit_zero: bool


def parse_equation(text: str) -> ParsedEquation:
    p = ExpressionParser()
    left, op, right = p.parse_relation(text)
    implicit = op is None
    if implicit:
        op, right = "=", Num(Fraction(0))
    vs = variables(left) | variables(right)
    if len(vs) > 1:
        raise ParseError("Solo se admite una incógnita; encontradas: "
                         + ", ".join(sorted(vs))
                         + ". (¿función desconocida? Soportadas: sin, cos, tan, exp, abs, ln, log, sqrt, root_n)")
    var = next(iter(vs)) if vs else "x"
    return ParsedEquation(text, op, left, right, sub(left, right), var, implicit)
