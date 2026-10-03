"""Registro de funciones. Añadir una función = añadir una entrada aquí.

Cada FunctionSpec define:
  evaluate(backend, v)   -> valor (mp / intervalo / exacto); puede exigir dominio
  constraints(arg_node)  -> restricciones [(nodo, relación)] con relación en {'>', '>=', '!='}
"""
from dataclasses import dataclass
from typing import Callable, Dict, List, Tuple
from .nodes import Func


@dataclass(frozen=True)
class FunctionSpec:
    name: str
    evaluate: Callable
    constraints: Callable
    doc: str = ""


def _tan(b, v):
    c = b.cos(v)
    b.require_nonzero(c, "tan(x) no está definida donde cos(x)=0")
    return b.div(b.sin(v), c)


def _ln(b, v):
    b.require_positive(v, "ln(x) exige x > 0")
    return b.log(v)


REGISTRY: Dict[str, FunctionSpec] = {
    "sin": FunctionSpec("sin", lambda b, v: b.sin(v), lambda a: [], "seno"),
    "cos": FunctionSpec("cos", lambda b, v: b.cos(v), lambda a: [], "coseno"),
    "tan": FunctionSpec("tan", _tan, lambda a: [(Func("cos", a), "!=")], "tangente"),
    "exp": FunctionSpec("exp", lambda b, v: b.exp(v), lambda a: [], "e^x"),
    "abs": FunctionSpec("abs", lambda b, v: b.abs(v), lambda a: [], "valor absoluto"),
    "ln": FunctionSpec("ln", _ln, lambda a: [(a, ">")], "logaritmo natural"),
    # 'log' (base 10 por defecto, o log_b) se trata en el evaluador por tener base opcional
    "log": FunctionSpec("log", None, lambda a: [(a, ">")], "logaritmo (base 10 por defecto)"),
}

FUNCTION_NAMES = tuple(REGISTRY)
