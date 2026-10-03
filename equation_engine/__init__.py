"""Motor matemático de ecuaciones (Bolzano + exactitud verificada). Sin interfaz web."""
from .engine import EquationEngine, EngineConfig
from .parser import parse_equation, ParseError
__all__ = ["EquationEngine", "EngineConfig", "parse_equation", "ParseError"]
