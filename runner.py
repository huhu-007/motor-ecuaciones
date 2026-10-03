"""Ejecución del motor con tiempo máximo duro (sin dependencias web).

Cada petición se ejecuta en un proceso aparte. Si excede el tiempo, el proceso se
termina de verdad (un hilo de Python no se puede matar, un proceso sí). Así una
ecuación pesada no puede bloquear el servidor.
"""
from __future__ import annotations

import json
import multiprocessing as mp
import threading
from typing import Any, Dict

# Configuración "ligera" para planes gratuitos (poca CPU).
LIGHT_CONFIG: Dict[str, Any] = {
    "analysis_budget": 3.0,     # segundos máximos por raíz en el reconocimiento exacto
    "max_analyzed_roots": 5,    # solo se analizan las 5 raíces más cercanas a 0
}


def _worker(conn, kind: str, text: str, digits: int, overrides: Dict[str, Any]) -> None:
    try:
        from equation_engine import EquationEngine

        if kind == "solve":
            eng = EquationEngine(digits=digits, **overrides)
            result = eng.solve(text)
        elif kind == "analyze":
            eng = EquationEngine(**overrides)
            result = eng.analyze_expression(text, digits)
        else:  # pragma: no cover
            raise ValueError(f"operación desconocida: {kind}")
        # to_json usa default=str: garantiza que todo es serializable
        conn.send(("ok", json.dumps(result, ensure_ascii=False, default=str)))
    except BaseException as exc:  # noqa: BLE001 - no queremos que el hijo muera sin avisar
        conn.send(("error", f"{type(exc).__name__}: {exc}"))
    finally:
        conn.close()


class Busy(Exception):
    """Todas las plazas de cálculo están ocupadas."""


class Runner:
    def __init__(self, max_concurrent: int = 1, queue_wait: float = 5.0):
        self._slots = threading.BoundedSemaphore(max_concurrent)
        self._queue_wait = queue_wait
        # "spawn" (proceso limpio) en vez de "fork": el servidor usa hilos y hacer fork
        # con hilos activos puede provocar bloqueos. Cuesta ~0,3 s por petición.
        self._ctx = mp.get_context("spawn")

    def run(self, kind: str, text: str, digits: int, timeout: float,
            overrides: Dict[str, Any] | None = None) -> Dict[str, Any]:
        """Devuelve el dict del motor, o un dict con status 'timeout' / 'internal_error'."""
        if not self._slots.acquire(timeout=self._queue_wait):
            raise Busy()
        try:
            overrides = dict(LIGHT_CONFIG, **(overrides or {}))
            parent, child = self._ctx.Pipe(duplex=False)
            proc = self._ctx.Process(target=_worker, args=(child, kind, text, digits, overrides),
                                     daemon=True)
            proc.start()
            child.close()
            try:
                if not parent.poll(timeout):
                    return {
                        "success": False,
                        "status": "timeout",
                        "messages": [f"El cálculo ha superado el tiempo máximo ({timeout:g} s). "
                                     "Prueba con una ecuación más sencilla."],
                    }
                tag, payload = parent.recv()
            except EOFError:
                return {"success": False, "status": "internal_error",
                        "messages": ["El proceso de cálculo terminó de forma inesperada."]}
            finally:
                if proc.is_alive():
                    proc.terminate()
                    proc.join(1)
                    if proc.is_alive():
                        proc.kill()
                proc.join(1)
                parent.close()
            if tag == "ok":
                return json.loads(payload)
            return {"success": False, "status": "internal_error",
                    "messages": ["Error interno del motor."], "error": payload}
        finally:
            self._slots.release()
