"""API HTTP del motor de ecuaciones (FastAPI).

Arranque local:   uvicorn api:app --reload
Render (start):   uvicorn api:app --host 0.0.0.0 --port $PORT

Variables de entorno (todas opcionales):
  ALLOWED_ORIGINS   orígenes CORS separados por comas (por defecto "*")
  SOLVE_TIMEOUT     segundos máximos por /solve   (por defecto 25)
  ANALYZE_TIMEOUT   segundos máximos por /analyze (por defecto 15)
  RATE_LIMIT        peticiones por minuto y por IP (por defecto 20)
  MAX_CONCURRENT    cálculos simultáneos           (por defecto 1)
"""
from __future__ import annotations

import os
import time
from collections import defaultdict, deque
from threading import Lock
from typing import Deque, Dict, List

from fastapi import FastAPI, Request
from fastapi.exceptions import RequestValidationError
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse
from pydantic import BaseModel, Field, field_validator

from runner import Busy, Runner

# ----------------------------------------------------------------- configuración
MAX_INPUT_CHARS = 200
SOLVE_TIMEOUT = float(os.getenv("SOLVE_TIMEOUT", "25"))
ANALYZE_TIMEOUT = float(os.getenv("ANALYZE_TIMEOUT", "15"))
RATE_LIMIT = int(os.getenv("RATE_LIMIT", "20"))
MAX_CONCURRENT = int(os.getenv("MAX_CONCURRENT", "1"))
ORIGINS = [o.strip() for o in os.getenv("ALLOWED_ORIGINS", "*").split(",") if o.strip()]

runner = Runner(max_concurrent=MAX_CONCURRENT, queue_wait=5.0)

app = FastAPI(title="Motor de ecuaciones", version="1.0.0",
              description="Resuelve ecuaciones con Bolzano y exactitud verificada.")

app.add_middleware(
    CORSMiddleware,
    allow_origins=ORIGINS,
    allow_credentials=False,   # con "*" no se pueden permitir credenciales
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["*"],
)


# ----------------------------------------------------------------- modelos de entrada
class SolveRequest(BaseModel):
    equation: str = Field(..., min_length=1, max_length=MAX_INPUT_CHARS,
                          description="Ecuación, p. ej. 'x^2 - 2 = 0'")
    digits: int = Field(30, ge=1, le=50, description="Decimales mostrados por raíz")
    detail: bool = Field(False, description="true = incluye 'recognition' e 'interval' de cada raíz")

    @field_validator("equation")
    @classmethod
    def _strip(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("La ecuación está vacía.")
        return v


class AnalyzeRequest(BaseModel):
    expression: str = Field(..., min_length=1, max_length=MAX_INPUT_CHARS,
                            description="Expresión constante, p. ej. '√3 + √2'")
    digits: int = Field(60, ge=10, le=100, description="Dígitos usados en el análisis")

    @field_validator("expression")
    @classmethod
    def _strip(cls, v: str) -> str:
        v = v.strip()
        if not v:
            raise ValueError("La expresión está vacía.")
        return v


# ----------------------------------------------------------------- limitación de tasa
_hits: Dict[str, Deque[float]] = defaultdict(deque)
_hits_lock = Lock()


def _client_ip(request: Request) -> str:
    fwd = request.headers.get("x-forwarded-for")
    if fwd:
        return fwd.split(",")[0].strip()
    return request.client.host if request.client else "desconocida"


def _rate_limited(ip: str) -> bool:
    now = time.monotonic()
    with _hits_lock:
        q = _hits[ip]
        while q and now - q[0] > 60:
            q.popleft()
        if len(q) >= RATE_LIMIT:
            return True
        q.append(now)
        # limpieza ocasional para que el diccionario no crezca sin límite
        if len(_hits) > 5000:
            for k in [k for k, v in _hits.items() if not v or now - v[-1] > 60]:
                _hits.pop(k, None)
        return False


def _error(status_code: int, status: str, message: str, headers: Dict[str, str] | None = None):
    return JSONResponse(status_code=status_code, headers=headers,
                        content={"success": False, "status": status, "messages": [message]})


# ----------------------------------------------------------------- manejo de errores
@app.exception_handler(RequestValidationError)
async def _validation_error(_: Request, exc: RequestValidationError):
    msgs: List[str] = []
    for e in exc.errors():
        campo = ".".join(str(p) for p in e.get("loc", []) if p != "body")
        msgs.append(f"{campo}: {e.get('msg', 'valor no válido')}" if campo else e.get("msg", "valor no válido"))
    return JSONResponse(status_code=422,
                        content={"success": False, "status": "invalid_input", "messages": msgs})


# ----------------------------------------------------------------- endpoints
@app.get("/health")
def health():
    """Para despertar el servidor (plan gratuito) y comprobar que responde."""
    return {"ok": True}


def _compact(result: dict) -> dict:
    for r in result.get("roots") or []:
        for k in ("recognition", "interval", "_sortkey"):
            r.pop(k, None)
    return result


@app.post("/solve")
def solve(req: SolveRequest, request: Request):
    if _rate_limited(_client_ip(request)):
        return _error(429, "rate_limited", "Demasiadas peticiones. Espera un minuto.", {"Retry-After": "60"})
    try:
        result = runner.run("solve", req.equation, req.digits, SOLVE_TIMEOUT)
    except Busy:
        return _error(503, "busy", "El servidor está ocupado con otro cálculo. Inténtalo de nuevo en unos segundos.",
                      {"Retry-After": "5"})
    return result if req.detail else _compact(result)


@app.post("/analyze")
def analyze(req: AnalyzeRequest, request: Request):
    if _rate_limited(_client_ip(request)):
        return _error(429, "rate_limited", "Demasiadas peticiones. Espera un minuto.", {"Retry-After": "60"})
    try:
        return runner.run("analyze", req.expression, req.digits, ANALYZE_TIMEOUT)
    except Busy:
        return _error(503, "busy", "El servidor está ocupado con otro cálculo. Inténtalo de nuevo en unos segundos.",
                      {"Retry-After": "5"})
