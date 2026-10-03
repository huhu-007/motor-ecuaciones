"""Interfaz de consola (flujo original en español). La lógica matemática vive en engine.py."""
import argparse
import threading
from fractions import Fraction
from pathlib import Path

from .engine import EquationEngine
from .parser import ParseError
from .bisection import BisectionDigits, texto_periodo, fraccion_periodica, PrecisionExhausted, UndefinedInside
from .formatter import ResultFormatter

ARCHIVO = Path.cwd() / "datos.txt"      # siempre el mismo archivo (se actualiza, no se crean más)
DELAY = 0.001                              # segundos entre decimales (el original ponía 0.001 aunque decía 1 s)

INSTRUCCIONES = """
==============================================================
        CALCULADORA DE RAÍCES POR BISECCIÓN (Bolzano)
==============================================================
QUÉ HACE:
  1. Te pide una ecuación y busca intervalos donde la función cambia de signo
     (teorema de Bolzano) entre -100 y 100, comprobando la CONTINUIDAD con
     aritmética de intervalos (los polos y puntos no definidos no cuentan).
  2. Reduce el intervalo a la mitad una y otra vez y muestra un decimal nuevo
     cada pocos segundos, sin límite de decimales.
  3. Si detecta periodicidad la muestra, da la fracción y se detiene.
  4. Con ENTER se para el cálculo y se guarda todo en datos.txt.
  5. Después eliges:
        R -> analiza el número (racional, radicales de grado 2-100, sumas de
             radicales, radicales anidados, π, e, ln, polinomio mínimo...)
        T -> termina

CÓMO ESCRIBIR LA ECUACIÓN:
  x^3+3x^2-8   2x^2-3   x^3 = 2x+5   e^x-2   ln(x)-1   sin(x)   cos(x)-1
  √x-2   ∛x+2   ⁴√x   root_5(x)   log_2(x)   (x+1)/(x-2)   2π   x^x
  Se puede omitir el * (3x, 2(x+1), (x+1)(x-2)).  Sin '=' se entiende '= 0'.
  Argumentos sin paréntesis solo si son un átomo: sin x; en otro caso sin(2x).
==============================================================
"""


def descripcion(iv):
    if iv.lo == iv.hi:
        return f"raíz exacta en x = {float(iv.lo):g}"
    return f"entre {float(iv.lo):g} y {float(iv.hi):g}" + (f"  (multiplicidad {iv.multiplicity})" if iv.multiplicity and iv.multiplicity > 1 else "")


def elegir_intervalo(intervalos):
    if len(intervalos) == 1:
        return intervalos[0]
    print("\nHe encontrado varias raíces:")
    for i, iv in enumerate(intervalos, 1):
        print(f"  {i}) {descripcion(iv)}")
    while True:
        s = input("¿Cuál quieres calcular? (número): ").strip()
        if s.isdigit() and 1 <= int(s) <= len(intervalos):
            return intervalos[int(s) - 1]
        print("Número no válido.")


def calcular(iv, delay):
    """Un decimal cada `delay` s; ENTER para parar. Devuelve el BisectionDigits (o None)."""
    parar = threading.Event()

    def esperar_enter():
        input()
        parar.set()

    hilo = threading.Thread(target=esperar_enter, daemon=True)
    hilo.start()
    print("\nCalculando... pulsa ENTER para parar.\n")
    try:
        bd = BisectionDigits(iv.sign_fn, iv.lo, iv.hi) if iv.lo != iv.hi else BisectionDigits(iv.sign_fn, iv.lo, iv.lo)
        bd.next_digit()
        print(f"{bd.signo}{bd.entero}.", end="", flush=True)
        while not parar.is_set() and bd.periodo is None:
            d = bd.next_digit()
            print(d, end="", flush=True)
            parar.wait(delay)
    except (PrecisionExhausted, UndefinedInside) as e:
        print(f"\n✘ {e}")
        return None
    if not parar.is_set():
        print(f"\n\nPeriodicidad detectada: {texto_periodo(bd.signo, bd.entero, bd.periodo)}")
        print(f"Es un número racional (evidencia numérica): {fraccion_periodica(bd.signo, bd.entero, bd.periodo)}")
        print("\nPulsa ENTER para continuar...")
        hilo.join()
    return bd


def guardar(ecuacion, bd):
    valor = bd.value_text()
    lineas = [f"ecuacion={ecuacion}", f"valor={valor}", f"decimales={len(bd.dec)}"]
    if bd.periodo:
        lineas.append(f"periodo={texto_periodo(bd.signo, bd.entero, bd.periodo)}")
        lineas.append(f"fraccion={fraccion_periodica(bd.signo, bd.entero, bd.periodo)}")
    ARCHIVO.write_text("\n".join(lineas) + "\n", encoding="utf-8")


def analizar(engine):
    if not ARCHIVO.exists():
        print(f"No encuentro {ARCHIVO.name}.")
        return
    datos = {}
    for linea in ARCHIVO.read_text(encoding="utf-8").splitlines():
        if "=" in linea:
            k, v = linea.split("=", 1)
            datos[k.strip()] = v.strip()
    print(f"\n--- Análisis de {ARCHIVO.name} ---\nEcuación: {datos.get('ecuacion', '?')}")
    if "fraccion" in datos:
        print(f"Periodicidad detectada: RACIONAL = {datos['fraccion']} (evidencia numérica; "
              "el motor la verifica con aritmética exacta al resolver la ecuación).")
        return
    dec = int(datos.get("decimales", 0))
    if dec < 12:
        print(f"Solo hay {dec} decimales: son muy pocos. Déjalo calcular más (mejor 40 o más).")
        return
    print(f"x ≈ {datos['valor']}  ({dec} decimales)\n")
    rep = engine.analyze_number(datos["valor"], dec)
    print(rep["statement"])
    best = rep.get("best")
    if best:
        otras = [f for f in best["equivalent_forms"] if f != best["exact_form"]]
        print(f"  forma: {best['exact_form']}" + (f"   (≡ {', '.join(otras)})" if otras else ""))
    if rep.get("rational_powers"):
        print("  x^n racional para n =", ", ".join(map(str, list(rep["rational_powers"])[:10])))
    pol = rep.get("polynomial")
    if pol and pol.get("found"):
        print(f"  polinomio mínimo probable: {pol['polynomial']}\n  {pol.get('radical_remark', '')}")
    for n in rep["notes"]:
        print("  ·", n)
    print(f"\nPrecisión: {rep['digits_available']} decimales; se usan {rep['precision']['digits_used_to_identify']} "
          f"para identificar y {rep['precision']['holdout_digits']} para validar.")
    print("Ojo: es evidencia numérica, no una demostración. Más decimales = más fiabilidad.")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--delay", type=float, default=DELAY, help="segundos entre decimales")
    args = ap.parse_args(argv)
    engine = EquationEngine()
    print(INSTRUCCIONES)
    while True:
        ecuacion = input("Ecuación: ").strip()
        if not ecuacion:
            continue
        try:
            eq, ev, res = engine.prepare(ecuacion)
        except ParseError as e:
            print(f"✘ {e} Inténtalo de nuevo.\n")
            continue
        if eq.relation != "=":
            print("✘ Las inecuaciones se entienden pero todavía no se resuelven.\n")
            continue
        break
    print("\nBuscando intervalos con cambio de signo...")
    intervalos = res.intervals
    if not intervalos:
        print(ResultFormatter.to_text(engine.solve(ecuacion)))
        return
    iv = elegir_intervalo(intervalos)
    print(f"Raíz elegida: {descripcion(iv)}")
    bd = calcular(iv, args.delay)
    if bd is None or bd.entero is None:
        print("\nNo hay datos que guardar.")
        return
    guardar(ecuacion, bd)
    print(f"\nCálculo detenido. Datos guardados en {ARCHIVO.name} ({len(bd.dec)} decimales).")
    while True:
        op = input("\n[R] Analizar si es una raíz    [T] Terminar\n> ").strip().lower()
        if op == "r":
            analizar(engine)
            break
        if op == "t":
            print("Programa terminado.")
            break
        print("Pulsa R o T.")


if __name__ == "__main__":
    main()
