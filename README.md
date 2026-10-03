# Motor matemático de ecuaciones (Bolzano + exactitud verificada)

Solo motor: parser, validación, dominio, búsqueda de raíces, bisección, análisis simbólico y resultado
estructurado (dict/JSON) listo para el futuro frontend Laravel. No hay interfaz web.

    python main.py                      # consola interactiva (flujo original, en español)
    python -m unittest -v tests.test_engine

    from equation_engine import EquationEngine
    eng = EquationEngine(digits=30)
    eng.solve("x^2 - 2 = 0")            # dict JSON-serializable
    eng.analyze_expression("√3 + √2")   # expresiones constantes (exacto si es posible)
    eng.analyze_number("1.41421356…", 60)   # reconocimiento numérico (solo evidencia)

Dependencias: Python ≥ 3.9 y `mpmath`.

## Arquitectura (módulos y responsabilidad)
| Módulo | Responsabilidad |
|---|---|
| `parser.py`, `nodes.py` | **ExpressionParser**: tokenizador + descenso recursivo → AST. Sin `eval`/`ast.parse`. |
| `functions.py` | Registro de funciones (sin, cos, tan, exp, abs, ln, log). Añadir una función = una entrada. |
| `evaluator.py` | **FunctionEvaluator**: un solo recorrido del AST con 3 backends (abajo). |
| `backends.py` | `MPBackend` (alta precisión), `IntervalBackend` (encierros rigurosos), `BackendPool` (sin estado global). |
| `exact.py` | `ExactBackend`: álgebra exacta de radicales (Besicovitch), múltiplos de π y potencias de e. |
| `domain.py` | **DomainAnalyzer**: restricciones (`x ≠ 2`, `x > 0`, `x ≥ 3`, `x ≠ π/2 + kπ`…). |
| `bolzano.py` | **BolzanoAnalyzer**: f(a), f(b) definidas, continuidad CERTIFICADA (intervalos), f(a)·f(b) ≤ 0. |
| `rootfinder.py` | **RootFinder**: rejilla configurable, anti-polo, fronteras de dominio, raíces tangentes, ruta polinómica (Yun + Sturm). |
| `bisection.py` | Bisección con decimales confirmados y periodicidad (algoritmo original conservado). |
| `polynomial.py` | Polinomios exactos sobre Q, Yun, Sturm, funciones racionales. |
| `recognizers.py` | **ExactFormRecognizer / RadicalRecognizer / PolynomialAnalyzer / ConstantRecognizer**. |
| `engine.py` | Fachada `EquationEngine` (orquesta todo y verifica de forma exacta). |
| `formatter.py` | **ResultFormatter**: dict → JSON / texto. La lógica nunca usa `print()`. |
| `cli.py` | Interfaz de consola original (ENTER, `datos.txt`, [R]/[T]). |

## Niveles de certeza (nunca se mezclan)
* `exact_proven`   – f(candidato)=0 en aritmética exacta **y** la raíz es única en el intervalo (Sturm / punto exacto).
* `exact_verified` – f(candidato)=0 exacto, pero la unicidad no está demostrada (funciones no polinómicas).
* `numerical_evidence` – ajuste de alta precisión (con cota de coincidencia fortuita). **No es una demostración**; el campo `exact_form` queda a `null` y el candidato va en `candidate_form`.
* `none` – sin forma exacta identificada (se informa del polinomio mínimo si se halla).

Estados: `solved`, `no_real_solution` (solo si está DEMOSTRADO: constante falsa o Sturm sin raíces),
`no_root_found_in_range`, `bolzano_not_applicable` (con el mensaje
«No se puede aplicar el teorema de Bolzano en este tipo de ecuaciones.»), `undefined_in_range`,
`identity`, `unsupported_relation`, `invalid_input`. Son afirmaciones distintas y no se confunden.

## Operaciones soportadas
+ − * / ^, multiplicación implícita (`2x`, `3(x+1)`, `(x+1)(x-1)`, `2π`, `2sin(x)`), decimales (`1.5`, `1,5`),
fracciones anidadas, potencias con exponente arbitrario (`x^2`, `x^(x+1)`, `2^x`, `e^x`, `x^x`),
`√`, `∛`, `∜`, `ⁿ√` (superíndices) y `root_n(x)` con n = 2…100 (pares exigen radicando ≥ 0; impares aceptan negativos),
`log(x)` (base 10), `log_b(x)`, `ln(x)`, `sin`, `cos`, `tan`, `exp`, `abs`, `sin^2(x)`, constantes simbólicas π y e,
`=` (sin `=` se entiende `= 0`). Inecuaciones: se interpretan (no se resuelven).
Reconocimiento exacto: racionales, `a·ⁿ√b/c`, `√a/√b`, `a/√b`, productos y cocientes (n = 2…100),
sumas/diferencias (`√a+√b`, `a√b+c√d`, `a+∛b+…`), radicales anidados `ⁿ√(a±b√c)`, q·π, q·e, ln(q), e^q,
polinomio mínimo (grados hasta donde alcance la precisión; configurable ≤ 100).

## No soportado (a propósito o todavía)
Integrales, derivadas, límites, sumatorios; más de una incógnita; inecuaciones (solo parseo); números complejos;
funciones inversas/hiperbólicas (el registro lo permite); notación científica `3e5`; factoriales; barras `|x|` (usa `abs`).

## Limitaciones matemáticas que permanecen
* Búsqueda de raíces solo en el rango configurado (por defecto [−100, 100]); fuera del rango solo se **cuentan** para polinomios (Sturm).
* Dos raíces de f no polinómica dentro de una misma celda sin cambio de signo pueden pasarse por alto (se avisa si hay mínimos |f|≈0).
* Raíces tangentes no polinómicas solo se aceptan si se verifican exactamente (si no, `unverified_candidates`).
* Unicidad de la raíz en el intervalo solo se demuestra para polinomios/funciones racionales.
* El reconocimiento de sumas de radicales usa familias acotadas (hasta 3 primos, índices mixtos pequeños); una forma fuera de ellas no se detecta (falso negativo, nunca falso positivo con la protección de altura 10^(D/4) y dígitos de reserva).
* Irreducibilidad del polinomio mínimo solo demostrada hasta grado 3; sin fórmulas de Cardano/Ferrari.
* Evaluación exacta: no hay `ln(2)` ni `sin(1)` simbólicos; raíz de una suma de radicales (√(1+√2) evaluado exactamente) no implementada.
* Convención real: `x^(p/q)` con q impar usa la raíz real (`(-8)^(1/3) = -2`); base negativa con exponente variable/irracional está indefinida; `0^0` indefinido; `log(x)` = log₁₀.
* Rendimiento: el reconocimiento puede consumir `analysis_budget` s por raíz (8 s por defecto); con muchas raíces (p. ej. `sin(x)`) solo se analizan las `max_analyzed_roots` más cercanas a 0.
