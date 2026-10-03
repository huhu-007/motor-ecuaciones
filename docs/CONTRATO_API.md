# Contrato de la API

Base: `https://TU-SERVICIO.onrender.com`. Todo es JSON. Los textos de mensajes (`messages`, `warnings`) vienen en español.

## GET /health
`200 {"ok": true}`. Úsalo al cargar la web para despertar el servidor.

## POST /solve
Petición:

    { "equation": "x^2 - 2 = 0", "digits": 30, "detail": false }

| campo | tipo | reglas |
|---|---|---|
| `equation` | string | 1–200 caracteres. Sin `=` se entiende `= 0` |
| `digits` | entero | 1–50, por defecto 30 (decimales por raíz) |
| `detail` | booleano | `true` añade `recognition` e `interval` en cada raíz (respuesta mucho más grande) |

Respuesta (campos principales):

| campo | significado |
|---|---|
| `success` | `true` si el motor llegó a una conclusión (incluye «no hay solución real» demostrada) |
| `status` | ver tabla de estados |
| `messages` | lista de textos para mostrar al usuario |
| `warnings` | avisos (p. ej. «Solo se busca en [-100, 100]…») |
| `equation`, `normalized_function`, `variable`, `relation` | lo que se interpretó |
| `domain.text` | dominio legible (p. ej. `x > 0`, `x ≠ π/2 + kπ (k entero)`); `domain.constraints` lo detalla |
| `search_range` | `[-100, 100]` |
| `bolzano_applicable` / `bolzano_message` | si Bolzano aplica y, si no, el mensaje |
| `roots` | lista de raíces (abajo) |

### Estados (`status`)
| valor | qué mostrar |
|---|---|
| `solved` | lista de raíces |
| `no_real_solution` | «sin solución real» (está **demostrado**) — usar `messages` |
| `no_root_found_in_range` | no se halló raíz en el rango (NO es una demostración) |
| `bolzano_not_applicable` | «No se puede aplicar el teorema de Bolzano en este tipo de ecuaciones.» + explicación en `messages` |
| `undefined_in_range` | la función no está definida en el rango |
| `identity` | la ecuación se cumple siempre en su dominio |
| `unsupported_relation` | inecuación: se interpreta pero no se resuelve |
| `invalid_input` | error de sintaxis; el texto está en `messages` |
| `timeout` | el cálculo tardó demasiado |
| `internal_error` | fallo interno |

### Cada raíz (`roots[i]`)
| campo | significado |
|---|---|
| `index` | número de raíz |
| `decimal_form` | valor decimal (texto, con `digits` decimales) |
| `exact_form` | forma exacta con símbolos (`-√2`, `π/2`…) o `null` |
| `exact_form_text` | misma forma en ASCII (`-sqrt(2)`) o `null` |
| `exactness` | nivel de certeza (abajo) |
| `candidate_form` | candidato de forma exacta (solo evidencia si `exactness` es `numerical_evidence`) |
| `multiplicity` | multiplicidad (o `null`) |
| `bolzano` | datos del intervalo: `a`, `b`, `sign_a`, `sign_b`, `applicable`, `reason`… |
| `proof` | detalle de la verificación exacta (o `null`) |
| `minimal_polynomial` | polinomio mínimo cuando se halla (`polynomial`, `degree`, `status`) |
| `analysis_skipped` | `true` si hay muchas raíces y esta no se analizó a fondo |

### Niveles de certeza (`exactness`) — mostrarlos siempre, nunca mezclarlos
| valor | significado | cómo presentarlo |
|---|---|---|
| `exact_proven` | f(raíz)=0 comprobado en aritmética exacta y raíz única en su intervalo | «Exacta (demostrada)» |
| `exact_verified` | f(raíz)=0 exacto, pero la unicidad no está demostrada | «Exacta (verificada)» |
| `numerical_evidence` | ajuste numérico de alta precisión; **no es demostración** | «Evidencia numérica» (usar `candidate_form`, nunca `exact_form`) |
| `none` | sin forma exacta identificada | solo el decimal |

## POST /analyze
Para expresiones **constantes** (sin variable), p. ej. `√3 + √2`.

    { "expression": "√3 + √2", "digits": 60 }

`expression`: 1–200 caracteres. `digits`: 10–100. Respuesta: `success`, `status` (`exact`, `numerical_evidence`,
`numerical_only`, `invalid_input`, `undefined`, `timeout`…), `exact_form`, `exact_form_text`, `exactness`,
`decimal_form`, `candidate_form`, `warnings`.

## Errores HTTP
Siempre con el mismo formato `{"success": false, "status": "...", "messages": ["..."]}`:

| HTTP | `status` | significado |
|---|---|---|
| 422 | `invalid_input` | petición mal formada o fuera de límites (p. ej. más de 200 caracteres) |
| 429 | `rate_limited` | demasiadas peticiones (cabecera `Retry-After`) |
| 503 | `busy` | el servidor está haciendo otro cálculo; reintentar en unos segundos |

Los errores de sintaxis de la ecuación NO son errores HTTP: llegan con 200 y `status: "invalid_input"`.
