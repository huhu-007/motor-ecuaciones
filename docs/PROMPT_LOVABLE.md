# Prompt para pegar en Lovable

(Sustituye `https://TU-SERVICIO.onrender.com` por la URL real de Render.)

---

Crea una web en español, de una sola página, llamada «Resolutor de ecuaciones». Es solo el frontend: toda la
matemática la hace una API externa que ya existe. NO reimplementes ningún cálculo en el navegador ni en el backend.

**API** (base: `https://TU-SERVICIO.onrender.com`, guárdala en una constante configurable):
- `GET /health` → `{"ok": true}`
- `POST /solve` con `{"equation": string, "digits": number}` → resultado de la ecuación
- `POST /analyze` con `{"expression": string, "digits": number}` → análisis de una expresión constante (p. ej. `√3+√2`)

**Arranque del servidor.** Está en un plan gratuito que se duerme: al cargar la página llama a `/health` y, mientras
no responda, muestra «Despertando el servidor, puede tardar hasta un minuto…». Haz la primera petición real con un
timeout de 90 s; las siguientes, de 40 s.

**Interfaz**
1. Campo de texto grande para la ecuación (máximo 200 caracteres), con ejemplos que se pueden pulsar:
   `x^2 - 2 = 0`, `x^3 - 6x^2 + 11x - 6 = 0`, `cos(x) = x`, `ln(x) = 0`.
2. Teclado matemático con botones que insertan en el cursor: `√(` `∛(` `ⁿ√` (inserta `root_n(`) `^` `π` `e` `sin(` `cos(` `tan(`
   `ln(` `log(` `exp(` `abs(` `( )` `=`.
3. Selector de modo: «Ecuación» (usa /solve) y «Expresión constante» (usa /analyze).
4. Botón «Resolver» (deshabilitado mientras carga, con indicador de progreso).
5. Selector de decimales (10–50, por defecto 30).

**Resultado de /solve.** Una tarjeta por cada elemento de `roots`, ordenadas como llegan, con:
- `decimal_form` en fuente monoespaciada.
- `exact_form` renderizada con KaTeX o con una fuente que soporte √, ⁿ√, π y superíndices; si es `null`, no mostrar forma exacta.
- Una insignia según `exactness`:
  - `exact_proven` → verde, «Exacta (demostrada)»
  - `exact_verified` → verde claro, «Exacta (verificada)»
  - `numerical_evidence` → ámbar, «Evidencia numérica, no demostrada»; mostrar `candidate_form` como «Candidato», NUNCA como solución exacta
  - `none` → gris, «Solo valor decimal»
- Un desplegable «Detalles» con `bolzano.reason`, el intervalo `bolzano.a` / `bolzano.b`, `multiplicity` y
  `minimal_polynomial.polynomial` si existen.
Encima de las tarjetas muestra `domain.text` como «Dominio: …», y todos los `warnings` en un aviso amarillo.
Si hay más de 12 raíces, muestra las 12 primeras y un botón «Ver todas».

**Según `status`** (muestra siempre los textos de `messages` tal cual vienen):
- `solved`: tarjetas de raíces.
- `no_real_solution`: mensaje claro «No tiene solución real» con los `messages`.
- `no_root_found_in_range`: «No se encontró ninguna raíz en [-100, 100]» (aclarando que NO significa que no exista).
- `bolzano_not_applicable`: mostrar exactamente «No se puede aplicar el teorema de Bolzano en este tipo de ecuaciones.» y debajo el resto de `messages`.
- `undefined_in_range`, `identity`, `unsupported_relation`: mostrar `messages` en un cuadro informativo.
- `invalid_input`: cuadro de error con `messages` (error de sintaxis); no borrar lo que el usuario escribió.
- `timeout`: «El cálculo ha tardado demasiado. Prueba con una ecuación más sencilla.»

**Resultado de /analyze.** Muestra `exact_form`, `decimal_form` y la misma insignia de `exactness`; si hay
`candidate_form` sin forma exacta, preséntalo como candidato. Muestra `warnings`.

**Errores HTTP.** El cuerpo siempre trae `messages`: 422 → mostrarlos; 429 → «Demasiadas peticiones, espera un minuto»;
503 → «El servidor está ocupado, reintentando…» y reintenta una vez a los 5 s; fallo de red → «No se pudo conectar con el servidor».

**Estilo.** Limpio y sobrio, tema claro/oscuro, adaptable a móvil, todo en español. Guarda un historial de las
últimas 10 consultas en el navegador (localStorage) con opción de borrarlo.
