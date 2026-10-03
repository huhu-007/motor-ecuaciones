# Desplegar la API en Render (gratis)

Necesitas: cuenta de GitHub y cuenta de Render (render.com). Nada se instala en tu ordenador.

## 1. Subir el código a GitHub
1. En github.com: **New repository** → nombre `motor-ecuaciones` → *Private* o *Public* → **Create**.
2. En la página del repositorio vacío pulsa **uploading an existing file**.
3. Descomprime el zip en tu ordenador y arrastra TODO el contenido de la carpeta
   (`equation_engine/`, `tests/`, `api.py`, `runner.py`, `requirements.txt`, `.python-version`, `main.py`…).
   Deben quedar `api.py` y `requirements.txt` en la raíz del repositorio, no dentro de otra carpeta.
4. **Commit changes**.

## 2. Crear el servicio en Render
1. Render → **New +** → **Web Service** → conecta tu cuenta de GitHub → elige `motor-ecuaciones`.
2. Rellena:
   - **Language**: Python 3
   - **Build Command**: `pip install -r requirements.txt`
   - **Start Command**: `uvicorn api:app --host 0.0.0.0 --port $PORT`
   - **Instance Type**: Free (comprueba que sigue apareciendo; si no, ve al punto "Si no hay plan gratuito")
3. (Opcional) En *Environment Variables* añade `PYTHON_VERSION` = `3.12.3` si Render no lee `.python-version`.
4. **Create Web Service**. El primer despliegue tarda unos minutos.

## 3. Comprobar que funciona
Render te da una URL tipo `https://motor-ecuaciones-xxxx.onrender.com`. Ábrela en el navegador con `/health` al final:

    https://motor-ecuaciones-xxxx.onrender.com/health      →  {"ok": true}
    https://motor-ecuaciones-xxxx.onrender.com/docs        →  página interactiva para probar /solve

En `/docs` abre **POST /solve** → *Try it out* → `{"equation": "x^2 - 2 = 0"}` → *Execute*.

## 4. Limitar quién puede llamarla (cuando tengas la web de Lovable)
En Render → tu servicio → *Environment* → añade `ALLOWED_ORIGINS` con la URL de tu web de Lovable
(por ejemplo `https://mi-app.lovable.app`; varias, separadas por comas). Mientras no la pongas, acepta cualquier origen.

## Cosas a saber del plan gratuito
- Se apaga tras 15 min sin uso; la primera petición tarda 30–60 s en despertarlo. La web debe llamar a `/health`
  al cargar y mostrar «Despertando el servidor…».
- Tiene muy poca CPU: las ecuaciones pesadas (p. ej. `sin(x)` con muchas raíces) irán más lentas y pueden devolver
  `status: "timeout"` (25 s en /solve, 15 s en /analyze).
- Solo se hace un cálculo a la vez (`MAX_CONCURRENT=1`); si llega otro, espera hasta 5 s y, si no hay hueco, devuelve HTTP 503.
- Límite por defecto: 20 peticiones por minuto y por IP (`RATE_LIMIT`).

## Si no hay plan gratuito
Las condiciones de los planes gratuitos cambian. Alternativa: Google Cloud Run (nivel gratuito, pide tarjeta). Dímelo y preparo un `Dockerfile`.

## Prueba local (opcional, con Python instalado)
    pip install -r requirements.txt
    uvicorn api:app --reload
    # abrir http://127.0.0.1:8000/docs
