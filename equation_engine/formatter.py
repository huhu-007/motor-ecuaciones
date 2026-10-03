"""ResultFormatter: convierte el dict del motor en texto legible (es) sin acoplar la lógica al print()."""
from __future__ import annotations
import json
from typing import Any, Dict, List

_STATUS = {
    "solved": "Resuelta",
    "no_real_solution": "La ecuación NO tiene solución real (demostrado)",
    "no_root_found_in_range": "No se encontró ninguna raíz en el rango de búsqueda",
    "bolzano_not_applicable": "No se pudo aplicar Bolzano",
    "undefined_in_range": "La función no está definida en el rango requerido",
    "identity": "Identidad (se cumple en todo el dominio)",
    "invalid_input": "Entrada no válida",
    "unsupported_relation": "Relación no soportada todavía",
}
_LEVEL = {
    "exact_proven": "EXACTA (demostrada)",
    "exact_verified": "EXACTA (f(candidato)=0 verificado; unicidad no demostrada)",
    "numerical_evidence": "evidencia numérica (NO demostrada)",
    "none": "sin forma exacta identificada",
}


class ResultFormatter:
    @staticmethod
    def to_json(result: Dict[str, Any], indent=2) -> str:
        return json.dumps(ResultFormatter.public(result), ensure_ascii=False, indent=indent, default=str)

    @staticmethod
    def public(result: Dict[str, Any]) -> Dict[str, Any]:
        """Quita claves internas (_sortkey...) para exponer el dict a la API."""
        def clean(o):
            if isinstance(o, dict):
                return {str(k): clean(v) for k, v in o.items() if not str(k).startswith("_")}
            if isinstance(o, list):
                return [clean(v) for v in o]
            return o
        return clean(result)

    @staticmethod
    def to_text(result: Dict[str, Any], max_roots: int = 20) -> str:
        r = result
        L: List[str] = [f"Ecuación: {r['equation']}", f"Estado: {_STATUS.get(r['status'], r['status'])}"]
        if r.get("domain"):
            L.append(f"Dominio: {r['domain']['text']}")
        if r.get("method"):
            L.append(f"Método: {r['method']}")
        L.append(f"Bolzano aplicable: {'sí' if r.get('bolzano_applicable') else 'no'}")
        for i, root in enumerate(r.get("roots", [])[:max_roots], 1):
            L.append(f"\nRaíz {i}: x ≈ {root['decimal_form']}  ({root['digits']} decimales)")
            if root.get("multiplicity") and root["multiplicity"] > 1:
                L.append(f"  multiplicidad: {root['multiplicity']}")
            L.append(f"  forma exacta: {root.get('exact_form') or root.get('candidate_form') or '—'}"
                     f"  → {_LEVEL.get(root.get('exactness', 'none'))}")
            mp = root.get("minimal_polynomial")
            if mp:
                L.append(f"  polinomio: {mp['polynomial']}  ({mp['status']})")
                if mp.get("radical_remark"):
                    L.append(f"  nota: {mp['radical_remark']}")
        if len(r.get("roots", [])) > max_roots:
            L.append(f"\n… y {len(r['roots']) - max_roots} raíces más.")
        for m in r.get("messages", []):
            L.append(f"\n{m}")
        for w in r.get("warnings", []):
            L.append(f"⚠ {w}")
        return "\n".join(L)
