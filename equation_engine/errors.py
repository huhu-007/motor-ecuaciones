class DomainError(Exception):
    """f no está definida en ese punto / la restricción de dominio se viola con certeza."""


class Uncertified(DomainError):
    """Con aritmética de intervalos no se puede garantizar el dominio (posible violación)."""


class NotExact(Exception):
    """El valor no es representable de forma exacta con el backend exacto."""
