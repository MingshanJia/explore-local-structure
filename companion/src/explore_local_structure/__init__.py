"""Local network structure without replacing NetworkX itself."""

from .closure import (
    average_closure,
    average_directed_closure,
    closure,
    directed_closure,
)
from .quadrangle import (
    average_inner_quadrangle_coefficient,
    average_outer_quadrangle_coefficient,
    global_quadrangle_coefficient,
    inner_quadrangle_coefficient,
    number_of_quadrangles,
    outer_quadrangle_coefficient,
    quadrangle_coefficients,
)

__version__ = "0.1.0"
__all__ = [
    "closure",
    "average_closure",
    "directed_closure",
    "average_directed_closure",
    "inner_quadrangle_coefficient",
    "outer_quadrangle_coefficient",
    "quadrangle_coefficients",
    "average_inner_quadrangle_coefficient",
    "average_outer_quadrangle_coefficient",
    "global_quadrangle_coefficient",
    "number_of_quadrangles",
]
