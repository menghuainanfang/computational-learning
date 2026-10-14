"""Optimize f: R^n -> R with vector inputs and structured results.

Example: gradient_descent(lambda w: (w[0] - 3)**2, [0.0]).x
Exports the solver, gradient utilities, and result containers.
"""
from .differentiation import make_gradient, finite_difference_gradient
from .optimizers import gradient_descent
from .results import OptimizationResult, IterationHistory

__all__ = ["gradient_descent", "make_gradient", "finite_difference_gradient",
           "OptimizationResult", "IterationHistory"]
