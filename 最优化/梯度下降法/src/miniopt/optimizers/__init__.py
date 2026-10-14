"""Unconstrained solvers mapping scalar objectives and vectors to results.

Gradient descent follows the negative gradient.
Example: gradient_descent(lambda w: w[0]**2, [1.0]).
"""
from .gradient_descent import gradient_descent

__all__ = ["gradient_descent"]
