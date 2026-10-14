"""Differentiate scalar objectives f: R^n -> R.

make_gradient(f) returns a callable mapping vectors to matching gradients.
Central differences approximate partial derivatives for verification.
Example: make_gradient(lambda w: w[0]**2)([2.0]) returns [4.0].
"""
import numpy as np
from autograd import grad


def validate_point(w):
    """Convert a nonempty finite real vector to an independent float array."""
    if np.iscomplexobj(w):
        raise ValueError("The parameter vector must be real-valued")
    point = np.asarray(w, dtype=float).copy()
    if point.ndim != 1 or point.size == 0:
        raise ValueError("The parameter vector must be a nonempty 1-D array")
    if not np.all(np.isfinite(point)):
        raise ValueError("The parameter vector must contain finite values")
    return point


def objective_value(f, w):
    """Evaluate a real scalar objective; reject vector-valued outputs."""
    value = np.asarray(f(w))
    if value.shape != () or np.iscomplexobj(value):
        raise ValueError("The objective must return a real scalar")
    return float(value)


def gradient_value(gradient, w):
    """Evaluate a real gradient and check its dimension."""
    value = np.asarray(gradient(w))
    if value.shape != w.shape or np.iscomplexobj(value):
        raise ValueError("The gradient must be real and match the parameter shape")
    return np.asarray(value, dtype=float)


def make_gradient(f):
    """Build an Autograd gradient. Use autograd.numpy inside the objective.

    Failed tracing raises ValueError with the underlying error as its cause.
    Nonfinite derivatives are returned for the solver to report numerically.
    """
    derivative = grad(f)

    def gradient(w):
        point = validate_point(w)
        objective_value(f, point)
        try:
            value = derivative(point)
        except Exception as exc:
            raise ValueError(
                "Automatic differentiation failed; use autograd.numpy operations "
                "or supply an analytical gradient"
            ) from exc
        return gradient_value(lambda _: value, point)

    return gradient


def finite_difference_gradient(f, w, epsilon=1e-6):
    """Approximate each partial derivative with symmetric differences."""
    point = validate_point(w)
    if not np.isfinite(epsilon) or epsilon <= 0:
        raise ValueError("epsilon must be finite and positive")
    result = np.empty_like(point)
    for i in range(point.size):
        offset = np.zeros_like(point)
        offset[i] = epsilon
        plus = objective_value(f, point + offset)
        minus = objective_value(f, point - offset)
        if not np.isfinite(plus) or not np.isfinite(minus):
            raise ValueError("Finite differences require finite objective values")
        result[i] = (plus - minus) / (2 * epsilon)
    return result
