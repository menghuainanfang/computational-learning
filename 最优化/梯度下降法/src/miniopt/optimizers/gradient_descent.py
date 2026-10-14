"""Gradient descent for real, differentiable f: R^n -> R.

Iterate w_(k+1) = w_k - alpha_k grad f(w_k), stopping at approximate
stationarity or a numerical failure. Inputs are a callable and a 1-D point;
the output is an OptimizationResult with aligned iteration history.
Example: gradient_descent(lambda w: (w[0]-3)**2, [0.0], alpha=0.1).
"""
import numpy as np
from ..differentiation import (
    validate_point, objective_value, gradient_value, make_gradient,
)
from ..results import IterationHistory, OptimizationResult
from ..step_size import choose_step, validate_strategy, StepSizeError


def gradient_descent(f, w0, step="fixed", alpha=0.1, max_iter=1000,
                     tol=1e-8, alpha_max=6.0, gradient=None):
    """Minimize f using Autograd or a supplied gradient(w).

    Invalid interfaces raise ValueError. Numerical failures return success=False
    at the last recorded iterate. nit counts accepted updates, including the
    last allowed update; its gradient is checked for stationarity as well.
    A custom step callable receives (f, w, direction, zero_based_iteration).
    """
    validate_strategy(step, alpha, alpha_max)
    if not isinstance(max_iter, (int, np.integer)) or max_iter < 0:
        raise ValueError("max_iter must be a nonnegative integer")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be finite and positive")
    w = validate_point(w0)
    derivative = make_gradient(f) if gradient is None else gradient
    value = objective_value(f, w)
    if not np.isfinite(value):
        raise ValueError("The initial objective value must be finite")
    g = gradient_value(derivative, w)
    norm = float(np.linalg.norm(g))
    points, values, norms, steps = [w.copy()], [value], [norm], [0.0]
    success = False
    message = "Maximum iterations reached"

    for k in range(max_iter + 1):
        if not np.isfinite(norm):
            message = "Nonfinite gradient at the recorded iterate"
            break
        if norm <= tol:
            success = True
            message = "Gradient norm meets tolerance (approximate stationarity)"
            break
        if k == max_iter:
            break
        try:
            length = choose_step(step, f, w, -g, k, alpha, alpha_max)
        except StepSizeError as exc:
            message = str(exc)
            break
        if length == 0:
            message = "No positive improving step found"
            break
        candidate = w - length * g
        if not np.all(np.isfinite(candidate)):
            message = "Nonfinite trial iterate"
            break
        if np.array_equal(candidate, w):
            message = "Step produced no change in the iterate"
            break
        candidate_value = objective_value(f, candidate)
        if not np.isfinite(candidate_value):
            message = "Nonfinite trial objective value"
            break
        candidate_gradient = gradient_value(derivative, candidate)
        w, value, g = candidate, candidate_value, candidate_gradient
        norm = float(np.linalg.norm(g))
        points.append(w.copy())
        values.append(value)
        norms.append(norm)
        steps.append(length)

    history = IterationHistory(np.asarray(points), np.asarray(values),
                               np.asarray(norms), np.asarray(steps))
    return OptimizationResult(w.copy(), value, len(steps) - 1,
                              success, message, history)
