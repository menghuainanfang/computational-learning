"""Choose positive step lengths along a descent direction.

Fixed steps use alpha; diminishing steps use alpha/(k+1). Bounded numerical
search approximately minimizes f(w + a*d) over [0, alpha_max].
choose_step returns a scalar; StepSizeError signals numerical search failure.
Example: choose_step('fixed', f, w, -g, 0, alpha=0.1).
"""
import numpy as np
from scipy.optimize import minimize_scalar
from .differentiation import objective_value


class StepSizeError(RuntimeError):
    """A step strategy could not produce a usable numerical step."""


def validate_strategy(step, alpha, alpha_max):
    """Check strategy configuration before optimization begins."""
    if not callable(step) and step not in ("fixed", "diminishing", "bounded_search"):
        raise ValueError("step must be fixed, diminishing, bounded_search, or callable")
    for name, value in (("alpha", alpha), ("alpha_max", alpha_max)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")


def choose_step(step, f, w, direction, k, alpha=0.1, alpha_max=6.0):
    """Return a step length; custom strategies receive (f, w, direction, k).

    Search checks both interval endpoints as well as SciPy's interior result.
    A zero return means no improving step was found.
    """
    validate_strategy(step, alpha, alpha_max)
    if not isinstance(k, (int, np.integer)) or k < 0:
        raise ValueError("k must be a nonnegative integer")
    if callable(step):
        candidate = float(step(f, w, direction, k))
    elif step == "fixed":
        candidate = float(alpha)
    elif step == "diminishing":
        candidate = float(alpha / (k + 1))
    else:
        def phi(a):
            value = objective_value(f, w + a * direction)
            if not np.isfinite(value):
                raise StepSizeError("Bounded search encountered a nonfinite objective")
            return value

        try:
            result = minimize_scalar(phi, bounds=(0.0, alpha_max), method="bounded")
            if (not result.success or not np.isfinite(result.x)
                    or not np.isfinite(result.fun) or not 0 <= result.x <= alpha_max):
                raise StepSizeError("Bounded numerical line search failed")
            candidates = [0.0, float(result.x), float(alpha_max)]
            values = [phi(a) for a in candidates]
            candidate = candidates[int(np.argmin(values))]
        except (RuntimeError, FloatingPointError, OverflowError) as exc:
            raise StepSizeError(f"Bounded numerical line search failed: {exc}") from exc
    if not np.isfinite(candidate) or candidate < 0:
        raise StepSizeError("The step size must be finite and nonnegative")
    return candidate
