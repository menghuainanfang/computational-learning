"""Armijo 回溯搜索：从初始步长开始，逐次缩小直到满足充分下降。

输入：f、当前向量 w、梯度 gradient、下降方向 direction，以及搜索参数。
输出：满足 f(w + alpha*d) <= f(w) + c*alpha*(gradient @ d) 的正步长。
例：armijo_backtracking(f, w, grad_f(w), -grad_f(w), alpha0=1.0)。
"""
import numpy as np
from .自动微分 import objective_value


class StepSizeError(RuntimeError):
    """在有限的搜索次数内没有找到可用步长。"""


def _check_inputs(w, gradient, direction, alpha0, shrink, c, max_backtracks):
    """检查搜索参数，并计算方向导数 gradient @ direction。"""
    w = np.asarray(w, dtype=float)
    gradient = np.asarray(gradient, dtype=float)
    direction = np.asarray(direction, dtype=float)
    if w.ndim != 1 or not w.size or w.shape != gradient.shape or w.shape != direction.shape:
        raise ValueError("w, gradient, and direction must be equal-sized 1-D vectors")
    if not (np.isfinite(w).all() and np.isfinite(gradient).all() and np.isfinite(direction).all()):
        raise ValueError("w, gradient, and direction must be finite")
    if not np.isfinite(alpha0) or alpha0 <= 0:
        raise ValueError("alpha0 must be finite and positive")
    if not (0 < shrink < 1 and 0 < c < 1):
        raise ValueError("0 < shrink < 1 and 0 < c < 1 are required")
    if not isinstance(max_backtracks, (int, np.integer)) or max_backtracks < 0:
        raise ValueError("max_backtracks must be a nonnegative integer")
    slope = float(np.dot(gradient, direction))
    if not np.isfinite(slope) or slope >= 0:
        raise ValueError("direction must be a strict descent direction")
    return w, direction, slope


def _trial(f, w, direction, f0, slope, c, alpha):
    """试一个步长；非有限函数值视为本次不满足 Armijo 条件。"""
    with np.errstate(over="ignore", divide="ignore", invalid="ignore"):
        value = objective_value(f, w + alpha * direction)
    threshold = f0 + c * alpha * slope
    return value, np.isfinite(value) and value <= threshold


def armijo_backtracking(f, w, gradient, direction, *, alpha0=1.0,
                        shrink=0.5, c=1e-4, max_backtracks=50):
    """只缩小步长，找到第一个满足 Armijo 条件的步长就返回。

    max_backtracks 是允许缩小的次数，初始步长也会试一次。
    该条件用于可微函数的严格下降方向。
    """
    w, direction, slope = _check_inputs(
        w, gradient, direction, alpha0, shrink, c, max_backtracks
    )
    f0 = objective_value(f, w)
    if not np.isfinite(f0):
        raise ValueError("f(w) must be finite")

    alpha = float(alpha0)
    for attempt in range(max_backtracks + 1):
        _, accepted = _trial(f, w, direction, f0, slope, c, alpha)
        if accepted:
            return alpha
        if attempt < max_backtracks:
            alpha *= shrink
            if alpha == 0:
                break
    raise StepSizeError("Armijo backtracking exhausted its trial budget")
