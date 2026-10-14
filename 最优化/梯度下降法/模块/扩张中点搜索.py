"""扩张＋单次中点搜索：先回溯，再扩张，最后至多一次中点试探。

输入：f、当前向量、梯度、下降方向，以及初始步长和扩张参数。
输出：一个满足 Armijo 条件的正步长；扩张还要求函数值继续降低。
例：expand_midpoint_search(f, w, g, -g, alpha0=0.18)。
它在首次扩张失败后只试一次中点，不继续二分，也不求精确线最小点。
"""
import numpy as np
from .自动微分 import objective_value
from .回溯搜索 import _check_inputs, _trial, StepSizeError


def expand_midpoint_search(f, w, gradient, direction, *, alpha0=0.1,
                           growth=2.0, shrink=0.5, c=1e-4, alpha_max=10.0,
                           max_backtracks=50, max_expansions=20):
    """先找可接受步长，再扩张；首次失败时只试一次两步长的中点。

    如果初始步长需要回溯，回溯成功后仍会进入扩张阶段。
    达到扩张次数上限或步长上界时，直接返回最后的成功步长。
    """
    w, direction, slope = _check_inputs(
        w, gradient, direction, alpha0, shrink, c, max_backtracks
    )
    if not np.isfinite(growth) or growth <= 1:
        raise ValueError("growth must be finite and > 1")
    if not np.isfinite(alpha_max) or alpha_max < alpha0:
        raise ValueError("alpha_max must be finite and >= alpha0")
    if not isinstance(max_expansions, (int, np.integer)) or max_expansions < 0:
        raise ValueError("max_expansions must be a nonnegative integer")
    f0 = objective_value(f, w)
    if not np.isfinite(f0):
        raise ValueError("f(w) must be finite")

    # 1. 如果初始步长不合适，先用回溯找到一个满足 Armijo 的步长。
    alpha = float(alpha0)
    accepted = False
    for attempt in range(max_backtracks + 1):
        value, accepted = _trial(f, w, direction, f0, slope, c, alpha)
        if accepted:
            break
        if attempt < max_backtracks:
            alpha *= shrink
            if alpha == 0:
                break
    if not accepted:
        raise StepSizeError("No Armijo-acceptable initial step found")

    best_alpha, best_value = alpha, value

    # 2. 扩张步长；同时满足 Armijo 和“比上次更低”才接受。
    for _ in range(max_expansions):
        proposed = min(alpha_max, best_alpha * growth)
        if proposed <= best_alpha:
            break
        proposed_value, proposed_ok = _trial(
            f, w, direction, f0, slope, c, proposed
        )
        if proposed_ok and proposed_value < best_value:
            best_alpha, best_value = proposed, proposed_value
            continue

        # 3. 首次失败后，只试一次成功步长与失败步长的算术中点。
        midpoint = best_alpha + 0.5 * (proposed - best_alpha)
        if best_alpha < midpoint < proposed:
            mid_value, mid_ok = _trial(f, w, direction, f0, slope, c, midpoint)
            if mid_ok and mid_value < best_value:
                best_alpha = midpoint
        break

    return best_alpha
