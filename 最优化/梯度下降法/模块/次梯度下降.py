"""次梯度迭代：用用户给定的次梯度更新不可微目标函数。

输入：f、初始向量、次梯度函数 subgradient 和固定或递减步长配置。
输出：IterationHistory，保存完整 max_iter 次更新及初始点。
公式：w_next = w - alpha_k * subgradient(w)。
例：subgradient_descent(f, [2.0], subgradient=np.sign, alpha=0.5)。
"""
import numpy as np
from .自动微分 import validate_point, objective_value, gradient_value
from .步长策略 import choose_step
from .结果记录 import IterationHistory


def subgradient_descent(f, w0, *, subgradient, step="fixed", alpha=0.5, max_iter=20):
    """按指定次数更新，返回历史；不使用普通梯度求解器的驻点停止规则。

    f 接收向量，subgradient 返回同形向量。次梯度由用户明确提供。
    即使某步选定次梯度为零，也继续记录后续迭代，方便比较步长行为。
    """
    if step not in ("fixed", "diminishing"):
        raise ValueError("Subgradient iteration supports fixed or diminishing steps")
    if not np.isfinite(alpha) or alpha <= 0:
        raise ValueError("alpha must be finite and positive")
    if not isinstance(max_iter, (int, np.integer)) or max_iter < 0:
        raise ValueError("max_iter must be a nonnegative integer")
    w = validate_point(w0)
    points, values, norms, steps = [], [], [], []

    for k in range(max_iter + 1):
        value = objective_value(f, w)
        g = gradient_value(subgradient, w)
        norm = float(np.linalg.norm(g))
        if not np.isfinite(value) or not np.isfinite(norm):
            raise ValueError("Objective and subgradient must be finite")
        points.append(w.copy())
        values.append(value)
        norms.append(norm)
        if k == 0:
            steps.append(0.0)
        if k == max_iter:
            break

        length = choose_step(step, f, w, -g, k, alpha=alpha)
        w = validate_point(w - length * g)
        steps.append(length)

    return IterationHistory(
        np.asarray(points), np.asarray(values), np.asarray(norms),
        np.asarray(steps), np.arange(1, len(values) + 1),
    )
