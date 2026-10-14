"""梯度下降：把求梯度、选步长和更新参数连接成一个循环。

输入：目标函数 f、初始向量 w0、步长策略和停止条件。
输出：OptimizationResult，包含最终点、停止原因和整个迭代记录。
每一步的数学公式：w_next = w - alpha_k * grad_f(w)。
例：gradient_descent(lambda w: (w[0] - 3)**2, [0.0], alpha=0.1)
"""
import numpy as np
from .自动微分 import (
    validate_point, objective_value, gradient_value, make_gradient,
)
from .结果记录 import IterationHistory, OptimizationResult
from .步长策略 import choose_step, validate_strategy, StepSizeError


def gradient_descent(f, w0, step="fixed", alpha=0.1, max_iter=1000,
                     tol=1e-8, alpha_max=6.0, gradient=None, step_options=None):
    """运行梯度下降，直到梯度足够小、达到迭代上限或遇到数值问题。

    f：目标函数；w0：初始点，例如 [2.0] 或 [2.0, -1.0]。
    step：fixed / diminishing / bounded_search / armijo / expand_midpoint。
    也可以传入步长函数；step_options 配置回溯与扩张策略的额外参数。
    alpha：固定或递减策略的初始步长；alpha_max：搜索区间上界。
    max_iter：最多更新几次；tol：梯度长度的容差。
    gradient：可选的梯度函数；省略时自动调用 make_gradient(f)。
    """
    # 1. 检查输入，准备目标函数的梯度。
    validate_strategy(step, alpha, alpha_max)
    if not isinstance(max_iter, (int, np.integer)) or max_iter < 0:
        raise ValueError("max_iter must be a nonnegative integer")
    if not np.isfinite(tol) or tol <= 0:
        raise ValueError("tol must be finite and positive")
    w = validate_point(w0)

    # 统计经求解器包装的调用，包括搜索和自动构造梯度中的函数调用。
    function_evaluations = 0
    gradient_evaluations = 0

    def counted_f(point):
        nonlocal function_evaluations
        function_evaluations += 1
        return f(point)

    derivative = make_gradient(counted_f) if gradient is None else gradient

    def evaluate_gradient(point):
        nonlocal gradient_evaluations
        gradient_evaluations += 1
        return gradient_value(derivative, point)

    # 2. 计算初始点的函数值和梯度，并记录初始状态。
    value = objective_value(counted_f, w)
    if not np.isfinite(value):
        raise ValueError("The initial objective value must be finite")
    g = evaluate_gradient(w)
    norm = float(np.linalg.norm(g))
    points = [w.copy()]
    values = [value]
    norms = [norm]
    steps = [0.0]
    calls = [function_evaluations]
    success = False
    message = "Maximum iterations reached"

    # 3. 每轮先判断是否停止，再完成一次更新。
    # 多检查一次，是为了判断最后一次更新后的点是否已经满足容差。
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
        # 下降方向是负梯度；步长由“步长策略”模块决定。
        direction = -g
        try:
            length = choose_step(
                step, counted_f, w, direction, k, alpha, alpha_max,
                gradient=g, step_options=step_options,
            )
        except StepSizeError as exc:
            message = str(exc)
            break
        if length == 0:
            message = "No positive improving step found"
            break
        # 核心更新：w_next = w + alpha_k * direction。
        candidate = w + length * direction
        if not np.all(np.isfinite(candidate)):
            message = "Nonfinite trial iterate"
            break
        if np.array_equal(candidate, w):
            message = "Step produced no change in the iterate"
            break
        candidate_value = objective_value(counted_f, candidate)
        if not np.isfinite(candidate_value):
            message = "Nonfinite trial objective value"
            break
        candidate_gradient = evaluate_gradient(candidate)
        # 接受新点，计算新梯度的长度，并将各项记录同步追加。
        w = candidate
        value = candidate_value
        g = candidate_gradient
        norm = float(np.linalg.norm(g))
        points.append(w.copy())
        values.append(value)
        norms.append(norm)
        steps.append(length)
        calls.append(function_evaluations)

    # 4. 把各步记录整理成数组，再返回一个结果对象。
    history = IterationHistory(np.asarray(points), np.asarray(values),
                               np.asarray(norms), np.asarray(steps), np.asarray(calls))
    return OptimizationResult(w.copy(), value, len(steps) - 1,
                              success, message, history,
                              function_evaluations, gradient_evaluations)
