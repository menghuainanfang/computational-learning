"""可复用的梯度下降模块。

输入：目标函数 f 和参数向量 w0。输出：最终点以及每一步的记录。
核心公式：w_next = w - alpha * grad_f(w)。
例：from 模块 import gradient_descent
    result = gradient_descent(lambda w: (w[0] - 3)**2, [0.0])
"""
from .自动微分 import make_gradient, finite_difference_gradient
from .梯度下降 import gradient_descent
from .结果记录 import OptimizationResult, IterationHistory

__all__ = ["gradient_descent", "make_gradient", "finite_difference_gradient",
           "OptimizationResult", "IterationHistory"]
