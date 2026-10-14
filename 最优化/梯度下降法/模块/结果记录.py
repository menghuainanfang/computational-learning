"""结果记录：保存最终答案和迭代过程。

梯度下降每走一步，就保存参数、函数值、梯度长度和所用步长。
OptimizationResult 保存最终点和停止原因；IterationHistory 保存数组记录。
例：result.x 是最终点，result.history.fun 是每一步的函数值。
历史包含初始点，所以有 nit + 1 条记录；初始点的 step_size 记为 0。
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class IterationHistory:
    x: np.ndarray              # 每一个参数向量，形状为 (记录数, 参数个数)
    fun: np.ndarray            # 对应的函数值
    gradient_norm: np.ndarray  # 对应的梯度长度
    step_size: np.ndarray      # 走到该点时使用的步长；初始记录为 0


@dataclass
class OptimizationResult:
    x: np.ndarray              # 最终参数向量
    fun: float                 # 最终函数值
    nit: int                   # 实际完成的更新次数
    success: bool              # 是否达到梯度容差（不保证全局最优）
    message: str               # 停止原因
    history: IterationHistory  # 全部迭代记录
