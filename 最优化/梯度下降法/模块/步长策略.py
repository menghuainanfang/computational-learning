"""步长策略：决定沿下降方向走多远。

输入：目标函数 f、当前点 w、方向 direction、迭代编号 k 和步长配置。
输出：一个步长 alpha_k；固定步长为 alpha，递减步长为 alpha/(k+1)。
有界搜索近似求解 min f(w + a * direction)，其中 0 <= a <= alpha_max。
例：choose_step("fixed", f, w, -grad_f(w), k=0, alpha=0.1)
"""
import numpy as np
from scipy.optimize import minimize_scalar
from .自动微分 import objective_value


class StepSizeError(RuntimeError):
    """步长计算遇到数值问题，供求解器记录停止原因。"""


def validate_strategy(step, alpha, alpha_max):
    """检查策略名和步长参数。"""
    if not callable(step) and step not in ("fixed", "diminishing", "bounded_search"):
        raise ValueError("step must be fixed, diminishing, bounded_search, or callable")
    for name, value in (("alpha", alpha), ("alpha_max", alpha_max)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{name} must be finite and positive")


def choose_step(step, f, w, direction, k, alpha=0.1, alpha_max=6.0):
    """根据策略返回步长；k 从 0 开始。返回 0 表示没有可用的正步长。"""
    validate_strategy(step, alpha, alpha_max)
    if not isinstance(k, (int, np.integer)) or k < 0:
        raise ValueError("k must be a nonnegative integer")
    if callable(step):
        candidate = float(step(f, w, direction, k))
    elif step == "fixed":
        # 策略一：每次都走相同的步长。
        candidate = float(alpha)
    elif step == "diminishing":
        # 策略二：迭代次数越大，步长越小。
        candidate = float(alpha / (k + 1))
    else:
        # 策略三：把多变量函数限制到下降方向，变成单变量函数 phi(a)。
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
            # 同时比较区间两端和 SciPy 找到的内部候选点。
            candidates = [0.0, float(result.x), float(alpha_max)]
            values = [phi(a) for a in candidates]
            candidate = candidates[int(np.argmin(values))]
        except (RuntimeError, FloatingPointError, OverflowError) as exc:
            raise StepSizeError(f"Bounded numerical line search failed: {exc}") from exc
    if not np.isfinite(candidate) or candidate < 0:
        raise StepSizeError("The step size must be finite and nonnegative")
    return candidate
