"""自动微分：把目标函数 f 变成梯度函数 grad_f。

f 接收一维参数向量，返回一个实数；grad_f 接收同样的向量，返回梯度向量。
数学含义：梯度的第 i 个分量是 f 对 w[i] 的偏导数。
例：grad_f = make_gradient(lambda w: w[0]**2)
    grad_f([2.0])  # 得到 [4.0]
"""
import numpy as np
from autograd import grad


def validate_point(w):
    """检查参数是否为有限实数构成的一维向量，并复制一份。"""
    if np.iscomplexobj(w):
        raise ValueError("The parameter vector must be real-valued")
    point = np.asarray(w, dtype=float).copy()
    if point.ndim != 1 or point.size == 0:
        raise ValueError("The parameter vector must be a nonempty 1-D array")
    if not np.all(np.isfinite(point)):
        raise ValueError("The parameter vector must contain finite values")
    return point


def objective_value(f, w):
    """计算 f(w)，检查返回的是一个实数，而不是数组。"""
    value = np.asarray(f(w))
    if value.shape != () or np.iscomplexobj(value):
        raise ValueError("The objective must return a real scalar")
    return float(value)


def gradient_value(gradient, w):
    """计算梯度，检查它的形状是否与参数向量相同。"""
    value = np.asarray(gradient(w))
    if value.shape != w.shape or np.iscomplexobj(value):
        raise ValueError("The gradient must be real and match the parameter shape")
    return np.asarray(value, dtype=float)


def make_gradient(f):
    """返回可调用的梯度函数；目标函数中的数组运算请使用 autograd.numpy。"""
    # Autograd 的 grad(f) 本身返回一个函数，并不是某一点的梯度值。
    derivative = grad(f)

    def gradient(w):
        point = validate_point(w)
        objective_value(f, point)
        try:
            return gradient_value(derivative, point)
        except Exception as exc:
            raise ValueError(
                "Automatic differentiation failed; use autograd.numpy operations "
                "or supply an analytical gradient"
            ) from exc

    return gradient


def finite_difference_gradient(f, w, epsilon=1e-6):
    """用 [f(w+h*e_i) - f(w-h*e_i)] / (2h) 近似偏导数，供核对使用。"""
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
