"""二维实例：最小化 f(w) = 0.5*w^T*A*w - b^T*w。

输入：初始点 [2, 3]；输出：最终点和二维等高线轨迹。
A 对称正定，解析解为 A*w=b。运行：python -m 实例.二维函数。
"""
import autograd.numpy as np
from 模块.自动微分 import make_gradient
from 模块.梯度下降 import gradient_descent
from 模块.绘图 import plot_trajectory_2d


A = np.array([[3.0, 1.0], [1.0, 2.0]])
b = np.array([1.0, 2.0])


def objective(w):
    return 0.5 * np.dot(w, np.dot(A, w)) - np.dot(b, w)


if __name__ == "__main__":
    grad_f = make_gradient(objective)
    result = gradient_descent(
        objective,
        w0=[2.0, 3.0],
        gradient=grad_f,
        step="bounded_search",
        alpha_max=6.0,
    )
    print("数值结果：", result.x)
    print("解析解：", np.linalg.solve(A, b))
    print("停止原因：", result.message)
    plot_trajectory_2d(objective, result, show=True)
