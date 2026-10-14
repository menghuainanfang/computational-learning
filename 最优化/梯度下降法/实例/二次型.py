"""五维实例：比较梯度下降结果与线性方程 A*w=b 的解析解。

输入：正定对角矩阵 A 和向量 b；输出：解的误差和收敛曲线。
高维问题沿用相同接口。运行：python -m 实例.二次型。
"""
import autograd.numpy as np
from 模块.自动微分 import make_gradient
from 模块.梯度下降 import gradient_descent
from 模块.绘图 import plot_convergence


A = np.diag(np.arange(1.0, 6.0))
b = np.ones(5)


def objective(w):
    return 0.5 * np.dot(w, np.dot(A, w)) - np.dot(b, w)


if __name__ == "__main__":
    result = gradient_descent(
        objective,
        w0=np.zeros(5),
        gradient=make_gradient(objective),
        step="fixed",
        alpha=0.2,
    )
    exact_solution = np.linalg.solve(A, b)
    error = np.linalg.norm(result.x - exact_solution)
    print("数值结果：", result.x)
    print("解析解：", exact_solution)
    print("与解析解的距离：", error)
    print("停止原因：", result.message)
    plot_convergence(result, show=True)
