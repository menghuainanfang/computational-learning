"""一维完整实例：求 f(w) = (w^4 + w^2 + 10w)/50 的驻点。

来源：原型 miniopt_v1.py 标注为 Machine Learning Refined 的 Exercise 3.5。
输入：初始点 [2.0]；输出：最终点、函数值、停止原因和两种图。
运行：在“梯度下降法”目录执行 python -m 实例.一维函数。
"""
import autograd.numpy as np

# 这些函数分别来自“模块”文件夹里的三个文件。
from 模块.自动微分 import make_gradient
from 模块.梯度下降 import gradient_descent
from 模块.绘图 import plot_convergence, plot_trajectory_1d


# 第一步：定义要优化的函数。即使只有一个未知数，也用向量表示。
def objective(w):
    x = w[0]
    return (x**4 + x**2 + 10*x) / 50.0


if __name__ == "__main__":
    # 第二步：设定初始点，以及“从函数得到梯度函数”。
    w0 = np.array([2.0])
    grad_f = make_gradient(objective)

    # 第三步：调用求解器。求解器内部会调用“步长策略”模块。
    result = gradient_descent(
        objective,
        w0=w0,
        gradient=grad_f,
        step="fixed",       # 可换成 "diminishing" 或 "bounded_search"
        alpha=1.0,          # 固定步长；递减策略则使用 1.0/(k+1)
        max_iter=1000,      # 最多更新 1000 次
        tol=1e-8,           # 梯度长度小于等于这个数时停止
    )

    # 第四步：读结果。x 是参数向量，fun 是目标函数的值。
    print("最终参数：", result.x)
    print("最终函数值：", result.fun)
    print("更新次数：", result.nit)
    print("是否达到梯度容差：", result.success)
    print("停止原因：", result.message)

    # 第五步：把同一个结果交给绘图模块。
    plot_convergence(result, show=False)
    plot_trajectory_1d(objective, result, show=True)
