"""在同一个一维函数上比较四种步长配置。

函数来自原型中的 Machine Learning Refined Exercise 3.5。
输入：相同的初始点 [2]；输出：各次结果和收敛曲线。
运行：python -m 实例.步长比较。
"""
from 模块.自动微分 import make_gradient
from 模块.梯度下降 import gradient_descent
from 模块.绘图 import compare_results


def objective(w):
    x = w[0]
    return (x**4 + x**2 + 10*x) / 50.0


if __name__ == "__main__":
    grad_f = make_gradient(objective)
    configurations = {
        "fixed alpha=1": ("fixed", 1.0),
        "fixed alpha=0.1": ("fixed", 0.1),
        "diminishing alpha=1": ("diminishing", 1.0),
        "bounded search": ("bounded_search", 1.0),
    }

    results = {}
    for name, (step, alpha) in configurations.items():
        result = gradient_descent(
            objective,
            w0=[2.0],
            gradient=grad_f,
            step=step,
            alpha=alpha,
        )
        results[name] = result
        print(name)
        print("  最终参数：", result.x, "函数值：", result.fun)
        print("  更新次数：", result.nit, "停止原因：", result.message)

    compare_results(results, show=True)
