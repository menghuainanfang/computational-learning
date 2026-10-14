"""线搜索实例：比较 Armijo 回溯与扩张＋单次中点策略。

输入：光滑二次函数 f(w)=0.5*(w[0]^2+10*w[1]^2)，初始点 [4,3]。
输出：函数值、更新次数、函数/梯度调用量和四种对比曲线。
运行：python -m 实例.线搜索比较。
"""
from pathlib import Path
from 模块.梯度下降 import gradient_descent
from 模块.绘图 import compare_line_search


# 第一步：定义目标函数；自动微分由求解器构造，以统计其中的函数调用。
def objective(w):
    return 0.5 * (w[0]**2 + 10.0 * w[1]**2)


if __name__ == "__main__":
    # 第二步：给两种策略设置相同的初始点、初始步长和更新次数。
    configurations = {
        "Armijo backtracking": ("armijo", {"c": 1e-4, "shrink": 0.5}),
        "Expansion + one midpoint": (
            "expand_midpoint",
            {"c": 1e-4, "shrink": 0.5, "growth": 2.0, "max_expansions": 12},
        ),
    }

    # 第三步：选择策略后，求解器内部调用对应线搜索模块。
    results = {}
    for name, (step, options) in configurations.items():
        result = gradient_descent(
            objective,
            w0=[4.0, 3.0],
            step=step,
            alpha=0.03,
            alpha_max=1.0,
            max_iter=60,
            tol=1e-7,
            step_options=options,
        )
        results[name] = result

        # 第四步：函数调用次数反映求值工作量，不等同于运行时间。
        print(name)
        print("  最终参数：", result.x, "函数值：", result.fun)
        print("  更新次数：", result.nit, "停止原因：", result.message)
        print("  函数调用：", result.function_evaluations,
              "梯度调用：", result.gradient_evaluations)

    # 第五步：比较按迭代次数和按函数调用量衡量的收敛过程。
    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)
    compare_line_search(results, show=True, save=output_dir / "compare_line_search.png")
