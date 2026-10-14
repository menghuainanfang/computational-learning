"""绝对值实例：比较固定步长和递减步长的 20 次次梯度更新。

来源：Machine Learning Refined Exercise 3.6。
f(w)=|w[0]| 在零点不可微；明确选择次梯度 sign(0)=0。
输入：初始点 [2.0] 与 [2.1]；输出：最终值和两组实验曲线。
运行：python -m 实例.绝对值函数。
"""
from pathlib import Path
import numpy as np
from 模块.次梯度下降 import subgradient_descent
from 模块.绘图 import compare_subgradient_runs


# 第一步：定义目标函数和选定的次梯度。
def objective(w):
    return abs(w[0])


def subgradient(w):
    return np.sign(w)


if __name__ == "__main__":
    # 第二步：保留原始初始点，并增加不会恰好落到零点的对照。
    cases = {}
    for w0 in (2.0, 2.1):
        # 第三步：同一函数使用两个步长策略，各完成 20 次更新。
        fixed = subgradient_descent(
            objective, [w0], subgradient=subgradient,
            step="fixed", alpha=0.5, max_iter=20,
        )
        diminishing = subgradient_descent(
            objective, [w0], subgradient=subgradient,
            step="diminishing", alpha=1.0, max_iter=20,
        )
        runs = {"fixed alpha=0.5": fixed, "diminishing alpha=1": diminishing}
        cases[f"w0={w0}"] = runs

        # 第四步：读取历史的最后一条记录。
        print("初始点：", w0)
        for name, history in runs.items():
            print(name, "最终参数：", history.x[-1], "函数值：", history.fun[-1])

    # 第五步：调用绘图模块；可重新生成的图片统一放在 outputs。
    output_dir = Path("outputs")
    output_dir.mkdir(exist_ok=True)
    compare_subgradient_runs(cases, show=True, save=output_dir / "exercise_3_6.png")
