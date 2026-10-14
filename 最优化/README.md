# 最优化

通过实现优化算法和运行数值实验，学习梯度、步长选择、驻点判定及收敛行为。

## 专题

| 专题 | 内容 | 语言 | 状态 |
| --- | --- | --- | --- |
| [梯度下降法](梯度下降法/README.md) | MiniOpt：通用梯度下降、自动微分、步长策略、轨迹绘图与测试 | Python | 已验证 |

梯度下降法专题源于 *Machine Learning Refined: Foundations, Algorithms, and Applications* 的练习。可复用算法与教材例题分别放在 `src/miniopt` 和 `examples` 中；数学说明见专题下的 `docs/mathematical_background.md`。

## 运行

从仓库根目录进入专题，安装独立环境后运行：

```powershell
cd 最优化/梯度下降法
python -m venv .venv
.venv\Scripts\python.exe -m pip install -e ".[dev]"
.venv\Scripts\python.exe -m pytest
.venv\Scripts\python.exe examples/gradient_descent_1d.py --show
```

依赖声明在专题的 `pyproject.toml` 中。实验使用 Python 3.13.2 验证；详细结果、数值限制和后续方向见 [实现报告](梯度下降法/docs/implementation_report.md)。
