# 最优化

通过实现优化算法和运行数值实验，学习梯度、步长选择、驻点判定及收敛行为。

## 专题

| 专题 | 内容 | 语言 | 状态 |
| --- | --- | --- | --- |
| [梯度下降法](梯度下降法/README.md) | 通用模块与完整实例：梯度下降、回溯、扩张中点和次梯度实验 | Python | 已验证 |

梯度下降法专题来源于 *Machine Learning Refined: Foundations, Algorithms, and Applications*
相关练习。通用代码放在 `梯度下降法/模块/`，具体问题放在 `梯度下降法/实例/`。
先从 [一维完整实例](梯度下降法/实例/一维函数.py) 开始阅读。

## 运行

从仓库根目录进入专题：

```powershell
cd 最优化/梯度下降法
python -m venv .venv
.venv\Scripts\python.exe -m pip install -e ".[dev]"
.venv\Scripts\python.exe -m 实例.一维函数
.venv\Scripts\python.exe -m pytest
```

依赖在 `pyproject.toml` 中声明。数学说明、验证结果和后续方向见专题的 `docs/`。
