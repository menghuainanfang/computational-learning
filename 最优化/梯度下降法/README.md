# 梯度下降法

这一专题分成两部分：**模块保存通用方法，实例定义具体问题并把模块连接起来。**
先读 `实例/一维函数.py`：它按五个步骤完整展示了一次求解。

## 文件放在哪里

```text
梯度下降法/
├── 模块/                   # 换一道题时可以继续使用
│   ├── 自动微分.py         # 从 f 得到梯度函数 grad_f
│   ├── 步长策略.py         # 决定每一步走多远
│   ├── 回溯搜索.py         # Armijo：不满足充分下降就缩小步长
│   ├── 扩张中点搜索.py     # 扩张，首次失败后只试一次中点
│   ├── 梯度下降.py         # 重复：求梯度 → 选步长 → 更新参数
│   ├── 次梯度下降.py       # 给不可微函数使用明确指定的次梯度
│   ├── 结果记录.py         # 保存答案和每一步的记录
│   ├── 绘图.py             # 把结果画成曲线或轨迹
│   └── __init__.py         # 方便导入常用函数
├── 实例/                   # 一道问题一个文件
│   ├── 一维函数.py         # 从这里开始读
│   ├── 二维函数.py         # 二维二次函数及等高线轨迹
│   ├── 二次型.py           # 五维问题，与解析解比较
│   ├── 步长比较.py         # 同一道题换不同步长
│   ├── 绝对值函数.py       # 20 次次梯度更新，比较固定和递减步长
│   ├── 线搜索比较.py       # 比较回溯与扩张中点的效果和求值工作量
│   └── __init__.py
├── tests/                  # 自动检查程序正确性
├── docs/                   # 数学说明与验证结果
└── pyproject.toml          # 安装依赖和模块
```

## 模块怎样拼起来

1. 实例定义 `objective(w)` 和初始向量 `w0`。
2. `make_gradient(objective)` 返回梯度函数 `grad_f`。
3. 实例把 `objective`、`w0`、`grad_f` 和步长配置交给 `gradient_descent`。
4. 求解器每轮调用 `choose_step`，然后执行 `w_next = w - alpha_k * grad_f(w)`。
5. 求解器返回 `result`；实例打印它，并交给绘图模块。

```mermaid
flowchart LR
    E[实例：目标函数与初始点] --> A[自动微分]
    A -->|梯度函数| G[梯度下降]
    E -->|目标函数、初始点、配置| G
    G -->|请求步长| S[步长策略]
    S -->|步长| G
    G -->|结果记录| P[打印与绘图]
```

建议阅读顺序：`实例/一维函数.py` → `模块/梯度下降.py` →
`模块/步长策略.py` → `模块/自动微分.py`。结果和绘图模块可以在需要时查看。

## 安装和运行

Python 3.10 或以上。在仓库根目录进入专题：

```powershell
cd 最优化/梯度下降法
python -m venv .venv
.venv\Scripts\python.exe -m pip install -e ".[dev]"
.venv\Scripts\python.exe -m 实例.一维函数
```

运行后会打印答案并显示图。其他实例的运行命令是：

```powershell
.venv\Scripts\python.exe -m 实例.二维函数
.venv\Scripts\python.exe -m 实例.二次型
.venv\Scripts\python.exe -m 实例.步长比较
.venv\Scripts\python.exe -m 实例.绝对值函数
.venv\Scripts\python.exe -m 实例.线搜索比较
```

在 macOS / Linux 上，将 `.venv\Scripts\python.exe` 换成 `.venv/bin/python`。
NumPy 用于数组，Autograd 用于求梯度，SciPy 用于有界步长搜索，Matplotlib 用于绘图。

## 换成自己的题目

复制 `实例/一维函数.py`，修改目标函数、初始点和求解配置。例如：

```python
import autograd.numpy as np
from 模块.自动微分 import make_gradient
from 模块.梯度下降 import gradient_descent

# f: R^n → R；w 是一维数组，函数返回一个实数。
def objective(w):
    return np.sum((w - 3.0)**2)

result = gradient_descent(
    objective,
    w0=[2.0, -1.0],
    gradient=make_gradient(objective),
    step="fixed",
    alpha=0.1,
)
print(result.x)
```

目标函数内的数组运算使用 `autograd.numpy`，避免把参与求导的量强制转换成
`float`。如果已经有解析梯度，可以把自己的 `grad_f(w)` 直接传给 `gradient`。
省略 `gradient` 时，求解器也会自动调用 `make_gradient`。

改变参数个数时，轨迹图也要选择对应的维数：一维用 `plot_trajectory_1d`，
二维用 `plot_trajectory_2d`；任意维数都可以使用 `plot_convergence`。

| step 配置 | 步长 | 主要参数 |
| --- | --- | --- |
| `"fixed"` | 每步相同 | `alpha` |
| `"diminishing"` | alpha / (k + 1) | `alpha` |
| `"bounded_search"` | 在 [0, alpha_max] 内数值搜索 | `alpha_max` |
| `"armijo"` | 不满足充分下降就缩小步长 | `alpha`、`step_options` |
| `"expand_midpoint"` | 扩张，首次失败后只试一次中点 | `alpha`、`alpha_max`、`step_options` |

## 回溯与扩张中点怎样使用

`alpha` 是每次线搜索的初始步长。`step_options` 设置缩小比例 `shrink`、
Armijo 系数 `c`、回溯次数 `max_backtracks`；扩张策略还接受扩张倍数
`growth` 和次数上限 `max_expansions`。`alpha_max` 直接传给求解器。

```python
result = gradient_descent(
    objective, w0=[4.0, 3.0],
    step="expand_midpoint", alpha=0.03, alpha_max=1.0,
    max_iter=60, tol=1e-7,
    step_options={"growth": 2.0, "shrink": 0.5, "c": 1e-4},
)
```

扩张中点策略先找到满足 Armijo 条件的步长；然后增大步长，要求同时满足
Armijo 条件和“函数值比上次更低”。第一次失败后，只试一次成功步长与
失败步长之间的中点。若扩张达到上界或次数上限，则返回最后的成功步长。
详细流程见 [线搜索与次梯度](docs/线搜索与次梯度.md)。

## 绝对值实例怎样使用模块

`实例/绝对值函数.py` 定义 `f(w)=|w[0]|`，并明确选择 `np.sign(w)` 为次梯度。
零点的普通导数不存在；取 `sign(0)=0` 是次梯度的一种合法选择。
实例调用 `模块/次梯度下降.py`，分别从 2.0 和 2.1 出发，使用固定步长 0.5
和递减步长 1/(k+1)，每组完整执行 20 次更新，再把历史交给绘图模块。
两组新实例生成的图片保存在 `outputs/`。

## 怎样读结果

| 字段 | 含义 |
| --- | --- |
| `result.x` | 最终参数向量 |
| `result.fun` | 最终函数值 |
| `result.nit` | 实际更新次数 |
| `result.success` | 是否达到梯度容差 |
| `result.message` | 停止原因 |
| `result.history` | 全部迭代记录 |
| `result.function_evaluations` | 求解器计数的目标函数调用次数 |
| `result.gradient_evaluations` | 总梯度函数调用次数 |

历史包括初始点，所以每个记录数组都有 `nit + 1` 项。`history.x` 保存参数，
`history.fun` 保存函数值，`history.gradient_norm` 保存梯度长度，
`history.step_size` 保存走到各点的步长；初始点的步长记为 0。
`history.function_evaluations` 保存记录各点时的累计函数调用次数，与其他
数组长度相同。计数也包括未接受的搜索试探；搜索失败时，总数可能
大于最后一个记录的累计数。自动构造梯度时，计数包含自动微分及其输入检查
中的目标函数调用。函数调用次数衡量求值工作量，不等同于运行时间。

绘图函数默认不显示图，设置 `show=True` 就会显示；设置 `save="图.png"`
可以保存图。`绘图.py` 与求解过程独立。

## 验证

```powershell
.venv\Scripts\python.exe -m pytest
.venv\Scripts\python.exe -m build
```

验证包括解析解对照、梯度核对、五种步长、单次中点规则、次梯度实例、
停止条件、求值计数、历史记录和绘图。
详细结果见 [实现报告](docs/implementation_report.md)，公式说明见
[数学背景](docs/mathematical_background.md)。

## 适用范围与限制

- 普通梯度下降与两种 Armijo 线搜索处理实数、无约束、可微的标量目标函数。
- 次梯度迭代由用户提供次梯度，仅支持固定或递减步长；按指定次数返回历史，
  不使用普通梯度求解器的驻点停止规则。
- 梯度很小表示接近驻点，不保证找到全局最小值。
- 固定或递减步长需要合适的尺度；有界数值搜索也不保证全局最优步长。
- 扩张中点是启发式策略，不保证得到线上全局最佳步长。
- 递减策略在当前一维例题中，1000 次更新后仍未达到梯度容差。
- 一维、二维问题可以画参数轨迹；高维问题可以画收敛曲线。
- 保存全部历史所需的内存随维数和迭代次数增长。
- 传入外部构造的梯度时，其中直接调用目标函数的次数不计入求解器的计数。

## 来源与后续方向

原型来自 *Machine Learning Refined: Foundations, Algorithms, and Applications*
相关练习；其中一维函数在原型中标注为 Exercise 3.5。
自动微分与有界搜索分别使用 Autograd 和 SciPy。
扩张中点是本项目的自定义策略；绝对值实例对应 Exercise 3.6。
后续可以研究上一轮步长的复用和相同求值预算下的比较。
代码采用 MIT 许可。
