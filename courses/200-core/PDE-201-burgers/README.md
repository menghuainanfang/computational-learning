# PDE-201｜Burgers方程的有限元方法

从自输运与黏性模型出发，建立弱形式、非线性装配与半离散ODE，再比较时间层选择、Picard/Newton及数值误差。先修为泛函分析、数学物理方程、团队有限元培训与FEM-101。

## 阅读

[认识Burgers方程](notes/01-认识Burgers方程.md)先给出一维方程，解释未知量、已知数据、非线性、自对流和黏性各有什么含义。通过曲线传播、变陡和平滑的图像认识解的变化，再给出完整初边值问题，区分二维标量与向量模型，并明确本课采用哪一个模型。

[物质导数与非线性有限元](notes/02-物质导数与非线性有限元.md)在已有模型认识的基础上，解释对流项的来源，建立弱形式与能量估计，随后推导质量矩阵、非线性张量、P1/P2装配与求积。

[时间推进与非线性求解](notes/03-时间推进与非线性求解.md)从半离散ODE构造下一时间层的方程，用完整一步算例连接新系数与有限元函数，再讨论一次冻结、Picard、Newton和离散能量。

[误差与计算证据](notes/04-误差与计算证据.md)区分空间、时间、代数误差，结合理论条件、制造解和固定半离散参考解释实验结果。

先读第一篇，确认能够说明方程描述什么、给定哪些数据、要求解哪个场，再进入推导与数值方法。程序入口见[实验说明](labs/README.md)。第四篇末尾提供自主讨论，第三篇的一步算例可用于对照公式与代码。

## 符号约定

| 符号 | 含义 |
|---|---|
| $d,\boldsymbol{\beta},D_\beta,\nu$ | 空间维数、固定输运方向、方向导数 $\boldsymbol{\beta}\cdot\nabla$、正黏性系数 |
| $u,u_h,\boldsymbol{U}$ | 连续场、有限元场、该场的系数向量 |
| $N_I,N$ | 理论内部自由度数、代码完整全局自由度数（含边界） |
| $\boldsymbol{M},\boldsymbol{K},\boldsymbol{C}(\boldsymbol{U})$ | 质量矩阵、刚度矩阵、非线性对流向量 |
| $\boldsymbol{B}(\boldsymbol{W}),\boldsymbol{G}(\boldsymbol{W})$ | 分别冻结速度因子、方向导数因子的矩阵 |
| $E,\boldsymbol{P}_E,\boldsymbol{U}^E$ | 单元、单元系数抽取映射、局部系数向量 |
| $t_n,\Delta t,k$ | 物理时间层、时间步长、新层内非线性迭代编号 |
| $m,p$ | 实验每方向网格分段数、有限元多项式次数 |
| $(v,w),\|v\|_0,\lvert v\rvert_1$ | $L^2$ 内积、$L^2$ 范数、梯度半范数 |
| $\|\boldsymbol{Z}\|_{\boldsymbol{M}},E^n$ | 质量范数、$\tfrac12(\boldsymbol{U}^n)^T\boldsymbol{M}\boldsymbol{U}^n$ |

正文用 $\varphi_i$ 表示标量基函数，代码中的 `phi` 为其求积点数组；`U`、`M`、`K`、`C`、`B`、`G` 分别对应 $\boldsymbol{U}$、$\boldsymbol{M}$、$\boldsymbol{K}$、$\boldsymbol{C}$、$\boldsymbol{B}$、$\boldsymbol{G}$。粗体对象的分量如 $U_i,M_{ij}$ 仍用非粗体。

单元积分中的 $i,j,k$ 表示基函数索引，作用在测试函数、函数值因子、导数因子上；该局部索引 $k$ 与时间迭代的 $k$ 按段落语境区分。程序的网格分段参数名为 `n`，实验正文改用 $m$ 避免与时间层混淆。矩阵公式默认按齐次边界消元后的内部系统理解；代码保留完整系数布局时会明确说明。

## 运行

从本课程根目录执行：

```text
python tests/验证非线性.py
python labs/运行时间实验.py
```

依赖FEALPy、NumPy、SciPy、Matplotlib，采用NumPy后端。验证环境采用 Python 3.13、FEALPy 3.4.0（源码安装）和 NumPy 2.5.3。配套数据由本目录程序生成。

## 状态与范围

**In Progress**。教材正文与程序、验证结果齐备，真实读者试学与语言审阅继续进行。实际覆盖一、二维标量P1/P2和Euler族；BDF2/RK为理论接口，未实现。空间实验为联合加密，未声称独立测量空间阶。

[实验说明](labs/README.md) · [验证说明](tests/README.md) · [资料来源](references.md)


## 符号规范依据

数学记号遵循[Meshflow 数学符号与领域记号约定](https://github.com/suanhaitech/meshflow_group/blob/main/standards/meshflow_mathematical_notation_standard.md)，向量、矩阵和张量采用 `\boldsymbol{}`，分量使用普通数学斜体。
