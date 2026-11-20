# FEM-101｜稳态反应–对流–扩散方程

本课从控制体守恒建立模型，连接弱形式、适定性、有限元装配、FEALPy实现和收敛验证。先修为泛函分析、数学物理方程及团队有限元基础培训。

## 学习目标与阅读入口

读者应能解释守恒与非守恒形式的关系，推导带Dirichlet边界的弱形式，说明适定性与误差估计所需条件，并将单元积分、自由度映射和边界消元对应到FEALPy实现。

1. [稳态输运模型与弱问题](notes/01-稳态输运模型与弱问题.md)：从守恒推导模型，解释变系数与材料界面，处理边界与弱形式，证明适定性。
2. [三角形有限元与全局装配](notes/02-三角形有限元与全局装配.md)：构造离散空间，连接参考单元、一般次数、变系数求积、全局编号和FEALPy。
3. [有限元误差与收敛验证](notes/03-有限元误差与收敛验证.md)：推导Céa与伴随误差估计，区分代数误差和离散误差，解释强对流振荡并核对实际收敛数据。
4. [FEALPy求解流程](notes/04-FEALPy求解流程.md)：按步骤运行数据、空间、装配、边界、求解和误差检查。

- [实验与自主问题](labs/README.md)
- [验证内容与结果](tests/README.md)
- [资料与实现来源](references.md)

前三篇按模型、离散、误差的依赖关系阅读，第四篇串联完整程序流程。公式沿用全课连续编号，代码与图像邻近相关推导。

## 主要符号

| 符号 | 约定 |
|---|---|
| $a,\boldsymbol{\beta},c;f,g$ | 扩散、速度、反应；源项与Dirichlet数据 |
| $\kappa$ | 守恒模型中的物理反应率，$c=\kappa+\nabla\cdot\boldsymbol{\beta}$ |
| $V,V_g$ | $H_0^1(\Omega)$；具有边界迹$g$的仿射空间 |
| $S_h,V_h$ | 连续分片一次空间；其零迹子空间 |
| $\mathcal{A},L$ | 连续双线性型与载荷泛函；离散全矩阵和向量另记$\boldsymbol{A},\boldsymbol{F}$ |
| $K,\boldsymbol{J}_K$ | 三角单元；参考到物理坐标映射的雅可比矩阵 |
| $\boldsymbol{D}^K,\boldsymbol{C}^K,\boldsymbol{R}^K$ | 单元扩散、对流、反应矩阵 |
| $i,j;I,J$ | 局部测试、试探编号；对应全局行、列编号 |
| $\mathcal I,\mathcal{D}$ | 内部与边界自由度集合 |
| $\|v\|_V,\lvert v\rvert_{H^1}$ | 均为零迹空间上的梯度范数$\|\nabla v\|_{L^2}$；在一般$H^1$中后者为半范数 |

完整范数$\|v\|_{H^1}^2=\|v\|_{L^2}^2+\lvert v\rvert_{H^1}^2$；$L^2$内积记为$(w,v)_{L^2}=\int_\Omega wv$。正则性、求积与几何条件在对应结论附近给出。各专题沿用全课公式编号1—29，新增解释中的未编号公式在本节就近引用。

## 运行

从本课程目录运行：

```text
python labs/稳态输运实验.py
```

依赖FEALPy、NumPy、SciPy和Matplotlib，使用NumPy后端。验证环境采用 Python 3.13、FEALPy 3.4.0（源码安装）、NumPy、SciPy 和 Matplotlib。其他 FEALPy 版本需核对接口。

## 范围与状态

状态：**In Progress**。教材正文和配套实验已具备可复现材料，语言修订与真实读者试学继续进行。

实现范围为二维、常系数、连续P1元及全Dirichlet边界。连续理论介绍一般维数，数值验证不代表一般维数与高阶实现。多分量和时间推进分别衔接FEM-102与TIM-101/FEM-103。

课程目录采用 README、notes、labs、tests、references.md 结构。对应团队目录为 `courses/100-foundations/FEM-101-reaction-advection-diffusion/`。


## 符号规范依据

数学记号遵循[Meshflow 数学符号与领域记号约定](https://github.com/suanhaitech/meshflow_group/blob/main/standards/meshflow_mathematical_notation_standard.md)，向量、矩阵和张量采用 `\boldsymbol{}`，分量使用普通数学斜体。

## 综合实践

[FEALPy求解流程](notes/04-FEALPy求解流程.md)将模型数据、网格与空间、积分装配、边界处理、线性求解和误差检查连成完整的七步流程，可在阅读第二篇后对照运行。

[通用模型与求解器](tests/收束求解器/README.md)提供维数、次数与方程数据的独立练习；配套验证及运行方法见[验证说明](tests/README.md)。
