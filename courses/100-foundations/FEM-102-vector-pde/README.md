# FEM-102｜多分量耦合系统

在同一张网格上求解多个互相影响的标量场，理解乘积空间、分量排列、非对角块和整体边界消元。先修为泛函分析、数学物理方程和团队有限元训练，以及FEM-101的标量装配。

## 专题与阅读顺序

| 专题 | 所解决的问题 | 与实现的联系 |
|---|---|---|
| [耦合模型与向量弱问题](notes/01-耦合模型与向量弱问题.md) | 交换反应怎样形成系统，边界提升与适定性为何成立 | 确定反应矩阵、制造源项和边界数据 |
| [乘积空间与块装配](notes/02-乘积空间与块装配.md) | 分量基函数怎样产生块矩阵，编号和边界怎样作用于全局系统 | TensorFunctionSpace、BlockForm和整体消元 |
| [求解实现与误差验证](notes/03-求解实现与误差验证.md) | 怎样阅读完整程序，误差阶从何而来，哪些证据支持结果 | 求值与积分、结构核对、加密表和收敛图 |

三篇共用同一个非对称双分量模型。第一篇建立连续问题，第二篇构造离散系统，第三篇评价计算结果；图像、程序和数据随相关主题连接。
- [实验与自主问题](labs/README.md)
- [验证内容](tests/README.md)
- [资料来源](references.md)

从课程根目录执行 `python labs/双分量耦合实验.py`。依赖FEALPy、NumPy、SciPy、Matplotlib，采用NumPy后端。验证环境采用 Python 3.13、FEALPy 3.4.0（源码安装）。

## 实际结果

n=8、16、32、64四档网格，双分量总自由度为162、578、2178、8450。最后一级整体L²误差观测阶1.99866，H¹半范数误差观测阶0.99934；边界误差均为0，内部归一化残差最大1.95e-13。各分量误差见CSV，装配、置换与解耦核对见JSON。

## 范围与状态

**In Progress**。完整正文、三幅图与实际可复现实验已形成，真实读者试学与语言审阅继续进行。本章程序仅覆盖二维、双分量、常系数、P1元和全Dirichlet条件；一般分量数理论不等于任意维数或高阶实现均已完成。对应团队目录为 `courses/100-foundations/FEM-102-vector-pde/`。


## 符号约定

| 符号 | 含义 |
|---|---|
| $d,m,N$ | 空间维数、分量数、每分量标量全局自由度数 |
| $r,s=1,\ldots,m$ | 数学分量编号；代码索引为 `r0=r-1`、`s0=s-1` |
| $i,j$ | 单元局部基函数编号；标量全局编号由 $m_K(i)$ 给出 |
| $\boldsymbol{C},\boldsymbol{C}_{\mathrm{sym}}$ | 反应矩阵及其对称部分，$\boldsymbol{C}_{\mathrm{sym}}=(\boldsymbol{C}+\boldsymbol{C}^T)/2$ |
| $\boldsymbol{A}_{rs},\boldsymbol{A}_K^{rs}$ | 全局分量块、单元分量块，块行测试、块列试探 |
| $\boldsymbol{K}_a,\boldsymbol{B},\boldsymbol{M},c_A$ | 标量扩散矩阵、对流矩阵、质量矩阵、双线性型连续性常数 |
| $\boldsymbol{u},\nabla\boldsymbol{u}$ | 向量场、分量×空间方向的梯度矩阵 |
| $\sigma,\tau$ | 制造解中的正弦乘积及其给定速度方向导数 |
| $I,D$ | 内部与Dirichlet自由度集合，下标出现在矩阵上时表示相应子矩阵 |


## 符号规范依据

数学记号遵循[Meshflow 数学符号与领域记号约定](https://github.com/suanhaitech/meshflow_group/blob/main/standards/meshflow_mathematical_notation_standard.md)，向量、矩阵和张量采用 `\boldsymbol{}`，分量使用普通数学斜体。
