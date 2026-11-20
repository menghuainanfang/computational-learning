# FEALPy求解流程

[课程入口](../README.md) · [单元与装配](02-三角形有限元与全局装配.md) · [误差理论与验证](03-有限元误差与收敛验证.md)

前面的推导已经给出弱问题和离散矩阵。本篇将它们按程序实际运行的顺序连起来：定义数据，建立网格与空间，计算各项积分，处理边界，求解系数，最后评估得到的有限元函数。完整实验见[稳态输运实验.py](../labs/稳态输运实验.py)。

## 第一步：准备模型数据

区域为单位正方形，扩散和反应系数均为1，给定速度为 $\boldsymbol\beta=(1,2)^T$。用

$$
u_\star(x,y)=1+\sin(\pi x)\sin(\pi y)
$$

制造右端。于是边界值为1，而不是0。将它代入 $-\Delta u+\boldsymbol\beta\cdot\nabla u+u=f$，得到

$$
f=2\pi^2\sin(\pi x)\sin(\pi y)
+\pi\cos(\pi x)\sin(\pi y)
+2\pi\sin(\pi x)\cos(\pi y)+u_\star.
$$

`exact`、`gradient`、`source` 分别给出真解、真解梯度及源项；`velocity` 在每个输入点返回二维速度。点数组的末轴是空间坐标，标量数据去掉这一轴，速度与梯度保留这一轴。`@cartesian` 仅声明输入是物理坐标，并不自动计算导数。

## 第二步：建立网格与函数空间

`TriangleMesh.from_box` 划分区域；`LagrangeFESpace(mesh,p=1)` 在这张网格上定义连续分片一次函数。网格提供几何和单元连接关系，空间提供基函数、梯度、自由度及 `cell_to_dof()`。

`p=1` 决定近似空间，积分器的 `q=3` 决定装配使用的求积规则。它们不能互相替代。当前常系数P1双线性项的最高多项式次数为2，三阶规则足以积分；正弦源项仍属于数值积分近似。

## 第三步：将弱形式的各项交给积分器

| 弱形式中的项 | 对应积分器 | 装配对象 |
|---|---|---|
| $\int a\nabla u_h\cdot\nabla v_h$ | `ScalarDiffusionIntegrator` | 扩散矩阵 |
| $\int(\boldsymbol\beta\cdot\nabla u_h)v_h$ | `ScalarConvectionIntegrator` | 对流矩阵 |
| $\int cu_hv_h$ | `ScalarMassIntegrator` | 反应矩阵 |
| $\int fv_h$ | `ScalarSourceIntegrator` | 载荷向量 |

扩散需要两侧基函数的梯度；对流只对试探侧求梯度，因此一般得到非对称矩阵；反应直接使用两侧基函数值。`ScalarMassIntegrator` 计算质量型积分结构，本课以它表达反应，并没有时间导数。

积分器先计算单元贡献。`BilinearForm` 管理前三项，`LinearForm` 管理源项；调用 `assembly()` 后，才按局部到全局编号累加成 $\boldsymbol A,\boldsymbol F$。此时尚未施加Dirichlet条件。

## 第四步：先处理边界，再求系数

把自由度分成内部 $I$ 与边界 $D$，数学上需要求解

$$
\boldsymbol A_{II}\boldsymbol U_I
=\boldsymbol F_I-\boldsymbol A_{ID}\boldsymbol g_D.
$$

`DirichletBC(space,gd=exact).apply(A,F)` 完成相应的右端修正，并在完整布局中用边界方程约束已知系数，返回 `Ab,Fb`。只在求解完成后把边界节点改成1，会漏掉边界对内部方程的贡献。

代码保留原来的 `A,F`，以便之后检验原系统的内部行。被边界条件替换后的边界行不再要求满足原弱方程。

## 第五步：恢复函数并选择线性求解方法

`spsolve(Ab,Fb,solver='scipy')` 求得系数向量；`uh[:]` 保存这些系数。空间和系数一起表示

$$
u_h=\sum_jU_j\varphi_j.
$$

`uh.value(bc)` 与 `uh.grad_value(bc)` 分别在求积点计算函数和梯度。节点系数与物理空间中的函数不同：误差积分、通量和导数都需要通过基函数重建。

本例对流矩阵非对称，使用稀疏直接法。对称正定系统可用CG，一般非对称系统可考虑GMRES及预条件；稀疏本身并不是使用CG的充分条件。程序不显式求逆，因为目标是解一个线性系统，形成逆矩阵会增加存储和计算。

## 第六步：检验方程、边界与函数误差

内部归一化残差检验 $(\boldsymbol A\boldsymbol U-\boldsymbol F)_I$；边界误差检验 $\boldsymbol U_D-\boldsymbol g_D$。它们回答代数方程及约束是否满足。与连续解的差异则用 $L^2$ 范数和梯度半范数积分评估。

误差评价采用独立的 `quadrature_formula(5)`，与装配用的 `q=3` 分开。权重归一化为和1后，只需乘物理单元面积，不能重复乘参考三角形面积。小残差不等于小离散误差；应结合逐级加密与理论阶判断。

## 一份可以直接运行的实现

以下代码可保存为独立Python文件运行，也与实验中的 `求解` 函数逐行对应。它使用同一个非零边界制造解，不依赖前文未定义的 `space` 或 `source`。

```python
"""稳态标量 CDR：FEALPy 装配、Dirichlet 消元及积分误差。"""
import numpy as np
from fealpy.decorator import cartesian
from fealpy.mesh import TriangleMesh
from fealpy.functionspace import LagrangeFESpace
from fealpy.fem import (BilinearForm, LinearForm, ScalarDiffusionIntegrator,
                        ScalarConvectionIntegrator, ScalarMassIntegrator,
                        ScalarSourceIntegrator, DirichletBC)
from fealpy.solver import spsolve
from fealpy.backend import backend_manager as bm


def 求解(n=8):
    """单位正方形上的稳态 CDR，a=c=1、b=(1,2)、边界值为1。"""
    bm.set_backend('numpy')
    mesh = TriangleMesh.from_box([0, 1, 0, 1], nx=n, ny=n)
    space = LagrangeFESpace(mesh, p=1)

    @cartesian
    def exact(p):
        x, y = p[..., 0], p[..., 1]
        return np.sin(np.pi*x) * np.sin(np.pi*y) + 1.

    @cartesian
    def gradient(p):
        x, y = p[..., 0], p[..., 1]
        return np.stack((np.pi*np.cos(np.pi*x)*np.sin(np.pi*y),
                         np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)), axis=-1)

    @cartesian
    def source(p):
        u = exact(p)
        g = gradient(p)
        return 2*np.pi**2*(u-1.) + g[..., 0] + 2*g[..., 1] + u

    form = BilinearForm(space)
    form.add_integrator(ScalarDiffusionIntegrator(coef=1., q=3))
    @cartesian
    def velocity(p):
        return np.broadcast_to(np.array([1., 2.]), p.shape)
    form.add_integrator(ScalarConvectionIntegrator(coef=velocity, q=3))
    form.add_integrator(ScalarMassIntegrator(coef=1., q=3))
    A = form.assembly()
    load = LinearForm(space)
    load.add_integrator(ScalarSourceIntegrator(source=source, q=3))
    F = load.assembly()
    # 保存原始系统；强制边界后，边界行不再代表原来的平衡方程。
    boundary_mask = space.is_boundary_dof()
    Ab, Fb = DirichletBC(space, gd=exact).apply(A, F)
    uh = space.function()
    uh[:] = spsolve(Ab, Fb, solver="scipy")
    interior = ~boundary_mask
    residual = np.linalg.norm((A @ uh[:] - F)[interior]) / max(1., np.linalg.norm(F[interior]))
    points = space.interpolation_points()
    boundary_error = np.max(np.abs(uh[:][boundary_mask] - exact(points)[boundary_mask]))
    # 误差积分独立采用更高精度，以免影响观测阶。
    bc, weights = mesh.quadrature_formula(5).get_quadrature_points_and_weights()
    points_q = mesh.bc_to_point(bc)
    error = uh.value(bc) - exact(points_q)
    grad_error = uh.grad_value(bc) - gradient(points_q)
    measures = mesh.entity_measure("cell")
    l2 = np.sqrt(np.einsum('cq,q,c->', error**2, weights, measures))
    h1 = np.sqrt(np.einsum('cqd,q,c->', grad_error**2, weights, measures))
    metrics = dict(n=n, dofs=space.number_of_global_dofs(), L2=float(l2),
                   H1_semi=float(h1), residual=float(residual), boundary_error=float(boundary_error))
    return mesh, uh, exact, metrics

if __name__ == '__main__':
    mesh, uh, exact, metrics = 求解(n=8)
    print(metrics)
```

## 第七步：换成实际问题的数据

实际求解只需区域、系数、源项和边界数据，不要求内部真解。没有解析解时，将 `gd=exact` 换成独立的边界函数，并移去依赖 `exact`、`gradient` 的误差计算；仍保留内部残差和边界核对。空间误差可以通过网格比较或可信参考解评估。

改变扩散或反应系数后，制造源项也要相应重算。例如扩散随位置变化时，必须从 $-
abla\cdot(a
abla u_\star)$ 重新推导，不能沿用常系数源项。对于当前物理题，改系数而保持源项不变也是允许的，但此时解已经改变，旧真解不再用于误差比较。

FEALPy负责网格、空间、积分器、装配和边界工具；模型数据和验证目的由使用者定义。已有CDR项目中的模型类与求解器是建立在库之上的项目接口，不能将其自定义方法直接当成FEALPy原生API。阅读扩展程序时，可沿“模型数据→空间→装配→边界→求解→误差”追踪同一数学对象，而不必从头学习另一套方程。

通用维数、次数与模型接口的完整实现见[收束求解器](../tests/收束求解器/README.md)，对应验证见[运行收束验证.py](../tests/运行收束验证.py)。它们与本篇二维P1主线共享装配流程，扩展能力以该实现及其验证为准。
