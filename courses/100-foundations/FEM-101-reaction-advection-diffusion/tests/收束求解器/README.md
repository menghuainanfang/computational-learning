# 课程收束｜稳态标量CDR求解器

读完专题正文并运行配套实验后，使用这份完整源码，把“数据—网格—空间—装配—边界—求解—验证”连成一次独立实践。

模型为

$$-\nabla\cdot(a\nabla u)+\boldsymbol b\cdot\nabla u+cu=f,\qquad u=g\ \text{on }\partial\Omega.$$

## 与原CDR项目的关系

来源为 [原CDR模型与求解器](https://github.com/menghuainanfang/research-workspace/tree/main/学习练手项目/cdr-solver/01_程序/cdr)。保留其CDRModel/SymbolicCDR数据组织、box_mesh网格入口、CDRSolver装配/求解分工、线性求解回调和制造解思路。

这里采用标量空间和稳态接口：`assemble()` 返回A、F，`solve_steady()`求解，`errors()`计算积分误差。数据函数只接收物理坐标，不携带时间、容量、初值或物质轴。使用标准Galerkin；SUPG的数学与源项目阅读入口见误差专题。

## 代码组织与阅读顺序

| 文件 | 先看什么 | 对应课程 |
|---|---|---|
| [模型.py](模型.py) | 系数、源项、边界，及从真解生成数据 | notes 01 |
| [求解器.py](求解器.py) 的 box_mesh 和初始化 | 线段/三角形/四面体，标量拉格朗日空间 | notes 02 |
| assemble | K、B、R三个积分器与载荷积分器 | 装配专题与求解流程 |
| solve_steady | 全Dirichlet边界、稀疏求解、内部残差 | 求解流程 |
| errors | 函数值及梯度的单元求积 | 误差专题与配套实验 |
| [运行收束验证.py](../运行收束验证.py) | 同一流程改变维数与次数 | 课程综合实践 |

`CDRModel`约定盒区域和边界识别。需要自定义数值模型时，继承它并实现 diffusion_coef、convection_coef、reaction_coef、source、dirichlet；计算真解误差时另提供solution与gradient。标量函数返回点数组去掉末尾坐标轴后的形状，速度与梯度再加维数轴。

## 如何运行

在课程目录执行，环境与课程README相同，另需SymPy：

```powershell
python tests/运行收束验证.py
```

结果写入 [收束验证结果.json](../收束验证结果.json)，不生成图像。运行入口显式选择NumPy后端，导入求解器包不切换后端。

下面片段在tests目录的脚本中运行，或将tests目录加入Python模块搜索路径：

```python
from fealpy.backend import backend_manager as bm
from 收束求解器 import SymbolicCDR, CDRSolver, box_mesh

bm.set_backend("numpy")
pde = SymbolicCDR(2, exact="1+sin(pi*x)*sin(pi*y)",
                  a="1+x", b=[1, 2], reaction=2)
solver = CDRSolver(pde, box_mesh(pde, n=8), p=1)
uh = solver.solve_steady()
print(solver.diagnostics)
print(solver.errors())
```

`exact`是一个标量表达式；速度是长度dim的列表。扩散系数为标量，符号源项按整个散度算子生成，包含扩散系数的空间导数。上例为非零边界。区域默认单位盒，也可以传 `box=[xmin,xmax,ymin,ymax]`。

实际输入不必知道真解，例如：

```python
pde = SymbolicCDR(1, source=2, boundary=1)
solver = CDRSolver(pde, box_mesh(pde, n=8))
uh = solver.solve_steady()
print(solver.diagnostics)
```

这段直接给出杆内均匀热源和两端边界。不能调用真解误差来评价未知真解的问题；应采用网格加密或合适的参考证据。

## 综合实践

1. 先自己画出每个方法的输入输出关系，再运行验证入口，定位三个积分器与每个检查量。
2. 在一维复现三节点手算；换成非零边界，解释内部右端怎样改变。
3. 把二维改成三维，补全速度分量和真解表达式。说明改变的是几何、基函数与积分，哪些求解步骤没有改变。
4. 改变a、b或c时，独立推导一份源项再与符号结果核对；不要只以求解成功作为正确证据。
5. 比较P1与P2的自由度和误差。分别说明h加密与p提高在空间上做了什么，并提交自己的参数、误差表和解释。

## 实测结果

2026年11月10日，FEALPy元数据3.4.0（本机源码）、NumPy 2.5.3、NumPy/CPU。18组为三个维数×两个次数×三档网格；各组使用非零边界、a=1+x、b=(1,…,dim)、c=2的正弦制造解。下面列出最后两档之间的观测阶，完整数据保存在JSON。

| 空间维数 | p | 最细n | 自由度 | L2误差 | L2阶 | H1半范数阶 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 16 | 17 | 2.171900e-03 | 2.001 | 0.998 |
| 1 | 2 | 16 | 33 | 3.076458e-05 | 2.998 | 1.998 |
| 2 | 1 | 16 | 289 | 5.040391e-03 | 1.978 | 0.989 |
| 2 | 2 | 16 | 1089 | 6.868538e-05 | 2.992 | 1.988 |
| 3 | 1 | 8 | 729 | 2.780200e-02 | 1.785 | 0.894 |
| 3 | 2 | 8 | 4913 | 7.278394e-04 | 3.047 | 1.910 |


此外，通过无真解输入与制造解输入的离散解一致性核对，以及一维变系数源项的独立手算核对。每组边界误差小于1e-11、内部相对残差小于1e-9。最后两档L2阶需大于p+0.65、梯度误差阶需大于p−0.35；这是所选光滑算例的检查范围，不要求每个粗网格区间都达到渐近阶。

## 适用范围

自动网格覆盖1—3维盒区域的全维单纯形网格；并非四维以上网格生成器。空间次数交给FEALPy的LagrangeFESpace，当前实际验证P1/P2。采用正标量扩散、标量线性反应、给定速度与全Dirichlet边界；强对流稳定化、混合边界、高维几何及其他后端不在本次实测范围。三维较粗网格的观测阶尚未达到渐近值，完整数据保留这一现象。
