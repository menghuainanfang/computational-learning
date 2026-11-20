# 课程收束｜稳态多分量CDR求解器

读完专题正文并运行配套实验后，阅读这套完整源码。模型为共享扩散a与给定速度b、矩阵C线性耦合、全Dirichlet的稳态系统。阅读顺序为 [模型.py](模型.py) → [求解器.py](求解器.py) → [验证入口](../验证耦合系统.py)。

## 原项目与本课程的关系

源码依据 [原CDR模型和块求解器](https://github.com/menghuainanfang/research-workspace/tree/main/学习练手项目/cdr-solver/01_程序/cdr)，保留CDRModel、SymbolicCDR、box_mesh、CDRSolver、TensorFunctionSpace、BlockForm与求解回调的组织。采用稳态标准Galerkin，数据只含空间坐标；原项目含时、SUPG与导出材料保持源位置。

## 运行及使用

从课程根目录执行 `python tests/验证耦合系统.py`。下面片段可在tests目录中的脚本运行：

```python
from fealpy.backend import backend_manager as bm
from 收束求解器 import SymbolicCDR, CDRSolver, box_mesh
bm.set_backend("numpy")
model = SymbolicCDR(2, components=2, source=[1, 0], boundary=[0, 0],
                    reaction=[[2, -1], [-1, 2]])
solver = CDRSolver(model, box_mesh(model, 16), p=1)
uh = solver.solve_steady()
print(uh[:].reshape(2, -1))
print(solver.diagnostics)
```

需要真解验证时使用 `exact=["sin(pi*x)*sin(pi*y)", "2*sin(pi*x)*sin(pi*y)"]`，程序生成匹配的源项和边界，随后可调用 `solver.errors()`。无真解输入不能调用这一误差入口。

`assemble()`返回未施加边界的A、F；`solve_steady()`保存原系统和修改后系统，并记录内部残差与边界误差。`component_field`提取给标量积分器使用的分量。导入包不切换后端，运行入口显式设置NumPy。

## 收束实践

先解释一个非对角积分项如何进入 `forms[r][s]`，再依次改变空间维数、分量数及次数。用不同的非零边界、非对称反应矩阵重新制造源项，记录每个分量的误差。最后将反应矩阵改为对角，说明为什么退化为分别求解标量问题。

继承CDRModel编写显式数组模型时，需提供diffusion_coef、convection_coef、reaction_coef、source、dirichlet；误差比较再提供solution、gradient。具体返回形状见第三篇note，尤其不能省略m=1的物质轴。

## 本次跨维数结果

FEALPy元数据3.4.0、NumPy 2.5.3、NumPy/CPU。下表列最后两档加密的向量整体误差阶；各分量误差及全部18组记录保存在JSON。

| dim | p | 最细n | 自由度 | L2误差 | L2阶 | H1半范数阶 |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 16 | 34 | 4.757162e-03 | 2.001 | 0.998 |
| 1 | 2 | 16 | 66 | 6.879011e-05 | 2.998 | 1.998 |
| 2 | 1 | 16 | 578 | 1.114007e-02 | 1.979 | 0.989 |
| 2 | 2 | 16 | 2178 | 1.535619e-04 | 2.991 | 1.988 |
| 3 | 1 | 8 | 1458 | 6.167121e-02 | 1.787 | 0.894 |
| 3 | 2 | 8 | 9826 | 1.625723e-03 | 3.045 | 1.910 |


## 适用范围

自动网格支持1—3维盒区域的线段、三角形、四面体，已验证P1/P2与双分量空间收敛，并增加三分量多项式检查。更多分量可按同一接口组织，但不同模型需另行验证。共享正标量扩散、给定速度、线性反应、全Dirichlet；非线性、混合边界、交叉扩散及其他后端不在本次实测范围。粗三维网格不一定进入渐近收敛区，完整数据保留各档实际阶。
