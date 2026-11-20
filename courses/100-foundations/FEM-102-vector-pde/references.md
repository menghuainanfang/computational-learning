# 资料与实现来源

## 研究程序与讲义

- [数学原理讲义§8、§9、§11：反应块、强制性与向量误差](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/02_%E5%AD%A6%E4%B9%A0%E8%AE%B2%E4%B9%89/01_%E6%95%B0%E5%AD%A6%E5%8E%9F%E7%90%86%E8%AE%B2%E4%B9%89.md)
- [程序入门讲义§5—8：分量轴、系数适配与整体边界](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/02_%E5%AD%A6%E4%B9%A0%E8%AE%B2%E4%B9%89/02_%E7%A8%8B%E5%BA%8F%E5%85%A5%E9%97%A8%E8%AE%B2%E4%B9%89.md)
- [网格与函数空间：乘积空间与自由度](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/03_FEALPy%E5%AD%A6%E4%B9%A0/01_%E7%BD%91%E6%A0%BC%E4%B8%8E%E5%87%BD%E6%95%B0%E7%A9%BA%E9%97%B4.md)
- [积分器与装配：单元贡献、标量子块与BlockForm](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/03_FEALPy%E5%AD%A6%E4%B9%A0/02_%E7%A7%AF%E5%88%86%E5%99%A8%E4%B8%8E%E8%A3%85%E9%85%8D.md)
- [边界、求解与源码阅读：已知边界的右端修正](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/03_FEALPy%E5%AD%A6%E4%B9%A0/03_%E8%BE%B9%E7%95%8C%E6%B1%82%E8%A7%A3%E4%B8%8E%E6%BA%90%E7%A0%81%E9%98%85%E8%AF%BB.md)
- [原模型源码：CDRModel、SymbolicCDR及at_time](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/01_%E7%A8%8B%E5%BA%8F/cdr/pde.py)
- [原求解器源码：assemble与_solve](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/01_%E7%A8%8B%E5%BA%8F/cdr/solver.py)

## 导师课程与本课对应

- [模型与向量弱形式](https://github.com/suanhaitech/meshflow_group/blob/4a77fa266253f8cd3f93220c9f8a422572dfa9a9/courses/100-foundations/FEM-102-vector-pde/notes/01-多分量模型与弱形式.md)
- [分量排列、变系数非对角块及边界](https://github.com/suanhaitech/meshflow_group/blob/4a77fa266253f8cd3f93220c9f8a422572dfa9a9/courses/100-foundations/FEM-102-vector-pde/notes/02-向量空间与块矩阵.md)
- [数据接口、积分器与块装配](https://github.com/suanhaitech/meshflow_group/blob/4a77fa266253f8cd3f93220c9f8a422572dfa9a9/courses/100-foundations/FEM-102-vector-pde/notes/03-FEALPy块装配与验证.md)

本课专题一吸收交换反应与模型适定性；专题二将乘积空间、点值表排列、变反应系数和跨分量边界移项展开到实际公式与代码；专题三区分模型数据、项目组织与FEALPy对象，并保留向量误差理论及当前制造解结果。源项目的含时与SUPG实现由原项目维护，本课标准Galerkin程序未引入这些能力。

[FEALPy源码](https://github.com/weihuayi/fealpy)用于网格、空间、积分器与求解接口。当前实验结果仍来自随附二维双分量P1程序；变系数线段积分是分析例，不表述为本次新增数值实验。未声明未读取的教材页码。
