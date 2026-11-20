# 有限元与偏微分方程课程

课程围绕模型认识、弱形式、离散、求解和数值验证展开。每门课的专题正文、配套程序、验证数据及来源分别放在 `notes/`、`labs/`、`tests/` 和 `references.md` 中。

| 课程 | 主要内容 | 专题数 |
|---|---|---:|
| [FEM-101：稳态反应–对流–扩散方程](100-foundations/FEM-101-reaction-advection-diffusion/README.md) | 局部守恒、弱问题、单元映射、全局装配与误差估计 | 4 |
| [FEM-102：多分量耦合系统](100-foundations/FEM-102-vector-pde/README.md) | 交换反应、乘积空间、块矩阵、整体边界与向量误差 | 3 |
| [PDE-201：Burgers 方程](200-core/PDE-201-burgers/README.md) | 方程认识、非线性有限元、时间推进、误差与能量 | 4 |

先学习101的标量有限元，再按需要阅读102或201。各课README说明先修知识、运行命令和实验实际覆盖范围。

本目录与[Meshflow课程目录](https://github.com/suanhaitech/meshflow_group/tree/main/courses)保持对应。图片由本仓库承载，正文使用可供团队课程引用的图片链接。实验保留源程序和小型CSV／JSON结果，运行环境与缓存不入库。

数学记号采用[Meshflow项目规范](https://github.com/suanhaitech/meshflow_group/blob/main/standards/meshflow_mathematical_notation_standard.md)。项目当前主写 `\boldsymbol{}`，与其链接的SuanhaiOS上位标准中 `\bm{}` 的差异由项目规范记录，仍待协调。
