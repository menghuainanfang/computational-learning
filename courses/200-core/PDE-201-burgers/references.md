# 资料与实现来源

## 模型、空间装配与时间方法

用户提供的[《Burgers方程学习笔记》](references/Burgers原版笔记.md)覆盖模型与物质导数、P1/P2非线性装配、专题A/B及完整时间步。第一篇的方程认识、第二篇的物质导数与质量矩阵说明，以及第三篇的一维三单元冻结梯度算例和函数重建，均参考了该笔记的相关推导。课程统一采用 $D_\beta$、$\boldsymbol{B}(\boldsymbol{W})$、$\boldsymbol{G}(\boldsymbol{W})$ 与 $\boldsymbol{C}(\boldsymbol{U})$ 的记号。

团队PDE-201的材料按专题提供详细推导：

- [模型、弱形式与空间半离散](https://github.com/suanhaitech/meshflow_group/blob/main/courses/200-core/PDE-201-burgers/notes/01-模型与空间半离散.md)：物质导数、标量与向量模型、弱形式和半离散系统。
- [非线性空间装配与高阶求积](https://github.com/suanhaitech/meshflow_group/blob/main/courses/200-core/PDE-201-burgers/notes/02-非线性空间装配与高阶求积.md)：P2节点与梯度、非线性张量、B/G及同冻结函数的矩阵积分比较。
- [半离散ODE与时间推进](https://github.com/suanhaitech/meshflow_group/blob/main/courses/200-core/PDE-201-burgers/notes/03-半离散ODE与时间推进.md)：初值、Euler族、BDF2/RK与基函数重建。
- [全隐式残差与非线性求解](https://github.com/suanhaitech/meshflow_group/blob/main/courses/200-core/PDE-201-burgers/notes/04-全隐式残差与非线性求解.md)：冻结矩阵性质、Picard/Newton、三种误差及计算成本。
- [误差、能量与数值验证](https://github.com/suanhaitech/meshflow_group/blob/main/courses/200-core/PDE-201-burgers/notes/05-误差能量与数值验证.md)：固定半离散时间参考、联合加密及能量平衡。

## 程序与演示

- [团队非线性有限元实现](https://github.com/suanhaitech/meshflow_group/blob/main/courses/200-core/PDE-201-burgers/labs/非线性有限元.py)：本课程序的主要来源，涵盖一二维P1/P2、C/B/G装配及Euler族。
- [团队时间实验](https://github.com/suanhaitech/meshflow_group/blob/main/courses/200-core/PDE-201-burgers/labs/运行时间实验.py)：固定网格上的DOP853参考与黏性对照。
- [二维演示模型与验证说明](https://github.com/menghuainanfang/research-workspace/blob/main/学习练手项目/burgers-lab/README.md)：交互页面的物理帧、时间格式、节点方程与能量面板。
- [一维求解与页面生成](https://github.com/menghuainanfang/research-workspace/blob/main/学习练手项目/burgers-lab/solve.py)：守恒通量形式、精确P1单元积分及Cole–Hopf案例。
- [FEM-101教材](https://github.com/menghuainanfang/computational-learning/tree/main/courses/100-foundations/FEM-101-reaction-advection-diffusion)：参考单元、单元到全局装配及线性有限元基础。
- [FEALPy源码](https://github.com/weihuayi/fealpy)：网格、Lagrange空间、基函数梯度、求积和全局编号。

正文采用统一的非线性向量记号C(U)、冻结速度矩阵B(W)和冻结梯度矩阵G(W)。演示说明中的对流矩阵C(U)对应本课B(U)，其中R(U)对应本课G(U)，需按具体积分定义辨认。程序实际覆盖范围及实验条件见课程README与第四篇。
