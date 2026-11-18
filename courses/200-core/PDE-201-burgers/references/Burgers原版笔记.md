# Burgers 方程学习笔记

# 一、Burgers 方程初印象

## 1\.1 一维黏性 Burgers 方程

我们先从最熟悉的一维形式开始：

$\displaystyle \boxed{u_t+uu_x-\nu u_{xx}=f,\qquad \nu>0.}$

其中：

- $u=u(x,t)$：未知标量函数；在 Burgers 的经典解释中，可以把它看成一维速度。

- $u_t$：固定空间位置的时间变化率。

- $uu_x$：非线性自对流项。

- $-\nu u_{xx}$：写在左侧的扩散项；若移到右侧就是 $+\nu u_{xx}$。

- $f=f(x,t)$：外部源项。

- $\nu$：运动黏度或扩散系数，若空间与时间具有常规物理量纲，则其量纲为 $\mathrm L^2/\mathrm T$。

判断是否线性的关键在于**未知函数进入方程的方式**。例如：

$\displaystyle \nu u_{xx}-u_t-c(x,t)u_x=-f$

在给定系数 $c(x,t)$ 的情况下仍然是线性的；但是 $uu_x$ 中的乘数 $u$ 本身也是未知函数，所以是非线性的。

## 1\.2 二维标量 Burgers 方程

现在让 $u$ 依赖二维位置：

$\displaystyle u=u(x,y,t)\in\mathbb R.$

我们学习过的一种二维标量 Burgers 型模型是：

$\displaystyle \boxed{u_t+u(u_x+u_y)-\nu\Delta u=f,\qquad \Delta u=u_{xx}+u_{yy}.}$

这是一种具体的二维标量模型，并不是一维 Burgers 在高维的唯一推广；它与后面“二维向量 Burgers”也不是同一个未知量结构。

我曾经想当然地把一维中的 $uu_x$ 变成 $u\nabla u$，但这会出现类型不一致：

$\displaystyle u\nabla u=(uu_x,uu_y)^T$

是向量，而 $u_t$、$\Delta u$、$f$ 都是标量，不能直接相加。

一种正确写法是引入依赖于解的对流速度：

$\displaystyle \mathbf b=(u,u)^T.$

于是：

$\displaystyle \mathbf b\cdot\nabla u=uu_x+uu_y.$

方程写成：

$\displaystyle \boxed{u_t+\mathbf b\cdot\nabla u-\nu\Delta u=f,\qquad\mathbf b=(u,u)^T.}$

## 1\.3 二维向量 Burgers 方程

若未知函数本身就是一个二维速度向量：

$\displaystyle \mathbf u=(u,v)^T,$

则二维向量 Burgers 可写为：

$\displaystyle \boxed{\mathbf u_t+(\mathbf u\cdot\nabla)\mathbf u-\nu\Delta\mathbf u=\mathbf f.}$

这里：

$\displaystyle \mathbf u\cdot\nabla=u\partial_x+v\partial_y.$

**括号中的东西首先是微分算子**，然后它才作用于向量场 $\mathbf u$：

$\displaystyle (\mathbf u\cdot\nabla)\mathbf u= \begin{pmatrix}  uu_x+vu_y\\uv_x+vv_y \end{pmatrix}.$

因此原来的向量方程是以下两条相互耦合的方程：

$\displaystyle \begin{cases}  u_t+uu_x+vu_y-\nu(u_{xx}+u_{yy})=f_1,\\  v_t+uv_x+vv_y-\nu(v_{xx}+v_{yy})=f_2. \end{cases}$

这里有两个未知标量函数 $u,v$，而上一节的二维标量模型只有一个未知函数 $u$。这一区别会直接影响有限元的自由度数、矩阵结构及求解方法。

## 1\.4 三维

三维向量速度场为 $\mathbf u=(u,v,w)^T$，对应的对流项为：

$\displaystyle (\mathbf u\cdot\nabla)\mathbf u= \begin{pmatrix}  uu_x+vu_y+wu_z\\  uv_x+vv_y+wv_z\\  uw_x+vw_y+ww_z \end{pmatrix}.$

## 1\.5 与 Navier–Stokes 方程的关系

在密度 $\rho$、动力黏度 $\mu$ 为常数，且流体为不可压缩牛顿流体的假设下：

$\displaystyle \begin{cases} \rho\bigl[\mathbf u_t+(\mathbf u\cdot\nabla)\mathbf u\bigr] =-\nabla p+\mu\Delta\mathbf u+\rho\mathbf b,\\ \nabla\cdot\mathbf u=0. \end{cases}$

令 $\nu=\mu/\rho$，$\pi=p/\rho$，得到：

$\displaystyle \begin{cases} \mathbf u_t+(\mathbf u\cdot\nabla)\mathbf u +\nabla\pi-\nu\Delta\mathbf u=\mathbf b,\\ \nabla\cdot\mathbf u=0. \end{cases}$

与向量 Burgers 比较，两者有相同形式的非线性自对流项和黏性扩散项。但是：

1. NS 还包含压力及其梯度，压力和速度通常需要耦合求解。

2. 不可压缩 NS 还要求 $\nabla\cdot\mathbf u=0$。

3. 一维 Burgers **不是**“一维不可压缩 NS 直接去掉压力”：一维不可压缩条件会给出 $u_x=0$，而 Burgers 可以有 $u_x\ne0$。

# 二、方程的推导

## 2\.1 两种观察流体的方式

经典力学中的牛顿第二定律是 $\mathbf F=m\mathbf a$，而质点的加速度是速度对时间的导数。

但流体中的速度是一个场：

$\displaystyle \mathbf u=\mathbf u(\mathbf x,t).$

我们需要区别两种观察方式：

- **欧拉描述**：固定空间位置 $\mathbf x$，观察那里的速度如何随时间改变。

- **拉格朗日描述**：跟随同一个流体质点，观察它在运动过程中经历的速度变化。

这不是说一定要更换参考系，而是观察位置是否随质点运动。

## 2\.2 从加速度差商开始

设某个流体质点在时刻 $t$ 的位置为 $\mathbf X(t)=\mathbf x$，经过 $\Delta t$ 后位移为 $\Delta\mathbf x$。

那么它的加速度为：

$\displaystyle \mathbf a= \lim_{\Delta t\to0} \frac{\mathbf u(\mathbf x+\Delta\mathbf x,t+\Delta t)-\mathbf u(\mathbf x,t)}{\Delta t}.$

这里比较的是**同一个质点**在两个时刻的速度。按照此前看视频时使用的技巧，在分子中加减 $\mathbf u(\mathbf x,t+\Delta t)$：

$\displaystyle \begin{aligned} \mathbf a=\lim_{\Delta t\to0}\Bigg[ &\frac{\mathbf u(\mathbf x,t+\Delta t)-\mathbf u(\mathbf x,t)}{\Delta t}\\ +&\frac{\mathbf u(\mathbf x+\Delta\mathbf x,t+\Delta t)-\mathbf u(\mathbf x,t+\Delta t)}{\Delta t} \Bigg]. \end{aligned}$

这样就把时间变化和位置变化拆开了。

## 2\.3 第一部分：局部加速度

第一项的空间位置没有改变，因此：

$\displaystyle \lim_{\Delta t\to0} \frac{\mathbf u(\mathbf x,t+\Delta t)-\mathbf u(\mathbf x,t)}{\Delta t} =\boxed{\mathbf u_t}.$

它反映的是固定地点的速度变化。

## 2\.4 第二部分：对流加速度

以二维为例，$\mathbf x=(x,y)$，$\Delta\mathbf x=(\Delta x,\Delta y)$。

在固定时刻 $t+\Delta t$ 对空间变量进行一阶 Taylor 展开：

$\displaystyle \begin{aligned} &\mathbf u(x+\Delta x,y+\Delta y,t+\Delta t)\\ ={}&\mathbf u(x,y,t+\Delta t) +\mathbf u_x\,\Delta x+\mathbf u_y\,\Delta y +o(\|\Delta\mathbf x\|). \end{aligned}$

把它代回第二项并除以 $\Delta t$：

$\displaystyle \mathbf u_x\frac{\Delta x}{\Delta t} +\mathbf u_y\frac{\Delta y}{\Delta t} +\frac{o(\|\Delta\mathbf x\|)}{\Delta t}.$

由于质点随速度场运动，$d\mathbf X/dt=\mathbf u(\mathbf X(t),t)$。因此：

$\displaystyle \frac{\Delta x}{\Delta t}\longrightarrow u, \qquad \frac{\Delta y}{\Delta t}\longrightarrow v.$

若速度场足够光滑、质点速度有界，则 $\|\Delta\mathbf x\|=O(\Delta t)$，余项除以 $\Delta t$ 仍趋于零。

所以第二项的极限为：

$\displaystyle u\mathbf u_x+v\mathbf u_y.$

用微分算子写就是：

$\displaystyle \boxed{(\mathbf u\cdot\nabla)\mathbf u.}$

## 2\.5 物质导数与链式法则

把两部分合并：

$\displaystyle \boxed{ \frac{D\mathbf u}{Dt}=\mathbf u_t+(\mathbf u\cdot\nabla)\mathbf u. }$

这里的大写 $D/Dt$ 叫作**物质导数**或**随体导数**。

如果把质点位置看成函数 $\mathbf X(t)$，也可以直接使用多元复合函数链式法则：

$\displaystyle \begin{aligned} \frac{d}{dt}\mathbf u(\mathbf X(t),t) &=\mathbf u_t+(D\mathbf u)\frac{d\mathbf X}{dt}\\ &=\mathbf u_t+(D\mathbf u)\mathbf u\\ &=\mathbf u_t+(\mathbf u\cdot\nabla)\mathbf u. \end{aligned}$

# 三、从强形式到弱形式

## 3\.1 初边值问题

我们先考虑一维问题：

$\displaystyle \begin{cases}  u_t+uu_x-\nu u_{xx}=f,&0<x<1,\ 0<t\le T,\\  u(x,0)=u_0(x),&0\le x\le1,\\  u(0,t)=u(1,t)=0,&0<t\le T. \end{cases}$

我们暂时选择齐次 Dirichlet 条件，因为它让分部积分的边界项消失，便于先理解核心方法。初值也应与边界条件相容，例如 $u_0(0)=u_0(1)=0$。

## 3\.2 弱形式

1. 先乘上测试函数 $v$ 并在区域上积分；

2. 对扩散项做分部积分，把一个空间导数转移给测试函数；

3. 再使用连续的分片线性函数进行离散。

取 $v\in H_0^1(0,1)$，有：

$\displaystyle \int_0^1u_tv\,dx +\int_0^1uu_xv\,dx -\nu\int_0^1u_{xx}v\,dx =\int_0^1fv\,dx.$

分部积分：

$\displaystyle -\nu\int_0^1u_{xx}v\,dx =-\nu[u_xv]_0^1+\nu\int_0^1u_xv_x\,dx.$

因为 $v(0)=v(1)=0$，边界项为零。因此弱形式为：

$\displaystyle \boxed{ \int_0^1u_tv\,dx +\int_0^1uu_xv\,dx +\nu\int_0^1u_xv_x\,dx =\int_0^1fv\,dx, \quad\forall v\in H_0^1(0,1). }$

这里的 $H_0^1(0,1)$ 表示具有平方可积的一阶弱导数、且边界迹为零的函数空间。现阶段只需理解：这样的函数不一定具有经典二阶导数，却足以使弱形式中的一阶导数积分有意义。

## 3\.3 二维标量 Burgers 的弱形式

取二维区域 $\Omega\subset\mathbb R^2$，考虑：

$\displaystyle \begin{cases}  u_t+u(u_x+u_y)-\nu\Delta u=f,&(x,y)\in\Omega,\ 0<t\le T,\\  u(\mathbf x,0)=u_0(\mathbf x),\\  u|_{\partial\Omega}=0. \end{cases}$

测试函数取 $v\in H_0^1(\Omega)$。乘上 $v$ 并积分后，扩散项运用多维分部积分（Green 公式）：

$\displaystyle -\int_\Omega(\Delta u)v\,d\mathbf x =\int_\Omega\nabla u\cdot\nabla v\,d\mathbf x -\int_{\partial\Omega}\partial_nu\,v\,ds.$

由于 $v$ 在边界为零，得到：

$\displaystyle \boxed{ \begin{aligned} \int_\Omega u_tv\,d\mathbf x &+\int_\Omega u(u_x+u_y)v\,d\mathbf x\\ &+\nu\int_\Omega\nabla u\cdot\nabla v\,d\mathbf x =\int_\Omega fv\,d\mathbf x. \end{aligned} }$

与一维比较，变化的是：积分区域由区间变成二维区域；$u_xv_x$ 变成 $\nabla u\cdot\nabla v$；对流项则按照我们选择的具体二维标量 Burgers 模型书写。

# 四、一维 P1 有限元离散

## 4\.1 从未知函数变成未知系数

将区间划分成节点：

$\displaystyle 0=x_0<x_1<\cdots<x_N=1.$

每两个相邻节点组成一个单元 $K_e=[x_e,x_{e+1}]$，长度 $h_e=x_{e+1}-x_e$。

我们选择 **P1 一次 Lagrange 元**：每个单元上的函数为一次多项式，整个区间上连续。

在一个单元 $K=[x_L,x_R]$，定义两个局部基函数：

$\displaystyle \varphi_L(x)=\frac{x_R-x}{h},\qquad \varphi_R(x)=\frac{x-x_L}{h}.$

它们在各自的节点上为1，在另一个节点上为0。于是：

$\displaystyle u_h(x,t)=U_L(t)\varphi_L(x)+U_R(t)\varphi_R(x).$

全局上写成：

$\displaystyle \boxed{u_h(x,t)=\sum_{j=1}^mU_j(t)\varphi_j(x).}$

这里 $m$ 是施加齐次 Dirichlet 条件后保留的自由度数。固定网格时，$\varphi_j$ 不随时间改变，因此：

$\displaystyle \partial_tu_h=\sum_{j=1}^m\dot U_j\varphi_j.$

## 4\.2 质量矩阵

在弱形式中取测试函数 $v_h=\varphi_i$，时间项为：

$\displaystyle \begin{aligned} \int_0^1(\partial_tu_h)\varphi_i\,dx &=\int_0^1\left(\sum_j\dot U_j\varphi_j\right)\varphi_i\,dx\\ &=\sum_j\dot U_j\underbrace{\int_0^1\varphi_j\varphi_i\,dx}_{M_{ij}}. \end{aligned}$

所以定义质量矩阵：

$\displaystyle \boxed{M_{ij}=\int_0^1\varphi_i\varphi_j\,dx.}$

时间项就变成 $M\dot U$。

对单元 $K=[x_L,x_R]$，直接积分得到：

$\displaystyle \boxed{ M^K=\frac h6\begin{pmatrix}2&1\\1&2\end{pmatrix}. }$

其中一个具体积分是：

$\displaystyle \int_K\varphi_L^2\,dx=\frac h3,$

而交叉项是：

$\displaystyle \int_K\varphi_L\varphi_R\,dx=\frac h6.$

## 4\.3 刚度矩阵

扩散项为：

$\displaystyle \nu\int_0^1(u_h)_x(\varphi_i)_x\,dx =\nu\sum_jU_j\int_0^1\varphi_j'\varphi_i'\,dx.$

定义：

$\displaystyle \boxed{K_{ij}=\int_0^1\varphi_i'\varphi_j'\,dx.}$

于是扩散项变成 $\nu KU$。

在一个 P1 线段单元上：

$\displaystyle \varphi_L'=-\frac1h,\qquad\varphi_R'=\frac1h.$

所以：

$\displaystyle \boxed{ K^K=\frac1h\begin{pmatrix}1&-1\\-1&1\end{pmatrix}. }$

**注意**：这里的 $K$ 是刚度矩阵，不要和表示单元的 $K$ 混淆；后者通常用上标 $M^K,K^K$ 表示“这个单元的贡献”。

## 4\.4 非线性对流项

我们需要计算：

$\displaystyle C_i(U)=\int_0^1u_h(u_h)_x\varphi_i\,dx.$

代入基函数展开：

$\displaystyle \begin{aligned} C_i(U) &=\int_0^1\left(\sum_jU_j\varphi_j\right)                  \left(\sum_kU_k\varphi_k'\right)\varphi_i\,dx\\ &=\sum_{j,k}U_jU_k\int_0^1\varphi_j\varphi_k'\varphi_i\,dx. \end{aligned}$

出现了 $U_jU_k$，所以它是关于未知系数的二次非线性函数。

在局部单元内：

$\displaystyle u_h=U_L\varphi_L+U_R\varphi_R,$

$\displaystyle \boxed{(u_h)_x=\frac{U_R-U_L}{h}.}$

这个导数在当前 P1 单元内是常数，因此可以从积分中提出去。

分别计算两个局部分量，得到：

$\displaystyle \boxed{ C^K(U)=\frac{U_R-U_L}{6} \begin{pmatrix}2U_L+U_R\\U_L+2U_R\end{pmatrix}. }$

这一步是整个非线性项离散的起点：**把导数换成 P1 插值函数的常导数后，剩余部分就由基函数积分决定**。它不是在另外做一遍有限差分，而是有限元插值函数求导的结果。

## 4\.5 一维半离散系统

将时间项、扩散项、非线性项、源项合在一起：

$\displaystyle \boxed{M\dot U+\nu KU+C(U)=F(t),}$

其中 $F_i(t)=\int_0^1 f(x,t)\varphi_i(x)\,dx$。

它叫**空间半离散系统**：空间已经变成有限维矩阵，但时间仍然连续，所以本质上是一个有限维非线性常微分方程组。

# 五、二维三角形 P1 单元

## 5\.1 二维 P1 单元可以用一个平面近似

考虑真实标量函数 $u(x,y)$。若把 $z=u(x,y)$ 画出来，它一般是弯曲的曲面。

有限元选择一个有限维近似空间。在每个三角形 $K$ 上，选择：

$\displaystyle \boxed{P_1(K)=\operatorname{span}\{1,x,y\}.}$

所以单元内近似函数为：

$\displaystyle \boxed{u_h|_K=a+bx+cy.}$

三个不共线的三角形顶点，加上它们对应的三个函数值，能够唯一确定这样的一个平面。

## 5\.2 标准三角形的一个算例

取：

$\displaystyle P_0=(0,0),\qquad P_1=(1,0),\qquad P_2=(0,1).$

三个形函数为：

$\displaystyle \boxed{ \lambda_0=1-x-y,\qquad\lambda_1=x,\qquad\lambda_2=y. }$

这里使用 $\lambda_i$ 表示三角形的重心坐标形函数，与前面使用的 $\varphi_i$ 是同一种 P1 基函数，只是局部记号不同。

插值函数：

$\displaystyle u_h=U_0\lambda_0+U_1\lambda_1+U_2\lambda_2.$

代入得到：

$\displaystyle u_h=U_0+(U_1-U_0)x+(U_2-U_0)y.$

于是：

$\displaystyle \boxed{ \nabla u_h= \begin{pmatrix}U_1-U_0\\U_2-U_0\end{pmatrix}. }$

而各形函数的梯度是：

$\displaystyle \nabla\lambda_0=\begin{pmatrix}-1\\-1\end{pmatrix},\quad \nabla\lambda_1=\begin{pmatrix}1\\0\end{pmatrix},\quad \nabla\lambda_2=\begin{pmatrix}0\\1\end{pmatrix}.$

所以也可以写成：

$\displaystyle \boxed{\nabla u_h=\sum_{j=0}^2U_j\nabla\lambda_j.}$

**这是二维 P1 的核心性质**：每个 $\lambda_j$ 是一次函数，因此它的梯度是常向量；由它们线性组合而来的 $u_h$，在该单元上的梯度也是常向量。

如果改用 P2 二次元，函数可能包含 $x^2,xy,y^2$，梯度就通常不再是常向量。

## 5\.3 任意三角形算

现实中的网格不可能都是单位直角三角形。设任意非退化三角形的三个顶点为：

$\displaystyle P_0=(x_0,y_0),\quad P_1=(x_1,y_1),\quad P_2=(x_2,y_2).$

单元内仍有：

$\displaystyle u_h=a+bx+cy.$

将三个节点值代入：

$\displaystyle \begin{cases} U_0=a+bx_0+cy_0,\\ U_1=a+bx_1+cy_1,\\ U_2=a+bx_2+cy_2. \end{cases}$

用后两个方程分别减去第一个：

$\displaystyle \begin{pmatrix} x_1-x_0&y_1-y_0\\ x_2-x_0&y_2-y_0 \end{pmatrix} \begin{pmatrix}b\\c\end{pmatrix} = \begin{pmatrix}U_1-U_0\\U_2-U_0\end{pmatrix}.$

因此：

$\displaystyle \boxed{ \nabla u_h= \begin{pmatrix} x_1-x_0&y_1-y_0\\ x_2-x_0&y_2-y_0 \end{pmatrix}^{-1} \begin{pmatrix}U_1-U_0\\U_2-U_0\end{pmatrix}. }$

只要三个顶点不共线，上述矩阵可逆。

## 5\.4 参考单元与雅可比矩阵

程序中通常不会每遇到一个三角形就重新解一次 $a,b,c$，而是利用参考单元与物理单元之间的几何变换。

取参考三角形：

$\displaystyle \widehat K=\{(\xi,\eta):\xi\ge0,\eta\ge0,\xi+\eta\le1\}.$

定义映射：

$\displaystyle \begin{pmatrix}x\\y\end{pmatrix} = \begin{pmatrix}x_0\\y_0\end{pmatrix} +J_K\begin{pmatrix}\xi\\\eta\end{pmatrix},$

其中：

$\displaystyle \boxed{ J_K= \begin{pmatrix} x_1-x_0&x_2-x_0\\ y_1-y_0&y_2-y_0 \end{pmatrix}. }$

注意，$J_K$ 的**每一列**是一条从 $P_0$ 出发的边向量。它是上一节系数方程中矩阵的转置，因此前面的梯度公式等价于：

$\displaystyle \boxed{\nabla u_h=J_K^{-T} \begin{pmatrix}U_1-U_0\\U_2-U_0\end{pmatrix}.}$

对形函数也有：

$\displaystyle \boxed{ \nabla\lambda_1=J_K^{-T}\begin{pmatrix}1\\0\end{pmatrix},\qquad \nabla\lambda_2=J_K^{-T}\begin{pmatrix}0\\1\end{pmatrix},\qquad \nabla\lambda_0=-\nabla\lambda_1-\nabla\lambda_2. }$

出现逆转置的原因是**多元函数链式法则**。这一步未来在三维四面体中会原样出现。

同时，面积变换满足：

$\displaystyle \boxed{|K|=\frac12|\det J_K|.}$

## 5\.5 一个一般三角形的算例

我们不直接停在抽象公式上。取：

$\displaystyle P_0=(0,0),\quad P_1=(2,0),\quad P_2=(0,1),$

节点值为：

$\displaystyle U_0=1,\qquad U_1=5,\qquad U_2=4.$

在这个三角形上：

$\displaystyle \lambda_0=1-\frac x2-y,\qquad\lambda_1=\frac x2,\qquad\lambda_2=y.$

所以：

$\displaystyle \begin{aligned} u_h&=1\cdot\lambda_0+5\lambda_1+4\lambda_2\\ &=1+2x+3y. \end{aligned}$

直接求导：

$\displaystyle \boxed{u_{h,x}=2,\qquad u_{h,y}=3.}$

这次梯度不是简单的 $(U_1-U_0,U_2-U_0)$，因为 $x$ 方向边长是2：

$\displaystyle \frac{U_1-U_0}{2}=2,\qquad\frac{U_2-U_0}{1}=3.$

这只是一个直角边恰好与坐标轴平行的例子；一般斜三角形要用上一节的 $J_K^{-T}$。

# 六、专题A：非线性有限元矩阵与高阶装配

## 6\.1 从P1弱形式得到局部矩阵

对二维标量模型：

$\displaystyle u_t+u(u_x+u_y)-\nu\Delta u=f,$

弱形式中的三个空间对象为：

$\displaystyle M_{ij}=\int_\Omega\varphi_i\varphi_j\,d\mathbf x, \qquad K_{ij}=\int_\Omega\nabla\varphi_i\cdot\nabla\varphi_j\,d\mathbf x,$

$\displaystyle C_i(U)=\int_\Omega u_h(u_{h,x}+u_{h,y})\varphi_i\,d\mathbf x.$

二维 P1 三角形每个单元有三个局部形函数，故局部质量矩阵、刚度矩阵均为 $3\times3$，局部非线性项是三维向量。

在任意三角形 $K$ 上：

$\displaystyle M_{ij}^K=\int_K\lambda_i\lambda_j\,d\mathbf x.$

P1 三角形上有一个标准积分公式：

$\displaystyle \int_K\lambda_i^2\,d\mathbf x=\frac{|K|}{6},\qquad \int_K\lambda_i\lambda_j\,d\mathbf x=\frac{|K|}{12}\quad(i\ne j).$

因此：

$\displaystyle \boxed{ M^K=\frac{|K|}{12} \begin{pmatrix} 2&1&1\\ 1&2&1\\ 1&1&2 \end{pmatrix}. }$



因为 $\nabla\lambda_i$ 在单元内是常向量：

$\displaystyle \begin{aligned} K^K_{ij} &=\int_K\nabla\lambda_i\cdot\nabla\lambda_j\,d\mathbf x\\ &=\boxed{|K|\nabla\lambda_i\cdot\nabla\lambda_j}. \end{aligned}$



这样一来，二维刚度矩阵的计算也不需要对复杂函数做积分，只需要三个梯度向量之间的点积。

对于第5\.5节的三角形，面积为 $|K|=1$，三个形函数梯度分别为：

$\displaystyle \nabla\lambda_0=(-1/2,-1)^T,\quad \nabla\lambda_1=(1/2,0)^T,\quad \nabla\lambda_2=(0,1)^T.$

所以：

$\displaystyle \boxed{ K^K= \begin{pmatrix} 5/4&-1/4&-1\\ -1/4&1/4&0\\ -1&0&1 \end{pmatrix}. }$

可以检查矩阵的每一行和都为零：这是因为 $\lambda_0+\lambda_1+\lambda_2=1$，所以三个梯度之和为零。

现在真正计算非线性部分：

$\displaystyle C_i^K(U)=\int_Ku_h(u_{h,x}+u_{h,y})\lambda_i\,d\mathbf x.$

定义：

$\displaystyle q_K:=u_{h,x}+u_{h,y}.$

因为 $u_h$ 是 P1 函数，$q_K$ 在单元内是常数（但依赖节点值 $U$，且可随时间变化），因此：

$\displaystyle C_i^K(U)=q_K\int_Ku_h\lambda_i\,d\mathbf x.$

又因为：

$\displaystyle u_h=\sum_{j=0}^2U_j\lambda_j,$

故：

$\displaystyle \begin{aligned} C_i^K(U) &=q_K\sum_{j=0}^2U_j\int_K\lambda_i\lambda_j\,d\mathbf x\\ &=q_K\sum_{j=0}^2M_{ij}^KU_j. \end{aligned}$

于是：

$\displaystyle \boxed{C^K(U)=q_KM^KU.}$

装配后得到如下空间半离散系统，时间仍然连续：

$\displaystyle \boxed{M\dot U+\nu KU+C(U)=F(t).}$

## 6\.2 先把几个容易混淆的对象分开

前面我们已经推导出有限元半离散方程：

$\displaystyle \boxed{M\dot U+\nu KU+C(U)=F(t).}$

其中：

- $M$：质量矩阵，由两个基函数相乘积分得到。

- $K$：刚度矩阵，由两个基函数的梯度做内积后积分得到。

- $C(U)$：**非线性对流向量**，不是一个固定矩阵。

- $U(t)$：所有自由度上、仍随时间变化的未知系数。

为简化方向导数，记：

$\displaystyle \mathbf s=(1,1)^T,\qquad D_su:=\mathbf s\cdot\nabla u=u_x+u_y.$

于是对流项可以简写为 $uD_su$。

注意区分四件事情：

1. **有限元展开**：$u_h=\sum_jU_j\varphi_j$，将未知函数表示成有限个系数。

2. **计算局部积分**：在一个三角形 $K$ 内求出 $M^K,K^K,C^K$ 或对流矩阵。

3. **全局装配**：按共享节点、共享边自由度的位置，将不同单元的贡献相加。

4. **时间离散及代数求解**：决定哪个时间层的 $U$ 已知、哪个未知，并解方程。这是后两个专题的主要内容。

换句话说，**“把非线性项写成矩阵乘向量”与“把非线性方程线性化”不是同一件事**。

---

## 6\.3 回顾 P1：为什么有一个特别简单的公式？

前一章已经推导了几何映射。这里沿用P1形函数 $\lambda_j$，令 $g_j=\mathbf s\cdot\nabla\lambda_j$，则 $q_K=g^TU^K$。它在当前单元内是常数，但依赖本单元的节点值，并非全局常数。

### 从非线性积分回看P1的简化

第 $i$ 个局部对流分量为

$\displaystyle \begin{aligned}C_i^K(U^K)&=\int_K u_hD_su_h\,\lambda_i\,d\mathbf x\\&=q_K\int_K u_h\lambda_i\,d\mathbf x\\&=q_K\sum_jU_j^K\int_K\lambda_i\lambda_j\,d\mathbf x.\end{aligned}$

最后的积分恰好是质量矩阵：

$\displaystyle M_{ij}^K=\int_K\lambda_i\lambda_j\,d\mathbf x=\frac{|K|}{12}\begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix}_{ij}.$

因此

$\displaystyle \boxed{C^K(U^K)=(g^TU^K)M^KU^K.}$

这与本章前面的P1推导一致。它看上去很简洁，是因为 **P1 的梯度为常数**，并不是对流项不再非线性。

同样地，局部刚度矩阵为

$\displaystyle K_{ij}^K=|K|\,\nabla\lambda_i\cdot\nabla\lambda_j.$

---

## 6\.4 局部结果为什么能够拼成全局结果？

考虑正方形 $[0,1]^2$，顶点编号为

$\displaystyle P_0=(0,0),\quad P_1=(1,0),\quad P_2=(0,1),\quad P_3=(1,1).$

将它分成两个三角形：

$\displaystyle K_1=(0,1,2),\qquad K_2=(1,3,2).$

两个单元的 P1 局部质量矩阵都是

$\displaystyle M^{K_e}=\frac1{24}\begin{pmatrix}2&1&1\\1&2&1\\1&1&2\end{pmatrix},\qquad e=1,2.$

按照单元的局部节点到全局节点的编号，将相同位置贡献累加，得到

$\displaystyle \boxed{M=\frac1{24}\begin{pmatrix}2&1&1&0\\1&4&2&1\\1&2&4&1\\0&1&1&2\end{pmatrix}.}$

例如 $M_{11}=2/24+2/24=4/24$，是因为节点 $P_1$ 出现在两个单元中。

非线性向量也遵循**同一个装配规则**：先在每个单元代入本单元的节点值计算 $C^{K_e}(U^{K_e})$，再把局部向量按全局自由度编号累加。

为了写得紧凑，令 $P_e$ 为从全局向量中选取本单元自由度的提取矩阵，则

$\displaystyle U^{K_e}=P_eU,\qquad\boxed{C(U)=\sum_eP_e^T C^{K_e}(P_eU).}$

线性矩阵则是

$\displaystyle M=\sum_eP_e^TM^{K_e}P_e,\qquad K=\sum_eP_e^TK^{K_e}P_e.$

可以发现：**全局装配并不要求局部项是线性的。** 非线性向量照样可以先计算、再累加。

> 上面的小网格只是演示装配，尚未施加 Dirichlet 边界条件。正式求解时还要处理边界自由度。
>
>

---

## 6\.5 从 P1 进入 P2：六节点元如何构造？

### 6\.5\.1 六个自由度从哪里来？

在一个直边三角形中，P2 的局部函数空间为

$\displaystyle P_2(K)=\operatorname{span}\{1,x,y,x^2,xy,y^2\}.$

其维数是 6。标准的拉格朗日 P2 三角形通常取：三个顶点，以及三条边各自的中点，共六个插值节点。

为什么特别选边中点？不仅因为容易计算，还因为**相邻三角形在公共边上的函数都是一元二次多项式**。只要边两端及中点的三个节点值一致，这两个一元二次多项式就在整条公共边上相等，从而保证拼接后 $u_h$ 连续。

需要注意：并不是任意摆放的六个节点都能唯一确定一个二维二次多项式；标准 P2 节点的几何布置满足唯一插值条件。

### 6\.5\.2 用重心坐标写出六个 P2 形函数

沿用 $\lambda_0,\lambda_1,\lambda_2$，满足

$\displaystyle \lambda_0+\lambda_1+\lambda_2=1.$

取局部自由度顺序为三个顶点 $0,1,2$，再依次为边中点 $01,12,20$。对应的六个 P2 形函数为

$\displaystyle \boxed{\begin{aligned}\varphi_0&=\lambda_0(2\lambda_0-1),\\\varphi_1&=\lambda_1(2\lambda_1-1),\\\varphi_2&=\lambda_2(2\lambda_2-1),\\\varphi_{01}&=4\lambda_0\lambda_1,\\\varphi_{12}&=4\lambda_1\lambda_2,\\\varphi_{20}&=4\lambda_2\lambda_0.\end{aligned}}$

这些函数的共同特征是：**在自己的节点取 1，在其余五个插值节点取 0**。因此局部近似解可以写成

$\displaystyle \boxed{u_h\big|_K=U_0^K\varphi_0+U_1^K\varphi_1+U_2^K\varphi_2+U_{01}^K\varphi_{01}+U_{12}^K\varphi_{12}+U_{20}^K\varphi_{20}.}$

如果是节点型 P2 基函数，这些系数正是六个插值节点上的函数值。

### 6\.5\.3 P2 的梯度为什么不能直接提出积分？

对顶点基函数和边中点基函数分别求梯度：

$\displaystyle \boxed{\nabla\varphi_i=(4\lambda_i-1)\nabla\lambda_i,\quad i=0,1,2,}$

$\displaystyle \boxed{\nabla\varphi_{ij}=4(\lambda_i\nabla\lambda_j+\lambda_j\nabla\lambda_i).}$

由于 $\lambda_i$ 本身随位置线性变化，这些梯度通常也是**单元内的一次函数**，不再是常向量。

例如，在参考三角形 $\widehat K$ 中：

$\displaystyle \lambda_0=1-x-y,\quad\lambda_1=x,\quad\lambda_2=y,$

所以

$\displaystyle \varphi_{01}=4x(1-x-y),\qquad\nabla\varphi_{01}=\begin{pmatrix}4-8x-4y\\-4x\end{pmatrix}.$

这次不能把整个 $D_su_h$ 当成一个常数提出积分号。**但它仍是已知次数的多项式，只要计算积分就可以了。**

---

## 6\.6 为什么非线性项会产生三阶张量？

这里重新从最普通的有限元展开出发。令一个单元有 $m$ 个自由度，对 P2 而言 $m=6$，则

$\displaystyle u_h=\sum_{j=1}^mU_j\varphi_j,\qquad D_su_h=\sum_{k=1}^mU_kD_s\varphi_k.$

对于测试函数 $\varphi_i$：

$\displaystyle \begin{aligned} C_i^K(U) &=\int_Ku_hD_su_h\,\varphi_i\,d\mathbf x\\ &=\int_K\left(\sum_jU_j\varphi_j\right) \left(\sum_kU_kD_s\varphi_k\right)\varphi_i\,d\mathbf x\\ &=\sum_{j,k}U_jU_k \underbrace{\int_K\varphi_i\varphi_jD_s\varphi_k\,d\mathbf x}_{T_{ijk}^K}.\end{aligned}$

定义

$\displaystyle \boxed{T_{ijk}^K:=\int_K\varphi_i\varphi_jD_s\varphi_k\,d\mathbf x,\qquad C_i^K(U)=\sum_{j,k}T_{ijk}^KU_jU_k.}$

**三个下标的来源：**

- $i$：当前测试函数（也就是结果向量的第几个分量）。

- $j$：前面的 $u_h$ 提供的基函数与系数。

- $k$：后面的 $D_su_h$ 提供的导数基函数与系数。

所以三阶张量并不神秘：**它就是一个需要三个下标才能查找的积分系数表**。P2 中的形状是 $6\times6\times6$，有 216 个数组位置，但不意味着要人工计算 216 次。

这个公式对于任意次数的拉格朗日多项式有限元都成立。二维三角形 P$p$ 元的局部自由度数为

$\displaystyle \boxed{m=\dim P_p(K)=\frac{(p+1)(p+2)}2.}$

例如 P1 为 3、P2 为 6、P3 为 10。对于直边仿射三角形及常方向 $\mathbf s$，被积函数 $\varphi_i\varphi_jD_s\varphi_k$ 的最高次数为

$\displaystyle \boxed{p+p+(p-1)=3p-1.}$

因此 P1 对应最高 2 次，P2 为 5 次，P3 为 8 次。若希望这些多项式积分精确，求积规则的代数精度至少应覆盖相应次数。

---

## 6\.7 为什么同一个非线性项能写成两种对流矩阵？

三阶张量给出

$\displaystyle C_i^K(U)=\sum_{j,k}T_{ijk}^KU_jU_k.$

这里有两个 $U$，我们可以先把其中一个与 $T$ 求和，留下另一个构成矩阵与向量的乘积。

### 6\.7\.1 把前面的函数 $u$ 吸收到矩阵里：$B(W)$

给定已知系数 $W$，对应 $w_h=\sum_jW_j\varphi_j$。定义

$\displaystyle \boxed{B_{ik}^K(W)=\sum_jT_{ijk}^KW_j=\int_Kw_h\varphi_iD_s\varphi_k\,d\mathbf x.}$

如果未知向量是 $U$，则

$\displaystyle \boxed{B^K(W)U\quad\longleftrightarrow\quad w_hD_su_h.}$

注意这里**冻结了前面的函数 **$u$，但新解的导数仍然未知，所以右侧对 $U$ 是线性的。

### 6\.7\.2 把后面的导数 $D_su$ 吸收到矩阵里：$G(W)$

这次定义

$\displaystyle \boxed{G_{ij}^K(W)=\sum_kT_{ijk}^KW_k=\int_K(D_sw_h)\varphi_i\varphi_j\,d\mathbf x.}$

于是

$\displaystyle \boxed{G^K(W)U\quad\longleftrightarrow\quad u_hD_sw_h.}$

这次**冻结的是后面的导数**。

### 6\.7\.3 为什么两种写法并不矛盾？

若两个因子使用的是**同一个向量** $U$，

$\displaystyle \boxed{B^K(U)U=G^K(U)U=C^K(U).}$

因为它们只是同一个乘积中的两个求和次序。

但一旦把矩阵所依赖的系数固定为另一个已知向量 $W$，一般就有

$\displaystyle \boxed{B^K(W)U\ne G^K(W)U.}$

这才是两种不同的线性化方式。**仅写成 **$C(U)=B(U)U$** 并没有消除非线性**，因为 $B(U)$ 自身还依赖未知量；真正形成线性系统要先固定 $W$。

|处理方式|连续表达|矩阵表达|冻结的对象|
|---|---|---|---|
|前因子冻结|$w_hD_su_h$|$B(W)U$|$u$|
|后因子冻结|$u_hD_sw_h$|$G(W)U$|$D_su$|

第7\.3节采用的 $B_{\mathrm{grad}}(W)$ 就是这里的 $G(W)$。这里用 $B(W)$ 表示冻结前因子，用 $G(W)$ 表示冻结导数，矩阵行对应测试函数，列对应待求系数。

---

## 6\.8 亲自计算一个 P2 矩阵元素：从节点值到积分结果

下面用一个可以手算到底的例子，把节点值、导数和积分连起来。我们不计算完整 $6\times6$ 矩阵，只计算其中**一个元素**，但每一步都走完。

### 第一步：选参考三角形和六节点函数

取

$\displaystyle K=\{(x,y):x\ge0,\ y\ge0,\ x+y\le1\}.$

选边 $12$ 中点的 P2 形函数

$\displaystyle \boxed{\varphi_{12}=4xy.}$

假设已知函数是

$\displaystyle \boxed{w_h(x,y)=x^2.}$

这确实属于 P2 空间。按“顶点 $0,1,2$；边中点 $01,12,20$”的顺序，它的六个节点值为

$\displaystyle \boxed{W=(0,\ 1,\ 0,\ 1/4,\ 1/4,\ 0)^T.}$

也就是说，这个已知函数可以完全由 P2 节点系数表示，并非随意假设一个与有限元无关的函数。

### 第二步：计算冻结梯度对应的矩阵 $G(W)$

先求已知函数的导数：

$\displaystyle D_sw_h=\partial_x(x^2)+\partial_y(x^2)=2x.$

注意：$2x$ **不是常数**，但因为 $w_h$ 已经给定，所以它是**已知的空间函数**。

我们计算 $G(W)$ 中行、列都对应 $\varphi_{12}$ 的那个元素：

$\displaystyle \begin{aligned}G_{12,12}^K(W)&=\int_K(D_sw_h)\varphi_{12}^2\,dx\,dy\\&=\int_K(2x)(4xy)^2\,dx\,dy\\&=32\int_Kx^3y^2\,dx\,dy.\end{aligned}$

对于标准三角形，单项式积分满足

$\displaystyle \boxed{\int_Kx^ay^b\,dx\,dy=\frac{a!b!}{(a+b+2)!},\qquad a,b\in\mathbb N_0.}$

所以

$\displaystyle \begin{aligned}G_{12,12}^K(W)&=32\frac{3!\,2!}{7!}=32\cdot\frac1{420}\\&=\boxed{\frac8{105}\approx0.07619047619}.\end{aligned}$

这已经是真正算出来的一个 **P2 局部对流矩阵元素**。

### 第三步：同一个 $w_h$，换成冻结前因子的矩阵 $B(W)$

现在再计算同一行列位置的 $B$：

$\displaystyle \begin{aligned}B_{12,12}^K(W)&=\int_Kw_h\varphi_{12}D_s\varphi_{12}\,dx\,dy\\&=\int_Kx^2(4xy)\,4(x+y)\,dx\,dy\\&=16\int_K(x^4y+x^3y^2)\,dx\,dy\\&=\boxed{\frac4{35}\approx0.11428571429}.\end{aligned}$

结果是

$\displaystyle \boxed{G_{12,12}^K(W)=\frac8{105}\ne\frac4{35}=B_{12,12}^K(W).}$

注意这里比较的是**两个矩阵的一个元素**，并不与 $B(W)W=G(W)W=C(W)$ 矛盾；整个矩阵与同一个 $W$ 相乘后，两种表示仍然给出相同的非线性向量。

这个例子说明：即便 P2 梯度不是常数，已知 $W$ 后，$B(W)$、$G(W)$ 的每个元素仍然都是可以计算的普通积分。

---

## 6\.9 实际计算中为什么使用数值求积？

刚才的 $x^3y^2$ 很容易直接积分，但一个 P2 单元有 36 个矩阵位置，整个网格还有许多单元。实际计算通常将积分转化为若干求积点的加权求和。

### 6\.9\.1 参考单元与物理单元的求积公式

在参考三角形 $\widehat K$ 上，选求积点 $\widehat{\mathbf x}_q$ 和权重 $\widehat\omega_q$，则

$\displaystyle \int_Kh(\mathbf x)\,d\mathbf x\approx\sum_{q=1}^Q\widehat\omega_q|\det J_K|\,h(\mathbf x_q),\qquad\mathbf x_q=F_K(\widehat{\mathbf x}_q).$

于是冻结梯度矩阵为

$\displaystyle \boxed{G_{ij}^K(W)\approx\sum_{q=1}^Q\widehat\omega_q|\det J_K|\,[D_sw_h](\mathbf x_q)\varphi_i(\mathbf x_q)\varphi_j(\mathbf x_q).}$

冻结函数矩阵则为

$\displaystyle \boxed{B_{ik}^K(W)\approx\sum_{q=1}^Q\widehat\omega_q|\det J_K|\,w_h(\mathbf x_q)\varphi_i(\mathbf x_q)[D_s\varphi_k](\mathbf x_q).}$

如果已将物理单元缩放吸收到积分权重 $\omega_q^K$ 中，就不能再额外乘一遍 $|\det J_K|$。

这些数值计算只涉及**已知的节点系数、形函数及其梯度、积分点和权重**。

### 6\.9\.2 使用太低阶的积分规则会怎样？

我们还用上节的例子：

$\displaystyle G_{12,12}^K(W)=32\int_Kx^3y^2\,dx\,dy=\frac8{105}.$

先试一个常见的三点规则：

$\displaystyle (x_q,y_q)=\left(\frac16,\frac16\right),\quad\left(\frac23,\frac16\right),\quad\left(\frac16,\frac23\right),$

每个参考三角形求积权重为 $1/6$。代入得

$\displaystyle \begin{aligned}G_{12,12}^{K,\,3\mathrm{pt}}&=\frac16\sum_{q=1}^3 32x_q^3y_q^2\\&=\boxed{\frac1{18}\approx0.05555555556}.\end{aligned}$

与精确值 $8/105\approx0.07619$ 不同。不是有限元本身算不出来，而是**三点规则仅能精确积分到二次多项式，无法精确计算这个五次被积函数**。

### 6\.9\.3 选用五次精度的七点规则

在参考三角形中，可以使用标准的七点五次精度求积规则。下表将对称的三个点归为一组；权重已经是适用于面积 $1/2$ 的参考三角形的权重。

|点组|$(x,y)$|每点权重|点数|
|---|---|---|---|
|中心|$(1/3,1/3)$|$0.112500000000000$|1|
|第 1 组|$(a,a),(a,b),(b,a)$，$a=0.470142064105115,\ b=0.059715871789770$|$0.066197076394253$|3|
|第 2 组|$(c,c),(c,d),(d,c)$，$c=0.101286507323456,\ d=0.797426985353087$|$0.062969590272414$|3|

将 $h(x,y)=32x^3y^2$ 代入，三组分别贡献约

$\displaystyle 0.01481481481,\qquad0.04954044803,\qquad0.01183521334.$

相加得到

$\displaystyle \boxed{G_{12,12}^{K,\,7\mathrm{pt}}\approx0.07619047619=\frac8{105}\ \text{（到舍入误差）}.}$

因此，对于直边仿射三角形的 P2 Burgers 对流项，使用能精确积分五次多项式的求积规则，便可以把这些局部矩阵积分精确计算出来（除数值舍入）。

**为什么说 P2 最高五次？** 因为对流积分含两个 P2 基函数（各最高二次）和一个 P2 基函数的方向导数（最高一次），次数总和 $2+2+1=5$。

> 本节使用了“恒定对流方向 $\mathbf s=(1,1)$、直边仿射单元、无额外非多项式系数”这些条件。如果以后变成曲边单元或更复杂的系数，积分未必仍是五次多项式，不能机械套用同一个精确度结论。
>
>

---

## 6\.10 P2 局部矩阵如何装配到整个网格？

现在我们已经知道如何在一个三角形上计算

$\displaystyle M^K,\quad K^K,\quad C^K(U^K),\quad B^K(W^K),\quad G^K(W^K).$

P2 的局部矩阵为 $6\times6$，局部对流向量为 6 维。全局装配仍然只看**每个局部自由度在整个网格中的全局编号**。

局部到全局的对应必须同时处理：

- 顶点自由度：相邻单元共享同一个顶点编号。

- 边中点自由度：相邻单元共享同一条边，必须共享该边中点对应的全局自由度。

- 单元内部自由度：标准 P2 拉格朗日三角形没有内部节点；P3 开始会出现三角形内部自由度。

举例：把一个正方形划成两个三角形，有 4 个顶点、5 条不同的边（4 条外边和 1 条内部对角线），因此标准 P2 元共有

$\displaystyle \boxed{4+5=9\text{ 个全局自由度},}$

**不是** $2\times6=12$，因为两个单元共享顶点与公共边中点。

在齐次 Dirichlet 边界条件下，外边界顶点及外边界边中点需要施加指定边界值。上面这个极小示例会只剩下内部对角线中点一个自由度；它适合说明装配，不适合作为精度研究的正式网格。

如果用 $P_e$ 提取第 $e$ 个单元的 6 个全局自由度，仍然有

$\displaystyle \boxed{\begin{aligned}C(U)&=\sum_eP_e^T C^{K_e}(P_eU),\\B(W)&=\sum_eP_e^TB^{K_e}(P_eW)P_e,\\G(W)&=\sum_eP_e^TG^{K_e}(P_eW)P_e.\end{aligned}}$

**计算局部积分的方式随 P1、P2、P3 变化，但“按照全局自由度编号累加”的装配规则不变。**

---

## 6\.11 从空间装配走向时间推进

把本章内容连起来，可以得到这样一条计算顺序：

$\displaystyle \boxed{\begin{gathered}\text{网格与有限元次数}\ \longrightarrow\ \text{形函数 }\varphi_j\\\longrightarrow\ u_h=\sum_jU_j\varphi_j\\\longrightarrow\ C_i^K(U)=\int_Ku_hD_su_h\,\varphi_i\\\longrightarrow\ T_{ijk}^K=\int_K\varphi_i\varphi_jD_s\varphi_k\\\longrightarrow\ C^K(U),\ B^K(W),\ G^K(W)\\\longrightarrow\ \text{按全局自由度装配}.\end{gathered}}$

回顾四个最重要的结论：

**第一，P1 和 P2 不存在原理上的断裂。** P1 的常梯度使 $C^K(U)=(g^TU)M^KU$ 成立；P2 不能直接提出梯度，但积分仍可通过多项式积分或数值求积计算。

**第二，三阶张量只是系数表。** 它的三个下标分别来自测试函数、前面的未知函数和后面的未知导数；不必把它视为新的物理量。

**第三，两种矩阵来自同一个乘积。** 冻结前面的 $u$ 使用 $B(W)$；冻结后面的 $D_su$ 使用 $G(W)$。当 $W=U$ 时 $B(U)U=G(U)U=C(U)$，当 $W$ 与 $U$ 不同时两种线性形式一般不同。

**第四，线性化和时间离散要分开。** 只有当 $W$ 已知时，$B(W)U$ 或 $G(W)U$ 才关于未知 $U$ 线性。选 $W=U^n$ 只是常见做法之一；它与同一时间步内部用 Picard 迭代更新 $W$ 的意义不同。

### 回到半离散系统

本章解决的是“**给定一个有限元函数，我怎样计算空间离散中的各项？**”

第七章接下来要解决的是“**空间离散后得到 **$M\dot U+\nu KU+C(U)=F(t)$**，这个有限维非线性常微分方程组如何从 **$t_n$** 推进到 **$t_{n+1}$**？**”

理解空间装配后，再比较不同时间层的取法、黏性系数的影响，以及稳定性与收敛性，才有明确的计算对象。

### 可以自检的三个问题

1. 为什么 P2 的 $G_{ij}(W)$ 虽然含位置变化的梯度，但在给定 $W$ 后仍然可以计算成一个具体数字？

2. 为什么 $B(U)U=G(U)U$，却不能断言 $B(W)=G(W)$？

3. 如果两个相邻 P2 三角形共用一条边，那么这条边的中点为什么必须使用同一个全局自由度编号？

能解释这三点后，就可以带着明确的空间矩阵进入第七章。

# 七、时间离散与参考文档的求解格式

## 7\.1 模型、时间网格与格式依据

本章对应[《基于 FEALPy 调和映射移动网格》编程数学文档的“算例一：Burgers 方程”](https://wnesm678i4.feishu.cn/wiki/VLg4wtFGpiN2GzkqsYOcwL3wnHc#WAqWdeoPqoE3tqx4NiDcCJ9Vnu5)。该节给出二维标量模型和两种 Euler 显隐组合格式；后面的求解代码使用第一种，测试部分也主要展示第一种。

参考文档取无源方程，黏性系数为 $\nu=1/Re$：

$\displaystyle u_t+u(u_x+u_y)-\nu\Delta u=0.$

令 $t_n=n\Delta t$，$U^n\approx U(t_n)$。空间离散后，$M\dot U+\nu KU+C(U)=F(t)$ 仍包含连续时间导数，需要进一步离散。以下保留 $F^{n+1}$ 以支持第八章的制造解验证；复现参考无源模型时取 $F^{n+1}=0$。

## 7\.2 主格式：对流显式、扩散隐式

参考文档的第一种格式是：

$\displaystyle \frac{u^{n+1}-u^n}{\Delta t}+u^n(u_x^n+u_y^n)-\nu\Delta u^{n+1}=f^{n+1}.$

时间导数采用后向差分，扩散放在新时间层，对流全部由旧时间层计算。这是隐显混合（IMEX）Euler 格式。对每个测试函数积分后得到：

$\displaystyle M\frac{U^{n+1}-U^n}{\Delta t}+\nu KU^{n+1}+C(U^n)=F^{n+1}.$

乘以 $\Delta t$ 并移项：

$\displaystyle \boxed{(M+\Delta t\nu K)U^{n+1}=MU^n-\Delta t C(U^n)+\Delta t F^{n+1}.}$

于是每个时间步只需解一次线性系统

当前时间步中，右端的非线性向量 $C(U^n)$ 已由旧解确定，不需要在这个时间步内进行非线性迭代。

## 7\.3 对照格式：冻结旧时间层梯度

参考文档的第二种格式为：

$\displaystyle \frac{u^{n+1}-u^n}{\Delta t}+u^{n+1}(u_x^n+u_y^n)-\nu\Delta u^{n+1}=f^{n+1}.$

旧解的梯度已经确定，因此该式仍然关于 $u^{n+1}$ 线性。定义：

$\displaystyle B_{\mathrm{grad},ij}(U^n)=\int_\Omega(u_{h,x}^n+u_{h,y}^n)\varphi_j\varphi_i\,d\mathbf x.$

P1 单元内 $q_K^n$ 为常数，所以：

$\displaystyle \boxed{B_{\mathrm{grad}}^K(U^n)=q_K^nM^K.}$

装配后得到：

$\displaystyle \boxed{\bigl(M+\Delta t\nu K+\Delta t B_{\mathrm{grad}}(U^n)\bigr)U^{n+1}=MU^n+\Delta t F^{n+1}.}$

参考文档将 $\Delta t B_{\mathrm{grad}}$ 记为 $M'$，将刚度矩阵记为 $A$。这只是符号区别。

# 八、专题B：从半离散系统到完整数值求解

## 8\.1 先确认空间离散交给了我们什么

我们研究的连续方程是

$\displaystyle \boxed{u_t+u(\boldsymbol\beta\cdot\nabla u)-\nu\Delta u=f, \qquad \boldsymbol\beta=(1,1)^T.}$

在区域 $\Omega\subset\mathbb R^2$、时间区间 $[0,T]$ 上，还需要给定初始条件 $u(\mathbf x,0)=u_0(\mathbf x)$ 和边界条件。只有 PDE 的形式而没有初始、边界条件，通常并不能确定一个具体的演化问题。

第六章的专题A已经说明：选定空间有限元空间 $V_h$，将 $u$ 近似为

$\displaystyle \boxed{u_h(\mathbf x,t)=\sum_{j=1}^{N}U_j(t)\varphi_j(\mathbf x).}$

这里 $\varphi_j(\mathbf x)$ 由**固定网格**确定；$U_j(t)$ 是**随时间变化的系数**。将其代入弱形式，得到

$\displaystyle \boxed{M\dot U(t)+\nu KU(t)+C(U(t))=F(t).}\quad(B.0)$

质量矩阵、刚度矩阵、非线性对流向量与载荷向量分别为

$\displaystyle \begin{aligned} M_{ij}&=\int_\Omega\varphi_i\varphi_j\,d\mathbf x,\\ K_{ij}&=\int_\Omega\nabla\varphi_i\cdot\nabla\varphi_j\,d\mathbf x,\\ C_i(U)&=\int_\Omega u_h(\boldsymbol\beta\cdot\nabla u_h)\varphi_i\,d\mathbf x,\\ F_i(t)&=\int_\Omega f(\mathbf x,t)\varphi_i\,d\mathbf x. \end{aligned}$

注意 $(\mathrm{B}.0)$ **还不是**一个固定的 $AU=b$，因为 $\dot U=dU/dt$ 仍是时间导数。

---

## 8\.2 半离散系统为什么是一个 ODE 组？

从 PDE 到 $(\mathrm{B}.0)$，我们只离散了**空间**：无限多个空间位置处的未知函数变成 $N$ 个有限元系数。但时间 $t$ 仍连续。

因此 $(\mathrm{B}.0)$ 是一个**有限维、非线性常微分方程组**，又称空间半离散方法或 *method of lines* 的结果。在 $M$ 可逆的情况下，可形式地写为

$\displaystyle \boxed{\dot U=\mathcal F_h(t,U) :=M^{-1}\bigl[F(t)-\nu KU-C(U)\bigr].}\quad(B.1)$

这和常微分方程课程里的 $y'=g(t,y)$ 是同一类型，只是未知量从一个数变成向量。

**为什么实际计算中不建议显式计算 **$M^{-1}$**？** 因为形成逆矩阵通常没有必要，也可能破坏稀疏性。若某个计算步骤需要 $M^{-1}r$，我们直接解 $Mz=r$ 即可。

### 有限元初始向量 $U^0$ 怎么确定？

从初始条件 $u_0$ 得到 $U^0$，常见方式有：

1. **节点插值**：对于通常的 Lagrange 有限元，将自由节点系数取为 $u_0$ 在相应节点处的值（前提是 $u_0$ 的点值有意义）。

2. $L^2$** 投影**：求 $U^0$ 使 $\int_\Omega(u_h^0-u_0)v_h\,d\mathbf x=0$ 对任意 $v_h\in V_h$ 成立，即求解质量矩阵方程

$\displaystyle \boxed{MU^0=b^0,\qquad b_i^0=\int_\Omega u_0\varphi_i\,d\mathbf x.}$

以上按齐次边界且已消去受限自由度理解；一般边界需要一致处理。**算法不能从凭空出现的 **$U^0$** 开始，它来自连续问题的初始数据。**

## 8\.3 时间差商如何把 ODE 变成代数方程？

令 $t_n=n\Delta t$，记 $U^n\approx U(t_n)$。目标是由已知 $U^n$ 求 $U^{n+1}$。

先看一元例子：

$\displaystyle y'=-y^2,\qquad y(0)=1,\qquad\Delta t=0.1.$

### 8\.3\.1 显式 Euler：在旧时间层计算右端

$\displaystyle \frac{y^{n+1}-y^n}{\Delta t}=-(y^n)^2 \quad\Longrightarrow\quad \boxed{y^{n+1}=y^n-\Delta t(y^n)^2.}$

因此 $y^1=0.9$。右端全已知，无需求解非线性方程。

### 8\.3\.2 后向 Euler：在新时间层计算右端

$\displaystyle \frac{y^{n+1}-y^n}{\Delta t}=-(y^{n+1})^2 \quad\Longrightarrow\quad \boxed{y^{n+1}+\Delta t(y^{n+1})^2=y^n.}$

代入给定数据后，选择靠近初值的正根，得到 $y^1\approx0.9160797831$。此时必须求解一个**非线性代数方程**。

### 8\.3\.3 冻结一个因子：仍使用新未知量，但线性化

将 $(y^{n+1})^2$ 改为 $y^ny^{n+1}$：

$\displaystyle \frac{y^{n+1}-y^n}{\Delta t}=-y^ny^{n+1} \quad\Longrightarrow\quad \boxed{y^{n+1}=\frac{y^n}{1+\Delta t\,y^n}.}$

于是 $y^1=1/1.1\approx0.9090909091$。

三种方法数值不同，因为它们求解的是**不同的时间离散方程**。冻结一个因子并不只是“求解器优化”，如果从头到尾固定旧时间层，它已经改变了时间格式本身。

## 8\.4 Burgers 中的一步时间推进怎样写？

对 $(\mathrm{B}.0)$ 采用后向 Euler 时间差商，且黏性扩散使用新时间层。定义

$\displaystyle \boxed{A_\nu=\frac{M}{\Delta t}+\nu K, \qquad r^n=\frac{M}{\Delta t}U^n+F^{n+1}.}\quad(B.2)$

于是当前时间步的基本形式为

$\displaystyle \boxed{A_\nu U^{n+1}+\text{对流项的时间离散}=r^n.}\quad(B.3)$

这里最重要的是：**时间导数如何近似**与**对流项用哪个时间层**是两个相关但不等同的设计选择。

第六章定义过两种矩阵。设已知 $W$ 对应函数 $w_h$，则

$\displaystyle \begin{aligned} B_{ij}(W)&=\int_\Omega w_h\varphi_i(\boldsymbol\beta\cdot\nabla\varphi_j)\,d\mathbf x &&\text{（冻结前面的函数）},\\ G_{ij}(W)&=\int_\Omega(\boldsymbol\beta\cdot\nabla w_h)\varphi_i\varphi_j\,d\mathbf x &&\text{（冻结后面的梯度）}. \end{aligned}$

并且

$\displaystyle \boxed{C(U)=B(U)U=G(U)U.}\quad(B.4)$

专题A的 $D_su$ 与这里的 $\boldsymbol\beta\cdot\nabla u$ 是同一个方向导数，$\boldsymbol\beta=(1,1)^T$ 是固定方向。本文用B表示冻结前因子，用G表示冻结导数，其含义由积分定义确定。第7\.3节的梯度矩阵就是这里的G。

将两个对流因子分别放在旧、新时间层，可以得到四种处理：

|当前时间步采用的对流表达|代数方程|性质|
|---|---|---|
|$u^n(\boldsymbol\beta\cdot\nabla u^n)$|$A_\nu U^{n+1}=r^n-C(U^n)$|对流显式、扩散隐式（IMEX）|
|$u^n(\boldsymbol\beta\cdot\nabla u^{n+1})$|$[A_\nu+B(U^n)]U^{n+1}=r^n$|关于新解线性|
|$u^{n+1}(\boldsymbol\beta\cdot\nabla u^n)$|$[A_\nu+G(U^n)]U^{n+1}=r^n$|关于新解线性|
|$u^{n+1}(\boldsymbol\beta\cdot\nabla u^{n+1})$|$A_\nu U^{n+1}+C(U^{n+1})=r^n$|非线性，需进一步求解|

这里只说明**它们怎样形成可求解的方程**；各格式的稳定性、误差阶、能量性质留给专题 C。不能认为“包含更多 $n+1$ 项”就必然更高阶。

## 8\.5 是否只有 Euler？BDF2 和 Runge–Kutta 又处在什么位置？

答案是：并不只有 Euler。BDF 与 RK 都是**时间积分方法**，不是与 Picard/Newton 相互排斥的求解器。

### 8\.5\.1 BDF2：用三个时间层近似时间导数

对等时间步长，BDF2 使用

$\displaystyle \boxed{\dot U(t_{n+1})\approx\frac{3U^{n+1}-4U^n+U^{n-1}}{2\Delta t}.}$

若其他项完全隐式，则

$\displaystyle \boxed{ M\frac{3U^{n+1}-4U^n+U^{n-1}}{2\Delta t} +\nu KU^{n+1}+C(U^{n+1})=F^{n+1}.}$

它对充分光滑的解具有**二阶时间精度的理论潜力**，但非线性项仍然是非线性的，可能仍需要 Newton/Picard。BDF2 需要 $U^{n-1},U^n$，因此第一步必须用合适的启动方法获得 $U^1$；若希望保持整体二阶，启动误差也需要匹配。

可以用旧时间层外推新时刻的对流系数，例如 $\widehat U^{n+1}=2U^n-U^{n-1}$，但这只是设计二阶线性化的一种思路；精度和稳定性还需要单独证明，不能将“冻结旧解”直接套到 BDF2 就宣称二阶。

### 8\.5\.2 Runge–Kutta：一个时间步内计算多个阶段

把半离散方程写成 $\dot U=\mathcal F_h(t,U)$。例如经典显式 RK4 每步需要计算四个阶段右端，然后加权组合，因而是**四阶时间积分方法**（满足相应光滑性条件）。

注意两点：

- 这里的“阶段”不是 Picard 或 Newton 的“非线性迭代”，即使都是一步内多次计算，其数学角色不同。

- 显式 RK4 虽具有高阶精度，却仍有稳定步长限制，尤其是扩散刚性较强时；某些隐式 RK 格式则需要在每个阶段求解代数/非线性系统。

因此**高阶 ≠ 一定更快，也 ≠ 一定更稳定**。比较时应该固定目标误差，测量总成本。

## 8\.6 怎样从 $U^{n+1}$ 重新得到有限元函数？

每做完一步计算，就得到了新时间层的系数向量。利用**同一组固定空间基函数**重建：

$\displaystyle \boxed{u_h(\mathbf x,t_{n+1})=\sum_{j=1}^N U_j^{n+1}\varphi_j(\mathbf x).}$

因而必须分清：

- 网格和基函数固定 $\Longrightarrow M,K$ 在本专题的假设下固定；

- 系数 $U^n$ 随时间改变 $\Longrightarrow u_h$ 改变；

- 对流矩阵 $B(U^n)$ 或 $G(U^n)$ 依赖解，下一物理时间步**一般需要重新装配/更新**；

- 若以后采用**移动网格**，基函数、质量矩阵、刚度矩阵乃至时间导数的处理都可能改变；那不属于当前固定网格问题。

---

## 8\.7 区分三层“求解”：时间积分、非线性求根与线性代数

当我们说“做一步”，至少可能指三件不同的事：

|层次|解决的问题|典型方法|指标|
|---|---|---|---|
|物理时间推进|$t_n\to t_{n+1}$|Euler、BDF、RK|$n$|
|同一时间层内的非线性求根|使离散残差 $R(U)=0$|Picard、Newton|$k$|
|当前线性子问题的求解|求 $Ax=b$ 或 $J\delta U=-R$|LU、GMRES、适用时的 CG|线性迭代指标（可另记 $\ell$）|

因此“后向 Euler \+ Newton \+ LU”是完全合法的组合；“后向 Euler 时间差分 \+ 冻结旧梯度 \+ LU”也是一条路线，但求解的是**另一种时间离散方程**。

对于一致质量矩阵、固定网格，$M$ 对称正定；标准扩散 $K$ 对称半正定（在适当齐次 Dirichlet 约束后通常正定），故 $A_\nu=M/\Delta t+\nu K$ 对称正定。**加入对流矩阵后整体矩阵不一定仍有相同性质**：$B(W)$ 一般非对称；$G(W)$ 虽是对称的加权质量型矩阵，却可能不定。因此选择 CG 等迭代线性求解器前要检查相应前提，不能一律采用。

## 8\.8 算一遍真正的时间步：从有限元系数到新系数

为了专注于线性代数，暂取一维 Burgers 方程

$\displaystyle u_t+uu_x-\nu u_{xx}=0, \qquad x\in(0,1),\quad u(0,t)=u(1,t)=0.$

将 $[0,1]$ 均分成三个 P1 区间单元，两端为零 Dirichlet 自由度，只保留两个内部节点。对三个长度为1/3的P1单元，按第四章的方法装配并消去端点，得到

$\displaystyle M=\begin{pmatrix}2/9&1/18\\1/18&2/9\end{pmatrix}, \qquad K=\begin{pmatrix}6&-3\\-3&6\end{pmatrix}.$

取已知旧解和参数

$\displaystyle U^n=\begin{pmatrix}1\\0\end{pmatrix}, \qquad\Delta t=0.1,\qquad\nu=0.1, \qquad F=0.$

采用**冻结旧梯度**的 $G$ 格式：

$\displaystyle \left[\frac M{\Delta t}+\nu K+G(U^n)\right]U^{n+1}=\frac M{\Delta t}U^n.$

三个单元内旧解的导数依次是 $3,-3,0$，由局部矩阵装配得到

$\displaystyle \boxed{G(U^n)=\begin{pmatrix}0&-1/6\\-1/6&-1/3\end{pmatrix}.}$

因而

$\displaystyle A=\frac M{\Delta t}+\nu K+G(U^n) =\begin{pmatrix}127/45&4/45\\4/45&112/45\end{pmatrix},$

$\displaystyle b=\frac M{\Delta t}U^n =\begin{pmatrix}20/9\\5/9\end{pmatrix}.$

解 $AU^{n+1}=b$，得到（可代回矩阵方程核对）

$\displaystyle \boxed{U^{n+1}=\begin{pmatrix}25/32\\25/128\end{pmatrix} \approx\begin{pmatrix}0.78125\\0.1953125\end{pmatrix}.}$

**究竟算出了什么？** 我们算出的不是一个“新网格”，而是**同一有限元空间中下一时刻的两个内部节点系数**。把它们代回 P1 基函数展开，就得到整个区间上的下一时刻近似函数。

下一步 $n+1\to n+2$ 时，$M,K$ 仍然是原矩阵；但 $G(U^{n+1})$ 一般与 $G(U^n)$ 不同，所以整体矩阵通常需要更新。特别注意：若采用显式对流、隐式扩散 $A_\nu U^{n+1}=r^n-C(U^n)$，且 $\Delta t,\nu$ 固定，**左侧 **$A_\nu$** 可以复用**；是否可复用矩阵依赖所选格式。

## 8\.9 若采用全隐式：先定义非线性残差

保持后向 Euler 的时间差分和新时间层扩散，但对流也取新时间层：

$\displaystyle A_\nu U^{n+1}+C(U^{n+1})=r^n.$

为简化书写，把当前时间步待求变量临时记为 $W$，定义

$\displaystyle \boxed{R(W)=A_\nu W+C(W)-r^n.}\quad(B.5)$

目标是求

$\displaystyle \boxed{R(W)=0.}\quad(B.6)$

这是一组**非线性代数方程**。$W$ 是未知系数向量，$R(W)$ 是各个离散方程的残差。求出满足 $(\mathrm{B}.6)$ 的向量后，才把它接受为 $U^{n+1}$。

## 8\.10 Picard 迭代：反复求线性系统，逼近同一个新时间层

回忆 $C(W)=B(W)W$。令 $W^{(0)}=U^n$ 为新时间层的初始猜测；在第 $k$ 次 Picard 迭代中，**仅冻结矩阵系数**，求解

$\displaystyle \boxed{[A_\nu+B(W^{(k)})]W^{(k+1)}=r^n.}\quad(B.7)$

然后用得到的 $W^{(k+1)}$ 更新系数，反复迭代，直到满足非线性收敛判据。若它收敛到 $W^*$ 且相关运算连续，则由 $B(W^*)W^*=C(W^*)$ 可知其固定点满足原来的全隐式残差方程。

也可以定义以 $G(W^{(k)})$ 冻结的 Picard 型迭代；两种迭代映射不同，收敛性质可能不同，但若收敛到一致的固定点，仍满足相同的全隐式方程。

**必须区分下面两件事：**

- **旧物理时间层一次冻结**：只解一次 $[A_\nu+B(U^n)]U^{n+1}=r^n$，得到的是一个**线性化时间格式**，不会自动满足完全隐式残差方程。

- **新时间层内反复 Picard**：使用 $B(W^{(0)}),B(W^{(1)}),\ldots$ 连续更新，目的是求解**原完全隐式方程**。

迭代 $W^{(0)}\to W^{(1)}\to W^{(2)}$ 的过程中，**物理时间始终是同一个 **$t_{n+1}$，并没有推进到 $t_{n+2}$。

## 8\.11 Newton 迭代：对整个残差求导，再求修正量

### 8\.11\.1 从一元方程理解

考虑上一节的一元后向 Euler 方程：

$\displaystyle R(z)=z+0.1z^2-1=0.$

在当前猜测 $z_k$ 附近，使用一阶 Taylor 近似：

$\displaystyle R(z_k+\delta z)\approx R(z_k)+R'(z_k)\delta z.$

令线性近似残差为零，得到

$\displaystyle \boxed{R'(z_k)\delta z=-R(z_k),\qquad z_{k+1}=z_k+\delta z.}$

取 $z_0=1$，则 $R(z_0)=0.1$、$R'(z_0)=1.2$，从而 $z_1=11/12\approx0.9166667$。第二轮后 $z_2\approx0.9160798122$，已经非常接近方程的正根 $z^*\approx0.9160797831$。

Newton 不是简单“把一个因子固定成上次迭代值”，而是对**整个残差函数**求导并修正当前猜测。

### 8\.11\.2 推广到有限维方程组

令 $W^{(k)}$ 为当前猜测，Jacobian（雅可比）矩阵为 $J_R(W^{(k)})$，求解

$\displaystyle \boxed{J_R(W^{(k)})\,\delta W^{(k)}=-R(W^{(k)}),}\quad(B.8)$

然后

$\displaystyle \boxed{W^{(k+1)}=W^{(k)}+\delta W^{(k)}.}\quad(B.9)$

如果直接更新不稳定，还可以研究带阻尼的 $W^{(k+1)}=W^{(k)}+\lambda_k\delta W^{(k)}$（$0<\lambda_k\le1$）；阻尼是否收敛仍需相应条件。

### 8\.11\.3 为什么雅可比中同时出现 $B$ 和 $G$？

定义**连续层面**的非线性算子

$\displaystyle \boxed{\mathcal N(u)=u(\boldsymbol\beta\cdot\nabla u).}$

对 $u$ 加入扰动 $\varepsilon w$：

$\displaystyle \begin{aligned} \mathcal N(u+\varepsilon w) &=(u+\varepsilon w)\bigl[\boldsymbol\beta\cdot\nabla(u+\varepsilon w)\bigr]\\ &=\mathcal N(u) +\varepsilon\Bigl[w(\boldsymbol\beta\cdot\nabla u)+u(\boldsymbol\beta\cdot\nabla w)\Bigr] +\varepsilon^2w(\boldsymbol\beta\cdot\nabla w). \end{aligned}$

因此其方向导数（例如在有界二维Lipschitz区域的H₀¹到H⁻¹映射下）是

$\displaystyle \boxed{\mathrm D\mathcal N(u)[w] =w(\boldsymbol\beta\cdot\nabla u) +u(\boldsymbol\beta\cdot\nabla w).}\quad(B.10)$

这里 $w$ 是**扰动方向**，不是弱形式中的测试函数 $v$。将该式代入有限元弱形式，并对离散系数求导，可得

$\displaystyle \boxed{J_C(W)=G(W)+B(W),}\quad(B.11)$

所以完整残差的雅可比矩阵为

$\displaystyle \boxed{J_R(W)=\frac M{\Delta t}+\nu K+G(W)+B(W).}\quad(B.12)$

这说明专题 A 的两种矩阵还有一个更深的作用：它们的**和**构成了完全隐式 Burgers 对流向量的雅可比矩阵。需要强调，这个结论基于当前的二次对流非线性及其一致离散；对于其他 PDE，不能不经推导就照搬。

## 8\.12 如何判断迭代收敛？怎样区分三种误差？

### 8\.12\.1 非线性代数迭代的停止判据

在固定时间层，求解 $R(W)=0$，一种常见的实用做法是**同时检查增量和残差**：

$\displaystyle \boxed{\frac{\|W^{(k+1)}-W^{(k)}\|}{1+\|W^{(k+1)}\|} <\varepsilon_{\mathrm{step}},}$

$\displaystyle \boxed{\frac{\|R(W^{(k+1)})\|}{1+\|r^n\|} <\varepsilon_{\mathrm{res}}.}$

同时设置最大迭代次数，检查线性子问题的求解成功与否；实际应用还可以使用绝对/相对容差及合适的尺度。增量小并不一定说明原方程的残差小；残差小在病态问题中也不必然意味着解的实际误差小。

**特别注意**：无论 Picard 还是 Newton，判断收敛时应将近似解代入**原始全隐式残差 **$R$，而非只检查当前被冻结的线性方程残差。

### 8\.12\.2 证明“会收敛”需要什么条件？

不能无条件保证任意格式、任意网格、任意黏性、任意初值都收敛。

- **Picard**：若迭代映射 $\Phi$ 在某个完备集合上映回自身且满足压缩性 $\|\Phi(X)-\Phi(Y)\|\le q\|X-Y\|$（$0<q<1$），可用 Banach 不动点定理证明固定点存在、唯一且迭代收敛。对于具体 Burgers 离散系统，必须证明这些条件，不是写出 Picard 公式就自动成立。

- **Newton**：在根附近，若雅可比可逆、足够光滑且初值足够接近根，经典定理给出局部快速收敛，适当条件下可达二次收敛 $\|e_{k+1}\|\le C\|e_k\|^2$；但它不是任意起点都收敛的全局保证。

- **线性子求解**：如果用迭代线性求解器，还要考虑其自身是否收敛、容差是否足够严格；这与 Picard/Newton 的非线性收敛是另一个层次。

### 8\.12\.3 三种“正确”不是同一件事

1. **代数求解误差**：当前时间离散方程 $R(U^{n+1})=0$ 是否求得足够准确？

2. **时间离散误差**：$U^n$ 与连续时间的半离散 ODE 解 $U(t_n)$ 相差多大？希望 $\Delta t\to0$ 时趋零。

3. **空间离散误差**：有限元半离散解 $u_h$ 与真实 PDE 解 $u$ 相差多大？希望 $h\to0$ 时趋零。

Newton 残差很小，只能说明**代数问题**较准确地求解了，不能单独证明 PDE 近似已经足够准确。专题 D 将通过制造解、控制变量、网格与时间加密分别验证空间和时间收敛。

### 8\.12\.4 怎样评价速度与成本？

直观上，计算成本包括：

$\displaystyle \boxed{\text{总成本}\approx\text{时间步数}\times \text{平均非线性迭代次数}\times\text{平均线性求解成本} +\text{装配等开销}.}$

这是一个粗略的**成本分解**，不是通用复杂度定理。高阶时间方法每步可能更贵，但若达到同样误差只需更少步数，总成本反而可能更低。比较算法应记录：**同等误差下的总耗时、线性/非线性迭代次数、矩阵更新次数、存储开销和失败率**。

## 8\.13 一套完整的数学求解流程（固定网格、后向 Euler、Newton）

以下不涉及任何特定程序库，只整理计算次序。

**输入**：$\Omega$、终止时间 $T$、$u_0$、$f$、边界数据、$\nu$、有限元次数、网格尺寸 $h$、时间步长 $\Delta t$、非线性/线性容差。

**阶段一：空间准备**

1. 网格剖分，确定有限元自由度和基函数 $\varphi_j$；施加边界条件（齐次 Dirichlet 时，可将边界自由度消元）。

2. 按单元装配质量矩阵 $M$ 和刚度矩阵 $K$；必要时准备数值求积规则。

3. 从 $u_0$ 插值或投影得到 $U^0$。当前物理时间层 $n=0$。

**阶段二：一个物理时间步 **$t_n\to t_{n+1}$

1. 根据 $t_{n+1}$ 计算载荷 $F^{n+1}$；定义 $A_\nu=M/\Delta t+\nu K$ 和 $r^n=MU^n/\Delta t+F^{n+1}$。

2. 若选**完全隐式后向 Euler**，定义 $R(W)=A_\nu W+C(W)-r^n$，初值取 $W^{(0)}=U^n$。

3. 在同一 $t_{n+1}$ 内，对 $k=0,1,\dots$ 重复：

    - 计算当前对流向量 $C(W^{(k)})$，形成**原残差**$R(W^{(k)})$；

    - 装配雅可比 $J_R(W^{(k)})=A_\nu+B(W^{(k)})+G(W^{(k)})$；

    - 求解线性修正方程 $J_R(W^{(k)})\delta W=-R(W^{(k)})$；

    - 更新 $W^{(k+1)}=W^{(k)}+\delta W$（必要时阻尼）；

    - 重新计算原残差并检查停止判据。

4. 若迭代成功，则接受 $U^{n+1}=W^{(k+1)}$；若失败，报告失败并采取更小步长、改进初值等应对措施，**不能悄悄把不收敛的近似值当作正确结果**。

**阶段三：推进与检验**

1. 重建 $u_h^{n+1}=\sum_jU_j^{n+1}\varphi_j$，记录时间、残差、迭代次数与需要的诊断量。

2. 若 $t_{n+1}<T$，令 $n\leftarrow n+1$，重复阶段二，直至到达终止时间。

3. 用制造解和网格/时间加密实验区分空间误差、时间误差、代数误差；如需研究稳定性，还应监测适当的范数、能量及非物理振荡。

**如果选择一次冻结旧时间层的线性化格式**，则阶段二中的第 5—7 步被一次已知矩阵的线性求解替代。例如冻结旧梯度：

$\displaystyle \boxed{\left[A_\nu+G(U^n)\right]U^{n+1}=r^n.}$

它更便宜，但这时应明确记录：求解的是**线性化时间格式**，而不是通过 Newton 精确求解了原完全隐式方程。

---

## 8\.14 连起计算流程与下一步问题

### 从连续问题到计算结果

**空间有限元离散**先将 PDE 转为有限维非线性 ODE：

$\displaystyle \boxed{M\dot U+\nu KU+C(U)=F.}$

**时间离散**再把 ODE 变成某一时间层的代数方程。若对流完全隐式，则代数方程是非线性的；若显式处理或冻结旧时间层某个因子，则得到关于新解的线性系统。

**Picard 和 Newton**不是另一种独立的时间积分法，而是在同一个物理时间层内部求解非线性代数方程的方法。它们内部可以反复调用 LU、GMRES 等线性求解器。

因此从连续问题到计算结果，应始终知道自己正处于哪一层：

$\displaystyle \boxed{ \begin{gathered} \text{PDE 初边值问题}\\ \Downarrow\ \text{弱形式与空间有限元}\\ M\dot U+\nu KU+C(U)=F\\ \Downarrow\ \text{Euler / BDF / RK 等时间积分}\\ \text{线性代数方程或非线性残差方程}\\ \Downarrow\ \text{必要时 Picard / Newton}\\ \text{线性方程求解}\\ \Downarrow\\ U^{n+1}\ \Longrightarrow\ u_h(\mathbf x,t_{n+1}). \end{gathered}}$

### 比较格式时先固定什么

为了公平比较四种对流时间层处理，专题 C 第一轮建议**固定**：网格、P1/P2 次数、边界条件、$\Delta t$、$\nu$、扩散项隐式处理、载荷与初始条件；**只改变**

$\displaystyle u^n(\boldsymbol\beta\cdot\nabla u^n),\quad u^n(\boldsymbol\beta\cdot\nabla u^{n+1}),\quad u^{n+1}(\boldsymbol\beta\cdot\nabla u^n),\quad u^{n+1}(\boldsymbol\beta\cdot\nabla u^{n+1}).$

**这里的 **$u^n$** 表示解的时间层，不是黏性系数 **$\nu$**；若排版容易混淆，应在讲义中始终将黏性写作 **$\nu$**，并为 **$u$** 保留明确上标。**

四种格式产生不同的矩阵结构和离散方程，其计算成本、稳定性与误差可能不同。专题 C 需要在明确理论假设的前提下分析，再由专题 D 的实验验证，而不是先断言哪个格式更好。

### 复习问题

1. 为什么得到 $M\dot U+\nu KU+C(U)=F$ 后还不能直接把它当作静态 $AU=b$？

2. 对于固定网格，$M,K$ 为什么通常不变，而 $B(U^n),G(U^n)$ 一般会变？

3. “只冻结一次 $U^n$”与“在新时间层内不断迭代 $W^{(k)}$”有什么不同？

4. 为什么 Newton 的雅可比中同时出现 $B(W)$ 和 $G(W)$？

5. 线性系统的残差、非线性迭代误差、时间误差、空间误差分别代表什么？

6. 高阶 BDF2 或 RK4 是否一定比一阶 Euler 更省时间？需要固定什么比较指标？

---

**下一步**：进入专题 C，系统研究对流项两因子的旧/新时间层组合、扩散项的隐显选择及黏性系数变化造成的稳定性与误差差异。专题 D 再使用制造解、网格加密和时间步加密进行定量验证。
