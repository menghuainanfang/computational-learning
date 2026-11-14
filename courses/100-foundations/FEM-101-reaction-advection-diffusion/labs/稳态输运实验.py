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
    form.add_integrator(ScalarDiffusionIntegrator(coef=1., q=5))
    @cartesian
    def velocity(p):
        return np.broadcast_to(np.array([1., 2.]), p.shape)
    form.add_integrator(ScalarConvectionIntegrator(coef=velocity, q=5))
    form.add_integrator(ScalarMassIntegrator(coef=1., q=5))
    A = form.assembly()
    load = LinearForm(space)
    load.add_integrator(ScalarSourceIntegrator(source=source, q=5))
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

def 核对装配():
    """分别核对几何梯度、局部对流矩阵及重复坐标累加。"""
    bm.set_backend('numpy')
    from scipy.sparse import coo_matrix
    mesh = TriangleMesh.from_box([0, 1, 0, 1], nx=2, ny=2)
    space = LagrangeFESpace(mesh, p=1)
    @cartesian
    def velocity(p):
        return np.broadcast_to(np.array([1., 2.]), p.shape)
    integrator = ScalarConvectionIntegrator(coef=velocity, q=5)
    local = integrator.assembly(space)
    dofs = space.cell_to_dof()
    rows = np.broadcast_to(dofs[:, :, None], local.shape)
    cols = np.broadcast_to(dofs[:, None, :], local.shape)
    count = space.number_of_global_dofs()
    manual = coo_matrix((local.ravel(), (rows.ravel(), cols.ravel())),
                        shape=(count, count)).tocsr()
    form = BilinearForm(space)
    form.add_integrator(ScalarConvectionIntegrator(coef=velocity, q=5))
    assembled = form.assembly().to_scipy()
    difference = manual - assembled
    assembly_error = np.max(np.abs(difference.data), initial=0.)
    bc = np.array([[1/3, 1/3, 1/3]])
    grad = space.grad_basis(bc)[:, 0]
    areas = mesh.entity_measure('cell')
    expected = np.broadcast_to(
        areas[:, None, None] / 3 * (grad @ np.array([1., 2.]))[:, None, :],
        local.shape)
    convection_error = np.max(np.abs(local - expected))
    nodes = mesh.entity('node')[mesh.entity('cell')]
    jacobian = np.stack((nodes[:, 1]-nodes[:, 0], nodes[:, 2]-nodes[:, 0]), axis=-1)
    reference_grad = np.array([[-1., -1.], [1., 0.], [0., 1.]])
    transformed = np.einsum('ij,cjk->cik', reference_grad, np.linalg.inv(jacobian))
    gradient_error = np.max(np.abs(transformed-grad))
    # 斜三角形独立检验：u=x+2y 的梯度应为 (1,2)。
    skew_j = np.array([[1., 1.], [0., 1.]])
    skew_gradient = np.linalg.solve(skew_j.T, np.array([1., 3.]))
    assert np.allclose(skew_gradient, [1., 2.])
    assert max(assembly_error, convection_error, gradient_error) < 1e-12
    return dict(local_shape=list(local.shape), cell_to_dof_shape=list(dofs.shape),
                gradient_shape=list(grad.shape), assembly_error=float(assembly_error),
                convection_error=float(convection_error), gradient_error=float(gradient_error),
                skew_gradient=skew_gradient.tolist())


def 绘制几何图(folder):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8))
    triangles = [np.array([[0,0],[1,0],[0,1]]),
                 np.array([[0.2,0.1],[2.2,0.4],[0.8,1.6]])]
    names = [['0: (0,0)','1: (1,0)','2: (0,1)'],
             [r'$P_0$',r'$P_1$',r'$P_2$']]
    for ax, points, labels in zip(axes,triangles,names):
        ax.add_patch(Polygon(points,facecolor='#dcecf5',edgecolor='#185572',lw=2))
        for p,label in zip(points,labels):
            ax.plot(*p,'o',color='#185572')
            ax.annotate(label,p,xytext=(4,8),textcoords='offset points',fontsize=12)
        for j in [1,2]:
            ax.annotate('',points[j],points[0],arrowprops=dict(arrowstyle='->',color='#bb5a26',lw=2))
        ax.set_aspect('equal');ax.set_xlim(-.25,2.6 if ax is axes[1] else 1.35)
        ax.set_ylim(-.25,1.95 if ax is axes[1] else 1.4)
        ax.set_xlabel('x' if ax is axes[1] else r'$\xi$')
        ax.set_ylabel('y' if ax is axes[1] else r'$\eta$')
    axes[0].set_title('参考单元：基函数统一定义')
    axes[1].set_title('物理单元：两条边组成雅可比矩阵')
    fig.suptitle(r'$x=P_0+J_K\hat{x},\quad \nabla_x\phi=J_K^{-T}\nabla_{\hat{x}}\hat\phi$',fontsize=14)
    fig.tight_layout();fig.savefig(folder/'参考单元映射.png',dpi=160);plt.close(fig)
    fig, axes=plt.subplots(1,2,figsize=(10,4))
    ax=axes[0]
    pts=np.array([[0,0],[1,0],[0,1],[1,1]])
    for indices,color,label in [([0,1,2],'#cfe5f0',r'$K_0$'),([1,3,2],'#f5dfc6',r'$K_1$')]:
        poly=pts[indices];ax.add_patch(Polygon(poly,facecolor=color,edgecolor='#444',lw=1.5))
        ax.text(*poly.mean(axis=0),label,fontsize=16,ha='center')
    for i,p in enumerate(pts):
        ax.plot(*p,'o',color='#222');ax.annotate(str(i),p,xytext=(6,5),textcoords='offset points',fontsize=14)
    ax.plot([1,0],[0,1],color='#bb5a26',lw=3)
    ax.set_xlim(-.2,1.25);ax.set_ylim(-.2,1.25);ax.set_aspect('equal');ax.axis('off')
    ax.set_title('两个单元共享全局自由度 1、2')
    axes[1].axis('off')
    axes[1].text(.02,.88,'局部编号       0     1     2',fontsize=13)
    axes[1].text(.02,.69,'K0 → 全局     0     1     2',fontsize=13)
    axes[1].text(.02,.50,'K1 → 全局     1     3     2',fontsize=13)
    axes[1].text(.02,.25,r'$A_{12}=A^{K_0}_{12}+A^{K_1}_{02}$',fontsize=17,color='#9d4420')
    axes[1].text(.02,.05,'同一全局位置的贡献相加',fontsize=13)
    fig.suptitle('局部矩阵如何装配为全局矩阵',fontsize=14)
    fig.tight_layout();fig.savefig(folder/'局部到全局装配.png',dpi=160);plt.close(fig)


def 绘制收支图(folder):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    fig,ax=plt.subplots(figsize=(8,4))
    ax.add_patch(Rectangle((.3,.25),.4,.5,facecolor='#e6f0f5',edgecolor='#185572',lw=2))
    ax.annotate('',(.45,.5),(.05,.5),arrowprops=dict(arrowstyle='->',lw=3,color='#185572'))
    ax.annotate('',(.95,.5),(.55,.5),arrowprops=dict(arrowstyle='->',lw=3,color='#185572'))
    ax.text(.1,.58,'输入通量',ha='center',fontsize=13)
    ax.text(.89,.58,'输出通量',ha='center',fontsize=13)
    ax.text(.5,.68,'控制体 W',ha='center',fontsize=15)
    ax.text(.5,.38,'源项补入 f　／　反应消耗 κu',ha='center',fontsize=13)
    ax.text(.5,.12,r'$\mathbf{J}=\mathbf{b}u-a\nabla u$',ha='center',fontsize=19)
    ax.text(.5,.02,'通量包含对流与扩散；箭头表示所示边界的输入、输出',ha='center',fontsize=12)
    ax.set_xlim(0,1);ax.set_ylim(-.04,.9);ax.axis('off')
    fig.tight_layout();fig.savefig(folder/'控制体收支.png',dpi=160);plt.close(fig)


def main():
    import csv, json
    from pathlib import Path
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    bm.set_backend('numpy')
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    base=Path(__file__).resolve().parents[1]
    figures=base/'notes'/'figures';data=base/'labs'/'results'
    figures.mkdir(exist_ok=True);data.mkdir(exist_ok=True)
    checks=核对装配()
    results=[]
    for n in (8,16,32,64):
        _,_,_,m=求解(n=n)
        m['h']=np.sqrt(2.)/n
        m['L2_order']=np.log2(results[-1]['L2']/m['L2']) if results else None
        m['H1_order']=np.log2(results[-1]['H1_semi']/m['H1_semi']) if results else None
        assert m['boundary_error']<1e-12 and m['residual']<1e-10
        results.append(m)
    assert 1.9<results[-1]['L2_order']<2.1 and .95<results[-1]['H1_order']<1.05
    with (data/'收敛结果.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(results[0]))
        writer.writeheader();writer.writerows(results)
    (data/'装配核对.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(1,2,figsize=(9,3.8))
    h=np.array([m['h'] for m in results])
    for ax,key,power,title in zip(axes,['L2','H1_semi'],[2,1],[r'$L^2$ 误差',r'$H^1$ 半范误差']):
        e=np.array([m[key] for m in results])
        ax.loglog(h,e,'o-',color='#185572',label='实测积分误差')
        ax.loglog(h,e[-1]*(h/h[-1])**power,'--',color='#bb5a26',label=f'{power} 阶参考线')
        ax.set_xlabel('最大单元直径 h');ax.set_ylabel('误差');ax.set_title(title)
        ax.set_xticks(h)
        ax.set_xticklabels([f'{v:.3f}' for v in h])
        from matplotlib.ticker import NullLocator
        ax.xaxis.set_minor_locator(NullLocator())
        ax.grid(True,which='both',alpha=.25);ax.legend()
    fig.tight_layout();fig.savefig(figures/'收敛曲线.png',dpi=160);plt.close(fig)
    绘制几何图(figures)
    绘制收支图(figures)
    print(json.dumps({'checks':checks,'convergence':results},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()

