"""双分量稳态输运：反应块、分量排列、边界与积分误差。"""
from pathlib import Path
import csv,json
import numpy as np
from scipy.sparse import kron,eye
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian
from fealpy.mesh import TriangleMesh
from fealpy.functionspace import LagrangeFESpace,TensorFunctionSpace
from fealpy.fem import BilinearForm,LinearForm,BlockForm,DirichletBC,ScalarDiffusionIntegrator,ScalarConvectionIntegrator,ScalarMassIntegrator,ScalarSourceIntegrator
from fealpy.solver import spsolve
C=np.array([[2.,-1.],[-.5,2.]])
@cartesian
def exact(p):
    s=np.sin(np.pi*p[...,0])*np.sin(np.pi*p[...,1])
    return np.stack((1+s,2+2*s),axis=-1)
@cartesian
def gradient(p):
    x,y=p[...,0],p[...,1]
    g=np.stack((np.pi*np.cos(np.pi*x)*np.sin(np.pi*y),np.pi*np.sin(np.pi*x)*np.cos(np.pi*y)),axis=-1)
    return np.stack((g,2*g),axis=-2)
@cartesian
def velocity(p):return np.broadcast_to(np.array([1.,2.]),p.shape)
def component_source(r,coef):
    @cartesian
    def source(p):
        u=exact(p);g=gradient(p)
        transport=g[...,r,0]+2*g[...,r,1]
        return 2*np.pi**2*(u[...,r]-(1 if r==0 else 2))+transport+np.einsum('s,...s->...',coef[r],u)
    return source

def solve(n,coef=C):
    bm.set_backend('numpy')
    mesh=TriangleMesh.from_box([0,1,0,1],nx=n,ny=n)
    scalar=LagrangeFESpace(mesh,p=1)
    space=TensorFunctionSpace(scalar,shape=(2,-1))
    forms=[];rhs=[]
    for r in range(2):
        row=[]
        for s in range(2):
            f=BilinearForm(scalar)
            f.add_integrator(ScalarMassIntegrator(coef=float(coef[r,s]),q=5))
            if r==s:
                f.add_integrator(ScalarDiffusionIntegrator(coef=1.,q=5))
                f.add_integrator(ScalarConvectionIntegrator(coef=velocity,q=5))
            row.append(f)
        forms.append(row)
        load=LinearForm(scalar)
        load.add_integrator(ScalarSourceIntegrator(source=component_source(r,coef),q=7))
        rhs.append(load.assembly())
    A=BlockForm(forms).assembly();F=np.concatenate(rhs)
    bd=np.tile(scalar.is_boundary_dof(),2)
    Ab,Fb=DirichletBC(space,gd=exact,threshold=bd).apply(A,F)
    uh=space.function();uh[:]=spsolve(Ab,Fb,solver='scipy')
    bc,w=mesh.quadrature_formula(7).get_quadrature_points_and_weights()
    pts=mesh.bc_to_point(bc);measure=mesh.entity_measure('cell')
    e=uh.value(bc)-exact(pts);g=uh.grad_value(bc)-gradient(pts)
    l2=np.sqrt(np.einsum('cqr,q,c->r',e**2,w,measure));h1=np.sqrt(np.einsum('cqrd,q,c->r',g**2,w,measure))
    data=exact(scalar.interpolation_points()).T.reshape(-1)
    row=dict(n=n,dofs=2*scalar.number_of_global_dofs(),L2_1=float(l2[0]),L2_2=float(l2[1]),H1_1=float(h1[0]),H1_2=float(h1[1]),L2=float(np.linalg.norm(l2)),H1=float(np.linalg.norm(h1)),boundary=float(np.max(np.abs(uh[:][bd]-data[bd]))),residual=float(np.linalg.norm((A@uh[:]-F)[~bd])/max(1.,np.linalg.norm(F[~bd]))))
    return row,A,scalar,uh

def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams['font.sans-serif']=['Microsoft YaHei'];plt.rcParams['axes.unicode_minus']=False
    root=Path(__file__).resolve().parents[1];out=root/'labs/results';figs=root/'notes/figures'
    out.mkdir(parents=True,exist_ok=True)
    figs.mkdir(parents=True,exist_ok=True)
    rows=[]
    for n in (8,16,32,64):
        row,_,_,_=solve(n)
        row['L2_order']=np.log2(rows[-1]['L2']/row['L2']) if rows else None
        row['H1_order']=np.log2(rows[-1]['H1']/row['H1']) if rows else None
        assert row['boundary']<1e-12 and row['residual']<1e-10
        rows.append(row)
    assert 1.9<rows[-1]['L2_order']<2.1 and .95<rows[-1]['H1_order']<1.05
    row,A,scalar,uh=solve(2)
    mass=BilinearForm(scalar);mass.add_integrator(ScalarMassIntegrator(coef=1.,q=5));M=mass.assembly().to_scipy()
    diag=BilinearForm(scalar);diag.add_integrator(ScalarDiffusionIntegrator(coef=1.,q=5));diag.add_integrator(ScalarConvectionIntegrator(coef=velocity,q=5));D=diag.assembly().to_scipy()
    expected=kron(eye(2),D)+kron(C,M)
    err=np.max(np.abs((A.to_scipy()-expected).data),initial=0.)
    N=scalar.number_of_global_dofs();perm=np.arange(2*N).reshape(2,N).T.reshape(-1)
    U=uh[:];nodeU=U[perm]
    permutation_error=np.max(np.abs((A.to_scipy()[perm][:,perm]@nodeU)-(A@U)[perm]))
    _,Az,_,_=solve(2,np.diag([2.,2.]))
    Az_scipy=Az.to_scipy()
    decoupling=max(np.max(np.abs(Az_scipy[:N,N:].data),initial=0.),
                   np.max(np.abs(Az_scipy[N:,:N].data),initial=0.))
    assert max(err,permutation_error,decoupling)<1e-12
    checks=dict(kronecker_error=float(err),permutation_error=float(permutation_error),zero_offdiag=float(decoupling),reaction_symmetric_eigenvalues=np.linalg.eigvalsh((C+C.T)/2).tolist())
    with (out/'收敛结果.csv').open('w',encoding='utf-8-sig',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=rows[0]);wr.writeheader();wr.writerows(rows)
    (out/'结构核对.json').write_text(json.dumps(checks,ensure_ascii=False,indent=2),encoding='utf8')
    fig,axes=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for ax,key,power in zip(axes,['L2','H1'],[2,1]):
        h=np.array([1/r['n'] for r in rows]);e=np.array([r[key] for r in rows]);ax.loglog(h,e,'o-',label='实测向量积分误差');ax.loglog(h,e[-1]*(h/h[-1])**power,'--',label=f'{power}阶参考线');ax.set_xlabel('h=1/n');ax.set_ylabel(r'$L^2$ 误差' if key=='L2' else r'$H^1$ 半范数误差');ax.legend();ax.grid(alpha=.3)
    fig.savefig(figs/'向量误差收敛.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(7,6));im=ax.imshow(A.to_scipy().toarray(),cmap='coolwarm');ax.axhline(N-.5,color='k');ax.axvline(N-.5,color='k');ax.set_xlabel('未知分量与空间自由度（列）');ax.set_ylabel('测试分量与空间自由度（行）');ax.set_title('四块矩阵：非对角块来自反应耦合');fig.colorbar(im,ax=ax);fig.tight_layout();fig.savefig(figs/'四块矩阵.png',dpi=150);plt.close(fig)
    fig,ax=plt.subplots(figsize=(9,3));ax.axis('off');ax.text(.05,.75,'分量优先： [u1[0], u1[1], u1[2] | u2[0], u2[1], u2[2]]',fontsize=13);ax.text(.05,.4,'节点优先： [u1[0], u2[0] | u1[1], u2[1] | u1[2], u2[2]]',fontsize=13);ax.text(.05,.05,'同一问题的两种编号：向量、矩阵行列、载荷和边界一起置换',fontsize=12);fig.tight_layout();fig.savefig(figs/'分量排列.png',dpi=150);plt.close(fig)
    print(json.dumps(dict(results=rows,checks=checks),ensure_ascii=False,indent=2))
if __name__=='__main__':main()
