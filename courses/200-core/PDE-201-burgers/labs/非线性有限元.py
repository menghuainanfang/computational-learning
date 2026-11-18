"""固定网格标量黏性Burgers：FEALPy空间、求积装配与四种Euler对流处理。"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve
from fealpy.backend import backend_manager as bm
from fealpy.mesh import IntervalMesh, TriangleMesh
from fealpy.functionspace import LagrangeFESpace


class BurgersFEM:
    """单位区间/正方形，全齐次Dirichlet；二维方向beta=(1,1)。"""
    def __init__(self, n=8, p=1, dim=2):
        bm.set_backend('numpy')
        if dim == 1:
            cells = np.column_stack((np.arange(n), np.arange(1, n+1)))
            self.mesh = IntervalMesh(np.linspace(0, 1, n+1)[:, None], cells)
        elif dim == 2:
            self.mesh = TriangleMesh.from_box([0, 1, 0, 1], nx=n, ny=n)
        else:
            raise ValueError('本课程实现单位区间和单位正方形')
        self.space = LagrangeFESpace(self.mesh, p=p)
        self.dim, self.p = dim, p
        bc, w = self.mesh.quadrature_formula(3*p+2).get_quadrature_points_and_weights()
        self.points = self.mesh.bc_to_point(bc)
        self.phi = np.broadcast_to(self.space.basis(bc), (self.mesh.number_of_cells(), len(w), self.space.number_of_local_dofs()))
        self.grad = self.space.grad_basis(bc)
        self.ds = self.grad.sum(axis=-1)
        self.weights = self.mesh.entity_measure('cell')[:, None]*w[None, :]
        self.dofs = self.space.cell_to_dof()
        self.N = self.space.number_of_global_dofs()
        self.nodes = self.space.interpolation_points()
        self.boundary = np.any(np.isclose(self.nodes, 0) | np.isclose(self.nodes, 1), axis=-1)
        self.I = np.flatnonzero(~self.boundary)
        self.rows = np.broadcast_to(self.dofs[:, :, None], (len(self.dofs), self.dofs.shape[1], self.dofs.shape[1])).ravel()
        self.cols = np.broadcast_to(self.dofs[:, None, :], (len(self.dofs), self.dofs.shape[1], self.dofs.shape[1])).ravel()
        self.M = self.matrix(np.einsum('cq,cqi,cqj->cij', self.weights, self.phi, self.phi))
        self.K = self.matrix(np.einsum('cq,cqid,cqjd->cij', self.weights, self.grad, self.grad))

    def matrix(self, local):
        return coo_matrix((local.ravel(), (self.rows, self.cols)), shape=(self.N, self.N)).tocsr()

    def vector(self, local):
        result = np.zeros(self.N)
        np.add.at(result, self.dofs.ravel(), local.ravel())
        return result

    def nonlinear(self, U):
        """C(U), B(U), G(U)使用同一求积，满足C=B U=G U及J_C=B+G。"""
        local = U[self.dofs]
        value = np.einsum('cqi,ci->cq', self.phi, local)
        derivative = np.einsum('cqi,ci->cq', self.ds, local)
        B = self.matrix(np.einsum('cq,cq,cqi,cqj->cij', self.weights, value, self.phi, self.ds))
        G = self.matrix(np.einsum('cq,cq,cqi,cqj->cij', self.weights, derivative, self.phi, self.phi))
        C = self.vector(np.einsum('cq,cq,cq,cqi->ci', self.weights, value, derivative, self.phi))
        return C, B, G

    def load(self, source, t):
        return self.vector(np.einsum('cq,cq,cqi->ci', self.weights, source(self.points, t), self.phi))

    def norm(self, U):
        return float(np.sqrt(max(0., U @ (self.M @ U))))

    def error(self, U, exact, t):
        values = np.einsum('cqi,ci->cq', self.phi, U[self.dofs])
        return float(np.sqrt(np.sum(self.weights*(values-exact(self.points, t))**2)))

    def step(self, U, dt, nu, load, scheme='newton', tol=1e-11, maxiter=40):
        A = self.M/dt+nu*self.K
        rhs = self.M@U/dt+load
        C, B, G = self.nonlinear(U)
        I = self.I
        if scheme in ('imex', 'front', 'gradient'):
            lhs, target = (A, rhs-C) if scheme == 'imex' else (A+(B if scheme == 'front' else G), rhs)
            V = np.zeros(self.N)
            V[I] = spsolve(lhs[I][:, I], target[I])
            algebraic = np.linalg.norm((lhs@V-target)[I])/(1+np.linalg.norm(target[I]))
            iterations = 1
        else:
            if scheme not in ('newton', 'picard'):
                raise ValueError('未知时间层或非线性求解方式')
            V = U.copy()
            for k in range(maxiter):
                C, B, G = self.nonlinear(V)
                residual = A@V+C-rhs
                delta = np.zeros(self.N)
                if scheme == 'newton':
                    delta[I] = spsolve((A+B+G)[I][:, I], -residual[I])
                    trial = V+delta
                else:
                    trial = np.zeros(self.N)
                    trial[I] = spsolve((A+B)[I][:, I], rhs[I])
                    delta = trial-V
                V = trial
                C, _, _ = self.nonlinear(V)
                algebraic = np.linalg.norm((A@V+C-rhs)[I])/(1+np.linalg.norm(rhs[I]))
                if algebraic < tol and np.linalg.norm(delta)/(1+np.linalg.norm(V)) < tol:
                    iterations = k+1
                    break
            else:
                raise RuntimeError('当前时间步非线性残差/增量未收敛')
        if not np.all(np.isfinite(V)):
            raise RuntimeError('时间步产生非有限系数')
        full = np.linalg.norm((A@V+self.nonlinear(V)[0]-rhs)[I])/(1+np.linalg.norm(rhs[I]))
        return V, dict(iterations=iterations, residual=float(algebraic), full_residual=float(full), boundary=float(np.max(np.abs(V[self.boundary]))))

    def run(self, initial, source, dt=.01, end=.1, nu=.1, scheme='newton'):
        steps = int(round(end/dt))
        if not np.isclose(steps*dt, end):
            raise ValueError('终止时间应为时间步的整数倍')
        U = initial(self.nodes, 0.).copy()
        U[self.boundary] = 0.
        history = []
        for n in range(steps):
            old = U.copy()
            t = (n+1)*dt
            U, info = self.step(old, dt, nu, self.load(source, t), scheme)
            C, _, _ = self.nonlinear(U if scheme in ('newton', 'picard') else old)
            if scheme == 'front': C = self.nonlinear(old)[1]@U
            if scheme == 'gradient': C = self.nonlinear(old)[2]@U
            load = self.load(source, t)
            info.update(t=t, energy=float(U@(self.M@U)/2), minimum=float(U.min()), maximum=float(U.max()),
                balance=float((U@(self.M@U)-old@(self.M@old))/2+(U-old)@(self.M@(U-old))/2+dt*nu*U@(self.K@U)+dt*U@C-dt*U@load))
            history.append(info)
        return U, history


def manufactured(dim, nu):
    def exact(x, t):
        return np.exp(-t)*np.prod(np.sin(np.pi*x), axis=-1)
    def source(x, t):
        u = exact(x, t)
        derivatives = np.zeros_like(u)
        for d in range(dim):
            factors = np.sin(np.pi*x).copy()
            factors[..., d] = np.pi*np.cos(np.pi*x[..., d])
            derivatives += np.exp(-t)*np.prod(factors, axis=-1)
        return (-1+nu*dim*np.pi**2)*u+u*derivatives
    return exact, source
