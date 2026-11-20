"""FEALPy标量子块装配与分量优先的稳态联立求解。"""
import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.mesh import IntervalMesh, TriangleMesh, TetrahedronMesh
from fealpy.functionspace import LagrangeFESpace, TensorFunctionSpace
from fealpy.fem import (BilinearForm, LinearForm, BlockForm, DirichletBC,
                       ScalarDiffusionIntegrator, ScalarConvectionIntegrator,
                       ScalarMassIntegrator, ScalarSourceIntegrator)
from fealpy.solver import spsolve
from .模型 import component_field


def box_mesh(pde, n):
    if not isinstance(n, int) or n < 1:
        raise ValueError("n应为每方向的正整数段数")
    if pde.dim == 1:
        return IntervalMesh(bm.linspace(*pde.domain(), n+1).reshape(-1, 1),
                            bm.stack((bm.arange(n), bm.arange(1, n+1)), axis=-1))
    if pde.dim == 2:
        return TriangleMesh.from_box(pde.domain(), nx=n, ny=n)
    return TetrahedronMesh.from_box(pde.domain(), nx=n, ny=n, nz=n)


class CDRSolver:
    """共享a、b，矩阵C耦合；标准Galerkin、全Dirichlet、NumPy/CPU。"""
    def __init__(self, pde, mesh, p=1, *, q=None, linear_solver=None):
        if bm.get_current_backend().backend_name != "numpy":
            raise RuntimeError("当前实测求解器使用NumPy/CPU")
        if not isinstance(p, int) or p < 1:
            raise ValueError("p应为正整数")
        if mesh.geo_dimension() != pde.dim or mesh.top_dimension() != pde.dim:
            raise ValueError("方程与全维网格应一致")
        self.pde, self.mesh = pde, mesh
        self.scalar_space = LagrangeFESpace(mesh, p=p)
        self.space = TensorFunctionSpace(self.scalar_space, shape=(pde.components, -1))
        self.q, self.linear_solver = p+4 if q is None else q, linear_solver
        boundary = pde.is_dirichlet_boundary(self.scalar_space.interpolation_points())
        self.boundary = bm.tile(boundary, pde.components)

    def assemble(self):
        m, q = self.pde.components, self.q
        bcs, _ = self.mesh.quadrature_formula(q).get_quadrature_points_and_weights()
        a = self.pde.diffusion_coef(self.mesh.bc_to_point(bcs))
        if not bm.all(bm.isfinite(a)) or bm.any(a <= 0):
            raise ValueError("积分点的扩散系数应有限且为正")
        forms, rhs = [], []
        for r in range(m):
            row = []
            for s in range(m):
                form = BilinearForm(self.scalar_space)
                form.add_integrator(ScalarMassIntegrator(
                    coef=component_field(self.pde.reaction_coef, (r, s)), q=q))
                if r == s:
                    form.add_integrator(ScalarDiffusionIntegrator(coef=self.pde.diffusion_coef, q=q))
                    form.add_integrator(ScalarConvectionIntegrator(coef=self.pde.convection_coef, q=q))
                row.append(form)
            forms.append(row)
            load = LinearForm(self.scalar_space)
            load.add_integrator(ScalarSourceIntegrator(
                source=component_field(self.pde.source, (r,)), q=q))
            rhs.append(load.assembly())
        return BlockForm(forms).assembly(), bm.concatenate(rhs)

    def solve_steady(self):
        self.A, self.F = self.assemble()
        self.Ab, self.Fb = DirichletBC(self.space, gd=self.pde.dirichlet,
                                     threshold=self.boundary).apply(self.A, self.F)
        values = (spsolve(self.Ab.copy(), self.Fb.copy(), solver="scipy") if self.linear_solver is None
                  else self.linear_solver(self.Ab.copy(), self.Fb.copy()))
        if values.shape != self.Fb.shape or not bm.all(bm.isfinite(values)):
            raise ValueError("线性求解应返回有限、同形状的向量")
        self.uh = self.space.function()
        self.uh[:] = values
        boundary_data = self.pde.dirichlet(self.scalar_space.interpolation_points()).T.reshape(-1)
        interior = ~self.boundary
        self.diagnostics = dict(dofs=self.space.number_of_global_dofs(),
            residual=float(bm.linalg.norm((self.A @ values-self.F)[interior]))
                     / max(1., float(bm.linalg.norm(self.F[interior]))),
            boundary_error=float(bm.max(bm.abs(values[self.boundary]-boundary_data[self.boundary]))))
        return self.uh

    def errors(self):
        """返回分量误差以及由分量平方和构成的向量误差。"""
        bc, weights = self.mesh.quadrature_formula(self.q+2).get_quadrature_points_and_weights()
        points = self.mesh.bc_to_point(bc)
        e = self.uh.value(bc)-self.pde.solution(points)
        g = self.uh.grad_value(bc)-self.pde.gradient(points)
        measure = self.mesh.entity_measure("cell")
        l2 = np.sqrt(np.einsum('cqr,q,c->r', e**2, weights, measure))
        h1 = np.sqrt(np.einsum('cqrd,q,c->r', g**2, weights, measure))
        return dict(L2_components=l2.tolist(), H1_components=h1.tolist(),
                    L2=float(np.linalg.norm(l2)), H1_semi=float(np.linalg.norm(h1)))
