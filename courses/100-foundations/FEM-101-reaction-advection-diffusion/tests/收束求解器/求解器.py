"""标量拉格朗日有限元：空间装配、边界消元、稀疏求解、积分误差。"""
import numpy as np
from fealpy.backend import backend_manager as bm
from fealpy.mesh import IntervalMesh, TriangleMesh, TetrahedronMesh
from fealpy.functionspace import LagrangeFESpace
from fealpy.fem import (BilinearForm, LinearForm, DirichletBC,
                       ScalarDiffusionIntegrator, ScalarConvectionIntegrator,
                       ScalarMassIntegrator, ScalarSourceIntegrator)
from fealpy.solver import spsolve


def box_mesh(pde, n):
    """盒区域的单纯形网格：线段、三角形或四面体。"""
    if not isinstance(n, int) or n < 1:
        raise ValueError("每方向剖分段数n必须为正整数")
    if pde.dim == 1:
        node = bm.linspace(*pde.domain(), n+1).reshape(-1, 1)
        cell = bm.stack((bm.arange(n), bm.arange(1, n+1)), axis=-1)
        return IntervalMesh(node, cell)
    if pde.dim == 2:
        return TriangleMesh.from_box(pde.domain(), nx=n, ny=n)
    return TetrahedronMesh.from_box(pde.domain(), nx=n, ny=n, nz=n)


class CDRSolver:
    """NumPy/CPU，标准Galerkin，全维单纯形网格上的标量稳态系统。"""
    def __init__(self, pde, mesh, p=1, *, q=None, linear_solver=None):
        if bm.get_current_backend().backend_name != "numpy":
            raise RuntimeError("当前求解与验证使用NumPy/CPU")
        if not isinstance(p, int) or p < 1:
            raise ValueError("拉格朗日空间次数p应为正整数")
        if mesh.geo_dimension() != pde.dim or mesh.top_dimension() != pde.dim:
            raise ValueError("模型与全维网格的维数应相同")
        self.pde, self.mesh = pde, mesh
        self.space = LagrangeFESpace(mesh, p=p)
        self.q = p+4 if q is None else q
        self.linear_solver = linear_solver
        self.boundary = pde.is_dirichlet_boundary(self.space.interpolation_points())

    def assemble(self):
        """对应弱形式三项K、B、R及载荷F，返回尚未施加边界的A、F。"""
        bc, _ = self.mesh.quadrature_formula(self.q).get_quadrature_points_and_weights()
        a = self.pde.diffusion_coef(self.mesh.bc_to_point(bc))
        if not bm.all(bm.isfinite(a)) or bm.any(a <= 0):
            raise ValueError("本课扩散系数应在积分点有限且为正")
        form = BilinearForm(self.space)
        form.add_integrator(ScalarDiffusionIntegrator(coef=self.pde.diffusion_coef, q=self.q))
        form.add_integrator(ScalarConvectionIntegrator(coef=self.pde.convection_coef, q=self.q))
        form.add_integrator(ScalarMassIntegrator(coef=self.pde.reaction_coef, q=self.q))
        load = LinearForm(self.space)
        load.add_integrator(ScalarSourceIntegrator(source=self.pde.source, q=self.q))
        return form.assembly(), load.assembly()

    def solve_steady(self):
        self.A, self.F = self.assemble()
        self.Ab, self.Fb = DirichletBC(self.space, gd=self.pde.dirichlet,
                                     threshold=self.boundary).apply(self.A, self.F)
        values = (spsolve(self.Ab.copy(), self.Fb.copy(), solver="scipy") if self.linear_solver is None
                  else self.linear_solver(self.Ab.copy(), self.Fb.copy()))
        if values.shape != self.Fb.shape or not bm.all(bm.isfinite(values)):
            raise ValueError("线性求解应返回有限、形状正确的系数向量")
        self.uh = self.space.function()
        self.uh[:] = values
        interior = ~self.boundary
        self.diagnostics = {
            "dofs": self.space.number_of_global_dofs(),
            "residual": float(bm.linalg.norm((self.A @ values-self.F)[interior]))
                        / max(1., float(bm.linalg.norm(self.F[interior]))),
            "boundary_error": float(bm.max(bm.abs(values[self.boundary]
                -self.pde.dirichlet(self.space.interpolation_points())[self.boundary])))}
        return self.uh

    def errors(self):
        """以单元求积计算L2误差、H1半范数误差，适用于所有支持维数。"""
        bc, weights = self.mesh.quadrature_formula(self.q+2).get_quadrature_points_and_weights()
        points = self.mesh.bc_to_point(bc)
        value_error = self.uh.value(bc)-self.pde.solution(points)
        grad_error = self.uh.grad_value(bc)-self.pde.gradient(points)
        measure = self.mesh.entity_measure("cell")
        l2 = np.sqrt(np.einsum('cq,q,c->', value_error**2, weights, measure))
        h1 = np.sqrt(np.einsum('cqd,q,c->', grad_error**2, weights, measure))
        return float(l2), float(h1)
