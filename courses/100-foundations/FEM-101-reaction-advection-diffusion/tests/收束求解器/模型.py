"""方程数据：-div(a grad u) + b.grad u + c u = f，全 Dirichlet。"""
import sympy as sp
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian


class CDRModel:
    """标量模型接口；标量场返回点数组前缀形状，向量场另加空间轴。"""
    def __init__(self, dim, box=None):
        if dim not in (1, 2, 3):
            raise ValueError("当前网格入口支持空间1、2、3维")
        self.dim = dim
        self.box = list(box) if box is not None else [0., 1.] * dim
        if len(self.box) != 2*dim or any(self.box[2*i] >= self.box[2*i+1] for i in range(dim)):
            raise ValueError("区域应按各方向的左、右端点给出")

    def domain(self):
        return list(self.box)

    @cartesian
    def is_dirichlet_boundary(self, p):
        mask = bm.zeros(p.shape[:-1], dtype=bm.bool)
        for axis in range(self.dim):
            mask |= bm.abs(p[..., axis] - self.box[2*axis]) < 1e-12
            mask |= bm.abs(p[..., axis] - self.box[2*axis+1]) < 1e-12
        return mask


class SymbolicCDR(CDRModel):
    """用x、y、z表达数据；可指定真解制造源项，也可直接给源项与边界。"""
    def __init__(self, dim, *, exact=None, source=None, boundary=None,
                 a=1, b=None, reaction=0, box=None):
        super().__init__(dim, box)
        coordinates = sp.symbols("x y z")[:dim]
        allowed = set(coordinates)

        def expr(v):
            value = sp.sympify(v)
            if not value.free_symbols <= allowed:
                raise ValueError("表达式只能含当前空间坐标，不包含时间或未赋值参数")
            return value

        diffusion = expr(a)
        velocity = [expr(v) for v in (b if b is not None else [0]*dim)]
        if len(velocity) != dim:
            raise ValueError("速度分量数应等于空间维数")
        reaction = expr(reaction)
        u = None if exact is None else expr(exact)
        if source is None:
            if u is None:
                raise ValueError("需要提供源项，或提供真解以生成源项")
            source = sum(-sp.diff(diffusion*sp.diff(u, x), x)
                         + velocity[i]*sp.diff(u, x) for i, x in enumerate(coordinates)) + reaction*u
        if boundary is None:
            if u is None:
                raise ValueError("无真解模型必须提供边界值")
            boundary = u
        self._functions = {}

        def compile_field(name, expressions, shape):
            self._functions[name] = ([sp.lambdify(coordinates, expr(v), "numpy")
                                     for v in expressions], shape)

        compile_field("a", [diffusion], ())
        compile_field("b", velocity, (dim,))
        compile_field("c", [reaction], ())
        compile_field("f", [source], ())
        compile_field("g", [boundary], ())
        if u is not None:
            compile_field("u", [u], ())
            compile_field("gradient", [sp.diff(u, x) for x in coordinates], (dim,))

    def _evaluate(self, name, p):
        if name not in self._functions:
            raise ValueError("当前模型未提供精确解或梯度，不能计算真解误差")
        functions, shape = self._functions[name]
        args = [p[..., i] for i in range(self.dim)]
        values = [bm.broadcast_to(bm.asarray(f(*args), dtype=p.dtype), p.shape[:-1]) for f in functions]
        return bm.stack(values, axis=-1).reshape(p.shape[:-1]+shape) if shape else values[0]

    @cartesian
    def diffusion_coef(self, p):
        return self._evaluate("a", p)

    @cartesian
    def convection_coef(self, p):
        return self._evaluate("b", p)

    @cartesian
    def reaction_coef(self, p):
        return self._evaluate("c", p)

    @cartesian
    def source(self, p):
        return self._evaluate("f", p)

    @cartesian
    def dirichlet(self, p):
        return self._evaluate("g", p)

    @cartesian
    def solution(self, p):
        return self._evaluate("u", p)

    @cartesian
    def gradient(self, p):
        return self._evaluate("gradient", p)
