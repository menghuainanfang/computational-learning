"""共享标量扩散与给定速度，反应矩阵耦合的稳态方程数据。"""
import sympy as sp
from fealpy.backend import backend_manager as bm
from fealpy.decorator import cartesian


class CDRModel:
    """点数组末轴为空间轴；解、源项与边界另带物质轴。"""
    def __init__(self, dim, components=2, box=None):
        if dim not in (1, 2, 3) or not isinstance(components, int) or components < 1:
            raise ValueError("使用1—3维、正整数分量数")
        self.dim, self.components = dim, components
        self.box = list(box) if box is not None else [0., 1.] * dim
        if len(self.box) != 2*dim or any(self.box[2*i] >= self.box[2*i+1] for i in range(dim)):
            raise ValueError("区域按各方向左、右端点排列")

    def domain(self):
        return list(self.box)

    @cartesian
    def is_dirichlet_boundary(self, p):
        mask = bm.zeros(p.shape[:-1], dtype=bm.bool)
        for i in range(self.dim):
            mask |= bm.abs(p[..., i]-self.box[2*i]) < 1e-12
            mask |= bm.abs(p[..., i]-self.box[2*i+1]) < 1e-12
        return mask


def component_field(function, component):
    """从多分量数据中取标量分量，保留Cartesian坐标约定。"""
    @cartesian
    def evaluate(p):
        return function(p)[(...,)+component]
    return evaluate


class SymbolicCDR(CDRModel):
    """使用x、y、z表达稳态数据；exact/source/boundary均为分量列表。"""
    def __init__(self, dim, *, exact=None, source=None, boundary=None,
                 a=1, b=None, reaction=0, components=None, box=None):
        count = len(exact) if exact is not None else components
        if count is None:
            raise ValueError("无真解时需指定components")
        super().__init__(dim, count, box)
        coordinates = sp.symbols("x y z")[:dim]

        def expr(v):
            value = sp.sympify(v)
            if not value.free_symbols <= set(coordinates):
                raise ValueError("表达式仅含当前空间坐标")
            return value

        a = expr(a)
        b = [expr(v) for v in (b if b is not None else [0]*dim)]
        C = (sp.Matrix(reaction) if isinstance(reaction, (list, tuple, sp.MatrixBase))
             else sp.eye(count)*expr(reaction)).applyfunc(expr)
        if len(b) != dim or C.shape != (count, count):
            raise ValueError("速度或反应矩阵的维数不符")
        u = [expr(v) for v in exact] if exact is not None else None
        if source is None:
            if u is None:
                raise ValueError("需给出源项或制造解")
            source = [sum(-sp.diff(a*sp.diff(u[r], x), x)+b[j]*sp.diff(u[r], x)
                          for j, x in enumerate(coordinates))
                      +sum(C[r, s]*u[s] for s in range(count)) for r in range(count)]
        if boundary is None:
            boundary = u
        if boundary is None:
            raise ValueError("无真解时需给出边界值列表")
        self._functions = {}

        def compile_field(name, values, shape):
            self._functions[name] = ([sp.lambdify(coordinates, expr(v), "numpy") for v in values], shape)

        compile_field("a", [a], ())
        compile_field("b", b, (dim,))
        compile_field("C", list(C), (count, count))
        for name, values in (("f", source), ("g", boundary), ("u", u)):
            if values is not None:
                if len(values) != count:
                    raise ValueError("源项、边界与真解列表长度应等于分量数")
                compile_field(name, values, (count,))
        if u is not None:
            compile_field("gradient", [sp.diff(v, x) for v in u for x in coordinates], (count, dim))

    def _evaluate(self, name, p):
        if name not in self._functions:
            raise ValueError("本模型未给精确解，不能计算真解误差")
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
        return self._evaluate("C", p)

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
