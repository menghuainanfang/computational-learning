"""课程收束：跨维数、跨次数、变扩散系数及无真解输入验证。"""
from pathlib import Path
from importlib.metadata import version
import json
import numpy as np
from fealpy.backend import backend_manager as bm
from 收束求解器 import SymbolicCDR, CDRSolver, box_mesh


def main():
    bm.set_backend("numpy")
    rows = []
    for dim in (1, 2, 3):
        xyz = ('x', 'y', 'z')[:dim]
        # 非零边界、变系数扩散和非零对流/反应共同检验装配。
        exact = '1+'+'*'.join(f'sin(pi*{x})' for x in xyz)
        model = SymbolicCDR(dim, exact=exact, a='1+x', b=list(range(1, dim+1)), reaction=2)
        for p in (1, 2):
            previous = None
            for n in ((4, 8, 16) if dim < 3 else (2, 4, 8)):
                solver = CDRSolver(model, box_mesh(model, n), p=p)
                solver.solve_steady()
                l2, h1 = solver.errors()
                row = dict(dim=dim, p=p, n=n, L2=l2, H1_semi=h1, **solver.diagnostics)
                assert row['residual'] < 1e-9 and row['boundary_error'] < 1e-11, row
                if previous:
                    row['L2_order'] = float(np.log2(previous['L2']/l2))
                    row['H1_order'] = float(np.log2(previous['H1_semi']/h1))
                rows.append(row)
                previous = row
            assert row['L2_order'] > p+.65 and row['H1_order'] > p-.35, row
    # 同一物理问题不提供真解，也应得到同一个离散解。
    model = SymbolicCDR(1, exact='1+x*(1-x)')
    physical = SymbolicCDR(1, source=2, boundary=1)
    solvers = [CDRSolver(m, box_mesh(m, 8)) for m in (model, physical)]
    values = [s.solve_steady()[:] for s in solvers]
    assert np.max(np.abs(values[0]-values[1])) < 1e-12
    try:
        solvers[1].errors()
    except ValueError:
        pass
    else:
        raise AssertionError("无真解模型不应报告真解误差")
    # 独立可手算数据：-( (1+x)u' )' + 2u' + 3u，u=x²。
    polynomial = SymbolicCDR(1, exact='x**2', a='1+x', b=[2], reaction=3)
    points = np.array([[.2], [.7]])
    assert np.allclose(polynomial.source(points), -2+3*points[:, 0]**2)
    data = dict(FEALPy=version('fealpy'), NumPy=np.__version__, backend='numpy',
                cases=rows, 无真解输入一致性='通过', 变系数源项手算核对='通过')
    out = Path(__file__).resolve().parent/'收束验证结果.json'
    out.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf8')
    print(json.dumps(data, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
