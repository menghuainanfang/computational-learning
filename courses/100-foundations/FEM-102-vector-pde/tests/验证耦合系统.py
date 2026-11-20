"""核对反应块、分量顺序、解耦退化及1—3维P1/P2的空间误差。"""
from pathlib import Path
from importlib.metadata import version
import json
import numpy as np
from fealpy.backend import backend_manager as bm
from 收束求解器 import SymbolicCDR, CDRSolver, box_mesh


def main():
    bm.set_backend("numpy")
    checks = {}
    S = np.array([[2., 1.], [1., 2.]])/6
    R = np.block([[2*S, -S], [-S, 2*S]])
    assert np.allclose(R @ [1, 1, 2, 2], [0, 0, 1.5, 1.5])
    values = np.array([[10, 20], [11, 21], [12, 22]])
    flat = values.T.reshape(-1)
    assert np.array_equal(flat, [10, 11, 12, 20, 21, 22])
    assert np.array_equal(flat.reshape(2, 3).T, values)
    checks['手算反应块及排列往返'] = True
    # 用实际FEALPy装配核对两个非对角块，包含空间变系数。
    model = SymbolicCDR(1, source=[1, 0], boundary=[0, 0], components=2,
                        reaction=[[2, 'x'], [-.25, 3]])
    solver = CDRSolver(model, box_mesh(model, 1))
    A, _ = solver.assemble()
    a = A.to_scipy().toarray()  # 仅4×4教学矩阵用于逐项核对。
    expected = np.array([[1, 1], [1, 3]])/12
    points = solver.scalar_space.interpolation_points()[:, 0]
    order = np.argsort(points)
    assert np.allclose(a[:2, 2:][np.ix_(order, order)], expected)
    assert np.allclose(a[2:, :2][np.ix_(order, order)], -.25*S)
    checks['非对角块积分与行列方向'] = True
    # 对角系统与两个单分量模型分别求解应相同。
    common = dict(a='1+x', b=[1, 2])
    model = SymbolicCDR(2, source=[1, 'x'], boundary=[0, 1], components=2,
                        reaction=[[2, 0], [0, 3]], **common)
    solver = CDRSolver(model, box_mesh(model, 8))
    actual = solver.solve_steady()[:].reshape(2, -1)
    for r, (f, g, c) in enumerate([(1, 0, 2), ('x', 1, 3)]):
        scalar = SymbolicCDR(2, source=[f], boundary=[g], components=1, reaction=c, **common)
        scalar_solver = CDRSolver(scalar, box_mesh(scalar, 8))
        assert np.max(np.abs(actual[r]-scalar_solver.solve_steady()[:])) < 1e-11
    checks['解耦退化到独立标量问题'] = True
    polynomial = SymbolicCDR(1, exact=['x*x', '1+x'], a='1+x', b=[2],
                             reaction=[[2, -1], [-.5, 3]])
    x = np.array([[.2], [.7]])
    expected = np.stack((2*x[:, 0]**2-x[:, 0]-3, 4+3*x[:, 0]-.5*x[:, 0]**2), axis=-1)
    assert np.allclose(polynomial.source(x), expected)
    checks['源项独立手算'] = True
    rows = []
    for dim in (1, 2, 3):
        xyz = ('x', 'y', 'z')[:dim]
        w = '1+'+'+'.join(xyz)+'+'+'*'.join(f'sin(pi*{x})' for x in xyz)
        model = SymbolicCDR(dim, exact=[w, f'2*({w})+1'], a='1+x',
                            b=list(range(1, dim+1)), reaction=[[2, -.5], [-.25, 3]])
        for p in (1, 2):
            old = None
            for n in ((4, 8, 16) if dim < 3 else (2, 4, 8)):
                solver = CDRSolver(model, box_mesh(model, n), p=p)
                solver.solve_steady()
                row = dict(dim=dim, p=p, n=n, **solver.errors(), **solver.diagnostics)
                assert row['residual'] < 1e-9 and row['boundary_error'] < 1e-11
                if old:
                    row['L2_order'] = float(np.log2(old['L2']/row['L2']))
                    row['H1_order'] = float(np.log2(old['H1_semi']/row['H1_semi']))
                rows.append(row)
                old = row
            assert row['L2_order'] > p+.6 and row['H1_order'] > p-.35, row
    # 三分量、空间变反应矩阵、不同非零边界检验轴与边界排列。
    model = SymbolicCDR(2, exact=['1+x*x', '2+y*y', '3+x*y'],
                        reaction=[[3, '-x', .2], [-.5, 4, .3], [.1, -.2, 2]])
    solver = CDRSolver(model, box_mesh(model, 3), p=2)
    solver.solve_steady()
    assert solver.errors()['L2'] < 1e-11 and solver.diagnostics['boundary_error'] < 1e-11
    checks['三分量P2多项式与不同边界'] = True
    data = dict(FEALPy=version('fealpy'), NumPy=np.__version__, backend='numpy', checks=checks, spatial=rows)
    out = Path(__file__).resolve().parent/'验证结果.json'
    out.write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf8')
    print('通过：', checks, '；空间加密18组。')


if __name__ == '__main__':
    main()
