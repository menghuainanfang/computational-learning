"""独立手算、矩阵恒等式、雅可比、边界及空间加密验证。"""
from pathlib import Path
import sys, json
from importlib.metadata import version
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'labs'))
from 非线性有限元 import BurgersFEM, manufactured


def main():
    checks = {}
    fem = BurgersFEM(3, dim=1)
    U = np.zeros(fem.N)
    U[np.argmin(abs(fem.nodes[:, 0]-1/3))] = 1.
    I = fem.I
    _, _, G = fem.nonlinear(U)
    assert np.allclose(fem.M[I][:, I].toarray(), [[2/9, 1/18], [1/18, 2/9]])
    assert np.allclose(G[I][:, I].toarray(), [[0, -1/6], [-1/6, -1/3]])
    V, info = fem.step(U, .1, .1, np.zeros(fem.N), 'gradient')
    assert np.allclose(V[I], [25/32, 25/128])
    checks['一维单步完整手算'] = True
    for dim in [1, 2]:
        for p in [1, 2]:
            f = BurgersFEM(4, p, dim)
            rng = np.random.default_rng(12)
            u, direction = rng.normal(size=(2, f.N))
            C, B, G = f.nonlinear(u)
            assert np.linalg.norm(C-B@u) < 1e-12
            assert np.linalg.norm(C-G@u) < 1e-12
            eps = 1e-6
            fd = (f.nonlinear(u+eps*direction)[0]-f.nonlinear(u-eps*direction)[0])/(2*eps)
            assert np.linalg.norm(fd-(B+G)@direction)/(1+np.linalg.norm(fd)) < 1e-8
    checks['一二维P1P2矩阵恒等式与雅可比'] = True
    rows = []
    # 固定终止时间，dt随h^(p+1)缩小，使一阶时间误差不支配空间观测。
    for p in [1, 2]:
        previous = None
        for n in [4, 8, 16]:
            f = BurgersFEM(n, p, 2)
            exact, source = manufactured(2, .1)
            dt = .02/(n//4)**(p+1)
            u, history = f.run(exact, source, dt=dt, end=.02)
            row = dict(p=p, n=n, dt=dt, L2=f.error(u, exact, .02),
                residual=max(x['residual'] for x in history), boundary=max(x['boundary'] for x in history))
            if previous: row['order'] = float(np.log2(previous['L2']/row['L2']))
            assert row['residual'] < 1e-10 and row['boundary'] == 0
            rows.append(row)
            previous = row
        assert rows[-1]['order'] > p+.5
    checks['二维制造解空间加密'] = True
    output = dict(FEALPy=version('fealpy'), NumPy=np.__version__, checks=checks, spatial=rows)
    Path(__file__).with_name('验证结果.json').write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(output, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
