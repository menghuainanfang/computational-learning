"""固定半离散系统比较时间误差；固定条件比较黏性与输运功。"""
from pathlib import Path
import json, csv
from importlib.metadata import version
import numpy as np
from scipy.integrate import solve_ivp
from scipy.sparse.linalg import splu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from 非线性有限元 import BurgersFEM


def main():
    out = Path(__file__).with_name('results')
    figures = Path(__file__).resolve().parents[1]/'notes'/'figures'
    figures.mkdir(exist_ok=True)
    plt.rcParams['font.sans-serif']=['Microsoft YaHei','SimHei','DejaVu Sans']
    plt.rcParams['axes.unicode_minus']=False
    out.mkdir(exist_ok=True)
    fem = BurgersFEM(n=8, p=1, dim=2)
    initial = lambda x, t: np.prod(np.sin(np.pi*x), axis=-1)
    zero = lambda x, t: np.zeros(x.shape[:-1])
    U0 = initial(fem.nodes, 0)
    U0[fem.boundary] = 0
    I = fem.I
    mass = splu(fem.M[I][:, I].tocsc())
    def rhs(t, interior):
        U = np.zeros(fem.N)
        U[I] = interior
        return mass.solve((-.1*(fem.K@U)-fem.nonlinear(U)[0])[I])
    reference = solve_ivp(rhs, (0, .1), U0[I], method='DOP853', rtol=2e-11, atol=2e-13)
    assert reference.success
    ref = np.zeros(fem.N)
    ref[I] = reference.y[:, -1]
    rows = []
    for scheme in ['imex', 'front', 'gradient', 'newton', 'picard']:
        prev = None
        for dt in [.02, .01, .005]:
            U, history = fem.run(initial, zero, dt=dt, end=.1, scheme=scheme)
            row = dict(scheme=scheme, dt=dt, error=fem.norm(U-ref),
                residual=max(x['residual'] for x in history), full_residual=max(x['full_residual'] for x in history),
                iterations=sum(x['iterations'] for x in history), balance=max(abs(x['balance']) for x in history))
            if prev: row['order'] = float(np.log2(prev['error']/row['error']))
            assert row['residual'] < 1e-10 and row['balance'] < 1e-9
            rows.append(row)
            prev = row
        assert rows[-1]['order'] > .8
    assert abs(rows[11]['error']-rows[14]['error']) < 1e-9
    cases = []
    values = []
    for nu in [.03, .1, .3]:
        U, history = fem.run(initial, zero, dt=.005, end=.1, nu=nu)
        cases.append(dict(nu=nu, history=history))
        values.append(U)
    data = dict(FEALPy=version('fealpy'), NumPy=np.__version__, n=8, p=1, end=.1,
        reference=dict(method='DOP853', rtol=2e-11, atol=2e-13, nfev=reference.nfev), time=rows, viscosity=cases)
    (out/'时间与黏性结果.json').write_text(json.dumps(data, ensure_ascii=False, indent=2), encoding='utf-8')
    with (out/'时间误差.csv').open('w', encoding='utf-8-sig', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[-1])); writer.writeheader(); writer.writerows(rows)
    fig, ax = plt.subplots(figsize=(7, 4))
    for scheme in ['imex', 'front', 'gradient', 'newton', 'picard']:
        selected = [r for r in rows if r['scheme'] == scheme]
        ax.loglog([r['dt'] for r in selected], [r['error'] for r in selected], 'o-', label=scheme)
    ax.set(xlabel='时间步长', ylabel='相对半离散参考的质量范数误差')
    ax.legend(); ax.grid(True); fig.tight_layout(); fig.savefig(figures/'时间误差.png', dpi=160); plt.close(fig)
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.5), constrained_layout=True)
    for ax, U, nu in zip(axes, values, [.03, .1, .3]):
        image = ax.tripcolor(fem.nodes[:, 0], fem.nodes[:, 1], fem.mesh.entity('cell'), U, shading='gouraud', vmin=0, vmax=1)
        ax.set(title=f'nu={nu}, t=0.1', xlabel='x', ylabel='y', aspect='equal')
    fig.colorbar(image, ax=axes, label='u');fig.savefig(figures/'同一时刻黏性对照.png', dpi=160);plt.close(fig)
    print(json.dumps(rows, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
