# MiniOpt 0.1.0 implementation report

MiniOpt provides reusable gradient descent for real scalar objectives in any
parameter dimension, with independent differentiation, step selection, results,
and plotting. The release passes 48 automated tests and can be installed both
as an editable checkout and from its wheel.

Repository: [computational-learning](https://github.com/menghuainanfang/computational-learning).
Project directory: `最优化/梯度下降法`.
License: MIT. Verification date: October 14, 2026.

## Project structure

```text
miniopt/
  pyproject.toml
  MANIFEST.in
  README.md
  LICENSE
  .gitignore
  src/miniopt/
    __init__.py
    differentiation.py
    step_size.py
    results.py
    visualization.py
    optimizers/
      __init__.py
      gradient_descent.py
  examples/
    gradient_descent_1d.py
    gradient_descent_2d.py
    quadratic_optimization.py
    learning_rate_comparison.py
  tests/
    test_differentiation.py
    test_step_size.py
    test_gradient_descent.py
    test_visualization.py
  docs/
    mathematical_background.md
    implementation_report.md
```

## Design decisions

- A nonempty one-dimensional NumPy vector represents all optimization
  parameters. Objectives return a real scalar; gradients match the vector.
- Autograd supplies derivatives. Analytical gradients can be passed directly
  for objectives using ordinary NumPy or external numerical code.
- Step selection uses functions: fixed, diminishing, bounded numerical search,
  or a custom callable. No strategy class hierarchy is needed.
- SciPy supplies bounded search. Its candidate and both interval endpoints
  are compared. Invalid numerical results, nonfinite evaluations, and zero
  steps stop the run with an explicit unsuccessful result.
- Result and history containers are dataclasses. All history arrays include
  the initial point and have `nit + 1` entries; the initial step entry is zero.
  Gradients are evaluated at the final accepted point, including when the
  iteration budget has been exhausted.
- `success=True` means approximate stationarity under the gradient tolerance.
  Configuration and interface errors raise exceptions; numerical failures
  return their termination reason with recorded history.
- Matplotlib is an optional dependency. Core solver imports work without it.
- Setuptools builds the package from a `src` layout. The source distribution
  includes examples, tests, documentation, and the license.

## Environment and dependencies

The project lives in `最优化/梯度下降法` within the course repository, with its
development environment in `.venv`. The original prototype is preserved.
No global Python configuration was changed.

| Component | Verified version | Role |
| --- | --- | --- |
| Python | 3.13.2 | Windows interpreter |
| NumPy | 2.5.4 | Vectors and numerical operations |
| Autograd | 1.9.1 | Automatic differentiation |
| SciPy | 1.18.1 | Bounded scalar minimization |
| Matplotlib | 3.11.2 | Optional figures |
| pytest | 9.1.1 | Automated tests |
| build | 1.6.1 | Distribution frontend |
| setuptools | 84.0.0 | Isolated build backend |

Git uses the existing user identity. The authenticated GitHub account is
`menghuainanfang`. Only package sources, examples, tests, and documentation
are committed; environments, build output, plots, and logs are ignored.

## Verification results

| Check | Result |
| --- | --- |
| Automated suite | 48 passed in 9.58 seconds |
| Dependency consistency | No broken requirements |
| Editable installation | Passed |
| Import outside the source directory | Passed |
| Wheel installation in a separate environment | Passed |
| Wheel-only sample optimization | Converged to [3, 3] |
| Core import without Matplotlib | Passed |
| Source archive and wheel build | Passed |
| Four runnable examples | Passed |
| Saved convergence plot | Rendered and visually inspected |

Tests cover one-dimensional and positive definite multidimensional quadratics,
the original quartic, analytical and central-difference gradient comparisons,
all step strategies, and the callable extension interface. Failure cases cover
invalid input dimensions, gradient shape mismatches, invalid rates and settings,
nonfinite objectives, derivatives and trial iterates, stationary initial points,
zero steps, failed searches, and consistent history including final gradients.

The original quartic with fixed alpha=1 reaches approximate stationarity in
43 updates at w=-1.2347728 and f(w)=-0.169969284. The five-dimensional quadratic
example has solution error 9.046e-9 relative to the analytical linear solve.

## Installation and usage

```bash
git clone https://github.com/menghuainanfang/computational-learning.git
cd computational-learning/最优化/梯度下降法
python -m venv .venv
# Windows PowerShell:
.venv\Scripts\Activate.ps1
python -m pip install -e ".[dev]"
python -m pytest
python examples/gradient_descent_1d.py --show
python -m build
```

```python
import autograd.numpy as np
from miniopt import gradient_descent

result = gradient_descent(lambda w: np.sum((w - 3)**2),
                          [2., -1.], step="fixed", alpha=.1)
print(result.x, result.fun, result.message)
```

## Limitations

Success establishes approximate stationarity, not a global minimum.
Fixed and diminishing schedules depend on objective scaling, while bounded
numerical search has no global optimality guarantee on a nonconvex line.
On the original quartic, diminishing alpha=1 reaches the 1,000-update limit
at w=-0.39532823 and f(w)=-0.075451461; reducing steps slows progress on this
example. This remains an unsuccessful termination in the reported comparison.

Objectives require compatible differentiable operations or a supplied gradient.
Constraints and complex variables are outside this release. Full history uses
memory proportional to updates times dimension. Trajectory figures support
one or two parameters. Verification was performed on Windows with Python
3.13.2; the declared older Python versions were not tested locally.

## Next milestone

Implement Armijo backtracking through the callable step interface and verify
its sufficient-decrease condition. Evaluation counts and optional history
sampling would then support experiments on larger problems.
