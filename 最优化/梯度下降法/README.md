# MiniOpt

MiniOpt is a small educational Python library for gradient descent on real,
unconstrained, differentiable scalar objectives. It grew from exercises in
*Machine Learning Refined: Foundations, Algorithms, and Applications*.
The reusable solver accepts any parameter dimension through a single vector
interface. Example objectives live outside the library.

## Installation

Python 3.10 or later is required. From the repository root, first run
`cd 最优化/梯度下降法`. Then install this package:

```bash
python -m venv .venv
# Windows PowerShell
.venv\Scripts\Activate.ps1
# macOS / Linux: source .venv/bin/activate
python -m pip install -e ".[dev]"
```

Use `python -m pip install -e .` for the core library, or `-e ".[plot]"`
for plotting. NumPy, Autograd, and SciPy are core dependencies. Matplotlib
is imported only by the visualization module.

## Minimal example

```python
import autograd.numpy as np
from miniopt import gradient_descent

def f(w):
    return np.sum((w - 3.0) ** 2)

result = gradient_descent(f, w0=[2.0, -1.0], alpha=0.1,
                          max_iter=1000, tol=1e-8)
print(result.x, result.fun, result.success, result.message)
```

An objective takes a nonempty one-dimensional real vector and returns a scalar.
Use `autograd.numpy` for operations on differentiable variables. Ordinary NumPy
operations may fail during tracing, and conversions such as `float(w[0])` break
the derivative. Constants and matrices can be ordinary NumPy arrays.
For an objective that uses ordinary NumPy or external code, supply an analytical
gradient: `gradient_descent(f, w0, gradient=grad_f)`. Its output must match `w0`.
`make_gradient(f)` can also be used independently; `finite_difference_gradient`
provides a central-difference comparison for small verification problems.

## Step sizes

| Strategy | Configuration | Definition |
| --- | --- | --- |
| Fixed | `step="fixed", alpha=0.1` | alpha at every update |
| Diminishing | `step="diminishing", alpha=0.1` | alpha / (k + 1), starting at k = 0 |
| Bounded numerical search | `step="bounded_search", alpha_max=6.0` | SciPy bounded minimization along the descent direction |

Bounded search compares the numerical interior candidate with both endpoints
of `[0, alpha_max]`. A custom callable `step(f, w, direction, k)` may return a
nonnegative scalar step, allowing another strategy without changing the solver.
A zero step stops the solver with `success=False`.

## Results and plots

`result.x`, `fun`, `nit`, `success`, and `message` describe the final recorded
iterate. `nit` counts accepted updates. `result.history` contains `x`, `fun`,
`gradient_norm`, and `step_size`, all with `nit + 1` records. Initial values are
included; `step_size[0] = 0`, and subsequent entries are the steps used to reach
those records. `history.x` has shape `(nit + 1, n)`.

Interface mistakes raise `ValueError`. Numerical search failures and nonfinite
trial values return an unsuccessful result at the last recorded iterate.
A nonfinite derivative at a recorded point is retained in the history and
reported as a failure. Exceptions raised by user functions retain their context.
`max_iter=0` evaluates the initial point without taking an update.

```python
from miniopt.visualization import plot_convergence, compare_results

fig, axes = plot_convergence(result, show=True, save="convergence.png")
# compare_results({"alpha=0.1": result_a, "alpha=0.2": result_b})
```

`plot_trajectory_1d(f, result)` plots iterates over an objective curve;
`plot_trajectory_2d(f, result)` plots iterates over contours. Each plotting
function returns its figure and axes, defaults to `show=False`, and accepts
an optional save path.

## Project structure

```text
src/miniopt/
  differentiation.py       Autograd and central differences
  step_size.py             Functional step strategies
  results.py               Dataclass results and aligned history
  visualization.py         Optional plotting
  optimizers/              Gradient descent
examples/                  Four runnable optimization experiments
tests/                     Analytical, regression, and failure tests
docs/                      Mathematical background and implementation report
```

Run examples from the checkout after installation:

```bash
python examples/gradient_descent_1d.py
python examples/gradient_descent_2d.py
python examples/quadratic_optimization.py
python examples/learning_rate_comparison.py
python -m pytest
python -m build
```

Examples accept `--show` to open figures and `--save PATH` to save them. Create
the destination directory first. The quadratic example prints its comparison
with the analytical solution.

## Mathematical background

Gradient descent updates `w_next = w - alpha * grad f(w)`. A small gradient
indicates approximate stationarity. See [mathematical background](docs/mathematical_background.md)
for descent conditions and the role of the step size.

## Limitations

- Approximate stationarity does not imply a global minimum; a stationary point
  may be a maximum or saddle point.
- Fixed and diminishing steps need suitable scales. They can increase the
  objective or diverge, and a diminishing schedule can make slow progress.
- Bounded numerical search is approximate and does not guarantee a globally
  optimal step on a nonconvex line. Nonfinite evaluations stop the search.
- Objectives must be differentiable and Autograd-compatible unless an
  analytical gradient is supplied. Constraints and complex variables are
  outside the current interface.
- Full iterate history uses memory proportional to dimension times updates;
  trajectory plots are limited to one or two parameters.
- This release is verified locally on Windows with Python 3.13.2. Other
  supported Python versions have not been exercised in this environment.

## Next milestone

Add Armijo backtracking through the existing callable step interface, with
tests for sufficient decrease. Then consider evaluation counts and optional
history sampling for larger problems.

## License

MIT. Copyright 2026 menghuainanfang.
