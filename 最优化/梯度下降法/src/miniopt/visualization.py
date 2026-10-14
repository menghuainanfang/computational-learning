"""Plot scalar convergence measures and low-dimensional iterate paths.

Inputs are OptimizationResult objects and, for paths, f: R^n -> R.
Functions return (figure, axes); show and save are optional.
Example: plot_convergence(result, save='convergence.png').
"""
import numpy as np
import matplotlib.pyplot as plt
from .differentiation import objective_value


def _finish(fig, show, save):
    fig.tight_layout()
    if save is not None:
        fig.savefig(save, dpi=160)
    if show:
        plt.show()


def plot_convergence(result, *, show=False, save=None):
    """Plot objective values and gradient norms, in any dimension."""
    return compare_results({"run": result}, show=show, save=save)


def compare_results(results, *, show=False, save=None):
    """Compare named runs, for example runs with different learning rates."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for name, result in results.items():
        axes[0].plot(result.history.fun, label=name)
        axes[1].semilogy(np.maximum(result.history.gradient_norm,
                                    np.finfo(float).tiny), label=name)
    axes[0].set(xlabel="Iteration", ylabel="f(w)", title="Objective value")
    axes[1].set(xlabel="Iteration", ylabel="Gradient norm", title="Stationarity")
    for ax in axes:
        ax.grid(alpha=0.3)
        ax.legend()
    _finish(fig, show, save)
    return fig, axes


def plot_trajectory_1d(f, result, *, bounds=None, show=False, save=None):
    """Overlay recorded iterates on a one-dimensional objective curve."""
    if result.history.x.shape[1] != 1:
        raise ValueError("A 1-D trajectory requires one parameter")
    x = result.history.x[:, 0]
    if bounds is None:
        margin = max(float(np.ptp(x)) * 0.15, 0.5)
        bounds = (float(x.min() - margin), float(x.max() + margin))
    grid = np.linspace(*bounds, 300)
    fig, ax = plt.subplots()
    ax.plot(grid, [objective_value(f, np.array([v])) for v in grid], label="f(w)")
    ax.plot(x, result.history.fun, "o-", markersize=3, label="Iterates")
    ax.set(xlabel="w[0]", ylabel="f(w)")
    ax.legend()
    _finish(fig, show, save)
    return fig, ax


def plot_trajectory_2d(f, result, *, bounds=None, show=False, save=None):
    """Overlay a two-dimensional trajectory on objective contours.

    bounds is ((xmin, xmax), (ymin, ymax)); defaults cover the trajectory.
    """
    points = result.history.x
    if points.shape[1] != 2:
        raise ValueError("A 2-D trajectory requires two parameters")
    if bounds is None:
        margin = np.maximum(np.ptp(points, axis=0) * 0.15, 0.5)
        bounds = list(zip(points.min(axis=0) - margin, points.max(axis=0) + margin))
    xx, yy = np.meshgrid(np.linspace(*bounds[0], 100),
                         np.linspace(*bounds[1], 100))
    zz = np.array([objective_value(f, p) for p in
                   np.column_stack((xx.ravel(), yy.ravel()))]).reshape(xx.shape)
    fig, ax = plt.subplots()
    contours = ax.contour(xx, yy, zz, levels=20)
    ax.clabel(contours, fontsize=8)
    ax.plot(points[:, 0], points[:, 1], "o-", markersize=3)
    ax.set(xlabel="w[0]", ylabel="w[1]")
    _finish(fig, show, save)
    return fig, ax
