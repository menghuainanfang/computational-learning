"""绘图：将优化结果画成收敛曲线或参数轨迹。

输入：结果对象 result；轨迹图还需要目标函数 f。
输出：(figure, axes)，可以显示，也可以保存到指定路径。
函数值曲线看目标函数的变化，梯度长度曲线看是否接近驻点。
例：plot_convergence(result, show=True)
"""
import numpy as np
import matplotlib.pyplot as plt
from .自动微分 import objective_value


def _finish(fig, show, save):
    fig.tight_layout()
    if save is not None:
        fig.savefig(save, dpi=160)
    if show:
        plt.show()


def plot_convergence(result, *, show=False, save=None):
    """画函数值与梯度长度随迭代次数变化的曲线，适用于任意维数。"""
    return compare_results({"run": result}, show=show, save=save)


def compare_results(results, *, show=False, save=None):
    """把多次运行的收敛曲线画在一起，例如比较不同步长。"""
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
    """画一维目标函数，并在曲线上标出迭代经过的点。"""
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
    """画二维等高线和迭代轨迹；bounds 为 ((xmin, xmax), (ymin, ymax))。"""
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


def compare_line_search(results, *, show=False, save=None):
    """比较二维线搜索：函数值、函数调用量、步长以及参数轨迹。

    输入为 {名称: OptimizationResult}。函数值采用对称对数坐标，
    既能显示接近零的变化，也能保留零值和负值。
    """
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for name, result in results.items():
        h = result.history
        if h.x.shape[1] != 2:
            plt.close(fig)
            raise ValueError("Line-search trajectory comparison requires two parameters")
        axes[0, 0].plot(h.fun, label=name)
        axes[0, 1].plot(h.function_evaluations, h.fun, label=name)
        axes[1, 0].plot(np.arange(1, result.nit + 1), h.step_size[1:], label=name)
        axes[1, 1].plot(h.x[:, 0], h.x[:, 1], "o-", markersize=2.5, label=name)
    axes[0, 0].set(title="Objective vs iteration", xlabel="Iteration", ylabel="f(w)")
    axes[0, 1].set(title="Objective vs function calls", xlabel="Function calls", ylabel="f(w)")
    for ax in axes[0]:
        ax.set_yscale("symlog", linthresh=1e-12)
    axes[1, 0].set(title="Accepted steps", xlabel="Update", ylabel="Step size")
    axes[1, 1].set(title="Parameter trajectories", xlabel="w[0]", ylabel="w[1]")
    for ax in axes.flat:
        ax.grid(alpha=0.25)
        ax.legend()
    _finish(fig, show, save)
    return fig, axes


def compare_subgradient_runs(cases, *, show=False, save=None):
    """画一维次梯度实验，每个初始点一行，比较函数值和参数变化。

    cases 的形式为 {初始点说明: {策略名称: IterationHistory}}。
    """
    fig, axes = plt.subplots(len(cases), 2, figsize=(11, 4 * len(cases)), squeeze=False)
    for row, (title, runs) in enumerate(cases.items()):
        for name, history in runs.items():
            if history.x.shape[1] != 1:
                plt.close(fig)
                raise ValueError("Subgradient trajectory comparison requires one parameter")
            axes[row, 0].plot(history.fun, "o-", markersize=3, label=name)
            axes[row, 1].plot(history.x[:, 0], "o-", markersize=3, label=name)
        axes[row, 0].set(title=f"{title}: objective", ylabel="f(w)")
        axes[row, 1].set(title=f"{title}: iterates", ylabel="w[0]")
        for ax in axes[row]:
            ax.axhline(0, color="gray", linewidth=0.7)
            ax.grid(alpha=0.25)
            ax.legend()
    for ax in axes[-1]:
        ax.set_xlabel("Iteration")
    _finish(fig, show, save)
    return fig, axes
