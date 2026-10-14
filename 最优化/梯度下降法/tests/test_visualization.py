"""Exercise convergence and trajectory figures using a headless backend.

Results supply iterates and scalar histories; saved figures must be readable.
Example: python -m pytest tests/test_visualization.py.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest
from 模块 import gradient_descent
from 模块.绘图 import (
    plot_convergence, compare_results, plot_trajectory_1d, plot_trajectory_2d,
)


def test_figures(tmp_path):
    f = lambda w: w[0]**2
    result = gradient_descent(f, [1.])
    path = tmp_path / "plot.png"
    fig, _ = plot_convergence(result, save=path)
    assert plt.imread(path).size > 0
    plt.close(fig)
    fig, _ = compare_results({"a": result, "b": result})
    plt.close(fig)
    fig, _ = plot_trajectory_1d(f, result)
    plt.close(fig)
    f2 = lambda w: w[0]**2 + w[1]**2
    result2 = gradient_descent(f2, [1., 2.])
    fig, _ = plot_trajectory_2d(f2, result2)
    plt.close(fig)
    with pytest.raises(ValueError, match="one parameter"):
        plot_trajectory_1d(f2, result2)
    with pytest.raises(ValueError, match="two parameters"):
        plot_trajectory_2d(f, result)
