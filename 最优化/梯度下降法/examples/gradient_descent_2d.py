"""Optimize a coupled quadratic in R^2 using bounded numerical search.

Input: [2, 3]. Output: final point and optional contour trajectory.
Example: python examples/gradient_descent_2d.py --save contour.png.
"""
import argparse
import autograd.numpy as np
from miniopt import gradient_descent
from miniopt.visualization import plot_trajectory_2d


def objective(w):
    a = np.array([[3., 1.], [1., 2.]])
    b = np.array([1., 2.])
    return .5 * np.dot(w, np.dot(a, w)) - np.dot(b, w)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--save")
    args = parser.parse_args()
    result = gradient_descent(objective, [2., 3.], step="bounded_search")
    print(f"x={result.x}, f={result.fun:.9f}, {result.message}")
    plot_trajectory_2d(objective, result, show=args.show, save=args.save)


if __name__ == "__main__":
    main()
