"""Minimize the original quartic f(w) = (w^4 + w^2 + 10w)/50.

Input: a vector with one parameter. Output: stationarity report and optional
trajectory figure. Example: python examples/gradient_descent_1d.py --show.
"""
import argparse
from miniopt import gradient_descent
from miniopt.visualization import plot_trajectory_1d


def objective(w):
    return (w[0]**4 + w[0]**2 + 10*w[0]) / 50


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--save")
    args = parser.parse_args()
    result = gradient_descent(objective, [2.], alpha=1.)
    print(f"x={result.x}, f={result.fun:.9f}, updates={result.nit}, {result.message}")
    plot_trajectory_1d(objective, result, show=args.show, save=args.save)


if __name__ == "__main__":
    main()
