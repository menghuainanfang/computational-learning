"""Compare an R^5 quadratic minimizer with the solution of Aw=b.

Input: symmetric positive definite A and b. Output: final error and optional
convergence figure. Example: python examples/quadratic_optimization.py.
"""
import argparse
import autograd.numpy as np
from miniopt import gradient_descent
from miniopt.visualization import plot_convergence


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--save")
    args = parser.parse_args()
    a = np.diag(np.arange(1., 6.))
    b = np.ones(5)
    def objective(w):
        return .5 * np.dot(w, np.dot(a, w)) - np.dot(b, w)
    result = gradient_descent(objective, np.zeros(5), alpha=.2)
    error = np.linalg.norm(result.x - np.linalg.solve(a, b))
    print(f"x={result.x}, solution error={error:.3e}, {result.message}")
    plot_convergence(result, show=args.show, save=args.save)


if __name__ == "__main__":
    main()
