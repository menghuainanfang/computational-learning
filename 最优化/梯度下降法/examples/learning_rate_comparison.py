"""Compare fixed, diminishing, and bounded steps on the original quartic.

Input: initial vector [2]. Output: reports and optional convergence comparison.
Example: python examples/learning_rate_comparison.py --show.
"""
import argparse
from miniopt import gradient_descent
from miniopt.visualization import compare_results


def objective(w):
    return (w[0]**4 + w[0]**2 + 10*w[0]) / 50


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--save")
    args = parser.parse_args()
    configurations = {
        "fixed alpha=1": ("fixed", 1.),
        "fixed alpha=0.1": ("fixed", .1),
        "diminishing alpha=1": ("diminishing", 1.),
        "bounded search": ("bounded_search", 1.),
    }
    results = {}
    for name, (step, alpha) in configurations.items():
        result = gradient_descent(objective, [2.], step=step, alpha=alpha)
        results[name] = result
        print(f"{name}: x={result.x}, f={result.fun:.9f}, {result.message}")
    compare_results(results, show=args.show, save=args.save)


if __name__ == "__main__":
    main()
