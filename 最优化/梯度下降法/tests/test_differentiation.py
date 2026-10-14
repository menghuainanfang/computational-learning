"""Check gradient vectors against analytical and numerical derivatives.

Run with pytest; scalar objectives map finite vectors to matching gradients.
Example: python -m pytest tests/test_differentiation.py.
"""
import autograd.numpy as anp
import numpy as np
import pytest
from 模块 import make_gradient, finite_difference_gradient


def test_matrix_gradient():
    a = anp.array([[3., 1.], [1., 2.]])
    b = anp.array([1., 2.])
    f = lambda w: 0.5 * anp.dot(w, anp.dot(a, w)) - anp.dot(b, w)
    w = np.array([0.7, -1.2])
    np.testing.assert_allclose(make_gradient(f)(w), a @ w - b)
    np.testing.assert_allclose(make_gradient(f)(w),
                               finite_difference_gradient(f, w), atol=1e-8)


def test_failed_differentiation():
    with pytest.raises(ValueError, match="Automatic differentiation failed"):
        make_gradient(lambda w: float(w[0]) ** 2)([1.])


def test_vector_objective():
    with pytest.raises(ValueError, match="real scalar"):
        make_gradient(lambda w: w ** 2)([1., 2.])


@pytest.mark.parametrize("epsilon", [0, -1, np.nan, np.inf])
def test_invalid_difference_step(epsilon):
    with pytest.raises(ValueError, match="epsilon"):
        finite_difference_gradient(lambda w: w[0]**2, [1.], epsilon)
