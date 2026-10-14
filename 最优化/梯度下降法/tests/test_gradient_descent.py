"""Verify solutions, termination, and aligned records of gradient descent.

Quadratic minima are known analytically; the quartic tests the original
example. Example: python -m pytest tests/test_gradient_descent.py.
"""
import autograd.numpy as anp
import numpy as np
import pytest
from 模块 import gradient_descent
from 模块.步长策略 import StepSizeError


@pytest.mark.parametrize("step", ["fixed", "diminishing", "bounded_search"])
def test_1d_quadratic(step):
    result = gradient_descent(lambda w: (w[0]-3)**2, [0.], step=step,
                              alpha=.5, tol=1e-9)
    assert result.success
    np.testing.assert_allclose(result.x, [3.], atol=1e-8)


def test_multidimensional_quadratic():
    a = anp.array([[3., 1., 0.], [1., 4., 1.], [0., 1., 2.]])
    b = anp.array([1., 2., -1.])
    f = lambda w: .5 * anp.dot(w, anp.dot(a, w)) - anp.dot(b, w)
    result = gradient_descent(f, [2., -1., 3.], alpha=.1)
    assert result.success
    np.testing.assert_allclose(result.x, np.linalg.solve(a, b), atol=1e-8)
    h = result.history
    assert h.x.shape == (result.nit + 1, 3)
    assert all(v.shape == (result.nit + 1,) for v in
               (h.fun, h.gradient_norm, h.step_size))
    assert h.step_size[0] == 0
    np.testing.assert_allclose(h.fun, [f(w) for w in h.x])
    np.testing.assert_allclose(h.gradient_norm,
                               [np.linalg.norm(a @ w - b) for w in h.x])


def test_original_quartic():
    f = lambda w: (w[0]**4 + w[0]**2 + 10*w[0]) / 50
    result = gradient_descent(f, [2.], alpha=1.)
    optimum = np.roots([4., 0., 2., 10.])
    optimum = optimum[np.abs(optimum.imag) < 1e-12].real[0]
    assert result.success
    assert result.x[0] == pytest.approx(optimum, abs=1e-7)


@pytest.mark.parametrize("point", [1., [], [[1., 2.]], [np.nan], [np.inf], [1j]])
def test_invalid_point(point):
    with pytest.raises(ValueError):
        gradient_descent(lambda w: anp.sum(w**2), point)


def test_gradient_dimension():
    with pytest.raises(ValueError, match="shape"):
        gradient_descent(lambda w: anp.sum(w**2), [1., 2.], gradient=lambda w: [1.])


def test_initial_nonfinite_objective():
    with pytest.raises(ValueError, match="initial objective"):
        gradient_descent(lambda w: np.inf, [1.])


def test_nonfinite_gradient():
    result = gradient_descent(lambda w: w[0]**2, [1.], gradient=lambda w: [np.nan])
    assert not result.success and result.nit == 0
    assert "Nonfinite gradient" in result.message


def test_nonfinite_trial_objective():
    result = gradient_descent(lambda w: w[0]**2 if w[0] > 0 else np.inf,
                              [1.], alpha=1., gradient=lambda w: 2*w)
    assert not result.success and result.nit == 0
    np.testing.assert_array_equal(result.x, [1.])


def test_nonfinite_trial_iterate():
    with np.errstate(over="ignore"):
        result = gradient_descent(lambda w: w[0], [1e308], alpha=1e308,
                                  gradient=lambda w: [2.])
    assert not result.success and result.nit == 0
    assert "Nonfinite trial iterate" in result.message


def test_nonfinite_final_gradient_history():
    result = gradient_descent(lambda w: w[0]**2, [1.], alpha=.5,
                              gradient=lambda w: [2.] if w[0] > 0 else [np.nan])
    assert not result.success and result.nit == 1
    assert result.history.x.shape == (2, 1)
    assert np.isnan(result.history.gradient_norm[-1])


def test_diminishing_iterates():
    result = gradient_descent(lambda w: w[0]**2, [1.], step="diminishing",
                              alpha=.1, max_iter=20, tol=1e-14)
    steps = .1 / np.arange(1., 21.)
    expected = np.concatenate(([1.], np.cumprod(1 - 2 * steps)))
    np.testing.assert_allclose(result.history.x[:, 0], expected)
    np.testing.assert_allclose(result.history.step_size[1:], steps)


def test_stationary_and_final_iteration():
    f = lambda w: w[0]**2
    assert gradient_descent(f, [0.]).nit == 0
    assert gradient_descent(f, [0.], max_iter=0).success
    result = gradient_descent(f, [1.], alpha=.5, max_iter=1)
    assert result.success and result.nit == 1


def test_iteration_limit():
    result = gradient_descent(lambda w: w[0]**2, [1.], max_iter=1)
    assert not result.success and result.nit == 1
    assert result.history.gradient_norm.shape == (2,)


def test_custom_gradient_and_step():
    result = gradient_descent(lambda w: float(np.dot(w, w)), [1., 2.],
                              gradient=lambda w: 2*w, step=lambda f, w, d, k: .5)
    assert result.success and result.nit == 1


def test_search_failure_and_zero_step():
    def failed(*args):
        raise StepSizeError("Search failed")
    for step in (failed, lambda *args: 0., lambda *args: np.nan):
        result = gradient_descent(lambda w: w[0]**2, [1.], step=step)
        assert not result.success and result.nit == 0


@pytest.mark.parametrize("kwargs", [dict(max_iter=-1), dict(max_iter=1.5),
    dict(tol=0), dict(tol=np.nan), dict(step="unknown")])
def test_invalid_configuration(kwargs):
    with pytest.raises(ValueError):
        gradient_descent(lambda w: w[0]**2, [1.], **kwargs)
