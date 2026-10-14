"""Check fixed, diminishing, and bounded steps along descent directions.

Tests compare numerical steps with known scalar minima and inject search
failures. Example: python -m pytest tests/test_step_size.py.
"""
from types import SimpleNamespace
import numpy as np
import pytest
import miniopt.step_size as module
from miniopt.step_size import choose_step, StepSizeError


def test_fixed_and_diminishing():
    f = lambda w: w[0]**2
    assert choose_step("fixed", f, np.array([1.]), np.array([-2.]), 4, alpha=.2) == .2
    assert choose_step("diminishing", f, np.array([1.]), np.array([-2.]), 4, alpha=.2) == .04


def test_bounded_interior_and_endpoint():
    f = lambda w: w[0]**2
    args = (f, np.array([1.]), np.array([-2.]), 0)
    assert choose_step("bounded_search", *args, alpha_max=2.) == pytest.approx(.5)
    assert choose_step("bounded_search", *args, alpha_max=.1) == .1


def test_zero_step():
    assert choose_step("bounded_search", lambda w: w[0]**2,
                       np.array([0.]), np.array([1.]), 0) == 0


@pytest.mark.parametrize("answer", [
    SimpleNamespace(success=False, x=.5, fun=0.),
    SimpleNamespace(success=True, x=np.nan, fun=0.),
    SimpleNamespace(success=True, x=.5, fun=np.inf),
    SimpleNamespace(success=True, x=7., fun=0.),
])
def test_invalid_search(monkeypatch, answer):
    monkeypatch.setattr(module, "minimize_scalar", lambda *a, **k: answer)
    with pytest.raises(StepSizeError):
        choose_step("bounded_search", lambda w: w[0]**2,
                    np.array([1.]), np.array([-2.]), 0)


def test_nonfinite_search_objective():
    with pytest.raises(StepSizeError, match="nonfinite"):
        choose_step("bounded_search", lambda w: np.nan,
                    np.array([1.]), np.array([-2.]), 0)


def test_search_exception(monkeypatch):
    def failed(*args, **kwargs):
        raise RuntimeError("Numerical failure")
    monkeypatch.setattr(module, "minimize_scalar", failed)
    with pytest.raises(StepSizeError, match="Numerical failure"):
        choose_step("bounded_search", lambda w: w[0]**2,
                    np.array([1.]), np.array([-2.]), 0)


@pytest.mark.parametrize("value", [0, -1, np.nan, np.inf])
def test_invalid_learning_rates(value):
    with pytest.raises(ValueError):
        choose_step("fixed", lambda w: w[0]**2,
                    np.array([1.]), np.array([-2.]), 0, alpha=value)
    with pytest.raises(ValueError):
        choose_step("bounded_search", lambda w: w[0]**2,
                    np.array([1.]), np.array([-2.]), 0, alpha_max=value)
