"""核对回溯和扩张中点的接受规则，以及完整次梯度实例。

输入为已知解析行为的二次函数和绝对值函数，输出由断言核对。
运行：python -m pytest tests/test_line_search.py。
"""
import numpy as np
import pytest

from 模块 import armijo_backtracking, expand_midpoint_search, gradient_descent
from 模块.回溯搜索 import StepSizeError


def squared(w):
    return float(np.dot(w, w))


def test_armijo_rejects_alpha_one_and_accepts_half():
    w = np.array([2.0]); g = 2*w; d = -g
    a = armijo_backtracking(squared, w, g, d, alpha0=1.0)
    assert a == 0.5


def test_expand_midpoint_uses_only_one_refinement():
    # From alpha=0.1, accepted increases 0.1->0.2->0.4.
    # 0.8 is worse than 0.4, midpoint 0.6 is also worse than 0.4.
    w = np.array([2.0]); g = 2*w
    probes = []
    def tracked(w_new):
        a = (2.0 - float(w_new[0])) / 4.0
        probes.append(a)
        return squared(w_new)
    a = expand_midpoint_search(tracked, w, g, -g, alpha0=0.1)
    assert a == pytest.approx(0.4)
    assert sum(abs(x-0.6) < 1e-9 for x in probes) == 1
    assert len(probes) == 6  # initial f, 0.1, 0.2, 0.4, 0.8, midpoint


def test_midpoint_can_improve_result():
    # start=0.2: 0.4 accepted, 0.8 worse, midpoint=0.6 worse => retain 0.4
    w = np.array([2.0]); g = 2*w
    assert expand_midpoint_search(squared, w, g, -g, alpha0=0.2) == pytest.approx(0.4)
    # start=0.3: 0.6 improves; 1.2 Armijo fails; midpoint 0.9 worse => 0.6
    assert expand_midpoint_search(squared, w, g, -g, alpha0=0.3, alpha_max=1.5) == pytest.approx(0.6)
    # start=0.15: 0.3 improves, 0.6 improves, 1.2 fails, midpoint 0.9 worse
    assert expand_midpoint_search(squared, w, g, -g, alpha0=0.15, alpha_max=2.0) == pytest.approx(0.6)


def test_one_midpoint_acceptance():
    # phi(alpha) = (2-4*alpha)**2, minimum at alpha=.5.
    # .18 -> .36 (better); .72 (worse); midpoint .54 (better than .36).
    w = np.array([2.0]); g = 2*w
    accepted = expand_midpoint_search(squared, w, g, -g, alpha0=0.18)
    assert accepted == pytest.approx(0.54)


def test_invalid_direction():
    w=np.array([1.0]); g=2*w
    with pytest.raises(ValueError, match="descent"):
        armijo_backtracking(squared,w,g,g)


def test_backtrack_when_initial_step_bad():
    w=np.array([2.0]);g=2*w
    got=expand_midpoint_search(squared,w,g,-g,alpha0=2.0,alpha_max=3.0)
    assert got > 0
    assert squared(w-got*g) < squared(w)


def test_solver_accepts_new_modes_without_automatic_differentiation():
    def grad_f(w):return 2*w
    for mode in ("armijo", "expand_midpoint"):
        result=gradient_descent(squared,[2.0],gradient=grad_f,step=mode,alpha=.1,max_iter=160)
        assert result.fun < 1e-10
        assert len(result.history.fun) == result.nit+1
        assert len(result.history.step_size) == result.nit+1


def test_no_acceptance_when_trial_budget_zero():
    w=np.array([2.0]);g=2*w
    with pytest.raises(StepSizeError):
        armijo_backtracking(squared,w,g,-g,alpha0=10.0,max_backtracks=0)


def test_exercise_3_6_original_and_perturbed_starts():
    from 模块.次梯度下降 import subgradient_descent
    from 实例.绝对值函数 import objective, subgradient
    def run(w0, step, alpha):
        return subgradient_descent(objective, [w0], subgradient=subgradient,
                                   step=step, alpha=alpha, max_iter=20)
    original = run(2.0, "fixed", 0.5)
    assert original.x[4, 0] == pytest.approx(0.0)
    assert original.fun[-1] == pytest.approx(0.0)
    assert original.x.shape == (21, 1)
    perturbed = run(2.1, "fixed", 0.5)
    assert perturbed.x[-1, 0] == pytest.approx(0.1)
    assert perturbed.x[-2, 0] == pytest.approx(-0.4)
    diminishing = run(2.1, "diminishing", 1.0)
    assert diminishing.fun[-1] < perturbed.fun[-1]
    np.testing.assert_allclose(diminishing.step_size[1:], 1 / np.arange(1., 21.))


def test_expansion_reaches_upper_bound_without_midpoint():
    probes = []
    def tracked(w):
        probes.append((2.0 - w[0]) / 4.0)
        return squared(w)
    w = np.array([2.0])
    assert expand_midpoint_search(tracked, w, 2*w, -2*w,
                                  alpha0=.1, alpha_max=.3) == pytest.approx(.3)
    np.testing.assert_allclose(probes, [0, .1, .2, .3])


def test_nonfinite_trial_is_backtracked():
    def positive_domain(w):
        return w[0]**2 if w[0] > 0 else np.inf
    w = np.array([1.0])
    assert armijo_backtracking(positive_domain, w, 2*w, -2*w) == .25
    assert expand_midpoint_search(positive_domain, w, 2*w, -2*w,
                                  alpha0=1.) == .375


@pytest.mark.parametrize("step", ["armijo", "expand_midpoint"])
def test_search_steps_meet_armijo_condition(step):
    result = gradient_descent(squared, [2.0, -1.0], gradient=lambda w: 2*w,
                              step=step, alpha=.1, max_iter=30)
    history = result.history
    for i in range(result.nit):
        g = 2 * history.x[i]
        threshold = history.fun[i] - 1e-4 * history.step_size[i + 1] * np.dot(g, g)
        assert history.fun[i + 1] <= threshold


@pytest.mark.parametrize("automatic", [False, True])
def test_evaluation_counts_match_actual_calls(automatic):
    calls = []
    def objective(w):
        calls.append(1)
        return w[0]**2
    gradient = None if automatic else lambda w: 2*w
    result = gradient_descent(objective, [2.0], gradient=gradient,
                              step="armijo", alpha=1., max_iter=5)
    assert result.success
    assert result.function_evaluations == len(calls)
    assert result.gradient_evaluations == result.nit + 1
    history = result.history
    assert history.function_evaluations.shape == (result.nit + 1,)
    assert history.function_evaluations[-1] == result.function_evaluations
    assert np.all(np.diff(history.function_evaluations) > 0)


def test_failed_search_counts_trials_and_keeps_last_point():
    result = gradient_descent(squared, [2.0], gradient=lambda w: 2*w,
                              step="armijo", alpha=10.,
                              step_options={"max_backtracks": 0})
    assert not result.success and result.nit == 0
    assert "trial budget" in result.message
    assert result.function_evaluations == 3
    assert result.history.function_evaluations[0] == 1
    np.testing.assert_array_equal(result.x, [2.0])


@pytest.mark.parametrize("options", [dict(growth=1), dict(max_expansions=-1),
                                       dict(shrink=1), dict(c=0)])
def test_invalid_expansion_configuration(options):
    w = np.array([2.0])
    with pytest.raises(ValueError):
        expand_midpoint_search(squared, w, 2*w, -2*w, **options)


def test_subgradient_vector_history():
    from 模块.次梯度下降 import subgradient_descent
    history = subgradient_descent(lambda w: np.sum(np.abs(w)), [2.1, -2.1],
                                  subgradient=np.sign, max_iter=20)
    assert history.x.shape == (21, 2)
    assert history.fun.shape == history.gradient_norm.shape == history.step_size.shape == (21,)
    np.testing.assert_allclose(history.x[-1], [.1, -.1])
    np.testing.assert_array_equal(history.function_evaluations, np.arange(1, 22))
