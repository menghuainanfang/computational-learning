"""Containers for iterates of w_(k+1) = w_k - alpha_k grad f(w_k).

OptimizationResult holds the final point, objective, and termination status.
History arrays include the initial point; step_size[0] is zero. Every array
has nit + 1 records, and x has shape (nit + 1, n).
Example: result.history.fun plots objective values at recorded iterates.
"""
from dataclasses import dataclass
import numpy as np


@dataclass
class IterationHistory:
    x: np.ndarray
    fun: np.ndarray
    gradient_norm: np.ndarray
    step_size: np.ndarray


@dataclass
class OptimizationResult:
    x: np.ndarray
    fun: float
    nit: int
    success: bool
    message: str
    history: IterationHistory
