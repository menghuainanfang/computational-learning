"""稳态标量 CDR 的模型、装配与求解。"""
from .模型 import CDRModel, SymbolicCDR
from .求解器 import CDRSolver, box_mesh

__all__ = ["CDRModel", "SymbolicCDR", "CDRSolver", "box_mesh"]
