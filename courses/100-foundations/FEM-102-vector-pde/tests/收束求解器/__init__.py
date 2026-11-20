"""稳态多分量线性反应耦合的课程收束源码。"""
from .模型 import CDRModel, SymbolicCDR, component_field
from .求解器 import CDRSolver, box_mesh

__all__ = ["CDRModel", "SymbolicCDR", "CDRSolver", "box_mesh", "component_field"]
