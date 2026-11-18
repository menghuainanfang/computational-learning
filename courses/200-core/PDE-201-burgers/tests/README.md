# 验证

从课程根目录执行 `python tests/验证非线性.py`：一维独立手算；一、二维P1/P2的 $\boldsymbol{C}(\boldsymbol{U})=\boldsymbol{B}(\boldsymbol{U})\boldsymbol{U}=\boldsymbol{G}(\boldsymbol{U})\boldsymbol{U}$ 及中心差商Jacobian；六组制造解联合加密。具体结果见[验证结果](验证结果.json)。

时间实验另运行 `python labs/运行时间实验.py`，使用固定半离散DOP853参考及三档时间步；要求实际方程残差小于1e-10、能量平衡缺口小于1e-9、末级观测阶大于0.8。Picard与Newton的终值时间误差差小于1e-9。

这些判据针对当前参数与光滑算例。代码核对、时间阶、联合误差与能量解释各有对象，不能互相替代。
