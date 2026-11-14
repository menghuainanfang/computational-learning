# 资料与实现来源

- [团队课程结构与内容规范](https://github.com/suanhaitech/meshflow_group/blob/main/courses/README.md)：课程编号、目录职责及建设状态。
- [团队FEM-101课程](https://github.com/suanhaitech/meshflow_group/tree/main/courses/100-foundations/FEM-101-reaction-advection-diffusion)：已有模型、有限元推导与程序资料。
- [FEALPy](https://github.com/weihuayi/fealpy)：网格、Lagrange空间、积分器、装配与边界处理。

本章接口已与本机FEALPy源码核对，涉及lagrange_fe_space.py、scalar_convection_integrator.py、bilinear_form.py与dirichlet_bc.py。Céa估计与伴随误差论证在正文展开；凸多边形上的椭圆正则性作为已有理论结论使用，未声明具体教材页码。
