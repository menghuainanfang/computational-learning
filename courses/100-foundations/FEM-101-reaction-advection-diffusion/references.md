# 资料与实现来源

## 研究项目资料

- [CDR数学原理讲义](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/02_%E5%AD%A6%E4%B9%A0%E8%AE%B2%E4%B9%89/01_%E6%95%B0%E5%AD%A6%E5%8E%9F%E7%90%86%E8%AE%B2%E4%B9%89.md)
- [CDR程序入门讲义](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/02_%E5%AD%A6%E4%B9%A0%E8%AE%B2%E4%B9%89/02_%E7%A8%8B%E5%BA%8F%E5%85%A5%E9%97%A8%E8%AE%B2%E4%B9%89.md)
- [FEALPy库与项目的边界](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/03_FEALPy%E5%AD%A6%E4%B9%A0/00_%E5%BA%93%E4%B8%8E%E9%A1%B9%E7%9B%AE%E7%9A%84%E8%BE%B9%E7%95%8C.md)
- [网格与函数空间](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/03_FEALPy%E5%AD%A6%E4%B9%A0/01_%E7%BD%91%E6%A0%BC%E4%B8%8E%E5%87%BD%E6%95%B0%E7%A9%BA%E9%97%B4.md)
- [积分器与装配](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/03_FEALPy%E5%AD%A6%E4%B9%A0/02_%E7%A7%AF%E5%88%86%E5%99%A8%E4%B8%8E%E8%A3%85%E9%85%8D.md)
- [边界、求解与源码阅读](https://github.com/menghuainanfang/research-workspace/blob/main/%E5%AD%A6%E4%B9%A0%E7%BB%83%E6%89%8B%E9%A1%B9%E7%9B%AE/cdr-solver/03_FEALPy%E5%AD%A6%E4%B9%A0/03_%E8%BE%B9%E7%95%8C%E6%B1%82%E8%A7%A3%E4%B8%8E%E6%BA%90%E7%A0%81%E9%98%85%E8%AF%BB.md)

数学原理讲义用于核对变系数扩散、空间与代数误差、Péclet数和一致SUPG残差；程序讲义及四篇FEALPy阅读材料用于核对维数与次数、回调数组约定、局部到全局装配和库／项目分工。零基础编程、含时推进与多物质接口按课程分工保留在原材料中。

## 团队课程资料

- [模型与弱形式](https://github.com/suanhaitech/meshflow_group/blob/main/courses/100-foundations/FEM-101-reaction-advection-diffusion/notes/01-%E7%A8%B3%E6%80%81%E6%A8%A1%E5%9E%8B%E4%B8%8E%E5%BC%B1%E5%BD%A2%E5%BC%8F.md)
- [局部几何、自由度与装配](https://github.com/suanhaitech/meshflow_group/blob/main/courses/100-foundations/FEM-101-reaction-advection-diffusion/notes/02-%E4%BD%8E%E7%BB%B4%E6%89%8B%E7%AE%97%E4%B8%8E%E4%BA%8C%E7%BB%B4%E8%A3%85%E9%85%8D.md)
- [实现、代数求解与强对流](https://github.com/suanhaitech/meshflow_group/blob/main/courses/100-foundations/FEM-101-reaction-advection-diffusion/notes/03-FEALPy%E6%B1%82%E8%A7%A3%E4%B8%8E%E6%95%B0%E5%80%BC%E9%AA%8C%E8%AF%81.md)

三篇专题沿用团队的稳态标量定位；连续证明与教材组织在本课正文展开。新增的界面通量解释由分域分部积分得到，高阶积分次数由被积多项式次数推导，不作为新增数值实验。

## FEALPy实现

- [FEALPy项目](https://github.com/weihuayi/fealpy)

接口按本机实际源码核对：`lagrange_fe_space.py` 的自由度、基函数与物理梯度；`scalar_diffusion_integrator.py` 默认装配路径中的系数求值、求积与`bilinear_integral`。脚本仍采用本课已验证的二维常系数P1模型，未声称本轮验证变系数、材料界面、高阶或SUPG。凸多边形椭圆正则性作为已有理论结论使用，未声明具体教材页码。
