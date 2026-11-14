# FEM-101｜稳态反应–对流–扩散方程

本课从控制体守恒建立模型，连接弱形式、适定性、有限元装配、FEALPy实现和收敛验证。先修为泛函分析、数学物理方程及团队有限元基础培训。

## 学习目标与阅读入口

读者应能解释守恒与非守恒形式的关系，推导带Dirichlet边界的弱形式，说明适定性与误差估计所需条件，并将单元积分、自由度映射和边界消元对应到FEALPy实现。

- [教材正文：稳态输运问题的有限元方法](notes/01-稳态输运问题的有限元方法.md)
- [实验与自主问题](labs/README.md)
- [验证内容与结果](tests/README.md)
- [资料与实现来源](references.md)

正文作为一个完整主题单元连续阅读，代码与图像紧邻相应推导。

## 运行

从本课程目录运行：

```text
python labs/稳态输运实验.py
```

依赖FEALPy、NumPy、SciPy和Matplotlib，使用NumPy后端。本机已验证解释器为 `D:/桌面/movingmesh/移动网格/.venv/Scripts/python.exe`，FEALPy源码为 `C:/Users/10251/fealpy/fealpy`。其他安装版本需核对接口。

## 范围与状态

状态：**In Progress**。教材正文和配套实验已具备可复现材料，语言修订与真实读者试学继续进行。

实现范围为二维、常系数、连续P1元及全Dirichlet边界。连续理论介绍一般维数，数值验证不代表一般维数与高阶实现。多分量和时间推进分别衔接FEM-102与TIM-101/FEM-103。

课程目录采用 README、notes、labs、tests、references.md 结构。对应团队目录为 `courses/100-foundations/FEM-101-reaction-advection-diffusion/`。

## 图片与资源

正文图片通过本个人仓库的原始文件链接引用，源图保存在 `notes/figures/`。复制正文到团队仓库时可保留这些图片链接。课程程序可在本地再生成图片与精简数据。
