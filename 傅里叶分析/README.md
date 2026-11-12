# 第24周傅里叶分析实验

[任务原文](任务要求.md) · [学习笔记](learning_note.md) · [辅助讲义](傅里叶分析学习讲义.md)

## 运行

本机验证环境：Python 3.13.2，NumPy 2.3.4，Matplotlib 3.10.7。不依赖FEALPy。

```powershell
python -m pip install numpy==2.3.4 matplotlib==3.10.7
python py/傅里叶实验.py
```

在本目录执行。程序从自身路径定位results，运行其他工作目录时也可传入程序的绝对路径。图像使用Agg后端，无需弹窗或图形桌面。

## 输出与验收

- 实验A：固定种子的复随机序列N=8、64、256；手写DFT和基2FFT对照NumPy，预热一次后重复10次平均计时。
- 性质：长度8实序列的共轭对称、Parseval、重构、线性、循环右移；检查基2FFT拒绝空序列、非2幂和非一维输入。
- 实验B：任务指定的5 Hz与30 Hz信号、10 Hz理想低通、长度4补零卷积。
- 误差统一使用绝对容差1e-10；频率和幅值核对显式设置rtol=0。运行通过表示这些指定输入通过检查，不代表完成任意信号的算法验证。
- results包含原始JSON、三份CSV和同一对比图的PNG、SVG版本。计时随硬件和负载变化，重跑后应同步笔记中的计时表。

学习笔记按四组问题和两组实验组织，不重复粘贴完整程序。辅助讲义提供可选推导，不增加必做实验。

## 图像与团队提交

完整材料及图像保存于[个人仓库的傅里叶分析目录](https://github.com/menghuainanfang/computational-learning/tree/main/傅里叶分析)。个人版保留PNG、SVG和数值数据。

培训仓库目录为 `weekly_outputs/official_training_stage_2/week_24/qijing/`。培训版笔记使用个人仓库图片的在线地址；该目录只提交文档、程序和数值JSON/CSV，不提交图像。培训版复现命令：

```powershell
python py/傅里叶实验.py --no-plots
```

默认运行生成图像；`--no-plots` 只生成数据。可用 `--output 路径` 指定输出位置。教材PDF保存在本机，不随作业提交。
