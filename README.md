# 课程代码仓库

本科课程学习中凡是需要写代码的内容——作业、实验、课程项目、课堂练习——统一存放在本仓库,按课程分类,方便回溯与迁移。

## 目录结构

一门课一个顶层文件夹,文件夹名用英文小写 + 连字符,各课程互不干扰:

```
.
├── README.md                    # 本文件:仓库总览与使用规范
├── .gitignore
├── numerical-analysis/          # 课程 1:数值分析
│   ├── README.md                # 该课程的说明(用下方模板)
│   ├── lab1-interpolation/      # 实验 / 作业:一个任务一个文件夹
│   └── lab2-numerical-integration/
└── python-programming/          # 课程 2:Python 程序设计
    ├── README.md
    └── ...
```

## 加一门新课

1. 新建顶层文件夹,命名为课程英文缩写,如 `numerical-analysis`;
2. 课程内按「实验 / 作业 / 专题」再分子文件夹,一个任务一个文件夹;
3. 把下面的课程 README 模板填好后,放到该课程根目录。

> 文件夹名建议用英文:命令行和跨平台环境对中文路径支持不一致,容易踩编码乱码的坑,排查成本很高。

### 每门课程的 README 模板

```markdown
# <课程名>

- 学期:2026-2027 第一学期
- 语言 / 工具:Python 3.13
- 依赖:见 requirements.txt(如有)

## 实验列表

| 编号 | 内容 | 状态 |
| ---- | ---- | ---- |
| lab1 | 三次样条插值 | 已完成 |
```

## 规范约定

### 环境与依赖

- Python 课程:在课程根目录放 `requirements.txt`,并注明所用的 Python 版本;
- C/C++ 课程:在课程 README 里注明编译器与构建方式(如 CMake 版本、编译命令);
- 不提交虚拟环境(`venv/`、`.venv/`)和运行产物(图片、日志、数据)——它们都能由代码重新生成。

### Git 提交

- 粒度:完成一个实验 / 一次作业,提交一次;
- message 格式:「课程 —— 做了什么」,例如:
  - `数值分析:完成三次样条插值实验`
  - `Python程序设计:实现命令行通讯录`
- push 之前先 `git status`、`git diff` 检查一遍,确认没有混入无关文件。

### 代码来源

- 尽量独立完成;
- 参考了教材、他人代码,或直接引用教材例题的,在文件开头注明出处。

## 课程清单(随学习进度更新)

| 文件夹 | 课程 | 语言 | 状态 |
| ------ | ---- | ---- | ---- |
| `numerical-analysis/` | 数值分析 | Python | 学习中 |
| `python-programming/` | Python 程序设计 | Python | 待开始 |

## .gitignore(仓库根目录,可直接复制使用)

```gitignore
# Python
__pycache__/
*.py[cod]
.venv/
venv/

# IDE
.idea/
.vscode/

# Jupyter
.ipynb_checkpoints/

# 运行输出与数据(按需取舍)
output/
data/

# 系统文件
.DS_Store
Thumbs.db
```
