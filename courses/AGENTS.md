# 课程编辑约定

本文件适用于 `courses/` 下的教材、程序说明和数值验证材料。

- 修改数学内容前读取并应用 [Meshflow数学符号规范](https://github.com/suanhaitech/meshflow_group/blob/main/standards/meshflow_mathematical_notation_standard.md)，按对象类型区分场、系数向量、矩阵和分量。
- 项目当前采用 `\boldsymbol{}`；与SuanhaiOS的 `\bm{}` 约定差异尚待协调，不宣称已修订上位标准。
- 保持控制方程、假设、边界条件及证明的含义；代码变量与数学记号的关系应明确。
- 正文使用Markdown。生成HTML时读取[HTML版式规范](https://github.com/suanhaitech/meshflow_group/blob/main/standards/html/meshflow_html_layout_standard.md)，按连续阅读或逐页放映分别选择模板。
- 图片保存在各课程 `notes/figures/`，正文通过本仓库的图片URL引用。保留复现程序与必要的小型结果。
- 区分理论证明、数值观察与计划；原始引用材料保留来源与其自身记号。
