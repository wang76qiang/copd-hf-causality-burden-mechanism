# T1 — GBD–MR 桥接分析（任务 B4）

**日期**: 2026-07-26 · **脚本**: `scripts/t1_bridge.py` · **数据表**: `results/t1_bridge_table.csv`

## 1. 两套口径的量化对比

**前提假设**：2021 年全球 HF 总患病数 ≈ 6,000 万（文献常用量级，GBD 2019 约 5,600–6,400 万；本分析显式标注为假设值，敏感性可线性缩放）。GBD 归因数 3,613,136 已经任务 A5a 从原始 IHME csv 验证。

| 框架 | COPD 患病率假设 | 效应量 | PAF | 2021 归因 HF 病例 |
|---|---|---|---|---|
| **GBD 比较风险评估** | 3% / 5% / 10% | **隐含 RR = 3.14 / 2.28 / 1.64** | **隐含 PAF = 6.02%** | **3,613,136** |
| **MR + Levin PAF** | 3% / 5% / 10% | OR = 1.15 (IVW FE) | 0.45% / 0.75% / 1.49% | 270,208 / 448,998 / 891,326 |

- GBD 框架隐含 PAF（6.02%）是 MR-Levin 中心 PAF（0.75%）的 **8.0 倍**；同一 6,000 万基数下，GBD 归因 361 万 vs MR 外推 44.9 万（中心）。
- 反推：要让 MR 的 OR=1.15 解释 GBD 的 361 万归因病例，COPD 人群患病率需高达 ~35%（远超实际）；要让 GBD 隐含 RR≈2.28 与 MR 一致，则 GBD 归因数需下修约 8 倍。
- 注：V8 中 "≈22.5 万例（中心）" 系将 PAF 0.748% 乘以错误的 3,003 万基数（即 A5a 裁定的聚合错误数），且该基数本应为"HF 总患病数"而非"COPD 归因数"。按 6,000 万 HF 总基数重算，中心值为 **44.9 万例**。

## 2. 方法学差异（为何不必一致）

1. **估计对象不同**：GBD 比较风险评估（CRA）估计的是**当期人群水平**的归因负担——基于 COPD 现患的 RR（来自观察性队列的 meta 分析）与暴露分布，包含反向因果（HF 加重 COPD 症状/诊断）、共同病因（吸烟、衰老、炎症）与检出偏倚的全部表型关联；MR 估计的是**终身遗传工具化的 COPD 易感性**对 HF 的因果效应，不受反向因果与后天混杂影响，但只反映"遗传预测的 COPD"这一子集。
2. **RR 来源不同**：GBD 隐含 RR≈2.28（p=5% 时）明显高于 MR 的 OR=1.15，差异方向与"观察性关联被混杂/反向因果放大"一致。
3. **暴露定义不同**：GBD 的暴露是临床/肺功能确诊的 COPD 现患状态；MR 的暴露是遗传倾向（lifetime liability），效应按"每单位对数优势"标度，通常小于临床发病效应。
4. **PAF 的种群语义不同**：GBD PAF 回答"若消除 COPD，当期 HF 负担减少多少"（含不可消除的共享病因路径）；Levin-MR PAF 回答"若 COPD 的遗传因果路径被消除，HF 减少多少"。两者不应直接互换。

## 3. 稿件讨论段落建议（英文，148 词）

> The GBD comparative risk framework attributed 3.61 million HF cases to COPD in 2021, implying a population-attributable fraction of ~6.0% against an assumed 60 million total HF cases, whereas our Mendelian randomization estimate (OR 1.15) yields a Levin PAF of only 0.45–1.49% (~0.27–0.89 million cases). This eight-fold gap is expected rather than contradictory: GBD attribution captures the full phenotypic association, including shared causes (smoking, ageing, systemic inflammation), reverse causation, and detection bias, whereas MR isolates the lifelong causal effect of genetically proxied COPD liability, free from such amplification. The two frameworks thus answer complementary questions—current attributable burden versus causal preventability—and we report both explicitly, cautioning against interpreting the GBD figure as the caseload removable by COPD prevention alone.

## 4. 局限

- 6,000 万 HF 总基数为假设值；若采用 5,600 万–6,400 万区间，GBD 隐含 PAF 在 5.6%–6.5% 之间，结论不变。
- GBD 隐含 RR 的反推依赖 Levin 公式（假设 PAF 与 RR 的单一路径关系），仅用于量级比较。
