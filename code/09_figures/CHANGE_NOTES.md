# MedComm 修图变更说明（Figure 5 / 6 / 7）

日期：2026-08-15
工作副本目录：`C:\Users\fhj\.scratch\ECM1_Medcomm_preaudit\figure_fix_scripts\`
源数据/原始脚本（未改动）：`C:\Users\fhj\.scratch\ECM1_Medcomm_preaudit\zenodo\extracted\copd-hf-causality-burden-mechanism-1.0.2\`
输出目录：`E:\投稿版 ECM1_submission\Medcomm_revised\MedComm_figures\`

## 生成方法

```bash
cd C:\Users\fhj\.scratch\ECM1_Medcomm_preaudit\figure_fix_scripts
set REPRO_ROOT=C:\Users\fhj\.scratch\ECM1_Medcomm_preaudit\zenodo\extracted\copd-hf-causality-burden-mechanism-1.0.2
py make_fig5.py   # -> output\Fig5_horizon_2050.{png,pdf}
py make_fig6.py   # -> output\Fig6_molecular_convergence.{png,pdf}
py make_fig7.py   # -> output\Fig7_cellular_executioners.{png,pdf}
```

- 修改后的脚本：`make_fig5.py`、`make_fig6.py`、`make_fig7.py`（本目录；`figstyle.py` 为原样拷贝）。
- 三个脚本仅改动输出路径（写到本目录 `output\`，不污染源仓库）及下述面板；数据、统计量、其余面板完全未动。
- PNG 复制为 `Figure_5/6/7.png`，PDF 复制为 `Figure_5/6/7.pdf`（矢量版，投稿备用）。

## 各图改动

### Figure 5（panel C）
- `10%: PAF 1.48% ≈ 0.89 M` 注释原位于 (10.9, 1010)，压在上升曲线上；现移到 (6.0, 2950) 空白区并加引线。
- 同步把 `3%`、`5% central` 两条注释也移到曲线与 GBD 虚线之间的空白带（0.9,1250 / 2.6,2100），统一加引线、加粗、字号 8→8.5。
- 更新陈旧代码注释：B 面板瀑布图旧注释 "10.93-2.5-1.8=6.63" 数值错误，改为已验证的 "10.93 − 2.93 − 1.83 = 6.17 M"（图中 6.2 M 不变）。
- 图注无需改动（数值、面板内容均未变）。

### Figure 6（panel F）
- 热图由全部 41 基因（行高约 7pt 字号，180 mm 下不可读）改为 15 个正文点名/关键基因：
  F13A1, S100A8, S100A12, CD163, SMOC2, SFRP4, SERPINE1, THBS1, SPP1, IL1R2, EGR1, CDKN1A, AREG, SOCS3, IL1RL1（均存在于 t1_crosscohort_validation.csv）。
  LUM/ASPN/ACE 不在 41 基因交集表内，仍由 panel E 箱线图展示，未纳入热图。
- 行基因名字号 7→11（hub 基因加粗，与 panel C 的 amber hub 集合一致），队列名字号 8.5→10.5，FDR 显著性标记 "•"→"●"（8→12pt），新增白色单元格分隔线，色阶条标签 8→10。
- 完整 41 基因矩阵数据仍在 `results/t1_crosscohort_validation.csv`，可整体移至补充材料。
- **图注需同步修改**（见下）。

### Figure 7（panel E + 标签字号）
- panel E 原 `t[:20]+"…"` 截断 endpoint 名称（如 "Asthma (only as main…"）。现改为显示 8 个最小 P 值 endpoint 的**完整名称**（最长两条手动折成两行，无省略号）：
  Asthma (only as main-diagnosis) (more control exclusions) / Asthma (more control exclusions) / Anorexia / Internal derangement of knee / Chronic rhinitis, nasopharyngitis and pharyngitis / Other primary thrombocytopenia / Cardiomyopathy in other diseases / Behçet disease。
- panel E 顶注补充 "(8 smallest-P endpoints shown)"，保证"2,470 endpoints · none significant"表述仍准确。
- 字号增大：A 刻度 8→9 / 8.5→9.5、图例注 7.2→8；B 刻度 8.5→9 / 9.5；C y 标签 7.2→8.2、图例 7.5→8.2；D y 标签 8.5→9.2；E y 标签 6.8→8、Bonferroni 标注 7→8、散点 55→65。
- D 面板 x 轴上限 0.335→0.42，避免 E 面板加宽的长标签与 D 的 "+0.23" 等数值标签相碰（已目视确认无重叠）。
- F 面板为用户提供的机制示意图（Fig7F.jpg），未改动。

## 输出文件规格（已核验）

| 文件 | 像素 | DPI | 大小 |
|---|---|---|---|
| Figure_5.png | 3866×2366 | ≈300 | 371 KB |
| Figure_6.png | 4439×3108 | ≈300 | 737 KB |
| Figure_7.png | 4436×3759 | ≈300 | 3.8 MB |

白底、尺寸与原图基本一致；另附同名矢量 PDF。

## 视觉检查结果（逐图读图确认）

- Figure 5C：三条 PAF 注释均在空白区、带引线，不与曲线/散点/GBD 虚线重叠；无截断。✔
- Figure 6F：15 行基因名、4 列队列名、色阶、显著性圆点全部清晰可读；边缘行（SMOC2）圆点经局部放大确认未被裁剪。✔
- Figure 7E：8 条 endpoint 全称完整无省略号、无截断；与 D 面板数值标签无重叠；A/B/C/D 标签增大后无碰撞；F 插图完整。✔

## 需要主文稿图注同步修改的建议

1. **Figure 6 图注（panel F）**：原文若写 "heat map of the 41-gene intersection / 41-gene program across four cohorts"，建议改为：
   > "(F) Cross-cohort direction consistency (log₂FC) of 15 representative genes from the 41-gene COPD–HF intersection across bulk COPD lung (GSE57148), single-cell COPD lung (CELLxGENE), bulk failing heart (GSE57338) and single-nucleus failing heart (Reichart et al.); dots mark FDR < 0.05 within cohort, white cells = not measured. Hub genes (STRING/cytoHubba) in bold. The full 41-gene matrix is provided in the Supplementary Material (Table Sx / Fig. Sx)."
2. **Figure 7 图注（panel E）**：若原文写 "top 10 endpoints"，建议改为：
   > "(E) PheWAS of the SERPINE1 cis variant rs7860931 across 2,470 FinnGen R12 endpoints; the eight smallest-P endpoints are shown with full names; none passes the Bonferroni threshold (2.0×10⁻⁵), indicating a clean safety window."
3. Figure 5 图注无需改动。
4. 补充材料建议新增：完整 41 基因 × 4 队列热图（可由原脚本取消 KEY_GENES 过滤生成）或直接引用 `t1_crosscohort_validation.csv`。
