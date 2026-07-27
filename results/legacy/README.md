# legacy_archive — V8 留档证据归档说明

归档日期：2026-07-26（T1 任务 X）
来源目录：`E:/COPD+HF/V8_TopJournal_RealData/R_scripts/results/`
原始文件生成日期：2026-07-05（文件系统 mtime）
生成脚本：`E:/COPD+HF/V8_TopJournal_RealData/R_scripts/01–08_*.R`
背景：暴露 COPD ebi-a-GCST90018807（10 工具 SNP），结局 HF ebi-a-GCST009541；
亚洲复制 BBJ bbj-a-103 → bbj-a-109。归档目的：为 T1 顶刊版提供可追溯的既有证据，
避免重复计算；凡 T1 阶段因 OpenGWAS 封锁无法重跑的分析均引用本目录结果。

## 文件清单与含义

### CAUSE 分析（01_cause_analysis.R）
- `cause_results_cause_summary.json`：CAUSE 拟合摘要。注意 status=fallback_causal_model：
  完整 CAUSE 模型比较需要全基因组汇总统计（>100k 变异），本地仅 10 个工具，
  完整 LOO-CV 失败，退回单独拟合 causal model。因果效应 gamma(log OR)=0.138，
  95%CI [0.031, 0.235]，OR=1.15 [1.03, 1.27]（带明确 caveat）。
- `cause_results_cause_fit_fallback.rds`：fallback 拟合对象（R 序列化）。
- `01_cause_analysis.log`：运行日志。

### MR-Clust 聚类（02_mrclust_analysis.R）
- `mrclust_results_clusters.csv`：10 个 SNP 的聚类归属与后验概率
  （cluster 1 均值 0.215 含 6 SNP；rs11525583 等归入 Null/其他簇 → 提示多效性分层）。
- `mrclust_results_cluster_table.csv`：各簇均值/大小汇总表。
- `mrclust_results_mrclust_summary.json`：机器可读摘要（n_snp=10, n_clusters=2）。
- `mrclust_results_mrclust_fit.rds`：mrclust 拟合对象。
- `02_mrclust_analysis.log`：运行日志。

### BAPC 疾病负担预测（03/08_bapc_projection.R）
- `bapc_results_aggregated_cases.csv`：按年龄段聚合的 GBD 病例数。
- `bapc_results_apc_predictions.csv`：APC 模型逐年预测。
- `bapc_results_yearly_projections.csv`：1990–2050 逐年投影（含情景）。
- `bapc_results_bapc_summary.json`：BAPC 拟合摘要。
- `03_bapc_projection.log`：运行日志。

### SERPINE1(PAI-1) pQTL 药物靶点 MR（04_serpine1_pqtl_mr.R）——阴性
- `serpine1_pqtl_mr_summary.json`：循环 PAI-1（prot-c-2925_9_1, SomaLogic, 5 工具）
  → HF：IVW beta=0.0122 (se 0.0302) p=0.686；Weighted median p=0.648；
  MR-Egger p=0.711；Egger 截距 p=0.759；IVW Q_p=0.031（轻度异质性）。
- `serpine1_pqtl_mr_mr_results.csv` / `_pleiotropy.csv` / `_harmonised.csv`：
  MR 结果、多效性检验、协调后数据。
- `04_serpine1_pqtl_mr.log`：运行日志。

### SERPINE1 lead SNP IEU PheWAS（05_phewas_serpine1.R）
- `phewas_serpine1_ieu_all.csv`：rs7860931 在 IEU OpenGWAS 的全表型关联（全部记录）。
- `phewas_serpine1_ieu_p1e-3.csv`：p<1e-3 子集。
- `phewas_serpine1_ieu_significant.csv`：显著子集。
- `phewas_serpine1_summary.json`：摘要。
- `05_phewas_serpine1.log`：运行日志。
- （T1 阶段已补充 FinnGen R12 PheWAS，见 `../t1_c1_finngen_r12_phewas_rs7860931.csv`。）

### 亚洲复制 MR（06_asian_mr_update.R）
- `06_asian_harmonised_data.csv`：BBJ bbj-a-103→bbj-a-109 协调数据（6 SNP）。
- `06_asian_mr_bbj_congestive_hf_result.csv`：IVW OR=0.933 (0.848–1.027) p=0.155；WM p=0.396。
- `06_asian_mr_heterogeneity.csv`：IVW Q=8.78 (df=5) p=0.118。
- `06_asian_mr_pleiotropy.csv`：Egger 截距检验。
- `06_asian_mr_wald_ratios.csv`：逐 SNP Wald ratio。
- `06_asian_mr_update.log`：运行日志。
- `_search_asian_hf_gwas.log` / `_search_asian_hf_gwas2.log`：当时 OpenGWAS 检索记录，
  结论：bbj-a-109 为平台上最大东亚 HF GWAS。

### 其他
- `hf_gwas_candidate_info.csv`：HF 结局 GWAS 候选数据集信息（选择 ebi-a-GCST009541 的依据）。

## 引用注意
CAUSE 结果为 fallback 模式（非完整模型比较），论文中引用时必须保留其 caveat；
BAPC 目录中另有 Python 代理投影（real_data_results.json 内 bapc_proxy），与正式
BAPC 结果区分使用。
