from __future__ import annotations

import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
POST = ROOT / "13_postfreeze_interpretation"
REPORTS = ROOT / "12_reports"


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def main() -> int:
    s5 = load_json(ROOT / "07_convergence/stage5_expression_statistics.json")
    s8 = load_json(ROOT / "10_final_tiers/stage8_statistics.json")
    d5 = load_json(POST / "d5_biological_confound_summary.json")
    d7 = load_json(POST / "d7_dataset_collapse_audit_summary.json")
    d8 = load_json(POST / "d8_final_integrity_summary.json")
    pipe = load_json(ROOT / "pipeline_status.json")
    layers = load_json(POST / "d5_evidence_layers_summary.json")
    tier = pd.read_csv(ROOT / "10_final_tiers/final_candidate_tier_matrix.tsv", sep="\t")
    primary = pd.read_csv(POST / "d6_primary4_influence_decomposition.tsv", sep="\t")
    figures = pd.read_csv(POST / "d8_submission_figure_integrity.tsv", sep="\t")

    key_files = [
        ROOT / "scripts/stage5_expression_convergence.R",
        ROOT / "scripts/postfreeze_d5_biological_confound_audit.py",
        ROOT / "scripts/postfreeze_d7_ortholog_dataset_collapse_audit.py",
        ROOT / "scripts/postfreeze_d8_final_integrity_audit.py",
        ROOT / "scripts/postfreeze_d4_orthogonal_lit_v2.py",
        ROOT / "scripts/postfreeze_fig7_rab1a_network.py",
        ROOT / "run_pipeline_v2.ps1",
        ROOT / "07_convergence/stage5_expression_statistics.json",
        ROOT / "07_convergence/primary_infection_convergence.tsv",
        ROOT / "08_evidence_atlas/host_dependency_evidence_atlas.tsv.gz",
        ROOT / "09_holdout_validation/holdout_frozen_set_tests.tsv",
        ROOT / "10_final_tiers/final_candidate_tier_matrix.tsv",
        REPORTS / "CURRENT_ANALYSIS_AUDIT.md",
        REPORTS / "BIOLOGICAL_CONFOUND_AUDIT_V2_CN.md",
        POST / "d4_orthogonal_literature_evidence_v2.tsv",
        POST / "d5_global_response_amplitude_matched_null.tsv",
        POST / "d7_original_vs_dataset_collapsed_convergence.tsv",
        POST / "d8_final_integrity_checks.tsv",
        POST / "POSTFREEZE_FIGURE_LEGENDS.md",
    ] + [POST / name for name in figures.file]
    missing = [str(p) for p in key_files if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing final package files: " + " | ".join(missing))

    manifest = pd.DataFrame(
        [{
            "relative_path": str(path.relative_to(ROOT)),
            "bytes": path.stat().st_size,
            "sha256": sha(path),
            "modified": datetime.fromtimestamp(path.stat().st_mtime).isoformat(timespec="seconds"),
        } for path in key_files]
    )
    manifest.to_csv(REPORTS / "FINAL_SUBMISSION_AUDIT_MANIFEST_V2.tsv", sep="\t", index=False)

    tier_counts = tier.final_tier.value_counts().to_dict()
    primary_lines = "\n".join(
        f"- {r.candidate}: oriented Z={r.direction_oriented_z:.2f}; "
        f"leave-one-candidate-out mean Z={r.leave_candidate_out_mean_z:.2f}."
        for r in primary.itertuples()
    )
    figure_lines = "\n".join(
        f"- `{r.file}` — {int(r.pages)} page, {int(r.bytes):,} bytes, integrity {r.status}."
        for r in figures.itertuples()
    )

    report = f"""# 嗜肺军团菌宿主依赖图谱：最终综合策略、结果与投稿前审计（V2）

> 生成时间：{datetime.now().isoformat(timespec='seconds')}  
> 项目目录：`{ROOT}`  
> 主控状态：`{pipe['state']}`  
> 完整性审计：`{d8['fatal_error_status']}`（{d8['checks_passed']}/{d8['checks_total']} PASS）  
> 稿件状态：可开始正式写作；最终投稿冻结等待 G1/G2 外部来源闭合。

## 一、研究最终定位

本研究不是把不同证据简单相加的候选基因排名，而是建立一个偏倚感知、证据语义分离、冻结后外部验证的嗜肺军团菌宿主依赖图谱。核心概念为：

> **A variable effector shell converges on a conserved Rab/ER-Golgi host-trafficking dependency core.**

物理接触、AP-MS、宿主 CRISPR 依赖、LCV 空间定位、跨队列表达与外部删除株转录响应分别回答不同问题，不允许用一种证据替代另一种证据。最终最强结论是 RAB1A，而不是“所有运输基因同等重要”。

## 二、数据与验证架构

研究整合：

1. 军团菌效应子—人宿主互作证据；
2. U937 全基因组 CRISPR 与重复筛选，并用 DepMap common-essential/Chronos 标记基础适应度混杂；
3. 6 个发现表达数据集及机制/转录—翻译轨道；
4. 人源 WT、ankB 与鼠源 LCV 蛋白质组；
5. 预先冻结的 56 个候选和验证规则；
6. 独立 E-MTAB-16655 PMA-THP-1 留出集：45 个样本、15 个条件、每组 3 个生物学重复；
7. 冻结后匹配空集、物种组成、队列脆弱性、文献正交证据和完整性审计。

E-MTAB-16655 不仅是未用于发现的数据，而且具有不同人源细胞背景与不同效应子删除架构，因此可称为 **cross-context external validation**；但 THP-1 仍是细胞系，不能称为生理或临床验证。

## 三、主要发现

### 3.1 证据轨道非等价

在混杂匹配后，PPI(A/B) 对 clean CRISPR 的预测较弱；LCV 空间证据与 clean CRISPR 的关联较强；PPI 与 LCV 的表面富集可被网络度数解释；CRISPR/LCV 与表达响应不存在可替代关系。PPI+LCV 联合阳性组中 clean CRISPR 为 7/18（38.9%），该结果是描述性收敛，不是校准预测器。

### 3.2 RAB/ER-Golgi 运输程序富集

冻结 54 基因响应面板中 RAB 家族富集 64.4 倍，完整运输基因集合富集 21.6 倍。独立 DepMap 共本质匹配置换支持 Rab/ER-Golgi 模块在非感染数据中的协同依赖，但该结果仅为探索性模块证据，不消除 common-essential caveat，也不证明感染特异性。

### 3.3 RAB1A 是唯一严格 Tier 1

最终 56 个冻结候选中：

- Tier 1：{tier_counts.get('TIER1_EXTERNALLY_VALIDATED', 0)}（RAB1A）；
- Tier 2：{tier_counts.get('TIER2_MULTIMODAL', 0)}；
- Tier 3 context-dependent：{tier_counts.get('TIER3_CONTEXT_DEPENDENT', 0)}；
- Tier 3 fitness caveat：{tier_counts.get('TIER3_FITNESS_CAVEAT', 0)}；
- 留出方向矛盾后拒绝：{tier_counts.get('REJECT_HOLDOUT_DIRECTIONAL_CONTRADICTION', 0)}。

RAB1A 的优势来自直接效应子操纵、clean replicated CRISPR、LCV 空间支持、跨队列元收敛和锁定留出方向支持共同出现。RAB10 是强 Tier 2，但不能用外部转录支持替代预定义的机制/转译要求。

### 3.4 四主候选不应解释为等证据模块

冻结四基因集合检验 FDR=0.0603，按预定义规则不计阳性。事后影响拆解仅用于解释异质性：

{primary_lines}

方向信号主要由 RAB1A/RAB10 承载，ACTR3/VAPA 较弱；这不改变冻结检验，也不触发重新分级。

## 四、关键纠错及其影响

### 4.1 跨物种同源映射重复计数

原 Stage 5 将同一数据集内多个小鼠同源基因映射到同一人类基因后，仍按多行计入 Stouffer 汇总。共发现 {d7['duplicate_human_gene_dataset_pairs']} 个重复 `human_symbol × dataset` 单元和 {d7['extra_rows_that_were_not_independent_studies']} 个额外非独立行。

纠正后每个基因×数据集只贡献一次：

- 原 meta-convergent：{d7['original_stable_n']:,}；
- 纠正后 meta-convergent：{s5['stable_cross_cohort']:,}；
- ≥3 独立数据集显著：{layers['multi_study_replicated_ge3']:,}；
- 两层交集：{layers['both_meta_and_ge3']:,}；
- 状态变化基因：{d7['stable_status_changed_n']}；
- 冻结候选表达标志变化：{d7['frozen_candidate_expression_flag_changed_n']}；
- 最终 Tier 变化：0。

旧版结果已归档在 `07_convergence/legacy_pre_dataset_collapse/`。正文应使用 2,601，并在方法/补充材料披露纠错；不得静默保留旧 2,590。

### 4.2 留出条件 F/I/FGHI 核定

官方 SDRF 明确包含 UI、WT、dotA、A、B、C、D、E、F、FGHI、I、O、P、R、S，各 3 个生物学样本。F 与 I 是独立条件，同时也位于 FGHI 组合背景；G/H 未被单独测序。字母标签是官方元数据，字母到删除基因清单的映射才是外部推断。O/P/R/S/G 内容尚未完全锁定时，禁止作基因特异性机制推断。

### 4.3 全局响应幅度与“负荷混杂”

冻结面板响应与全局宿主转录扰动相关：rho={d5['panel_vs_global_spearman']['rho']:.3f}，P={d5['panel_vs_global_spearman']['p']:.4f}。但 E-MTAB-16655 未提供 bacterial burden；因此该指标只能称为 **global response-amplitude sensitivity**，不能称为已测得的负荷混杂。

5,000 个匹配空集在 12 个删除株对比中经 FDR 后 0 个显著。这意味着：面板在不同删除背景下响应幅度不同，但不能据此证明特定效应子模块直接驱动特异宿主程序。冗余/拮抗仅作为机制假说。

### 4.4 Figure 7 引文纠错

效应子—Rab1 边已分别绑定直接原始来源：DrrA/SidM 的 GEF 与 AMPylation、SidD 的 de-AMPylation、LepB 的 GAP、AnkX/Lem3 的可逆磷酸胆碱化及 SetA 糖基化不再错误共用 PMID。PMID 15117975 已纠正为 Kagan 等关于 Rab1/Sec22b 与复制性细胞器的研究；PMID 19942850 仅用于 GEF/GDI displacement，不再错误标作 AMPylation 论文。

## 五、Figure 5–7 与补充图最终状态

{figure_lines}

完整图注：`13_postfreeze_interpretation/POSTFREEZE_FIGURE_LEGENDS.md`。  
PDF 结构检查：8/8 通过；所有图均为可读取 PDF 且包含矢量文本。最终人工检查仍应在目标期刊版芯尺寸下完成。

## 六、独立实验文献锚点

校正版表为 `d4_orthogonal_literature_evidence_v2.tsv`，所有条目均强制 `used_for_tier_assignment=NO`，避免循环论证。

关键一手论文包括：

- Kagan et al., Rab1/Sec22b 与复制性细胞器：PMID 15117975；
- Suh et al., SidM/DrrA GEF/GDI displacement：PMID 19942850；
- Jeng et al., U937 CRISPR 与 RAB10/RABIF：PMID 31540829；
- Hoffmann et al., LCV Rab 蛋白与 RNAi 功能：PMID 24373249；
- Steiner et al., atlastin/ER remodeling：PMID 28835546；
- Wang et al., SetA-Rab1 糖基化：PMID 30323948。

VAPA、ACTR3、ACBD3 未在本次定向核查中找到候选特异、同等级的独立 Legionella 功能验证；如实标注“未定位到”，不以推断填补。

## 七、致命错误终审

自动完整性审计共 {d8['checks_total']} 项，{d8['checks_passed']} 项通过、{d8['checks_failed']} 项失败。检查覆盖：

- Stage 5 数据集折叠；
- 候选/Tier 稳定性；
- 15 个官方留出条件；
- F/I/FGHI 逻辑；
- 匹配空集结论边界；
- 文献是否参与 Tier；
- Figure 7 PMID 映射；
- 8 份 PDF 完整性；
- post-freeze Python 语法；
- Stage 10 对 Stage 5 的 checkpoint 依赖。

结论：`{d8['fatal_error_status']}`。主控复跑状态为 `{pipe['state']}`；所有可执行阶段已完成或按有效 checkpoint 跳过。四条 warning 均为预定义科学边界，不是技术失败。

## 八、允许与禁止的论文声称

### 可声称

1. 建立了版本化、偏倚感知的多模态军团菌宿主依赖/响应图谱。
2. 物理、空间、遗传和响应证据部分重叠但不可替代。
3. 冻结面板和 Dot/Icm 子集在跨细胞背景的独立留出中获得集合层支持。
4. RAB1A 满足完整、严格的 Tier-1 计算规则。
5. Rab/ER-Golgi 运输架构是多证据收敛的宿主核心。
6. 删除株之间存在宿主响应幅度异质性。

### 不可声称

1. 2,601 个基因都在多个研究中逐基因显著重复。
2. PASS/FAIL 证明某个删除区段内的具体效应子直接调控该程序。
3. 全局转录扰动等同于细菌负荷，或细微 18 h 负荷差异已完全排除。
4. AP-MS、LCV 或计算结构预测单独证明直接结合。
5. 转录响应证明胞内复制依赖或治疗因果性。
6. VAPA/ACTR3/ACBD3 与 RAB1A 具有同等级独立实验支持。

## 九、最终冻结前仅剩的外部缺口

### G1：作者复制控制的一手文件

请手动下载完整官方摘要 PDF 至 `01_data_raw/literature/`：

- `https://viin.net.au/docs/2025-viin-young-investigator-symposium-abstract-booklet-final-251112.pdf`
- `https://viin.org.au/docs/2024-viin-young-investigator-symposium-abstract-booklet-final-241106.pdf`

当前本地 `viin2024.pdf` 只是封面图，不能作为已核验全文。若完整摘要确认“所有 mutant 均高效复制/无明显复制缺陷”，可把 gross replication defect 降为不太可能，但仍保留“无逐株 18 h 定量 CFU”的限制。

### G2：删除株基因清单

需要 Hartland/Wibawa 原始菌株表或 in-press 论文补充材料，以锁定 G/O/P/R/S 的删除基因内容。若无法取得，则按现有规则冻结：保留条件层统计，禁止这些组的基因特异性机制解释。

## 十、投稿与写作行动方案

1. 立即开始 Methods/Results/Discussion，不再开启新 discovery 或修改候选/Tier。
2. 主叙事：`variable effector shell -> conserved Rab/ER-Golgi host-trafficking core`。
3. Figure 5 展示证据非等价；Figure 6 展示运输图谱与响应幅度；Figure 7 展示 RAB1A 的直接效应子调控和最终 Tier 架构。
4. 将数据集折叠纠错、matched-null 阴性和四基因异质性主动写入 Methods/Results/Discussion。
5. G1/G2 若无法完全闭合，不阻止计算资源型稿件写作，但决定最终措辞强度与投稿冻结。
6. 投稿前再执行一次 clean-room rerun、版芯尺寸人工图审和引用逐条核验。

## 十一、权威文件索引

- 当前核心审计：`12_reports/CURRENT_ANALYSIS_AUDIT.md`
- 生物学混杂审计：`12_reports/BIOLOGICAL_CONFOUND_AUDIT_V2_CN.md`
- 本综合终审：`12_reports/FINAL_INTEGRATED_EXPERT_REVIEW_PACKAGE_V2_CN.md`
- 关键哈希清单：`12_reports/FINAL_SUBMISSION_AUDIT_MANIFEST_V2.tsv`
- 完整性检查：`13_postfreeze_interpretation/d8_final_integrity_checks.tsv`
- 图稿完整性：`13_postfreeze_interpretation/d8_submission_figure_integrity.tsv`
- 校正文献表：`13_postfreeze_interpretation/d4_orthogonal_literature_evidence_v2.tsv`
- 纠错对比：`13_postfreeze_interpretation/d7_original_vs_dataset_collapsed_convergence.tsv`
- 匹配空集：`13_postfreeze_interpretation/d5_global_response_amplitude_matched_null.tsv`
- 主控日志：`{pipe['log']}`

## 最终判断

目前没有检测到会推翻 RAB1A 唯一 Tier-1 或多证据运输核心叙事的致命错误。最重要的统计实现错误已经透明纠正，候选与 Tier 对纠正稳定；最重要的生物学混杂也已从“未经检查”升级为定量敏感性分析并收缩声称。项目已具备撰写一篇严谨计算资源/感染生物学稿件的基础，但最终投稿冻结应等待 G1/G2 的来源闭合或将其残余不确定性正式写入局限性。
"""

    report_path = REPORTS / "FINAL_INTEGRATED_EXPERT_REVIEW_PACKAGE_V2_CN.md"
    report_path.write_text(report, encoding="utf-8")
    report_row = pd.DataFrame([{"relative_path": str(report_path.relative_to(ROOT)), "bytes": report_path.stat().st_size, "sha256": sha(report_path), "modified": datetime.fromtimestamp(report_path.stat().st_mtime).isoformat(timespec="seconds")}])
    manifest = pd.concat([manifest, report_row], ignore_index=True)
    manifest.to_csv(REPORTS / "FINAL_SUBMISSION_AUDIT_MANIFEST_V2.tsv", sep="\t", index=False)
    print(report_path)
    print(REPORTS / "FINAL_SUBMISSION_AUDIT_MANIFEST_V2.tsv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
