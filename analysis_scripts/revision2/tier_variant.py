"""Revision-2 discovery-only, parameterised re-implementation of the frozen evidence-tier rule.

Purpose
-------
Provide a *config-driven* version of ``tier()`` from
``scripts/stage6c_integrated_evidence_atlas.py`` so that the candidate rule can be
relaxed/tightened one requirement at a time (WP4) without ever touching the frozen
original script.

Hard boundaries (enforced, not merely documented)
-------------------------------------------------
1. ``select_candidates()`` is discovery-only. It reads the pre-holdout evidence
   atlas and recomputes candidate eligibility from raw discovery fields. It does
   **not** take ``evidence_tier`` as a decision input; that column is used only by
   the parity auditor, as a reference.
2. Every holdout-derived column is rejected at load time (see ``FORBIDDEN_COLUMNS``).
   In particular ``final_tier`` / ``final_candidate_tier_matrix.tsv`` must never be
   read here.
3. Post-holdout classification (``stage8_final_tier_matrix.assign``) is deliberately
   **not** implemented in this module.

Note on ``domain_count``
------------------------
The frozen ``domain_count`` in ``stage8_final_tier_matrix.py`` sums seven columns,
two of which are holdout-derived (``holdout_directional_support``,
``holdout_multieffector_support``). It is therefore holdout-contaminated and is not
reused. ``discovery_domain_count()`` below counts five discovery domains only.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Iterable, Sequence

import pandas as pd

# --------------------------------------------------------------------------- #
# Column contracts
# --------------------------------------------------------------------------- #

REQUIRED_DISCOVERY_COLUMNS: tuple[str, ...] = (
    "human_symbol",
    "functional_tier",
    "ppi_AB_flag",
    "apms_high_confidence",
    "fitness_artifact_compatible",
    "stable_cross_cohort",
    "any_mechanism_convergent",
    "lcv_spatial_supported",
)

#: Any column containing one of these substrings is holdout-derived and forbidden.
FORBIDDEN_COLUMN_SUBSTRINGS: tuple[str, ...] = (
    "holdout",
    "final_tier",
    "final_claim_boundary",
    "domain_count",
)

PANEL_TIER_PREFIXES: tuple[str, ...] = ("TIER1", "TIER2")

TIER1_MULTIMODAL = "TIER1_MULTIMODAL"
TIER2_DIRECT_FUNCTIONAL = "TIER2_DIRECT_FUNCTIONAL"
TIER2_DIRECT_ORTHOGONAL = "TIER2_DIRECT_ORTHOGONAL"
TIER2_FUNCTIONAL_RESPONSE_MECHANISM = "TIER2_FUNCTIONAL_RESPONSE_MECHANISM"
TIER3_FITNESS_CAVEAT = "TIER3_FITNESS_CAVEAT"
TIER3_SINGLE_ROUTE = "TIER3_SINGLE_ROUTE"
BACKGROUND = "BACKGROUND"

F1_F2: tuple[str, ...] = (
    "F1_GW_AND_BOTH_BATCH_CONCORDANT",
    "F2_MULTILAYER_CONCORDANT",
)
F3 = "F3_SINGLE_LAYER_OR_DISCORDANT"


# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #

#: Orthogonal-evidence modes accepted by the Tier-1 and Tier-2 branches.
ORTHOGONAL_MODES: tuple[str, ...] = (
    "any",                       # response OR mechanism OR spatial  (frozen behaviour)
    "response_only",             # response
    "no_spatial",                # response OR mechanism
    "response_and_mechanism",    # response AND mechanism
    "response_mechanism_spatial",  # response AND mechanism AND spatial
)


@dataclass(frozen=True)
class TierConfig:
    """One complete, explicit candidate-rule configuration.

    ``ORIGINAL`` reproduces ``stage6c_integrated_evidence_atlas.tier`` exactly.
    """

    name: str = "ORIGINAL"

    # --- functional branch -------------------------------------------------
    functional_set: tuple[str, ...] = F1_F2
    fitness_blocks_functional: bool = True

    # --- orthogonal branches ----------------------------------------------
    spatial_in_orthogonal: bool = True
    t1_orthogonal: str = "any"
    t2_direct_orthogonal: str = "any"
    t1_min_discovery_domains: int = 0

    # --- tier-structure switches ------------------------------------------
    t2_functional_response_branch: bool = True
    tier3_single_route: bool = True
    direct_only: bool = False

    # --- panel admission override ------------------------------------------
    #: Global minimum number of discovery domains required for *panel admission*.
    #: Applied in ``select_candidates`` after tier assignment, so tier labels stay
    #: interpretable. This is an admission rule, not a tier relabelling, and it is
    #: the quantity manipulated by the domain-count sensitivity series.
    panel_min_discovery_domains: int = 0

    # --- bookkeeping ------------------------------------------------------
    notes: str = ""
    changed_requirement: str = "none"

    def __post_init__(self) -> None:
        for mode_field in ("t1_orthogonal", "t2_direct_orthogonal"):
            value = getattr(self, mode_field)
            if value not in ORTHOGONAL_MODES:
                raise ValueError(
                    f"{mode_field}={value!r} not in {ORTHOGONAL_MODES}"
                )
        if self.t1_min_discovery_domains < 0 or self.t1_min_discovery_domains > 5:
            raise ValueError("t1_min_discovery_domains must be within 0..5")
        if self.panel_min_discovery_domains < 0 or self.panel_min_discovery_domains > 5:
            raise ValueError("panel_min_discovery_domains must be within 0..5")


ORIGINAL = TierConfig()

#: Group A — single-condition perturbations (one requirement changed at a time).
GROUP_A: dict[str, TierConfig] = {
    "A1_admit_F3": replace(
        ORIGINAL,
        name="A1_admit_F3",
        functional_set=F1_F2 + (F3,),
        changed_requirement="functional set includes F3 single-layer/discordant",
    ),
    "A2_fitness_does_not_block": replace(
        ORIGINAL,
        name="A2_fitness_does_not_block",
        fitness_blocks_functional=False,
        changed_requirement="baseline-fitness caveat no longer blocks functional evidence",
    ),
    "A3_T1_orthogonal_response_only": replace(
        ORIGINAL,
        name="A3_T1_orthogonal_response_only",
        t1_orthogonal="response_only",
        changed_requirement=(
            "Tier 1 orthogonal evidence narrowed to expression response only; "
            "mechanism and LCV spatial evidence no longer satisfy the requirement"
        ),
    ),
    "A4_T2_orthogonal_no_spatial": replace(
        ORIGINAL,
        name="A4_T2_orthogonal_no_spatial",
        t2_direct_orthogonal="no_spatial",
        changed_requirement="Tier 2 direct-orthogonal ignores spatial evidence",
    ),
    "A5_T1_require_response_and_mechanism": replace(
        ORIGINAL,
        name="A5_T1_require_response_and_mechanism",
        t1_orthogonal="response_and_mechanism",
        changed_requirement="Tier 1 requires response AND mechanism",
    ),
    "A6_T2_orthogonal_require_response_and_mechanism": replace(
        ORIGINAL,
        name="A6_T2_orthogonal_require_response_and_mechanism",
        t2_direct_orthogonal="response_and_mechanism",
        changed_requirement="Tier 2 direct-orthogonal requires response AND mechanism",
    ),
    "A7_drop_T2_functional_response_branch": replace(
        ORIGINAL,
        name="A7_drop_T2_functional_response_branch",
        t2_functional_response_branch=False,
        changed_requirement="Tier 2 functional-response-mechanism branch removed",
    ),
    "A8_T1_require_four_domains": replace(
        ORIGINAL,
        name="A8_T1_require_four_domains",
        t1_min_discovery_domains=4,
        changed_requirement="Tier 1 additionally requires >=4 discovery domains",
    ),
}

#: Variants that were specified, then measured and found to be mathematical no-ops.
#: They are recorded here for the audit trail but are **not** part of the reported
#: variant set, because a no-op carries no information about rule sensitivity.
#:
#: ``A8_T1_require_three_domains``: Tier 1 already requires direct evidence,
#: fitness-filtered functional evidence and at least one orthogonal domain, i.e. a
#: minimum of three discovery domains. Requiring ``popcount(domains) >= 3`` is
#: therefore implied by the rule and changes nothing: 0 rows relabelled, panel
#: 56 -> 56. Detected by ``run_variant_effect_audit.py`` and removed.
REMOVED_NOOP: dict[str, TierConfig] = {
    "A8_T1_require_three_domains": replace(
        ORIGINAL,
        name="A8_T1_require_three_domains",
        t1_min_discovery_domains=3,
        changed_requirement="Tier 1 additionally requires >=3 discovery domains (implied by Tier 1; no-op)",
        notes="REMOVED FROM FORMAL WP4: measured no-op, 0 rows relabelled, panel 56 -> 56",
    ),
}

#: Group B — structural ablations (reported separately, never pooled with Group A).
GROUP_B: dict[str, TierConfig] = {
    "B1_direct_only": replace(
        ORIGINAL,
        name="B1_direct_only",
        direct_only=True,
        changed_requirement="panel restricted to direct physical evidence",
    ),
    "B2_drop_T3_single_route": replace(
        ORIGINAL,
        name="B2_drop_T3_single_route",
        tier3_single_route=False,
        changed_requirement="Tier 3 single-route class removed",
    ),
}

ALL_VARIANTS: dict[str, TierConfig] = {"ORIGINAL": ORIGINAL, **GROUP_A, **GROUP_B}


#: Domain-count sensitivity series.
#:
#: This is **not** a symmetric threshold perturbation of a single rule. It imposes a
#: global minimum number of discovery domains on *panel admission* and reports the
#: whole series, because picking any single cut-off would invite the objection that
#: the cut-off was chosen to protect the result. Report as:
#: ``Domain-count sensitivity analysis, not a symmetric threshold perturbation``.
DOMAIN_SENSITIVITY: dict[str, TierConfig] = {
    f"D{cut}_panel_min_discovery_domains": replace(
        ORIGINAL,
        name=f"D{cut}_panel_min_discovery_domains",
        panel_min_discovery_domains=cut,
        changed_requirement=f"panel admission requires >={cut} discovery domains",
    )
    for cut in (3, 4, 5)
}


#: Declared effect of each variant, stated **before** measurement from the rule text
#: alone. ``run_variant_effect_audit.py`` compares this against the measured effect and
#: fails on any mismatch, so that no variant can be re-described after the fact.
#:
#: ``panel`` is one of ``expand`` / ``contract`` / ``unchanged``.
#: ``rows`` is one of ``some`` / ``none``.
DECLARED_EFFECT: dict[str, dict[str, str]] = {
    "A1_admit_F3": {
        "panel": "expand",
        "mode": "tier_relabel",
        "reason": "F3 single-layer/discordant functional evidence newly qualifies genes",
    },
    "A2_fitness_does_not_block": {
        "panel": "expand",
        "mode": "tier_relabel",
        "reason": "fitness-caveated CRISPR evidence starts counting as clean functional evidence",
    },
    "A3_T1_orthogonal_response_only": {
        "panel": "unchanged",
        "mode": "tier_relabel",
        "reason": (
            "demoted Tier 1 genes satisfy direct + functional_clean and therefore land in "
            "TIER2_DIRECT_FUNCTIONAL; panel membership cannot change by construction"
        ),
    },
    "A4_T2_orthogonal_no_spatial": {
        "panel": "contract",
        "mode": "tier_relabel",
        "reason": "Tier 2 direct-orthogonal genes supported only by LCV spatial evidence lose their only panel branch",
    },
    "A5_T1_require_response_and_mechanism": {
        "panel": "unchanged",
        "mode": "tier_relabel",
        "reason": "demoted Tier 1 genes land in TIER2_DIRECT_FUNCTIONAL",
    },
    "A6_T2_orthogonal_require_response_and_mechanism": {
        "panel": "contract",
        "mode": "tier_relabel",
        "reason": "Tier 2 direct-orthogonal genes with only one orthogonal domain lose their only panel branch",
    },
    "A7_drop_T2_functional_response_branch": {
        "panel": "contract",
        "mode": "tier_relabel",
        "reason": "genes admitted only by the functional-response-mechanism branch have no direct evidence and drop out",
    },
    "A8_T1_require_four_domains": {
        "panel": "unchanged",
        "mode": "tier_relabel",
        "reason": "demoted Tier 1 genes land in TIER2_DIRECT_FUNCTIONAL; this is a tier-architecture check",
    },
    "B1_direct_only": {
        "panel": "contract",
        "mode": "admission_filter",
        "reason": "panel restricted to genes with direct physical evidence; admission filter, no tier relabelling",
    },
    "B2_drop_T3_single_route": {
        "panel": "unchanged",
        "mode": "tier_relabel",
        "reason": "TIER3_SINGLE_ROUTE is outside the panel, so relabelling it to BACKGROUND cannot move panel members",
    },
    "D3_panel_min_discovery_domains": {
        "panel": "contract",
        "mode": "admission_filter",
        "reason": "admission gate removes panel genes with fewer than 3 discovery domains; no tier relabelling",
    },
    "D4_panel_min_discovery_domains": {
        "panel": "contract",
        "mode": "admission_filter",
        "reason": "admission gate removes panel genes with fewer than 4 discovery domains; no tier relabelling",
    },
    "D5_panel_min_discovery_domains": {
        "panel": "contract",
        "mode": "admission_filter",
        "reason": "admission gate removes panel genes with fewer than 5 discovery domains; no tier relabelling",
    },
}

ALL_CONFIGS: dict[str, TierConfig] = {**ALL_VARIANTS, **DOMAIN_SENSITIVITY}


#: Declared reporting class of each variant. WP4 has three classes and they must not be
#: pooled: ``panel_effect`` asks whether panel *membership* moves, ``tier_only`` asks
#: whether tier and primary assignment move while membership is fixed, and ``noop``
#: carries no information and is excluded from formal results.
VARIANT_CLASS: dict[str, str] = {
    "A1_admit_F3": "panel_effect",
    "A2_fitness_does_not_block": "panel_effect",
    "A3_T1_orthogonal_response_only": "tier_only",
    "A4_T2_orthogonal_no_spatial": "panel_effect",
    "A5_T1_require_response_and_mechanism": "tier_only",
    "A6_T2_orthogonal_require_response_and_mechanism": "panel_effect",
    "A7_drop_T2_functional_response_branch": "panel_effect",
    "A8_T1_require_four_domains": "tier_only",
    "B1_direct_only": "panel_effect",
    "B2_drop_T3_single_route": "tier_only",
    "D3_panel_min_discovery_domains": "panel_effect",
    "D4_panel_min_discovery_domains": "panel_effect",
    "D5_panel_min_discovery_domains": "panel_effect",
    "A8_T1_require_three_domains": "noop",
}

PANEL_EFFECT_VARIANTS: tuple[str, ...] = tuple(
    name for name, klass in VARIANT_CLASS.items() if klass == "panel_effect"
)
TIER_ONLY_VARIANTS: tuple[str, ...] = tuple(
    name for name, klass in VARIANT_CLASS.items() if klass == "tier_only"
)


# --------------------------------------------------------------------------- #
# Row-level evidence
# --------------------------------------------------------------------------- #

def as_bool(value: object) -> bool:
    """Parse the many boolean spellings present in the atlas TSV."""
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    text = str(value).strip().lower()
    if text in {"true", "1", "1.0", "yes", "y", "t"}:
        return True
    if text in {"false", "0", "0.0", "no", "n", "f", "", "nan", "none", "nat"}:
        return False
    return False


@dataclass(frozen=True)
class DiscoveryEvidence:
    """Exactly the discovery fields the rule is allowed to see."""

    gene: str
    functional_tier: str
    direct: bool
    fitness_artifact_compatible: bool
    response: bool
    mechanism: bool
    spatial: bool
    functional: bool
    functional_clean: bool
    ppi_AB_only: bool = False
    translatome: bool = False

    def discovery_domain_count(self) -> int:
        """Discovery domains only (0..5). Never includes holdout evidence."""
        return int(
            self.direct
            + self.functional_clean
            + self.response
            + self.mechanism
            + self.spatial
        )


def build_evidence(row: pd.Series, cfg: TierConfig = ORIGINAL) -> DiscoveryEvidence:
    functional_tier = str(row.get("functional_tier", "")).strip()
    functional = functional_tier in cfg.functional_set
    fitness = as_bool(row.get("fitness_artifact_compatible"))
    functional_clean = functional and not (cfg.fitness_blocks_functional and fitness)
    return DiscoveryEvidence(
        gene=str(row.get("human_symbol", "")),
        functional_tier=functional_tier,
        direct=as_bool(row.get("ppi_AB_flag")) or as_bool(row.get("apms_high_confidence")),
        fitness_artifact_compatible=fitness,
        response=as_bool(row.get("stable_cross_cohort")),
        mechanism=as_bool(row.get("any_mechanism_convergent")),
        spatial=as_bool(row.get("lcv_spatial_supported")),
        functional=functional,
        functional_clean=functional_clean,
        ppi_AB_only=as_bool(row.get("ppi_AB_flag")),
        translatome=(
            as_bool(row.get("mechanism_translation_effector_module_convergent"))
            or as_bool(row.get("mechanism_rna_rpf_interaction_convergent"))
        ),
    )


def _orthogonal(ev: DiscoveryEvidence, mode: str, spatial_allowed: bool) -> bool:
    response, mechanism = ev.response, ev.mechanism
    spatial = ev.spatial and spatial_allowed
    if mode == "any":
        return bool(response or mechanism or spatial)
    if mode == "response_only":
        return bool(response)
    if mode == "no_spatial":
        return bool(response or mechanism)
    if mode == "response_and_mechanism":
        return bool(response and mechanism)
    if mode == "response_mechanism_spatial":
        return bool(response and mechanism and spatial)
    raise ValueError(f"unknown orthogonal mode {mode!r}")


# --------------------------------------------------------------------------- #
# The parameterised rule
# --------------------------------------------------------------------------- #

def tier_variant(ev: DiscoveryEvidence, cfg: TierConfig = ORIGINAL) -> str:
    """Return the discovery-side evidence tier for one gene under ``cfg``.

    With ``cfg is ORIGINAL`` this is a literal re-implementation of the frozen
    ``tier()`` function. No holdout information is used at any point.
    """
    direct = ev.direct
    functional = ev.functional
    functional_clean = ev.functional_clean
    response = ev.response
    mechanism = ev.mechanism

    t1_orth = _orthogonal(ev, cfg.t1_orthogonal, cfg.spatial_in_orthogonal)
    t2_orth = _orthogonal(ev, cfg.t2_direct_orthogonal, cfg.spatial_in_orthogonal)
    t1_domains_ok = ev.discovery_domain_count() >= cfg.t1_min_discovery_domains

    if direct and functional_clean and t1_orth and t1_domains_ok:
        return TIER1_MULTIMODAL
    if direct and functional_clean:
        return TIER2_DIRECT_FUNCTIONAL
    if cfg.t2_functional_response_branch and functional_clean and response and mechanism:
        return TIER2_FUNCTIONAL_RESPONSE_MECHANISM
    if direct and t2_orth:
        return TIER2_DIRECT_ORTHOGONAL
    if functional and ev.fitness_artifact_compatible:
        return TIER3_FITNESS_CAVEAT
    single_route = (
        functional_clean
        or (response and mechanism)
        or direct
        or (ev.spatial and cfg.spatial_in_orthogonal and (response or mechanism))
    )
    if cfg.tier3_single_route and single_route:
        return TIER3_SINGLE_ROUTE
    return BACKGROUND


def entry_route_flags(ev: DiscoveryEvidence, cfg: TierConfig = ORIGINAL) -> dict[str, bool]:
    """Code-level route components.

    The frozen pipeline stores **no per-gene route label**. ``stage6d`` selects the panel
    purely by ``evidence_tier ∈ {TIER1*, TIER2*}``; the "two frozen routes" in the
    manuscript are described as requirements, but no script evaluates them as a filter.
    Route must therefore be *derived*, and the derivation must be stated explicitly.

    Three operational definitions are reported separately so that the manuscript's prose
    can be checked against what the code actually does:

    ``route1_strict``
        Class A/B physical evidence (curated PPI A/B only) + fitness-filtered CRISPR
        functional evidence + at least one expression/translation response. This is the
        closest operational reading of the Methods sentence for Route 1.
    ``route1_broad``
        The branch actually implemented by ``tier()`` for Tier 1: any direct physical
        evidence (curated A/B or high-confidence AP-MS) + fitness-filtered CRISPR
        functional evidence + any orthogonal evidence (response, mechanism or LCV
        spatial).
    ``route2``
        Fitness-filtered CRISPR functional evidence + response + mechanism, i.e. the
        ``TIER2_FUNCTIONAL_RESPONSE_MECHANISM`` branch.
    """
    translatome = as_bool(ev.__dict__.get("translatome", False))
    return {
        "route1_strict": bool(
            as_bool(ev.__dict__.get("ppi_AB_only", False))
            and ev.functional_clean
            and (ev.response or translatome)
        ),
        "route1_broad": bool(ev.direct and ev.functional_clean and (ev.response or ev.mechanism or ev.spatial)),
        "route2": bool(ev.functional_clean and ev.response and ev.mechanism),
        "direct_evidence": bool(ev.direct),
        "functional_clean": bool(ev.functional_clean),
        "fitness_artifact_compatible": bool(ev.fitness_artifact_compatible),
        "discovery_domain_count": ev.discovery_domain_count(),
    }


# --------------------------------------------------------------------------- #
# Table-level entry point
# --------------------------------------------------------------------------- #

def assert_discovery_only(columns: Iterable[str]) -> None:
    """Raise if any holdout-derived column is present in the input table."""
    offending = sorted(
        {
            str(c)
            for c in columns
            for token in FORBIDDEN_COLUMN_SUBSTRINGS
            if token in str(c)
        }
    )
    if offending:
        raise ValueError(
            "Holdout-derived column(s) present in a discovery-only table: "
            + ", ".join(offending)
        )


def load_discovery_atlas(path: str, columns: Sequence[str] | None = None) -> pd.DataFrame:
    """Load the pre-holdout atlas, refusing any holdout-derived column."""
    frame = pd.read_csv(path, sep="\t", low_memory=False,
                        usecols=list(columns) if columns else None)
    assert_discovery_only(frame.columns)
    missing = [c for c in REQUIRED_DISCOVERY_COLUMNS if c not in frame.columns]
    if missing:
        raise KeyError("Missing required discovery column(s): " + ", ".join(missing))
    return frame


#: All columns added by :func:`annotate` carry this prefix so that they can never be
#: confused with, or silently alias, an atlas column of the same name.
OUTPUT_PREFIX = "tv_"


def annotate(frame: pd.DataFrame, cfg: TierConfig = ORIGINAL) -> pd.DataFrame:
    """Attach ``variant_tier`` and prefixed route flags for every row of ``frame``."""
    rows = []
    for _, row in frame.iterrows():
        ev = build_evidence(row, cfg)
        item = {
            "variant_tier": tier_variant(ev, cfg),
            **entry_route_flags(ev, cfg),
        }
        rows.append(item)
    out = pd.DataFrame(rows, index=frame.index)
    out = out.rename(columns={c: f"{OUTPUT_PREFIX}{c}" for c in out.columns if c != "variant_tier"})
    clash = [c for c in out.columns if c in frame.columns]
    if clash:
        raise ValueError(
            "Variant output would overwrite discovery column(s): " + ", ".join(clash)
        )
    return pd.concat([frame, out], axis=1)


def select_candidates(frame: pd.DataFrame, cfg: TierConfig = ORIGINAL) -> pd.DataFrame:
    """Return the candidate panel implied by ``cfg`` (TIER1* ∪ TIER2*).

    Eligibility is recomputed from discovery fields; the frozen ``evidence_tier``
    column is never consulted.

    ``cfg.panel_min_discovery_domains`` is applied **after** tier assignment, as an
    admission gate. Tier labels themselves remain those of ``tier_variant`` so that the
    domain-count sensitivity series stays interpretable as "who is admitted", not as a
    second, silent tier rule.
    """
    annotated = annotate(frame, cfg)
    panel = annotated.loc[
        annotated.variant_tier.str.startswith(PANEL_TIER_PREFIXES, na=False)
    ].copy()
    if cfg.panel_min_discovery_domains > 0:
        panel = panel.loc[
            panel[f"{OUTPUT_PREFIX}discovery_domain_count"].astype(int)
            >= cfg.panel_min_discovery_domains
        ].copy()
    if cfg.direct_only:
        panel = panel.loc[panel[f"{OUTPUT_PREFIX}direct_evidence"].astype(bool)].copy()
    return panel


def candidate_genes(frame: pd.DataFrame, cfg: TierConfig = ORIGINAL) -> set[str]:
    panel = select_candidates(frame, cfg)
    return {str(g).upper() for g in panel.human_symbol}
