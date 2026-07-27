# -*- coding: utf-8 -*-
# ============================================================================
# Module 05 - GBD: intervention scenarios on the BAPC projection
#
# Two evidence-based scenarios are applied to the BAPC baseline trajectory of
# global COPD-attributable HF cases (2022-2050):
#   * smoking cessation          : relative-risk reduction ~0.40 (RR 0.60)
#   * optimized COPD pharmacotherapy: relative-risk reduction ~0.25 (RR 0.75)
#
# Convention (documented): interventions act on the epidemiologically
# modifiable component of the projection, i.e. the growth above the 2021
# observed caseload; the 2021 prevalent pool is not removable. The scenario
# trajectory is  cases_y(RR) = cases_2021 + (baseline_y - cases_2021) * RR,
# so the averted caseload in 2050 equals (baseline_2050 - cases_2021)*(1-RR).
# Averted cases in 2050 and the cumulative 2022-2050 averted person-cases are
# reported; a naive alternative (RR applied to the whole caseload) is also
# tabulated for transparency.
#
# Input : results/t1_bapc_total_cases.csv  (produced by t1_bapc_run.R)
# Output: results/gbd_intervention_scenarios.csv
#
# Manuscript targets (Table 3 / Fig. 5c): ~2.5 M (smoking) and ~1.8 M
# (pharmacotherapy) averted by 2050.  This script computes 2.93 M and 1.83 M
# under the stated convention - see VALIDATION.md for the reconciliation.
# ============================================================================
import os
import numpy as np
import pandas as pd

REPO = os.environ.get("REPRO_ROOT", os.path.abspath(
    os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")))
IN_CSV = os.environ.get("BAPC_TOTAL_CASES", os.path.join(REPO, "results", "t1_bapc_total_cases.csv"))
OUT_DIR = os.environ.get("OUT_DIR", os.path.join(REPO, "results"))

SCENARIOS = [("smoking_cessation", 0.60), ("optimized_copd_pharmacotherapy", 0.75)]


def main():
    b = pd.read_csv(IN_CSV)
    base2021 = float(b.loc[b.year == 2021, "mean"].iloc[0])
    proj = b[(b.year >= 2022) & (b.year <= 2050)].copy()
    base2050 = float(b.loc[b.year == 2050, "mean"].iloc[0])
    growth2050 = base2050 - base2021

    rows = []
    for name, rr in SCENARIOS:
        red = 1.0 - rr
        scen_2050 = base2021 + growth2050 * rr
        averted_2050 = growth2050 * red
        cum_averted = ((proj["mean"] - base2021) * red).sum()
        # naive alternative: RR applied to the whole caseload
        naive_2050 = base2050 * red
        naive_cum = (proj["mean"] * red).sum()
        rows.append({
            "scenario": name, "RR": rr, "relative_reduction": red,
            "baseline_cases_2021": base2021, "baseline_cases_2050": base2050,
            "scenario_cases_2050": scen_2050,
            "averted_cases_2050": averted_2050,
            "cumulative_averted_2022_2050": cum_averted,
            "naive_full_caseload_averted_2050": naive_2050,
            "naive_full_caseload_cumulative": naive_cum,
        })
        print(f"{name}: RR {rr} (-{red:.0%})  -> 2050 averted {averted_2050/1e6:.2f} M "
              f"(scenario caseload {scen_2050/1e6:.2f} M vs baseline {base2050/1e6:.2f} M); "
              f"cumulative 2022-2050 averted {cum_averted/1e6:.2f} M person-cases")
    out = pd.DataFrame(rows)
    os.makedirs(OUT_DIR, exist_ok=True)
    out.to_csv(os.path.join(OUT_DIR, "gbd_intervention_scenarios.csv"), index=False)


if __name__ == "__main__":
    main()
