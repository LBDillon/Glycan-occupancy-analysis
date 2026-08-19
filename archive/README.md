# Archive

Superseded code and outputs, kept rather than deleted because they are the
record of how withdrawn numbers were produced. Nothing here is part of the
current pipeline, and nothing here should be quoted.

## runners/

| File | Superseded by | Why |
|---|---|---|
| `build_manifest.py` | `pipeline/04_build_candidate_manifest.py` | built the manifest from `matched_pairs.csv`, so a site had to be matched before anyone asked whether it could be scored — the ordering that let unscoreable sites into balanced matched sets |
| `phase4_analysis.py` | `pipeline/09_analyse_scores.py` | clustered only on the occupied site's ortholog cluster, ignoring control-protein reuse |
| `retention_primary.py` | `pipeline/08_design.py` | targeted run over a matched set that no longer exists; overwrote its own output on resume |
| `retention_analysis.py` | `pipeline/10_analyse_retention_by_class.py` and `10b_analyse_retention_paired.py` | single-set analysis, predating the eukaryotic secretory class |
| `control_features.py` | `pipeline/03b_secretory_features.py` | still current in spirit for the first two control sets, but uninstrumented and superseded by the chunked version |

`score_unmatched.py` was deleted rather than archived: it read a file from
`/tmp`, so it could not be reproduced and has no provenance value.

## Superseded result files

These remain in `results/` under their original names and should be ignored:

- `contrasts_vs_*.csv`, `phase4_primary_analysis.json`, `matching_balance.json`
  — first-round outputs, computed with the defective scorer
- `mpnn_conditional_scores*.csv` — pre-correction scores; the valid rows are
  still used by the retention bridge, the invalid ones are what figure 6 plots

## Superseded documents

- `docs/archive/phase4_primary_result_SUPERSEDED.md` — the first result document.
  Its headline estimate of −0.057 SD is withdrawn: it was computed against a
  reference SD inflated by invalid rows.
