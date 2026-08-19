# results/

Rebuilt from `pipeline/`; only this file and `.gitignore` are tracked. Folders
run in pipeline order, so what a stage produces tells you where to look.

| Folder | What it holds | Produced by |
|---|---|---|
| `datasets/` | the site tables themselves — the 4,307-site universe, the 922 occupied sites, each control set, and the structural features computed for them | stages 01–03 |
| `manifests/` | one row per site with the model indices of its three sequon residues, plus the scoreability verdict decided **before** matching | stages 04–05 |
| `matching/` | the matched pairs each comparison rests on, and the balance reports that say whether matching worked | stage 06 |
| `scores/` | **sequon scores** — one row per site, per model. `scores_<set>.csv` | stage 07 |
| `designs/` | **generated sequences summarised as retention** — one row per site, per model | stage 08 |
| `analysis/` | contrasts, confidence intervals, verdicts, significance tests. The end products | stages 09–14 |
| `figures/` | every figure, explained in `docs/figures.md` | stages 20–23 |
| `exploratory/` | superseded outputs and sensitivity branches. **Nothing here should be quoted** | various |

## The four files to look at first

| File | What it says |
|---|---|
| `analysis/analysis_optimal.json` | the primary comparison: 16 pairs, +0.458 SD, inconclusive |
| `analysis/analysis_secretory.json` | the best-powered comparison: 262 pairs, +0.073 SD |
| `analysis/significance.csv` | all eight tests, corrected — none survives |
| `analysis/retention_paired.json` | design retention as a paired contrast, per control set |

## What is in `exploratory/` and why

- **`*_primary.*`** — an earlier label for the internal-control comparison, before
  matching became deterministic. Superseded by `*_optimal.*`.
- **`*_greedy_seed0.*`, `*_k5.*`** — matching sensitivity branches. Real analyses,
  kept for the record, but not the reported result.
- **`contrasts_vs_*.csv`, `phase4_primary_analysis.json`, `matching_balance.json`**
  — first-round outputs computed with the defective scorer. Withdrawn.
- **`mpnn_retention_primary.csv`** — a targeted design run over matched pairs that
  no longer exist after rematching.
- **`retention_analysis.json`** — single-set retention analysis, predating the
  eukaryotic secretory class.
- **`mpnn_retention.csv`** — the working file the main sweep appended to; the
  frozen snapshot in `designs/` is the one analysed.

Two files in `scores/` deserve a note. `mpnn_conditional_scores*.csv` are
**pre-correction** and contain the 105 invalid rows: they are kept because
figure 6 plots those rows, and because their valid rows still feed the
score-versus-retention bridge in figure 4. They are not used for any contrast.
