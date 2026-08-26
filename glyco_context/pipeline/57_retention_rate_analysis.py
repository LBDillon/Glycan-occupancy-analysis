"""Report the sequon retention rate against its controls, with provenance.

Separate from stage 56 so the statistics can be recomputed without regenerating
designs, and so they live in tested code rather than in whatever was typed at a
prompt. The aggregation is in `glyco_context.retention_rate`.

Usage:
    57_retention_rate_analysis.py [--rows ...] [--boot 4000]
"""
import argparse, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from glyco_context.retention_rate import QUANTITIES, summarise
from experimental_glycosylation_sites.provenance import hash_file, _git_state

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--rows", default="glyco_context/results/analysis/sequon_retention_rate.csv")
parser.add_argument("--dropped", default="glyco_context/results/analysis/sequon_retention_rate_dropped.csv")
parser.add_argument("--boot", type=int, default=4000)
parser.add_argument("--out", default="glyco_context/results/analysis/sequon_retention_rate_summary.json")
args = parser.parse_args()

source = Path(args.rows)
rows = pd.read_csv(source, low_memory=False)
summary = summarise(rows, n_boot=args.boot)

print(f"{summary['design_rows']} design-site rows | {summary['sites']} sites | "
      f"{summary['proteins']} proteins\n")
print(f"{'quantity':40}{'mean':>9}{'95% CI':>21}")
for column, label in QUANTITIES.items():
    if column not in summary:
        continue
    q = summary[column]
    print(f"{label:40}{100*q['mean']:8.1f}%  [{100*q['ci_low']:6.1f}%,{100*q['ci_high']:6.1f}%]")

if "control_minus_sequon" in summary:
    q = summary["control_minus_sequon"]
    print(f"\n{'control minus sequon (exact)':40}{100*q['mean']:8.1f}pp"
          f"  [{100*q['ci_low']:6.1f} ,{100*q['ci_high']:6.1f} ]")
    print("\nAn interval spanning zero is no detectable excess loss. It is not a\n"
          "demonstration that the rates are equal: this one remains compatible with\n"
          f"an excess as large as {100*q['ci_high']:.1f} percentage points.")

# Coverage: what was dropped matters as much as what was kept.
dropped = Path(args.dropped)
summary["dropped_chains"] = (
    pd.read_csv(dropped).reason.value_counts().to_dict() if dropped.exists() else {})
summary["source"] = {str(source.resolve()): hash_file(source)}
summary["n_boot"] = args.boot
summary["git"] = _git_state()
Path(args.out).write_text(json.dumps(summary, indent=2, default=str))
print(f"\nwrote {Path(args.out).resolve()}")
