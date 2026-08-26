"""What survives hiding the motif: the masking comparison, for any model.

Both arms of every comparison contain the sequon, so motif recognition alone
cannot separate them. This asks what the surroundings say once the motif is
hidden, as a paired change in contrast between the two masking schemes.

This existed nowhere in the repository: the ESMC masking result was computed by
hand and only the answer was kept, so the benchmark's sharpest finding could not
be regenerated.

**It is reported in log odds, the scale the rest of the analysis uses.** The
figure quoted in the earlier summaries -- +0.12 falling to +0.007, a change of
+0.113 -- is in predicted-probability units, which appear nowhere else. Both
describe the same collapse; they are not the same number and should not be
compared. On this scale the same ESMC comparison reads +0.404 falling to -0.155.

Usage:
    15_masking_comparison.py --visible esmc_single --hidden esmc_joint
                             [--label secretory] [--boot 20000]
"""
import argparse, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites import analysis_paths as paths
from experimental_glycosylation_sites.masking import masking_change
from experimental_glycosylation_sites.provenance import hash_file, _git_state

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--visible", required=True, help="variant scored with the motif visible")
parser.add_argument("--hidden", required=True, help="variant scored with the motif hidden")
parser.add_argument("--label", default="secretory")
parser.add_argument("--boot", type=int, default=20000)
parser.add_argument("--out", default=None)
args = parser.parse_args()

for role, variant in (("visible", args.visible), ("hidden", args.hidden)):
    path = paths.contrasts(args.label, variant)
    if not path.exists():
        raise SystemExit(
            f"{role} contrasts not found: {path}\n"
            f"Run 09_analyse_scores.py {args.label} --variant {variant} first.")

visible_path = paths.contrasts(args.label, args.visible)
hidden_path = paths.contrasts(args.label, args.hidden)
result = masking_change(pd.read_csv(visible_path), pd.read_csv(hidden_path),
                        n_boot=args.boot)

print(f"masking comparison, {args.label}: {args.visible} against {args.hidden}")
print(f"  pairs {result['n_pairs']}   resample units {result['n_units']}\n")
print(f"  contrast, motif visible : {result['mean_visible']:+.4f}")
print(f"  contrast, motif hidden  : {result['mean_hidden']:+.4f}")
print(f"  change                  : {result['mean']:+.4f}"
      f"  [{result['ci_low']:+.4f}, {result['ci_high']:+.4f}]   p={result['p']:.4f}")
print("\n  A large positive change means the preference depended on seeing the motif"
      "\n  in context. A change near zero means it did not.")

out = Path(args.out) if args.out else Path(
    f"results/analysis/masking_{args.label}_{args.visible}_vs_{args.hidden}.json")
out.parent.mkdir(parents=True, exist_ok=True)
result.update({
    "label": args.label, "visible": args.visible, "hidden": args.hidden,
    "n_boot": args.boot,
    "inputs": {str(visible_path.resolve()): hash_file(visible_path),
               str(hidden_path.resolve()): hash_file(hidden_path)},
    "git": _git_state()})
out.write_text(json.dumps(result, indent=2, default=str))
print(f"\nwrote {out.resolve()}")
