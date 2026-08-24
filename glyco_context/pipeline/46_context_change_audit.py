"""Account for every difference between the original and corrected tables.

Not a summary of the diff -- a proof of it. Each changed row is assigned to the
named correction that explains it, and anything left over is reported as
UNEXPLAINED rather than absorbed, because an unattributable change means
something moved that we did not intend to move.

Usage:
    46_context_change_audit.py --old <path to the original context_features.csv>
                               [--new results/datasets/context_features.csv]
                               [--out results/datasets/context_change_audit.csv]
"""
import argparse, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from glyco_context.change_audit import KEY, attribute_changes
from experimental_glycosylation_sites.provenance import hash_file, _git_state

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--old", required=True)
parser.add_argument("--new", default="glyco_context/results/datasets/context_features.csv")
parser.add_argument("--out", default="glyco_context/results/datasets/context_change_audit.csv")
args = parser.parse_args()

old_path, new_path = Path(args.old), Path(args.new)
old = pd.read_csv(old_path, low_memory=False)
new = pd.read_csv(new_path, low_memory=False)
print(f"old: {len(old)} rows, {len(old.columns)} columns  {old_path.resolve()}")
print(f"new: {len(new)} rows, {len(new.columns)} columns  {new_path.resolve()}")

shared_keys = old[KEY].merge(new[KEY], on=KEY, how="inner")
print(f"sites in both tables: {len(shared_keys)}")
only_old = len(old) - len(shared_keys)
only_new = len(new) - len(shared_keys)
if only_old or only_new:
    print(f"  only in old: {only_old}   only in new: {only_new}")

audit = attribute_changes(old, new)
print(f"\nchanged sites: {len(audit)} of {len(shared_keys)}")
print("\nby correction:")
counts = audit.category.value_counts().to_dict() if len(audit) else {}
for name, count in sorted(counts.items(), key=lambda kv: -kv[1]):
    print(f"  {name:34}{count:6d}")

# The silent case, called out because it is the one the triplet check could not
# see: these matched the expected triplet in the original table.
invisible = audit[audit.category == "invisible_gap_jump"] if len(audit) else audit
print(f"\npreviously invisible gap-jumps (matched the triplet, wrong residues): "
      f"{len(invisible)}")
if len(invisible):
    print(invisible[KEY].head(12).to_string(index=False))

unexplained = audit[audit.category == "UNEXPLAINED"] if len(audit) else audit
if len(unexplained):
    print(f"\nUNEXPLAINED changes: {len(unexplained)}")
    print(unexplained.head(12).to_string(index=False))

if len(audit):
    audit.to_csv(args.out, index=False)
    print(f"\naudit -> {Path(args.out).resolve()}")

summary = {
    "old": {str(old_path.resolve()): {"rows": int(len(old)), "sha256": hash_file(old_path)}},
    "new": {str(new_path.resolve()): {"rows": int(len(new)), "sha256": hash_file(new_path)}},
    "sites_in_both": int(len(shared_keys)),
    "changed_sites": int(len(audit)),
    "by_correction": counts,
    "unexplained": int(len(unexplained)),
    "git": _git_state(),
}
record = Path(args.out).with_suffix(".json")
record.write_text(json.dumps(summary, indent=2, default=str))
print(f"summary -> {record.resolve()}")
raise SystemExit(1 if len(unexplained) else 0)
