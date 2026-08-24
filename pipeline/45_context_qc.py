"""QC report for the context table: coverage, continuity, and invariants.

Coverage numbers alone pass a table whose features describe the wrong residue,
so this reports the mapping quality that makes the features interpretable and
then asserts the invariants. It exits non-zero when an invariant is broken --
a QC report that only prints cannot stop a bad table being used.

Usage:
    45_context_qc.py [--features results/datasets/context_features.csv]
                     [--out results/datasets/context_qc.json]
"""
import argparse, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.context_qc import check_invariants
from experimental_glycosylation_sites.context_views import asn_matches, is_core
from experimental_glycosylation_sites.provenance import hash_file, _git_state

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--features", default="results/datasets/context_features.csv")
parser.add_argument("--out", default="results/datasets/context_qc.json")
args = parser.parse_args()

source = Path(args.features)
frame = pd.read_csv(source, low_memory=False)
positions = ("n", "plus1", "plus2")
complete_ss = frame[[f"{p}_dssp_ok" for p in positions]].fillna(False).astype(bool).all(axis=1)
continuous = frame.mapping_continuous.fillna(False).astype(bool)

report = {"rows": int(len(frame)), "by_population": {}}
print(f"{'population':24}{'sites':>7}{'continuous':>12}{'all-3 SS':>10}{'core':>8}{'asn-only':>10}")
core, asn = is_core(frame), asn_matches(frame)
for name, group in frame.groupby("population"):
    idx = group.index
    stats = {
        "sites": int(len(group)),
        "mapping_continuous": int(continuous[idx].sum()),
        "dssp_all_three_positions": int(complete_ss[idx].sum()),
        "triplet_core": int(core[idx].sum()),
        "asn_centred": int(asn[idx].sum()),
    }
    report["by_population"][name] = stats
    print(f"{name:24}{stats['sites']:7d}{stats['mapping_continuous']:12d}"
          f"{stats['dssp_all_three_positions']:10d}{stats['triplet_core']:8d}"
          f"{stats['asn_centred']:10d}")

icodes = frame[[f"{p}_icode" for p in positions]].fillna("").astype(str)
report["sites_with_an_insertion_code"] = int((icodes != "").any(axis=1).sum())
print(f"\nsites with an insertion code at any position: "
      f"{report['sites_with_an_insertion_code']}")

violations = check_invariants(frame)
report["invariant_violations"] = violations
report["source"] = {str(source): {"rows": int(len(frame)), "sha256": hash_file(source)}}
report["git"] = _git_state()

print("\ninvariants:")
if not violations:
    print("  all satisfied")
for violation in violations:
    print(f"  BROKEN {violation['invariant']}: {violation['rows']} rows "
          f"(e.g. {violation['example']})\n         {violation['description']}")

Path(args.out).write_text(json.dumps(report, indent=2, default=str))
print(f"\nreport -> {args.out}")
raise SystemExit(1 if violations else 0)
