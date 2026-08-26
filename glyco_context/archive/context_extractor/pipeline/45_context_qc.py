"""QC report for the corrected context table. Gates, rather than describes.

Coverage numbers alone pass a table whose features describe the wrong residue,
so this reports mapping quality alongside coverage and then asserts the
invariants. It exits non-zero when a required condition is violated: a report
that only prints cannot stop a bad table being used.

Usage:
    45_context_qc.py [--features results/datasets/context_features.csv]
                     [--manifest results/datasets/context_manifest.csv]
                     [--views results/datasets] [--out results/datasets/context_qc.json]
"""
import argparse, glob, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from glyco_context.context_features import FEATURE_COLUMNS
from glyco_context.context_merge import KEY
from glyco_context.context_qc import check_invariants
from glyco_context.context_views import (asn_matches,
                                                            exclusion_reason,
                                                            is_core)
from glyco_context.sequence_qc import sequence_context_failures
from experimental_glycosylation_sites.provenance import hash_file, _git_state

POSITIONS = ("n", "plus1", "plus2")


def _json_keys(series) -> dict:
    """groupby(...).size() gives tuple keys; JSON object keys must be strings."""
    return {" | ".join(str(part) for part in key) if isinstance(key, tuple) else str(key): int(value)
            for key, value in series.items()}

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--features", default="glyco_context/results/datasets/context_features.csv")
parser.add_argument("--manifest", default="glyco_context/results/datasets/context_manifest.csv")
parser.add_argument("--views", default="glyco_context/results/datasets")
parser.add_argument("--shard-pattern",
                    default="glyco_context/results/datasets/context_features.shard*.csv")
parser.add_argument("--out", default="glyco_context/results/datasets/context_qc.json")
args = parser.parse_args()

source = Path(args.features)
frame = pd.read_csv(source, low_memory=False)
manifest = pd.read_csv(args.manifest, low_memory=False)
report: dict = {"rows": int(len(frame))}
problems: list = []


def section(title):
    print(f"\n=== {title} ===")


# --- 1. counts -------------------------------------------------------------
section("merged table")
report["population_counts"] = frame.population.value_counts().to_dict()
print(f"rows: {len(frame)}")
print(frame.population.value_counts().to_string())

# --- 2. manifest coverage --------------------------------------------------
section("manifest coverage")
eligible = manifest[manifest.features_available.astype(bool)]
eligible = eligible[eligible.population.isin(frame.population.unique())]
eligible = eligible.dropna(subset=["structure_pdb_id", "structure_chain_id",
                                   "structure_resseq"])[KEY].drop_duplicates()
covered = eligible.merge(frame[KEY].assign(_seen=1), on=KEY, how="left")
absent = int(covered._seen.isna().sum())
report["manifest_eligible_rows"] = int(len(eligible))
report["manifest_rows_absent"] = absent
print(f"eligible manifest rows: {len(eligible)}  absent from table: {absent}")
if absent:
    problems.append(f"{absent} eligible manifest rows missing from the table")

# --- 3. duplicates and failures -------------------------------------------
section("duplicates and extraction failures")
duplicates = int(frame.duplicated(KEY).sum())
failure_rows = 0
for path in glob.glob(args.shard_pattern.replace(".csv", "_failures.csv")):
    if Path(path).stat().st_size > 5:
        failure_rows += len(pd.read_csv(path, low_memory=False))
report["duplicate_keys"] = duplicates
report["extraction_failures"] = failure_rows
print(f"duplicate context keys: {duplicates}\nrecorded failures: {failure_rows}")
if duplicates:
    problems.append(f"{duplicates} duplicate context keys")
if failure_rows:
    problems.append(f"{failure_rows} recorded extraction failures")

# --- 4. sequence completeness ---------------------------------------------
section("expected sequence completeness")
joined = frame.merge(
    manifest[KEY + [c for c in ("uniprot_length", "sequon_triplet") if c in manifest]]
    .drop_duplicates(KEY), on=KEY, how="left", suffixes=("", "_m"))
if "sequon_triplet" not in joined.columns and "triplet_expected" in joined.columns:
    joined["sequon_triplet"] = joined.triplet_expected
sequence_failures = sequence_context_failures(joined)
report["sequence_failures"] = (
    _json_keys(sequence_failures.groupby(["population", "reason"]).size())
    if len(sequence_failures) else {})
if len(sequence_failures):
    print(sequence_failures.groupby(["population", "reason"]).size().to_string())
    problems.append(f"{len(sequence_failures)} rows without usable sequence context")
else:
    print("every row: sequence present, complete N-X-S/T triplet (X != P), valid coordinates")

# --- 5. continuity, triplet agreement, per-position coverage --------------
section("mapping quality and coverage by population")
continuous = frame.mapping_continuous.fillna(False).astype(bool)
matches = frame.triplet_matches.fillna(False).astype(bool)
header = (f"{'population':24}{'sites':>7}{'contin.':>9}{'triplet=':>9}"
          + "".join(f"{p+' RSA':>10}{p+' SS':>9}" for p in POSITIONS))
print(header)
by_population = {}
for name, group in frame.groupby("population"):
    idx = group.index
    stats = {"sites": int(len(group)),
             "mapping_continuous": int(continuous[idx].sum()),
             "triplet_agrees": int(matches[idx].sum())}
    line = f"{name:24}{stats['sites']:7d}{stats['mapping_continuous']:9d}{stats['triplet_agrees']:9d}"
    for position in POSITIONS:
        rsa = int(group[f"{position}_rsa"].notna().sum())
        ss = int(group[f"{position}_ss"].notna().sum())
        stats[f"{position}_rsa"] = rsa
        stats[f"{position}_ss"] = ss
        line += f"{rsa:10d}{ss:9d}"
    by_population[name] = stats
    print(line)
report["by_population"] = by_population

# --- 6. primary panel coverage --------------------------------------------
section("primary panel coverage by population (non-null %)")
panel = [c for c in FEATURE_COLUMNS if c in frame.columns]
coverage = {}
for name, group in frame.groupby("population"):
    coverage[name] = {c: round(100 * float(group[c].notna().mean()), 1) for c in panel}
report["panel_coverage_percent"] = coverage
thin = {(p, c): v for p, columns in coverage.items()
        for c, v in columns.items() if v < 50.0}
for name in coverage:
    worst = sorted(coverage[name].items(), key=lambda kv: kv[1])[:4]
    print(f"  {name:24} lowest: " + ", ".join(f"{c}={v}%" for c, v in worst))
report["panel_features_below_50_percent"] = {f"{p}|{c}": v for (p, c), v in thin.items()}

# --- 7. invariants ---------------------------------------------------------
section("invariants")
violations = check_invariants(frame)
report["invariant_violations"] = violations
if not violations:
    print("all satisfied")
for violation in violations:
    print(f"  BROKEN {violation['invariant']}: {violation['rows']} rows "
          f"(e.g. {violation['example']})")
    problems.append(f"invariant {violation['invariant']} broken on {violation['rows']} rows")

# --- 8. views and exclusions ----------------------------------------------
section("derived views")
view_sizes = {}
for name in ("triplet_core", "asn_core", "construct_review"):
    path = Path(args.views) / f"context_{name}.csv"
    if path.exists():
        view = pd.read_csv(path, low_memory=False)
        view_sizes[name] = {"rows": int(len(view)),
                            "sha256": hash_file(path),
                            "path": str(path.resolve())}
        print(f"  {name:18}{len(view):7d}  {path.resolve()}")
report["views"] = view_sizes
computed = {"triplet_core": int(is_core(frame).sum()),
            "asn_core": int(asn_matches(frame).sum())}
for name, expected in computed.items():
    actual = view_sizes.get(name, {}).get("rows")
    if actual is not None and actual != expected:
        problems.append(f"view {name} has {actual} rows, expected {expected}")

section("exclusion reasons from triplet_core")
reasons = exclusion_reason(frame)
table = frame.assign(exclusion_reason=reasons)
excluded = table[table.exclusion_reason.ne("")]
report["exclusion_reasons"] = (
    _json_keys(excluded.groupby(["exclusion_reason", "population"]).size())
    if len(excluded) else {})
print(excluded.groupby(["exclusion_reason", "population"]).size().to_string()
      if len(excluded) else "none")

# --- 9. provenance ---------------------------------------------------------
report["inputs"] = {
    str(source.resolve()): {"rows": int(len(frame)), "sha256": hash_file(source)},
    str(Path(args.manifest).resolve()): {"sha256": hash_file(Path(args.manifest))},
}
report["git"] = _git_state()
report["problems"] = problems

Path(args.out).write_text(json.dumps(report, indent=2, default=str))
section("result")
print(f"report -> {Path(args.out).resolve()}")
if problems:
    for problem in problems:
        print(f"  FAIL {problem}")
    raise SystemExit(1)
print("all gates passed")
