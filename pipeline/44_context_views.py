"""Build the three analysis-ready views of the context table.

One table cannot serve the atlas: a crystallographer's N->Q knockout, a +1
differing between isoform and construct, and a +2 that was never resolved are
three different facts. See context_views for what each view means and, for
`asn_core`, what it may not be used for.

Sequence context is joined from the manifest rather than recomputed here,
because it is a property of the UniProt sequence: a structure that stops short
of the C-terminus must not shorten it. Row-level sequence checks run *before*
the split, so a partially loaded cache fails here rather than producing views
whose sequence columns are quietly empty.

Usage:
    44_context_views.py [--features results/datasets/context_features.csv]
                        [--manifest results/datasets/context_manifest.csv]
                        [--outdir results/datasets]
"""
import argparse, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.context_views import (exclusion_reason,
                                                            split_views)
from experimental_glycosylation_sites.sequence_qc import (add_sequence_distances,
                                                          sequence_context_failures)
from experimental_glycosylation_sites.provenance import hash_file, _git_state

KEY = ["accession", "position", "population"]

# Sequence context and provenance the views carry so they stand alone.
MANIFEST_COLUMNS = ["uniprot_length", "sequon_triplet", "plus1_class",
                    "occupancy_status", "support_count", "structure_choice",
                    "n_structures_examined", "structure_icode",
                    "glycan_modelled_at_site"]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--features", default="results/datasets/context_features.csv")
parser.add_argument("--manifest", default="results/datasets/context_manifest.csv")
parser.add_argument("--outdir", default="results/datasets")
args = parser.parse_args()

source = Path(args.features)
frame = pd.read_csv(source, low_memory=False)
manifest = pd.read_csv(args.manifest, low_memory=False)

available = [c for c in MANIFEST_COLUMNS if c in manifest.columns]
frame = frame.merge(manifest[KEY + available].drop_duplicates(KEY),
                    on=KEY, how="left", suffixes=("", "_manifest"))
print(f"{len(frame)} rows; joined {len(available)} manifest columns")

# --- row-level sequence checks, before anything is split -------------------
failures = sequence_context_failures(frame)
if len(failures):
    print("\nsequence-context failures by population and reason:")
    print(failures.groupby(["population", "reason"]).size().to_string())
    raise SystemExit(
        f"SEQUENCE QC FAILED: {len(failures)} row(s) lack usable sequence context. "
        "This detects partial cache loss; resolve it rather than dropping rows.")
print("sequence context: every row has a sequence, a complete N-X-S/T triplet "
      "(X != P) and valid coordinates")

frame = add_sequence_distances(frame)
frame["exclusion_reason"] = exclusion_reason(frame)

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)
views = split_views(frame)
written = {}
print()
for name, view in views.items():
    path = outdir / f"context_{name}.csv"
    view.to_csv(path, index=False)
    written[str(path.resolve())] = {"rows": int(len(view)),
                                    "columns": int(len(view.columns)),
                                    "sha256": hash_file(path)}
    counts = view.population.value_counts().to_dict() if "population" in view else {}
    print(f"{name:18} {len(view):6d}  {counts}")

print("\nexclusion reasons from triplet_core:")
excluded = frame[frame.exclusion_reason.ne("")]
if len(excluded):
    print(excluded.groupby(["exclusion_reason", "population"]).size().to_string())

provenance = {
    "inputs": {str(source.resolve()): {"rows": int(len(frame)),
                                       "sha256": hash_file(source)},
               str(Path(args.manifest).resolve()): {"sha256": hash_file(Path(args.manifest))}},
    "outputs": written,
    "git": _git_state(),
}
record = outdir / "context_views_provenance.json"
record.write_text(json.dumps(provenance, indent=2, default=str))
print(f"\nprovenance -> {record.resolve()}")
