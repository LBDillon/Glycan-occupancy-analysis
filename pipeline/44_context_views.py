"""Split the context table into the three views an analysis may use.

One table cannot serve every question, because the reasons a site is imperfect
are not interchangeable -- see context_views for what each view means and, for
`asn_centred`, what it may not be used for.

Usage:
    44_context_views.py [--features results/datasets/context_features.csv]
                        [--outdir results/datasets]
"""
import argparse, json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.context_views import split_views
from experimental_glycosylation_sites.provenance import hash_file, _git_state

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--features", default="results/datasets/context_features.csv")
parser.add_argument("--outdir", default="results/datasets")
args = parser.parse_args()

source = Path(args.features)
frame = pd.read_csv(source, low_memory=False)
views = split_views(frame)

outdir = Path(args.outdir)
outdir.mkdir(parents=True, exist_ok=True)
written = {}
for name, view in views.items():
    path = outdir / f"context_{name}.csv"
    view.to_csv(path, index=False)
    written[str(path)] = {"rows": int(len(view)), "sha256": hash_file(path)}
    counts = view.population.value_counts().to_dict() if "population" in view else {}
    print(f"{name:18} {len(view):6d}  {counts}")

provenance = {
    "source": {str(source): {"rows": int(len(frame)), "sha256": hash_file(source)}},
    "outputs": written,
    "git": _git_state(),
}
(outdir / "context_views_provenance.json").write_text(json.dumps(provenance, indent=2))
print(f"\nprovenance -> {outdir / 'context_views_provenance.json'}")
