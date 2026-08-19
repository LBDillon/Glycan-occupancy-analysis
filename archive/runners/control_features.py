import sys, time, pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.structures import build_site_evidence
from experimental_glycosylation_sites.features import build_features

sites = pd.read_csv("results/negative_control_sites.csv", low_memory=False)
proteins = pd.read_csv("data/cache/negative_control_proteins.csv.gz", low_memory=False)
sequences = dict(zip(proteins.Entry, proteins.Sequence.fillna("")))

cache = Path("data/cache/pdb")
paths = {}
for p in list(cache.glob("*.pdb")) + list(cache.glob("*.cif")):
    paths.setdefault(p.stem.upper(), p)

manifest, keep = {}, []
for acc, group in sites.groupby("accession"):
    ids = [x.strip() for x in str(group.pdb_ids.iloc[0] or "").split(";") if x.strip()]
    hit = next((i for i in ids if i.upper() in paths), None)
    if hit and sequences.get(acc):
        manifest[acc] = {"accession": acc, "pdb_id": hit.upper(),
                         "output_path": str(paths[hit.upper()]), "all_pdb_ids": ""}
        keep.append(acc)

scoped = sites[sites.accession.isin(keep)].copy()
print(f"control sequons with a cached structure: {len(scoped)} across {len(keep)} proteins", flush=True)

chunks, n, t0 = [], len(scoped), time.time()
for start in range(0, n, 500):
    part = scoped.iloc[start:start + 500]
    chunks.append(build_site_evidence(part[["accession", "position"]], sequences, manifest))
    done = min(start + 500, n)
    rate = done / max(time.time() - t0, 1e-9)
    print(f"  mapped {done}/{n} ({rate:.0f}/s, eta {(n-done)/max(rate,1e-9)/60:.0f} min)", flush=True)
ev = pd.concat(chunks, ignore_index=True)

# key-based merge, never positional: build_features re-sorts its output
merged = scoped.merge(ev, on=["accession", "position"], how="left", validate="one_to_one")
merged["occupancy_status"] = "control_" + merged.control_set

feats = build_features(merged, paths, carry_columns=("control_set",))
feats.to_csv("results/negative_control_features.csv", index=False)

# provenance must agree with the source inventory, by key
check = feats.merge(sites[["accession", "position", "control_set"]]
                    .rename(columns={"control_set": "source_set"}),
                    on=["accession", "position"], how="left", validate="one_to_one")
bad = int((check.control_set != check.source_set).sum())
mismatch = int((check.occupancy_status != "control_" + check.source_set).sum())
print(f"\nprovenance disagreeing with source: {bad}")
print(f"occupancy_status disagreeing with source: {mismatch}")
assert bad == 0 and mismatch == 0, "provenance handoff still broken"

ok = feats[feats.features_available]
print(f"features computed: {len(ok)} of {len(feats)}")
print(ok.groupby("control_set").size().to_string())
