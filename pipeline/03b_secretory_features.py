"""Structural features for the secretory-eukaryotic-unannotated control set.

Same machinery as archive/runners/control_features.py, pointed at the new set. Kept as
its own runner because the set has a different provenance and a different, and
weaker, claim behind its negative label — mixing it into the same output file
would make it easy to forget which is which.
"""
import sys, time
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.features import build_features
from experimental_glycosylation_sites.structures import build_site_evidence

sites = pd.read_csv("results/secretory_unannotated_sites_raw.csv", low_memory=False)
proteins = pd.read_csv("data/cache/secretory_unannotated_proteins.csv.gz", low_memory=False)
sequences = dict(zip(proteins.Entry, proteins.Sequence.fillna("")))

cache = Path("data/cache/pdb")
paths = {}
for p in list(cache.glob("*.pdb")) + list(cache.glob("*.cif")):
    paths.setdefault(p.stem.upper(), p)

manifest, keep = {}, []
for accession, group in sites.groupby("accession"):
    ids = [x.strip() for x in str(group.pdb_ids.iloc[0] or "").split(";") if x.strip()]
    hit = next((i for i in ids if i.upper() in paths), None)
    if hit and sequences.get(accession):
        manifest[accession] = {"accession": accession, "pdb_id": hit.upper(),
                               "output_path": str(paths[hit.upper()]), "all_pdb_ids": ""}
        keep.append(accession)

scoped = sites[sites.accession.isin(keep)].copy()
print(f"sequons with a cached structure: {len(scoped)} across {len(keep)} proteins", flush=True)

chunks, total, t0 = [], len(scoped), time.time()
for start in range(0, total, 500):
    part = scoped.iloc[start:start + 500]
    chunks.append(build_site_evidence(part[["accession", "position"]], sequences, manifest))
    done = min(start + 500, total)
    rate = done / max(time.time() - t0, 1e-9)
    print(f"  mapped {done}/{total} ({rate:.0f}/s, eta {(total-done)/max(rate,1e-9)/60:.0f} min)",
          flush=True)
evidence = pd.concat(chunks, ignore_index=True)

# key-based merge, never positional: build_features re-sorts its output
merged = scoped.merge(evidence, on=["accession", "position"], how="left", validate="one_to_one")
merged["occupancy_status"] = "control_" + merged.control_set

# Chunked purely for visibility. build_features is one uninstrumented call, so
# a single invocation over 1,500 structures runs for half an hour with no way to
# tell progress from a stall. Solvent-accessibility is computed per structure, so
# chunking by structure changes nothing about the result.
chunks, t0 = [], time.time()
groups = list(merged.groupby("structure_pdb_id", dropna=False))
for index, (_, part) in enumerate(groups, 1):
    chunks.append(build_features(part, paths, carry_columns=("control_set",)))
    if index % 100 == 0 or index == len(groups):
        elapsed = time.time() - t0
        print(f"  featured {index}/{len(groups)} structures "
              f"({elapsed/index:.2f}s each, eta {(len(groups)-index)*elapsed/index/60:.0f} min)",
              flush=True)
feats = pd.concat(chunks, ignore_index=True)
feats.to_csv("results/secretory_unannotated_features.csv", index=False)

check = feats.merge(sites[["accession", "position", "control_set"]]
                    .rename(columns={"control_set": "source_set"}),
                    on=["accession", "position"], how="left", validate="one_to_one")
bad = int((check.control_set != check.source_set).sum())
mismatch = int((check.occupancy_status != "control_" + check.source_set).sum())
assert bad == 0 and mismatch == 0, "provenance handoff broken"
print(f"\nprovenance checks passed")

ok = feats[feats.features_available]
print(f"features computed: {len(ok)} of {len(feats)}")
