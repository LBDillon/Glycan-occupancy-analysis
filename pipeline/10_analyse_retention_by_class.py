"""Design retention across all five classes, on one footing.

Merges the original sweep (occupied, internal controls, bacterial, cytosolic)
with the later eukaryotic secretory sweep, maps every site to its class, and
reports mean retention with intervals bootstrapped over PROTEINS — several
sequons on one chain share a single set of 32 designs, so a site-level bootstrap
would be far too narrow.

This is the better-powered half of the project. The conditional score asks what
probability the model holds at a site; this asks what it actually writes. Both
are reported, and they should be read together.
"""
import argparse, json, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites import analysis_paths as paths

KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]
FULL = "std_frac_full_sequon_retained"
N_BOOT, BOOT_SEED = 4000, 11

CLASSES = [
    ("occupied_supported", "occupied (experimental glycan)"),
    ("observed_unmodified", "internal control"),
    ("control_secretory_eukaryotic_unannotated", "eukaryotic secretory control"),
    ("control_bacterial_extracytoplasmic", "bacterial control"),
    ("control_cytosolic_eukaryotic", "cytosolic control"),
]

_parser = argparse.ArgumentParser(description=__doc__)
paths.add_variant_argument(_parser)
VARIANT = _parser.parse_args().variant

SOURCES = paths.retention_sources(VARIANT)
print(f"variant: {VARIANT or '(legacy)'}")

frames = []
for retention_path, manifest_path in SOURCES:
    if not (Path(retention_path).exists() and Path(manifest_path).exists()):
        print(f"  skipping {retention_path} (not built yet)")
        continue
    retention = pd.read_csv(retention_path, low_memory=False)
    manifest = pd.read_csv(manifest_path, low_memory=False).drop_duplicates(KEY)
    for frame in (retention, manifest):
        for key in KEY:
            frame[key] = frame[key].astype(str)
    frames.append(retention.merge(manifest[KEY + ["occupancy_status"]], on=KEY, how="left"))

data = pd.concat(frames, ignore_index=True).drop_duplicates(KEY)

# Sites the model cannot decode are excluded: sample() writes their native
# residue back unchanged, so their retention describes the parser, not the model.
# Which sites those are is model-specific, so the variant's own scoreability
# pass is used when it exists -- the unsuffixed files are ProteinMPNN's, and
# filtering ESM-IF or CARBonAra by them would apply another model's exclusions.
sources = paths.scoreability_sources(VARIANT)
if VARIANT:
    borrowed = [s.name for s in sources if paths.suffix(VARIANT) not in s.stem]
    if borrowed:
        print(f"  NOTE: {borrowed} have no pass for variant {VARIANT!r}, so those "
              "exclusions are ProteinMPNN's rather than this model's")
able = pd.concat([pd.read_csv(p, low_memory=False) for p in sources],
                 ignore_index=True)
for key in KEY:
    able[key] = able[key].astype(str)
able = able.drop_duplicates(KEY)
before = len(data)
data = data.merge(able[KEY + ["scoreable"]], on=KEY, how="left")
data = data[data.scoreable == True].copy()
print(f"retention rows {before} -> {len(data)} after dropping undecodable sites\n")


def protein_bootstrap(frame, n_boot=N_BOOT, seed=BOOT_SEED):
    rng = np.random.default_rng(seed)
    groups = [g[FULL].to_numpy() for _, g in frame.groupby("accession")]
    if len(groups) < 2:
        return float("nan"), float("nan")
    draws = np.empty(n_boot)
    for i in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        draws[i] = np.concatenate([groups[j] for j in pick]).mean()
    return tuple(np.percentile(draws, [2.5, 97.5]))


summary = {}
print(f"{'class':34s} {'sites':>6s} {'proteins':>9s} {'retention':>10s} {'95% CI':>18s} {'never':>7s}")
print("-" * 92)
for key, label in CLASSES:
    group = data[data.occupancy_status == key]
    if group.empty:
        print(f"{label:34s} {'--':>6s}  (not present)")
        continue
    mean = float(group[FULL].mean())
    low, high = protein_bootstrap(group)
    never = float((group[FULL] == 0).mean())
    summary[key] = {"label": label, "n_sites": int(len(group)),
                    "n_proteins": int(group.accession.nunique()),
                    "mean_retention": round(mean, 4),
                    "ci95": [round(low, 4), round(high, 4)],
                    "frac_never_retained": round(never, 4)}
    print(f"{label:34s} {len(group):>6d} {group.accession.nunique():>9d} "
          f"{mean:>10.4f} [{low:>+7.4f},{high:>+7.4f}] {never:>7.1%}")

occ = summary.get("occupied_supported")
if occ:
    print("\ndifference from occupied (occupied minus control):")
    for key, entry in summary.items():
        if key == "occupied_supported":
            continue
        print(f"  {entry['label']:32s} {occ['mean_retention'] - entry['mean_retention']:+.4f}")

paths.retention_by_class_json(VARIANT).write_text(json.dumps({
    "unit": "one site; 32 unconstrained designs per chain at temperature 0.1",
    "bootstrap": "over proteins, not sites",
    "excluded": "sites ProteinMPNN cannot decode",
    "classes": summary,
}, indent=2))
data.to_csv(str(paths.retention_all_classes(VARIANT)), index=False)
print(f"\nwrote {paths.retention_by_class_json(VARIANT)} and "
      f"{paths.retention_all_classes(VARIANT)}")
