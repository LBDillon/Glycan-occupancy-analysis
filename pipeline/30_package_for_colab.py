"""Build the smallest bundle that lets stages 05/07/08 run somewhere else.

The structure cache is ~13 GB across two directories, and almost none of it is
needed: the model-dependent stages only ever open chains that appear in the
dataset manifest or in one of the matched sets. Those are ~950 chains in ~900
entries, about 1.2 GB uncompressed and roughly a fifth of that gzipped, which
fits comfortably in Drive and uploads once.

Two populations go in, and leaving either out breaks something quietly:

  * every site in `candidate_manifest_dataset` — `09_analyse_scores.py` derives
    its reference SD from ALL scoreable dataset sites, so a bundle holding only
    matched sites would silently rescale every effect size;
  * both members of every pair in each matched set — a paired analysis loses the
    whole pair when one side is missing.

Usage:  30_package_for_colab.py [--out DIR] [--comparison NAME ...] [--tar]
"""
import argparse, gzip, hashlib, json, shutil, sys, time
import pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.runner_support import structure_paths

DEFAULT_COMPARISONS = ("optimal", "secretory", "bacterial", "cytosolic")
KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--out", default="results/colab_bundle")
parser.add_argument("--comparison", action="append", default=[])
parser.add_argument("--tar", action="store_true", help="also write <out>.tar")
parser.add_argument("--no-compress", action="store_true",
                    help="copy structures verbatim instead of gzipping")
args = parser.parse_args()
comparisons = tuple(args.comparison) or DEFAULT_COMPARISONS

out = Path(args.out)
(out / "structures").mkdir(parents=True, exist_ok=True)
(out / "manifests").mkdir(parents=True, exist_ok=True)
(out / "matching").mkdir(parents=True, exist_ok=True)

manifests = {}
for which in ("dataset", "controls", "secretory"):
    path = Path(f"results/manifests/candidate_manifest_{which}.csv")
    if path.exists():
        manifests[which] = pd.read_csv(path, low_memory=False)

# The retention stage runs on its own manifests, not the candidate pools, so a
# bundle built from the candidate pools alone silently lacks ~800 chains.
EXTRA_MANIFESTS = ("scoring_manifest", "manifest_matched_controls",
                   "manifest_matched_secretory")
extra = {}
for stem in EXTRA_MANIFESTS:
    path = Path(f"results/manifests/{stem}.csv")
    if path.exists():
        extra[stem] = pd.read_csv(path, low_memory=False)

# --- decide which sites the bundle has to cover ---------------------------
wanted_sites = set()
dataset = manifests.get("dataset")
if dataset is not None:
    # the reference-SD population, in full
    for a, p in zip(dataset.accession, dataset.position):
        wanted_sites.add((str(a), int(p)))

# Every site in every manifest we ship, so no stage can open a structure the
# bundle lacks. Cheap insurance: the union adds ~800 chains over the matched
# sets, and a missing structure is a silently skipped chain rather than an error.
for frame in extra.values():
    for a, p in zip(frame.accession, frame.position):
        wanted_sites.add((str(a), int(p)))

covered = []
for name in comparisons:
    pairs_path = Path(f"results/matching/matched_pairs_{name}.csv")
    if not pairs_path.exists():
        print(f"  skipping {name}: {pairs_path} absent")
        continue
    pairs = pd.read_csv(pairs_path, low_memory=False)
    for prefix in ("case_", "control_"):
        for a, p in zip(pairs[f"{prefix}accession"], pairs[f"{prefix}position"]):
            wanted_sites.add((str(a), int(p)))
    covered.append(name)

all_manifest = pd.concat(list(manifests.values()) + list(extra.values()),
                         ignore_index=True).drop_duplicates(KEY)
needed = all_manifest[[(str(a), int(p)) in wanted_sites
                       for a, p in zip(all_manifest.accession, all_manifest.position)]]
pdb_ids = sorted({str(x).upper() for x in needed.structure_pdb_id})

print(f"comparisons covered : {', '.join(covered)}")
print(f"sites               : {len(needed)}")
print(f"chains              : {needed.groupby(['structure_pdb_id','structure_chain_id']).ngroups}")
print(f"structures to bundle: {len(pdb_ids)}")

# --- copy the structures --------------------------------------------------
paths = structure_paths()
copied = missing = 0
raw_bytes = out_bytes = 0
t0 = time.time()
for i, pdb_id in enumerate(pdb_ids, 1):
    source = paths.get(pdb_id)
    if source is None:
        missing += 1
        continue
    raw_bytes += source.stat().st_size
    if args.no_compress:
        target = out / "structures" / source.name
        shutil.copyfile(source, target)
    else:
        # Gzip in place. The loaders read .pdb/.cif by suffix, so the notebook
        # gunzips on arrival rather than teaching every parser about gzip.
        target = out / "structures" / (source.name + ".gz")
        with open(source, "rb") as fh, gzip.open(target, "wb", compresslevel=6) as gz:
            shutil.copyfileobj(fh, gz)
    out_bytes += target.stat().st_size
    copied += 1
    if i % 200 == 0:
        print(f"  {i}/{len(pdb_ids)} ({time.time()-t0:.0f}s)", flush=True)

# --- copy the tables ------------------------------------------------------
tables = []
for which in manifests:
    for stem in (f"candidate_manifest_{which}", f"scoreability_{which}"):
        src = Path(f"results/manifests/{stem}.csv")
        if src.exists():
            shutil.copyfile(src, out / "manifests" / src.name)
            tables.append(src.name)
for stem in EXTRA_MANIFESTS:
    src = Path(f"results/manifests/{stem}.csv")
    if src.exists():
        shutil.copyfile(src, out / "manifests" / src.name)
        tables.append(src.name)
for stem in ("scoreability", "scoreability_secretory"):
    src = Path(f"results/manifests/{stem}.csv")
    if src.exists():
        shutil.copyfile(src, out / "manifests" / src.name)
        tables.append(src.name)
for name in covered:
    src = Path(f"results/matching/matched_pairs_{name}.csv")
    shutil.copyfile(src, out / "matching" / src.name)
    tables.append(src.name)

(out / "BUNDLE.json").write_text(json.dumps({
    "built": time.strftime("%Y-%m-%dT%H:%M:%S"),
    "comparisons": covered,
    "n_sites": int(len(needed)),
    "n_chains": int(needed.groupby(["structure_pdb_id", "structure_chain_id"]).ngroups),
    "n_structures": copied,
    "n_structures_missing": missing,
    "compressed": not args.no_compress,
    "bytes_raw": raw_bytes,
    "bytes_bundled": out_bytes,
    "tables": sorted(tables),
    "note": "Structures cover every dataset-manifest site (the reference-SD "
            "population) plus both members of every matched pair.",
}, indent=2))

print(f"\nbundled {copied} structures ({missing} missing)")
print(f"  {raw_bytes/1e9:.2f} GB -> {out_bytes/1e9:.2f} GB")
print(f"  wrote {out}/")
if args.tar:
    archive = shutil.make_archive(str(out), "tar", root_dir=out)
    print(f"  wrote {archive} ({Path(archive).stat().st_size/1e9:.2f} GB)")
