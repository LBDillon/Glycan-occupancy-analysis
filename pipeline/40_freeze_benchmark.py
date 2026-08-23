"""Snapshot the benchmark state so later work cannot move it underneath.

`results/` is gitignored and rebuilds from `pipeline/` — which is right for
tables that are cheap to regenerate, and wrong for this week's: the corrected
ProteinMPNN scores, both models' retention and the significance tables cost
roughly twenty hours across ARC and a laptop, and several were produced by code
that has since changed.

So this copies the small, decisive artefacts into a dated, VERSIONED directory
and records a SHA-256 for everything else. The claims become reproducible from
git; the large inputs stay out of it but any drift in them is detectable.

What goes in, and why only this:

  * `analysis/` — contrasts, per-comparison JSONs, significance tables. These
    *are* the claims. Small, and nothing else can regenerate them once an input
    moves.
  * `matching/` — the frozen pairs every comparison rests on. Model-independent,
    so freezing them is what makes two models comparable at all.
  * a manifest hashing every score and design table, plus the code state.

Scores and designs are hashed rather than copied: ~10 MB, mechanically
reproducible from the manifests given the same code, and the hash is enough to
prove a later run used the same numbers.

Usage:  40_freeze_benchmark.py [--label 2026-08-23] [--check]
"""
import argparse, hashlib, json, shutil, subprocess, sys
from pathlib import Path

import pandas as pd

ROOT = Path("results")
FREEZE_ROOT = Path("benchmark_frozen")

COPY_DIRS = ("analysis", "matching")
HASH_GLOBS = ("scores/*.csv", "designs/*.csv", "manifests/*.csv")


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def git(*args) -> str:
    try:
        return subprocess.run(["git", *args], capture_output=True, text=True,
                              check=True).stdout.strip()
    except Exception:
        return "unknown"


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--label", default="2026-08-23")
parser.add_argument("--check", action="store_true",
                    help="verify an existing freeze instead of writing one")
args = parser.parse_args()

target = FREEZE_ROOT / args.label
manifest_path = target / "MANIFEST.json"

if args.check:
    if not manifest_path.exists():
        raise SystemExit(f"no freeze at {target}")
    manifest = json.loads(manifest_path.read_text())
    moved, missing = [], []
    for rel, expected in manifest["hashed"].items():
        path = ROOT / rel
        if not path.exists():
            missing.append(rel)
        elif sha256(path) != expected:
            moved.append(rel)
    print(f"freeze {args.label}: {len(manifest['hashed'])} hashed files")
    print(f"  missing : {len(missing)}")
    print(f"  CHANGED : {len(moved)}")
    for rel in moved[:10]:
        print(f"      {rel}")
    if missing[:5]:
        print("  missing examples:", missing[:5])
    raise SystemExit(1 if (moved or missing) else 0)

target.mkdir(parents=True, exist_ok=True)
copied = 0
for sub in COPY_DIRS:
    src = ROOT / sub
    if not src.is_dir():
        continue
    dest = target / sub
    dest.mkdir(parents=True, exist_ok=True)
    for path in sorted(src.glob("*")):
        if path.is_file():
            shutil.copyfile(path, dest / path.name)
            copied += 1

hashed = {}
for pattern in HASH_GLOBS:
    for path in sorted(ROOT.glob(pattern)):
        hashed[str(path.relative_to(ROOT))] = sha256(path)

# a human-readable index of what the freeze actually asserts
claims = {}
for path in sorted((ROOT / "analysis").glob("significance*.csv")):
    frame = pd.read_csv(path)
    variant = path.stem.replace("significance", "").lstrip("_") or "(legacy)"
    claims[variant] = [
        {"outcome": r.outcome, "comparison": r.comparison,
         "effect": round(float(r.effect), 6),
         "p_bh": round(float(r.p_bh), 6)}
        for r in frame.itertuples(index=False)
    ]

manifest_path.write_text(json.dumps({
    "label": args.label,
    "created": pd.Timestamp.utcnow().isoformat(),
    "git_commit": git("rev-parse", "HEAD"),
    "git_dirty": bool(git("status", "--porcelain")),
    "copied_files": copied,
    "hashed": hashed,
    "claims": claims,
    "note": "analysis/ and matching/ are copied verbatim; scores, designs and "
            "manifests are hashed only. Re-run with --check to detect drift.",
}, indent=2))

print(f"froze {copied} files into {target}")
print(f"hashed {len(hashed)} score/design/manifest tables")
print(f"recorded {sum(len(v) for v in claims.values())} significance rows "
      f"across {len(claims)} variants")
print(f"git {git('rev-parse', '--short', 'HEAD')}"
      f"{' (dirty)' if git('status', '--porcelain') else ''}")
