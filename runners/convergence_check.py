import sys, json, time
import numpy as np, pandas as pd
from pathlib import Path
sys.path.insert(0, "src")
from experimental_glycosylation_sites.mpnn_scoring import (
    load_model, conditional_probabilities, sequon_score)

KEY = ["accession", "position", "structure_pdb_id", "structure_chain_id"]
m = pd.read_csv("results/scoring_manifest.csv", low_memory=False)
unique = m.drop_duplicates(KEY).reset_index(drop=True)

# BLINDED: sampled with a fixed seed from all unique scoreable sites, without
# reference to occupancy_status or control_set.
sample = unique.sample(n=50, random_state=20260817).reset_index(drop=True)
assert "occupancy_status" in sample.columns  # present but deliberately unused here

paths = {}
for d in ("data/cache/pdb", "../ortholog_sequon_conservation/results/database_current/structures/pdb"):
    for p in list(Path(d).glob("*.pdb")) + list(Path(d).glob("*.cif")):
        paths.setdefault(p.stem.upper(), p)

model = load_model(Path("../../ProteinMPNN"))
rows, t0 = [], time.time()
for i, r in enumerate(sample.itertuples(index=False), 1):
    path = paths.get(str(r.structure_pdb_id).upper())
    if path is None:
        continue
    idx0 = (int(r.n_model_index), int(r.plus1_model_index), int(r.plus2_model_index))
    try:
        probs16, computed = conditional_probabilities(path, r.structure_chain_id, model,
                                            n_decoding_orders=16, seed=0, positions=list(idx0))
    except Exception as exc:
        print(f"  [{i}] {r.accession} {r.structure_pdb_id}: {type(exc).__name__}", flush=True)
        continue
    entry = {"accession": r.accession, "position": r.position}
    for n in (4, 8, 16):
        entry[f"score_{n}"] = sequon_score(probs16[:n], *idx0,
                                           computed=computed)["conditional_sequon_score"]
    rows.append(entry)
    if i % 10 == 0:
        print(f"  {i}/50 ({(time.time()-t0)/i:.1f}s per site)", flush=True)

df = pd.DataFrame(rows)
df.to_csv("results/convergence_check.csv", index=False)

# reference SD from the blinded sample, on the 16-order scores
ref_sd = float(df.score_16.std(ddof=1))
d8 = (df.score_8 - df.score_16).abs() / ref_sd
d4 = (df.score_4 - df.score_16).abs() / ref_sd
out = {
    "n_sites": len(df),
    "sample_seed": 20260817,
    "reference_sd_in_sample": round(ref_sd, 4),
    "8_vs_16": {"median_sd_units": round(float(d8.median()), 4),
                "p95_sd_units": round(float(d8.quantile(0.95)), 4),
                "correlation": round(float(df.score_8.corr(df.score_16)), 5)},
    "4_vs_16": {"median_sd_units": round(float(d4.median()), 4),
                "p95_sd_units": round(float(d4.quantile(0.95)), 4),
                "correlation": round(float(df.score_4.corr(df.score_16)), 5)},
}
passes = (out["8_vs_16"]["median_sd_units"] < 0.02
          and out["8_vs_16"]["p95_sd_units"] < 0.05
          and out["8_vs_16"]["correlation"] > 0.99)
out["decision"] = "adopt 8 decoding orders" if passes else "escalate to 16 decoding orders"
Path("results/convergence_check.json").write_text(json.dumps(out, indent=2))
print()
print(json.dumps(out, indent=2))
