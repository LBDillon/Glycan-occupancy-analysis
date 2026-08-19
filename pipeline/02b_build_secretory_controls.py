"""Build the secretory-eukaryotic-unannotated control set.

This set makes the opposite trade from the other two. It matches the occupied
sites on BOTH taxonomy and compartment — eukaryotic, secreted or membrane — so
neither of the confounds the other sets carry applies. It pays for that with a
weaker negative label: these sequons are not annotated as glycosylated, which is
not the same as being annotated unglycosylated.

About half of all eukaryotic secretory proteins with a solved structure carry a
glycoprotein keyword, so the unannotated half certainly contains real
glycosylation sites nobody has recorded. Those false negatives make the two
groups more similar than they truly are, which pushes any measured difference
towards zero. That is the conservative direction given the expected answer is
"no difference", but it means this set can support "no difference detected" far
better than it could support a positive finding.

Stage 1 here: pick one deposited structure per protein and fetch it. Features,
scoreability, matching and scoring follow in the usual runners.
"""
import json, sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "src")
from experimental_glycosylation_sites.config import load_config
from experimental_glycosylation_sites.controls import control_structure_targets
from experimental_glycosylation_sites.structures import fetch_structures

SITES = Path("results/datasets/secretory_unannotated_sites_raw.csv")
config = load_config(Path("config/default.toml"))

sites = pd.read_csv(SITES, low_memory=False)
print(f"{len(sites)} sequons in {sites.accession.nunique()} proteins")

wanted = control_structure_targets(sites, {"secretory_eukaryotic_unannotated": None}, seed=0)
print(f"{len(wanted)} proteins with at least one cross-referenced entry; one structure each")

stats = fetch_structures(
    wanted,
    Path(config.paths["cache_dir"]) / "pdb",
    delay=float(config.api.get("delay_seconds", 0.34)),
    timeout=int(config.api.get("timeout_seconds", 60)),
    per_accession_cap=1,
)
print(json.dumps(stats, indent=2))
