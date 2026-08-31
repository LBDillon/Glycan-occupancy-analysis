# Running the model stages on this laptop

`gca` is a standalone checkout: it holds no structures and no model
checkpoints. Both live in the SugarFix tree beside it, and the built-in
fallback paths only resolve when `gca` sits *inside* that tree. Standalone, the
two environment variables below are not optional — without them every chain
drops as `structure_not_cached` and a run produces an empty table while exiting
zero.

```bash
export GCA_STRUCTURE_DIRS=/Users/lauradillon/SugarFix-merge/analysis/experimental_glycosylation_sites/data/cache/pdb:/Users/lauradillon/SugarFix-merge/analysis/ortholog_sequon_conservation/results/database_current/structures/pdb
export PROTEINMPNN_DIR=/Users/lauradillon/SugarFix-merge/ProteinMPNN
```

That resolves 8,702 structures, which covers all 248 chains of
`candidate_manifest_dataset.csv`. `PROTEINMPNN_DIR` must be the checkout
containing `protein_mpnn_utils.py` and `vanilla_model_weights/v_48_020.pt`.

## Which models actually run here

Checked 2026-08-31 against every conda environment on the machine.

| Model | Environment | State | Run it on Colab |
|---|---|---|---|
| **ProteinMPNN** | `envs/esm3` | **Works.** ~8.5 s/chain on CPU, so ~35 min for 248 chains at 32 designs. | — |
| **ESM3** | `envs/esm3` | Works, but >2 min/chain for 32 designs on CPU — roughly 8+ hours for the same set. `mps` is available and untested. | [`cross_protein_esm3.ipynb`](../scripts/colab/cross_protein_esm3.ipynb) |
| **ESM-IF** | none | **Blocked.** Needs `fair-esm` (for `esm.inverse_folding`) plus `biotite`, `torch-geometric`, `torch-scatter`. No environment has them. | [`cross_protein_esm_if.ipynb`](../scripts/colab/cross_protein_esm_if.ipynb) |
| **CARBonAra** | none | **Blocked.** Needs `blosum`, `gemmi`, `h5py`, `scikit-learn`, `Bio`, *and* a CARBonAra checkout. No checkout exists on this machine. | [`cross_protein_carbonara.ipynb`](../scripts/colab/cross_protein_carbonara.ipynb) |

The blocker for ESM-IF is not a missing package but a name collision, and it is
the one documented in [`models.md`](models.md): `fair-esm` and
EvolutionaryScale's `esm` both install a top-level module called `esm`. The
`esm3` environment holds the second, so ESM-IF cannot join it. It needs its own
environment, as it has on ARC (`venv-if/`).

Nothing here is a code problem — the adapters are fine and the registry
lazy-imports, so an absent model is unavailable rather than fatal. It is purely
that two of the four environments have never been built on this machine.

## Sanity check before a long run

```bash
python -c "
import sys; sys.path.insert(0,'src')
from experimental_glycosylation_sites.runner_support import structure_paths, build_adapter
print('structures:', len(structure_paths(())))
print(build_adapter('proteinmpnn','cpu',max_batch=None).describe())
"
```

Expect 8,702 structures and `{'model': 'v_48_020', ...}`. A `FileNotFoundError`
naming `protein_mpnn_utils.py` means `PROTEINMPNN_DIR` is unset; a structure
count of 0 means `GCA_STRUCTURE_DIRS` is.

## The Colab notebooks

Three notebooks in [`../scripts/colab/`](../scripts/colab/) cover the models this
machine cannot run, for the cross-protein retention experiment specifically. Each
does the same five things — check out the branch, unpack the release bundle, run
a preflight, run `08_design.py --save-sequences` over
`candidate_manifest_dataset.csv`, then `57_cross_protein_sequon_retention.py` on
the sequences it saved.

They are three rather than one for the reason in [`models.md`](models.md):
`fair-esm` and EvolutionaryScale's `esm` both install a top-level module named
`esm`, so ESM-IF and ESM3 cannot share a runtime. CARBonAra is separate because
its dependency set is unrelated to both, following the split in
`scripts/arc/glyco_setup.sh`.

**They clone from GitHub, so anything not pushed does not exist for them.** Each
asserts that `--save-sequences` and stage 57 are present in the checkout and
stops with an explanation if they are not — which is what a stale
`/content/module` or an unpushed branch looks like.

The ESM-IF notebook carries an optional ProteinMPNN validation arm. It reruns the
model this laptop already did and asserts the answer comes back at 13.49% pattern
retention, which checks the whole path — bundle, structures, manifest, both
stages — against a known number before the other model's results are trusted.

→ [`running_on_arc.md`](running_on_arc.md) for the cluster, where the job arrays
and all four environments live.
