# experimental_glycosylation_sites

> **⚠ Correction, 2026-08-20 — ProteinMPNN results are being regenerated.**
> `mpnn_scoring.ALPHABET` held a three-letter lookup table from inside
> `parse_PDB_biounits` instead of ProteinMPNN's token alphabet. Consequence:
> `p_asn_at_n` was reading **P(aspartate)**, and designed sequences were decoded
> with the wrong letters. P(Ser) and P(Thr) were correct by coincidence, so the
> hydroxyl half of every sequon score is sound and the asparagine half is not.
> **Do not quote any ProteinMPNN score, retention figure or SD-standardised
> effect produced before this date** — including the numbers currently in
> `docs/primary_result.md`, `docs/significance.md` and `docs/figures.md`.
> Scoreability, matching and the manifests are unaffected. Full account:
> [`docs/correction_2026-08-20_alphabet.md`](docs/correction_2026-08-20_alphabet.md).

Site-level N-linked glycosylation occupancy evidence derived from the ortholog
sequon-conservation database.

- Evidence source guide: [`docs/evidence_sources.md`](docs/evidence_sources.md)
- Downstream analysis options: [`docs/analysis_options.md`](docs/analysis_options.md)

---

# Running the whole thing

Everything runs from this directory. Stages are numbered in execution order;
letters mark optional branches. All of it is reproducible from committed code —
`results/` is gitignored and rebuilds from these scripts.

```bash
cd analysis/experimental_glycosylation_sites
export KMP_DUPLICATE_LIB_OK=TRUE OMP_NUM_THREADS=1   # macOS duplicate libomp
```

### 1. Build the site tables

```bash
python -m experimental_glycosylation_sites run          # the 4,307-site universe
python -m experimental_glycosylation_sites fetch-controls
python pipeline/02b_build_secretory_controls.py         # eukaryotic secretory set
```

### 2. Structural features

```bash
python pipeline/03b_secretory_features.py               # ~35 min
```

### 3. Manifest, then scoreability — in that order

Scoreability is settled **before** matching, from coordinates alone. Doing it
afterwards lets unscoreable sites into matched sets and removes them later,
unbalancing sets that matching had just balanced.

```bash
python pipeline/04_build_candidate_manifest.py dataset
python pipeline/04_build_candidate_manifest.py controls
python pipeline/04_build_candidate_manifest.py secretory
python pipeline/05_scoreability.py results/manifests/candidate_manifest_dataset.csv \
                                   results/manifests/scoreability_dataset.csv
```

### 4. Match

```bash
python pipeline/06_match_primary.py        # internal controls: optimal + sensitivities
python pipeline/06b_match_diagnostics.py   # bacterial, cytosolic, secretory
```

### 5. Score and design — the two model-dependent stages

Both take `--model` (default `proteinmpnn`) and `--device` (default `cpu`), so a
second model needs no second copy of either file.

```bash
python pipeline/07_score.py  <manifest.csv> results/scores/scores_<set>.csv     # ~5s/chain
python pipeline/08_design.py <manifest.csv> results/retention_<set>.csv  # ~25s/chain

python pipeline/07_score.py  <manifest.csv> results/scores/scores_<set>_esm_if.csv --model esm_if
python pipeline/08_design.py <manifest.csv> results/designs/retention_<set>_esm_if.csv --model esm_if
```

ESM-IF also needs its own scoreability pass — its parser is biotite, not
Biopython, so which sites it can address is a different question:

```bash
python pipeline/05_scoreability.py <manifest.csv> <out.csv> --model esm_if
```

**On a GPU.** Bundle the ~900 structures these stages actually open (~1.2 GB, a
fifth of that gzipped, versus 13 GB for the whole cache) and run the notebook:

```bash
python pipeline/30_package_for_colab.py --out results/colab_bundle --tar
# upload results/colab_bundle.tar to Drive, then open
# notebooks/esm_if_and_mpnn_gpu.ipynb in Colab
```

### 6. Analyse

Every analysis stage takes `--variant`, naming which run's numbers to read and
write. Omitting it reproduces the original ProteinMPNN filenames and numbers
exactly.

```bash
python pipeline/09_analyse_scores.py optimal        # PRIMARY
python pipeline/09_analyse_scores.py secretory      # parallel, best powered
python pipeline/10_analyse_retention_by_class.py
python pipeline/10b_analyse_retention_paired.py
python pipeline/11_significance.py                  # all 8 tests, corrected
```

A second model, or the corrected ProteinMPNN run, is the same commands with a tag:

```bash
python pipeline/09_analyse_scores.py optimal --variant alphabet_corrected
python pipeline/09_analyse_scores.py optimal --variant esm_if
python pipeline/11_significance.py               --variant esm_if
```

The tag suffixes every input and output — `scores_dataset_esm_if.csv` in,
`analysis_optimal_esm_if.json` out — so two models' results never overwrite each
other. **A variant whose score file is missing stops with an error rather than
falling back to the default**, because reading the wrong score file is otherwise
silent. Model provenance in the output JSON is read from the score file's own
`model` / `conditioning` / `n_orders` columns, never restated by the analysis.

Matching carries no variant: it is built from RSA, neighbour counts and
hydrophobic fraction, never from model output, so every model shares one set of
pairs.

### 7. Figures

```bash
python pipeline/20_figures_summary.py
python pipeline/21_figures_all_classes.py
python pipeline/22_figures_control_provenance.py
python pipeline/23_figures_primary.py optimal
```

Supporting checks: `12_matching_sensitivity.py` (200-seed sweep),
`13_name_audit.py` (are the "unannotated" controls really unannotated),
`14_convergence_check.py` (how many decoding orders are enough).

## Where things live

| Path | Contents |
|---|---|
| `data/raw/`, `data/cache/` | inputs and API/structure caches; never written by analysis |
| `results/manifests/candidate_manifest_*.csv` | one row per site, with its three model indices |
| `results/manifests/scoreability_*.csv` | which sites the model can evaluate, decided pre-matching |
| `results/matching/matched_pairs_*.csv` | the pairs each comparison rests on |
| `results/scores/scores_*.csv` | **sequon scores** — one row per site per model |
| `results/retention_*.csv`, `mpnn_retention*.csv` | **designs** — retention per site |
| `results/analysis/analysis_*.json`, `contrasts_*.csv` | contrasts, intervals, verdicts |
| `results/analysis/significance.csv` | all eight tests with corrections |
| `results/figures/` | every figure |
| `archive/` | superseded runners and outputs, kept for provenance |

Datasets used: a dated UniProt snapshot, GlyGen and GlyConnect API caches, and
deposited PDB/mmCIF structures. All read-only and referenced by path — see
**Inputs are read-only** below.

---

# Adding another model

The analysis downstream of scoring is model-agnostic: matching, contrasts, the
cluster bootstrap, significance testing and the figures all work on tables keyed
by `(accession, position)`. Adding ESM-IF, ESM3 or anything else means writing
**one adapter** and touching nothing else.

Three models are registered:

| Name | Conditions on | Implements | Notes |
|---|---|---|---|
| `proteinmpnn` | backbone + all other native residues, 8 decoding orders | scorer + designer | [correction](docs/correction_2026-08-20_alphabet.md) |
| `esm_if` | backbone + native prefix (autoregressive) | scorer + designer | [doc](docs/second_model_esm_if.md) |
| `esmc` | **sequence only** (masked LM) | scorer | [doc](docs/third_model_esmc.md) |

**`esmc` cannot be installed alongside `esm_if`.** `fair-esm` and
EvolutionaryScale's `esm` both claim the top-level import name `esm`. Use a
separate environment for each; the registry lazy-imports, so the absent model is
unavailable rather than breaking the package.

`src/experimental_glycosylation_sites/adapters/` defines two protocols. A model
may implement either or both:

| Protocol | Question it answers | Feeds |
|---|---|---|
| `SequonScorer` | what probability does the model hold at the three sequon residues? | `07_score.py` |
| `SequenceDesigner` | what does the model write when redesigning the chain? | `08_design.py` |

`SequonScorer` splits its work in two so that neither model wastes effort:
`prepare_chain(path, chain, positions)` does the once-per-chain computation and
returns an opaque context, `score_from(context, indices, expected_triplet)` reads
one sequon out of it, and `score_site` is the two composed. ProteinMPNN pays per
position, so it wants every position on a chain at once; ESM-IF decodes the whole
chain in one pass, so its second sequon should cost nothing.

**Steps**

1. Copy `adapters/proteinmpnn.py` to `adapters/<model>.py`.
2. Implement `decodable_positions`, `prepare_chain`/`score_from`/`score_site`,
   and/or `design`, plus `describe()` for the provenance columns.
3. Register it in `adapters/__init__.py`.
4. Run stages 5–7 with `--model <name>`; the comparison figures gain a series.

**Two invariants, both learned the hard way**

- **Never score a residue the model did not actually evaluate.** ProteinMPNN
  returns a zero row for residues with incomplete backbones, which exponentiates
  to twenty-one ones and scores +13.8. That defect inverted the sign of the first
  result. Report such sites as unscoreable; the scorer raises rather than
  returning a value.
- **`decodable_positions` must not need a model pass.** Scoreability has to be
  answerable before matching, or matched sets lose members afterwards.
- **Never trust the manifest's index across a parser boundary.** `model_index`
  is an ordinal into the chain as Biopython reads it. A model whose own parser
  disagrees must translate, and must verify that the mapped residues reproduce
  the manifest's triplet — ESM-IF disagreed on ~5% of sites before its mapping
  was built. Related: the alphabet a model returns probabilities in is an
  assumption until you round-trip it against that model's own output. Assuming
  it is what produced the 2026-08-20 correction.

Verify conformance:

```python
from experimental_glycosylation_sites import adapters
from experimental_glycosylation_sites.adapters.base import SequonScorer, SequenceDesigner
a = adapters.load("proteinmpnn")
isinstance(a, SequonScorer), isinstance(a, SequenceDesigner)   # (True, True)
```

---

# Where to start reading

| Document | For |
|---|---|
| [`docs/primary_result.md`](docs/primary_result.md) | the result and its limits |
| [`docs/figures.md`](docs/figures.md) | all nine figures explained in prose |
| [`docs/concepts.md`](docs/concepts.md) | what the terms mean, in plain language |
| [`docs/significance.md`](docs/significance.md) | the eight tests and why none survives correction |
| [`docs/negative_controls.md`](docs/negative_controls.md) | the four control sets and what evidence stands behind each |
| [`docs/correction_2026-08-18.md`](docs/correction_2026-08-18.md) | what was corrected and why |
| [`docs/correction_2026-08-20_alphabet.md`](docs/correction_2026-08-20_alphabet.md) | **the alphabet defect — read before quoting any ProteinMPNN number** |
| [`docs/second_model_esm_if.md`](docs/second_model_esm_if.md) | ESM-IF: what its conditional is, the index mapping, how to run it |
| [`docs/third_model_esmc.md`](docs/third_model_esmc.md) | ESMC: sequence-only baseline, the two masking schemes, the environment split |
| [`docs/adding_models_explainer.md`](docs/adding_models_explainer.md) | **how the benchmark went from one model to three, and what went wrong on the way** |
| `config/scoring_frozen.toml` | the frozen configuration and both amendments |

---

## Purpose

The ortholog sequon-conservation database is built around **pairs**. Each row is
one orthologous comparison in which one protein carries an N-X-S/T sequon and
the other has lost it. That unit is right for evolutionary questions and wrong
for claims about glycosylation, for two reasons.

**A sequon is a motif, not a modification.** N-X-S/T is necessary but not
sufficient for N-linked glycosylation. Occupancy depends on local structure,
topology, and whether the protein ever meets an oligosaccharyltransferase.
"This pair lost a sequon" and "this pair lost a glycan" are different claims,
and pair rows only support the first.

**Pair rows repeat sites.** One asparagine on one protein appears once per
ortholog it was compared against. In the current analysis-ready set:

```
13,816 analysis-ready pair rows  →  4,307 unique (accession, position) sites
```

Counting rows would turn one biochemical fact into up to dozens of apparent
independent observations. That collapse is the reason this module exists.

So the unit of analysis here is **one UniProt accession plus one 1-indexed
residue position**. Ortholog pairs are retained as context in a separate table
(`site_pair_associations.csv`) where they cannot inflate site counts.

### What the module asks of each site

Four independent evidence layers, each individually toggleable in
`config/default.toml`:

| Layer | Source | Nature | Default |
|---|---|---|---|
| `uniprot` | Exact CARBOHYD N-linked feature plus ECO codes | Curated, auditable per feature | on |
| `glygen` | GlyGen protein-detail API `site_category` | Independent aggregation, largely mass spectrometry | on |
| `structure` | ASN-glycan `LINK` records in cached PDB files | Direct physical observation | on |
| `glyconnect` | GlyConnect REST API site strings | Corroboration only, low coverage | off |

Layers are independent. Any one of them can support a site on its own, and a
site that fails UniProt's evidence policy can still be supported by a structural
glycan linkage.

### Absence of annotation is not evidence of absence

A site with no supporting layer is `unknown` — **never** a biological negative.
Well-studied proteins accumulate annotations; obscure ones do not, so treating
unannotated sites as negatives would measure curation effort and report it as
glycobiology.

`occupancy_status` therefore takes three values, and the third is deliberately
narrow. From annotation alone a site is `occupied_supported` or `unknown`; no
annotation source can establish that a site was *examined and found bare*.

`observed_unmodified` is populated from structures instead, and only where
absence is informative: the entry models a glycan at some other residue, so
sugars demonstrably survived preparation and this depositor demonstrably
modelled them, and the protein was expressed in a host that can glycosylate.
**32 sites across 25 proteins** meet both conditions.

They are best described as *sequons with no modelled glycan under
internal-control conditions*. That is strong evidence of absence by structural
standards while still being a statement about the deposited model rather than
the molecule — not a definitive biochemical negative, but the most informative
internal control available, and the one the occupancy comparison runs on.

The same caution applies to structures. `structure_residue_resolved` means the
residue is modelled with no glycan linkage attached. It is not evidence the site
was unoccupied — glycans are routinely removed before crystallisation, proteins
are frequently expressed in bacterial systems that do not glycosylate, and
glycans that are present are often left unmodelled through disorder.

## Inputs are read-only and referenced by path

Nothing under `analysis/ortholog_sequon_conservation/` or `data/raw/` is ever
written, moved, or copied. Every input is referenced by a configurable path in
`config/default.toml`, resolved **relative to the config file's own directory**:

| Config key | Points at | Consumed by |
|---|---|---|
| `pairs_master` | Canonical pair table | Candidate universe, associations |
| `homology_qc` | Homology QC buckets per pair | Strict / plausible subsets |
| `uniprot_tsv` | Dated UniProt snapshot (gzipped TSV) | CARBOHYD features, sequences |
| `structure_manifest` | Cached PDB download manifest | Structure layer |
| `proteins_master` | Canonical protein table | Declared and existence-checked; not yet read |
| `structure_dir` | Cached PDB directory | Declared and existence-checked; the manifest's `output_path` is used instead |
| `existing_structural_context` | Prior per-site RSA / SS table | Declared and existence-checked; consumed by downstream analysis, not by this pipeline |

`validate` (below) checks every non-output path exists, so a wrong path fails
loudly instead of producing a silently empty join.

Outputs go only to `results/` and `data/cache/`, both inside this module. Their
contents are gitignored; only the `.gitignore` and `README.md` are tracked.

## Setup

Requires Python >= 3.12. Use the project interpreter:

```bash
cd analysis/experimental_glycosylation_sites
/Users/lauradillon/miniforge3/bin/python3.12 -m pip install -e .
```

That installs `pandas>=2.0`, `biopython>=1.80`, `numpy>=1.24` and puts the `src/`
layout package on the path. Run the tests to confirm the install:

```bash
/Users/lauradillon/miniforge3/bin/python3.12 -m pytest -q
```

## Commands

Four commands, all through `python -m experimental_glycosylation_sites`. Each
accepts `--config PATH` (default: `config/default.toml`) and validates every
configured input before doing anything else.

```bash
cd analysis/experimental_glycosylation_sites
PY=/Users/lauradillon/miniforge3/bin/python3.12
```

**`validate`** — check that every configured input path resolves. Does no work
and touches no output. Run it first after moving the module or editing paths.

```bash
$PY -m experimental_glycosylation_sites validate
# All configured inputs resolve.
```

**`fetch-glygen`** — populate the GlyGen cache for every candidate accession.
Resumable: accessions already in the cache are skipped, so an interrupted fetch
can simply be rerun. Rate-limited by `api.delay_seconds`. Prints the cache path.

```bash
$PY -m experimental_glycosylation_sites fetch-glygen
# cache: .../data/cache/glygen_protein_detail.jsonl
```

**`fetch-glyconnect`** — same, for GlyConnect. Two API calls per accession
(search by UniProt accession for a protein id, then fetch that protein's
structures). Only needed if you enable the `glyconnect` layer.

```bash
$PY -m experimental_glycosylation_sites fetch-glyconnect
```

**`run`** — run every enabled layer, write all result tables, print the summary
counts as JSON. Offline by default; pass `--fetch` to refresh the API caches
first.

```bash
$PY -m experimental_glycosylation_sites run
$PY -m experimental_glycosylation_sites run --fetch          # refresh caches first
$PY -m experimental_glycosylation_sites run --config config/uniprot_only.toml
```

## Outputs

Written to `results/` (configurable via `paths.results_dir`). Every candidate
site appears exactly once across `experimental_sites_all.csv` and
`excluded_sites.csv`.

| File | Rows | Meaning |
|---|---|---|
| `candidate_sites.csv` | one per candidate site | The deduplicated `(accession, position)` universe. The denominator for everything else. |
| `site_pair_associations.csv` | one per (site, pair) | Ortholog context: `pair_id`, `cluster_id`, paired accession, source, `homology_qc_bucket`. Kept separate so pair rows can never inflate site counts. |
| `uniprot_exact_n_linked_sites.csv` | one per parsed CARBOHYD feature | Full inventory of what the UniProt snapshot actually contained, including features no candidate site uses, and including range / uncertain / non-N-linked features with their `parse_status`. Provenance, not a headline result — it lets exclusions be audited against the source. |
| `experimental_sites_uniprot_baseline.csv` | one per UniProt-qualifying site | UniProt-only pass. This is the frozen regression baseline the snapshot test guards. |
| `experimental_sites_all.csv` | one per supported site | All enabled layers combined. `experimental_positive == True`. |
| `experimental_sites_strict_plus_plausible.csv` | subset of the above | Sites with at least one association in `strict_ortholog_like` or `plausible_ortholog_like`. |
| `experimental_sites_strict.csv` | subset of the above | Sites with at least one `strict_ortholog_like` association. |
| `excluded_sites.csv` | one per unsupported site | Every candidate no layer supports, each with a machine-readable `exclusion_reason`. These are `unknown`, not negatives. |
| `curator_inferred_sensitivity_sites.csv` | one per `manual_curator_inference` site | Sensitivity set: sites whose strongest UniProt tier is ECO:0000305. Written when `policy.curator_inferred_sensitivity` is true. |
| `glygen_site_evidence.csv` | one per candidate site | GlyGen layer detail: tier, raw categories, evidence databases, PubMed ids, GlyTouCan ids. Written when the `glygen` layer is on. |
| `glyconnect_site_evidence.csv` | one per candidate site | GlyConnect layer detail: supported flag, glycan count, compositions. Written when the `glyconnect` layer is on. |
| `structure_site_evidence.csv` | one per candidate site | Structure layer detail: `structure_tier`, PDB id, chain, resseq, insertion code, the residue actually observed at the mapped position, and `structure_detail`. Read `structure_tier` **and** `structure_detail` together — see the limitations in `docs/analysis_options.md`. |
| `summary.json` | — | All headline counts: candidates, exclusions, per-subset sites and proteins, baseline and enriched, and per-layer support totals. |
| `provenance.json` | — | Everything needed to explain a result set: every input path with SHA-256, size and mtime; the effective config and its hash; git commit and dirty flag; all counts; and per-layer cache sizes. |

### Key site-level columns

| Column | Meaning |
|---|---|
| `uniprot_tier` | Strongest ECO tier at the exact position, or `exact_feature_absent` |
| `uniprot_evidence_codes` | Pipe-delimited ECO codes retained verbatim |
| `glygen_tier`, `structure_tier` | Per-layer classification (see `docs/evidence_sources.md`) |
| `support_sources` | Pipe-delimited supporting layers in fixed order, e.g. `uniprot\|structure` |
| `support_count` | How many independent layers support the site |
| `experimental_positive` | Boolean, derived from the configured policy |
| `occupancy_status` | `occupied_supported` or `unknown` |
| `in_strict`, `in_strict_plus_plausible` | Ortholog subset membership |
| `n_associations` | How many pair rows this one site collapsed from |

### Current build

From the latest full run (`results/datasets/summary.json`):

| Quantity | Value |
|---|---|
| Analysis-ready pair rows | 13,816 |
| Unique candidate sites | 4,307 |
| UniProt baseline, all | 505 sites / 401 proteins |
| UniProt baseline, strict + plausible | 357 sites / 303 proteins |
| UniProt baseline, strict only | 321 sites / 278 proteins |
| Enriched (all layers), all | 922 sites / 703 proteins |
| Enriched, strict + plausible | 460 sites / 381 proteins |
| Enriched, strict only | 396 sites / 333 proteins |
| Excluded | 3,385 (922 + 3,385 = 4,307) |
| Sites with UniProt support | 505 |
| Sites with GlyGen support | 774 |
| Sites with structural glycan support | 172 |
| Sites failing UniProt policy but carrying a structural glycan | 32 |

Those 32 sites are the clearest demonstration that the layers are independent:
each has a modelled glycan covalently linked to the asparagine in a cached PDB
entry, yet its UniProt annotation rests on a sequence model, a curator
inference, or no evidence code at all. Twenty-four of them are still supported
by structure alone once GlyGen is included.

Enrichment adds 417 sites to the 505-site UniProt baseline without moving the
baseline itself: 383 of the additions carry GlyGen support, 32 carry a
structural glycan, and 8 carry both. The GlyGen layer was populated from 1,714
cached accessions — the candidates that carry a GlyGen cross-reference, which
are the only ones worth requesting (`provenance.json` records
`glygen_cached_accessions: 1714`).

## Caching, offline reruns, and determinism

API responses are cached as JSONL under `data/cache/`, one line per accession
recording `{accession, detail, error}`:

- `data/cache/glygen_protein_detail.jsonl`
- `data/cache/glyconnect_protein_detail.jsonl`

Three consequences.

**Reruns are offline.** `run` without `--fetch` never opens a socket. It reads
the caches, the UniProt TSV, the canonical tables, and the cached PDB files.
Once the caches are populated the whole module runs with no network.

**Fetches are resumable.** Each fetch loads the existing cache first and skips
accessions already present, appending and flushing line by line. An interrupted
or rate-limited fetch is fixed by rerunning the same command. Failures are
recorded as a row with a null `detail` and a non-empty `error`, so a failed
accession is visible rather than silently absent — and because only rows with a
non-null `detail` are loaded, a failure is retried on the next fetch.

**Results are deterministic.** Every table is written through `table_io.write_table`,
which sorts on `(accession, position)` and resets the index before writing, so
two runs over the same inputs produce byte-identical CSVs and diff cleanly.
`provenance.json` hashes every input so a changed result can always be traced to
a changed input.

## Reading a snapshot test failure

`tests/test_snapshot.py` compares the UniProt-only counts against the frozen
fixture `tests/snapshots/uniprot_baseline_2026-04-27.json`. It runs against the
canonical inputs when they are present and skips with an explicit reason when
they are not.

On mismatch it prints the count deltas *and* a short SHA-256 fingerprint of every
configured input file:

```
Baseline counts drifted from the 2026-04-27 snapshot.
Deltas: { "all_sites": { "expected": 505, "actual": 511 } }
Current input fingerprints: { "pairs_master": "3dabc958457eebbe", ... }
Investigate whether UniProt or the canonical tables changed before editing
the snapshot file.
```

**Do not edit the fixture to make the test pass.** The numbers in it are
snapshot expectations, not biological constants, and they never appear in
production filtering logic. Their whole job is to notice that something upstream
moved. Editing them to match discards exactly the signal the test exists to
raise.

What to do instead:

1. **Compare fingerprints.** Which input changed? Match the printed fingerprints
   against the `sha256` values recorded in the last known-good
   `results/datasets/provenance.json`. The first sixteen characters are what the test
   prints.
2. **A new UniProt snapshot** is the most common cause. Sites gain and lose
   evidence codes between releases, and counts move legitimately. That is a data
   event, not a bug.
3. **A changed canonical table** (`pairs_master.csv` or the homology QC table)
   changes the candidate universe or the subsets. Check `candidate_sites` in the
   deltas: if it moved, the input universe changed, not the evidence policy.
4. **No input changed** means the code changed. That is a genuine regression:
   fix the code, not the fixture.
5. **Only once the change is understood and intended**, freeze a *new* dated
   fixture alongside the old one, with a `_comment` recording how the numbers
   were obtained and why they moved. The dated filename is the point — the old
   snapshot stays as the record of the previous data state.

## Layout

```
analysis/experimental_glycosylation_sites/
├── README.md
├── pyproject.toml
├── config/default.toml            # paths, layer toggles, evidence policy, API settings
├── docs/
│   ├── evidence_sources.md        # what each source can and cannot establish
│   └── analysis_options.md        # downstream analyses, controls, limitations
├── src/experimental_glycosylation_sites/
│   ├── cli.py                     # four commands
│   ├── config.py                  # load TOML, resolve paths, validate existence
│   ├── models.py                  # SiteKey, UniProtFeature
│   ├── uniprot.py                 # CARBOHYD parsing from the gzipped TSV
│   ├── evidence.py                # ECO tier policy, layer combination
│   ├── glygen.py                  # protein-detail fetch, cache, category mapping
│   ├── glyconnect.py              # two-step API, site-string parsing
│   ├── orthologs.py               # candidate universe, associations, subsets
│   ├── structures.py              # chain mapping, LINK parsing, resolution ladder
│   ├── pipeline.py                # orchestration
│   ├── provenance.py              # input hashes, config hash, git state
│   └── table_io.py                # deterministic CSV writing
├── tests/                         # unit tests, fixtures, frozen snapshot
├── data/{raw,cache}/              # gitignored contents
└── results/                       # gitignored contents
```
