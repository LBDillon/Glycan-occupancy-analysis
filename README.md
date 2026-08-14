# experimental_glycosylation_sites

Site-level N-linked glycosylation occupancy evidence derived from the ortholog
sequon-conservation database.

- Evidence source guide: [`docs/evidence_sources.md`](docs/evidence_sources.md)
- Downstream analysis options: [`docs/analysis_options.md`](docs/analysis_options.md)

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

`occupancy_status` therefore takes only two values in this phase:
`occupied_supported` and `unknown`. A third value, `observed_unmodified`, is
defined in `evidence.py` and deliberately left unpopulated: no current source
can establish that a site was *examined and found bare*. A test asserts it is
never emitted.

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

From the latest full run (`results/summary.json`):

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
   `results/provenance.json`. The first sixteen characters are what the test
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
