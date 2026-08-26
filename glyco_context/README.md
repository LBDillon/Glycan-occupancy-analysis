# Glyco-site context analysis

What the sequence and structural environments of experimentally supported
N-linked glycosylation sites actually look like, and whether they differ from
matched unoccupied sequons.

Separate from the occupancy benchmark in the parent directory, which asks a
different question — *what do design models distinguish?* This asks *what does
the biology look like?* The two meet at Step 4, where context differences are
tested against the frozen model contrasts.

## Layout

    pipeline/   41  build the context manifest (one row per biological site)
                43  extract structural features   43b  merge shards, gated
                44  the three analysis views      45   enforcing QC report
                46  old-vs-new change audit
                47  describe the occupied distribution        (Step 2)
                48  population comparisons                    (Step 3)
                49  matched-pair comparisons                  (Amendment 1)

    src/glyco_context/   context_features  extraction
                         context_merge     shard merging, fatal on loss
                         context_views     the three views and exclusion reasons
                         context_qc        invariants
                         context_stats     effect sizes, cluster bootstrap, BH
                         sequence_qc       row-level sequence checks
                         change_audit      attribute every changed row

    docs/       why this exists, the mapping correction, the pre-specification,
                the findings, plus ../docs/glossary.md (every term defined) and
                figure_captions.md (figures carry no explanatory text)
    results/    generated tables, reports and figures (git-ignored)
    tests/      run from the repository root with the rest of the suite

## Running it

Stages read the shared inputs by configuration, so run from the **repository
root** with the input roots set:

    export GCA_DATA_ROOTS=/path/to/data:/path/to/analysis/data
    export GCA_DATASETS_DIR=/path/to/results/datasets
    export GCA_STRUCTURE_DIRS=/path/to/pdb:/path/to/more/pdb
    python glyco_context/pipeline/44_context_views.py

## What it depends on

Not standalone, by design. It consumes the occupancy resource's evidence tables
and shares its structure handling — `experimental_glycosylation_sites.features`
for the accessibility model (so RSA stays on the same scale as the frozen
matching), plus `runner_support`, `provenance` and `input_paths`. Keeping one
accessibility implementation matters more than package independence: two would
drift, and the benchmark's matching is built on this one.

## Read this before using the results

`docs/findings_2026-08-26_context_retention.md` — the live result. Designs drift
away from natural occupied context with the sequon fully protected (+0.071
[0.025, 0.122]), and ProteinMPNN is not measurably better at preserving it than
changing the same number of residues at random.

The earlier comparative analysis is archived in
`archive/comparative_analysis/`. Once composition is controlled by matching,
nearly every context difference between occupied and unoccupied sequons
disappears, and the population-level comparison is confounded with protein
identity. Do not quote those population numbers as occupancy effects.
