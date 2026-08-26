# Glyco-site context analysis

What the sequence and structural environments of experimentally supported
N-linked glycosylation sites actually look like, and whether they differ from
matched unoccupied sequons.

Separate from the occupancy benchmark in the parent directory, which asks a
different question — *what do design models distinguish?* This asks *what does
the biology look like?* The two meet at Step 4, where context differences are
tested against the frozen model contrasts.

## Layout

    pipeline/   50  the natural reference picture
                51  redesign with the sequon held fixed, and measure local chemistry
                52  distance from natural occupied context, and the random control
                53  result figures          54  feature distributions
                55  composition control: is the shift local or global?

    src/glyco_context/   local_chemistry   the fifteen-feature panel
                         fixed_design      holding positions fixed, and verifying it
                         context_distance  distance from the reference, protein held out
                         context_stats     effect sizes, cluster bootstrap, BH

    docs/       why this exists, the pre-specification, the findings
    archive/    context_extractor/     the 97-column structural description
                comparative_analysis/  occupied against unoccupied, a negative result
    results/    generated tables, reports and figures (git-ignored)
    tests/      run from the repository root with the rest of the suite

Two things are archived rather than deleted, each with a README explaining why.
The extractor's output is retained and still read; only its code is archived.

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
