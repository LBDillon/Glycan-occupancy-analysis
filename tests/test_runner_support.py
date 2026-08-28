"""Resumable runner plumbing that must survive interrupted output files."""

from __future__ import annotations

from experimental_glycosylation_sites.runner_support import read_resumable_csv


def test_headerless_output_from_a_failed_run_resumes_as_empty(tmp_path):
    """Stage 07 used to write a lone newline when every site failed.

    A retry then raised pandas.EmptyDataError before it could score anything,
    turning a recoverable model failure into a permanently poisoned output.
    """
    output = tmp_path / "scores.csv"
    output.write_text("\n")

    resumed = read_resumable_csv(output)

    assert resumed.empty
    assert list(resumed.columns) == []


def test_empty_resume_can_supply_the_schema_needed_for_finalisation(tmp_path):
    """An empty retry must still be safe to deduplicate and write again."""
    output = tmp_path / "scores.csv"
    output.write_text("\n")
    keys = ["accession", "position"]

    resumed = read_resumable_csv(output, empty_columns=keys)
    finalised = resumed.drop_duplicates(keys)

    assert finalised.empty
    assert list(finalised.columns) == keys
