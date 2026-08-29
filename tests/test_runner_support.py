"""Resumable runner plumbing that must survive interrupted output files."""

from __future__ import annotations

from experimental_glycosylation_sites.runner_support import (
    needs_header, read_resumable_csv)


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


def test_a_zero_byte_output_still_asks_for_a_header(tmp_path):
    """The root cause of the headerless append.

    A run killed before its first flush leaves a zero-byte file. `not
    path.exists()` is False for it, so the append wrote no header and the rows
    became unreadable.
    """
    output = tmp_path / "scores.csv"
    output.touch()

    assert needs_header(output)
    assert needs_header(tmp_path / "absent.csv")


def test_a_written_output_does_not_repeat_its_header(tmp_path):
    output = tmp_path / "scores.csv"
    output.write_text("accession,position\nP01861,42\n")

    assert not needs_header(output)


def test_a_headerless_append_is_refused_rather_than_silently_mistyped(tmp_path):
    """The failure this pairing exists to prevent.

    Appending without a header leaves a file that parses cleanly and takes its
    first DATA row as the column names. It surfaced as a bare KeyError from
    drop_duplicates at the very end of a stage -- after 225 of 248 chain groups
    had been scored -- losing the whole run's compute to a naming accident.
    """
    output = tmp_path / "scores.csv"
    output.write_text("P01861,42\nP02671,17\n")   # no header line

    try:
        read_resumable_csv(output, empty_columns=["accession", "position"])
    except SystemExit as exc:
        assert "accession" in str(exc) and "header" in str(exc)
    else:
        raise AssertionError("a headerless file was accepted")


def test_the_append_cycle_round_trips_from_a_zero_byte_file(tmp_path):
    """What stages 07, 08 and 08b actually do, end to end."""
    import pandas as pd

    output = tmp_path / "scores.csv"
    output.touch()                                  # crashed run's leftover
    keys = ["accession", "position"]

    for rows in ([{"accession": "P01861", "position": 42}],
                 [{"accession": "P02671", "position": 17}]):
        pd.DataFrame(rows).to_csv(output, mode="a",
                                  header=needs_header(output), index=False)

    resumed = read_resumable_csv(output, empty_columns=keys)

    assert list(resumed.columns) == keys
    assert len(resumed.drop_duplicates(keys)) == 2
