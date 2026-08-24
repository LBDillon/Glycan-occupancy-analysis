"""External input locations resolved by configuration, not by relative guesses.

The repository holds code, tests and generated tables; the bulk inputs (the
UniProt release, the structure cache) live outside it and are read only. Hard
coding `../../data/...` tied every stage to one checkout layout and made the
repository non-portable -- the exact failure that stopped the corrected rerun.
"""
from __future__ import annotations

import pytest

from experimental_glycosylation_sites.input_paths import (
    MissingInput,
    resolve_input,
    structure_dirs,
)


def test_resolve_input_finds_a_file_under_a_configured_root(tmp_path, monkeypatch):
    release = tmp_path / "raw" / "uniprot" / "release.tsv.gz"
    release.parent.mkdir(parents=True)
    release.write_text("")
    monkeypatch.setenv("GCA_DATA_ROOTS", str(tmp_path))
    assert resolve_input("raw/uniprot/release.tsv.gz") == release


def test_resolve_input_prefers_the_first_configured_root(tmp_path, monkeypatch):
    """Roots are searched in order, so an override can shadow the default."""
    first, second = tmp_path / "a", tmp_path / "b"
    for root in (first, second):
        (root / "raw").mkdir(parents=True)
        (root / "raw" / "x.tsv").write_text(root.name)
    monkeypatch.setenv("GCA_DATA_ROOTS", f"{first}:{second}")
    assert resolve_input("raw/x.tsv").read_text() == "a"


def test_missing_input_names_every_location_searched(tmp_path, monkeypatch):
    """A stage that cannot find its input must say where it looked.

    The alternative is what happened before: a relative path that silently
    resolves to nothing outside its original working directory.
    """
    monkeypatch.setenv("GCA_DATA_ROOTS", str(tmp_path))
    with pytest.raises(MissingInput) as excinfo:
        resolve_input("raw/uniprot/absent.tsv.gz")
    message = str(excinfo.value)
    assert "absent.tsv.gz" in message
    assert str(tmp_path) in message
    assert "GCA_DATA_ROOTS" in message, "the error should say how to configure it"


def test_structure_dirs_puts_configured_directories_first(tmp_path, monkeypatch):
    configured = tmp_path / "cache"
    configured.mkdir()
    monkeypatch.setenv("GCA_STRUCTURE_DIRS", str(configured))
    assert structure_dirs()[0] == configured
