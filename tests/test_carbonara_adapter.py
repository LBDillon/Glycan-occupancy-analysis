"""The CARBonAra adapter's contract with the pipeline.

Hermetic throughout: the model is a fake and the upstream parser a stub, so none
of this needs a checkout, weights, gemmi, blosum, CUDA or the network.

What these tests are for is the seam rather than the arithmetic — that stage 05
gets manifest-space booleans, that stage 07 gets the same row schema every other
model produces, and that the package still imports when CARBonAra is absent,
which is the normal state of every environment except the one running it.
"""
from __future__ import annotations

import ast
import json
import subprocess
import sys

import numpy as np
import pytest

from experimental_glycosylation_sites import adapters, carbonara_scoring as cs
from experimental_glycosylation_sites.adapters import carbonara as carbonara_adapter
from experimental_glycosylation_sites.adapters.base import (SequenceDesigner,
                                                            SequonScorer)
from test_carbonara_scoring import (BACKBONE, FakeCARBonAra, SEQUON_CHAIN,
                                    SEQUON_SEQ, build_pdb, simple_chain,
                                    stub_structure, write)


@pytest.fixture
def fake_backend(monkeypatch):
    """A CARBonAra whose parse and model are both ours."""
    monkeypatch.setattr(cs, "_load_structure",
                        lambda pdb_text, carbonara_dir=None: stub_structure(pdb_text))
    model = FakeCARBonAra(len(SEQUON_CHAIN))
    monkeypatch.setattr(carbonara_adapter, "load_model",
                        lambda *a, **k: model)
    return model


@pytest.fixture
def chain(tmp_path):
    return write(tmp_path, build_pdb(simple_chain(SEQUON_CHAIN)))


# --------------------------------------------------------------------------
# 1. Registration and laziness.
# --------------------------------------------------------------------------

def test_carbonara_is_registered():
    assert "carbonara" in adapters.available()


def test_importing_the_package_does_not_import_carbonara():
    """A registry that imported every model would break on the absent one.

    Run in a subprocess because another test in this session may already have
    imported the adapter, which would make an in-process check pass for the
    wrong reason.
    """
    script = (
        "import sys; sys.path.insert(0, 'src');"
        "import experimental_glycosylation_sites as p;"
        "from experimental_glycosylation_sites import adapters;"
        "names = adapters.available();"
        "assert 'carbonara' in names, names;"
        "leaked = [m for m in ('carbonara', 'gemmi', 'blosum') if m in sys.modules];"
        "assert not leaked, leaked;"
        "assert 'experimental_glycosylation_sites.adapters.carbonara'"
        " not in sys.modules;"
        "print('ok')"
    )
    result = subprocess.run([sys.executable, "-c", script],
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "ok" in result.stdout


def test_constructing_the_adapter_does_not_need_a_checkout():
    """Discovery is deferred to first use, so a missing checkout surfaces when
    the model is actually wanted rather than at import."""
    adapter = carbonara_adapter.CARBonAraAdapter()
    assert adapter.name == "carbonara"
    assert adapter._model is None


# --------------------------------------------------------------------------
# 2. Protocol conformance.
# --------------------------------------------------------------------------

def test_the_adapter_is_a_scorer_and_not_a_designer():
    """Scorer only, by design: upstream generation samples stochastically from
    raw confidences, so there is no retention number this benchmark can use."""
    adapter = carbonara_adapter.CARBonAraAdapter()
    assert isinstance(adapter, SequonScorer)
    assert not isinstance(adapter, SequenceDesigner)
    assert not hasattr(adapter, "design")


def test_the_adapter_implements_the_whole_scorer_protocol():
    adapter = carbonara_adapter.CARBonAraAdapter()
    for method in ("describe", "decodable_positions", "prepare_chain",
                   "score_from", "score_site"):
        assert callable(getattr(adapter, method)), method


def test_describe_names_the_conditioning_precisely():
    """It must never be poolable with ProteinMPNN's average over eight sampled
    decoding orders, nor with ESM-IF's causal prefix."""
    described = carbonara_adapter.CARBonAraAdapter().describe()
    assert described["conditioning"] == "conditional_all_other_native"
    assert described["n_orders"] == 1
    assert described["model"] == cs.DEFAULT_MODEL


# --------------------------------------------------------------------------
# Checkout discovery.
# --------------------------------------------------------------------------

def test_checkout_discovery_honours_the_environment_override(tmp_path, monkeypatch):
    from experimental_glycosylation_sites.runner_support import carbonara_dir

    (tmp_path / "carbonara.py").write_text("# stand-in\n")
    monkeypatch.setenv("CARBONARA_DIR", str(tmp_path))
    assert carbonara_dir() == tmp_path


def test_a_half_cloned_checkout_is_rejected(tmp_path, monkeypatch):
    """Directory existence is not evidence of a usable checkout.

    The fallback candidates are emptied for the duration: whether this machine
    happens to have a checkout at `~/CARBonAra` is not what is under test, and
    leaving them in place made this pass or fail depending on the developer's
    home directory.
    """
    from experimental_glycosylation_sites import runner_support

    monkeypatch.setattr(runner_support, "CARBONARA_CANDIDATES", ())
    monkeypatch.setenv("CARBONARA_DIR", str(tmp_path))
    with pytest.raises(FileNotFoundError, match="carbonara.py"):
        runner_support.carbonara_dir()


def test_discovery_reports_every_place_it_looked(tmp_path, monkeypatch):
    """The error has to be actionable without reading the source."""
    from experimental_glycosylation_sites import runner_support

    monkeypatch.delenv("CARBONARA_DIR", raising=False)
    monkeypatch.setattr(runner_support, "CARBONARA_CANDIDATES",
                        (str(tmp_path / "nowhere"),))
    with pytest.raises(FileNotFoundError) as caught:
        runner_support.carbonara_dir()
    message = str(caught.value)
    assert "nowhere" in message
    assert "CARBONARA_DIR" in message


# --------------------------------------------------------------------------
# Against a real checkout, when one happens to be present.
# --------------------------------------------------------------------------

def _checkout():
    from experimental_glycosylation_sites.runner_support import carbonara_dir

    try:
        return carbonara_dir()
    except FileNotFoundError:
        return None


def test_the_assumed_alphabet_matches_the_checkout():
    """The regression test for the defect this module is most exposed to.

    `CARBONARA_RESNAMES` is a constant transcribed from upstream, and every
    probability column depends on its order. Compared against the checkout's own
    `std_aminoacids` whenever there is a checkout to compare against — no model
    load, no weights, just the source of truth.
    """
    checkout = _checkout()
    if checkout is None:
        pytest.skip("no CARBonAra checkout")
    encoding = checkout / "src" / "data_encoding.py"
    if not encoding.exists():
        pytest.skip("checkout has no src/data_encoding.py")

    # Parsed rather than imported: importing the module would pull in torch and
    # gemmi, and the literal is what needs checking, not the import machinery.
    tree = ast.parse(encoding.read_text())
    upstream = None
    for node in ast.walk(tree):
        if not isinstance(node, ast.Assign):
            continue
        if not any(getattr(t, "id", None) == "std_aminoacids" for t in node.targets):
            continue
        # np.array([...]) -> the list literal inside the call
        value = node.value
        if isinstance(value, ast.Call) and value.args:
            value = value.args[0]
        upstream = tuple(ast.literal_eval(value))
        break

    assert upstream is not None, f"no std_aminoacids assignment found in {encoding}"
    assert upstream == cs.CARBONARA_RESNAMES, (
        f"checkout says {upstream}, this module assumes {cs.CARBONARA_RESNAMES}")
    cs.verify_alphabet(upstream)


# --------------------------------------------------------------------------
# 15. Stage 05 — scoreability in manifest space.
# --------------------------------------------------------------------------

def test_decodable_positions_is_a_manifest_space_boolean_array(chain, fake_backend):
    """`05_scoreability.py` does `decodable[int(r.n_model_index)]`."""
    adapter = carbonara_adapter.CARBonAraAdapter()
    decodable = adapter.decodable_positions(chain, "A")
    assert decodable.dtype == bool
    assert len(decodable) == len(SEQUON_SEQ)
    assert decodable.all()


def test_stage_05_indexing_selects_the_expected_site(tmp_path, fake_backend):
    """The exact expression stage 05 evaluates, against a chain with one
    unscoreable residue in the middle of the sequon."""
    residues = [("ALA", 1, " ", BACKBONE, "ATOM  "),
                ("GLY", 2, " ", BACKBONE, "ATOM  "),
                ("ASN", 3, " ", BACKBONE, "ATOM  "),
                ("LYS", 4, " ", ("N", "CA", "C"), "ATOM  "),
                ("SER", 5, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    decodable = carbonara_adapter.CARBonAraAdapter().decodable_positions(path, "A")

    indices = (2, 3, 4)
    bad = [i for i in indices if i >= len(decodable) or not bool(decodable[i])]
    assert bad == [3]


def test_an_unreadable_chain_is_all_false_not_an_exception(chain, fake_backend):
    adapter = carbonara_adapter.CARBonAraAdapter()
    assert not adapter.decodable_positions(chain, "Q").any()


# --------------------------------------------------------------------------
# 16. Stage 07 — the row schema.
# --------------------------------------------------------------------------

EXPECTED_SCORE_KEYS = {
    "conditional_sequon_score", "conditional_sequon_score_sd",
    "n_decoding_orders", "p_asn_at_n", "p_ser_at_plus2", "p_thr_at_plus2",
    "p_ser_or_thr_at_plus2", "p_pro_at_plus1", "logit_p_asn",
    "logit_p_ser_or_thr", "probs_n", "probs_plus1", "probs_plus2",
}


def test_score_site_returns_the_shared_schema(chain, fake_backend):
    adapter = carbonara_adapter.CARBonAraAdapter()
    scored = adapter.score_site(chain, "A", (2, 3, 4), expected_triplet="NKS")
    assert set(scored) == EXPECTED_SCORE_KEYS


def test_the_schema_matches_the_other_models(chain, fake_backend):
    """Downstream is model-agnostic, so a column the others lack would either
    break a merge or quietly become a column of NaN."""
    from experimental_glycosylation_sites.esmif_scoring import sequon_score as esm

    class Alphabet:
        @staticmethod
        def get_idx(letter):
            return "ACDEFGHIKLMNPQRSTVWY".index(letter)

    rows = np.full((3, 20), 1.0 / 20)
    reference = esm(rows, Alphabet(), 0, 1, 2)

    adapter = carbonara_adapter.CARBonAraAdapter()
    scored = adapter.score_site(chain, "A", (2, 3, 4))
    assert set(scored) == set(reference)


def test_stage_07_can_serialise_the_row(chain, fake_backend):
    """The exact transformation `07_score.py` applies before writing."""
    adapter = carbonara_adapter.CARBonAraAdapter()
    provenance = adapter.describe()
    scored = adapter.score_site(chain, "A", (2, 3, 4), expected_triplet="NKS")

    vectors = {k: json.dumps([round(x, 6) for x in scored.pop(k)])
               for k in ("probs_n", "probs_plus1", "probs_plus2")}
    row = {"accession": "P00000", "position": 3, "structure_pdb_id": "TEST",
           "structure_chain_id": "A", "triplet": "NKS", "subtype": "NXS",
           **provenance, **scored, **vectors}

    for key in ("probs_n", "probs_plus1", "probs_plus2"):
        assert len(json.loads(row[key])) == 20
    assert row["conditioning"] == "conditional_all_other_native"
    assert np.isfinite(row["conditional_sequon_score"])
    json.dumps(row)          # the whole row must survive a round trip


def test_probability_vectors_are_in_carbonaras_documented_order(chain, fake_backend):
    """Twenty entries, abundance-sorted, N at 13 — not twenty-one and not
    alphabetical. Anything reading these columns needs the order to be stated."""
    adapter = carbonara_adapter.CARBonAraAdapter()
    scored = adapter.score_site(chain, "A", (2, 3, 4))
    assert len(scored["probs_n"]) == 20
    assert scored["probs_n"][cs.AA_INDEX["N"]] == pytest.approx(scored["p_asn_at_n"])
    assert scored["probs_plus2"][cs.AA_INDEX["S"]] == pytest.approx(
        scored["p_ser_at_plus2"])
    assert scored["probs_plus1"][cs.AA_INDEX["P"]] == pytest.approx(
        scored["p_pro_at_plus1"])


# --------------------------------------------------------------------------
# prepare_chain / score_from, as the runners use them.
# --------------------------------------------------------------------------

def test_prepare_chain_computes_only_the_requested_positions(chain, fake_backend):
    """CARBonAra pays one forward pass per position, so this is the difference
    between three passes and one per residue on the chain."""
    adapter = carbonara_adapter.CARBonAraAdapter()
    adapter.prepare_chain(chain, "A", [2, 3, 4])
    assert len(fake_backend.calls) == 3


def test_a_second_sequon_on_a_chain_reuses_the_prepared_context(chain, fake_backend):
    adapter = carbonara_adapter.CARBonAraAdapter()
    context = adapter.prepare_chain(chain, "A", [2, 3, 4])
    before = len(fake_backend.calls)
    adapter.score_from(context, (2, 3, 4), expected_triplet="NKS")
    adapter.score_from(context, (2, 3, 4), expected_triplet="NKS")
    assert len(fake_backend.calls) == before


def test_score_from_refuses_a_triplet_the_model_did_not_read(chain, fake_backend):
    adapter = carbonara_adapter.CARBonAraAdapter()
    context = adapter.prepare_chain(chain, "A", [2, 3, 4])
    with pytest.raises(cs.SequonMismatchError):
        adapter.score_from(context, (2, 3, 4), expected_triplet="NQT")


def test_score_from_refuses_a_position_that_was_never_evaluated(chain, fake_backend):
    """An unscoreable site must fail closed, not borrow a neighbour's row."""
    adapter = carbonara_adapter.CARBonAraAdapter()
    context = adapter.prepare_chain(chain, "A", [2])
    with pytest.raises(cs.IncompleteBackboneError):
        adapter.score_from(context, (2, 3, 4))


def test_an_incomplete_backbone_site_is_dropped_before_the_model_runs(
        tmp_path, fake_backend):
    residues = [("ALA", 1, " ", BACKBONE, "ATOM  "),
                ("GLY", 2, " ", BACKBONE, "ATOM  "),
                ("ASN", 3, " ", ("N", "CA"), "ATOM  "),
                ("LYS", 4, " ", BACKBONE, "ATOM  "),
                ("SER", 5, " ", BACKBONE, "ATOM  ")]
    path = write(tmp_path, build_pdb(residues))
    adapter = carbonara_adapter.CARBonAraAdapter()
    context = adapter.prepare_chain(path, "A", [2, 3, 4])
    assert 2 not in context[1]
    with pytest.raises(cs.IncompleteBackboneError):
        adapter.score_from(context, (2, 3, 4))
