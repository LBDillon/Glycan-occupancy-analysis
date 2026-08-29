"""ESM3's two knobs: the structure track, and what is masked at the sequon.

Hermetic — a fake model and tokeniser, so none of this needs EvolutionaryScale's
`esm`, a gated checkpoint, a GPU or the network.

The tests that matter are the ones asserting the manipulation really happens:
that `seq_only` withholds the structure tokens rather than merely being labelled
that way, and that `single` and `joint` mask different things. A mode that is
named but not applied would produce a full set of plausible numbers.
"""
from __future__ import annotations

import numpy as np
import pytest

from experimental_glycosylation_sites import esm3_scoring as e3

VOCAB = ["<cls>", "<pad>", "<eos>", "<unk>", "L", "A", "G", "V", "S", "E",
         "R", "T", "I", "D", "P", "K", "Q", "N", "F", "Y", "M", "H", "W",
         "C", "X", "B", "U", "Z", "O", ".", "-", "|", "<mask>"]
MASK_ID = VOCAB.index("<mask>")


class FakeTokenizer:
    mask_token_id = MASK_ID

    def encode(self, text):
        return ([VOCAB.index("<cls>")] + [VOCAB.index(c) for c in text]
                + [VOCAB.index("<eos>")])

    def convert_ids_to_tokens(self, ids):
        return [VOCAB[i] for i in ids]

    def convert_tokens_to_ids(self, token):
        return VOCAB.index(token)


class FakeModel:
    """Records every forward call so the manipulation can be inspected."""

    def __init__(self):
        self.calls = []

    def forward(self, sequence_tokens=None, structure_tokens=None):
        import torch

        self.calls.append({"sequence": sequence_tokens.clone(),
                           "structure": None if structure_tokens is None
                           else structure_tokens.clone()})
        b, t = sequence_tokens.shape
        logits = torch.zeros((b, t, len(VOCAB)))
        return type("Out", (), {"sequence_logits": logits})()


def context_for(sequence="MANKSTV", with_structure=True):
    import torch

    tokenizer = FakeTokenizer()
    tokens = torch.tensor(tokenizer.encode(sequence))
    structure = torch.arange(len(tokens)) if with_structure else None
    return (sequence, tokens, structure), tokenizer


# --------------------------------------------------------------------------
# Provenance: the four modes must be distinguishable in the output.
# --------------------------------------------------------------------------

def test_every_mode_combination_has_its_own_conditioning_string():
    seen = {e3.conditioning(s, m)
            for s in e3.STRUCTURE_MODES for m in e3.MASK_MODES}
    assert len(seen) == 4
    assert e3.conditioning("struct_cond", "single") == \
        "masked_structure_conditioned_single"
    assert e3.conditioning("seq_only", "joint") == "masked_sequence_only_joint"


def test_the_structure_track_is_named_in_the_conditioning():
    """A structure-conditioned and a sequence-only run of the SAME model would
    otherwise be indistinguishable once written to disk."""
    assert "structure_conditioned" in e3.conditioning("struct_cond", "single")
    assert "sequence_only" in e3.conditioning("seq_only", "single")


@pytest.mark.parametrize("structure_mode,mask_mode",
                         [("nonsense", "single"), ("struct_cond", "nonsense")])
def test_an_unknown_mode_is_refused(structure_mode, mask_mode):
    with pytest.raises(ValueError):
        e3._check(structure_mode, mask_mode)


# --------------------------------------------------------------------------
# The manipulation actually happens.
# --------------------------------------------------------------------------

def test_seq_only_withholds_the_structure_tokens():
    """The whole value of this model is that the track can be switched off."""
    context, tokenizer = context_for()
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (2, 3, 4),
                                 structure_mode="seq_only", mask_mode="single")
    assert all(c["structure"] is None for c in model.calls)


def test_struct_cond_passes_the_structure_tokens():
    context, tokenizer = context_for()
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (2, 3, 4),
                                 structure_mode="struct_cond", mask_mode="single")
    assert all(c["structure"] is not None for c in model.calls)
    # one structure row per masked sequence row
    assert model.calls[0]["structure"].shape[0] == model.calls[0]["sequence"].shape[0]


def test_single_masks_one_position_per_row():
    context, tokenizer = context_for()
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (2, 3, 4),
                                 mask_mode="single")
    batch = model.calls[0]["sequence"]
    assert batch.shape[0] == 3, "one row per scored position"
    for row, index in enumerate((2, 3, 4)):
        masked = (batch[row] == MASK_ID).nonzero().flatten().tolist()
        assert masked == [index + e3.TOKEN_OFFSET]


def test_joint_masks_all_three_at_once_in_a_single_row():
    context, tokenizer = context_for()
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (2, 3, 4),
                                 mask_mode="joint")
    batch = model.calls[0]["sequence"]
    assert batch.shape[0] == 1, "joint needs one pass, not three"
    masked = (batch[0] == MASK_ID).nonzero().flatten().tolist()
    assert masked == [2 + e3.TOKEN_OFFSET, 3 + e3.TOKEN_OFFSET,
                      4 + e3.TOKEN_OFFSET]


def test_the_token_offset_is_where_the_masks_land():
    """Position i is token i+1 because <cls> is prepended. An error here masks
    a neighbouring residue and reads a plausible number for the wrong site."""
    context, tokenizer = context_for("MANKSTV")
    model = FakeModel()
    e3.conditional_probabilities(context, model, tokenizer, (0,), mask_mode="single")
    batch = model.calls[0]["sequence"]
    assert int(batch[0][0]) == VOCAB.index("<cls>")
    assert int(batch[0][1]) == MASK_ID, "position 0 must mask token 1"


def test_a_wrong_token_offset_is_caught_at_load():
    class Shifted(FakeTokenizer):
        def encode(self, text):
            return [VOCAB.index(c) for c in text]      # no <cls>
    with pytest.raises(RuntimeError, match="token offset"):
        e3._assert_token_offset(Shifted())


def test_the_correct_token_offset_is_accepted():
    e3._assert_token_offset(FakeTokenizer())


# --------------------------------------------------------------------------
# Score and guards.
# --------------------------------------------------------------------------

def rows(**overrides):
    base = {i: np.full(len(VOCAB), 1.0 / len(VOCAB)) for i in (0, 1, 2)}
    base.update(overrides)
    return base


def test_the_score_matches_a_hand_calculation():
    tokenizer = FakeTokenizer()
    n_row = np.zeros(len(VOCAB)); n_row[VOCAB.index("N")] = 0.5
    n_row[VOCAB.index("A")] = 0.5
    p1 = np.zeros(len(VOCAB)); p1[VOCAB.index("P")] = 0.1; p1[VOCAB.index("G")] = 0.9
    p2 = np.zeros(len(VOCAB)); p2[VOCAB.index("S")] = 0.5
    p2[VOCAB.index("T")] = 0.25; p2[VOCAB.index("A")] = 0.25

    score = e3.sequon_score({0: n_row, 1: p1, 2: p2}, tokenizer, 0, 1, 2)
    assert score["p_asn_at_n"] == pytest.approx(0.5)
    assert score["p_ser_or_thr_at_plus2"] == pytest.approx(0.75)
    assert score["p_pro_at_plus1"] == pytest.approx(0.1)
    assert score["conditional_sequon_score"] == pytest.approx(np.log(3.0) / 2)
    assert score["conditional_sequon_score_sd"] == 0.0
    assert score["n_decoding_orders"] == 1


def test_an_unnormalised_row_is_refused():
    bad = rows(); bad[1] = bad[1] * 3
    with pytest.raises(e3.InvalidProbabilityVector):
        e3.check_scoreable(bad, (0, 1, 2))


def test_a_missing_position_is_refused():
    with pytest.raises(e3.InvalidProbabilityVector):
        e3.check_scoreable({0: np.full(len(VOCAB), 1.0 / len(VOCAB))}, (0, 1, 2))


def test_a_triplet_mismatch_is_refused():
    e3.check_triplet("MANKSTV", (2, 3, 4), "NKS")
    with pytest.raises(e3.SequonMismatchError):
        e3.check_triplet("MANKSTV", (2, 3, 4), "NQT")


# --- design ---------------------------------------------------------------
# ESM3 can redesign a backbone, unlike ESMC. These check the plumbing that can
# be checked without the gated checkpoint: the schedule, the refusals, and that
# the sequence track really is cleared before generation rather than the model
# being handed the native sequence back.

class FakeChain:
    def __init__(self, sequence):
        self.sequence = sequence


class FakeProtein:
    """Stands in for ESMProtein, recording what generation was handed."""

    def __init__(self, chain):
        self.from_chain = chain
        self.sequence = chain.sequence
        self.coordinates = object()
        self.function_annotations = object()

    @classmethod
    def from_protein_chain(cls, chain):
        return cls(chain)


class FakeGenerationConfig:
    def __init__(self, track=None, num_steps=None, temperature=None):
        self.track, self.num_steps, self.temperature = track, num_steps, temperature


class FakeGenerator:
    """Returns a sequence that depends on the torch seed, so draws are visible."""

    def __init__(self, length, sequence=None):
        self.length, self.forced, self.seen = length, sequence, []

    def generate(self, protein, config):
        import torch

        self.seen.append({"sequence_in": protein.sequence,
                          "num_steps": config.num_steps,
                          "temperature": config.temperature,
                          "track": config.track})
        if self.forced is not None:
            return type("Out", (), {"sequence": self.forced})()
        letters = "ACDEFGHIKLMNPQRSTVWY"
        draw = torch.randint(0, len(letters), (self.length,))
        return type("Out", (), {"sequence": "".join(letters[i] for i in draw)})()


@pytest.fixture
def fake_esm(monkeypatch):
    """Inject the two esm names design_sequences imports at call time."""
    import sys
    import types

    api = types.ModuleType("esm.sdk.api")
    api.ESMProtein = FakeProtein
    api.GenerationConfig = FakeGenerationConfig
    for name, module in (("esm", types.ModuleType("esm")),
                         ("esm.sdk", types.ModuleType("esm.sdk")),
                         ("esm.sdk.api", api)):
        monkeypatch.setitem(sys.modules, name, module)
    return api


def _patch_chain(monkeypatch, sequence):
    monkeypatch.setattr(e3, "checked_chain",
                        lambda *a, **k: (sequence, FakeChain(sequence)))


def test_design_steps_are_bounded_at_both_ends():
    """One pass per residue is unaffordable; one pass total is not a joint draw."""
    assert e3.design_steps(16) == e3.DESIGN_MIN_STEPS          # short chain
    assert e3.design_steps(100_000) == e3.DESIGN_MAX_STEPS     # very long chain
    assert e3.DESIGN_MIN_STEPS < e3.design_steps(400) < e3.DESIGN_MAX_STEPS


def test_an_explicit_step_count_never_exceeds_the_chain():
    assert e3.design_steps(30, num_steps=12) == 12
    assert e3.design_steps(10, num_steps=999) == 10


def test_design_clears_the_sequence_track_before_generating(monkeypatch, fake_esm):
    """The native sequence must not be handed back to the model.

    If it were, the model would be completing a sequence it can already see and
    retention would measure copying rather than redesign.
    """
    sequence = "MANKSTVQW"
    _patch_chain(monkeypatch, sequence)
    model = FakeGenerator(len(sequence))

    e3.design_sequences("x.pdb", "A", model, n_designs=2, temperature=0.1)

    assert len(model.seen) == 2
    assert all(call["sequence_in"] is None for call in model.seen)
    assert all(call["track"] == "sequence" for call in model.seen)
    assert all(call["temperature"] == 0.1 for call in model.seen)


def test_designs_come_back_full_length_and_one_per_request(monkeypatch, fake_esm):
    sequence = "MANKSTVQW"
    _patch_chain(monkeypatch, sequence)

    designs = e3.design_sequences("x.pdb", "A", FakeGenerator(len(sequence)),
                                  n_designs=5, temperature=0.1)

    assert len(designs) == 5
    assert all(len(d) == len(sequence) for d in designs)


def test_the_same_seed_reproduces_the_same_designs(monkeypatch, fake_esm):
    """A resumed run must not silently redraw a site it already has."""
    sequence = "MANKSTVQW"
    _patch_chain(monkeypatch, sequence)
    kwargs = dict(n_designs=4, temperature=0.5, seed=7)

    first = e3.design_sequences("x.pdb", "A", FakeGenerator(len(sequence)), **kwargs)
    again = e3.design_sequences("x.pdb", "A", FakeGenerator(len(sequence)), **kwargs)

    assert first == again
    assert len(set(first)) > 1, "every design identical: the seed is not advancing"


def test_a_short_sequence_from_the_model_is_refused(monkeypatch, fake_esm):
    """Silently accepting it would shift every index in the retention read-out."""
    sequence = "MANKSTVQW"
    _patch_chain(monkeypatch, sequence)
    truncated = FakeGenerator(len(sequence), sequence="MANK")

    with pytest.raises(e3.DesignFailedError, match="4 residues"):
        e3.design_sequences("x.pdb", "A", truncated, n_designs=1, temperature=0.1)


def test_an_empty_sequence_from_the_model_is_refused(monkeypatch, fake_esm):
    _patch_chain(monkeypatch, "MANKSTVQW")

    with pytest.raises(e3.DesignFailedError, match="no sequence"):
        e3.design_sequences("x.pdb", "A", FakeGenerator(9, sequence=""),
                            n_designs=1, temperature=0.1)


def test_a_generation_failure_is_not_a_parse_failure():
    """They are counted separately, so they must not share a type."""
    assert not issubclass(e3.DesignFailedError, e3.ChainUnreadableError)
    assert not issubclass(e3.ChainUnreadableError, e3.DesignFailedError)


def test_the_adapter_refuses_to_design_without_the_structure_track():
    """seq_only has no backbone, so a 'redesign' would be an unconditional draw."""
    from experimental_glycosylation_sites.adapters.esm3 import ESM3Adapter

    adapter = ESM3Adapter(structure_mode="seq_only")

    with pytest.raises(ValueError, match="structure track"):
        adapter.design("x.pdb", "A", n_designs=4, temperature=0.1)


def test_the_adapter_records_a_non_default_schedule_in_its_provenance():
    """Designs made under different schedules are not the same estimand."""
    from experimental_glycosylation_sites.adapters.esm3 import ESM3Adapter

    assert "design_num_steps" not in ESM3Adapter().describe()
    assert ESM3Adapter(num_steps=12).describe()["design_num_steps"] == 12


def test_designs_are_not_labelled_autoregressive():
    """Stage 08's default label would be wrong for a masked diffusion model.

    Retention measured by left-to-right decoding and by iterative unmasking are
    different quantities; this column is what stops them being pooled.
    """
    from experimental_glycosylation_sites.adapters.esm3 import ESM3Adapter

    generation = ESM3Adapter().describe_generation()

    assert generation["generation"] == "masked_diffusion_unmasking"
    assert "autoregressive" not in generation["generation"]
    assert generation["native_procedure"] is True


def test_the_adapter_accepts_max_batch_without_claiming_to_use_it():
    """Stage 08 passes it to every adapter; ESM3 generates one design at a time."""
    from experimental_glycosylation_sites.adapters.esm3 import ESM3Adapter

    adapter = ESM3Adapter(max_batch=8)

    assert adapter.max_batch == 8
    assert "max_batch" not in adapter.describe_generation()


# --- batched design -------------------------------------------------------
# One design at a time left the GPU idle: a pilot measured ~13 unmasking steps
# per second whatever the chain length, so cost was steps x designs and barely
# touched by length. These check the batching without needing a real model.

class FakeBatchGenerator(FakeGenerator):
    """Exposes batch_generate, and records the size of every batch it is given."""

    def __init__(self, length, sequence=None, short_by=0):
        super().__init__(length, sequence)
        self.batches, self.short_by = [], short_by

    def batch_generate(self, proteins, configs):
        assert len(proteins) == len(configs), "one config per protein"
        self.batches.append(len(proteins))
        outs = [self.generate(p, c) for p, c in zip(proteins, configs)]
        return outs[:len(outs) - self.short_by] if self.short_by else outs


def test_batching_is_used_when_the_model_offers_it(monkeypatch, fake_esm):
    sequence = "MANKSTVQW"
    _patch_chain(monkeypatch, sequence)
    model = FakeBatchGenerator(len(sequence))

    designs = e3.design_sequences("x.pdb", "A", model, n_designs=32,
                                  temperature=0.1)

    assert len(designs) == 32
    assert model.batches, "batch_generate was never called"
    assert sum(model.batches) == 32


def test_a_long_chain_is_batched_within_the_slot_budget(monkeypatch, fake_esm):
    """The same bound that stopped ProteinMPNN exhausting memory on long chains."""
    from experimental_glycosylation_sites.retention import (DESIGN_SLOT_BUDGET,
                                                            batch_for_length)

    long_sequence = "A" * 1200
    _patch_chain(monkeypatch, long_sequence)
    model = FakeBatchGenerator(len(long_sequence))

    e3.design_sequences("x.pdb", "A", model, n_designs=32, temperature=0.1)

    expected = batch_for_length(len(long_sequence), 32)
    assert expected < 32, "a 1200-residue chain should not run 32 at once"
    assert max(model.batches) <= expected
    assert max(model.batches) * len(long_sequence) <= DESIGN_SLOT_BUDGET


def test_a_short_chain_uses_one_batch(monkeypatch, fake_esm):
    _patch_chain(monkeypatch, "MANKSTVQW")
    model = FakeBatchGenerator(9)

    e3.design_sequences("x.pdb", "A", model, n_designs=32, temperature=0.1)

    assert model.batches == [32]


def test_max_batch_caps_the_batch(monkeypatch, fake_esm):
    _patch_chain(monkeypatch, "MANKSTVQW")
    model = FakeBatchGenerator(9)

    e3.design_sequences("x.pdb", "A", model, n_designs=32, temperature=0.1,
                        max_batch=8)

    assert max(model.batches) <= 8
    assert sum(model.batches) == 32


def test_a_batch_returning_too_few_designs_is_refused(monkeypatch, fake_esm):
    """Silently accepting it would score retention over fewer designs than 32."""
    _patch_chain(monkeypatch, "MANKSTVQW")
    model = FakeBatchGenerator(9, short_by=1)

    with pytest.raises(e3.DesignFailedError, match="designs for a batch"):
        e3.design_sequences("x.pdb", "A", model, n_designs=4, temperature=0.1)


def test_batching_can_be_turned_off(monkeypatch, fake_esm):
    """The sequential path must stay reachable: it needs only `generate`."""
    _patch_chain(monkeypatch, "MANKSTVQW")
    model = FakeBatchGenerator(9)

    designs = e3.design_sequences("x.pdb", "A", model, n_designs=4,
                                  temperature=0.1, use_batch=False)

    assert len(designs) == 4
    assert model.batches == [], "batch_generate was used despite use_batch=False"


def test_a_model_without_batch_generate_still_works(monkeypatch, fake_esm):
    _patch_chain(monkeypatch, "MANKSTVQW")

    designs = e3.design_sequences("x.pdb", "A", FakeGenerator(9),
                                  n_designs=4, temperature=0.1)

    assert len(designs) == 4


def test_which_path_was_used_is_recorded(monkeypatch, fake_esm):
    """The two paths draw the same distribution but consume randomness
    differently, so a run is reproducible against itself and not across them."""
    from experimental_glycosylation_sites.adapters.esm3 import ESM3Adapter

    assert ESM3Adapter().describe_generation()["batched"] is True
    assert ESM3Adapter(use_batch=False).describe_generation()["batched"] is False
