"""Chunking the marginalisation batch must be arithmetic, not approximation.

The whole structure is tiled once per sample, so a batch of 16 is thousands of
residue-slots on a long chain — which took out seven of eight ARC tasks with
OUT_OF_MEMORY. Splitting the batch fixes that only if the split changes nothing
about the answer: same rows, same order, same values.

The sampling that makes marginalisation stochastic happens outside the forward
pass, so chunking cannot perturb it. These cover the part that could go wrong.
"""
from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")


def chunked_forward(token_batch, chunk, run):
    """The pattern used in `marginalised_probabilities`."""
    total = token_batch.shape[0]
    size = total if chunk is None else max(1, min(total, chunk))
    pieces = [run(token_batch[start:start + size])
              for start in range(0, total, size)]
    return torch.cat(pieces, dim=0) if len(pieces) > 1 else pieces[0]


def _run(batch):
    """A stand-in whose output depends on the row, so reordering would show."""
    return batch.float().unsqueeze(-1) * torch.tensor([1.0, 2.0])


def test_chunking_reproduces_the_unchunked_result_exactly():
    batch = torch.arange(16 * 5).reshape(16, 5)
    whole = chunked_forward(batch, None, _run)
    for chunk in (1, 2, 3, 4, 8, 16, 32):
        assert torch.equal(chunked_forward(batch, chunk, _run), whole), chunk


def test_row_order_is_preserved():
    batch = torch.arange(9 * 2).reshape(9, 2)
    out = chunked_forward(batch, 4, _run)
    assert torch.equal(out[0], _run(batch[:1])[0])
    assert torch.equal(out[-1], _run(batch[-1:])[0])


def test_a_chunk_larger_than_the_batch_is_one_pass():
    batch = torch.arange(3 * 2).reshape(3, 2)
    assert torch.equal(chunked_forward(batch, 99, _run), _run(batch))


def test_uneven_final_chunk_is_not_dropped():
    batch = torch.arange(7 * 2).reshape(7, 2)
    assert chunked_forward(batch, 3, _run).shape[0] == 7
