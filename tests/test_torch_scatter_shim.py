"""The shim must agree with the real torch_scatter, where that is installed."""
from __future__ import annotations

import importlib.util

import pytest

torch = pytest.importorskip("torch")

from experimental_glycosylation_sites._torch_scatter_shim import (  # noqa: E402
    scatter, scatter_add)

HAVE_REAL = importlib.util.find_spec("torch_scatter") is not None
requires_real = pytest.mark.skipif(not HAVE_REAL, reason="torch_scatter not installed")


def test_counts_edges_like_esm_if_does():
    """The one call ESM-IF actually makes: count arrivals at each node."""
    dst = torch.tensor([0, 0, 1, 3, 3, 3])
    counts = scatter_add(torch.ones_like(dst), dst, dim_size=5)
    assert counts.tolist() == [2, 1, 0, 3, 0]


def test_dim_size_controls_the_output_length():
    src = torch.ones(4)
    index = torch.tensor([0, 0, 1, 1])
    assert scatter_add(src, index, dim_size=7).shape == (7,)


def test_empty_index_gives_zeros():
    out = scatter_add(torch.empty(0), torch.empty(0, dtype=torch.long), dim_size=3)
    assert out.tolist() == [0.0, 0.0, 0.0]


def test_broadcasts_a_1d_index_over_2d_source():
    src = torch.arange(8, dtype=torch.float).reshape(4, 2)
    index = torch.tensor([0, 0, 1, 1])
    out = scatter_add(src, index, dim=0, dim_size=2)
    assert out.tolist() == [[2.0, 4.0], [10.0, 12.0]]


@requires_real
@pytest.mark.parametrize("dim_size", [5, 8])
def test_matches_real_torch_scatter_add(dim_size):
    import torch_scatter

    torch.manual_seed(0)
    src = torch.randn(12, 3)
    index = torch.randint(0, 5, (12,))
    mine = scatter_add(src, index, dim=0, dim_size=dim_size)
    theirs = torch_scatter.scatter_add(src, index, dim=0, dim_size=dim_size)
    assert torch.equal(mine, theirs)


@requires_real
def test_matches_real_torch_scatter_on_the_esm_if_call():
    import torch_scatter

    dst = torch.randint(0, 20, (100,))
    mine = scatter_add(torch.ones_like(dst), dst, dim_size=20)
    theirs = torch_scatter.scatter_add(torch.ones_like(dst), dst, dim_size=20)
    assert torch.equal(mine, theirs)


@requires_real
@pytest.mark.parametrize("reduce", ["sum", "mean", "max"])
def test_scatter_matches_real(reduce):
    import torch_scatter

    torch.manual_seed(1)
    src = torch.randn(10, 2)
    index = torch.randint(0, 4, (10,))
    mine = scatter(src, index, dim=0, dim_size=4, reduce=reduce)
    theirs = torch_scatter.scatter(src, index, dim=0, dim_size=4, reduce=reduce)
    assert torch.allclose(mine, theirs, atol=1e-6)
