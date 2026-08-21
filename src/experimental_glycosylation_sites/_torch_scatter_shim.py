"""A native-torch stand-in for the two `torch_scatter` names ESM-IF imports.

`torch_scatter` is a compiled extension that must match the exact torch build.
PyG stops publishing wheels for older torch releases, and building from source
needs `nvcc` — which ARC login nodes do not provide. So on a current torch the
package is simply unavailable, and ESM-IF cannot import at all.

It is also barely needed. `esm/inverse_folding/gvp_modules.py` imports
`scatter_add` and `scatter`, uses `scatter_add` exactly once — to count how many
edges arrive at each node — and never uses `scatter`. Both have had native torch
equivalents since 1.x.

So rather than pin the whole stack to whatever torch PyG last shipped wheels for,
this provides the two names and registers itself as `torch_scatter` when the real
package is absent. `install()` is a no-op when the real one is importable, so a
machine that has it keeps using it.

The implementations are checked against the real `torch_scatter` in
`tests/test_torch_scatter_shim.py`, which skips where it is not installed.
"""
from __future__ import annotations

import sys
import types


def _broadcast(index, src, dim: int):
    """Expand a 1-D index to `src`'s shape, which is what scatter_add_ wants."""
    if index.dim() == src.dim():
        return index
    shape = [1] * src.dim()
    shape[dim] = -1
    return index.view(shape).expand_as(src)


def scatter_add(src, index, dim: int = -1, out=None, dim_size: "int | None" = None):
    """`torch_scatter.scatter_add`, via `Tensor.scatter_add_`."""
    import torch

    dim = dim if dim >= 0 else src.dim() + dim
    index = index.to(torch.long)

    if out is None:
        size = list(src.shape)
        if dim_size is not None:
            size[dim] = dim_size
        else:
            size[dim] = int(index.max()) + 1 if index.numel() else 0
        out = torch.zeros(size, dtype=src.dtype, device=src.device)

    return out.scatter_add_(dim, _broadcast(index, src, dim), src)


def scatter(src, index, dim: int = -1, out=None, dim_size: "int | None" = None,
            reduce: str = "sum"):
    """`torch_scatter.scatter`, via `Tensor.scatter_reduce_`.

    Present for import compatibility: ESM-IF imports it but never calls it.
    """
    import torch

    if reduce in ("sum", "add"):
        return scatter_add(src, index, dim, out, dim_size)

    dim = dim if dim >= 0 else src.dim() + dim
    index = index.to(torch.long)
    native = {"mean": "mean", "min": "amin", "max": "amax", "mul": "prod"}
    if reduce not in native:
        raise ValueError(f"unsupported reduce {reduce!r}")

    if out is None:
        size = list(src.shape)
        size[dim] = dim_size if dim_size is not None else (
            int(index.max()) + 1 if index.numel() else 0)
        out = torch.zeros(size, dtype=src.dtype, device=src.device)

    return out.scatter_reduce_(dim, _broadcast(index, src, dim), src,
                               reduce=native[reduce], include_self=False)


def install() -> bool:
    """Register the shim as `torch_scatter` unless the real package is present.

    Returns True if the shim was installed. Must be called before anything
    imports `esm.inverse_folding`, which pulls torch_scatter in at import time.
    """
    try:
        import torch_scatter  # noqa: F401
        return False
    except ImportError:
        pass

    module = types.ModuleType("torch_scatter")
    module.scatter_add = scatter_add
    module.scatter = scatter
    module.__doc__ = ("Native-torch shim from experimental_glycosylation_sites; "
                      "the compiled torch_scatter is not installed.")
    module.__SHIM__ = True
    sys.modules["torch_scatter"] = module
    return True
