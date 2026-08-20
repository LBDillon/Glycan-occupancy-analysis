"""Model adapters. Register new models here."""
from .base import SequenceDesigner, SequonScorer

# Imported lazily: each adapter pulls in its model's dependencies, and a missing
# one should disable that model rather than break the package.
_REGISTRY = {
    "proteinmpnn": ("experimental_glycosylation_sites.adapters.proteinmpnn",
                    "ProteinMPNNAdapter"),
    "esm_if": ("experimental_glycosylation_sites.adapters.esm_if",
               "ESMIFAdapter"),
    # Sequence-only. Needs EvolutionaryScale's `esm`, which cannot be installed
    # alongside `fair-esm` -- both claim the import name `esm`.
    "esmc": ("experimental_glycosylation_sites.adapters.esmc", "ESMCAdapter"),
}


def available() -> list[str]:
    return sorted(_REGISTRY)


def load(name: str, **kwargs):
    import importlib

    if name not in _REGISTRY:
        raise KeyError(f"unknown model {name!r}; registered: {available()}")
    module_path, class_name = _REGISTRY[name]
    return getattr(importlib.import_module(module_path), class_name)(**kwargs)


__all__ = ["SequonScorer", "SequenceDesigner", "available", "load"]
