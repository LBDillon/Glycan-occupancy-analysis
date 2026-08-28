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
    # Scorer and designer, though design() uses a custom
    # independent_calibrated_sampling path rather than upstream's stochastic
    # imprint_sampling. Needs a CARBonAra checkout plus gemmi, blosum,
    # scikit-learn and h5py, none of them core dependencies.
    "carbonara": ("experimental_glycosylation_sites.adapters.carbonara",
                  "CARBonAraAdapter"),
    # Causal and sequence-only. Scorer only: generation is unconditioned by any
    # backbone, so there is no redesign of a chain to measure retention on.
    # Needs `transformers`, which is not a core dependency.
    "progen2": ("experimental_glycosylation_sites.adapters.progen2",
                "ProGen2Adapter"),
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
