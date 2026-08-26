"""Redesign glycoprotein chains with the sequon held fixed, and score the context.

The experiment: force ProteinMPNN to keep a naturally occupied N-X-S/T motif,
then ask what it does to the residues around it. See
docs/prespecification_fixed_sequon_context_retention.md.

Three variants per site, all on the same backbone so only sequence differs:

    wild_type   the deposited sequence
    design      ProteinMPNN with the sequon positions masked out of design
    random      the same number of positions altered, residues drawn from the
                chain's own composition -- so the control is arbitrary change of
                matched size, not a composition-destroying shuffle

Sites whose ProteinMPNN parse and BioPython geometry do not align are dropped
and counted rather than scored on indices that cannot be trusted.

Usage:
    51_design_fixed_sequon.py [--proteins 40] [--designs 32] [--policy full_sequon]
"""
import argparse, sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")
sys.path.insert(0, "glyco_context/src")
from experimental_glycosylation_sites.retention import design_sequences
from experimental_glycosylation_sites.runner_support import (proteinmpnn_dir,
                                                             structure_paths)
from experimental_glycosylation_sites.mpnn_scoring import (chain_mapping, load_model,
                                                           to_manifest_space)
from experimental_glycosylation_sites.structures import _parse_chains
from glyco_context.fixed_design import sequon_positions, verify_sequon_index
from glyco_context.local_chemistry import (chemistry_panel, disulfide_indices,
                                           shell_indices_from_structure)

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--core", default="glyco_context/results/datasets/context_triplet_core.csv")
parser.add_argument("--indices", default="results/manifests/candidate_manifest_dataset.csv")
parser.add_argument("--proteins", type=int, default=40)
parser.add_argument("--designs", type=int, default=32)
parser.add_argument("--temperature", type=float, default=0.1)
parser.add_argument("--policy", default="full_sequon")
parser.add_argument("--fix-disulfides", action="store_true",
                    help="also hold cysteines in detected disulfides fixed, as the pre-specification states. Off by default so existing results reproduce; turning it on requires regenerating the designs.")
parser.add_argument("--seed", type=int, default=0)
parser.add_argument("--reuse-designs", default=None,
                    help="a sequences CSV from an earlier run. Designs are read from it rather than regenerated, so a control or panel fix costs no model time.")
parser.add_argument("--panels-only", action="store_true",
                    help="wild-type panels for every site, no design. This builds the natural reference distribution, which needs all the sites but none of the designs.")
parser.add_argument("--out", default="glyco_context/results/analysis/fixed_sequon_panels.csv")
args = parser.parse_args()

core = pd.read_csv(args.core, low_memory=False)
occ = core[core.population == "occupied"].copy()
indices = pd.read_csv(args.indices, low_memory=False)
index_columns = ["n_model_index", "plus1_model_index", "plus2_model_index"]
occ = occ.merge(indices[["accession", "position"] + index_columns].drop_duplicates(
    ["accession", "position"]), on=["accession", "position"], how="inner")
occ = occ.dropna(subset=index_columns)
print(f"occupied sites with validated model indices: {len(occ)}")

# Balanced across NXS/NXT, one site per chain group taken together so a chain is
# designed once and every site on it is measured from the same designs.
rng = np.random.default_rng(args.seed)
chains = occ.groupby(["structure_pdb_id", "structure_chain_id"])
groups = list(chains)
subtype_of = {key: group.subtype.iloc[0] for key, group in groups}
by_subtype = {"NXS": [], "NXT": []}
for key, _ in groups:
    by_subtype.setdefault(subtype_of[key], []).append(key)
selected = []
for subtype in ("NXS", "NXT"):
    pool = by_subtype.get(subtype, [])
    rng.shuffle(pool)
    selected.extend(pool[: args.proteins // 2])
if args.panels_only:
    selected = {key for key, _ in groups}
selected = set(selected)
print(f"selected {len(selected)} chains "
      f"({sum(1 for k in selected if subtype_of[k]=='NXS')} NXS, "
      f"{sum(1 for k in selected if subtype_of[k]=='NXT')} NXT)")

paths = structure_paths(())
stored = (pd.read_csv(args.reuse_designs, low_memory=False)
          if args.reuse_designs else None)
if stored is not None:
    print(f"reusing designs from {args.reuse_designs}")
model = (None if args.panels_only or stored is not None
         else load_model(proteinmpnn_dir(), device="cpu"))

rows, sequences, dropped, t0 = [], [], [], time.time()
protected_cys = set()
for n, (key, group) in enumerate([g for g in groups if g[0] in selected], start=1):
    pdb_id, chain_id = key
    path = paths.get(str(pdb_id).upper())
    if path is None:
        dropped.append({"pdb": pdb_id, "chain": chain_id, "reason": "structure_not_cached"})
        continue
    # Everything below is in the manifest's index space: the wild type comes
    # from the same parse the manifest counts, designs are projected back into
    # it, and the shell is built from the same residue list. One index space
    # throughout is the only reason this is checkable at all.
    try:
        mapping, _ = chain_mapping(path, str(chain_id), pdb_id)
        native = next((c for c in _parse_chains(path, str(pdb_id))
                       if c.chain_id == str(chain_id)), None)
        if not mapping or native is None:
            raise KeyError("chain cannot be mapped onto ProteinMPNN's parse")
        wild_type = native.sequence
    except Exception as exc:
        dropped.append({"pdb": pdb_id, "chain": chain_id,
                        "reason": f"{type(exc).__name__}: {str(exc)[:60]}"})
        continue

    fixed, site_shells = [], {}
    for site in group.itertuples():
        n_index = int(site.n_model_index)
        # Index check first, against the decoded sequence itself. If this fails
        # the shell would be centred on the wrong residue and every number for
        # this site would be quietly wrong.
        if not verify_sequon_index(wild_type, n_index, str(site.triplet_expected)):
            dropped.append({"pdb": pdb_id, "chain": chain_id,
                            "accession": site.accession, "position": site.position,
                            "expected": str(site.triplet_expected),
                            "observed": wild_type[n_index:n_index + 3]
                                        if n_index + 3 <= len(wild_type) else "",
                            "reason": "sequon_index_mismatch"})
            continue
        shell = shell_indices_from_structure(
            # n_resseq / n_icode are the QC record of the residue actually
            # measured, which is what the shell must be centred on.
            path, str(chain_id), int(site.n_resseq),
            getattr(site, "n_icode", ""), wild_type, n_index)
        if shell is None:
            dropped.append({"pdb": pdb_id, "chain": chain_id,
                            "accession": site.accession, "position": site.position,
                            "reason": "sequence_structure_misalignment"})
            continue
        site_shells[(site.accession, site.position)] = (n_index, shell, site)
        # The mask ProteinMPNN reads is in its own index space, so the sequon
        # positions are translated on the way in and the designs on the way out.
        fixed.extend(mapping[i] for i in sequon_positions(n_index, policy=args.policy)
                     if i in mapping)
    if not site_shells:
        continue

    if args.fix_disulfides:
        cys = disulfide_indices(path, str(chain_id), wild_type)
        fixed.extend(mapping[i] for i in cys if i in mapping)
        protected_cys.update(cys)

    if args.panels_only:
        designs = []
    elif stored is not None:
        prior = stored[(stored.structure_pdb_id == pdb_id)
                       & (stored.structure_chain_id == chain_id)
                       & (stored.variant == "design")].sort_values("replicate")
        designs = prior.sequence.tolist()
        if not designs:
            dropped.append({"pdb": pdb_id, "chain": chain_id,
                            "reason": "no stored designs for this chain"})
            continue
    else:
        designs = [to_manifest_space(d, mapping) for d in
                   design_sequences(path, str(chain_id), model, n_designs=args.designs,
                                    temperature=args.temperature, seed=args.seed,
                                    fixed_positions=sorted(set(fixed)))]

    # Random control: the same number of altered positions, drawn from this
    # chain's own composition, and never touching the fixed sequon positions.
    # `fixed` is in ProteinMPNN's index space, for the sampler's mask. The random
    # control edits the manifest-space wild type, so it needs the sequon
    # positions in *that* space -- using the model-space set here let the control
    # mutate the very sequon it is supposed to hold constant.
    protected = {i for (n_index, _, _) in site_shells.values()
                 for i in sequon_positions(n_index, policy=args.policy)}
    protected |= protected_cys
    designable = [i for i in range(len(wild_type)) if i not in protected]
    background = np.array(list(wild_type))
    randoms = []
    for design in designs:
        changed = sum(1 for i in designable if design[i] != wild_type[i])
        seq = list(wild_type)
        if changed and designable:
            spots = rng.choice(designable, size=min(changed, len(designable)),
                               replace=False)
            for spot in spots:
                # Draw until the residue actually differs. Sampling the chain's
                # composition without this leaves the original residue in the
                # draw, so a fraction of "mutations" are no-ops -- about
                # sum(freq^2), which is ~6.6% here -- and the control ends up
                # less perturbed than the design it is matched to.
                current = wild_type[spot]
                alternatives = background[background != current]
                if len(alternatives):
                    seq[spot] = str(rng.choice(alternatives))
        randoms.append("".join(seq))

    # Full sequences kept alongside the panels: the panel is a summary, and a
    # structure prediction later needs the sequence itself.
    sequences.append({"structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                      "variant": "wild_type", "replicate": 0, "sequence": wild_type})
    for replicate, (design, control) in enumerate(zip(designs, randoms)):
        sequences.append({"structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                          "variant": "design", "replicate": replicate, "sequence": design})
        sequences.append({"structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                          "variant": "random", "replicate": replicate, "sequence": control})

    for (accession, position), (n_index, shell, site) in site_shells.items():
        common = {"accession": accession, "position": position,
                  "structure_pdb_id": pdb_id, "structure_chain_id": chain_id,
                  "subtype": site.subtype, "n_rsa": site.n_rsa,
                  "n_ss_coarse": site.n_ss_coarse,
                  "shell_size": len(shell)}
        rows.append({**common, "variant": "wild_type", "replicate": 0,
                     "n_mutations": 0,
                     **chemistry_panel(wild_type, n_index, shell)})
        for replicate, (design, control) in enumerate(zip(designs, randoms)):
            mutations = sum(1 for i in designable if design[i] != wild_type[i])
            rows.append({**common, "variant": "design", "replicate": replicate,
                         "n_mutations": mutations,
                         **chemistry_panel(design, n_index, shell)})
            rows.append({**common, "variant": "random", "replicate": replicate,
                         "n_mutations": mutations,
                         **chemistry_panel(control, n_index, shell)})
    if n % 5 == 0:
        print(f"  {n}/{len(selected)} chains, {len(rows)} rows, "
              f"{(time.time()-t0)/60:.1f} min", flush=True)

frame = pd.DataFrame(rows)
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
frame.to_csv(args.out, index=False)
print(f"\n{len(frame)} rows over {frame.groupby(['accession','position']).ngroups} sites "
      f"-> {Path(args.out).resolve()}")
if sequences:
    seq_path = Path(args.out).with_name("fixed_sequon_sequences.csv")
    pd.DataFrame(sequences).drop_duplicates(
        ["structure_pdb_id", "structure_chain_id", "variant", "replicate"]).to_csv(
        seq_path, index=False)
    print(f"sequences -> {seq_path.resolve()}")

if dropped:
    pd.DataFrame(dropped).to_csv(Path(args.out).with_name("fixed_sequon_dropped.csv"), index=False)
    print(f"dropped: {len(dropped)}")
    print(pd.DataFrame(dropped).reason.value_counts().to_string())
