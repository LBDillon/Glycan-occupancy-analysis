# Methods — locating a sequon in a model's own view of a chain

*How a site in the manifest is turned into an index a structure model can be
read at, why two indices are involved rather than one, and the checks that
confirm the mapping. Written 2026-08-24.*

Every model-facing measurement in this repository — the conditional sequon
score, sequon retention under redesign — is read at a numeric position in a
sequence the model produced. Getting that position right is the whole
measurement: an index that is off by one reads a different residue and returns a
perfectly plausible number for the wrong thing.

Two parsers are involved, and they do not index a chain the same way.

## The two parses

**This repository's parse** (`structures._parse_chains`) lists the residues that
were actually observed, in order. The manifest's `model_index` is an ordinal
into that list:

```python
"model_index": index - 1,          # 0-based, as a model reads the chain
```

**ProteinMPNN's parse** (`protein_mpnn_utils.parse_PDB_biounits`) does something
different. It walks the author residue numbering from lowest to highest and
emits one slot per number, inserting a placeholder where a number is absent:

```python
for resn in range(min_resn,max_resn+1):
    if resn in seq:
        for k in sorted(seq[resn]): seq_.append(aa_3_N.get(seq[resn][k],20))
    else: seq_.append(20)
```

Token 20 is `X`. So a chain numbered 24–445 with twenty numbers missing produces
a sequence of 422 positions, not 402. Insertion codes are expanded in the inner
loop, so 36, 36A and 36B occupy three slots in both parses.

The two indices coincide **only when a chain is numbered without gaps**. Most
depositions have gaps.

## How each model is read

**ESM-IF** resolves this by alignment, comparing the two sequences and remapping
where they differ (`esmif_scoring.py`):

```python
if native.sequence == esm_seq:
    # The common case: the two parsers agree residue for residue.
    to_esm = {i: i for i in range(len(esm_seq))}
else:
    to_esm = {a - 1: b - 1 for a, b in _alignment_pairs(native.sequence, esm_seq)}
```

**ProteinMPNN** is read through `mpnn_scoring.build_index_map`, which
reconstructs the enumeration above rather than aligning. The difference between
the parses is not a sequence difference — the residues are the same — so it is
arithmetic on residue numbers, and arithmetic is exact where an aligner has to
choose:

```python
mapping, model_index, native_index = {}, 0, 0
for number in range(min(by_number), max(by_number) + 1):
    if number in by_number:
        for _ in sorted(by_number[number]):
            mapping[native_index] = model_index
            native_index += 1
            model_index += 1
    else:
        model_index += 1                     # ProteinMPNN's placeholder slot
```

The reconstruction is a hypothesis about what the other parser did, so it is
checked against that parser's own output before being used. Any disagreement
returns an empty mapping and the chain is dropped:

```python
if native_index != len(native_sequence):
    return {}
for source, target in mapping.items():
    if target >= len(model_sequence) or model_sequence[target] != native_sequence[source]:
        return {}
```

## Worked examples

Three chains, chosen to cover a wide gap, a narrow one and none at all. In each
case the site is a known N-X-S/T sequon and the manifest records the triplet it
should be.

### 9G3Q chain A — A0A1S4F0I0 position 225, expected `NES`

```
biopython observed residues : 402      resseq 24..445, span 422
numbering gaps              : 20 missing residue numbers
proteinmpnn sequence length : 422      = 402 observed + 20 placeholders
manifest n_model_index      : 181      (an ordinal into the 402)
raw   seq[181:184]          : 'LKN'    read without mapping
mapped index                : 201      shift +20
mapped seq[201:204]         : 'NES'    matches the manifest triplet
```

The shift is exactly the number of absent residue numbers preceding the site.

### 5H5Y chain A — A0A023YYV9 position 257, expected `NRS`

```
biopython observed residues : 286      resseq 21..316, span 296
numbering gaps              : 10
proteinmpnn sequence length : 296
manifest n_model_index      : 226
raw   seq[226:229]          : 'VNI'
mapped index                : 236      shift +10
mapped seq[236:239]         : 'NRS'    matches
```

Again the shift equals the gap count: ten missing numbers, ten positions.

### 4EBY chain A — A8R7E6 position 52, expected `NSS`

```
biopython observed residues : 200      resseq 25..224, span 200
numbering gaps              : 0
proteinmpnn sequence length : 200
manifest n_model_index      : 27
raw   seq[27:30]            : 'NSS'
mapped index                : 27       shift 0
mapped seq[27:30]           : 'NSS'    matches
```

A contiguously numbered chain is a fixed point of the mapping, which is the
behaviour the identity case requires.

### A chain the mapping declines — 5O2W chain A

```
biopython observed residues : 248
proteinmpnn sequence length : 247      shorter, not longer
```

Gap-filling can only lengthen the sequence, so a shorter one means ProteinMPNN
excluded a residue this parse kept — a difference the reconstruction does not
model. The verification fails and the chain is dropped rather than read at an
index that cannot be justified.

## Independent confirmation

The reconstruction is checked against the sequences, which could in principle
agree while the model still decodes something else. The stronger check runs the
model: positions are held fixed during a real ProteinMPNN design, and the
residues that come back fixed are compared with the manifest.

```
5H5Y:A 257   manifest 'NRS'   design at mapped index 'NRS'
9G3Q:A 225   manifest 'NES'   design at mapped index 'NES'
4EBY:A  52   manifest 'NSS'   design at mapped index 'NSS'
4B8V:A  89   manifest 'NCS'   design at mapped index 'NCS'
```

Whatever the model writes at a fixed position is the native residue at that
position, so this confirms the index end to end rather than through a
reconstruction of the parse.

## Checks that run in the pipeline

- `build_index_map` returns `{}` unless every mapped position matches the model
  parse residue for residue.
- `verify_sequon_index` confirms the decoded sequence really reads N-X-S/T with
  X not proline at the claimed index, independently of any structure parse.
- `05_scoreability` records a site as unscoreable when its indices fall outside
  the model's decodable positions.

A site failing any of these is dropped and counted, never read at a best guess.
