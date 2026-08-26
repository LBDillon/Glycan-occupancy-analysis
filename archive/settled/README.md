# Archived — settled questions

Each of these answered a question that is now closed. Kept because the answer
rests on them, removed from the live pipeline because rerunning them would
change nothing.

- **`12_matching_sensitivity.py`** — the 200-seed sweep that showed greedy
  matching's significance was being decided by an arbitrary seed while its
  direction was stable. That is why matching is deterministic and no seed is
  involved. The question does not reopen.

- **`14_convergence_check.py`** and **`42_marginal_k_check.py`** — how many
  samples ESM-IF's marginalised joint score needs before it stops moving. Fixed
  a free parameter that would otherwise have been chosen arbitrarily. Settled.

- **`30_package_for_colab.py`** — bundled structures for Colab, which was
  abandoned for ARC: sessions died at twelve hours and a dropped runtime lost
  the work. Job arrays solved it. The route is not in use.
