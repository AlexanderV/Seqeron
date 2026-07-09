---
type: source
title: "Evidence: COMPGEN-ANI-001 (Average Nucleotide Identity — ANIb)"
tags: [validation, comparative-genomics]
doc_path: docs/Evidence/COMPGEN-ANI-001-Evidence.md
sources:
  - docs/Evidence/COMPGEN-ANI-001-Evidence.md
source_commit: 908b2139cdb38b9897cafa7d07a0d525758cb281
ingested: 2026-07-09
created: 2026-07-09
updated: 2026-07-09
---

# Evidence: COMPGEN-ANI-001

The validation-evidence artifact for test unit **COMPGEN-ANI-001** — Average Nucleotide
Identity under the **ANIb** definition (Goris et al. 2007): fragment the query genome, place
each fragment against the reference, filter by identity/coverage cut-offs, and average the
per-fragment identities (`CalculateAni` / `CalculateReciprocalAni`). This is the **first
Comparative-genomics** family Evidence file and one instance of the templated per-algorithm
[[algorithm-validation-evidence|evidence artifact]] pattern; the algorithm, its parameters,
invariants, worked oracles, and corner cases are summarized in
[[average-nucleotide-identity]]. Its sibling COMPGEN unit is the shared synteny anchor
[[synteny-and-rearrangement-detection]]. See [[test-unit-registry]] for how units are tracked.

## What this file records

- **Online sources:**
  - **Goris et al. 2007** (IJSEM 57:81–91, DOI 10.1099/ijs.0.64483-0, authority rank 1) — the
    ANIb method: query "cut into consecutive 1020 nt fragments" (mirroring ~1 kb DDH
    fragmentation), each searched by BLASTN against the reference, best match kept; ANI =
    "mean identity of all BLASTN matches that showed more than 30 % overall sequence identity
    (recalculated to an identity along the entire sequence) over an alignable region of at
    least 70 % of their length". Species boundary: ANI ≈ 95 % ↔ 70 % DDH. Reverse searching
    performed "to provide reciprocal values".
  - **Konstantinidis & Tiedje 2005** (PNAS 102(7):2567–2572, rank 1) — ANI of shared genes as
    a robust relatedness measure; ANI ≈ 94 % ↔ 70 % DNA–DNA reassociation.
  - **pyani** (Pritchard et al. 2016, ANIb reference implementation, rank 3) — confirms 1020 nt
    fragments, the 70 % coverage / 30 % identity cut-offs, and the mean-of-qualifying-fragments
    average. The 2026-06-23 refresh fetched the pyani `anib.html` source: gapped BLASTN
    (`-xdrop_gap_final 150`, `ani_alnlen = blast_alnlen - blast_gaps`), identity/coverage
    recalculated over query length (`ani_pid = ani_alnids/qlen`, `ani_coverage = ani_alnlen/qlen`),
    and paired jobs in both directions.
- **Algorithm behaviour (from the artifact):** fragment → best-match placement → per-fragment
  identity/coverage recalculated over the **query-fragment length** → discard fragments below
  30 % identity or 70 % alignable → mean over survivors. Gapped Smith-Waterman placement
  recovers indels; reciprocal ANI = mean of both directions (symmetric).
- **Datasets (documented oracles):**
  - *Synthetic exact-arithmetic* (ref `AAAACCCCGGGGTTTT`, fragLen 4): identical → 1.0; one
    mismatch in last fragment → 0.9375; `AATT` last fragment → 0.875; identity cut-off
    (`CGTC` = 0 excluded) → 1.0; alignable cut-off (ref shorter than fragment) → 0; query
    shorter than fragment → 0.
  - *Gapped placement*: `AAAACCCC` vs `AAAATCCCC` → ungapped 0.875 < gapped 1.0.
  - *Reciprocal*: identical genomes → 1.0; symmetry ANI(A,B)=ANI(B,A); `AAAACGTC`/`AAAAAAAA`
    fragLen 4 → (1.0 + 1.0)/2 = 1.0.

## Deviations and assumptions

The artifact records **no unresolved correctness-affecting assumptions** — the earlier
"ungapped placement" assumption is **resolved** by the 2026-06-23 fix (gapped Smith-Waterman
placement + `CalculateReciprocalAni` reciprocal mean, per Goris/pyani). The single remaining
item is an explicit **engine decision, not a correctness gap**: the gapped path uses the
library's own Smith-Waterman aligner (`SequenceAligner.LocalAlign`, BLAST DNA scoring) rather
than the NCBI BLASTN engine — full DP, more sensitive than BLAST's heuristic seeding, with
identity/coverage on the same recalculated-over-fragment definition. Numeric ANI on real
genomes may differ slightly from NCBI-BLASTN pipelines because the alignment engine differs;
indel handling is correct (documented in algorithm doc §5.3). No contradictions among sources
— Goris, Konstantinidis & Tiedje, and pyani agree on the fragmentation, cut-offs, averaging,
gapped placement, and reciprocal computation.
