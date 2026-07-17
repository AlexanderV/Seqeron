---
type: concept
title: "Splice acceptor site prediction (3' splice site: AG, polypyrimidine tract, branch point)"
tags: [splicing, motif, algorithm]
mcp_tools:
  - find_branch_points
sources:
  - docs/algorithms/Splicing/Acceptor_Site_Detection.md
  - docs/Evidence/SPLICE-ACCEPTOR-001-Evidence.md
  - docs/Evidence/SPLICE-DONOR-001-Evidence.md
  - docs/Evidence/SPLICE-PREDICT-001-Evidence.md
source_commit: 3ae2cfc56ff0d4e03325af1ce7f3fcfd1313ab69
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: splice-acceptor-001-evidence
      evidence: "Test Unit ID: SPLICE-ACCEPTOR-001 — Algorithm: Acceptor Site Detection"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:known-motif-search
      source: splice-acceptor-001-evidence
      evidence: "AcceptorPwm position-weight-matrix scoring of the 3' splice site — the position-weight-matrix branch of the degenerate/consensus motif family that known-motif-search names as its non-exact counterpart"
      confidence: medium
      status: current
---

# Splice acceptor site prediction (3' splice site)

The **acceptor site** (**3' splice site**, 3'ss) is the intron end that the
[spliceosome](https://en.wikipedia.org/wiki/Spliceosome) cleaves and ligates to the
downstream exon. In the canonical **GU-AG rule** (>99% of introns) the intron begins
`GU` at the 5' **donor** site and ends `AG` at the 3' acceptor. Acceptor-site
prediction scores candidate `AG` dinucleotides for the surrounding signal that marks a
*real* 3'ss versus a **cryptic** intronic `AG`. This concept is the anchor for the
**splicing family** (acceptor / donor / branch point); it is the first ingested member,
and its 5' partner is [[splice-donor-site-prediction]] (SPLICE-DONOR-001 — the `GU`/GT
donor).
Validated under test unit **SPLICE-ACCEPTOR-001** — record
[[splice-acceptor-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern. Seqeron exposes
`FindAcceptorSites` (default) plus opt-in `FindAcceptorBranchPoint` and
`ScoreAcceptorMaxEnt`.

## The three 3'ss signals

A functional acceptor is a composite of three cis-elements, ordered 5'→3' inside the
intron and ending at the exon boundary:

1. **Branch point sequence (BPS)** — an adenosine (the lariat branch nucleotide)
   ~18–50 nt upstream of the AG. The human consensus is **`yUnAy`** (Gao et al. 2008):
   positions −3..+1 with the branch **A at position 0**, conserved as
   y@−3 ≈ 79%, U@−2 ≈ 74.6%, **A@0 ≈ 92.3%**, y@+1 ≈ 75.1% (n = 181 lariat-confirmed
   sites; median location −26 nt). Corroborated genome-wide by Mercer et al. (2015),
   59,359 branchpoints, ~19–35 nt from the 3'ss.
2. **Polypyrimidine tract (PPT)** — a C/U-rich stretch (usually 15–20 nt, ~5–40 nt
   before the intron 3' end) between the branch point and the AG. Recruits **U2AF65**;
   a **strong** (continuous) PPT gives robust splicing, a **weak** (purine-interrupted)
   PPT scores lower and can be skipped. Per Gao (2008) the PPT spans 4–24 nt downstream
   of the branch point.
3. **The AG dinucleotide + context** — consensus `(Yn)NYAG|G` (Burge et al. 1999).
   Shapiro & Senapathy (1987), from 3,700+ introns, put the acceptor AG at
   near-100% conservation (position −2 A, −1 G), a strong C at −3 (~65–70%), and
   G as the most frequent first exonic base (+1, ~50%).

## Three scoring surfaces

| Method | Model | Source |
|--------|-------|--------|
| `FindAcceptorSites` (default) | **AcceptorPwm** position-weight matrix + PPT-quality term, normalized to [0,1] via `(score/(count+1) + 2)/4` | Shapiro & Senapathy 1987 (PWM weights); normalization is a design-choice heuristic |
| `FindAcceptorBranchPoint` (opt-in) | scans the **18–40 nt** window upstream of the AG for `yUnAy`; conservation-weighted score `matched/maxScore`, `maxScore = 0.790+0.746+0.923+0.751 = 3.210` (perfect `yUnAy` = 1.0); reports the best candidate iff score ≥ `minScore` (**default 0.5**) | Gao et al. 2008 + Mercer et al. 2015 |
| `ScoreAcceptorMaxEnt` (opt-in) | **MaxEntScan score3ss**: max-entropy model over a 23-nt window (20 intronic + 3 exonic; AG fixed at 0-based 18–19), captures position dependencies beyond a PWM; score = `log2(P_maxent/P_bgd)` | Yeo & Burge 2004 |

The PWM+PPT default score does **not** include a branch-point term; branch-point
detection is a separate opt-in call. All three are additive — the default scorer is
unchanged by the opt-ins.

### MaxEntScan factorisation + provenance

`ScoreAcceptorMaxEnt` reproduces the maxentpy `score3` factorisation exactly: the
conserved AG is scored separately (`cons1_3[A]·cons2_3[G]/(bgd_3[A]·bgd_3[G])`) and the
remaining 21-nt "rest" is factorised over nine overlapping sub-sequences (five
multiplied, four divided). The precomputed probability tables (`score3_matrix.txt`,
**82,560 records**) were embedded as `Data/maxent_score3.txt` from the **MIT-licensed
maxentpy port** (kepbod/maxentpy). The canonical cross-check is
`score3('ttccaaacgaacttttgtAGgga')` → **2.89** (full precision 2.886773); two further
reference values (**8.19**, **−0.08**) reproduce exactly, so a wrong table or
factorisation would fail the 2.89 check. **Licence flag:** the bundled table is the
MIT port, not the original Burge-lab Perl models (which carry academic terms) —
recorded in `Data/maxent_score3.LICENSE.md`.

## Method contract (algorithm spec)

Per `docs/algorithms/Splicing/Acceptor_Site_Detection.md` (SPLICE-ACCEPTOR-001), the
`SpliceSitePredictor` public surface is:

| Signature | Role |
|-----------|------|
| `FindAcceptorSites(string sequence, double minScore = 0.5, bool includeNonCanonical = false)` | default scan; emits `SpliceSite` records |
| `FindAcceptorBranchPoint(string sequence, int agEnd, double minScore = 0.5)` | opt-in branch-point detection → `BranchPointResult` |
| `ScoreAcceptorMaxEnt(string window)` | opt-in Yeo & Burge 2004 MaxEntScan `score3ss`, in bits |

Private helpers: `ScoreAcceptorSite`, `ScoreU12AcceptorSite`, `ScoreBranchPointConsensus`.
The scan uppercases the sequence and maps `T`→`U`; canonical `AG` scanning always runs,
begins at index **15** (so the upstream PPT window exists), and `includeNonCanonical=true`
additionally scans U12 `AC`. A site is emitted only when its normalized score ≥ `minScore`.

**`SpliceSite` output record:** `Position` (for `AG`, the index `i+1` of the `G` — the
terminal *intronic* base, **not** the first exonic base; for `AC`, the index of the `C`),
`Type` (`Acceptor` for canonical `AG`, `U12Acceptor` for `AC`), `Motif`
(`GetMotifContext(sequence, position-1, 15, 2)` local context), `Score` (normalized), and
`Confidence` (clamped linear interpolation of the score).

**Canonical scoring internals** (`ScoreAcceptorSite`): (1) a **PPT term** — count
pyrimidines in the window `[position-15, position-3)`, divide by `12`, scale by `2`; plus
(2) a **sparse `AcceptorPwm`** applied at offsets `-15, -10, -5, -4, -3, -2, -1, 0` (indexed
as `position + 2 + offset`). The combined value is normalized by `(score/(count+1) + 2)/4`
and clamped to `[0,1]`. **U12 scoring** (`ScoreU12AcceptorSite`): `C` at `-1` and `-2`, a
pyrimidine at `-3`, plus an upstream PPT fraction, normalized by the maximum score **3.5**.

**Invariants (INV-01…INV-07):** `Score` and `Confidence` are in `[0,1]` (both scorers clamp);
canonical `AG` → `Type = Acceptor`, `AC` → `Type = U12Acceptor`; canonical `Position = i+1`
(the `G`); branch-point `Score ∈ [0,1]` with a perfect `yUnAy` = 1.0, reported only when the
branch A lies **18–40 nt inclusive** upstream of the `AG`; and branch-point detection is
**additive** — `FindAcceptorSites` output is unchanged whether or not it is called (the
default acceptor score carries no branch-point term).

## Non-canonical acceptors and corner cases

- **U12-type (minor spliceosome):** ~0.4% of human introns; the 3'ss uses **AC**
  instead of AG. Seqeron's opt-in `includeNonCanonical=true` scores these against the
  **YCCAC** consensus (Hall & Padgett 1994; Jackson 1991) rather than a fixed constant.
- **GC-AG introns:** ~1% of mammalian introns use GC at the 5' donor but keep AG at the
  acceptor (Burge et al. 1999) — the acceptor logic is unaffected.
- **Cryptic acceptors:** intronic `AG` dinucleotides are decoys; PPT quality + AG
  position distinguish real from cryptic.
- **Guards:** sequences **< 20 nt** return empty (insufficient PPT context); no AG in
  the scannable range returns empty; scan starts after position ~15. DNA `T`↔RNA `U`
  equivalence and case-insensitivity are honoured.

## Relation to other units

The 5' partner is [[splice-donor-site-prediction]] — the `GU`/GT donor at the *other* end
of the same intron under the GU-AG rule. The composite
[[gene-structure-prediction-intron-exon]] (SPLICE-PREDICT-001) **consumes this acceptor unit**,
pairing each `AG` 3'ss with a `GT/GU` donor into an intron and assembling the full
exon-intron gene structure (with exon typing, phase, and the spliced sequence). The donor is a single contiguous `MAG|GURAGU`
motif (no branch point / PPT), so its default scorer is a plainer **consensus match
fraction** and its MaxEntScan **score5ss** uses a single 16,384-entry table (vs this
acceptor's PWM+PPT default and score3ss's 82,560-record overlapping-sub-window model).

The default acceptor scorer is a **position-weight-matrix** scan — the PWM branch of
the degenerate/consensus motif family whose exact-match baseline is
[[known-motif-search]] and whose fixed-catalog DNA analogue is
[[regulatory-element-detection]]. The `yUnAy` and `AcceptorPwm` consensus strings are
IUPAC/consensus motifs in the sense of [[consensus-sequence]]. The MaxEntScan model is
a probabilistic sequence-motif scorer, distinct from the exact and PWM scanners.

## Reference sources

Wikipedia **RNA splicing / Polypyrimidine tract / Spliceosome** (rank 4);
**Shapiro & Senapathy 1987** (3'ss nucleotide statistics, PWM weights);
**Burge, Tuschl & Sharp 1999** (`(Yn)NYAG|G` structure, GC-AG / U12);
**Yeo & Burge 2004** (MaxEntScan score3ss) with the MIT-licensed maxentpy tables;
**Gao et al. 2008** (`yUnAy` human branch-point consensus, per-position conservation);
**Mercer et al. 2015** (genome-wide branchpoint distribution);
**Hall & Padgett 1994 / Jackson 1991** (U12 YCCAC). No source contradictions — the
encyclopedic, statistical, and max-entropy sources are mutually consistent; the only
non-source elements are the [0,1] normalization heuristic and API-shape guards.
