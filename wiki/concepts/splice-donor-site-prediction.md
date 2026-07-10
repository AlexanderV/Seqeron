---
type: concept
title: "Splice donor site prediction (5' splice site: GU/GT dinucleotide, MAG|GURAGU consensus, MaxEntScan score5ss)"
tags: [splicing, motif, algorithm]
sources:
  - docs/Evidence/SPLICE-DONOR-001-Evidence.md
source_commit: ce6f817f61151956d1e97909c1ccf5d70f0c333c
created: 2026-07-10
updated: 2026-07-10
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: splice-donor-001-evidence
      evidence: "Test Unit ID: SPLICE-DONOR-001 — Algorithm: Donor (5') Splice Site Detection"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:splice-acceptor-site-prediction
      source: splice-donor-001-evidence
      evidence: "The 5' donor (GU) and 3' acceptor (AG) are the two ends of the same intron under the canonical GU-AG rule; this donor unit is the sibling the acceptor anchor page anticipated"
      confidence: high
      status: current
---

# Splice donor site prediction (5' splice site)

The **donor site** (**5' splice site**, 5'ss) is the intron *start* — the
[spliceosome](https://en.wikipedia.org/wiki/Spliceosome) cleaves the exon-intron
boundary here and the intron begins with an almost invariant **GU** (GT in DNA)
dinucleotide. Under the canonical **GU-AG rule** (>99% of introns) the intron opens `GU`
at this donor and ends `AG` at the 3' [[splice-acceptor-site-prediction|acceptor]]; the
two are the paired ends of one intron. Donor-site prediction scores candidate `GT/GU`
dinucleotides for the surrounding signal that marks a *real* 5'ss versus a **cryptic**
`GT`. This is the **donor member of the splicing family** whose anchor is the
[[splice-acceptor-site-prediction]] page. Validated under test unit **SPLICE-DONOR-001** —
record [[splice-donor-001-evidence]]; [[test-unit-registry]] tracks the unit and
[[algorithm-validation-evidence]] describes the artifact pattern. Seqeron exposes
`FindDonorSites` / `ScoreDonorSite` (default) plus opt-in `ScoreDonorMaxEnt`.

## The 5'ss consensus

Unlike the 3' acceptor (a composite of branch point + polypyrimidine tract + AG), the
donor is a **single contiguous 9-nt motif** spanning the boundary, positions −3..+6
(the `|` is the exon-intron junction):

**`MAG | GU RAGU`** (M = A/C, R = A/G) — Alberts *Molecular Biology of the Cell* /
Shapiro & Senapathy 1987 `(C/A)AG|GU(A/G)AGU`. Per-position conservation from
Shapiro & Senapathy (thousands of sites):

| Position | −3 | −2 | −1 | 0 | +1 | +2 | +3 | +4 | +5 |
|----------|----|----|----|----|----|----|----|----|----|
| Consensus | M (A/C) | A | **G** | **G** | **U** | R (A/G) | A | G | U |
| Conservation | ~35% each | ~60% | ~80% | ~100% | ~100% | | | | |

Positions **0 (G)** and **+1 (U)** are near-invariant; **−1 (G)** is strongly conserved
(~80%) and critical for **U1 snRNP** base-pairing (U1 snRNP binds the GU during E-complex
formation). Because branch point and polypyrimidine tract live on the 3' side, the donor
page has **no branch-point surface** — the structural asymmetry with the acceptor.

## Two scoring surfaces

| Method | Model | Source |
|--------|-------|--------|
| `FindDonorSites` / `ScoreDonorSite` (default) | **IUPAC consensus match fraction** — score = matches / positions scored against the `MAG|GURAGU` binary weights (1.0 = match, 0.0 = mismatch), normalized to [0,1]; **no PWM, no ad-hoc formula** | `MAG|GURAGU` consensus (Wikipedia RNA splicing / Alberts; Shapiro & Senapathy 1987) |
| `ScoreDonorMaxEnt` (opt-in) | **MaxEntScan score5ss**: max-entropy model over a **9-nt** window (3 exon + 6 intron); score = `log2(P_maxent/P_bgd)` | Yeo & Burge 2004 |

Both are additive — the default scorer is unchanged by the opt-in. This is a deliberately
**simpler default** than the acceptor's `FindAcceptorSites` (which is a Shapiro/Senapathy
**PWM + PPT-quality** scorer with the `(score/(count+1)+2)/4` normalization): the donor
default is a plain consensus-match fraction, so **GC donors score lower for free** —
position +1 (C) mismatches the invariant `U`, capping them at 8/9 ≈ 0.889 vs 9/9 = 1.0 for
GT, with no explicit penalty constant.

### MaxEntScan score5 factorisation + provenance

`ScoreDonorMaxEnt` reproduces the maxentpy `score5` factorisation exactly. The conserved
`GT` at 0-based positions **3..4** is scored separately and removed; the remaining 7-nt
"**rest**" sequence is `window[0:3] + window[5:9]`, whose max-entropy probability is looked
up **directly in a single 4⁷ = 16,384-entry table** keyed by the 7-mer string. This is the
key contrast with **score3ss** (acceptor): score5 is **single-matrix with no overlapping
sub-windows**, where score3 factorises a 21-nt rest over nine overlapping sub-sequences
(82,560 records). Final score = `log2(GT_term · rest_term)` where
`GT_term = cons1_5[G]·cons2_5[T] / (bgd_5[G]·bgd_5[T])` and `rest_term = matrix[rest]`
(`bgd_5 = {A:.27,C:.23,G:.23,T:.27}`; `cons1_5[G]=.9896`; `cons2_5[T]=.9884`). The
precomputed table is embedded as `Data/maxent_score5.txt` (16,384 records) from the
**MIT-licensed maxentpy port** (kepbod/maxentpy). Canonical cross-checks:
`score5('cagGTAAGT')` → **10.86** (10.858313, the primary check),
`score5('gagGTAAGT')` → **11.08**, `score5('taaATAAGT')` → **−0.12** (a non-GT donor).
**Licence flag:** the bundled table is the MIT port, not the original Burge-lab Perl models
(academic terms) — recorded in `Data/maxent_score5.LICENSE.md`.

## Non-canonical donors and corner cases

- **GC-AG introns:** ~0.5–1% of U2-type introns use **GC** instead of GT at the 5'ss,
  still AG at the acceptor (Burge et al. 1999) — valid, processed by the major
  spliceosome, and naturally lower-scoring (see above).
- **U12-type (minor spliceosome):** ~0.3% of introns use **AT-AC** boundaries — an **AT**
  donor with the extended `/ATATCC/` motif, scored against that consensus when
  `includeNonCanonical=true`.
- **Cryptic donors:** intronic/exonic `GT` decoys resemble the consensus but are not used;
  point mutations can activate them.
- **Guards:** sequences **shorter than the 9-nt window** cannot be scored → empty; no
  `GT/GU` in the sequence → empty; empty/null → empty (never throws). DNA `T`↔RNA `U`
  equivalence and case-insensitivity (`ToUpperInvariant`) are honoured.

## Relation to other units

The donor and [[splice-acceptor-site-prediction|acceptor]] are the two halves of the same
intron under the GU-AG rule; together with a future branch-point-as-its-own-unit they form
the splicing family anchored on the acceptor page. The default donor scorer is a
**consensus / IUPAC-degenerate** motif matcher (`MAG|GURAGU` is a
[[consensus-sequence|consensus motif]] / [[iupac-degenerate-consensus|IUPAC-degenerate
pattern]]), the position-independent cousin of the acceptor's position-weight-matrix
scanner and of the broader [[known-motif-search]] / [[regulatory-element-detection]] motif
family. The MaxEntScan model is the probabilistic sequence-motif scorer shared with the
acceptor's score3ss.

## Reference sources

Wikipedia **RNA splicing / Spliceosome** (rank 4, `MAG|GURAGU` consensus, GU-AG rule, U1
snRNP GU binding); **Shapiro & Senapathy 1987** (rank 1, 5'ss position frequencies + PWM
approach); **Burge, Tuschl & Sharp 1999** (rank 1, GC-AG / U12 AT-AC / `ATATCC`);
**Yeo & Burge 2004** (rank 1, MaxEntScan score5ss) with the MIT-licensed maxentpy 16,384-record
table; **Alberts et al.** *Molecular Biology of the Cell* (donor consensus). No source
contradictions — the encyclopedic, statistical, and max-entropy sources are mutually
consistent; the only non-source elements are the [0,1] match-fraction normalization
(a design choice) and API-shape guards. All prior assumptions are **RESOLVED**: PWM values
replaced by IUPAC binary consensus weights, score normalization by a plain match fraction,
and the old GC-donor 0.7 penalty removed (GC donors self-penalize via the +1 mismatch).
