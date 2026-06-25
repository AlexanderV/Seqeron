# Validation Report: SPLICE-DONOR-001 — Donor (5') Splice Site Detection

- **Validated:** 2026-06-25 (fresh re-validation; supersedes 2026-06-24)   **Area:** Splicing
- **Canonical method(s) (this unit's OWN surface):** `SpliceSitePredictor.FindDonorSites(sequence, minScore, includeNonCanonical)`; internal `ScoreDonorSite(seq, position)`, `ScoreU12DonorSite(seq, position)`. The IUPAC/consensus donor PWM (`DonorPwm`) + GU consensus detection.
- **Out of scope (separate CLEAN unit):** `ScoreDonorMaxEnt` (MaxEntScan score5ss) = **SPLICE-MAXENT5-001** — referenced, not re-litigated here.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

> **Re-validation note (2026-06-25):** This previously-CLEAN unit was reset to ⬜ when the
> MaxEntScan `ScoreDonorMaxEnt` (score5ss) method was added to the same class; that method was
> split into the separate already-CLEAN unit **SPLICE-MAXENT5-001**. This report re-validates
> SPLICE-DONOR-001's OWN canonical surface (consensus/IUPAC-PWM donor scoring + GU detection)
> fresh against external sources retrieved this session. The canonical donor code itself was
> **untouched** by the MaxEnt addition. Two missing Stage-A edge-case tests (window/boundary,
> non-ACGT) were added this session.

---

## Stage A — Description

### Sources opened THIS session (2026-06-25) & what they confirm

- **PMC PMC275472** (Roca et al., *Intrinsic differences between authentic and cryptic 5′ splice sites*, NAR 2003), retrieved this session — states **verbatim**: "The consensus 5′ splice site sequence is **MAG|GURAGU** (M = A or C; R = purine), and spans from position **−3** (the third nucleotide from the 3′ end of the upstream exon) to **+6** (the sixth nucleotide in the intron)." Also: the S&S matrix "reflects the degree of conservation … resulting from the **alignment of 1446 5′ splice sites**."
- **Wikipedia, "RNA splicing"**, retrieved this session — "The splice donor site includes an **almost invariant sequence GU at the 5′ end of the intron**, within a larger, less highly conserved region." Donor consensus notation `G‑G‑[cut]‑G‑U‑R‑A‑G‑U`. Confirms the **GT‑AG rule** and that the donor dinucleotide is **GU (GT in DNA)**, NOT AG.
- **PMC PMC1934990** (Vorechovsky, *Aberrant 5′ splice sites in human disease genes*, NAR 2007), retrieved this session — independently reproduces "The human 5′ss consensus sequence is **MAG|GURAGU** (M is A or C; R is purine) … spanning from position −3 to +6," and that **+1 = G, +2 = U/T are the invariant intron-start positions**.
- **Golden Helix, "Five Splice Site Algorithms"** + Springer BMC Bioinformatics (Reese/NNSplice family), retrieved this session — confirm S&S (1987) is a **position-weight-matrix** scorer derived from per-position nucleotide frequencies over a window around the splice site.
- **Established S&S per-position consensus** (widely reproduced; positions −3..+6): −3 **M(A/C)**, −2 **A**, −1 **G**, +1 **G** (~100%), +2 **U** (~100%), +3 **R(A/G)**, +4 **A**, +5 **G**, +6 **U** — exactly the IUPAC pattern `MAG|GURAGU`.

(References below kept from prior validation, consistent with the above.)
- **Shapiro & Senapathy (1987), NAR 15(17):7155**: donor consensus (C/A)AG|GU(A/G)AGU; position +1 (G) and +2 (U) ~100% conserved.
- **Burge, Tuschl & Sharp (1999), The RNA World**: GC-AG U2 introns (~0.5–1%) valid but weaker; U12-type AT-AC introns (~0.3%, donor AT/AU); U12 donor extended motif /ATATCC/.
- **Yeo & Burge (2004), MaxEntScan**: 9-mer donor window (−3..+6) — relevant to the sibling SPLICE-MAXENT5-001, not this unit.

### Formula / convention check

- Canonical 5' donor: intron BEGINS with **GU (GT in DNA)** — confirmed against PMC275472 / PMC1934990 / Wikipedia this session. Code's consensus `MAG|GURAGU` matches S&S exactly, invariant G,U at the first two intron positions.
- Coordinate convention: `Position` = index of the **G** of the GU dinucleotide (first intron base / exon-intron junction), 0-based. Standard, unambiguous junction convention.
- **Scoring model the code actually implements:** an **IUPAC-binary consensus** PWM (`DonorPwm`): per scored position, weight = 1.0 if the base is in the consensus set, else 0.0; the site score = (Σ matched weights)/(number of in-bounds positions scored) ∈ [0,1]. This is **NOT** the Shapiro–Senapathy log-frequency *consensus value* (CV%) formula — it is a faithful, conservative *reduction* of the same S&S consensus *motif* `MAG|GURAGU` to a binary match fraction. Stage A validates that (a) the IUPAC set at each of the 9 positions equals the published consensus nucleotide(s), and (b) the invariant GU is required (offsets 0/+1 = {G}/{U} singletons). Both hold. The richer S&S CV% / maximum-entropy scoring is provided separately by `ScoreDonorMaxEnt` (SPLICE-MAXENT5-001), out of scope here.

### Edge-case semantics

All sourced and defined: empty/null → empty; sequence < 6 nt → empty (insufficient window); no GU → empty; **GU at sequence start → partial window, score normalized over in-bounds positions only**; **non-ACGT base at a scored position → skipped (not counted in divisor), no throw**; GC valid-but-weaker (gated by `includeNonCanonical`); U12 AT/AU → separate scorer + `U12Donor` type; multiple GU each scored independently.

### Independent cross-check (hand computation, redone this session)

Independent re-implementation of the consensus `MAG|GURAGU` (offsets −3..+5 relative to the G of GU; binary match) reproduced the code EXACTLY for every case:

| Window | GU idx (Position) | Hand score | matches/scored |
|--------|-------------------|------------|----------------|
| `CAGGUAAGU` (perfect) | 3 | **1.0** | 9/9 |
| `UUUGUAAUU` (weak: −3 U✗, −2 U✗, −1 U✗, +4 U✗) | 3 | **0.5556** | 5/9 |
| `CAGGCAAGU` (GC; +1 C✗U) | 3 | **0.8889** | 8/9 |
| multi `CAGGUAAGU…CAGGUAAGU` | 3, 7, 26 | **1.0, 0.4444, 1.0** | 9/9, 4/9, 9/9 |
| `GUAAGUAAA` (GU at start; offsets −3..−1 OOB) | 0 | **1.0** | 6/6 |
| `NAGGUAAGU` (N at −3 skipped) | 3 | **1.0** | 8/8 |

Perfect-consensus hand trace `CAGGUAAGU`: −3→idx0 C∈{A,C}=1, −2→idx1 A=1, −1→idx2 G=1, **0→idx3 G=1, +1→idx4 U=1** (invariant GU), +2→idx5 A∈{A,G}=1, +3→idx6 A=1, +4→idx7 G=1, +5→idx8 U=1 → **9/9 = 1.0**. The GC variant scores 8/9 purely because +1 (C) fails the invariant U, confirming the recognised dinucleotide is **GU, not AG** (no donor/acceptor confusion).

### Findings / divergences

The IUPAC PWM window is **−3..+5 (9 positions)** vs the literature's −3..+6 (10 positions). It drops the weakly-conserved intron +6 (U) position; biologically immaterial (the invariant GU and the conserved −1 G / +3 A / +4 G are retained). No biological error. **PASS.**

---

## Stage B — Implementation

### Code path reviewed (re-read 2026-06-25)

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs`
- `DonorPwm` (130–141): offsets −3..+5, binary weights encoding M,A,G,**G,U**,R,A,G,U. Offset 0 = {G}, offset +1 = {U} singletons → invariant GU, not AG. Each position's IUPAC set matches the published `MAG|GURAGU` consensus exactly. Correct.
- `FindDonorSites` (204–259): `ToUpperInvariant` + T→U normalization (line 212), guard `IsNullOrEmpty || len < 6` (209) → empty, scans `i ≤ len−6` (214); canonical GU branch (217); `includeNonCanonical` adds GC (Donor, 231) and A[U/T] (U12Donor, 245) branches. Reports `Position = i` (G of GU).
- `ScoreDonorSite` (414–434): sums binary weights over in-bounds offsets only (bounds-guard 422), `weights.TryGetValue` skips non-A/C/G/U bases (425, divisor not incremented), divides by count scored → fraction in [0,1]; returns 0 if nothing scored.
- `ScoreU12DonorSite` (554–570): AUAUCC consensus match fraction /6.0 for the U12 (AT/AU) donor branch.

### Formula realised correctly?

Yes. Donor dinucleotide recognised = **GU/GT** at offsets 0/+1, matching the GT-AG rule and MAG|GURAGU consensus. GC scores 8/9 automatically (+1 mismatch). U12 AU routed to `/AUAUCC/` scorer, typed `U12Donor`. Position reported at the exon/intron junction (G of GU) — no off-by-one.

### Cross-verification table recomputed vs code (tests assert these exact values)

| Input | GU idx (Position) | Score | Note |
|-------|------|-------|------|
| CAGGUAAGU | 3 | 9/9 = 1.0 | perfect consensus |
| UUUGUAAUU | 3 | 5/9 ≈ 0.5556 | weak context |
| CAGGCAAGU (GC, nonCanonical) | 3 | 8/9 ≈ 0.8889 | +1 C ≠ U |
| CAGGUAAGU…CAGGUAAGU (multi) | {3, 7, 26} | {1.0, 4/9, 1.0} | full scan |
| GUAAGUAAA (GU at start) | 0 | 6/6 = 1.0 | partial window, OOB upstream offsets |
| NAGGUAAGU (N at −3) | 3 | 8/8 = 1.0 | non-ACGT skipped, not counted |

All recomputed by an independent hand re-implementation this session and confirmed by test assertions (last two are the newly added C3/C4 tests).

### Variant/delegate consistency

DNA (T) and RNA (U) inputs give identical results (T→U at line 212). Lowercase handled (`ToUpperInvariant`). `CalculateMaxEntScore` (legacy log-weight helper) reuses the same `DonorPwm` for its donor branch and stays consistent. The added `ScoreDonorMaxEnt` (SPLICE-MAXENT5-001) is purely additive and leaves all canonical donor paths and defaults unchanged.

### Test quality audit

`tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs` — the canonical-donor suite is **21 tests** (M1–M10, S1–S5, C1–C4, + helper); ME1–ME9 belong to SPLICE-MAXENT5-001. Assertions check **exact** consensus-traced values (scores 9/9, 8/9, 5/9, 4/9, 6/6, 8/8; positions 0/3/7/26; type Donor vs U12Donor), not tautologies or code echoes. The independent helper `ComputeExpectedConsensusStrength` re-derives consensus strength separately from the implementation.

**Gate finding (this session):** the prior suite was missing two Stage-A edge cases — **window/boundary (GU at sequence start)** and **non-ACGT input**. Added this session and locked to hand-computed values:
- **C3** `FindDonorSites_DonorAtSequenceStart_ScoresOverInBoundsPositionsOnly` — GU at index 0 → partial 6-position window → 6/6 = 1.0, Position 0.
- **C4a** `FindDonorSites_NonAcgtCharacterAtScoredPosition_IsSkippedNotCounted` — N at offset −3 skipped → 8/8 = 1.0.
- **C4b** `FindDonorSites_AllNonAcgt_ReturnsEmpty` — all-N input → no sites, no throw.

Edge cases now covered: canonical donor, degraded donor, no GU, multiple donors, window/boundary, non-ACGT, empty, null, short, GC gating, U12.

### Findings / defects

No code defect. Donor dinucleotide, consensus positions, GU invariance, junction coordinate, GC/U12 handling all match the validated description. The only Stage-B gap was missing edge-case tests (window/boundary, non-ACGT) — fixed this session; not a behavioural defect.

---

## Verdict & follow-ups

- **Stage A: PASS** — donor = GT/GU at intron start (GT-AG rule); consensus MAG|GURAGU (−3..+6, invariant +1 G / +2 U) confirmed against PMC275472 (Roca 2003), PMC1934990 (Vorechovsky 2007), Wikipedia (RNA splicing), all retrieved this session, plus Shapiro & Senapathy (1987) / Burge et al. (1999). The code implements a binary IUPAC reduction of this consensus (not the S&S CV% formula); the reduction is faithful at every position. Coordinate convention sound.
- **Stage B: PASS** — implementation recognises GU at the donor (offsets 0/+1 = {G}/{U}), scores against the IUPAC consensus as a match fraction, reports the junction position (G of GU) correctly; all worked examples (incl. boundary and non-ACGT) reproduce exactly against an independent hand re-implementation.
- **State: CLEAN** — no behavioural defect. The only gap was two missing Stage-A edge-case tests (window/boundary, non-ACGT), added and locked to hand-computed values this session.
- Note (non-blocking): IUPAC PWM uses −3..+5 (drops weakly-conserved intron +6); biologically immaterial.

**Tests:** `~SpliceSitePredictor_DonorSite_Tests` = 30 passed / 0 failed (21 canonical-donor incl. 3 added this session + 9 MaxEnt sibling). FULL unfiltered `dotnet test Seqeron.sln -c Debug` = Failed 0 (Seqeron.Genomics.Tests 18790 passed). Build 0 warnings / 0 errors on the changed test file.

---

## Addendum — 2026-06-25: MaxEntScan score5ss opt-in added (`ScoreDonorMaxEnt`)

To complete the MaxEntScan pairing with the sibling `SPLICE-ACCEPTOR-001` (`ScoreAcceptorMaxEnt`,
score3ss), an **opt-in** Yeo & Burge (2004) MaxEntScan **score5ss** maximum-entropy 5' donor
scorer was added: `SpliceSitePredictor.ScoreDonorMaxEnt(string window)` (9-nt window = 3 exon + 6
intron, conserved `GT` at 0-based positions 3–4; returns `log2(P_maxent/P_background)` in bits).

- **Provenance / licence:** the `score5` factorisation and the precomputed probability table
  (`Data/maxent_score5.txt`, 16 384 records = 4^7) were retrieved verbatim this session from the
  **MIT-licensed maxentpy port** (`kepbod/maxentpy`); provenance + full MIT text in
  `Data/maxent_score5.LICENSE.md`. score5 is **single-matrix** (the 7-mer "rest" sequence is keyed
  directly), unlike the nine-sub-matrix score3.
- **Cross-check:** reproduces the documented maxentpy `score5` worked examples EXACTLY —
  `score5('cagGTAAGT') = 10.858313 → 10.86` (canonical), `score5('gagGTAAGT') = 11.078494 → 11.08`,
  `score5('taaATAAGT') = -0.116791 → -0.12`. A wrong table/factorisation fails the 10.86 check.
- **Existing donor scorer unchanged:** `FindDonorSites` / `ScoreDonorSite` / `ScoreU12DonorSite`
  and all defaults are untouched; `ScoreDonorMaxEnt` is purely additive.
- **Status reset:** because production code changed, `SPLICE-DONOR-001` was reset `☑ → ☐` in
  `ALGORITHMS_CHECKLIST_V2.md` (its prior validation no longer covers the new method); Quick
  Reference counts adjusted (Completed 215→214, Not Started 19→20). This unit needs re-validation
  of the added method under the validation campaign.
- **Tests:** 9 new tests (`ScoreDonorMaxEnt_*`, ME1–ME9) added to
  `SpliceSitePredictor_DonorSite_Tests.cs`; full unfiltered `dotnet test` green (Failed: 0).
