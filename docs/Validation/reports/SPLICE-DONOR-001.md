# Validation Report: SPLICE-DONOR-001 — Donor (5') Splice Site Detection

- **Validated:** 2026-06-12   **Area:** Splicing
- **Canonical method(s):** `SpliceSitePredictor.FindDonorSites(sequence, minScore, includeNonCanonical)`; internal `ScoreDonorSite(context)` (+ `ScoreU12DonorSite`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm

- **Wikipedia, "RNA splicing"** (https://en.wikipedia.org/wiki/RNA_splicing): the splice **donor** site includes an almost invariant **GU** at the 5' end of the intron (GT in DNA). Extended donor consensus `G‑G‑[cut]‑G‑U‑R‑A‑G‑U` (R = purine). The **GU‑AG rule** (GU at 5' end, AG at 3' end of intron) accounts for >99% of splicing (canonical). Non-canonical = minor spliceosome.
- **Wikipedia, "Minor spliceosome / U12-type introns"**: U12-type introns are not exclusively AT‑AC; donor consensus historically `ATATCC…`, 3' `YCCAC`. Distinguished mainly by 5' SS and branch site, not just terminal dinucleotides.
- **Shapiro & Senapathy (1987), NAR 15(17):7155** (confirmed via Springer/NAR cross-references and the Shapiro–Senapathy algorithm page): donor consensus **MAG|GURAGU** (M = A/C, R = A/G), positions **−3 to +6** relative to the exon/intron junction; the **GU at the first two intron bases is essentially invariant**; position −3 prefers C/A; PWM-based scoring.
- **Burge, Tuschl & Sharp (1999)** and **Yeo & Burge (2004)**: GC‑AG U2 introns (~0.5–1%) valid but weaker; U12 AT‑AC (~0.3%); higher PWM/log-odds score = stronger site.

### Formula / convention check

- Canonical 5' donor: intron BEGINS with **GT (GU in RNA)** — confirmed (GT‑AG rule). The spec's `MAG|GURAGU` matches Shapiro & Senapathy exactly, with the invariant G,U at the first two intron positions.
- Coordinate convention: the spec/code report the donor **Position = index of the G** of the GU dinucleotide (first intron base / exon‑intron junction), 0-based. This is a standard, unambiguous junction convention.
- Donor dinucleotide is **GT/GU**, NOT AG — no donor/acceptor confusion. GC‑AG and U12 AT/AU handling present and gated behind `includeNonCanonical`.

### Edge-case semantics

All sourced and defined: empty/null → empty; sequence < 6 nt → empty (insufficient window); no GU → empty; GC valid-but-weaker; multiple GU each scored independently.

### Independent cross-check (hand computation)

For `...exon] GU AAGU [intron...` worked example `CAGGUAAGU` (perfect consensus, GU at index 3):
offsets −3..+5 map to C,A,G,G,U,A,A,G,U →
−3 C∈{A,C}=1, −2 A=A=1, −1 G=G=1, **0 G=G=1, +1 U=U=1** (invariant GU), +2 A∈{A,G}=1, +3 A=A=1, +4 G=G=1, +5 U=U=1 → **9/9 = 1.0**.
Spec-expected donor position = 3, score = 1.0 — matches.

### Findings / divergences

The spec window is **−3..+5 (9 positions)** vs the literature's −3..+6 (10 positions, i.e. through intron +6). This drops the weakly-conserved +6 (T) position; it is biologically immaterial for donor recognition (the discriminating, invariant signal GU and the conserved −1 G / +3 A / +4 G / +5 G/T are retained). The binary IUPAC weighting (match=1, else 0) is a faithful, conservative reduction of the Shapiro–Senapathy frequency consensus, explicitly documented in the Assumption Register. No biological error. PASS.

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs`
- `DonorPwm` (lines 102–113): offsets −3..+5, binary weights encoding M,A,G,**G,U**,R,A,G,U. Invariant G (offset 0) and U (offset +1) are the only allowed bases at the dinucleotide — correct (GT, not AG).
- `FindDonorSites` (149–204): T→U normalization, ToUpperInvariant, guard `len<6`, scans `i ≤ len−6`; canonical GU branch; `includeNonCanonical` adds GC (Donor) and A[U/T] (U12Donor) branches. Reports `Position = i` (G of GU).
- `ScoreDonorSite` (286–306): sums binary weights over in-bounds offsets, divides by count scored → fraction in [0,1].

### Formula realised correctly?

Yes. The dinucleotide recognised at the donor is **GU/GT** at offsets 0/+1, matching the validated GT‑AG rule and the MAG|GURAGU consensus. GC scores lower automatically (offset +1 C mismatches invariant U → max 8/9). U12 AU routed to a separate `/AUAUCC/` consensus scorer and typed `U12Donor`. Position reported at the exon/intron junction (G of GU) — no off-by-one.

### Cross-verification table recomputed vs code (tests pass on these exact values)

| Input | GU idx (Position) | Score | Note |
|-------|------|-------|------|
| CAGGUAAGU | 3 | 9/9 = 1.0 | perfect consensus |
| UUUGUAAUU | 3 | 5/9 | weak context |
| CAGGCAAGU (GC, nonCanonical) | 3 | 8/9 | +1 C≠U |
| CAGG**U**AAGU…CAGG**U**AAGU (multi) | {3, 7, 26} | {1.0, 4/9, 1.0} | full scan |

All recomputed by hand and confirmed by the test assertions.

### Variant/delegate consistency

DNA (T) and RNA (U) inputs give identical results (T→U at line 157). Lowercase handled. `CalculateMaxEntScore` uses the same `DonorPwm` for the donor branch.

### Test quality audit

`tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs` — 18 spec tests + helper, 24 total executed. Assertions check **exact** sourced values (score 9/9, 8/9, 5/9, 4/9; position 3/7/26; type Donor vs U12Donor), not tautologies. An independent helper (`ComputeExpectedConsensusStrength`) re-derives the consensus strength separately from the implementation. Edge cases (null, empty, <6, no-GU, GC gating, U12) all covered.

### Findings / defects

None. Donor dinucleotide, consensus positions, GU invariance, junction coordinate, and GC/U12 handling all match the validated description.

---

## Verdict & follow-ups

- **Stage A: PASS** — donor = GT/GU at intron start (GT‑AG rule); consensus MAG|GURAGU confirmed against Shapiro & Senapathy (1987), Burge et al. (1999), Wikipedia; coordinate convention sound.
- **Stage B: PASS** — implementation recognises GU at the donor, scores against the IUPAC consensus, reports the junction position correctly; worked examples reproduce exactly.
- **State: CLEAN** — no defect found. No code changes required.
- Note (non-blocking): PWM uses −3..+5 (drops weakly-conserved intron +6); biologically immaterial, documented.

**Tests:** `~DonorSite` filter = 24 passed / 0 failed. Full suite `Seqeron.Genomics.Tests` = 4486 passed, 0 failed.
