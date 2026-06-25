# Validation Report: SPLICE-DONOR-001 — Donor (5') Splice Site Detection

- **Validated:** 2026-06-24   **Area:** Splicing
- **Canonical method(s):** `SpliceSitePredictor.FindDonorSites(sequence, minScore, includeNonCanonical)`; internal `ScoreDonorSite(seq, position)` (+ `ScoreU12DonorSite`)
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS

---

## Stage A — Description

### Sources opened & what they confirm

- **Wikipedia, "Shapiro–Senapathy algorithm"** (verified via WebSearch): the 5' splice site consensus is **MAG|GURAGU** (M = A/C, R = purine), spanning position **−3** (third nt from the 3' end of the upstream exon) to **+6** (sixth intron nt). The S&S consensus matrix is built from alignment of 1446 5' splice sites; the algorithm scores authentic and cryptic donor sites by consensus match.
- **GT-AG rule** (multiple sources): **GU is invariant at the 5' splice site** (intron start), AG at the 3' end, for the major (U2) intron class. The GU dinucleotide is required. Donor = GU/GT, NOT AG — no donor/acceptor confusion.
- **Shapiro & Senapathy (1987), NAR 15(17):7155** (Evidence doc + cross-ref): donor consensus (C/A)AG|GU(A/G)AGU; position 0 (G) and +1 (U) ~100% conserved; −1 (G) ~80%; −2 (A) ~60%; −3 (M) A/C ~35% each.
- **Burge, Tuschl & Sharp (1999), The RNA World**: GC-AG U2 introns (~0.5–1%) valid but weaker; U12-type AT-AC introns (~0.3%, donor AT/AU); U12 donor extended motif /ATATCC/. Consistent with Burset, Seledtsov & Solovyev (2000) GC-AG vs GT-AG frequency analysis.
- **Yeo & Burge (2004), MaxEntScan**: 9-mer donor window (−3..+6); higher score = stronger site.

### Formula / convention check

- Canonical 5' donor: intron BEGINS with **GU (GT in DNA)** — confirmed. Spec's `MAG|GURAGU` matches S&S exactly, with invariant G,U at the first two intron positions.
- Coordinate convention: `Position` = index of the **G** of the GU dinucleotide (first intron base / exon-intron junction), 0-based. Standard, unambiguous junction convention.
- IUPAC binary scoring (match = 1.0, no match = 0.0, score = matches / positions scored) — a faithful, conservative reduction of the S&S frequency consensus (S1, S3); documented in the Assumption Register (A1/A2/A3 all eliminated).

### Edge-case semantics

All sourced and defined: empty/null → empty; sequence < 6 nt → empty (insufficient window); no GU → empty; GC valid-but-weaker (gated by `includeNonCanonical`); U12 AT/AU → separate scorer + `U12Donor` type; multiple GU each scored independently.

### Independent cross-check (hand computation)

Worked example `CAGGUAAGU` (perfect consensus, GU at index 3). DonorPwm offsets −3..+5 evaluated at `position+offset`:
−3→idx0 C∈{A,C}=1, −2→idx1 A=1, −1→idx2 G=1, **0→idx3 G=1, +1→idx4 U=1** (invariant GU), +2→idx5 A∈{A,G}=1, +3→idx6 A=1, +4→idx7 G=1, +5→idx8 U=1 → **9/9 = 1.0**, Position = 3. Matches spec.

`CAGGCAAGU` (GC): identical except +1→idx4 C ≠ U → **8/9**, confirming GC donor < GT donor and that the recognised dinucleotide is GU (not AG).

### Findings / divergences

The PWM window is **−3..+5 (9 positions)** vs the literature's −3..+6 (10 positions). It drops the weakly-conserved intron +6 (T) position; biologically immaterial (invariant GU and the conserved −1 G / +3 A / +4 G are retained). No biological error. **PASS.**

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/SpliceSitePredictor.cs`
- `DonorPwm` (102–113): offsets −3..+5, binary weights encoding M,A,G,**G,U**,R,A,G,U. Offset 0 = G (only), offset +1 = U (only) → invariant GU, not AG. Correct.
- `FindDonorSites` (149–204): T→U normalization (line 157), `ToUpperInvariant`, guard `len < 6`, scans `i ≤ len−6`; canonical GU branch (162); `includeNonCanonical` adds GC (Donor) and A[U/T] (U12Donor) branches. Reports `Position = i` (G of GU).
- `ScoreDonorSite` (286–306): sums binary weights over in-bounds offsets, divides by count scored → fraction in [0,1].

### Formula realised correctly?

Yes. Donor dinucleotide recognised = **GU/GT** at offsets 0/+1, matching the GT-AG rule and MAG|GURAGU consensus. GC scores 8/9 automatically (+1 mismatch). U12 AU routed to `/AUAUCC/` scorer, typed `U12Donor`. Position reported at the exon/intron junction (G of GU) — no off-by-one.

### Cross-verification table recomputed vs code (tests assert these exact values)

| Input | GU idx (Position) | Score | Note |
|-------|------|-------|------|
| CAGGUAAGU | 3 | 9/9 = 1.0 | perfect consensus |
| UUUGUAAUU | 3 | 5/9 | weak context |
| CAGGCAAGU (GC, nonCanonical) | 3 | 8/9 | +1 C ≠ U |
| CAGGUAAGU…CAGGUAAGU (multi) | {3, 7, 26} | {1.0, 4/9, 1.0} | full scan |

All recomputed by hand and confirmed by test assertions.

### Variant/delegate consistency

DNA (T) and RNA (U) inputs give identical results (T→U at line 157). Lowercase handled (`ToUpperInvariant`). `CalculateMaxEntScore` reuses the same `DonorPwm` for the donor branch.

### Test quality audit

`tests/Seqeron/Seqeron.Genomics.Tests/SpliceSitePredictor_DonorSite_Tests.cs` — 18 tests (M1–M10, S1–S5, C1–C2, + helper). Assertions check **exact** sourced values (scores 9/9, 8/9, 5/9, 4/9; positions 3/7/26; type Donor vs U12Donor), not tautologies. Independent helper `ComputeExpectedConsensusStrength` re-derives consensus strength separately from the implementation. Edge cases (null, empty, <6, no-GU, GC gating, U12) all covered.

### Findings / defects

None. Donor dinucleotide, consensus positions, GU invariance, junction coordinate, GC/U12 handling all match the validated description.

---

## Verdict & follow-ups

- **Stage A: PASS** — donor = GT/GU at intron start (GT-AG rule); consensus MAG|GURAGU confirmed against Shapiro & Senapathy (1987), Burge et al. (1999), Wikipedia (Shapiro–Senapathy algorithm); coordinate convention sound.
- **Stage B: PASS** — implementation recognises GU at the donor, scores against the IUPAC consensus, reports the junction position correctly; worked examples reproduce exactly.
- **State: CLEAN** — no defect found. No code changes required.
- Note (non-blocking): PWM uses −3..+5 (drops weakly-conserved intron +6); biologically immaterial, documented in the Assumption Register.

**Tests:** `~SpliceSitePredictor_DonorSite_Tests` = 18 passed / 0 failed. Build succeeded, 0 warnings. No code touched, so no full-suite run required.

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
