# Validation Report: CHROM-CENT-001 — Centromere classification + α-satellite suprachromosomal-family assignment

- **Validated:** 2026-06-26   **Area:** Chromosome
- **Re-validation trigger:** limitation-fix commit `887a9945` ADDED suprachromosomal-family (SF) assignment
  (`AssignSuprachromosomalFamily` + `LoadBundledAlphaSatelliteReference` + new types). Prior validation
  SUPERSEDED; unit re-validated fresh this session with primary focus on the new surface and a confirmation
  pass over the pre-existing canonical surface.
- **Canonical / surface validated this session:**
  - **NEW (primary focus):** `ChromosomeAnalyzer.AssignSuprachromosomalFamily(sequence, reference=null)`,
    `LoadBundledAlphaSatelliteReference()`, types `SuprachromosomalFamily`, `AlphaSatelliteReferenceMonomer`,
    `SuprachromosomalFamilyResult`, `AlphaSatelliteBoxType`; bundled CC0 Dfam reference
    `Resources/AlphaSatelliteReference.fasta` (ALR / ALRa / ALRb).
  - **Pre-existing (confirmation):** Levan classification (`CalculateArmRatio`,
    `ClassifyChromosomeByArmRatio`, `DetermineCentromereType`), α-satellite / CENP-B detection
    (`DetectAlphaSatellite`, `FindCenpBBoxes`), HOR detection (`DetectHigherOrderRepeat`).
- **Source file:** `src/Seqeron/Algorithms/Seqeron.Genomics.Chromosome/ChromosomeAnalyzer.cs`
- **Test files:** `ChromosomeAnalyzer_SuprachromosomalFamily_Tests.cs` (13, NEW),
  `ChromosomeAnalyzer_Centromere_Tests.cs`, `ChromosomeAnalyzer_MutationKillers_Tests.cs`,
  `ChromosomeAnalyzerTests.cs` (275 total under `~ChromosomeAnalyzer`, all green).
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **State:** ✅ CLEAN (no defect found; no code changed this session)

---

## Stage A — Description

### Sources opened THIS session & what they confirm

- **Dfam REST API** (retrieved this session, verbatim `consensus_sequence`):
  - `DF000000029` **ALR** — length 171, classification `root;Tandem_Repeat;Satellite;Centromeric`.
  - `DF000000014` **ALRa** — length 172, same classification (≈83.7% to ALR).
  - `DF000000015` **ALRb** — length 169, same classification.
  All three consensus strings match the bundled `Resources/AlphaSatelliteReference.fasta` **byte-for-byte**
  (verified this session — see Stage B).
- **Dfam licence** (Storer, Hubley, Rosen et al. 2021, *Mobile DNA* 12:2; dfam.org) — Dfam data is
  dedicated to the **public domain under CC0**; freely redistributable. The bundled reference is genuinely CC0.
- **McNulty SM & Sullivan BA (2018)** "Alpha satellite DNA biology", *Chromosome Res* (PMC6121732) —
  confirmed **verbatim**:
  - SF1 dimeric (J1·J2; HSA1,3,5,6,7,10,12,16,19); SF2 dimeric (D1·D2; eleven chromosomes);
    SF3 pentameric (W1–W5; HSA1,11,17,X); SF4 monomeric (M1; acrocentric p-arms + Y); SF5 irregular (R1·R2).
  - A/B-box rule **verbatim:** "A-type monomers include J1, D2, W4, W5, M1, and R2 monomers, while B-type
    consist of J2, D1, W1–W3, and R1 monomers." and "B-type monomers contain CENP-B boxes, while A-type
    monomers contain a binding site for pJα."
- **Shepelev et al. (2009)** PLOS Genet 5:e1000641 — SF taxonomy origin confirmed via cross-source:
  SF1 = J1+J2, SF2 = D1+D2, SF3 = W1–W5, SF4 = M1, SF5 = R1+R2 (twelve SF monomer types in five families).
- **Masumoto et al. (1989)** CENP-B box 5′-`YTTCGTTGGAARCGGGA`-3′ (Y=C/T, R=A/G) — the 17-bp motif; used by
  `FindCenpBBoxes`.
- **Levan A, Fredga K, Sandberg AA (1964)** *Hereditas* 52(2):201–220 — arm-ratio r = L/S and centromeric
  index ci = 100·p/(p+q); cut-points 1.7 / 3.0 / 7.0 (ci 37.5 / 25.0 / 12.5) over m / sm / st / a / T.
  (Re-confirmed; unchanged from prior round.)

### A/B-box typing of the bundled monomers (sequence-verifiable, hand-checked this session)

Scanning the retrieved Dfam consensus strings for the IUPAC CENP-B box `[CT]TTCGTTGGAA[AG]CGGGA`:
- **ALRb** → one match at **0-based position 126** (`CTTCGTTGGAAACGGGA`) → **B-type**. ✔ (matches the
  bundled `boxtype=B` and the README claim "CENP-B box at consensus position 126").
- **ALR** and **ALRa** → **no match** → **A-type**. ✔
This reproduces the sourced rule (B-type carries the CENP-B box; A-type does not).

### SF assignment rule vs. the published families

The method has only two reproducible SF-determining signals available from the CC0 reference: the **HOR
period** (from `DetectHigherOrderRepeat`) and the **A/B-box composition** of one HOR unit. The mapping is a
faithful reduction of the Shepelev/McNulty taxonomy onto those two signals:

| Signal | Code rule | Published basis | Verdict |
|---|---|---|---|
| period multiple of 5 | → SF3 | SF3 is pentameric (W1–W5) | ✔ diagnostic |
| period 1, all A-type | → SF4 | SF4 = M1, monomeric, A-type | ✔ |
| period 2 (A+B dimer) | → {SF1,SF2} | SF1/SF2 both dimeric A→B | ✔ (cannot separate from CC0) |
| irregular A/B, both box types | → SF5 | SF5 = R1·R2 irregular | ✔ |
| otherwise | → Unknown | — | ✔ |

### Edge-case semantics (sourced / defined)
- Empty / null / shorter than one 171-bp monomer → not alpha-satellite, `Unknown`, period 0, empty pattern.
- A query below the **≥60% identity** gate to the closest reference monomer → not alpha-satellite. The gate
  is justified by alphoid divergence (~16% from consensus, 50–70% between monomer classes; Waye & Willard
  1987; PMC6121732) — random DNA sits well below it (measured ~52% this session), real monomers ~85–98%.

### Stage A verdict
**✅ PASS.** Every SF fact, the A/B-box rule, the CENP-B box motif/position, the three Dfam sequences, and the
CC0 licence are confirmed against authoritative external first-sources retrieved this session. The reduction
of the taxonomy onto {HOR period + A/B composition} is sound, and the residual (SF1-vs-SF2 and non-period-5
SF3) is an honest open boundary, not a hidden defect — separating them needs the SF-resolved consensus
monomer library (J1/J2/D1/D2/W1–W5/M1/R1/R2), which is not CC0/redistributable.

---

## Stage B — Implementation

### Reference data integrity (byte-for-byte vs Dfam)
- `AlphaSatelliteReference.fasta`: the ALR/ALRa/ALRb consensus lines are **identical** to the Dfam REST
  `consensus_sequence` fields retrieved this session (lengths 171/172/169; the `N`-containing positions
  preserved). `boxtype` headers (A/A/B) match the sequence-derived CENP-B typing.
- `LoadBundledAlphaSatelliteReference` parses the FASTA from the embedded resource (csproj
  `EmbeddedResource`), upper-cases, and returns three records with the right `Name`/`Accession`/`Sequence`/
  `BoxType` — confirmed by independent probe (`ref count=3`; ALR/ALRa/ALRb; acc/len/box all correct).

### Code path reviewed (`ChromosomeAnalyzer.cs`)
- `AssignSuprachromosomalFamily` (`:1090`): guards empty/null/<171; splits into `len/171` whole monomers
  (trailing partial ignored); best-matches each monomer to every reference via
  `SequenceAligner.GlobalAlign` → `CalculateStatistics().Identity`; accepts a monomer as alpha-satellite
  when best identity ≥ 60%; `matchedCount==0` → not-alpha/`Unknown`; takes period from
  `DetectHigherOrderRepeat`; builds the first-unit A/B pattern; delegates to `ClassifyFamily`.
- `ClassifyFamily` (`:1187`): period%5==0→SF3; period 1 → SF4 (all-A) else SF5; period 2 → {SF1,SF2};
  else both-box-types → SF5; else Unknown. Matches the Stage-A table exactly.

### Independent cross-check — SF assignments reproduced this session (own probe, not test echoes)

| Input (real Dfam-derived monomers) | IsAlpha | Family | period | best | meanId | pattern |
|---|---|---|---|---|---|---|
| ALRa ×8 (monomeric A) | true | **Sf4** | 1 | ALRa | 95.5% | [A] |
| (ALRa+ALRb) ×6 (dimer A·B) | true | **Sf1OrSf2Dimeric** | 2 | ALRa | 97.5% | [A,B] |
| pentamer 3B+2A ×6 (W1–W3 B, W4–W5 A) | true | **Sf3** | 5 | ALRb | 84.7% | [B,B,B,A,A] |
| irregular A/B array | true | **Sf5** | 1 | ALRa | 97.0% | (irregular) |
| random DNA (seed 1234, 400 bp) | **false** | **Unknown** | 0 | null | 51.9% | [] |
| ALRb ×8 (monomeric B-only) | true | **Sf5** | 1 | ALRb | 92.2% | [B] |

`FindCenpBBoxes`: ALRb→`[126]`, ALRa→`[]`, ALR→`[]`. All values match the sourced expectations and the
test assertions; the negative case (random DNA, 51.9% < 60% gate → not alpha-satellite) is genuine.

### Pre-existing surface (confirmation pass)
- **Levan** (`ClassifyChromosomeByArmRatio` / `DetermineCentromereType`): normalises to r = max(p/q, q/p),
  switch ≤1.7 m / ≤3.0 sm / <7.0 st / else a; p=0/ratio≤0 → Telocentric. Correct per Levan 1964 (this was
  the defect fixed in the prior round; it remains correct and the lock-tests trace to hand-computed r).
- **α-satellite / CENP-B / HOR detectors:** unchanged by commit `887a9945` (additive only). The new method's
  M-SF-8 test asserts `DetectAlphaSatellite`/`DetectHigherOrderRepeat` are byte-identical before/after a call
  — the additive contract holds.
- All **275** `~ChromosomeAnalyzer` tests green.

### Test-quality audit (HARD gate)
- All public surface covered: `LoadBundledAlphaSatelliteReference` (M-SF-1/2), `AssignSuprachromosomalFamily`
  default + caller-supplied reference overload (M-SF-3..8, S-SF-1), case-insensitivity (S-SF-2), identity
  sanity (C-SF-1), edge cases (empty/null/short, empty-reference→throws).
- Assertions are **sourced, not code echoes**: A/B typing traces to the PMC6121732 rule + the CENP-B box at
  position 126; SF families trace to the dimeric/pentameric/monomeric/irregular taxonomy; identity to the
  measured Dfam-derived values; the negative case to the 60% gate. Reference monomers are the REAL Dfam
  strings (period-5 case uses deterministic mild point-variants that stay above the gate and keep parent
  A/B type).
- Deterministic (fixed local RNG seeds; the shared-RNG hazard does not apply — RNG is local to each test).

### Stage B verdict
**✅ PASS.** The code faithfully realises the validated description; the bundled reference is byte-identical
to CC0 Dfam; A/B typing is correct (ALRb carries the CENP-B box at 126); the SF rule reproduces the published
family signals on worked examples plus a real negative case. No defect found; no code changed this session.

---

## Verdict & follow-ups
- **Stage A:** ✅ PASS — SF taxonomy, A/B-box rule, CENP-B box, the three Dfam sequences + CC0 licence, and
  Levan thresholds confirmed against external first-sources retrieved this session.
- **Stage B:** ✅ PASS — reference byte-verified vs Dfam; SF assignments + CENP-B positions reproduced
  independently; additive contract holds; tests sourced.
- **State:** ✅ **CLEAN.** Full unfiltered `dotnet test Seqeron.sln -c Debug` **Failed: 0**
  (Seqeron.Genomics.Tests 18860 passed); no code/test changed this session.
- **Documented open boundary (acceptable, not a defect):** SF1 vs SF2 (both dimeric A→B) and SF3 arrays whose
  period is not a multiple of 5 (e.g. dodecameric DXZ1) are not resolved from the CC0 reference; they need the
  SF-resolved consensus monomer library (not CC0/redistributable). Recorded in `Resources/README.md`,
  `LIMITATIONS.md` §2, and the TestSpec residual.

---

## Historical note
The prior round (2026-06-25, Levan focus) found and fully fixed a real defect in
`ClassifyChromosomeByArmRatio` (diverged from Levan 1964; omitted Subtelocentric) and de-green-washed its
tests. That fix remains in place and re-confirmed here. The present session re-validates the unit after the
SF-assignment addition (commit `887a9945`), which is additive and does not touch the Levan surface.

## Runtime enforcement (LimitationPolicy)

Under the default `LimitationPolicy.DefaultMode = Strict`, the unresolved `Sf1OrSf2Dimeric` (SF1-vs-SF2) call throws `Seqeron.Genomics.Core.SeqeronLimitationException` (named limitation + workaround; see [LIMITATIONS.md](../LIMITATIONS.md) › Runtime enforcement and `LimitationCatalog`). `Permissive` mode returns the historical best-effort value. This is an additive policy layer; the validated contract and `✅ CLEAN` verdict are unchanged.
