# Validation Report: ANNOT-GENE-001 — Gene prediction / Shine-Dalgarno RBS detection

- **Validated:** 2026-06-24   **Area:** Annotation (gene prediction / RBS)
- **Canonical method(s):**
  - `GenomeAnnotator.PredictGenes(dna, minOrfLength, prefix)` — ORF-based gene prediction, both strands
  - `GenomeAnnotator.FindRibosomeBindingSites(dna, upstreamWindow, minDistance, maxDistance)` — forward-strand SD (legacy)
  - `GenomeAnnotator.FindRibosomeBindingSitesBothStrands(...)` — both strands, reports `(position, sequence, score, strand)`
  - private `ScanStrandForShineDalgarno(...)` — shared single-strand scan
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** CLEAN (no defect found; independently re-derived and mutation-checked)

This is an **independent** re-validation (fresh context). Code and its tests were treated as
untrusted: all consensus/spacing values were re-retrieved this session from primary sources, and
all reverse-strand coordinates were re-derived by hand in Python (`/tmp/derive.py`) without lifting
the repo's expected values.

## Stage A — Description

### Sources opened this session (URLs + exact content)

1. **Wikipedia "Shine-Dalgarno sequence"** (WebFetched):
   - Six-base consensus: **"AGGAGG"**.
   - Anti-SD 3' tail of 16S rRNA: **"5'-YACCUCCUUA-3'"** (Y = pyrimidine).
   - Location: **"It is generally located around 8 bases upstream of the start codon AUG."**
   - Note: T4 early genes are GAGG-dominated; the shorter motifs in the code (GGAGG, AGGAG, GAGG,
     AGGA) are all substrings of the AGGAGG consensus, biologically justified.
2. **Chen et al. (1994) NAR 22(23):4953-4957, PubMed 7528374** (WebFetched abstract):
   - **"Measurements of protein expression demonstrated an optimal aligned spacing of 5 nt for both
     series."** Aligned spacing = SD 3' end → A of the AUG.

### Model (validated)
Prokaryotic ORF-based gene finding: scan all 3 frames on both strands for start→stop, filter by
minimum length, label CDS. RBS = AGGAGG (and substrings) found upstream of a start codon at an
**aligned spacing** ∈ [minDistance, maxDistance]. SD is a feature of the mRNA, which for a
reverse-strand gene is the reverse complement of the forward genomic strand; therefore the
reverse-strand SD is found by scanning the reverse complement, not the forward strand.

### Edge-case semantics (sourced/defined)
- Empty/null input → empty (guard present).
- No ORF / no upstream window → no RBS (SD is defined relative to a start codon).
- Spacing at `maxDistance` included, above it excluded (window semantics).
- `protein_length` excludes the stop codon: `(End-Start)/3 − 1`.

### Stage A findings
PASS. Consensus AGGAGG, anti-SD 5'-YACCUCCUUA-3', ~8 nt location, and the 5 nt optimal aligned
spacing all match authoritative sources. The default `[4,15]` window is an honestly-declared
permissive heuristic search window (the method reports candidates; it does not score optimality),
not a description error. No biological error.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`
- `FindOrfs` (91-120): forward + reverse passes; reverse coords mapped to forward via
  `adjStart = len − orf.End`, `adjEnd = len − orf.Start`.
- `FindRibosomeBindingSites` (267-284): forward-only via `.Where(o => !o.IsReverseComplement)`.
- `FindRibosomeBindingSitesBothStrands` (302-344): forward pass + reverse pass; reverse ORFs
  un-mapped to revComp space (`o with { Start = len−o.End, End = len−o.Start }`), scanned on
  `revComp`, motif mapped back via `forwardPosition = len − hit.position − motifLen` (341).
- `ScanStrandForShineDalgarno` (350-384): spacing test `distanceToStart = orf.Start − genomicPos −
  motif.Length` ∈ [minDistance, maxDistance]; score = `motif.Length / 6.0`.
- `PredictGenes` (389-421): `protein_length = ProteinSequence.TrimEnd('*').Length` (excludes stop).

### Reverse-strand coordinate mapping — independently re-derived (NOT from repo tests)
For `CreateSequenceWithReverseStrandSd("AGGAGG", spacing 8)`: source construct
S = `C×10 + AGGAGG + C×8 + (ATG + A×297 + TAA)`, |S| = 10+6+8+303 = **327**. Genomic input = RC(S).
Code scans RC(RC(S)) = S, AGGAGG at revPos = 10 → `forwardPos = 327 − 10 − 6 = `**`311`**.
Independently verified `RC(S)[311..317] = "CCTCCT"` (the anti-SD complement on the forward strand);
aligned spacing in revComp space = 24 − 10 − 6 = **8** ✓. Other positions re-derived:
R4 (spacing 15, |S|=334) → `334 − 10 − 6 = `**`318`**; R3 reverse part at forward offset 327 →
**638**. All match the test expectations and `/tmp/derive.py`.

### Cross-verification table (recomputed vs code)
| Case | Expected (independent) | Test asserts | Code output |
|------|------------------------|--------------|-------------|
| Fwd AGGAGG @ spacing 8 | pos 10, score 1.0 | 10, 1.0 | ✓ |
| Rev AGGAGG @ spacing 8 | fwd pos 311, '-', anti-SD CCTCCT | 311, '-' | ✓ |
| Rev @ maxDist 15 | fwd pos 318 | 318 | ✓ |
| Rev (combined) | fwd pos 638 | 638 | ✓ |
| Spacing 16 (>max) | rejected | empty | ✓ |
| protein_length, minimal_orf(100) | (303)/3−1 = 100 | 100 | ✓ |

### Variant / behaviour consistency
- R7 asserts byte-for-byte parity of `(position, sequence, score)` between legacy forward-only and
  the both-strands `+` hits → no forward-behaviour drift.
- Reverse pass reuses the identical `ScanStrandForShineDalgarno` (no logic duplication).

### Test quality audit + independent mutation check
39 tests (32 legacy + R1-R7 reverse). Assertions are exact positions/sequences/scores/strands or
exact emptiness — no bare "no throw". R1 pre-guards the construct geometry
(`sequence.Substring(311,6) == "CCTCCT"`) so the position assertion cannot silently track a helper
change.

**Independent mutation (by validator):** replaced the reverse mapping with
`forwardPosition = hit.position` → **3 tests failed** (R1/R3/R4), confirming the mapping is genuinely
locked. Source reverted; `git diff` clean; full filtered suite re-run green (39/39).

### Stage B findings
PASS. The code realises the validated biology; the reverse coordinate mapping
`forwardPos = len − revPos − motifLen` is correct (independently re-derived and mutation-falsified);
boundary, no-site, both-strand and protein-length cases are handled; tests are strict, exact, sourced.

## Verdict & follow-ups
- **Stage A:** PASS · **Stage B:** PASS · **End-state:** CLEAN.
- **Filtered suite:** `GenomeAnnotator_Gene_Tests` → 39 passed, 0 failed.
- No code or test change required. No follow-up.
