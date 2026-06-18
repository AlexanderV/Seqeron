# Validation Report: ANNOT-GENE-001 — Shine-Dalgarno / Ribosome-Binding-Site detection (reverse-strand enhancement)

- **Validated:** 2026-06-17   **Area:** Annotation (gene prediction / RBS)
- **Re-validated enhancement:** commit `2e29eff` — reverse-strand Shine-Dalgarno reporting
  (new `GenomeAnnotator.FindRibosomeBindingSitesBothStrands` overload + shared
  `ScanStrandForShineDalgarno` factoring).
- **Canonical method(s):**
  - `GenomeAnnotator.FindRibosomeBindingSites(dna, upstreamWindow, minDistance, maxDistance)` — forward only (legacy, unchanged behaviour)
  - `GenomeAnnotator.FindRibosomeBindingSitesBothStrands(dna, upstreamWindow, minDistance, maxDistance)` — both strands, `(int position, string sequence, double score, char strand)` (NEW)
  - private `ScanStrandForShineDalgarno(...)` — shared single-strand scan
- **Stage A verdict:** ✅ PASS
- **Stage B verdict:** ✅ PASS
- **End-state:** ✅ CLEAN (no defect found; existing tests confirmed strict and falsifiable)

This is an **independent** re-validation. The implementer's code and its tests were treated as
untrusted; all consensus/spacing numbers were re-retrieved this session, and all reverse-strand
positions were re-derived by hand (Python, `/tmp/derive.py`) without lifting the repo's expected
values.

## Stage A — Description

### Sources opened this session (URLs + exact content)

1. **Shine-Dalgarno consensus & anti-SD** — Wikipedia "Shine-Dalgarno sequence"
   (https://en.wikipedia.org/wiki/Shine-Dalgarno_sequence), WebFetched:
   - Six-base consensus: **"AGGAGG"**.
   - Anti-SD 3' tail of 16S rRNA: **"5'-YACCUCCUUA-3'"** (Y = pyrimidine).
   - Location: **"around 8 bases upstream of the start codon AUG"**.
   - The SD is described as a feature of the **mRNA** (it base-pairs with 16S rRNA), i.e. it lies
     on the transcribed strand. (Primary: Shine & Dalgarno (1975) *Nature* 254:34-38.)
2. **Optimal aligned spacing** — Chen et al. (1994) *NAR* 22(23):4953-4957,
   PubMed 7528374 (https://pubmed.ncbi.nlm.nih.gov/7528374/), WebFetched abstract:
   - **"an optimal aligned spacing of 5 nt for both series"** (aligned spacing = SD 3' end to the
     A of the AUG initiation codon). Functional window in the literature spans ~4-9 nt; the code's
     default `[minDistance=4, maxDistance=15]` is a permissive superset of this and is honestly a
     heuristic search window, not a claim that 15 nt is optimal.
3. **Reverse-strand transcription orientation** — Wikipedia "Sense (molecular biology)"
   (https://en.wikipedia.org/wiki/Sense_(molecular_biology)), WebFetched:
   - For a minus-strand gene the **mRNA is the reverse complement of the plus (forward) genomic
     strand** in that region ("the sequence of the RNA transcript will look identical to the
     positive-sense strand, apart from … uracil instead of thymine"; the minus strand is the
     template).

### Biology of the enhancement (validated)

A reverse-strand gene's mRNA = reverse complement of the forward genomic strand. The SD motif lies
on that mRNA, upstream (5') of the reverse-strand start codon. Therefore, to find a reverse-strand
SD one must scan the **reverse complement** of the genomic sequence for AGGAGG upstream of the
reverse-strand ORF starts — exactly what the enhancement does. The forward-strand genomic bases at
a reverse SD locus read the **anti-SD complement** (AGGAGG → CCTCCT), which is *not* what the
forward-only scan looks for; hence the two strands are correctly independent.

### Edge-case semantics

- **No ORF / no upstream window** → no hit (correct; SD is defined relative to a start codon).
- **Spacing at/over the window boundary** → included at `maxDistance`, excluded above it
  (aligned-spacing window semantics).
- **Empty input** → empty (guard present).

### Stage A findings

PASS. Consensus, anti-SD, spacing and reverse-strand orientation all match authoritative sources.
The only soft point — the wide `[4,15]` default window vs the literature's ~5 nt optimum — is an
honestly-declared heuristic search window (the method *reports candidates*, it does not score
optimality), not a description error. No biological error.

## Stage B — Implementation

### Code path reviewed
`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`
- `FindRibosomeBindingSites` (lines 262-279) — forward-only; `.Where(o => !o.IsReverseComplement)`.
- `FindRibosomeBindingSitesBothStrands` (lines 297-339) — forward pass + reverse pass.
- `ScanStrandForShineDalgarno` (lines 345-379) — shared scan; spacing test
  `distanceToStart = orf.Start - genomicPos - motif.Length` ∈ `[minDistance, maxDistance]`;
  score = `motif.Length / 6.0`.
- `FindOrfs` reverse-coordinate mapping (lines 109-119): reverse ORFs are found in `revComp` space
  then mapped to forward via `adjStart = len - end`, `adjEnd = len - start`.

### Forward→reverse coordinate mapping — independently re-derived (NOT taken from repo tests)

Reverse pass builds `revComp = RC(dna)` and **undoes** the FindOrfs forward mapping via
`o with { Start = len - o.End, End = len - o.Start }` to recover revComp-space ORF coords
(round-trip: `len - (len - origStart) = origStart` ✓). It then scans `revComp`, obtaining a motif at
revComp-space index `revPos`, and reports `forwardPos = len - revPos - motifLen`.

Hand derivation for the R1 construct (`CreateSequenceWithReverseStrandSd("AGGAGG", spacing 8)`):
the source SD+ORF construct S = `C×10 + AGGAGG + C×8 + (ATG + A×297 + TAA)`, |S| = 10+6+8+303 = **327**;
the genomic input is `RC(S)`, |RC| = 327. The code scans `RC(RC(S)) = S`, finds AGGAGG at revPos = 10,
so `forwardPos = 327 - 10 - 6 = `**`311`**. Independently I verified `RC(S)[311..317] = "CCTCCT"`
(the anti-SD complement on the forward strand). This matches the code's output and the test's
expectation — derived independently in `/tmp/derive.py`, not lifted from the repo.

Other re-derived positions: R4 (spacing 15, |S| = 334) → `334 - 10 - 6 = `**`318`**; R3 reverse part
at forward offset 327 → reverse hit at **638**. All confirmed.

The mapping is correct: `revComp[revPos, revPos+L)` corresponds to forward
`[len-revPos-L, len-revPos)`, whose 5' (lowest) forward index is `len - revPos - L`.

### Spacing semantics check
`distanceToStart = orf.Start - genomicPos - motif.Length` equals the Chen "aligned spacing"
(SD 3' end → AUG first base). Verified numerically: for the test constructs `distanceToStart`
parameter == computed aligned spacing for spacings 4/15/16. So R4 (15 → detected at the boundary)
and R5 (16 → rejected) correctly pin `maxDistance=15`.

### Variant / behaviour consistency
- Legacy `FindRibosomeBindingSites` signature & forward-only behaviour preserved (R7 asserts
  byte-for-byte parity of `(position, sequence, score)` with the both-strands `+` hits).
- Reverse pass is the *same* scan on the reverse complement (no logic duplication / drift).

### Test quality audit (HARD gate)

`tests/.../GenomeAnnotator_Gene_Tests.cs` reverse region (R1-R7) + 32 legacy tests = 39 total.

| Gate criterion | Result |
|---|---|
| Forward-strand hit covered | ✅ R7 (parity), R3 (forward hit @10 in combined construct) |
| Reverse-strand hit covered, exact position | ✅ R1 (311), R3 (638), R4 (318) — all match my independent derivation |
| No-site / no-ORF | ✅ R6 (no ORF → empty) |
| Boundary spacing | ✅ R4 (=maxDistance detected), R5 (>maxDistance rejected) |
| Assertions trace to sourced values, not code echoes | ✅ positions hand-derived from construct geometry; consensus AGGAGG, score 6/6=1.0 sourced |
| Wrong strand mapping must fail a test | ✅ verified by mutation (below) |
| No green-washing (no bare "no throw") | ✅ every test asserts exact position/sequence/score/strand or exact emptiness |

**Mutation testing (independent, by the validator):**
- Mutant A — `forwardPos = hit.position` (drop the reverse mapping): **3 tests fail** (R1, R3, R4).
- Mutant B — off-by-one `forwardPos = len - hit.position - motifLen + 1`: **3 tests fail** (R1, R3, R4).
- Both reverted; `git diff` on the source = clean.

The R1 anti-SD guard `sequence.Substring(311,6) == "CCTCCT"` independently fixes the construct
geometry before the method is even called, so the position assertion cannot silently track a code
change in the helper.

Minor note (not a defect): R2 ("forward-only ignores reverse SD") asserts the legacy method returns
empty on a reverse-only construct. The empty result follows from there being no forward ORF rather
than from explicit strand filtering, so R2 alone is a weak separation guard — but R3 (combined
construct: forward-only must report only the `+` hit @10 and the both-strands method must add the
`-` hit @638) and R7 (exact forward parity) together fully lock the forward/reverse separation.
Coverage is adequate; no test added.

### Stage B findings
PASS. The code faithfully realises the validated biology; the coordinate mapping is correct
(independently re-derived); the boundary and no-site cases are handled; tests are strict, exact,
sourced, and proven falsifiable by mutation. No code or test change required.

## Verdict & follow-ups

- **Stage A:** ✅ PASS · **Stage B:** ✅ PASS · **End-state:** ✅ CLEAN.
- **Full unfiltered suite:** `dotnet test` (no filter) → **6768 passed, 0 failed, 0 skipped-as-fail**,
  build 0 errors (4 pre-existing NUnit2007 warnings in unrelated `ApproximateMatcher_EditDistance_Tests.cs`).
- No defect logged. No follow-up required.
