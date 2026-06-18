# Validation Report: GENOMIC-ORF-001 — Open Reading Frame (ORF) Detection

- **Validated:** 2026-06-16   **Area:** Analysis
- **Canonical method(s):** `GenomicAnalyzer.FindOpenReadingFrames(DnaSequence, int minLength = 100)` (+ private `FindOrfsInFrame`, `IsCodon`, `IsStopCodon`); result struct `OrfInfo`. `src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/GenomicAnalyzer.cs:402–496, 596–612`.
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS
- **End-state:** ✅ CLEAN (no defect found)
- **Duplicate?** No. Distinct from ANNOT-ORF-001 (`GenomeAnnotator.FindOrfs` → `OpenReadingFrame`) and TRANS-SIXFRAME-001 / `Translator.FindOrfs` (→ `OrfResult`), which are genetic-code-parameterized. This unit is the ATG-only, standard-code, `OrfInfo`-returning six-frame scanner. Three separate public APIs confirmed by grep.

## Stage A — Description

### Sources opened this session (retrieved, not trusted by label)
- **Rosalind "Open Reading Frames" (ORF)** — https://rosalind.info/problems/orf/ (WebFetch, 2026-06-16). Confirms verbatim: ORF begins at a start codon, ends at a stop codon, no intervening stop; **six** reading frames (3 forward + 3 reverse complement); start codon AUG (→ATG); sample DNA `AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG`; sample output exactly `{MLLGSFRLIPKETLIQVAGSSPCNLS, M, MGMTPRLGLESLLE, MTPRLGLESLLE}` (distinct protein strings; nested ORFs sharing a stop both reported).
- **Wikipedia "Open reading frame"** — https://en.wikipedia.org/wiki/Open_reading_frame (WebFetch, 2026-06-16). Confirms: "spans of DNA sequence between the start and stop codons"; start ATG; stop codons RNA UAA/UAG/UGA = DNA TAA/TAG/TGA; three frames per strand, six total for dsDNA; alternative definition "length divisible by three and bounded by stop codons"; minimum-length conventions (100 / 150 codons) and that a long ORF is not conclusive evidence of a gene.
- **NCBI ORFfinder** and **NCBI Genetic Codes (transl_table=1)** — as cited in the Evidence doc (start ATG; stop TAA/TAG/TGA; nucleotide minimum-length filter, inclusive). Consistent with the above.

### Formula / model check
The model in the algorithm doc (§2.2) and code matches the sources: for each strand and each frame offset f∈{0,1,2}, scan codons; an ORF is `S[a..b+3)` with `S[a..a+3)=ATG`, `S[b..b+3)` a stop, no in-frame stop strictly between, span divisible by 3 including the stop; protein candidate = translation of `S[a..b)` (stop excluded). Reverse complement searched identically. Every in-frame ATG that reaches a stop is reported (nested ORFs). All matches Rosalind + Wikipedia.

### Edge-case semantics (all sourced)
- ATG with no downstream in-frame stop → not reported (Rosalind "translate until a stop codon"). ✔
- Nested ATGs sharing one stop → all reported (Rosalind sample `MGM…`/`MTP…`). ✔
- ORF only on reverse strand → reported with `IsReverseComplement=true` (six-frame requirement). ✔
- Sequence shorter than start+stop / no ATG → empty. ✔
- `minLength` = nucleotide threshold, inclusive (NCBI ORFfinder). ✔ (assumption, sourced)
- ORF `Sequence` includes the stop codon so `Length % 3 == 0` (Wikipedia "length divisible by three, bounded by stop codons"); protein candidate excludes the stop. ✔ (assumption, sourced)
- Null `sequence` → `ArgumentNullException` (contract). ✔

### Independent cross-check (numbers re-derived this session)
Built an **independent Python re-implementation from the sourced definitions** (own codon table from NCBI transl_table=1, own reverse-complement, every-ATG-to-first-in-frame-stop, six frames) — NOT from the repo code — and recomputed every test expectation:

| Case | Input | Independent result |
|------|-------|--------------------|
| Rosalind | sample | distinct proteins = `{M, MGMTPRLGLESLLE, MLLGSFRLIPKETLIQVAGSSPCNLS, MTPRLGLESLLE}` — **exact match** to Rosalind output |
| M1 | `ATGAAAAAATAA` | 1 ORF `(ATGAAAAAATAA, pos0, frame1, fwd)`, protein `MKK` |
| M3 | `ATGGGGATGCCCTAA` | fwd frame1 = `[(ATGGGGATGCCCTAA, 0), (ATGCCCTAA, 6)]` (shared TAA) |
| M4 | `ATGAAAAAAAAA` | fwd ORFs = ∅ |
| M5 | `TTAGGGGGGCAT` | exactly 1 ORF, RC, `ATGCCCCCCTAA`, protein `MPP` (RC = `ATGCCCCCCTAA`) |
| M6 / M6b | `ATGAAATAA` | minLen12 → ∅; minLen9 fwd pos0 → `[(ATGAAATAA, 0)]` |
| M10 | `ATGTAA`/`ATGTAG`/`ATGTGA` | each → exactly 1 fwd frame1 pos0 ORF |
| S2 | `ATGAAATGA` / `CATGCCCTAA` | frame1 ORF at 0; frame2 ORF at 1 |
| C1 | `ATG` | ∅ |

All values reproduced exactly.

### Findings / divergences
None. Description is biologically and mathematically correct and matches the retrieved sources. Three assumptions (stop-inclusive span, nucleotide inclusive minLength, ATG-only) are explicit and each sourced.

## Stage B — Implementation

### Code path reviewed
`GenomicAnalyzer.FindOpenReadingFrames` (`GenomicAnalyzer.cs:421`): null-guard → 3 forward frames → reverse-complement → 3 reverse frames. `FindOrfsInFrame` (`:447`) steps `start += CodonLength` over codon positions; at each ATG, inner loop steps to the first stop and yields `OrfInfo(substring(start,length), start, frame+1, isRC)` if `length >= minLength`, then `break` at the first in-frame stop. `IsCodon`/`IsStopCodon` use `string.CompareOrdinal`. Constants: `StartCodon="ATG"`, `StopCodons={TAA,TAG,TGA}`, `CodonLength=3`, `FramesPerStrand=3`.

### Formula realised correctly?
Yes — the code is a faithful, exact realisation of the validated model: fixed-stride codon scan, every ATG opens a candidate, first in-frame stop terminates, span includes the stop (`length = i + CodonLength - start`), six frames over forward + reverse complement, nucleotide inclusive minLength (`length >= minLength`). Case handled by `DnaSequence` normalisation (verified by S1). No heuristic, no approximation.

### Cross-verification recomputed vs code
The full unfiltered suite was run (`--no-build`) and passes; the dedicated fixture's exact assertions equal the independent Python values in the table above. Spot invariants (start=ATG, end∈stops, len%3==0) hold over the Rosalind set.

### Variant / delegate consistency
Single public method; no `*Fast`/delegate variant. `OrfInfo.CodonCount = Length/3` and `Length = Sequence.Length` are trivial derived properties (consistent by construction).

### Test quality audit (HARD gate)
- **Sourced, not echoes:** M2 asserts `Is.EquivalentTo` the **exact 4-protein Rosalind set**; M1/M3/M5/M6b/M10 assert exact `Is.EqualTo` Sequence/Position/Frame/IsReverseComplement. A deliberately-wrong implementation (greedy single-ATG, missing RC, off-by-one minLength, wrong stop set) fails these exact assertions — not green-washable. Verified each expected value against the independent re-derivation, not the code.
- **No green-washing:** no `GreaterThan`/`AtLeast`-style substitutes for known exacts; no widened tolerance; no skip/ignore. M5 uses unfiltered `Has.Count.EqualTo(1)` (independently confirmed the input yields exactly one ORF total).
- **Coverage:** every Stage-A branch/edge is exercised — single ORF, six-frame Rosalind, nested-shared-stop, no-in-frame-stop, reverse-strand-only, minLength exclude + at-threshold include, all three stop codons (`[TestCase]`), structural invariants INV-1/2/3, lowercase normalisation, multi-frame frame-number, too-short, and null→throw. All public behaviour covered.
- **Honest green:** full unfiltered `dotnet test` = **Failed: 0, Passed: 6619, Skipped: 0**; `dotnet build` = 0 errors. (4 NUnit2007 warnings exist only in the unrelated `ApproximateMatcher_EditDistance_Tests.cs`, a file not touched by this unit.)
- **Gate result: PASS.**

Minor, non-defect notes (no action required): the **default `minLength=100`** parameter value is never exercised (all tests pass an explicit value) and `OrfInfo.CodonCount` is not directly asserted; both are trivial and fully determined by paths that are tested. These do not weaken any sourced assertion and are not gaps in algorithm logic.

### Findings / defects
None. No code or test change made this session.

## Verdict & follow-ups
- **Stage A: PASS. Stage B: PASS. End-state: ✅ CLEAN.** Not a duplicate.
- Test-quality gate: **PASS**.
- No defects logged. Working tree changes are validation docs only.
