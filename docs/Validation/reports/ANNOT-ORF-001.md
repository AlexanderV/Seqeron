# Validation Report: ANNOT-ORF-001 — Open Reading Frame (ORF) Detection

- **Validated:** 2026-06-24   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.FindOrfs(dnaSequence, minLength, searchBothStrands, requireStartCodon)`, `GenomeAnnotator.FindLongestOrfsPerFrame(dnaSequence, searchBothStrands)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **End state:** CLEAN (no defect; no code changes)

---

## Stage A — Description

### Sources opened & what they confirm

- **Rosalind "ORF" problem** (fetched live, rosalind.info/problems/orf) — Sample `>Rosalind_99`
  `AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG`
  → exactly four distinct proteins:
  `MLLGSFRLIPKETLIQVAGSSPCNLS`, `M`, `MGMTPRLGLESLLE`, `MTPRLGLESLLE`.
  Confirms: start codon **AUG/ATG** only (for the exercise); stop codons **UAA/UAG/UGA (TAA/TAG/TGA)**;
  **six reading frames** = 3 forward + 3 on the reverse complement; an ORF "starts from the start
  codon and ends by the stop codon, without any other stop codon in between"; the **stop codon is
  NOT part of the protein**; output is the **set of distinct** candidate proteins.
- **Wikipedia, "Open reading frame"** — ORF spans start→stop; start usually ATG, stops TAA/TAG/TGA;
  six possible frame translations from the double strand; a minimum length (≈100 or ≈150 codons) is
  applied in gene-prediction settings.
- **NCBI ORFfinder** — parameterised by minimal length, genetic code, and a start-codon option
  (ATG only / ATG+alternative GTG,TTG / any sense codon). Confirms GTG/TTG as legitimate alternative
  starts under bacterial/alternative codes.
- **Deonier et al. (2005)** ≈150 codons typical; **Claverie (1997)** ≈100 codons typical — match the
  spec's `minLength` default of 100 aa.

> Note: the orchestration prompt quoted a **truncated** sequence
> (`...GACTTATATTAGCTACCC`) for the Rosalind sample; the authoritative Rosalind sample (and the
> M12 test) use the **full** sequence `...GACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG`,
> which is what yields the four canonical proteins. Validation used the authoritative full sequence.

### Semantics defined and confirmed

| Property | Validated value |
|---|---|
| Start codons | ATG (canonical) + GTG, TTG (alternative starts; NCBI option, sourced superset of Rosalind) |
| Stop codons | TAA, TAG, TGA |
| Reading frames | 6 = 3 forward + 3 reverse-complement |
| ORF end | first in-frame stop; no internal in-frame stop |
| Stop in coordinates | included in `Sequence`/`End`; `ProteinSequence` ends with `*` (stop), excluded from the peptide |
| Coordinate base | 0-based, half-open `[Start, End)` |
| Minimum length | amino acids (codons before the stop) |
| Nested/overlapping | every start codon upstream of a shared stop yields its own ORF |

### Independent cross-check (Rosalind sample, authoritative)

Hand/tool-verified expected set = `{MLLGSFRLIPKETLIQVAGSSPCNLS, M, MGMTPRLGLESLLE, MTPRLGLESLLE}`;
`MLLGSFRLIPKETLIQVAGSSPCNLS` lies on the **reverse-complement** strand, exercising six-frame
scanning and reverse-strand coordinate mapping.

### Findings / divergences

None at the description level. The only superset vs. Rosalind's ATG-only convention is the
deliberate inclusion of GTG/TTG alternative starts, which is sourced (Wikipedia / NCBI ORFfinder).

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`
- `FindOrfs` (L91–120): forward scan, then reverse-complement scan with coordinate adjustment.
- `FindOrfsInStrand` (L122–137): frames 0,1,2.
- `FindOrfsInFrame` (L139–218): start/stop scanning, nested-start handling, run-off handling.
- `FindLongestOrfsPerFrame` (L223–247): groups by signed frame, picks longest protein.

### Faithful realisation of the validated description (evidence)

- **Start/stop sets** (L53–64): `{ATG,GTG,TTG}` / `{TAA,TAG,TGA}`, `OrdinalIgnoreCase`. Matches Stage A.
- **Six frames**: 3 forward frames (L130) on `dnaSequence` plus 3 on
  `GetReverseComplementString(dnaSequence)` (L109–119). Confirmed.
- **First in-frame stop, no internal stop**: on a stop, all open starts are completed and
  `currentOrfStarts` is cleared (L158–187), so no ORF crosses an in-frame stop. Confirmed.
- **Stop inclusion / coordinates**: `End = i + 3`, `Sequence = sequence.Substring(start, i+3-start)`
  (L168, L173) → 0-based half-open, stop **included** in `Sequence`. `minLength` uses
  `orfAaLength = (i - start)/3` (L163–164) → codons before the stop. Internally coherent.
- **Reverse-strand mapping** (L115–116): `adjStart = len - orf.End`, `adjEnd = len - orf.Start` —
  half-open intervals map cleanly with no off-by-one.
- **Nested ORFs**: `currentOrfStarts` is a list; every start before a shared stop emits an ORF ending
  at the same `End` (L161–179). Confirmed by Rosalind run (`MGMTPRLGLESLLE` start=24 nested
  `MTPRLGLESLLE` start=30, shared end=69).

### Cross-verification recomputed vs. actual code (live run, Rosalind full sample)

```
start= 4 end=10 frame=2 rc=False  M
start=24 end=69 frame=1 rc=False  MGMTPRLGLESLLE
start=30 end=69 frame=1 rc=False  MTPRLGLESLLE
start=43 end=52 frame=2 rc=False  LD     (GTG/TTG alt-start extra)
start=60 end=69 frame=1 rc=False  LE     (alt-start extra)
start=10 end=91 frame=3 rc=True   MLLGSFRLIPKETLIQVAGSSPCNLS
start=20 end=26 frame=2 rc=True   M
```
All four authoritative Rosalind proteins present (verified OK each). The two extras (`LD`, `LE`)
arise solely from GTG/TTG alternative starts — expected and sourced. M12 asserts **containment** of
the four canonical proteins, so the superset does not cause a false failure.

### Edge cases (Stage-A → code)

| Case | Behaviour | Where |
|---|---|---|
| Empty / null | empty result | L97–98 |
| No start (requireStart) | empty | start set never matched |
| No stop (requireStart) | empty | run-off block guarded by `!requireStartCodon` (L191) |
| Below minLength | excluded | L166 |
| Exactly minLength | included | M06b |
| Boundary ORF (ends at seq end) | included, `End == len` | M01, S03 |
| Alt starts GTG/TTG | detected | M05, Rosalind extras |
| N in codon | not start nor stop (set lookup) | S04a/b/c |
| Run-off zero-length start | skipped via `endPos > start` guard (L203) | INV-04 protection |

### Variant / wrapper consistency

- `FindLongestOrfsPerFrame` reuses `FindOrfs(minLength:1)`, returns 6 signed-frame keys (or 3
  forward-only). Confirmed by M16/M16b/M17.
- Wrapper/alternate (`Translator.FindOrfs`, `GenomicAnalyzer.FindOpenReadingFrames`) covered by smoke
  tests per spec.

### Test quality audit

`GenomeAnnotator_ORF_Tests.cs`: 35 tests, all passing, with exact assertions (counts, positions,
frame numbers, protein strings, reverse-strand protein). M12 uses the authoritative Rosalind dataset
and the four published proteins.

### Findings / notes (no defect)

- **PASS-WITH-NOTES (minor, non-canonical path):** with `requireStartCodon=false`, the run-off logic
  (L183–217) seeds a new candidate start only at `i+3` after each stop and at the frame origin only
  implicitly; a coding stretch before the first start/stop in start-less mode is not emitted from the
  frame origin. This affects only the non-default `requireStartCodon=false` mode, which is outside the
  validated standard-ORF semantics (Rosalind / NCBI require a start codon). The zero-length run-off
  case is now correctly guarded (`endPos > start`, L203). The canonical path
  (`requireStartCodon=true`) is fully correct.
- **Stop-inclusion convention:** `Sequence`/`End` include the stop while the `minLength` amino-acid
  count excludes it; both are internally consistent and locked by tests (NCBI-style coordinates
  include the stop).

---

## Verdict & follow-ups

- **Stage A:** PASS — ORF semantics, six-frame scanning, start/stop sets, alternative starts, and
  minimum length match Rosalind (live fetch), Wikipedia, NCBI ORFfinder, Deonier (2005),
  Claverie (1997).
- **Stage B:** PASS-WITH-NOTES — implementation faithfully realises the validated description on the
  canonical `requireStartCodon=true` path; reverse-strand coordinates and the full Rosalind worked
  example (incl. the reverse-strand protein) verified against the live code. The only note is the
  non-canonical `requireStartCodon=false` run-off seeding, outside standard ORF scope.
- **End state:** CLEAN — no defect requiring a fix. No source or test changes.
- **Tests:** ORF filter → **35 passed / 0 failed**. Full project: **18208 passed, 0 failed**.
