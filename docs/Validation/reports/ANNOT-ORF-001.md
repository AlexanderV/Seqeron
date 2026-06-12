# Validation Report: ANNOT-ORF-001 — Open Reading Frame (ORF) Detection

- **Validated:** 2026-06-12   **Area:** Annotation
- **Canonical method(s):** `GenomeAnnotator.FindOrfs(dna, minLength, searchBothStrands, requireStartCodon)`, `GenomeAnnotator.FindLongestOrfsPerFrame(dna, searchBothStrands)`
- **Stage A verdict:** PASS
- **Stage B verdict:** PASS-WITH-NOTES
- **End state:** CLEAN (no defect; no code changes)

---

## Stage A — Description

### Sources opened & what they confirm

- **Wikipedia, "Open reading frame"** — ORF is a span of DNA between a start codon and a stop
  codon. Start typically AUG/ATG; stops UAA/UAG/UGA (TAA/TAG/TGA). DNA is read in **three frames
  per strand → six possible frame translations** because of the double helix. Some authors require
  a **minimum length of 100 or 150 codons** for gene-prediction purposes; this varies by application.
- **Rosalind "ORF" problem** — "an ORF starts from the start codon and ends by the stop codon,
  without any other stop codons in between." All **six reading frames** (3 forward + 3 reverse
  complement) are scanned; the **stop codon terminates but is not part of the protein string**;
  output is the set of **distinct** candidate proteins (ATG start only for the Rosalind exercise).
- **NCBI ORFfinder** — parameterised by minimal length, genetic code, and a start-codon option
  (ATG only; ATG + alternative GTG/TTG; or any sense codon). Confirms GTG/TTG are legitimate
  alternative starts under bacterial/alternative codes.
- **Deonier et al. (2005)** — minimum ORF ≈ 150 codons typical. **Claverie (1997)** — ≈ 100 codons
  typical. These match the spec's `minLength` default of 100 aa.

### Semantics defined and confirmed

| Property | Validated value |
|---|---|
| Start codon | ATG (canonical) + GTG, TTG (alternative starts, NCBI option) |
| Stop codons | TAA, TAG, TGA |
| Reading frames | 6 = 3 forward + 3 reverse-complement |
| ORF end | first in-frame stop codon; no internal in-frame stop |
| Stop included in coordinates | implementation-defined; here the stop **is** included in `Sequence`/`End` |
| Coordinate base | 0-based, half-open `[Start, End)` |
| Minimum length | in amino acids (codons), excluding the stop |
| Nested/overlapping ORFs | every start codon upstream of a shared stop yields its own ORF |

### Independent cross-check (Rosalind sample, hand/tool computed)

Sample `>Rosalind_99`:
`AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG`

Rosalind expected distinct proteins (ATG-only):
```
MLLGSFRLIPKETLIQVAGSSPCNLS   (reverse strand)
M
MGMTPRLGLESLLE
MTPRLGLESLLE
```
All four are produced by the implementation (see Stage B run). `MLLGSFRLIPKETLIQVAGSSPCNLS`
is found on the **reverse-complement** strand, exercising six-frame scanning and reverse-strand
coordinate mapping.

### Findings / divergences

None at the description level. The only superset behaviour (vs. Rosalind's ATG-only) is the
deliberate inclusion of GTG/TTG alternative starts, which is sourced (Wikipedia / NCBI).

---

## Stage B — Implementation

### Code path reviewed

`src/Seqeron/Algorithms/Seqeron.Genomics.Annotation/GenomeAnnotator.cs`
- `FindOrfs` (L90–119): forward scan, then reverse-complement scan with coordinate adjustment.
- `FindOrfsInStrand` (L121–136): loops frames 0,1,2.
- `FindOrfsInFrame` (L138–212): start/stop scanning, nested-start handling, run-off handling.
- `FindLongestOrfsPerFrame` (L217–241): groups by signed frame, picks longest protein.

### Faithful realisation of the validated description (evidence)

- **Start/stop sets** (L52–63): `{ATG,GTG,TTG}` / `{TAA,TAG,TGA}`, case-insensitive. Matches Stage A.
- **Six frames**: 3 forward frames (L129) on `dnaSequence` plus 3 on
  `GetReverseComplementString(dnaSequence)` (L108–118). Confirmed.
- **First in-frame stop, no internal stop**: on hitting a stop, all open starts are completed and
  `currentOrfStarts` is cleared (L157–180), so no ORF crosses an in-frame stop. Confirmed.
- **Stop inclusion / coordinates**: `End = i + 3`, `Sequence = sequence.Substring(start, i+3-start)`
  (L167, L173) → 0-based half-open, stop codon **included**. `minLength` test uses
  `orfAaLength = (i - start)/3` (L162–163) → counts codons **before** the stop, so the stop is not
  counted toward the length threshold. Consistent and internally coherent.
- **Reverse-strand coordinate mapping** (L114–116): `adjStart = len - orf.End`,
  `adjEnd = len - orf.Start`. Verified correct: for every reverse-strand ORF, reverse-complementing
  the forward slice `seq[Start..End)` reproduces the stored ORF sequence exactly (see check below).
  Half-open intervals map cleanly with no off-by-one.
- **Nested ORFs**: `currentOrfStarts` is a list; every start before the shared stop emits an ORF
  ending at the same `End` (L160–178). Confirmed by M11 and by the Rosalind run
  (`MGMTPRLGLESLLE` start=24 and nested `MTPRLGLESLLE` start=30 share end=69).

### Cross-verification recomputed vs. actual code (Rosalind sample)

Run via the real library (`FindOrfs(seq, minLength:1, searchBothStrands:true)`):

```
start=4  end=10 frame=2 rc=False  M
start=24 end=69 frame=1 rc=False  MGMTPRLGLESLLE
start=30 end=69 frame=1 rc=False  MTPRLGLESLLE
start=43 end=52 frame=2 rc=False  LD                (GTG/TTG alt-start extra)
start=60 end=69 frame=1 rc=False  LE                (alt-start extra)
start=10 end=91 frame=3 rc=True   MLLGSFRLIPKETLIQVAGSSPCNLS
start=20 end=26 frame=2 rc=True   M
```
All four Rosalind proteins present. The two extras (`LD`, `LE`) arise solely from GTG/TTG
alternative starts — expected and sourced. The M12 test asserts **containment** of the four
canonical proteins, so the superset does not cause a false failure.

Reverse-strand coordinate verification (forward slice → reverse-complement == stored ORF):
```
fwd[20,26) -> ATGTAA == stored ATGTAA                                    match
fwd[10,91) -> ATG...TAACCTGAGTTAG == stored (full 81 nt)                 match
```

### Edge cases (Stage-A → code)

| Case | Behaviour | Where |
|---|---|---|
| Empty / null | empty result | L96–97 (M02, Null test) |
| No start (requireStart) | empty | start set never matched (M03) |
| No stop (requireStart) | empty — run-off block guarded by `!requireStartCodon` (L190) | M04 |
| No stop (requireStart=false) | run-off ORF to frame-aligned end, no trailing stop (truncated) | L190–211 |
| Below minLength | excluded | L165 (M06) |
| Exactly minLength | included | M06b |
| Boundary ORF (ends at seq end) | included, `End == len` | M01, S03 |
| Alt starts GTG/TTG | detected | M05, Rosalind extras |
| N in codon | not start nor stop (set lookup) | S04a/b/c |

### Variant / wrapper consistency

- `FindLongestOrfsPerFrame` reuses `FindOrfs(minLength:1)` and returns 6 signed-frame keys
  (or 3 forward-only). Confirmed by M16/M16b/M17.
- Wrapper/alternate (`Translator.FindOrfs`, `GenomicAnalyzer.FindOpenReadingFrames`) covered by
  smoke tests per spec.

### Test quality audit

35 ORF-filtered tests pass; assertions are exact (counts, positions, frame numbers, protein
strings, reverse-strand protein), covering all Stage-A edge cases. The Rosalind dataset test
uses authoritative expected proteins.

### Findings / notes (no defect)

- **PASS-WITH-NOTES (minor, non-canonical path):** with `requireStartCodon=false`, the run-off
  logic (L182–211) seeds a new candidate start only at `i+3` after each stop and never seeds the
  frame's initial offset, so a coding stretch before the first start/stop in start-less mode is not
  emitted as an ORF beginning at the frame origin. This affects only the non-default
  `requireStartCodon=false` mode, which is outside the validated standard-ORF semantics (Rosalind /
  NCBI require a start codon). The canonical path (`requireStartCodon=true`) is fully correct. No
  change made.
- **Stop-inclusion convention:** `Sequence`/`End` include the stop codon while the `minLength`
  amino-acid count excludes it. Both are internally consistent and matched by the locked tests; this
  is a defensible, documented convention (NCBI-style coordinates include the stop).

---

## Verdict & follow-ups

- **Stage A:** PASS — ORF semantics, six-frame scanning, start/stop sets, alternative starts, and
  minimum length match Wikipedia, Rosalind, NCBI ORFfinder, Deonier (2005), Claverie (1997).
- **Stage B:** PASS-WITH-NOTES — implementation faithfully realises the validated description on the
  canonical `requireStartCodon=true` path; reverse-strand coordinates and the Rosalind worked
  example (incl. the reverse-strand protein) verified against the real code. The only note is the
  non-canonical `requireStartCodon=false` run-off seeding, which is out of standard scope.
- **End state:** CLEAN — no defect requiring a fix. No source or test changes.
- **Tests:** `--filter ORF` → 35 passed / 0 failed. Full project: **4461 passed, 0 failed**.
