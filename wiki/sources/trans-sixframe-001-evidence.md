---
type: source
title: "Evidence: TRANS-SIXFRAME-001 (Six-frame translation & START→STOP ORF finding)"
tags: [validation, annotation]
doc_path: docs/Evidence/TRANS-SIXFRAME-001-Evidence.md
sources:
  - docs/Evidence/TRANS-SIXFRAME-001-Evidence.md
source_commit: 950ce49428fde05020ff3c08e70ac1231947fc59
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: TRANS-SIXFRAME-001

Validation-evidence artifact for test unit **TRANS-SIXFRAME-001** — **six-frame translation**
(translate a nucleotide sequence in all 3 forward + 3 reverse-complement reading frames) and
**START→STOP ORF finding**. This is the **same whole-sequence `Translator` surface** validated
by [[trans-prot-001-evidence]] (TRANS-PROT-001), approached from the reference-implementation
angle (EMBOSS / Biopython) rather than the biology-oracle angle; both are synthesized in
[[genetic-code-translation]]. It sits above the per-codon [[trans-codon-001-evidence]]
(TRANS-CODON-001) table lookup, and is the genetic-code-parameterized sibling of the ATG-only
[[open-reading-frame-detection]] scanner. One instance of the templated
[[algorithm-validation-evidence|evidence artifact]] pattern; [[test-unit-registry]] tracks the unit.

## What this file records

- **Online sources (all reference implementations / official specs):**
  - **EMBOSS `transeq`** — six frame values `1,2,3,F / -1,-2,-3,R / 6` (exactly six frames);
    documents **both** reverse-frame numbering conventions (see below).
  - **Biopython `Bio.SeqUtils.six_frame_translations`** — the **governing algorithm**: reverse
    strand computed once (`anti = reverse_complement(seq)`); forward `frames[i+1] = translate(seq[i:…])`
    for offsets 0/1/2; reverse `frames[-(i+1)] = translate(anti[i:…])[::-1]` for offsets 0/1/2;
    dict keyed `+1,+2,+3 / -1,-2,-3`.
  - **EMBOSS `getorf`** — ORF definitions; `-find 1` = **START→STOP** (the model implemented);
    default `-minsize` 30 **nucleotides**; both strands searched, reverse ORFs marked `(REVERSE SENSE)`.
  - **NCBI Genetic Codes, transl_table=1** — Standard code `AAs`/`Starts` rows; starts TTG/CTG/ATG;
    stops TAA/TAG/TGA; initiator translated as Met by display convention.
  - **Wikipedia "Reading frame"** — six ways to read any DNA (3 each direction); reverse frames read
    5'→3' on the complementary strand = reverse-complement read forward.
- **Datasets (documented oracles):**
  - **`ACTGG` six-frame** — reverse complement `CCAGT`; **this repo returns frame −1 = `P`**
    (RC offset 0, `CCA`→P), whereas EMBOSS's phase-locked default returns `S`. Convention
    difference, not a biology error (see below).
  - **39-nt `ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG`** — full six-frame table:
    +1 `MAIVMGR*KGAR*`, +2 `WPL*WAAERVPD`, +3 `GHCNGPLKGCPI`, −1 `LSGTLSAAHYNGH`,
    −2 `YRAPFQRPITMA`, −3 `IGHPFSGPLQWP`.
  - **ORF `GGGATGAAACCCTAAGGG`** — START→STOP: StartPosition 3, EndPosition 14 (inclusive of stop),
    Frame 1, protein `MKP` (start included, stop excluded), 3 aa / 12 nt.
- **Corner cases / failure modes documented:** incomplete trailing codon ignored
  (`fragment_length = 3*((len−i)//3)`); reverse-frame numbering ambiguity (two conventions);
  no START under `-find 1` → no ORF; ORF running off the end (no downstream stop) is incomplete;
  `minLength` filtering; IUPAC-ambiguous codons absent from the 64-codon table.

## The reverse-frame numbering convention (the distinctive detail)

Two documented conventions disagree on what `−1/−2/−3` mean:

- **Phase-locked** (EMBOSS `transeq` default): frame `−1` uses the *same codon phase* as forward
  frame `1`, so its codons correspond position-by-position with frame `1`.
- **Independent offset** (Biopython / **this repo**): frame `−k` = translation of the
  **reverse complement read 5'→3' at offset `k−1`** (`−1`=offset 0, `−2`=offset 1, `−3`=offset 2).
  There is *no* codon correspondence between frame `1` and `−1`.

The repository follows the **Biopython independent-offset** convention — the dominant
reference-implementation behaviour, explicitly listed as an accepted "alternative" in the EMBOSS
transeq docs. Both are correct biology; only the −1/−2/−3 **labels** differ (hence the `ACTGG`
P-vs-S discrepancy is a labelling convention, not a bug). In Biopython the trailing `[::-1]`
only reverses the **display string** for alignment under the forward sequence — each reverse
frame's residue content is the RC translated 5'→3' at that offset.

## Assumptions recorded (source-anchored, not invented)

1. **Reverse-frame numbering = Biopython independent-offset** (`−k` = RC offset `k−1`).
2. **Stop → `*`, IUPAC-ambiguous codon → `X`** (matches the existing `GeneticCode.Translate`
   behaviour; the `X`-for-ambiguity divergence from the doc's exception promise is discussed in
   [[genetic-code-translation]]).
3. **ORF `minLength` counts amino acids** — an API-shape choice; getorf's `-minsize` is in
   nucleotides. Well-defined for any value.

## Verification methodology & deviations

Six-frame output cross-checked against the Biopython `six_frame_translations` algorithm applied to
NCBI table 1; the 39-nt table and the `MKP` ORF derived from those specs. **No deviation** beyond
the two deliberate convention choices (reverse-frame labelling; amino-acid `minLength`). The
`Translator` ORF finder here (START→STOP, both strands, alternative starts TTG/CTG per table 1)
is deliberately **not** contract-equivalent to `GenomicAnalyzer.FindOpenReadingFrames`
([[open-reading-frame-detection]], ATG-only / standard-code / nucleotide `minLength`) nor to the
annotation-layer `GenomeAnnotator.FindOrfs` — callers pick the entry point deliberately.
