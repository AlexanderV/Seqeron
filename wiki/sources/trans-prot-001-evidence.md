---
type: source
title: "Evidence: TRANS-PROT-001 (Whole-sequence protein translation — framed / six-frame / ORF)"
tags: [validation, annotation]
doc_path: docs/Evidence/TRANS-PROT-001-Evidence.md
sources:
  - docs/Evidence/TRANS-PROT-001-Evidence.md
source_commit: 7122c87a12f8c52d64ed7d5f5241ff5aa19879ef
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: TRANS-PROT-001

The validation-evidence artifact for test unit **TRANS-PROT-001** — **whole-sequence
nucleotide → protein translation** (the `Translator` class, area *Translation*): framed
translation in a chosen reading frame, **six-frame** translation, and genetic-code-parameterized
ORF finding. This is the **sequence-level layer above** the per-codon table lookup validated in
[[trans-codon-001-evidence]] (TRANS-CODON-001, `GeneticCode.Translate`); the two units together
are synthesized in [[genetic-code-translation]]. One instance of the templated
[[algorithm-validation-evidence|evidence artifact]] pattern; see [[test-unit-registry]] for unit
tracking, and [[open-reading-frame-detection]] for the sibling standard-code ORF scanner.

## What this file records

- **Online sources:**
  - **Wikipedia "Translation (biology)"** — proteins built from RNA templates; codons read as
    triplets 5'→3'; AUG start (Met); UAA/UAG/UGA stop; termination via release factors.
  - **Wikipedia "Reading frame"** — three forward frames (offset 0/1/2) per strand; six frames
    total on dsDNA (3 forward + 3 reverse-complement), each oriented 5'→3'; mRNA has three
    possible frames, one translated.
  - **Wikipedia "Open reading frame"** — ORF = start→stop span; six-frame translation; a stop is
    expected ~once per 21 codons in random sequence; short ORFs (sORFs) < 100 codons can be real.
  - **NCBI Genetic Codes** (`wprintgc.cgi`) — alternative translation tables (1 Standard,
    2 Vertebrate-Mito, 3 Yeast-Mito, 11 Bacterial/Plastid) that a translator should support.
  - **UniProt P01308** — human preproinsulin (110 aa) as a verified biological oracle.
- **Datasets (documented oracles):**
  - **Insulin B chain** (UniProt P01308, positions 25–54): DNA
    `TTCGTGAACCAGCACCTGTGCGGCTCCCACCTGGTGGAAGCTCTGTACCTGGTGTGTGGGGAGCGTGGCTTCTTCTACACACCCAAGACC`
    → protein `FVNQHLCGSHLVEALYLVCGERGFFYTPKT` (30 aa) — the correctness fixture.
  - Standard codon examples (ATG/AUG→M, GCN→A, TAA/TAG/TGA→`*`).
  - Six-frame layout: forward frames +1/+2/+3, reverse-complement −1/−2/−3.
- **Corner cases / failure modes documented:** empty/null sequence; length not divisible by 3
  (trailing partial codon); no start codon → ORF finder returns empty; stop before minimum length
  → ORF filtered; multiple ORFs → multiple results; alternative genetic codes; case-insensitive
  input; DNA vs RNA (automatic T→U).

## Implementation notes (from the doc's code review)

The **`Translator`** class:
1. Accepts **DNA or RNA**, converting `T→U` automatically.
2. `frame` parameter **0, 1, or 2** (throws for other values).
3. Optional **`toFirstStop`** — terminate translation at the first in-frame stop.
4. **Six-frame** translation returns a dictionary keyed **−3…+3 excluding 0** (three forward,
   three reverse-complement).
5. **ORF finding** with configurable minimum length and both-strand search.
6. Alternative genetic codes selectable via a **`GeneticCode`** parameter — so `Translator`
   composes the [[genetic-code-translation|codon table]] rather than duplicating it.

## Verification methodology & deviations

All four implemented genetic-code tables were verified **codon-by-codon** against the NCBI
Genetic Codes page (updated 2024-09-23) using NCBI's positional `AAs`/`Starts` encoding, and the
insulin B-chain DNA→protein fixture was verified against UniProt P01308. The doc records **no
deviation** for the translation itself. Note the `Translator` ORF finder is genetic-code-
parameterized and distinct from `GenomicAnalyzer.FindOpenReadingFrames`
([[open-reading-frame-detection]], ATG-only / standard-code) and the annotation-layer
`GenomeAnalyzer.FindOrfs` — callers pick the entry point deliberately.
