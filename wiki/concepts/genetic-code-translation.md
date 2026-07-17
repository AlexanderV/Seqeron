---
type: concept
title: "Genetic code (codon → amino-acid translation tables)"
tags: [annotation, algorithm]
mcp_tools:
  - translate_dna
  - translate_rna
sources:
  - docs/algorithms/Translation/Codon_Translation.md
  - docs/algorithms/Translation/Protein_Translation.md
  - docs/algorithms/Translation/Six_Frame_Translation.md
  - docs/Evidence/TRANS-CODON-001-Evidence.md
  - docs/Evidence/TRANS-PROT-001-Evidence.md
  - docs/Evidence/TRANS-SIXFRAME-001-Evidence.md
source_commit: c9d9954b23845153e6f530e3ebb74f54dfdc8753
created: 2026-07-10
updated: 2026-07-17
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: trans-codon-001-evidence
      evidence: "Test Unit ID: TRANS-CODON-001 — Area: Translation — codon → amino-acid lookup (GeneticCode.Translate)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: trans-prot-001-evidence
      evidence: "Test Unit ID: TRANS-PROT-001 — Area: Translation — whole-sequence framed / six-frame translation + ORF (Translator)"
      confidence: high
      status: current
    - predicate: relates_to
      object: concept:test-unit-registry
      source: trans-sixframe-001-evidence
      evidence: "Test Unit ID: TRANS-SIXFRAME-001 — Area: Translation — six-frame translation + START→STOP ORF finding (Translator), Biopython/EMBOSS reference-implementation anchors"
      confidence: high
      status: current
---

# Genetic code (codon → amino-acid translation tables)

The **genetic code** is the many-to-one map from **64 nucleotide triplets (codons)** to
**20 amino acids + 3 stop signals**. Seqeron models it as `GeneticCode` in
`Seqeron.Genomics.Core` (`GeneticCode.cs`): a single-codon lookup (`Translate(codon) → char`)
plus start/stop predicates and a reverse lookup. It is the **foundational table** that every
higher codon-level operation reads — [[open-reading-frame-detection|ORF detection]] scans it
for ATG→stop spans, [[codon-optimization]] substitutes within its synonymous families, and the
[[relative-synonymous-codon-usage|RSCU]] / [[codon-adaptation-index|CAI]] /
[[effective-number-of-codons|ENC]] measures all partition codons by the amino acid this table
assigns, and [[variant-effect-annotation-vep|VEP-style variant annotation]] translates a variant's
reference and alternate codons through it to decide synonymous / missense / stop_gained / stop_lost.
Validated as test unit **TRANS-CODON-001** ([[trans-codon-001-evidence]]); the
**whole-sequence** layer above it — framed / six-frame translation and ORF finding via the
`Translator` class — is validated as **TRANS-PROT-001** ([[trans-prot-001-evidence]]) and
described below. See [[test-unit-registry]] for how the units are tracked and
[[algorithm-validation-evidence]] for the evidence-artifact pattern.

## The four supported tables

Seqeron ships **four** of the NCBI genetic-code tables (out of 33), each a static singleton
(`GeneticCode.Standard`, `.VertebrateMitochondrial`, `.YeastMitochondrial`, `.BacterialPlastid`,
plus `GetByTableNumber(int)`):

| NCBI table | Property | Key deviations from Standard |
|-----------|----------|------------------------------|
| **1** Standard | `Standard` | baseline; starts AUG (primary) + UUG/CUG (alt) |
| **2** Vertebrate Mitochondrial | `VertebrateMitochondrial` | AGA/AGG → **Stop** (not Arg); AUA → **Met** (not Ile); UGA → **Trp** (not Stop) |
| **3** Yeast Mitochondrial | `YeastMitochondrial` | CUN (CUU/CUC/CUA/CUG) → **Thr** (not Leu); AUA → **Met**; UGA → **Trp** |
| **11** Bacterial/Archaeal/Plant-Plastid | `BacterialPlastid` | same codon table as Standard, **extra start codons** (GUG/UUG/CUG/AUU/AUC/AUA) |

Codon→amino-acid mappings and the start/stop sets are taken **directly from the NCBI `AAs` and
`Starts` strings** for tables 1/2/3/11 (NCBI Genetic Codes, updated 2024-09-23); the Evidence
records **no deviation** in the mappings themselves. Stop codons translate to `'*'`.

## Degeneracy (codon redundancy)

The code is degenerate: most amino acids have several synonymous codons. Only **Met (AUG)** and
**Trp (UGG)** are single-codon families — the same fact that makes them fixed points under
[[codon-optimization]] and gives them relative adaptiveness `w ≡ 1` in [[codon-adaptation-index|CAI]].
Six-fold families (Leu, Ser, Arg) sit at the other extreme. The full family-size distribution in
the Standard table: **1-fold** Met, Trp; **2-fold** Phe, Tyr, His, Gln, Asn, Lys, Asp, Glu, Cys;
**3-fold** Ile; **4-fold** Val, Pro, Thr, Ala, Gly; **6-fold** Leu, Ser, Arg.
`GetCodonsForAminoAcid(aa)` is the reverse lookup used to verify degeneracy.

`Translate`, `IsStartCodon`, and `IsStopCodon` are all **O(1)** — a fixed `T→U`/upper-case
normalization followed by a single dictionary or set lookup per codon.

## Input contract and normalization

- **DNA or RNA** input both accepted — `Translate` normalizes `T → U` internally (so `ATG`≡`AUG`).
- **Case-insensitive** (`AUG` = `aug` = `AuG`, upper-cased internally).
- Codon **must be exactly 3 characters**; empty/`null`/wrong-length → `ArgumentException`
  (`ArgumentNullException` for null).
- `IsStartCodon` / `IsStopCodon` test membership in the table's start/stop sets; a start codon
  used mid-sequence still translates to its amino acid (**M**, not fMet — the code does not model
  formyl-methionine).

## IUPAC ambiguity → 'X' (documented as an intentional simplification)

The implementation treats **valid IUPAC ambiguity codons** (any triplet over
`ACGURYMKSWBDHVN`, e.g. `NNN`, `ANN`, `RAY`) as **untranslatable but not invalid**: `Translate`
returns `'X'` (unknown amino acid) rather than throwing. Only a triplet containing a
**non-IUPAC** character (e.g. `Z`, or `XYZ`, `12G`) raises `ArgumentException`. The
**algorithm spec** (`Codon_Translation.md`, TRANS-CODON-001, §5.2/§5.3/§6.1) documents this
`'X'` return as the **intended, "intentionally simplified"**
behaviour — ambiguity is *collapsed* to unknown rather than resolved probabilistically or
expanded across concrete codons — so the implementation **matches its algorithm spec**. This
only ever **contradicted** the older *Evidence* doc's *Documented Corner Cases* / *Known Failure
Modes* tables, which had stated that an "Unknown codon (e.g., NNN)" should yield an
`ArgumentException`; the spec doc supersedes that corner-case expectation. The NCBI mapping
tables match exactly; the whole distinction lives in the **API-contract layer** (ambiguity
handling), not in the code tables. (`GeneticCode.cs`, `Translate`.)

## Whole-sequence framed translation (the `Translator` layer)

Above the single-codon lookup sits the **`Translator`** class (test unit **TRANS-PROT-001**,
[[trans-prot-001-evidence]]), which walks a whole nucleotide sequence in triplets and emits a
protein. It **composes** this genetic-code table rather than duplicating it — the code table is
selectable via a `GeneticCode` parameter, so all four tables (1/2/3/11) and their alternative
start/stop sets flow through unchanged.

- **Reading frame.** A `frame` parameter **0, 1, or 2** offsets the first codon; other values
  throw. Each single strand thus has three frames, and reverse-complementing the strand gives
  three more — **six frames** total (Wikipedia "Reading frame").
- **Six-frame translation** returns a dictionary keyed **−3…−1, +1…+3 (0 excluded)**: three
  forward frames and three on the reverse complement, each read 5'→3'. **Reverse-frame numbering
  follows the Biopython "independent-offset" convention** (validated as **TRANS-SIXFRAME-001**,
  [[trans-sixframe-001-evidence]]): frame **−k = the reverse complement translated 5'→3' at
  offset k−1** (−1 = RC offset 0, −2 = offset 1, −3 = offset 2), with *no* codon correspondence
  to forward frame +1. This is the dominant reference-implementation behaviour
  (Biopython `six_frame_translations`) and the documented "alternative" in EMBOSS `transeq`;
  the EMBOSS *default* is instead **phase-locked** (frame −1 shares frame +1's phase). Both are
  correct biology — only the −1/−2/−3 **labels** differ. (Oracle: `ACTGG` → frame −1 = `P` here
  vs `S` under EMBOSS phase-locked; and the 39-nt six-frame table in [[trans-sixframe-001-evidence]].)
- **`toFirstStop`** optionally terminates translation at the first in-frame stop (the codon-table
  `'*'`), matching the "translate until a stop" convention; otherwise stops render as `*` and
  reading continues to the sequence end.
- **Input contract** mirrors the codon layer: **DNA or RNA** (automatic `T→U`), case-insensitive;
  a trailing partial codon (length not divisible by 3) is simply not translated.
- **ORF finding.** `Translator.FindOrfs` is a **genetic-code-parameterized** ORF scanner —
  a **START→STOP** model (EMBOSS `getorf -find 1`): a span from a start codon to the next
  in-frame stop, protein **including** the start residue and **excluding** the stop, with
  configurable minimum length (counted in **amino acids**, not nucleotides) and both-strand
  search. It is deliberately **not** contract-equivalent to the ATG-only / standard-code
  `GenomicAnalyzer.FindOpenReadingFrames` ([[open-reading-frame-detection]]) nor the
  prokaryotic-start annotation-layer `GenomeAnnotator.FindOrfs`; callers pick the entry point
  deliberately. (Oracle: `GGGATGAAACCCTAAGGG` → one ORF, start 3, end 14 inclusive, `MKP`.)

**Correctness oracle:** the human insulin **B chain** (UniProt P01308, positions 25–54) DNA
`TTCGTG…AAGACC` (90 nt) → `FVNQHLCGSHLVEALYLVCGERGFFYTPKT` (30 aa). All four tables verified
codon-by-codon against NCBI (2024-09-23); **no deviation** recorded for translation itself.

### Method contract (algorithm spec)

From the TRANS-PROT-001 algorithm spec (`Protein_Translation.md`), the `Translator` entry points
(`Seqeron.Genomics.Core/Translator.cs`) — signatures, defaults, and I/O:

| Method | Signature | Returns |
|--------|-----------|---------|
| `Translate` | `(DnaSequence \| RnaSequence \| string, GeneticCode? = Standard, int frame = 0, bool toFirstStop = false)` | `ProteinSequence` |
| `TranslateSixFrames` | `(DnaSequence, GeneticCode? = Standard)` | `IReadOnlyDictionary<int, ProteinSequence>` (keys `±1…±3`) |
| `FindOrfs` | `(DnaSequence, GeneticCode? = Standard, int minLength = 100, bool searchBothStrands = true)` | `IEnumerable<OrfResult>` (start, end, frame, protein) |

- **Null / empty contract.** The `DnaSequence` and `RnaSequence` overloads **throw
  `ArgumentNullException`** on null; the **string** overload returns an **empty `ProteinSequence`**
  for null/empty input. An invalid `frame` (not `0/1/2`) throws **`ArgumentOutOfRangeException`**.
  `TranslateSixFrames` **always returns all six keys**, even for an empty sequence.
- **`FindOrfs` defaults.** `minLength = 100` (in **amino acids**), `searchBothStrands = true`;
  reverse-strand ORFs are reported with **negative frame** values.
- **Invariants (verified in `TranslatorTests.cs`).**
  **INV-01** `Translate(seq, frame:0)` == `TranslateSixFrames(seq)[1]` (same helper, forward frame 0).
  **INV-02** protein length ≤ `floor((sequenceLength − frame) / 3)` (only complete codons consumed).
  **INV-03** `TranslateSixFrames` returns exactly the six keys `1,2,3,−1,−2,−3`.
  **INV-04** `OrfResult.NucleotideLength == EndPosition − StartPosition + 1`.
- **Accepted deviations / not-implemented (spec §5.3–§5.4).** An ORF that **runs off the end** of
  the scanned sequence without hitting a stop is still emitted if it meets `minLength` (broader than
  a strict stop-terminated definition). **Reverse-strand ORF coordinates are in the
  reverse-complement scan's coordinate system** — callers needing original-strand genomic positions
  must **remap externally**. **Nested ORFs** (opening a new ORF in the same frame before the current
  one closes) are not reported. These are the practical-scanner limits: no gene-model scoring, splice
  awareness, or regulatory context.

### Six-frame translation contract (TRANS-SIXFRAME-001 algorithm spec)

`Six_Frame_Translation.md` (TRANS-SIXFRAME-001) is the **canonical primary spec** for the
six-frame surface and its START→STOP ORF scanner, cross-anchored to Biopython
`six_frame_translations`, EMBOSS `transeq`, and EMBOSS `getorf`. Beyond the shared
`Translator` contract above, it pins:

- **Frame construction & offsets.** Forward frame `+f` = translation of the input at offset
  `f−1`; reverse frame `−f` = translation of the **reverse complement** at offset `f−1` (the
  Biopython `frames[-(i+1)] = translate(anti[i:])` convention — `−1` = RC offset 0). Each frame
  consumes only **complete** codons, so its length is exactly **⌊(len−offset)/3⌋** and any
  trailing 1–2 nt are dropped (spec INV-04).
- **Never early-terminates.** `TranslateSixFrames` renders internal stop codons as `*` and reads
  to the end of every frame — it does **not** honour `toFirstStop` (contrast the single-frame
  `Translate(..., toFirstStop: true)`). It computes the reverse complement **once** and fills all
  six frames in one shared offset loop. Empty input → six empty frames and no ORFs;
  IUPAC-ambiguous codons → `X` (inherited from `GeneticCode.Translate`).
- **`OrfResult` fields.** `StartPosition` = 0-based first base of the START codon (in the scanned
  strand's coordinates); `EndPosition` = 0-based **last** base of the STOP codon, **inclusive**;
  `Frame` = ±1…±3; `Protein` = residues with START included, STOP excluded; derived
  `NucleotideLength = EndPosition − StartPosition + 1` and `AminoAcidLength = Protein.Length`
  (spec INV-05/INV-06).
- **Complexity.** Both `TranslateSixFrames` and `FindOrfs` are **O(n) time / O(n) space** — one
  reverse complement plus six linear codon passes; the repository suffix tree is not applicable
  (sequential triplet decoding, not substring search).
- **Worked oracle.** `TranslateSixFrames("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")` →
  `frames[+1] = "MAIVMGR*KGAR*"`, `frames[−1] = "LSGTLSAAHYNGH"`;
  `FindOrfs("GGGATGAAACCCTAAGGG", minLength: 1, searchBothStrands: false)` → one ORF Start=3,
  End=14 (inclusive), Frame=+1, `MKP`.

## Scope

Four tables only (1/2/3/11) of NCBI's 33. The single-codon lookup is `GeneticCode.Translate`;
whole-sequence framed / six-frame translation and ORF finding live in `Translator` (above), and
the MCP surface exposes translation as `TranslateDna`/`TranslateRna`.

## Reference sources

**NCBI Genetic Codes** (`transl_table` 1/2/3/11 `AAs`+`Starts` strings, official spec),
**Wikipedia "Genetic code"** (64→20+3 degeneracy), **Wikipedia "Start codon"** (AUG universal,
GUG/UUG alternatives), **Wikipedia "Stop codon"** (UAA ochre / UAG amber / UGA opal), plus the
historical Nirenberg & Matthaei 1961 and Crick 1968 citations. Full record in
[[trans-codon-001-evidence]]. For the whole-sequence layer: **Wikipedia "Translation (biology)"**
(triplet reading, start/stop, release factors), **Wikipedia "Reading frame"** (three/six frames),
**Wikipedia "Open reading frame"** (six-frame, sORFs), and **UniProt P01308** (insulin B-chain
oracle) — full record in [[trans-prot-001-evidence]]. The six-frame surface is independently
anchored to reference implementations — **Biopython `six_frame_translations`** (governing
reverse-frame convention), **EMBOSS `transeq`** (frame numbering / phase-locked alternative) and
**EMBOSS `getorf`** (START→STOP ORF model) — in [[trans-sixframe-001-evidence]].
