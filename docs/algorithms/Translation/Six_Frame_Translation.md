# Six-Frame Translation and ORF Finding

| Field | Value |
|-------|-------|
| Algorithm Group | Translation |
| Test Unit ID | TRANS-SIXFRAME-001 |
| Related Projects | Seqeron.Genomics.Core |
| Implementation Status | Production |
| Last Reviewed | 2026-09-28 |

## 1. Overview

Six-frame translation produces the six possible protein readings of a double-stranded
nucleotide sequence: three forward frames (reading offsets 0, 1, 2 of the given strand)
and three reverse frames (offsets 0, 1, 2 of its reverse complement) [1][5]. Open
Reading Frame (ORF) finding then locates coding-candidate regions that begin at a START
codon and end at a STOP codon [4]. Both operations are exact and specification-driven:
they apply a genetic-code table (NCBI translation tables) deterministically, with no
heuristic or probabilistic step. They are used in gene prediction, unannotated-sequence
exploration, and proteogenomics.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

A codon is a triplet of nucleotides read 5'→3'; the genetic code maps each of the 64
codons to one amino acid or a stop signal [3]. Because reading can start at any of three
offsets and the antiparallel complementary strand provides three more, a duplex sequence
has exactly six reading frames [5].

### 2.2 Core Model

Let `s` be the input sequence (T→U normalised) and `G` the genetic-code map.
Forward frame `f ∈ {1,2,3}` is the concatenation
`G(s[o]s[o+1]s[o+2]) · G(s[o+3]…) · …` with offset `o = f-1`, consuming only complete
codons [2]. Reverse frame `-f` is the same construction applied to `revcomp(s)` with
offset `o = f-1` [2]; this is Biopython's `frames[-(i+1)] = translate(anti[i:])`
convention (see §5.4). An ORF (getorf `-find 1`) is a maximal region beginning at a
START codon and ending at the next in-frame STOP codon, both selected from the active
genetic code [4]. The initiator residue is reported as Met whatever the start codon
(getorf `-methionine`, default Y [6]; Biopython `translate(cds=True)` [2]; NCBI [3]). The Standard code (NCBI table 1) has start codons TTG, CTG, ATG and
stop codons TAA, TAG, TGA [3].

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | `TranslateSixFrames` yields exactly 6 frames keyed {+1,+2,+3,−1,−2,−3} | three forward + three reverse frames of a duplex [1][5] |
| INV-02 | Forward frame `+f` = translation at offset `f−1` of the input | direct codon reading [2] |
| INV-03 | Reverse frame `−f` = translation at offset `f−1` of the reverse complement | Biopython reverse-frame loop [2] |
| INV-04 | Each frame length = ⌊(len−offset)/3⌋; trailing partial codon ignored | only complete codons consumed [2] |
| INV-05 | Every ORF starts at a START codon and its protein starts with `M`; if terminated, EndPosition is the STOP's last base (inclusive, INSDC CDS convention [7]) and the protein excludes the STOP; if open, EndPosition is the last base of the last complete codon | getorf START→STOP model [4][6]; INSDC [7] |
| INV-06 | `NucleotideLength = EndPosition − StartPosition + 1` = 3·(AminoAcidLength + 1) if terminated, 3·AminoAcidLength if open; `AminoAcidLength = Protein.Length` | inclusive coordinates [6][7] |

### 2.5 Comparison with Related Methods

| Aspect | This module (Biopython convention) | EMBOSS transeq default |
|--------|-----------------------------------|------------------------|
| Reverse-frame label | frame −1 = reverse complement at offset 0 | frame −1 = phase-locked to forward frame 1 |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| dna | DnaSequence | required | Sequence to translate / scan | non-null |
| geneticCode | GeneticCode? | Standard (table 1) | Codon→amino-acid table | NCBI translation table; `FindOrfs` rejects the dual-coding-stop tables 27, 28, 31 (`ArgumentException`) |
| minLength | int | 100 | Minimum ORF length, in **amino acids** | ≥ 0 |
| searchBothStrands | bool | true | Also scan the reverse complement for ORFs | — |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (six frames) | IReadOnlyDictionary<int, ProteinSequence> | Keys +1,+2,+3,−1,−2,−3 → protein per frame |
| OrfResult.StartPosition | int | 0-based first base of the START codon (in the scanned strand's coordinates) |
| OrfResult.EndPosition | int | 0-based, inclusive: last base of the STOP codon (INSDC CDS convention [7]; getorf prints this − 3); for an open ORF the last base of the last complete codon |
| OrfResult.Frame | int | +1..+3 forward, −1..−3 reverse |
| OrfResult.Protein | ProteinSequence | Translated residues; initiator reported as `M` (getorf `-methionine` [6]); STOP excluded |

### 3.3 Preconditions and Validation

Null `dna` throws `ArgumentNullException`. Input is upper-cased and T→U normalised
before codon lookup; indexing is 0-based; ORF EndPosition is inclusive. An empty
sequence yields six empty frames and no ORFs. IUPAC-ambiguous codons are resolved as in
Biopython (inherited from `GeneticCode.Translate`: e.g. `GCN`→`A`, `TAR`→`*`, `RAY`→`B`, `NNN`/`TAN`→`X`).

## 4. Algorithm

### 4.1 High-Level Steps

1. Validate input; resolve the genetic code (default Standard).
2. Compute the reverse complement once.
3. For offsets 0,1,2: translate the forward strand (→ +1,+2,+3) and the reverse
   complement (→ −1,−2,−3), consuming only complete codons.
4. For ORFs: in each frame, on entering a START codon begin accumulating residues;
   (the first residue is written as `M`); on the next in-frame STOP, emit an ORF if its
   protein length ≥ `minLength`; if the strand ends mid-ORF, emit the open ORF (ending at
   the last complete codon) if long enough. Starts inside an open ORF are not reported
   separately (getorf `ORF[frame]` flag [6]).
5. If `searchBothStrands`, repeat ORF scanning on the reverse complement with negative
   frame labels.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

Standard genetic code (NCBI table 1): start codons {TTG, CTG, ATG}; stop codons
{TAA, TAG, TGA} [3]. Codon length = 3; three reading frames per strand [5].

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| TranslateSixFrames | O(n) | O(n) | n = sequence length; one revcomp + six linear passes |
| FindOrfs | O(n) | O(n) | linear scan per frame on each searched strand |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [Translator.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Core/Translator.cs)

- `Translator.TranslateSixFrames(DnaSequence, GeneticCode?)`: six-frame translation.
- `Translator.FindOrfs(DnaSequence, GeneticCode?, int minLength, bool searchBothStrands)`: ORF finding.
- `OrfResult` (record struct): ORF coordinates, frame, protein, derived lengths.

### 5.2 Current Behavior

Translation uses a direct O(n) codon scan; the repository suffix tree is not applicable
(this is sequential triplet decoding, not substring search). Frames are computed in a
single offset loop with one shared reverse complement. `TranslateSixFrames` renders
internal stop codons as `*` (it does not terminate early), unlike `Translate(...,
toFirstStop: true)`.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Six frames keyed +1..+3 / −1..−3, forward offsets 0/1/2 and reverse-complement offsets 0/1/2 (Biopython six_frame_translations) [2].
- Trailing partial codon ignored (Biopython `fragment_length` truncation) [2].
- Standard code start/stop codons {TTG,CTG,ATG} / {TAA,TAG,TGA} (NCBI table 1) [3].
- ORF = START→STOP region (EMBOSS getorf `-find 1`, `getorf_FindORFs`) [4][6]; initiator → `M` (`-methionine` default Y) [6]; open ORF ends at last complete codon (`WriteORF(start, pos+2)`) [6].
- End coordinate includes the STOP codon (INSDC feature table CDS: "location includes stop codon") [7]; getorf's printed end excludes it (`WriteORF(start, pos-1)`) — a coordinate-convention difference of exactly 3.
- Tables 27/28/31: every stop codon is dual-coding (Biopython `CodonTable`: `stop_codons ⊂ forward_table`), so `FindOrfs` throws `ArgumentException` — as Biopython refuses `to_stop=True` for these tables [2].

**Intentionally simplified:**

- ORF `minLength` is measured in amino acids; **consequence:** users converting from EMBOSS getorf `-minsize` (nucleotides) must divide by 3.

**Not implemented:**

- EMBOSS phase-locked reverse-frame numbering; **users should rely on:** the Biopython convention used here (the EMBOSS-documented alternative) — no separate phase-locked mode is provided.
- getorf modes other than START→STOP (`-find` 0,2,3,4,5); **users should rely on:** `TranslateSixFrames` plus manual stop-to-stop scanning if needed.

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Reverse-frame numbering = Biopython independent-offset (not EMBOSS phase-locked) | Assumption | −1/−2/−3 labels differ from EMBOSS default | accepted | Documented alternative in EMBOSS transeq [1]; see INV-03 |
| 2 | `minLength` in amino acids, not nucleotides | Deviation | Different parameter unit than getorf | accepted | §3.1 |
| 3 | End coordinate includes the STOP (INSDC) | Convention | getorf's printed end = EndPosition − 3 | accepted | [7] |
| 4 | Reverse-strand ORF coordinates in reverse-complement coordinates | Convention | getorf prints forward-strand coordinates (`len−s+1`) | accepted | §6.2 |
| 5 | A START in the very last complete codon of a frame (no STOP) is reported as a 1-aa open ORF when `minLength ≤ 1` | Divergence | getorf never writes it (branch-order artefact of `getorf_FindORFs`: the START branch opens the ORF and the loop ends) | accepted | unobservable for `minLength ≥ 2` (getorf default minsize 30 nt = 10 aa) |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Null `dna` | `ArgumentNullException` | contract |
| Empty sequence | six empty frames; no ORFs | no complete codons [2] |
| Length not multiple of 3 | trailing 1–2 nt ignored | Biopython truncation [2] |
| No START codon | no ORF | getorf START→STOP [4] |
| ORF runs to sequence end without STOP | open ORF emitted if ≥ minLength; ends at last complete codon | getorf `WriteORF(start, pos+2)` [6] |
| Alternative start codon (TTG, CTG, GTG…) | initiator reported as `M` | getorf `-methionine` [6]; Biopython `cds=True` [2] |
| Genetic code 27, 28 or 31 in `FindOrfs` | `ArgumentException` | no unambiguous stop codon; Biopython `to_stop` refusal [2] |

### 6.2 Limitations

Reverse-frame labels follow the Biopython convention and will not match EMBOSS transeq's
phase-locked default. ORF finding implements only the START→STOP model (getorf `-find 1`)
and does not provide stop-to-stop or flanking-nucleotide modes. ORF coordinates for
reverse-strand hits are expressed in the reverse complement's coordinate frame.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var dna = new DnaSequence("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG");
var frames = Translator.TranslateSixFrames(dna);
// frames[1]  == "MAIVMGR*KGAR*"
// frames[-1] == "LSGTLSAAHYNGH"

var orfs = Translator.FindOrfs(new DnaSequence("GGGATGAAACCCTAAGGG"),
                               minLength: 1, searchBothStrands: false);
// one ORF: Start=3, End=14, Frame=1, Protein="MKP"
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [Translator_SixFrames_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Core/Translator_SixFrames_Tests.cs) — covers `INV-01`..`INV-06`
- Evidence: [TRANS-SIXFRAME-001-Evidence.md](../../../docs/Evidence/TRANS-SIXFRAME-001-Evidence.md)

## 8. References

1. Rice P, Longden I, Bleasby A. 2000. EMBOSS: transeq application documentation. EMBOSS. https://emboss.sourceforge.net/apps/cvs/emboss/apps/transeq.html
2. Cock PJA et al. 2009. Biopython. Bioinformatics 25(11):1422–1423. Source Bio/SeqUtils/__init__.py (six_frame_translations). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
3. NCBI. The Genetic Codes — Standard Code (transl_table=1). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
4. Rice P, Longden I, Bleasby A. 2000. EMBOSS: getorf application documentation. EMBOSS. https://emboss.sourceforge.net/apps/cvs/emboss/apps/getorf.html
5. Wikipedia contributors. Reading frame (cites Lodish 2007; Pierce 2012). https://en.wikipedia.org/wiki/Reading_frame
6. EMBOSS getorf source: `emboss/getorf.c` (`getorf_FindORFs`, `getorf_WriteORF`) and `emboss/acd/getorf.acd` (`-methionine` default Y). https://raw.githubusercontent.com/kimrutherford/EMBOSS/master/emboss/getorf.c
7. INSDC. The DDBJ/ENA/GenBank Feature Table Definition — CDS: "sequence of nucleotides that corresponds with the sequence of amino acids in a protein (location includes stop codon)". https://www.insdc.org/submitting-standards/feature-table/
