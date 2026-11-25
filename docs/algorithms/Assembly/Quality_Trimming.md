# Quality Trimming (BWA / cutadapt running-sum)

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-TRIM-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-11-25 |

## 1. Overview

Quality trimming removes low-quality bases from the ends of sequencing reads before assembly or alignment. This implementation uses the running-sum method introduced by BWA and reused by cutadapt: it subtracts a quality cutoff from each base's Phred score, computes partial sums toward each read end, and cuts at the index where the partial sum is minimal [1][2]. The method is specification-driven and exact (deterministic), and is refined relative to a naive per-base cutoff in that it tolerates isolated high-quality bases inside a low-quality tail [1]. Reads shorter than a minimum length after trimming are dropped.

## 2. Scientific / Formal Basis

### 2.1 Domain Context

Each base in a FASTQ read carries a Phred quality score `Q = -10·log10(P)`, where `P` is the estimated probability that the base call is wrong [3]. In the Sanger/Illumina-1.8+ encoding, qualities are stored as ASCII characters with an offset of 33: ASCII 33–126 encode Phred 0–93 [3]. Sequencing error rates rise toward read ends, so end trimming improves downstream accuracy.

### 2.2 Core Model

For a read with Phred qualities `q_0 … q_{l-1}` and cutoff `C`, define the per-base score `s_i = q_i − C`. To trim the 3' end, compute, for each index `x`, the partial sum from `x` to the end `Σ_{i≥x} s_i`, and cut the read at the index where this sum is minimal; the bases before that index are retained [1]. Equivalently (BWA), accumulate `C − q_i` from the 3' end and keep the position of the maximum accumulated value [2]. The same procedure applied to the other end trims the 5' end [1].

BWA's reference loop [2]:

```c
for (l = p->len - 1; l >= BWA_MIN_RDLEN; --l) {
    s += trim_qual - (p->qual[l] - 33);
    if (s < 0) break;
    if (s > max) max = s, max_l = l;
}
```

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | The trimmed sequence is a contiguous substring of the input read | A single cut index is chosen per end [1][2] |
| INV-02 | Trimmed length ≤ original length | Trimming only removes bases [1][2] |
| INV-03 | Cutoff < 1 ⇒ read is not trimmed (only the min-length filter may drop it) | BWA `trim_qual < 1` guard [2] |
| INV-04 | Every output read has length ≥ minLength | Post-trim min-length filter (§4.2) |
| INV-05 | Quality decoded as Phred = ASCII − offset (33 default, 64 opt-in), malformed qualities rejected | Sanger encoding [3]; BWA `qual - 33` [2]; cutadapt `--quality-base` [1] |

### 2.5 Comparison with Related Methods

| Aspect | Running-sum (this) | Naive per-base cutoff |
|--------|--------------------|------------------------|
| Boundary rule | argmin of partial sums up to the `s < 0` early break | first base ≥ cutoff from each end |
| Isolated good base in bad tail | may be retained [1] | always cut at first good base |
| Origin | BWA / cutadapt [1][2] | ad-hoc |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| `reads` | `IReadOnlyList<(string sequence, string quality)>` | required | Reads with Phred+offset quality strings | non-null; `quality.Length == sequence.Length`; characters in [offset, `~`] |
| `minQuality` | `int` | 20 | Quality cutoff `C` for the 3' end only (3-argument overload = cutadapt `-q C` = `-q 0,C`; BWA `aln -q C`) | `< 1` disables trimming [2] |
| `cutoffFront` / `cutoffBack` | `int` | — | Separate 5'/3' cutoffs (5-argument overload; cutadapt `-q 5',3'`) | `< 1` leaves that end untrimmed |
| `minLength` | `int` | 50 | Reads shorter than this after trimming are dropped (cutadapt `-m`) | any |
| `phredOffset` | `int` | 33 | Quality ASCII offset (cutadapt `--quality-base`): 33 Sanger, 64 Illumina 1.3–1.7 | [33, 126] |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IReadOnlyList<string>` | Surviving trimmed sequences, in input order |

### 3.3 Preconditions and Validation

`reads`, and each read's sequence and quality, must be non-null (`ArgumentNullException`). Qualities are decoded by the canonical Core FASTQ decoder `QualitySequence(sequence, quality, phredOffset)` (Phred = `ASCII − offset`): a length mismatch or a character outside [offset, `~`] throws `ArgumentException` (Cock et al. 2010; Biopython `QualityIO`); before 2026-11 a short quality string threw `IndexOutOfRangeException`, a long one was silently ignored and characters below `!` decoded to negative Phred. The 5' and 3' passes run independently over the full read. A read whose trimmed length is `< minLength` (including fully removed reads, length 0) is excluded from the output.

## 4. Algorithm

### 4.1 High-Level Steps

1. Decode the qualities with `QualitySequence` (validation above). A cutoff `< 1` needs no special case: the first step of that end's scan is already negative, so the end is not trimmed (cutadapt loop; BWA's explicit `trim_qual < 1` guard) [1][2].
2. 3' end: over the full read `quality[0..n)`, scan from the end accumulating `(Phred − cutoff)`; stop as soon as the running sum becomes positive (the cutadapt/BWA `s < 0` early break, sign-flipped); set `end` to the index of the minimal partial sum reached before the break [1][2].
3. 5' end: independently over the full read `quality[0..n)`, scan from the start accumulating `(Phred − cutoff)` with the same early break; set `start` to one past the index of the minimal partial sum [1][2].
4. If `start ≥ end`, the good-quality segment is empty (cutadapt `start >= stop ⇒ (0,0)`): drop the read.
5. Otherwise, if `end − start ≥ minLength`, emit `sequence[start..end)`; else drop the read.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Phred decoding offset = 33 by default (`PhredAsciiOffset`), configurable (`phredOffset`, cutadapt `--quality-base`) [3].
- Min-length filter drops survivors below `minLength` (cutadapt `--minimum-length` semantics).
- This implementation reproduces cutadapt's `quality_trim_index` exactly, including the `s < 0` early break (the accumulator stops once the high-quality suffix/prefix outweighs the low-quality run) and the `start >= stop ⇒ empty` rule [1][2]. The early break is essential: it is what "allows some good-quality bases among the bad ones" — a high-quality base near an end protects the bases on its interior side from being trimmed by that end's pass. A naive global-minimum scan *without* the break over-trims (it would cut at the deepest minimum anywhere in the read).
- BWA additionally floors the 3' scan at `BWA_MIN_RDLEN = 35` [2]; this is available as the opt-in `bwaMinLength` (default 0 = cutadapt).

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Trim all reads | O(n · r) | O(1) extra per read | n = reads, r = read length; two linear passes per read |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.QualityTrimReads(reads, minQuality, minLength)`: 3' end only with one cutoff (cutadapt `-q C`), Phred+33 (= `QualityTrimReads(reads, 0, minQuality, minLength)`).
- `SequenceAssembler.QualityTrimReads(reads, cutoffFront, cutoffBack, minLength, phredOffset = 33, bwaMinLength = 0)`: cutadapt `-q 5',3' --quality-base N -m minLength`; `bwaMinLength: SequenceAssembler.BwaMinReadLength` = BWA `aln -q`.
- `SequenceAssembler.QualityTrimIndex(phredScores, cutoffFront, cutoffBack)`: line-by-line port of cutadapt 5.2 `qualtrim.quality_trim_index`, returns `[start, stop)`.

### 5.2 Current Behavior

Trimming runs the 5' pass and the 3' pass independently over the full read (cutadapt `quality_trim_index`), each with the `s < 0` early break; if the resulting windows cross (`start >= stop`) the read is dropped. A cutoff `< 1` leaves that end unchanged. Cross-checked against the installed cutadapt 5.2 on 5000 random reads (seed 2026; Phred+33/+64, cutoffs −3..45 per end, minLength 0..10): 5000/5000 identical. Reads with length `< minLength` after trimming (or `0`) are omitted. This is a single short scan per read, not a substring-search problem, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Running-sum cut at the index of the minimal partial sum of `(q − cutoff)`, applied independently to both ends, **including the `s < 0` early break** [1][2].
- The `start >= stop ⇒ empty good-quality segment` rule [1].
- Phred+33 decoding (`ASCII − 33`) [2][3].
- Cutoff `< 1` disables trimming [2].

**BWA 35-base floor (opt-in, A8-2).** BWA's `bwa_trim_read` scans `for (l = p->len - 1; l >= BWA_MIN_RDLEN; --l)` with
`#define BWA_MIN_RDLEN 35` (`bwtaln.h`) [2], so the first 35 bases are never examined and a read keeps at least
min(35, length) bases. `QualityTrimReads(reads, cutoffFront, cutoffBack, minLength, phredOffset, bwaMinLength)` takes the
floor as `bwaMinLength` (`SequenceAssembler.BwaMinReadLength` = 35; default 0 = cutadapt's unrestricted scan, which is why
the default trims below 35 bases); `QualityTrimIndex(..., bwaMinLength)` applies it to the 3' scan (BWA has no 5' pass). The
MCP tool `quality_trim_reads` exposes it as `bwaMinLength`. Oracle: a verbatim transcription of `bwa_trim_read`
(`tools/review/oracles/B10/wp4_oracle.py`); 1500 random reads x 7 seeds (length 0..100, cutoff 0..30) give identical kept lengths.

**Intentionally simplified:** none for the trimming rule itself (cutadapt `quality_trim_index` and BWA `bwa_trim_read` are both reproduced).

**Not implemented:**

- cutadapt `--nextseq-trim` (two-colour G handling) and `--poly-a`; other trimming models (Trimmomatic SLIDINGWINDOW / LEADING / TRAILING) live in `QualityScoreAnalyzer` (IO).

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Both-end passes run independently over the full read | Resolved | Matches cutadapt `quality_trim_index` exactly | accepted | Cutadapt runs the 5' and 3' passes independently, then drops the read if `start >= stop` [1] |
| 2 | `minLength` post-trim filter | Assumption | Drops short survivors | accepted | cutadapt min-length semantics |
| 3 | Single cutoff (3-argument overload) trims the 3' end only | Resolved (audit A8-1) | Equals cutadapt `-q C` = `-q 0,C` ("If only one cutoff is given, it is used for the 3' end only") and BWA `aln -q C`; both ends need both cutoffs (`cutoffFront`, `cutoffBack`) | resolved | cutadapt 5.2 `quality_trim_index(q, 0, C)` oracle `tools/review/oracles/B10/trim_cutadapt_oracle.py`, 20000/20000 identical |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| All-high-quality read | unchanged | partial-sum minimum at the end [1] |
| All-low-quality read | dropped (length 0) | both end passes consume the read [1] |
| Cutoff ≤ 0 | unchanged (subject to min-length) | BWA `trim_qual < 1` [2] |
| Trimmed length < minLength | read dropped | min-length filter (INV-04) |
| Empty `reads` | empty result | nothing to process |
| `reads` null | `ArgumentNullException` | input validation |
| quality length ≠ sequence length, or character outside [offset, `~`] | `ArgumentException` | malformed FASTQ (Cock et al. 2010) via `QualitySequence` |

### 6.2 Limitations

The encoding is not auto-detected: Phred+64 input must be passed with `phredOffset: 64` (as cutadapt `--quality-base 64`). Does not perform adapter removal or per-window (sliding-window) trimming.

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
var reads = new[] { ("ACGTACGTAC", "KI;<)(,%#$") }; // Phred 42,40,26,27,8,7,11,4,2,3
IReadOnlyList<string> kept = SequenceAssembler.QualityTrimReads(reads, minQuality: 10, minLength: 1);
// kept[0] == "ACGT"  (trimmed to the first four bases)
```

**Numerical walk-through (cutadapt example [1]):**

Qualities `42,40,26,27,8,7,11,4,2,3`, cutoff 10. After subtraction: `32,30,16,17,-2,-3,1,-6,-8,-7`. Partial sums from the end: `(70),(38),8,-8,-25,-23,-20,-21,-15,-7`. Minimum `-25` at index 4 → keep the first four bases [1].

### 7.3 Related Tests, Evidence, or Documents

- Tests: [SequenceAssembler_QualityTrimReads_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/Unit/Alignment/SequenceAssembler_QualityTrimReads_Tests.cs) — covers `INV-01`–`INV-05`
- Evidence: [ASSEMBLY-TRIM-001-Evidence.md](../../../docs/Evidence/ASSEMBLY-TRIM-001-Evidence.md)

## 8. References

1. Cutadapt project. 2024. Source `src/cutadapt/qualtrim.pyx` v5.2 (https://raw.githubusercontent.com/marcelm/cutadapt/v5.2/src/cutadapt/qualtrim.pyx) and `cli.py` `parse_cutoffs` (single value → (0, C)).  Algorithm details — Cutadapt documentation (quality trimming). https://cutadapt.readthedocs.io/en/stable/algorithms.html
2. Li, H. BWA source: `bwa_trim_read` (bwaseqio.c) and `BWA_MIN_RDLEN` (bwtaln.h). https://github.com/lh3/bwa/blob/master/bwaseqio.c , https://github.com/lh3/bwa/blob/master/bwtaln.h
3. Cock, P.J.A., Fields, C.J., Goto, N., Heuer, M.L., Rice, P.M. 2010. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Research 38(6):1767–1771. https://doi.org/10.1093/nar/gkp1137
