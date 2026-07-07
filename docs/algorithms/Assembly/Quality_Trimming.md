# Quality Trimming (BWA / cutadapt running-sum)

| Field | Value |
|-------|-------|
| Algorithm Group | Assembly |
| Test Unit ID | ASSEMBLY-TRIM-001 |
| Related Projects | Seqeron.Genomics.Alignment |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

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
| INV-05 | Quality decoded as Phred = ASCII − 33 | Sanger encoding [3]; BWA `qual - 33` [2] |

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
| `reads` | `IReadOnlyList<(string sequence, string quality)>` | required | Reads with Phred+33 quality strings | non-null; `quality.Length == sequence.Length` per read |
| `minQuality` | `int` | 20 | Quality cutoff `C` subtracted from each Phred score | `< 1` disables trimming [2] |
| `minLength` | `int` | 50 | Reads shorter than this after trimming are dropped | ≥ 0 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| return | `IReadOnlyList<string>` | Surviving trimmed sequences, in input order |

### 3.3 Preconditions and Validation

`reads` must be non-null (`ArgumentNullException` otherwise). Quality is decoded as Phred = `ASCII − 33` (0-based indexing). Each read's quality string is assumed equal in length to its sequence. Trimming is applied to the 3' end then the 5' end of the surviving window. A read whose trimmed length is `< minLength` (including fully removed reads, length 0) is excluded from the output.

## 4. Algorithm

### 4.1 High-Level Steps

1. If `minQuality < 1`, skip trimming for the read (keep full sequence) [2].
2. 3' end: over the full read `quality[0..n)`, scan from the end accumulating `(Phred − cutoff)`; stop as soon as the running sum becomes positive (the cutadapt/BWA `s < 0` early break, sign-flipped); set `end` to the index of the minimal partial sum reached before the break [1][2].
3. 5' end: independently over the full read `quality[0..n)`, scan from the start accumulating `(Phred − cutoff)` with the same early break; set `start` to one past the index of the minimal partial sum [1][2].
4. If `start ≥ end`, the good-quality segment is empty (cutadapt `start >= stop ⇒ (0,0)`): drop the read.
5. Otherwise, if `end − start ≥ minLength`, emit `sequence[start..end)`; else drop the read.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

- Phred decoding offset = 33 (named constant `PhredAsciiOffset`) [3].
- Min-length filter drops survivors below `minLength` (cutadapt `--minimum-length` semantics).
- This implementation reproduces cutadapt's `quality_trim_index` exactly, including the `s < 0` early break (the accumulator stops once the high-quality suffix/prefix outweighs the low-quality run) and the `start >= stop ⇒ empty` rule [1][2]. The early break is essential: it is what "allows some good-quality bases among the bad ones" — a high-quality base near an end protects the bases on its interior side from being trimmed by that end's pass. A naive global-minimum scan *without* the break over-trims (it would cut at the deepest minimum anywhere in the read).
- BWA additionally floors trimmed length at `BWA_MIN_RDLEN = 35` [2]; this implementation omits that fixed floor (the configurable `minLength` filter plays the corresponding role), so reads can be trimmed below 35 bases unless `minLength` forbids it.

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Trim all reads | O(n · r) | O(1) extra per read | n = reads, r = read length; two linear passes per read |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [SequenceAssembler.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.Alignment/SequenceAssembler.cs)

- `SequenceAssembler.QualityTrimReads(reads, minQuality, minLength)`: public running-sum trimmer with min-length filter.
- `SequenceAssembler.TrimEnd / TrimStart` (private): single-end running-sum cut helpers.

### 5.2 Current Behavior

Trimming runs the 3' pass and the 5' pass independently over the full read (cutadapt `quality_trim_index`), each with the `s < 0` early break; if the resulting windows cross (`start >= stop`) the read is dropped. A cutoff `< 1` returns the read unchanged. Reads with length `< minLength` after trimming (or `0`) are omitted. This is a single short scan per read, not a substring-search problem, so the repository suffix tree is not applicable.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Running-sum cut at the index of the minimal partial sum of `(q − cutoff)`, applied independently to both ends, **including the `s < 0` early break** [1][2].
- The `start >= stop ⇒ empty good-quality segment` rule [1].
- Phred+33 decoding (`ASCII − 33`) [2][3].
- Cutoff `< 1` disables trimming [2].

**Intentionally simplified:**

- BWA's fixed `BWA_MIN_RDLEN = 35` floor is omitted; the configurable `minLength` filter plays the corresponding role; **consequence:** trimmed boundaries match cutadapt exactly and reads can be trimmed below 35 bases unless `minLength` forbids it [1][2].

**Not implemented:**

- Phred+64 (legacy Illumina) decoding; **users should rely on:** Phred+33 (Sanger) input, the current standard [3].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Both-end passes run independently over the full read | Resolved | Matches cutadapt `quality_trim_index` exactly | accepted | Cutadapt runs the 5' and 3' passes independently, then drops the read if `start >= stop` [1] |
| 2 | `minLength` post-trim filter | Assumption | Drops short survivors | accepted | cutadapt min-length semantics |

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

### 6.2 Limitations

Assumes Phred+33 (Sanger) encoding; legacy Phred+64 input would be mis-decoded. Assumes quality string length equals sequence length per read. Does not perform adapter removal or per-window (sliding-window) trimming.

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

1. Cutadapt project. 2024. Algorithm details — Cutadapt documentation (quality trimming). https://cutadapt.readthedocs.io/en/stable/algorithms.html
2. Li, H. BWA source: `bwa_trim_read` (bwaseqio.c) and `BWA_MIN_RDLEN` (bwtaln.h). https://github.com/lh3/bwa/blob/master/bwaseqio.c , https://github.com/lh3/bwa/blob/master/bwtaln.h
3. Cock, P.J.A., Fields, C.J., Goto, N., Heuer, M.L., Rice, P.M. 2010. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. Nucleic Acids Research 38(6):1767–1771. https://doi.org/10.1093/nar/gkp1137
