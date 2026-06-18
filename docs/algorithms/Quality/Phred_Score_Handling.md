# Phred Score Handling

| Field | Value |
|-------|-------|
| Algorithm Group | Quality |
| Test Unit ID | QUALITY-PHRED-001 |
| Related Projects | Seqeron.Genomics.IO |
| Implementation Status | Production |
| Last Reviewed | 2026-06-13 |

## 1. Overview

Phred score handling converts between the ASCII characters of a FASTQ quality line and the integer Phred quality scores they encode, and converts a quality line between the two Phred ASCII encodings (Phred+33 and Phred+64). It is a specification-driven, exact operation: each character maps to exactly one Phred score and back, with no approximation. It is used whenever sequencing quality data must be read, written, or normalized between platform encodings [1].

## 2. Scientific / Formal Basis

### 2.1 Domain Context

FASTQ files store, per base, a Phred quality score Q encoded as a single ASCII character. Two Phred ASCII encodings are in common use: the Sanger / Illumina 1.8+ "Phred+33" encoding and the older Illumina 1.3–1.7 "Phred+64" encoding. They differ only in the ASCII offset added to the score; the Phred score itself has the same meaning in both [1].

### 2.2 Core Model

The Phred quality score is defined as Q = −10 log₁₀(P), where P is the estimated probability that the base call is wrong [1].

Character ↔ score relations [1]:

- Decode: Q = ord(char) − offset
- Encode: char = chr(Q + offset)

Per Cock et al. (2010) [1]:

- **Phred+33 (Sanger / Illumina 1.8+):** offset 33; ASCII 33–126 encodes Phred scores 0–93.
- **Phred+64 (Illumina 1.3+):** offset 64; ASCII 64–126 encodes Phred scores 0–62.

Because the Phred score is identical across both variants, converting a quality line from one variant to the other is a pure re-offset: decode with the source offset, re-encode with the target offset (a constant byte shift of ±31) [1]. (Solexa scores, defined as Q_solexa = −10 log₁₀(P/(1−P)), are a different, non-Phred scale and are out of scope for this unit [1].)

### 2.4 Properties and Invariants

| ID | Invariant | Holds because |
|----|-----------|---------------|
| INV-01 | Phred+33: Q = ord(char) − 33, valid for Q ∈ [0, 93] (ASCII 33–126) | Cock et al. (2010) [1] |
| INV-02 | Phred+64: Q = ord(char) − 64, valid for Q ∈ [0, 62] (ASCII 64–126) | Cock et al. (2010) [1] |
| INV-03 | `ToQualityString(ParseQualityString(s, e), e) == s` for any valid quality string `s` under encoding `e` | Decode and encode are mutual inverses on the valid range [1] |
| INV-04 | `ConvertEncoding` preserves the Phred score across variants | Phred score is variant-invariant; conversion is a re-offset [1] |

### 2.5 Comparison with Related Methods

| Aspect | Phred encoding (this unit) | Solexa encoding |
|--------|----------------------------|-----------------|
| Score formula | Q = −10 log₁₀(P) [1] | Q = −10 log₁₀(P/(1−P)) [1] |
| Min score | 0 | −5 [1] |
| Cross-conversion | lossless re-offset between Phred+33/Phred+64 [1] | lossy numeric conversion to Phred [1] |

## 3. Contract

### 3.1 Inputs and Parameters

| Name | Type | Default | Description | Constraints |
|------|------|---------|-------------|-------------|
| qualityString | `string` | required | FASTQ quality line | non-null; each char must decode within the encoding's valid Phred range |
| scores | `IReadOnlyList<int>` | required | Phred scores to encode | non-null; each score within the encoding's valid range |
| encoding | `QualityEncoding` | `Phred33` | Phred+33 / Phred+64 / Auto | Auto resolved via `DetectEncoding` on parse |
| fromEncoding / toEncoding | `QualityEncoding` | required | Source/target encoding for conversion | Phred+33 or Phred+64 |

### 3.2 Output / Return Value

| Field | Type | Description |
|-------|------|-------------|
| (ParseQualityString) | `int[]` | Phred scores, one per input character; empty for empty input |
| (ToQualityString) | `string` | ASCII quality string; empty for empty input |
| (ConvertEncoding) | `string` | Quality string re-encoded under the target encoding |

### 3.3 Preconditions and Validation

Inputs are 0-based, ASCII. Null inputs raise `ArgumentNullException`. On parse, a character that decodes to a Phred score outside the encoding's valid range ([0,93] for Phred+33, [0,62] for Phred+64 [1]) raises `ArgumentOutOfRangeException`. On encode, a score outside that range raises `ArgumentOutOfRangeException` — this also covers the Phred+33→Phred+64 overflow case (a source score > 62 has no Phred+64 representation [1]). Empty input yields empty output. Encoding is case-sensitive ASCII; no T↔U or case normalization applies to quality lines.

## 4. Algorithm

### 4.1 High-Level Steps

1. Select offset and valid score range from the encoding (33/[0,93] or 64/[0,62]).
2. **Parse:** for each char, compute Q = char − offset; validate Q within range; collect.
3. **Encode:** for each score, validate within range; emit char = score + offset.
4. **Convert:** parse with source encoding, then encode with target encoding.

### 4.2 Decision Rules, Scoring, Reference Tables, or Data Structures

| Encoding | Offset | ASCII range | Phred range | Source |
|----------|--------|-------------|-------------|--------|
| Phred+33 (Sanger/Illumina 1.8+) | 33 | 33–126 | 0–93 | [1] |
| Phred+64 (Illumina 1.3+) | 64 | 64–126 | 0–62 | [1] |

### 4.3 Complexity

| Operation | Time | Space | Notes |
|-----------|------|-------|-------|
| Parse / Encode / Convert | O(n) | O(n) | single linear pass over the n characters/scores |

## 5. Implementation Notes

### 5.1 Location and Entry Points

**Implementation location:** [QualityScoreAnalyzer.cs](../../../src/Seqeron/Algorithms/Seqeron.Genomics.IO/QualityScoreAnalyzer.cs)

- `QualityScoreAnalyzer.ParseQualityString(qualityString, encoding)`: decodes a quality string to Phred scores with range validation.
- `QualityScoreAnalyzer.ToQualityString(scores, encoding)`: encodes Phred scores to a quality string with range validation.
- `QualityScoreAnalyzer.ConvertEncoding(qualityString, fromEncoding, toEncoding)`: re-offsets a quality string between Phred encodings.

### 5.2 Current Behavior

Offsets and valid ranges are named constants citing Cock et al. (2010). `ParseQualityString` resolves `Auto` via the existing `DetectEncoding` heuristic; `ToQualityString` treats `Auto` as the modern `Phred+33` default. This is not a search/matching unit, so the repository suffix tree is not applicable. The pre-existing `QualityStringToPhred`/`PhredToQualityString` helpers remain for other callers; the new canonical methods add explicit range validation that those helpers do not perform.

### 5.3 Conformance to Theory / Spec

**Implemented (verbatim from the cited theory/spec):**

- Phred+33 decode/encode with offset 33, valid scores 0–93 (ASCII 33–126) [1].
- Phred+64 decode/encode with offset 64, valid scores 0–62 (ASCII 64–126) [1].
- Lossless Phred+33 ↔ Phred+64 conversion as a score-preserving re-offset [1].

**Intentionally simplified:**

- (none)

**Not implemented:**

- Solexa (`fastq-solexa`) encoding and the lossy Solexa↔Phred numeric conversion; **users should rely on:** a dedicated Solexa converter (not provided) — only the two Phred variants are in scope [1].

### 5.4 Deviations and Assumptions

| # | Item | Type | Impact | Status | Notes |
|---|------|------|--------|--------|-------|
| 1 | Exception type for malformed/out-of-range chars and overflow is `ArgumentOutOfRangeException` | Assumption | Determines failure mode users catch | accepted | Range bounds are source-backed [1]; the .NET exception type is an API-shape choice |

## 6. Edge Cases and Limitations

### 6.1 Edge Cases

| Case | Expected Behavior | Rationale |
|------|-------------------|-----------|
| Empty string / empty scores | empty result | identity boundary |
| Null input | `ArgumentNullException` | public-API contract |
| Char below offset (negative Q) | `ArgumentOutOfRangeException` | Phred Q ≥ 0 [1] |
| Phred+33 score > 62 → Phred+64 | `ArgumentOutOfRangeException` | not representable in Phred+64 (max 62) [1] |
| Phred+64 → Phred+33 (any valid score) | always succeeds | Phred+33 range (0–93) ⊇ Phred+64 range (0–62) [1] |

### 6.2 Limitations

Handles only the Phred+33 and Phred+64 ASCII encodings; Solexa scores are not supported. Validation is range-based per the spec; it does not enforce platform-specific narrower ranges (e.g. raw Illumina 0–40) since those are expectations, not format limits [1].

## 7. Examples and Related Material

### 7.1 Worked Example

**API usage example:**

```csharp
int[] q = QualityScoreAnalyzer.ParseQualityString("!5?I~", QualityScoreAnalyzer.QualityEncoding.Phred33);
// q == [0, 20, 30, 40, 93]
string s = QualityScoreAnalyzer.ToQualityString(q, QualityScoreAnalyzer.QualityEncoding.Phred33);
// s == "!5?I~"
string p33 = QualityScoreAnalyzer.ConvertEncoding("@h~",
    QualityScoreAnalyzer.QualityEncoding.Phred64, QualityScoreAnalyzer.QualityEncoding.Phred33);
// p33 == "!I_"  (Phred scores 0, 40, 62 preserved)
```

### 7.3 Related Tests, Evidence, or Documents

- Tests: [QualityScoreAnalyzer_ParseQualityString_Tests.cs](../../../tests/Seqeron/Seqeron.Genomics.Tests/QualityScoreAnalyzer_ParseQualityString_Tests.cs) — covers `INV-01`–`INV-04`
- Evidence: [QUALITY-PHRED-001-Evidence.md](../../../docs/Evidence/QUALITY-PHRED-001-Evidence.md)

## 8. References

1. Cock, P.J.A., Fields, C.J., Goto, N., Heuer, M.L., Rice, P.M. 2010. The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. *Nucleic Acids Research* 38(6):1767–1771. https://doi.org/10.1093/nar/gkp1137
