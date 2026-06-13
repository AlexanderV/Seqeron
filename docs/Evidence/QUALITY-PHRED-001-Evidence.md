# Evidence Artifact: QUALITY-PHRED-001

**Test Unit ID:** QUALITY-PHRED-001
**Algorithm:** Phred Score Handling (FASTQ quality string parsing, encoding, and Phred+33 ↔ Phred+64 conversion)
**Date Collected:** 2026-06-13

---

## Online Sources

### Cock et al. (2010) — The Sanger FASTQ file format and Solexa/Illumina variants

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/
**Accessed:** 2026-06-13
**Retrieval:** WebFetch of the PMC article page (PMC2847217). The NCBI `/pmc/articles/` URL 301-redirected to the `pmc.ncbi.nlm.nih.gov` host, which was then fetched directly.
**Authority rank:** 1 (peer-reviewed paper; *Nucleic Acids Research*) — also the de-facto FASTQ format specification.

**Key Extracted Points:**

1. **Sanger / Phred+33 encoding:** Verbatim: "Sanger FASTQ files use ASCII 33–126 to encode PHRED qualities from 0 to 93 (i.e. PHRED scores with an ASCII offset of 33)." → offset = 33, valid Phred Q ∈ [0, 93], valid ASCII ∈ [33, 126].
2. **Illumina 1.3+ / Phred+64 encoding:** Verbatim: "The Illumina 1.3+ FASTQ variant encodes PHRED scores with an ASCII offset of 64, and so can hold PHRED scores from 0 to 62 (ASCII 64–126)." → offset = 64, valid Phred Q ∈ [0, 62], valid ASCII ∈ [64, 126].
3. **Phred quality score definition:** Q = −10 log₁₀(P), where P is the estimated probability of a base-call error.
4. **Char ↔ score relation (Sanger):** ASCII character = chr(Q + 33); i.e. Q = ord(char) − 33. The same relation with offset 64 applies to Illumina 1.3+.
5. **Cross-variant conversion:** The Phred quality score is invariant across the Sanger and Illumina 1.3+ variants; only the ASCII offset differs. Converting between "fastq-sanger" (offset 33) and "fastq-illumina" (offset 64) is a pure re-offset of the same Phred score. (Solexa scores, by contrast, require Equation (3) numeric conversion and are not Phred — out of scope for this unit.)
6. **Solexa note (out of scope):** Solexa quality scores are defined as Q_solexa = −10 log₁₀(P/(1−P)) with offset 64 (ASCII 59–126, scores −5 to 62). Conversion Solexa→Phred is lossy. This unit handles only the two Phred variants.

---

## Documented Corner Cases and Failure Modes

### From Cock et al. (2010)

1. **Out-of-range characters:** A Phred+33 file uses ASCII 33–126; a Phred+64 file uses ASCII 64–126. A character below the offset would decode to a negative Phred score, which is not valid for either Phred variant (Phred Q ≥ 0). Such input is malformed for the stated encoding.
2. **Phred+64 → Phred+33 always representable:** Since Phred+64 holds Q ∈ [0, 62] and Phred+33 holds Q ∈ [0, 93], every valid Phred+64 score is representable in Phred+33.
3. **Phred+33 → Phred+64 may overflow:** A Phred+33 score in (62, 93] exceeds the Phred+64 maximum (62) and cannot be encoded in Phred+64.

---

## Test Datasets

### Dataset: Sanger Phred+33 worked characters

**Source:** Cock et al. (2010), §"Sanger FASTQ format" (offset 33; ASCII 33–126 → Q 0–93).

| Character | ASCII | Phred Q (Phred+33) |
|-----------|-------|--------------------|
| `!`       | 33    | 0  |
| `5`       | 53    | 20 |
| `?`       | 63    | 30 |
| `I`       | 73    | 40 |
| `~`       | 126   | 93 |

### Dataset: Illumina 1.3+ Phred+64 worked characters

**Source:** Cock et al. (2010), §"Illumina 1.3+ FASTQ variant" (offset 64; ASCII 64–126 → Q 0–62).

| Character | ASCII | Phred Q (Phred+64) |
|-----------|-------|--------------------|
| `@`       | 64    | 0  |
| `h`       | 104   | 40 |
| `~`       | 126   | 62 |

### Dataset: Phred+64 → Phred+33 conversion (Phred score preserved)

**Source:** Derived from Cock et al. (2010) points 4–5 (Phred score invariant; re-offset by −64 then +33, i.e. shift each byte by −31).

| Input (Phred+64) | Q | Output (Phred+33) |
|------------------|---|-------------------|
| `@`  (64)        | 0  | `!`  (33) |
| `h`  (104)       | 40 | `I`  (73) |
| `~`  (126)       | 62 | `_`  (95) |

---

## Assumptions

1. **ASSUMPTION: Validation strictness on decode.** Cock et al. define the valid ASCII/score ranges per variant but do not prescribe an exact exception type for malformed bytes. This unit raises `ArgumentOutOfRangeException` when a character decodes to a Phred score outside the variant's documented valid range (Phred+33: [0,93]; Phred+64: [0,62]) — chosen to surface malformed input rather than silently producing negative/invalid scores. Range bounds themselves are source-backed; only the choice of exception type is an assumption (API-shape, non-correctness-affecting).
2. **ASSUMPTION: Phred+33 → Phred+64 overflow.** When converting a Phred+33 score > 62 to Phred+64 (not representable, per corner case 3), this unit raises `ArgumentOutOfRangeException`. The non-representability is source-backed; the exception type is the assumption.

---

## Recommendations for Test Coverage

1. **MUST Test:** `ParseQualityString` decodes Phred+33 boundary chars `!`→0, `~`→93 and interior `5`→20, `?`→30, `I`→40 — Evidence: Cock et al. (2010) point 1.
2. **MUST Test:** `ParseQualityString` decodes Phred+64 boundary chars `@`→0, `~`→62 and interior `h`→40 — Evidence: Cock et al. (2010) point 2.
3. **MUST Test:** `ToQualityString` encodes scores back to the exact characters (round-trip with parse) for both variants — Evidence: Cock et al. (2010) point 4.
4. **MUST Test:** `ConvertEncoding` Phred+64→Phred+33 maps `@h~` → `!I_` preserving Phred scores — Evidence: Cock et al. (2010) points 4–5.
5. **MUST Test:** `ConvertEncoding` Phred+33→Phred+64 maps `!I` → `@h` (Q 0,40) — Evidence: Cock et al. (2010) points 4–5.
6. **MUST Test:** Out-of-range decode (char below offset → negative Q) throws — Evidence: corner case 1.
7. **MUST Test:** Phred+33→Phred+64 overflow (Q>62, e.g. `~`=93) throws — Evidence: corner case 3.
8. **SHOULD Test:** empty string parses to empty / encodes to empty (identity) — Rationale: trivial boundary.
9. **SHOULD Test:** null input handling — Rationale: documented failure mode for public API.
10. **COULD Test:** round-trip property parse∘encode = identity over a random valid Phred+33 string with fixed seed — Rationale: invariant check.

---

## References

1. Cock, P.J.A., Fields, C.J., Goto, N., Heuer, M.L., Rice, P.M. (2010). The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants. *Nucleic Acids Research*, 38(6), 1767–1771. https://doi.org/10.1093/nar/gkp1137 (retrieved via https://pmc.ncbi.nlm.nih.gov/articles/PMC2847217/, accessed 2026-06-13)

---

## Change History

- **2026-06-13**: Initial documentation.
