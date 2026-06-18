# Evidence Artifact: ASSEMBLY-CORRECT-001

**Test Unit ID:** ASSEMBLY-CORRECT-001
**Algorithm:** K-mer spectrum (two-sided) read error correction
**Date Collected:** 2026-06-13

---

## Online Sources

### Musket: a multistage k-mer spectrum-based error corrector for Illumina sequence data

**URL:** https://academic.oup.com/bioinformatics/article/29/3/308/257257/Musket-a-multistage-k-mer-spectrum-based-error
**Accessed:** 2026-06-13 (retrieved via WebFetch of the OUP article page)
**Authority rank:** 1 (peer-reviewed paper, *Bioinformatics*) / 3 (reference tool)

**Key Extracted Points:**

1. **Trusted vs untrusted k-mer:** Quoted from the article — "the *k*-mers whose multiplicity exceeds a coverage *cut-off* are deemed to be trusted and otherwise, untrusted." The cut-off is chosen automatically from "the smallest density around the valley" of the k-mer coverage histogram.
2. **Trusted base rule:** "If a base is covered by any trusted *k*-mer, the base is deemed to be trusted and untrusted, otherwise."
3. **Single-substitution assumption:** "our two-sided correction conservatively assumes that there is at most one substitution error in any *k*-mer of a read."
4. **Two-sided correction goal:** "our two-sided correction aims to find a unique alternative base that makes all *k*-mers that cover position *i* trusted"; it evaluates "both the leftmost and the rightmost *k*-mers that cover position *i*."
5. **Ambiguity rule:** "For a single base, if more than one alternative is found to make both the leftmost and the rightmost *k*-mers trusted, the base will keep unchanged as a result of ambiguity."

### Quake: quality-aware detection and correction of sequencing errors

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3156955/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full text)
**Authority rank:** 1 (peer-reviewed paper, *Genome Biology*)

**Key Extracted Points:**

1. **Trusted/untrusted definition:** "high coverage *k*-mers as *trusted*, because they are highly likely to occur in the genome, and low coverage *k*-mers as *untrusted*."
2. **Trusted base / region:** corrections keep "every base covered by the right-most trusted *k*-mer to be correct and removing them from the error region" — bases within trusted k-mers are not modified.
3. **Correction search:** corrections are searched until "a set of corrections *C* makes all *k*-mers in the region trusted"; errors are localized first to the intersection, then the union, of the read's untrusted k-mers.
4. **Edit model:** corrections are "single-base edits (edit distance)" / "nucleotide edits" — substitutions over the four nucleotides.

### Mining statistically-solid k-mers for accurate NGS error correction

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6311904/
**Accessed:** 2026-06-13 (retrieved via WebFetch of the PMC full text)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Solid k-mer:** "A k-mer frequently occurring in NGS reads" that "often does not contain any sequencing error."
2. **Weak k-mer:** a k-mer not meeting the frequency threshold f0, which "often contains sequencing errors."
3. **Two-sided correction:** a correction is accepted only when forward and backward searches meet and "the k-mers contained in both of them are solid," ensuring unambiguous solutions.

---

## Documented Corner Cases and Failure Modes

### From Musket (Liu et al. 2013)

1. **Ambiguous position:** more than one alternative base makes the covering k-mers trusted → the base is left unchanged (no correction made).
2. **At-most-one-error assumption:** the two-sided rule assumes ≤1 substitution per k-mer; multiple close errors within one k-mer may not be correctable.

### From Quake (Kelley et al. 2010)

1. **No correcting set:** when no single set of corrections makes the region's k-mers trusted, the read is left uncorrected (and in Quake, may be trimmed/discarded — out of scope here).
2. **Solid k-mer containing an error / weak k-mer without an error:** the frequency cut-off cannot perfectly separate erroneous from correct k-mers (general limitation of k-spectrum methods).

---

## Test Datasets

### Dataset: Derived single-substitution worked example (k=3, cut-off=2)

**Source:** Constructed from the Musket/Quake definitions above; all intermediate counts derived by hand and confirmable from the formal rules.

| Parameter | Value |
|-----------|-------|
| Reads | `ACGTACGT` ×3 (true), `ACGTTCGT` ×1 (substitution A→T at index 4) |
| k | 3 |
| Trusted cut-off (minKmerFrequency) | 2 |
| k-mer spectrum | ACG=7, CGT=8, GTA=3, TAC=3, GTT=1, TTC=1, TCG=1 |
| Trusted k-mers (≥2) | ACG, CGT, GTA, TAC |
| Untrusted k-mers (<2) | GTT, TTC, TCG |
| Erroneous read positions covered only by untrusted k-mers | index 4 only |
| Unique correcting base at index 4 | A (C→GTC absent; G→GTG absent; only A→GTA,TAC,ACG all trusted) |
| Expected output | `ACGTACGT` ×4 (error read fully corrected; true reads unchanged) |

### Dataset: Ambiguity example (k=1, cut-off=2)

**Source:** Constructed to exercise the Musket ambiguity rule (>1 valid alternative → unchanged).

| Parameter | Value |
|-----------|-------|
| Reads | `A`, `A`, `C`, `C`, `G` |
| k | 1 |
| Trusted cut-off | 2 |
| 1-mer spectrum | A=2, C=2, G=1 |
| Trusted 1-mers (≥2) | A, C |
| Correction of read `T` | candidates A and C both yield a trusted 1-mer → ambiguous → unchanged |
| Expected output for `T` | `T` (unchanged) |

---

## Assumptions

1. **ASSUMPTION: Default parameter values** — Musket/Quake choose k and the coverage cut-off automatically from the data (coverage-histogram valley); the library exposes both as parameters and uses fixed defaults (`kmerSize=15`, `minKmerFrequency=2`). Defaults are non-behavioral for the tested contract because every behavioral test passes k and cut-off explicitly; the algorithm's correctness for given (k, cut-off) is fully source-defined.

---

## Recommendations for Test Coverage

1. **MUST Test:** single-substitution read is fully corrected to the trusted sequence (worked example k=3, cut-off=2) — Evidence: Musket two-sided rule (find unique base making covering k-mers trusted).
2. **MUST Test:** position covered by a trusted k-mer is never modified — Evidence: Musket/Quake "base covered by any trusted k-mer is trusted."
3. **MUST Test:** ambiguous position (>1 valid alternative) is left unchanged — Evidence: Musket ambiguity rule.
4. **MUST Test:** no valid correcting base → base unchanged — Evidence: Quake "no set of corrections makes the region trusted" → uncorrected.
5. **MUST Test:** output count and per-read length are preserved (substitution-only model) — Evidence: Quake/Musket nucleotide-edit (substitution) model.
6. **SHOULD Test:** correct reads pass through unchanged — Rationale: trusted reads must not be altered.
7. **SHOULD Test:** case-insensitive handling (lowercase input upper-cased) — Rationale: matches existing analyzer convention.
8. **COULD Test:** reads shorter than k contribute no k-mers and are returned unchanged — Rationale: boundary of the k-mer window.
9. **MUST Test:** null reads → ArgumentNullException; kmerSize < 1 → ArgumentOutOfRangeException — Rationale: documented failure modes.

---

## References

1. Liu Y, Schmidt B, Maskell DL (2013). Musket: a multistage k-mer spectrum-based error corrector for Illumina sequence data. *Bioinformatics* 29(3):308-315. https://doi.org/10.1093/bioinformatics/bts690
2. Kelley DR, Schatz MC, Salzberg SL (2010). Quake: quality-aware detection and correction of sequencing errors. *Genome Biology* 11:R116. https://doi.org/10.1186/gb-2010-11-11-r116
3. Song L, Florea L (2018, retrieved 2026-06-13). Mining statistically-solid k-mers for accurate NGS error correction. PMC6311904. https://pmc.ncbi.nlm.nih.gov/articles/PMC6311904/

---

## Change History

- **2026-06-13**: Initial documentation.
