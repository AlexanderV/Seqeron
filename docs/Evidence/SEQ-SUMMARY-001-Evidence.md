# Evidence Artifact: SEQ-SUMMARY-001

**Test Unit ID:** SEQ-SUMMARY-001
**Algorithm:** Sequence Summary (aggregation of nucleotide composition + Shannon entropy + linguistic complexity + melting temperature into one record)
**Date Collected:** 2026-06-14

---

## Online Sources

<!-- SEQ-SUMMARY aggregates already-validated per-metric methods. The summary's correctness
     is defined field-by-field: each summary field must equal the value produced by its
     canonical per-metric method on the same input. The sources below establish the canonical
     definition of each aggregated metric so the per-field expectations are evidence-derived. -->

### Biopython — `Bio.SeqUtils` source (`gc_fraction`)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
**Accessed:** 2026-06-14 (WebFetch of the raw source file in this session)
**Authority rank:** 3 (reference implementation in an established bioinformatics library)

**Key Extracted Points:**

1. **GC content count expression (verbatim):** `gc = sum(seq.count(x) for x in "CGScgs")`; the function returns `gc / length`, a float. GC counting includes lowercase, so it is case-insensitive.
2. **Length under the default `ambiguous="remove"` (verbatim):** `length = gc + sum(seq.count(x) for x in "ATWUatwu")` — GC fraction is GC over counted unambiguous bases.
3. **Empty sequence (verbatim):** `if length == 0: return 0` before the division, so an all-non-counted / empty sequence returns 0.

### Biopython — `Bio.SeqUtils.MeltingTemp` source (`Tm_Wallace`, `Tm_GC`)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py
**Accessed:** 2026-06-14 (WebFetch of the raw source file in this session)
**Authority rank:** 3 (reference implementation)

**Key Extracted Points:**

1. **Wallace rule (verbatim docstring):** "Tm = 4 degC * (G + C) + 2 degC * (A+T)". Docstring worked example `mt.Tm_Wallace('ACGTTGCAATGCCGTA')` returns `48.0`. Described as a rule of thumb for primers of ~14–20 nt.
2. **GC formula general form (verbatim):** "Tm = A + B(%GC) - C/N + salt correction - D(%mismatch)".
3. **Marmur & Doty 1962 valueset (verbatim):** "Tm = 69.3 + 0.41(%GC) - 650/N" (A=69.3, B=0.41, C=650, N = length). This is the basic GC/Marmur-Doty Tm family used for longer sequences.

### Wikipedia — "Entropy (information theory)" (citing the primary Shannon 1948)

**URL:** https://en.wikipedia.org/wiki/Entropy_(information_theory)
**Accessed:** 2026-06-14 (WebFetch in this session)
**Authority rank:** 4 (Wikipedia citing the primary Shannon 1948)

**Key Extracted Points:**

1. **Shannon entropy formula (verbatim):** `H(X) := − ∑ x∈X p(x) log p(x)`.
2. **Unit (verbatim):** "Base 2 gives the unit of bits (or 'shannons')." Hence `log₂` yields entropy in bits.
3. **Maximum entropy (verbatim):** "The maximal entropy of an event with n different outcomes is log_b(n): it is attained by the uniform probability distribution." So a sequence over k equally-frequent symbols has H = log₂(k) bits.
4. **Primary citation:** Claude Shannon (1948), "A Mathematical Theory of Communication".

### Wikipedia — "Linguistic sequence complexity" (citing the primary Trifonov 1990)

**URL:** https://en.wikipedia.org/wiki/Linguistic_sequence_complexity
**Accessed:** 2026-06-14 (WebFetch in this session)
**Authority rank:** 4 (Wikipedia citing the primary Trifonov 1990)

**Key Extracted Points:**

1. **Vocabulary usage U (verbatim):** for oligomers of a given size, U is "the ratio of the actual vocabulary size of a given sequence to the maximal possible vocabulary size for a sequence of that length."
2. **Combination across word sizes (verbatim):** complexity combines as a **product** `C = U₁U₂…Uᵢ…Uw`.
3. **Range (verbatim):** "0 < C < 1 for DNA sequence fragments."
4. **Primary citation:** Edward Trifonov (1990).

---

## Documented Corner Cases and Failure Modes

### From Biopython `gc_fraction`

1. **Empty sequence:** returns 0 (zero-length denominator guarded). The summary's GcContent inherits this: an empty input yields GcContent = 0.

### From Biopython `MeltingTemp`

1. **Length-dependent formula selection:** the Wallace rule is documented as suitable for short oligos; the GC/Marmur-Doty formula applies to longer sequences. The summary selects Wallace when `length < 14` and the GC formula otherwise (a length-dispatch consistent with the sibling SEQ-TM-001 unit and the `ThermoConstants.WallaceMaxLength = 14` boundary).

### Aggregation-specific (this unit)

1. **Field consistency:** because the summary is a pure aggregation, every field MUST equal the value its canonical per-metric method returns on the identical input. Any divergence (different rounding, different alphabet handling, different formula selection) is a defect of the summary, independent of whether each underlying metric is itself correct (that is the scope of the sibling units SEQ-COMPOSITION-001, SEQ-ENTROPY-PROFILE-001, SEQ-TM-001).

---

## Test Datasets

### Dataset: Balanced tetramer repeat "ATGCATGC"

**Source:** Derived directly from the cited formulas (composition counting; Shannon H = log₂(k); Wallace Tm = 2(A+T)+4(G+C)).

| Parameter | Value |
|-----------|-------|
| Sequence | ATGCATGC |
| Length | 8 |
| Composition A,T,G,C,U,N | 2,2,2,2,0,0 |
| GcContent | 0.5 (= 4/8) |
| Entropy (bits) | 2.0 (= log₂ 4, four symbols each freq 1/4) |
| MeltingTemperature (°C) | 24.0 (length<14 → Wallace = 2·4 + 4·4) |
| Complexity | equals `CalculateLinguisticComplexity("ATGCATGC")` (independently 0.83968253968…, mean of per-k vocabulary-usage ratios) |

### Dataset: 16-mer "ATGCATGCATGCATGC" (Marmur-Doty branch)

**Source:** Marmur-Doty GC Tm family; repo variant `64.9 + 41·(GC − 16.4)/N` (validated in SEQ-TM-001).

| Parameter | Value |
|-----------|-------|
| Sequence | ATGCATGCATGCATGC |
| Length | 16 |
| GC count | 8 |
| MeltingTemperature (°C) | 43.375 (length≥14 → GC formula = 64.9 + 41·(8−16.4)/16) |

---

## Assumptions

1. **ASSUMPTION: Tm formula-selection threshold (length < 14).** The summary passes `useWallaceRule: sequence.Length < 14` to the melting-temperature method. The 14 nt boundary is the sibling SEQ-TM-001 convention (`ThermoConstants.WallaceMaxLength`); Biopython documents Wallace as a rule of thumb for ~14–20 nt without fixing an exact switch point. This is a non-correctness-affecting choice *for the summary*, because the summary's contract is "MeltingTemperature equals `CalculateMeltingTemperature` with this flag" — the threshold belongs to the already-validated SEQ-TM-001 unit, and the summary is tested for equality with that canonical method on the same input.

---

## Recommendations for Test Coverage

1. **MUST Test:** Each summary field equals its canonical per-metric method's value on the same input (Length, GcContent, Entropy, Complexity, MeltingTemperature, Composition dictionary). — Evidence: per-metric source definitions above; aggregation-consistency requirement.
2. **MUST Test:** Exact values on the worked example "ATGCATGC" (GcContent 0.5, Entropy 2.0, Tm 24.0). — Evidence: derived from the cited formulas.
3. **MUST Test:** Tm uses the GC/Marmur-Doty branch for length ≥ 14 (e.g. 43.375 for the 16-mer). — Evidence: Biopython length-dependent formula selection + SEQ-TM-001 constants.
4. **SHOULD Test:** Composition dictionary keys/counts (A,T,G,C,U,N) match `CalculateNucleotideComposition`; U and N counted for RNA / ambiguous inputs. — Rationale: composition is the only nested structure in the record.
5. **SHOULD Test:** Empty and null input return the zero/degenerate summary (Length 0, GcContent 0, Entropy 0, Complexity 0, Tm 0). — Rationale: documented empty-sequence behavior.
6. **COULD Test:** Case-insensitivity (lowercase input gives identical summary). — Rationale: per-metric methods uppercase internally.

---

## References

1. Shannon, C. E. (1948). A Mathematical Theory of Communication. Bell System Technical Journal, 27(3):379–423. https://doi.org/10.1002/j.1538-7305.1948.tb01338.x (formula/units confirmed via https://en.wikipedia.org/wiki/Entropy_(information_theory), fetched 2026-06-14).
2. Trifonov, E. N. (1990). Making sense of the human genome. In: Structure & Methods (linguistic complexity / vocabulary usage). Definition confirmed via https://en.wikipedia.org/wiki/Linguistic_sequence_complexity (fetched 2026-06-14).
3. Marmur, J., Doty, P. (1962). Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. J Mol Biol 5:109–118. GC Tm form confirmed via Biopython `Bio.SeqUtils.MeltingTemp` source, https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py (fetched 2026-06-14).
4. Wallace, R. B. et al. (rule of thumb Tm = 2(A+T)+4(G+C)). Confirmed via Biopython `Bio.SeqUtils.MeltingTemp.Tm_Wallace` source/docstring, https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/MeltingTemp.py (fetched 2026-06-14).
5. Cock, P. J. A. et al. (2009). Biopython. Bioinformatics 25(11):1422–1423. https://doi.org/10.1093/bioinformatics/btp163 (gc_fraction source, https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py, fetched 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
