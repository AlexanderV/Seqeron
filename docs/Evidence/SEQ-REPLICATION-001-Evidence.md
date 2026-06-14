# Evidence Artifact: SEQ-REPLICATION-001

**Test Unit ID:** SEQ-REPLICATION-001
**Algorithm:** Replication Origin Prediction (cumulative GC-skew minimum)
**Date Collected:** 2026-06-14

---

## Online Sources

### Rosalind — Minimum Skew Problem (BA1F)

**URL:** https://rosalind.info/problems/ba1f/
**Accessed:** 2026-06-14
**Authority rank:** 3 (reference problem with an exact published worked example / canonical algorithm definition)

**Retrieved via:** WebFetch of the URL above.

**Key Extracted Points:**

1. **Skew definition:** "The skew of a DNA string Genome, denoted Skew(Genome), [is] the difference between the total number of occurrences of 'G' and 'C' in Genome." Computed by iterating each position: G contributes +1, C contributes −1, and A/T contribute 0, starting from 0 at position 0.
2. **Minimum Skew Problem:** "Find a position in a genome minimizing the skew. Given: A DNA string Genome. Return: All integer(s) i minimizing Skew(Prefix_i(Text)) over all values of i (from 0 to |Genome|)." So positions are 0-based prefix indices in [0, |Genome|].
3. **Sample input (verbatim):** `CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG` (length 100).
4. **Sample output (verbatim):** `53 97` — the positions minimizing the skew.

**Independent re-derivation (in this session):** running the per-nucleotide cumulative skew (G:+1, C:−1, A/T:0; Skew_0 = 0) over the sample input yields a global minimum value of −4 at prefix indices 53 and 97, reproducing the published output exactly.

---

### Grigoriev A (1998) — Analyzing genomes with cumulative skew diagrams

**URL:** https://academic.oup.com/nar/article/26/10/2286/1030593 (DOI: https://doi.org/10.1093/nar/26.10.2286)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary literature, Nucleic Acids Research)

**Retrieved via:** WebSearch ("Grigoriev 1998 Analyzing genomes with cumulative skew diagrams Nucleic Acids Research") then WebFetch of the Oxford Academic article page.

**Key Extracted Points:**

1. **Abstract (verbatim):** "A novel method of cumulative diagrams shows that the nucleotide composition of a microbial chromosome changes at two points separated by about a half of its length. These points coincide with sites of replication origin and terminus for all bacteria where such sites are known."
2. **Construction:** a running (cumulative) sum of the nucleotide skew along the sequence, integrated from a start point to each position.
3. **Origin/terminus location:** the cumulative GC-skew diagram reaches its global minimum at the replication origin and its global maximum near the terminus; the two extrema are separated by roughly half the chromosome length.
4. **Strand bias:** the leading strand contains more guanine than cytosine.

---

### Wikipedia — GC skew (used for its cited primaries)

**URL:** https://en.wikipedia.org/wiki/GC_skew
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia; relied on for the formula statement and its cited primaries Lobry 1996, Grigoriev 1998)

**Retrieved via:** WebSearch ("Lobry 1996 ... GC skew origin replication") then WebFetch of the article.

**Key Extracted Points:**

1. **Formula (verbatim):** "GC skew = (G − C)/(G + C)".
2. **Replication features (verbatim):** "the maximum value of the cumulative skew corresponds to the terminal, and the minimum value corresponds to the origin of replication."
3. **Strand bias (verbatim):** "the leading strand contains more guanine (G) and thymine (T), whereas the lagging strand contains more adenine (A) and cytosine (C)."
4. **Cited primaries:** Lobry, J. R. (1996) Mol Biol Evol 13:660–665; Grigoriev, A. (1998) Nucleic Acids Res 26:2286–2290.

---

### Lobry JR (1996) — Asymmetric substitution patterns in the two DNA strands of bacteria

**URL:** https://pubmed.ncbi.nlm.nih.gov/8676740/ (Mol Biol Evol 13(5):660–665)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary literature)

**Retrieved via:** WebSearch ("Lobry 1996 asymmetric substitution patterns ... GC skew origin replication"); PubMed record and the review summarizing it confirmed the points below.

**Key Extracted Points:**

1. Lobry first reported (1996) compositional asymmetry between the two strands in three bacteria (E. coli, B. subtilis, H. influenzae): departure from intrastrand A=T and C=G equifrequency.
2. GC/AT skews switch sign at the origin and terminus of replication; this sign switch is used to confirm/predict the origin of replication.

---

## Documented Corner Cases and Failure Modes

### From Rosalind BA1F

1. **Ties:** the problem asks for ALL positions minimizing the skew (e.g. `53 97`); an implementation returning a single position must define a deterministic tie-break (this unit returns the first/smallest minimizing index).
2. **Prefix indexing:** positions range over i ∈ [0, |Genome|], i.e. there are |Genome|+1 prefix values; Skew_0 = 0 before any base.

### From Grigoriev 1998 / Lobry 1996

3. **No asymmetry → no signal:** if a sequence has no net G/C strand bias the cumulative diagram is flat (amplitude 0) and the origin/terminus are not meaningfully resolved.

---

## Test Datasets

### Dataset: Rosalind BA1F sample

**Source:** Rosalind, Minimum Skew Problem (BA1F), https://rosalind.info/problems/ba1f/

| Parameter | Value |
|-----------|-------|
| Genome | `CCTATCGGTGGATTAGCATGTCCCTGTACGTTTCGCCGCGAACTAGTTCACACGGCTTGATGGCAAATGGTTTTTCCGGCGACCGTAATCGTCCACCGAG` |
| Length | 100 |
| Minimum skew value | −4 |
| Positions of minimum skew | 53, 97 |
| First minimizing position (this unit's tie-break) | 53 |

### Dataset: Tiny worked examples (derived from the BA1F definition)

**Source:** Definition derivation from Rosalind BA1F / Grigoriev 1998.

| Sequence | Skew diagram (Skew_0..Skew_n) | Min value @ pos | Max value @ pos |
|----------|-------------------------------|-----------------|-----------------|
| `CCGGGG` | 0,−1,−2,−1,0,+1,+2 | −2 @ 2 | +2 @ 6 |
| `GGGCCC` | 0,+1,+2,+3,+2,+1,0 | 0 @ 0 | +3 @ 3 |
| `AATT`   | 0,0,0,0,0 | 0 @ 0 | 0 @ 0 |

---

## Assumptions

1. **ASSUMPTION: IsSignificant semantics.** No authoritative source defines a numeric "significance" cutoff for an origin call. The previous implementation used an invented threshold (`amplitude > count × 0.01`); that constant is untraceable and is removed. `IsSignificant` is redefined as the threshold-free, evidence-neutral predicate `max > min` (the diagram has non-zero amplitude, i.e. a detectable strand-composition asymmetry exists per Lobry 1996 / Grigoriev 1998). This is the weakest non-invented definition; callers needing a quantitative confidence measure should inspect the skew amplitude directly.

---

## Recommendations for Test Coverage

1. **MUST Test:** BA1F sample genome returns PredictedOrigin = 53 (first minimizing prefix index). — Evidence: Rosalind BA1F sample output `53 97`.
2. **MUST Test:** Per-nucleotide skew uses G:+1, C:−1, A/T:0 with Skew_0 = 0 (small derived examples `CCGGGG`, `GGGCCC`). — Evidence: Rosalind BA1F definition.
3. **MUST Test:** PredictedTerminus = position of the global maximum. — Evidence: Grigoriev 1998 / Wikipedia (max = terminus).
4. **MUST Test:** First-occurrence tie-break when several positions share the extreme value. — Evidence: BA1F returns multiple minimizers.
5. **SHOULD Test:** Flat diagram (no G/C, or balanced) → origin = terminus = 0, IsSignificant = false. — Rationale: documented "no asymmetry" corner case.
6. **SHOULD Test:** Case-insensitive on the string overload; A/T bases do not move the diagram. — Rationale: counting convention.
7. **COULD Test:** Empty / null input returns zero prediction (string) / throws (DnaSequence). — Rationale: documented input handling.

---

## References

1. Rosalind. Minimum Skew Problem (BA1F). https://rosalind.info/problems/ba1f/
2. Grigoriev, A. (1998). Analyzing genomes with cumulative skew diagrams. Nucleic Acids Research 26(10):2286–2290. https://doi.org/10.1093/nar/26.10.2286
3. Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. Molecular Biology and Evolution 13(5):660–665. https://pubmed.ncbi.nlm.nih.gov/8676740/
4. Wikipedia. GC skew. https://en.wikipedia.org/wiki/GC_skew (accessed 2026-06-14; used for cited primaries 2 and 3).

---

## Change History

- **2026-06-14**: Initial documentation.
