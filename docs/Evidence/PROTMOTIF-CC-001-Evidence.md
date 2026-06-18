# Evidence Artifact: PROTMOTIF-CC-001

**Test Unit ID:** PROTMOTIF-CC-001
**Algorithm:** Coiled-Coil Prediction (heptad-repeat a/d hydrophobic-core detection)
**Date Collected:** 2026-06-14

---

## Online Sources

### Mason JM & Arndt KM (2004) — "Coiled coil domains: stability, specificity, and biological implications" (ChemBioChem 5(2):170–176)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC2567808/ (review building on this work; primary record: https://pubmed.ncbi.nlm.nih.gov/14760737/ , DOI https://doi.org/10.1002/cbic.200300781)
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed review) / 4 (Wikipedia citing this primary)

**Retrieval:** WebSearch "Mason Arndt 2004 coiled coil review ChemBioChem heptad repeat a d positions buried hydrophobic" → opened PMC2567808 with WebFetch.

**Key Extracted Points:**

1. **Heptad notation:** Verbatim from the fetched page — "The heptad repeat, denoted [abcdefg]n, typically has hydrophobic residues at a and d, and polar/charged residues at e and g."
2. **Hydrophobic core:** "Positions a and d contain predominantly hydrophobic amino acids that drive helix association through their burial in the coiled-coil core."
3. **Specificity:** "Interactions among a, d, e and g residues account for most structural specificity in coiled coils." (a/d = hydrophobic core; e/g = electrostatic edge.)

---

### Wikipedia — "Coiled coil" (citing Mason & Arndt 2004)

**URL:** https://en.wikipedia.org/wiki/Coiled_coil
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing primary literature — primary used, not Wikipedia itself)

**Retrieval:** WebFetch of the article.

**Key Extracted Points:**

1. **Heptad pattern:** Verbatim — "Coiled coils usually contain a repeated pattern, hxxhcxc, of hydrophobic (h) and charged (c) amino-acid residues, referred to as a heptad repeat." (h at positions a and d.)
2. **a/d residue identity:** Verbatim — "The positions in the heptad repeat are usually labeled abcdefg, where a and d are the hydrophobic positions, often being occupied by isoleucine, leucine, or valine." This fixes the hydrophobic-core residue set used by this implementation: **{I, L, V}**.
3. **Packing:** "The packing in a coiled-coil interface is exceptionally tight, with almost complete van der Waals contact between the side-chains of the a and d residues." (Burial of the a/d hydrophobic stripe is the driving force.)
4. **Cited primary:** Reference 10 = Mason JM, Arndt KM (Feb 2004), ChemBioChem 5(2):170–176.

---

### Wikipedia — "Heptad repeat"

**URL:** https://en.wikipedia.org/wiki/Heptad_repeat
**Accessed:** 2026-06-14
**Authority rank:** 4 (Wikipedia citing primaries: Chambers et al. 1990, J Gen Virol)

**Retrieval:** WebFetch of the article.

**Key Extracted Points:**

1. **Pattern:** A repeating pattern of seven amino acids at positions a b c d e f g following the motif H P P H C P C: positions a and d are hydrophobic (H), c/g charged (C), b/e/f polar (P).
2. **Leucine zipper:** Leucine zippers have "predominantly leucine in the d position of the heptad repeat" — confirms L is a canonical a/d residue.

---

### Lupas A, Van Dyke M, Stock J (1991) — "Predicting coiled coils from protein sequences" (Science 252:1162–1164)

**URL:** https://www.science.org/doi/10.1126/science.252.5009.1162 ; record https://pubmed.ncbi.nlm.nih.gov/2031185/
**Accessed:** 2026-06-14
**Authority rank:** 1 (peer-reviewed primary)

**Retrieval:** WebSearch "Lupas 1991 predicting coiled coils from protein sequences heptad ... canonical positions residues".

**Key Extracted Points:**

1. **Heptad register / window:** The COILS method scores a residue by sliding a window over the sequence; "Because a window can be assigned seven different heptad repeat frames and a residue can occupy 28 positions in a gliding window, there are 196 preliminary scores for each residue." This establishes (a) the **seven heptad registers** (frames) that must be tried, and (b) the canonical **window of 28 residues** (4 heptads). A residue's score is taken as the maximum over windows/registers.
2. **Heptad definition:** "Coiled coils arise from a characteristic heptad repeat sequence normally denoted as (abcdefg)n with generally hydrophobic residues in positions a and d."
3. **PSSM not used:** The full COILS position-specific residue-frequency table (the 21×20 probabilities) was **not** retrievable in this session. Per the unit's directive, the COILS PSSM is therefore deliberately NOT implemented; this unit implements the fully-specified heptad-register a/d hydrophobic-core detection instead (registers + 28-window are reused from Lupas).

---

## Documented Corner Cases and Failure Modes

### From Lupas et al. (1991) / Mason & Arndt (2004)

1. **Sub-window sequences:** A sequence shorter than the scoring window cannot be scored (no full window exists) → no prediction.
2. **Heptad register ambiguity:** The correct register is unknown a priori; all seven frames must be evaluated and the best taken (Lupas: "seven different heptad repeat frames").
3. **Minimum length:** Naturally occurring coiled coils are built from multiple heptads; isolated single hydrophobic positions are not coiled coils (Mason & Arndt: (abcdefg)1-(abcdefg)2-(abcdefg)3 …).

---

## Test Datasets

### Dataset: GCN4 leucine-zipper canonical heptad (illustrative, fully derivable)

**Source:** Mason & Arndt (2004); Wikipedia "Coiled coil" — a/d = {I,L,V}; leucine at d.

A synthetic peptide built from the repeat `LEIQAQK` (a=L, b=E, c=I, d=Q… ) is NOT used for exact values to avoid ambiguity. Instead exact expected values are derived from minimal, unambiguous constructions (below), since the a/d-occupancy fraction is a closed-form count, not a tabulated empirical value.

### Dataset: Closed-form a/d occupancy (derivation, not tabulated constant)

**Source:** Definition — score = (count of a/d positions whose residue ∈ {I,L,V}) / (count of a/d positions in window), maximized over 7 registers.

| Construction | a/d residues in best register | Score |
|--------------|-------------------------------|-------|
| `(LAALAAA)` repeated (L at a and d every heptad) | all a/d = L | 1.0 |
| All-glycine sequence | no a/d ∈ {I,L,V} | 0.0 |
| `(LAAAAAA)` repeated (L at a only) | half of a/d are L | 0.5 |

---

## Assumptions

1. **ASSUMPTION: Hydrophobic-core residue set = {I, L, V}** — Justification: the verbatim authoritative statement (Wikipedia citing Mason & Arndt 2004) is "a and d ... often being occupied by isoleucine, leucine, or valine." Other residues (A, M, F) appear in some coiled coils but are not part of the single, unambiguously stated set; restricting to {I,L,V} keeps every constant source-traceable. This is the only modeling choice; it is exactly the set named in the source.
2. **ASSUMPTION: Window=28, MinRegion=21, Threshold=0.5 defaults** — window 28 (4 heptads) and 7 registers are from Lupas (1991); MinRegion=21 (3 heptads) from Mason & Arndt's "(abcdefg)1-2-3" multi-heptad requirement; threshold 0.5 = "predominantly hydrophobic" (>half of a/d positions occupied). These are parameters with documented defaults; all are caller-overridable.

---

## Recommendations for Test Coverage

1. **MUST Test:** Perfect heptad (L at every a and d) yields one region with score 1.0 spanning the sequence. — Evidence: definition + Mason & Arndt a/d burial.
2. **MUST Test:** Sequence with no {I,L,V} returns no regions (score 0 < threshold). — Evidence: a/d residue set {I,L,V}.
3. **MUST Test:** Sequence shorter than window returns empty. — Evidence: Lupas window rule.
4. **MUST Test:** Best-register selection — a coiled coil offset from frame 0 is still found via the 7-register max. — Evidence: Lupas "seven different heptad repeat frames".
5. **MUST Test:** Region shorter than MinRegion (21) is rejected even if scoring. — Evidence: multi-heptad requirement.
6. **SHOULD Test:** Half-occupancy (L at a only) gives score exactly 0.5 (boundary). — Rationale: threshold boundary.
7. **SHOULD Test:** Null / empty input returns empty. — Rationale: standard validation.
8. **COULD Test:** Case-insensitivity (lowercase residues recognised). — Rationale: input normalization.

---

## References

1. Mason JM, Arndt KM. 2004. Coiled coil domains: stability, specificity, and biological implications. ChemBioChem 5(2):170–176. https://doi.org/10.1002/cbic.200300781 (PMID 14760737)
2. Lupas A, Van Dyke M, Stock J. 1991. Predicting coiled coils from protein sequences. Science 252(5009):1162–1164. https://doi.org/10.1126/science.252.5009.1162 (PMID 2031185)
3. Chambers P, Pringle CR, Easton AJ. 1990. Heptad repeat sequences are located adjacent to hydrophobic regions in several types of virus fusion glycoproteins. J Gen Virol 71(12):3075–3080. https://doi.org/10.1099/0022-1317-71-12-3075
4. Wikipedia. Coiled coil. https://en.wikipedia.org/wiki/Coiled_coil (accessed 2026-06-14)
5. Wikipedia. Heptad repeat. https://en.wikipedia.org/wiki/Heptad_repeat (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation.
