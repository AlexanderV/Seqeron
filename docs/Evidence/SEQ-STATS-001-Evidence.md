# Evidence Artifact: SEQ-STATS-001

**Test Unit ID:** SEQ-STATS-001
**Algorithm:** Sequence Composition Statistics (nucleotide composition, GC content, GC/AT skew)
**Date Collected:** 2026-06-13

---

## Online Sources

### Biopython — `Bio.SeqUtils` source (`gc_fraction`, `GC_skew`)

**URL:** https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
**Accessed:** 2026-06-13 (fetched the raw source file in this session)
**Authority rank:** 3 (reference implementation in an established bioinformatics library)

**Key Extracted Points:**

1. **GC content (`gc_fraction`):** Computes `gc = sum(seq.count(x) for x in "CGScgs")` and `length = gc + sum(seq.count(x) for x in "ATWUatwu")`; returns `gc / length` as a float in [0, 1]. Default `ambiguous="remove"` "will only count GCS and will only include ACTGSWU when calculating the sequence length."
2. **GC content empty handling:** "Note that this will return zero for an empty sequence."
3. **GC skew (`GC_skew`):** `skew = (g - c) / (g + c)` where `g = s.count("G") + s.count("g")` and `c = s.count("C") + s.count("c")`. Default `window=100`.
4. **GC skew zero-denominator handling:** "Returns 0 for windows without any G/C by handling zero division errors" — the code catches `ZeroDivisionError` and sets `skew = 0.0`.

### Wikipedia — "GC skew" (using cited primary, Lobry 1996)

**URL:** https://en.wikipedia.org/wiki/GC_skew
**Accessed:** 2026-06-13 (fetched in this session) · re-fetched 2026-10-14 (SRC citation audit: supports the claims above)
**Authority rank:** 4 (Wikipedia; the formula traces to the Lobry 1996 primary below)

**Key Extracted Points:**

1. **GC skew formula:** `GC skew = (G − C)/(G + C)`.
2. **AT skew formula:** `AT skew = (A − T)/(A + T)`.
3. **Interpretation:** "positive GC skew represents richness of G over C and the negative GC skew represents richness of C over G."
4. **Sign-convention note:** Lobry's original notation was `(C − G)/(C + G)`; modern implementations flip it to `(G − C)/(G + C)`.

### Lobry (1996) — primary source for strand compositional asymmetry

**URL:** https://academic.oup.com/mbe/article-lookup/doi/10.1093/oxfordjournals.molbev.a025626 (DOI: 10.1093/oxfordjournals.molbev.a025626)
**Accessed:** 2026-06-13 (resolved DOI and fetched the article-lookup page in this session)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Citation confirmed:** J. R. Lobry (1996). "Asymmetric substitution patterns in the two DNA strands of bacteria." *Molecular Biology and Evolution* 13(5):660–665.
2. **Finding:** "There was a departure from intrastrand equifrequency between A and T or between C and G, showing that the substitution patterns of the two strands of DNA were asymmetric." This is the basis for measuring G/C and A/T compositional skew.

---

## Documented Corner Cases and Failure Modes

### From Biopython `Bio.SeqUtils`

1. **Empty sequence:** `gc_fraction` returns 0 for an empty sequence.
2. **No G or C present (zero denominator):** `GC_skew` returns 0.0 when `g + c == 0` (caught `ZeroDivisionError`). By the same reasoning AT skew returns 0 when `a + t == 0`.
3. **Mixed case:** GC counting includes lowercase (`"CGScgs"`, `s.count("g")`), so the result is case-insensitive.

---

## Test Datasets

### Dataset: Hand-derived worked examples (from the formulas above)

**Source:** Formulas from Biopython `Bio.SeqUtils` and Wikipedia "GC skew" (cited above).

| Input | A | T | G | C | U | GC content (G+C)/total | AT content | GC skew (G−C)/(G+C) | AT skew (A−T)/(A+T) |
|-------|---|---|---|---|---|------------------------|------------|---------------------|---------------------|
| `ATGC` | 1 | 1 | 1 | 1 | 0 | 2/4 = 0.5 | 2/4 = 0.5 | 0/2 = 0 | 0/2 = 0 |
| `GGGC` | 0 | 0 | 3 | 1 | 0 | 4/4 = 1.0 | 0 | 2/4 = 0.5 | 0 (a+t=0) |
| `AAAT` | 3 | 1 | 0 | 0 | 0 | 0/4 = 0.0 | 4/4 = 1.0 | 0 (g+c=0) | 2/4 = 0.5 |
| `GCCC` | 0 | 0 | 1 | 3 | 0 | 4/4 = 1.0 | 0 | −2/4 = −0.5 | 0 (a+t=0) |
| `AAUUGGCC` | 2 | 0 | 2 | 2 | 2 | 4/8 = 0.5 | 4/8 = 0.5 | 0/4 = 0 | 2/2 = 1.0 (a−t)/(a+t) |

---

## Assumptions

1. **ASSUMPTION: Degenerate IUPAC codes (S, W, R, Y, …) are not counted toward composition totals.** Biopython's `gc_fraction` counts `S` toward GC and `W` toward the denominator. The repository implementation counts only A/T/G/C/U toward GC/AT totals and routes other letters to `CountN`/`CountOther`. For sequences over the standard {A,T,G,C,U} alphabet (the unit's scope) the two agree exactly; the difference only manifests on degenerate symbols. This is a documented repository convention (not a simplification): Biopython `gc_fraction` parity, including S/W, is available from the canonical `CalculateGcFraction(GcAmbiguityMode)`; the amino-acid composition's Biopython `len(seq)` denominator is opt-in since 2026-10-14 (B03 F35, `AminoAcidCompositionDenominator.SequenceLength`).

---

## Recommendations for Test Coverage

1. **MUST Test:** GC content = (G+C)/(A+T+G+C+U) on a known sequence — Evidence: Biopython `gc_fraction`.
2. **MUST Test:** GC skew = (G−C)/(G+C) with exact value including a negative case — Evidence: Wikipedia "GC skew" / Biopython `GC_skew`.
3. **MUST Test:** AT skew = (A−T)/(A+T) with exact value — Evidence: Wikipedia "GC skew".
4. **MUST Test:** Exact A/T/G/C/U/N/Other counts and Length — Evidence: definition of nucleotide composition.
5. **MUST Test:** Empty/null → all-zero composition (GC content 0) — Evidence: Biopython `gc_fraction` empty handling.
6. **SHOULD Test:** Case-insensitivity (upper == lower == mixed) — Rationale: Biopython counts lowercase.
7. **SHOULD Test:** Zero-denominator skew (no G/C → GC skew 0; no A/T → AT skew 0) — Rationale: Biopython `ZeroDivisionError` handling.
8. **COULD Test:** `SummarizeNucleotideSequence` aggregates the same GC content / counts — Rationale: thin summary wrapper.
9. **COULD Test:** `CalculateAminoAcidComposition` returns exact residue counts and length — Rationale: protein-variant wrapper (MW/pI/hydrophobicity belong to SEQ-MW/PI/HYDRO units).

---

## References

1. Lobry, J. R. (1996). Asymmetric substitution patterns in the two DNA strands of bacteria. *Molecular Biology and Evolution* 13(5):660–665. https://doi.org/10.1093/oxfordjournals.molbev.a025626
2. Cock, P. J. A. et al. Biopython, `Bio.SeqUtils` (`gc_fraction`, `GC_skew`). https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py (accessed 2026-06-13)
3. Wikipedia contributors. GC skew. https://en.wikipedia.org/wiki/GC_skew (accessed 2026-06-13)

---

## Addendum 2026-10-14 (FIN-B03 F35, F39, F50) — composition denominator, MCP options

Executed references (copied from B03.md F35/F50):
- Biopython 1.88 `ProteinAnalysis` (aromaticity / `amino_acids_percent`, denominator `len(seq)`): "MKW" 0.33333333333333337 / 33.333333333333336; "MKW*" 0.25 / 25.0; "MKW-A" 0.2 / 20.0; "mkwfy**" 0.4285714285714286 / 14.285714285714286; "FWY12K" 0.5 / 16.666666666666668; "ACDEFGHIKLMNPQRSTVWYX" 0.14285714285714285 / 4.761904761904762.
- EMBOSS 6.6.0 `pepstats` Mole% (denominator = residues): "MKW*" Residues 3, K 33.333; "MKW-A" Residues 4, K 25.000; "MKWX" Residues 4, K 25.000. The library default is the pepstats convention; the earlier wording "unlike Biopython/EMBOSS `len(seq)`" was wrong for EMBOSS.
- API: `CalculateAminoAcidComposition(seq, AminoAcidCompositionDenominator)`, `CalculateAminoAcidFractions(seq, denominator)`; MCP `amino_acid_composition(denominator, includeFractions)` (F50); MCP `nucleotide_composition` reports `n` (F39).
- Tests: `SequenceStatistics_AminoAcidCompositionDenominator_Tests` (18 cases); `Seqeron.Mcp.Sequence.Tests/B03SequenceStatisticsMcpOptionsTests` (F39 N counting, F50 +7).

## Change History

- **2026-06-13**: Initial documentation.
- **2026-09-28** (review campaign B03): added Biopython 1.88 cross-check (`gc_fraction`, `GC_skew`,
  `ProteinAnalysis.aromaticity`/`count_amino_acids`) and EMBOSS pepstats/Eamino.dat class definitions
  (Aromatic F+H+W+Y, Charged B+D+E+H+K+R+Z) for the amino-acid composition ratios; see
  docs/Validation/reports/SEQ-STATS-001.md "Review 2026-09".
- **2026-10-14**: FIN-B03 R2-WP5 — addendum F35/F39/F50; Assumption 1 reworded (convention, not simplification).
