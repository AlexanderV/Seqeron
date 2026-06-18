# Evidence Artifact: PROTMOTIF-PATTERN-001

**Test Unit ID:** PROTMOTIF-PATTERN-001
**Algorithm:** Protein Pattern Matching Methods (FindMotifByPattern, FindMotifByProsite, ConvertPrositeToRegex, FindDomains)
**Date Collected:** 2026-06-14

---

## Online Sources

### PROSITE Pattern Syntax (ScanProsite documentation / PROSITE User Manual §IV.E)

**URL:** https://prosite.expasy.org/scanprosite/scanprosite_doc.html
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 2 (official database specification)

**Key Extracted Points:**

1. **Amino acid codes:** "The standard IUPAC one letter code for the amino acids is used in PROSITE."
2. **Wildcard:** "The symbol 'x' is used for a position where any amino acid is accepted."
3. **Ambiguous positions (character class):** "[ALT] stands for Ala or Leu or Thr." → maps to regex `[ALT]`.
4. **Excluded positions (negated class):** "{AM} stands for all any amino acid except Ala and Met." → maps to regex `[^AM]`.
5. **Separation:** "Each element in a pattern is separated from its neighbor by a '-'." → separator is dropped in the regex.
6. **Exact repetition:** "x(3) corresponds to x-x-x" and "A(3) corresponds to A-A-A." → maps to `.{3}` / `A{3}`.
7. **Range repetition:** "x(2,4) corresponds to x-x or x-x-x or x-x-x-x." → maps to `.{2,4}`.
8. **Range restriction:** "Ranges can only be used with 'x', for instance 'A(2,4)' is not a valid pattern element." (fixed counts on a residue letter such as `A(3)` are still valid).
9. **Terminal anchors:** "When a pattern is restricted to either the N- or C-terminal of a sequence, that pattern respectively starts with a '<' symbol or ends with a '>' symbol." → `<` maps to `^`, `>` maps to `$`.
10. **Extended ScanProsite query metacharacter:** negative search uses `<{C}*>` (the `*` Kleene star is a ScanProsite query extension, not part of the standard PA-line grammar).

### PROSITE Entry PS00001 (worked example)

**URL:** https://prosite.expasy.org/PS00001
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 5 (curated database entry, with primary references)

**Key Extracted Points:**

1. **Accession/Name/Description:** PS00001, "ASN_GLYCOSYLATION", "N-glycosylation site."
2. **Pattern (PA line):** `N-{P}-[ST]-{P}.` → regex `N[^P][ST][^P]`.

### PROSITE Entry PS00005 (worked example)

**URL:** https://prosite.expasy.org/PS00005
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 5

1. **Accession/Name/Description:** PS00005, "PKC_PHOSPHO_SITE", "Protein kinase C phosphorylation site."
2. **Pattern:** `[ST]-x-[RK]` → regex `[ST].[RK]`.

### PROSITE Entry PS00016 (worked example)

**URL:** https://prosite.expasy.org/PS00016
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 5

1. **Accession/Name/Description:** PS00016, "RGD", "Cell attachment sequence."
2. **Pattern:** `R-G-D.` → regex `RGD` (period terminates the pattern).

### PROSITE Entry PS00017 (worked example)

**URL:** https://prosite.expasy.org/PS00017
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 5

1. **Accession/Name/Description:** PS00017, "ATP_GTP_A", "ATP/GTP-binding site motif A (P-loop)."
2. **Pattern:** `[AG]-x(4)-G-K-[ST]` → regex `[AG].{4}GK[ST]`.

### PROSITE Entry PS00029 (worked example)

**URL:** https://prosite.expasy.org/PS00029
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 5

1. **Accession/Name/Description:** PS00029, "LEUCINE_ZIPPER", "Leucine zipper pattern."
2. **Pattern:** `L-x(6)-L-x(6)-L-x(6)-L` → regex `L.{6}L.{6}L.{6}L`.

### De Castro et al. (2006) — ScanProsite

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC1538847/
**Accessed:** 2026-06-14 (fetched via WebFetch)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Citation:** de Castro E, Sigrist CJA, Gattiker A, Bulliard V, Langendijk-Genevaux PS, Gasteiger E, Bairoch A, Hulo N. "ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins." Nucleic Acids Research. 2006 Jul 14;34(Web Server issue):W362–W365. DOI: 10.1093/nar/gkl124. PMID: 16845026.
2. **Signature kinds:** PROSITE signatures are "patterns (regular expressions), used for short motif detection, or generalized profiles (weight matrices)." (confirms patterns are realized as regular expressions).

### Schneider & Stephens (1990) — Sequence logos / information content

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC332411/ (citation + abstract) and WebSearch result block summarizing the Rsequence formula
**Accessed:** 2026-06-14 (fetched via WebFetch + WebSearch)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Citation:** Schneider TD, Stephens RM. "Sequence logos: a new way to display consensus sequences." Nucleic Acids Res. 1990 Oct 25;18(20):6097–6100. DOI: 10.1093/nar/18.20.6097. PMID: 2172928.
2. **Per-position information content (Rsequence):** Rseq = Smax − Sobs = log2 N − Σ pn log2 pn, where N is the number of distinct symbols and pn is the observed frequency of symbol n. Information is "measured in bits."
3. **Maxima:** "The maximum sequence conservation per site is log2 4 = 2 bits for DNA/RNA and log2 20 ≈ 4.32 bits for proteins." For a position that allows k of the 20 protein residues with uniform probability, IC = log2(20) − log2(k) = log2(20/k) bits.

---

## Documented Corner Cases and Failure Modes

### From PROSITE Pattern Syntax (scanprosite_doc.html)

1. **Period terminator:** A PROSITE PA line ends with a period `.`; characters after the period are not part of the pattern.
2. **Ranges only on `x`:** `x(2,4)` is valid; a range on a residue letter such as `A(2,4)` is *not* a valid PROSITE pattern element (fixed counts like `A(3)` remain valid).
3. **Unsupported ScanProsite query metacharacters:** the `*` Kleene star (`<{C}*>`) belongs to the ScanProsite *query* extension, not the standard PA-line grammar; a PA-line→regex converter must not treat it as a residue.

### From Schneider & Stephens (1990)

1. **Single-residue position contributes maximum IC:** a fixed letter (k=1) contributes log2(20/1) = log2(20) ≈ 4.32 bits; a wildcard `x` (k=20) contributes log2(20/20) = 0 bits.

---

## Test Datasets

### Dataset: PS00016 RGD (literal pattern conversion + match)

**Source:** PROSITE PS00016, https://prosite.expasy.org/PS00016

| Parameter | Value |
|-----------|-------|
| PROSITE pattern | `R-G-D` |
| Expected regex | `RGD` |
| Test sequence | `AAARGDAAA` |
| Expected match start (0-based) | 3 |
| Expected matched substring | `RGD` |
| IC score (3 fixed residues) | 3 × log2(20) ≈ 12.965784284662087 bits |

### Dataset: PS00001 N-glycosylation (exclusion + class conversion + match)

**Source:** PROSITE PS00001, https://prosite.expasy.org/PS00001

| Parameter | Value |
|-----------|-------|
| PROSITE pattern | `N-{P}-[ST]-{P}` |
| Expected regex | `N[^P][ST][^P]` |
| Test sequence | `AANASAAANGTAAAA` |
| Expected match starts (0-based) | 2, 8 |
| Expected matched substrings | `NASA`, `NGTA` |

### Dataset: PS00017 P-loop (fixed-range x(4) conversion)

**Source:** PROSITE PS00017, https://prosite.expasy.org/PS00017

| Parameter | Value |
|-----------|-------|
| PROSITE pattern | `[AG]-x(4)-G-K-[ST]` |
| Expected regex | `[AG].{4}GK[ST]` |

### Dataset: Information content scoring (Schneider & Stephens 1990)

**Source:** Schneider & Stephens (1990), DOI 10.1093/nar/18.20.6097

| Pattern position | Allowed residues k | IC = log2(20/k) bits |
|------------------|--------------------|----------------------|
| Fixed letter (e.g. `R`) | 1 | log2(20) ≈ 4.321928094887363 |
| Class `[ST]` | 2 | log2(10) ≈ 3.321928094887362 |
| Wildcard `.` (`x`) | 20 | 0 |

---

## Assumptions

1. **ASSUMPTION: Overlapping-match enumeration.** `FindMotifByPattern` wraps the regex in a zero-width lookahead `(?=(...))` so it enumerates *all* start positions, including overlapping ones. The PROSITE/ScanProsite papers retrieved this session do not state whether overlapping occurrences are reported. This is a documented repository contract, not a PROSITE-mandated behavior; it does not change the set of *positions where the pattern can start*, only whether overlapping starts are all listed. Marked ASSUMPTION; it is verified by tests as a stable contract but not claimed to be PROSITE-specified.

2. **ASSUMPTION: E-value model.** `CalculateEValue` uses E = (N − L + 1) × 2^(−IC), the expected number of random matches under a uniform i.i.d. amino-acid background. This is a standard combinatorial expectation consistent with the IC definition (Schneider & Stephens 1990) but the exact E-value reported by ScanProsite (which uses Swiss-Prot residue frequencies) was not retrieved; the repository value is therefore a model-defined quantity, tested for its defining formula rather than against a ScanProsite number.

---

## Recommendations for Test Coverage

1. **MUST Test:** `ConvertPrositeToRegex` produces the exact regex for each retrieved worked example (PS00001, PS00005, PS00016, PS00017, PS00029) and for each syntax atom (x, x(n), x(n,m), [..], {..}, A(n), <, >, trailing `.`). — Evidence: scanprosite_doc.html; PS00001/05/16/17/29.
2. **MUST Test:** `FindMotifByPattern` returns exact 0-based start/end/substring for a literal and a class pattern, and `Score` equals the information content Σ log2(20/k_i). — Evidence: Schneider & Stephens (1990); PROSITE entries.
3. **MUST Test:** `FindMotifByProsite` end-to-end (PROSITE → regex → matches) on PS00001 and PS00016 yields the dataset positions. — Evidence: PS00001, PS00016.
4. **MUST Test:** `FindDomains` detects the P-loop / kinase ATP-binding (PS00017-equivalent) signature at the correct position in a synthetic sequence. — Evidence: PS00017.
5. **SHOULD Test:** overlapping enumeration via lookahead; unsupported `*` rejected with FormatException. — Rationale: documented repository contract / "reject don't silently drop".
6. **COULD Test:** null/empty inputs return empty (no exception) for every method. — Rationale: defensive contract.

---

## References

1. Hulo N, Bairoch A, Bulliard V, Cerutti L, Cuche BA, de Castro E, Lachaize C, Langendijk-Genevaux PS, Sigrist CJA. (2008). The 20 years of PROSITE. Nucleic Acids Res. 36(Database):D245–D249. https://doi.org/10.1093/nar/gkm977
2. de Castro E, Sigrist CJA, Gattiker A, Bulliard V, Langendijk-Genevaux PS, Gasteiger E, Bairoch A, Hulo N. (2006). ScanProsite: detection of PROSITE signature matches and ProRule-associated functional and structural residues in proteins. Nucleic Acids Res. 34(Web Server):W362–W365. https://doi.org/10.1093/nar/gkl124
3. Schneider TD, Stephens RM. (1990). Sequence logos: a new way to display consensus sequences. Nucleic Acids Res. 18(20):6097–6100. https://doi.org/10.1093/nar/18.20.6097
4. ExPASy / PROSITE. PROSITE pattern syntax (ScanProsite documentation). https://prosite.expasy.org/scanprosite/scanprosite_doc.html (accessed 2026-06-14)
5. PROSITE entries PS00001, PS00005, PS00016, PS00017, PS00029. https://prosite.expasy.org/PS00001 etc. (accessed 2026-06-14)

---

## Change History

- **2026-06-14**: Initial documentation for PROTMOTIF-PATTERN-001.
