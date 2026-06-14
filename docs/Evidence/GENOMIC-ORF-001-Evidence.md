# Evidence Artifact: GENOMIC-ORF-001

**Test Unit ID:** GENOMIC-ORF-001
**Algorithm:** Open Reading Frame (ORF) Detection — `GenomicAnalyzer.FindOpenReadingFrames`
**Date Collected:** 2026-06-14

---

## Online Sources

### Rosalind — "Open Reading Frames" (ORF) problem

**URL:** https://rosalind.info/problems/orf/
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 4 (educational reference with an exact worked dataset; used here for the worked example, cross-checked against ranks 1–2)

**Key Extracted Points:**

1. **ORF definition:** A sequence in DNA/RNA potentially able to encode a protein that begins with a start codon and ends with a stop codon, containing no internal stop codons.
2. **Reading frames:** A DNA string generates six total reading frames — three from the forward strand and three from its reverse complement.
3. **Start codon:** AUG (codes for methionine, M), indicating the beginning of translation.
4. **Stop codons:** one of three RNA codons indicating termination (UAA/UAG/UGA → DNA TAA/TAG/TGA).
5. **Return value:** all **distinct** candidate protein strings derived by translating ORFs into amino acids until a stop codon is encountered. (Nested/overlapping ORFs that share a stop are both reported — see Test Dataset, e.g. `MGMTPRLGLESLLE` and `MTPRLGLESLLE`.)

### Wikipedia — "Open reading frame" (using its cited statements)

**URL:** https://en.wikipedia.org/wiki/Open_reading_frame
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 4

**Key Extracted Points:**

1. **Definition:** "Reading frames are defined as spans of DNA sequence between the start and stop codons." An ORF may contain a start codon (usually AUG) and by definition cannot extend beyond a stop codon (usually UAA, UAG or UGA in RNA).
2. **Six-frame translation:** "a DNA strand has three distinct reading frames"; with double-stranded DNA "there are six possible frame translations" — three reading frames on each complementary strand.
3. **Minimum length:** "some authors say that an ORF should have a minimal length, e.g. 100 codons or 150 codons." Length divisible by three.

### NCBI ORFfinder (tool reference)

**URL:** https://www.ncbi.nlm.nih.gov/orffinder/
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 5 (NCBI tool)

**Key Extracted Points:**

1. **Start-codon options:** "ATG" only; "ATG" and alternative initiation codons; Any sense codon ("find all stop-to-stop ORFs").
2. **Minimal ORF length:** selectable nucleotide values 30 / 75 / 150 / 300 / 600; the tool restricts the search to ORFs with length **equal or more than** the selected value (nucleotide units).

### NCBI Genetic Codes — Standard code (transl_table=1)

**URL:** https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
**Accessed:** 2026-06-14 (retrieved with WebFetch)
**Authority rank:** 2 (official NCBI specification)

**Key Extracted Points:**

1. **Standard start codon:** ATG (Methionine, M).
2. **Stop codons:** TAA, TAG, TGA — all designated "*" (termination).

---

## Documented Corner Cases and Failure Modes

### From Rosalind / Wikipedia

1. **Nested ORFs sharing a stop:** A downstream ATG in the same frame opens a second, shorter ORF terminated by the same stop; both protein candidates are reported (Rosalind sample: `MGMTPRLGLESLLE` vs `MTPRLGLESLLE`).
2. **No in-frame stop:** A reading begun at an ATG with no subsequent in-frame stop is NOT a complete ORF and yields no protein candidate (Rosalind translates "until a stop codon is encountered").
3. **Both strands:** ORFs must be searched on the reverse complement as well as the forward strand (six-frame search).

### From NCBI ORFfinder

1. **Minimum length filter:** ORFs shorter than the selected nucleotide length are excluded (length ≥ threshold).

---

## Test Datasets

### Dataset: Rosalind ORF sample (`Rosalind_99`)

**Source:** Rosalind ORF problem, https://rosalind.info/problems/orf/ (retrieved 2026-06-14).

| Parameter | Value |
|-----------|-------|
| Input DNA | `AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG` |
| Expected distinct proteins (4) | `MLLGSFRLIPKETLIQVAGSSPCNLS`, `M`, `MGMTPRLGLESLLE`, `MTPRLGLESLLE` |
| Frames searched | 6 (forward + reverse complement, 3 each) |
| Start codon | ATG |
| Stop codons | TAA, TAG, TGA |

Independent re-derivation (Standard code, every-ATG-to-first-in-frame-stop, both strands, distinct proteins) reproduces exactly these 4 proteins.

### Dataset: Single canonical forward ORF (derivation)

**Source:** Standard genetic code (NCBI transl_table=1); construction.

| Parameter | Value |
|-----------|-------|
| Input DNA | `ATG` + `AAA` + `AAA` + `TAA` = `ATGAAAAAATAA` |
| ORF DNA span (incl. stop) | `ATGAAAAAATAA` (12 nt, divisible by 3) |
| Protein candidate (excl. stop) | `MKK` |
| Start position (0-based) | 0, frame 1 |

---

## Assumptions

1. **ASSUMPTION: ORF span includes the terminating stop codon in the reported DNA `Sequence`.** Sources define the ORF as start-codon to stop-codon; the reported nucleotide span includes the stop so that `Length % 3 == 0` and the boundaries are explicit (Wikipedia: "length divisible by three … bounded by stop codons"). The translated protein candidate excludes the stop (Rosalind translates "until a stop codon"). This affects the reported `Sequence`/`Length` only, not which ORFs are detected.
2. **ASSUMPTION: `minLength` is measured in nucleotides** (matching NCBI ORFfinder's nucleotide length options), inclusive lower bound (length ≥ minLength). Default 100 nt retained from the existing public API; it is a caller-supplied parameter and any value is honored.
3. **ASSUMPTION: Standard start codon ATG only** (NCBI ORFfinder default "ATG only"; alternative initiation codons are out of scope for this unit and noted as Not Implemented).

---

## Recommendations for Test Coverage

1. **MUST Test:** Rosalind sample — six-frame, both strands, returns exactly the 4 distinct proteins. — Evidence: Rosalind ORF dataset.
2. **MUST Test:** Single forward ORF `ATGAAAAAATAA` → one ORF, sequence `ATGAAAAAATAA`, position 0, frame 1, protein `MKK`. — Evidence: Standard code derivation.
3. **MUST Test:** Nested ORFs sharing a stop are both reported. — Evidence: Rosalind (`MGM…`/`MTP…`).
4. **MUST Test:** Reading begun at ATG with no in-frame stop yields no ORF. — Evidence: Rosalind / Wikipedia.
5. **MUST Test:** Reverse-complement-only ORF is detected with `IsReverseComplement=true`. — Evidence: Rosalind six-frame requirement.
6. **MUST Test:** minLength filtering excludes shorter ORFs and includes exactly-at-threshold ORFs. — Evidence: NCBI ORFfinder length filter.
7. **MUST Test:** Invariants — every ORF starts with ATG, ends with a stop codon, length divisible by 3. — Evidence: Wikipedia/NCBI definition.
8. **SHOULD Test:** Lowercase input handled (case-insensitive). — Rationale: real-data robustness.
9. **SHOULD Test:** Multiple distinct stop codons (TAA/TAG/TGA) all recognized. — Rationale: completeness of stop set.
10. **COULD Test:** Empty / too-short sequence → empty result. — Rationale: edge guard.

---

## References

1. Rosalind. "ORF — Open Reading Frames." https://rosalind.info/problems/orf/ (accessed 2026-06-14).
2. Wikipedia contributors. "Open reading frame." https://en.wikipedia.org/wiki/Open_reading_frame (accessed 2026-06-14).
3. NCBI. "ORFfinder (Open Reading Frame Finder)." https://www.ncbi.nlm.nih.gov/orffinder/ (accessed 2026-06-14).
4. NCBI. "The Genetic Codes (transl_table=1, Standard)." https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi (accessed 2026-06-14).

---

## Change History

- **2026-06-14**: Initial documentation.
