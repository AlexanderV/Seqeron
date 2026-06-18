# Evidence Artifact: EPIGEN-METHYL-001

**Test Unit ID:** EPIGEN-METHYL-001
**Algorithm:** Methylation Site Detection, Sequence-Context Classification (CpG/CHG/CHH), and Methylation Profile
**Date Collected:** 2026-06-13

---

## Online Sources

### IUPAC Nucleotide Ambiguity Codes (Cornish-Bowden 1985) — Los Alamos HIV DB table

**URL:** https://www.hiv.lanl.gov/content/sequence/HelpDocs/IUPAC.html
**Accessed:** 2026-06-13 (WebFetch of the URL above)
**Authority rank:** 2 (official IUPAC-IUB nomenclature, tabulated)

**Key Extracted Points:**

1. **Code H:** The table states `H = A or C or T` — i.e. "not G". Companion rows: `B = C or G or T`, `V = A or C or G`, `D = A or G or T`, `N = G or A or T or C`.
2. **Source citation:** The page cites "Cornish-Bowden (1985) IUPAC-IUB SYMBOLS FOR NUCLEOTIDE NOMENCLATURE. Nucl. Acids Res. 13: 3021-3030" as the origin of the codes.

### Krueger & Andrews (2011) — Bismark, Bioinformatics 27(11):1571–1572 (PMC3102221)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3102221/
**Accessed:** 2026-06-13 (WebFetch of the PMC article)
**Authority rank:** 3 (reference implementation / peer-reviewed methods paper)

**Key Extracted Points:**

1. **Context discrimination:** "methylation calls in Bismark take the surrounding sequence context into consideration and discriminate between cytosines in CpG, CHG and CHH context."
2. **Meaning of H:** "whereby H can be either A, T or C."
3. **Symmetry:** "In plants, methylation is quite common in both the symmetric CpG or CHG, and asymmetric CHH context." The strand-specific output is "very useful to study asymmetric methylation (hemi- or CHH methylation) in a strand-specific manner."

### Lister et al. (2009) — Human DNA methylomes at base resolution, Nature 462:315–322

**URL:** https://pubmed.ncbi.nlm.nih.gov/19829295/ (DOI 10.1038/nature08514)
**Accessed:** 2026-06-13 (WebSearch returning the abstract + PubMed WebFetch)
**Authority rank:** 1 (peer-reviewed primary literature)

**Key Extracted Points:**

1. **Non-CG contexts defined:** Methylation in H1 embryonic stem cells occurs "in non-CG contexts (mCHG and mCHH, where H = A, C or T)."
2. **Prevalence H1:** Non-CG methylation (mCHG + mCHH) comprises "almost 25% of all cytosines at which DNA methylation is identified" in H1 ES cells ("Nearly one-quarter of all methylation identified in embryonic stem cells was in a non-CG context").
3. **Prevalence IMR90:** "of the methylcytosines detected in the IMR90 genome, 99.98% were in the CG context" — i.e. somatic/differentiated cells are essentially CpG-only.

### Schultz, Schmitz & Ecker (2012) — 'Leveling' the playing field, Trends in Genetics 28(12):583–585

**URL:** https://www.cell.com/trends/genetics/abstract/S0168-9525(12)00171-0 (DOI 10.1016/j.tig.2012.10.012)
**Accessed:** 2026-06-13 (WebSearch of the paper text; methods reproduced in PMC6686912, fetched)
**Authority rank:** 1 (peer-reviewed primary literature defining the metric)

**Key Extracted Points:**

1. **Weighted methylation level:** "Weighted methylation is the sum of the methylated reads divided by the number of total reads (methylated and un-methylated) for a certain cytosine sequence context" — i.e. WML = (Σ methylated reads) / (Σ total reads) over all cytosines of the context.
2. **Per-context aggregation:** The metric is computed per sequence context (CG, CHG, CHH separately), motivating separate CpG/CHG/CHH summaries.
3. **Reproduced in methods (PMC6686912, fetched 2026-06-13):** "Weighted DNA methylation was calculated for CG sites by dividing the total number of aligned methylated reads by the total number of methylated plus unmethylated reads (Schultz et al. 2012)."

### Per-cytosine fractional methylation level (Schultz 2012 / Lister 2009 convention) — PMC6686912

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6686912/
**Accessed:** 2026-06-13 (WebFetch of the PMC article)
**Authority rank:** 3 (peer-reviewed methods reproducing Schultz 2012)

**Key Extracted Points:**

1. **Per-cytosine level:** Methylation level at a single cytosine = methylated reads / total (methylated + unmethylated) reads at that position — a fraction in [0, 1]. This is the per-site `MethylationLevel` value carried by each site.

---

## Documented Corner Cases and Failure Modes

### From IUPAC (Cornish-Bowden 1985) / Bismark (2011)

1. **H excludes G:** A trinucleotide `C x G` where x ∈ {A,C,T} is CHG, but `C G ...` (x = G) is CpG, never CHG. The middle/last H position must be A, C, or T; a G there changes the context class.
2. **CHH is asymmetric:** `C H H` has no complementary symmetric counterpart, unlike CpG and CHG which are palindromic/symmetric across the two strands.

### From sequence-context geometry (Bismark)

1. **Incomplete trailing context:** A cytosine at sequence end whose downstream context is truncated cannot be classified into CHG vs CHH (the third base is unknown); a `CG` at the end is still an unambiguous CpG because only two bases are needed.
2. **Non-ACGT / ambiguous bases:** If the base after C, or the two downstream bases, are not in {A,C,G,T}, the context is undetermined and the cytosine is not a classifiable methylation site.

---

## Test Datasets

### Dataset: Hand-derived context-classification sequence

**Source:** Derived directly from the context definitions (Bismark 2011; Lister 2009; IUPAC H = A/C/T).

| Parameter | Value |
|-----------|-------|
| Sequence | `CGACAGCAA` |
| C at index 0 | `CG…` → CpG |
| C at index 3 | `CAG` → CHG (H=A, then G) |
| C at index 6 | `CAA` → CHH (H=A, H=A) |
| Total classifiable C sites | 3 |

### Dataset: Weighted methylation level worked example

**Source:** Schultz et al. (2012) definition WML = Σ(methylated)/Σ(total).

| Parameter | Value |
|-----------|-------|
| Site A | methylated=8, total=10 → fraction 0.8 |
| Site B | methylated=2, total=10 → fraction 0.2 |
| Weighted CpG level | (8+2)/(10+10) = 10/20 = 0.5 |
| Unweighted mean of fractions | (0.8+0.2)/2 = 0.5 (equal here; differs when coverage unequal) |

### Dataset: Lister (2009) biological prevalence (context realism)

**Source:** Lister et al. (2009), Nature 462:315–322.

| Parameter | Value |
|-----------|-------|
| IMR90 (somatic) methylcytosines in CG context | 99.98% |
| H1 ES cells methylcytosines in non-CG (CHG+CHH) | ~25% |

---

## Assumptions

1. **ASSUMPTION: Per-cytosine level for sequence-only input is 0.** `FindMethylationSites(sequence)` has no bisulfite read data, so each site's `MethylationLevel` and `Coverage` are 0 until measured. This is a representational default (the site is *potentially* methylatable); it is not a claim from the sources that the level is zero. Used only by the sequence-only entry point; methylation values come from `CalculateMethylationFromBisulfite` (out of this unit's scope) or caller-supplied sites.
2. **ASSUMPTION: "Methylated CpG" summary threshold in the profile (level ≥ 0.5).** The `MethylatedCpGSites` count uses a 0.5 binary cutoff. Schultz (2012) recommends a binomial test rather than a fixed fraction cutoff; 0.5 is a descriptive convenience for the count field and does not affect the continuous methylation-level outputs, which follow Schultz (2012) exactly.

---

## Recommendations for Test Coverage

1. **MUST Test:** `GetMethylationContext` returns CpG for `CG`, CHG for `C`+H+`G`, CHH for `C`+H+`H`, with H ∈ {A,C,T} — Evidence: Bismark (2011), IUPAC H = A/C/T.
2. **MUST Test:** H position rejecting G — `CGG` middle is G → CpG not CHG; classification driven by exact bases — Evidence: IUPAC (G ∉ H).
3. **MUST Test:** `FindMethylationSites` enumerates one site per classifiable C with correct Type and 0-based Position on `CGACAGCAA` — Evidence: context geometry.
4. **MUST Test:** Trailing/incomplete context not classified; `CG` at end still CpG — Evidence: Bismark context geometry.
5. **MUST Test:** Profile weighted/mean methylation per context and counts on caller-supplied sites with known levels — Evidence: Schultz (2012) WML formula.
6. **SHOULD Test:** Case-insensitivity (lowercase input) — Rationale: sequences are commonly lowercase-masked.
7. **SHOULD Test:** Null/empty sequence → empty; empty site list → zeroed profile — Rationale: documented validation contract.
8. **COULD Test:** Non-ACGT base in H/third position → not classified — Rationale: undetermined context per IUPAC.

---

## References

1. Cornish-Bowden A. (1985). Nomenclature for incompletely specified bases in nucleic acid sequences: recommendations 1984 (IUPAC-IUB). Nucleic Acids Research 13(9):3021–3030. https://doi.org/10.1093/nar/13.9.3021 (table retrieved at https://www.hiv.lanl.gov/content/sequence/HelpDocs/IUPAC.html)
2. Krueger F, Andrews SR. (2011). Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics 27(11):1571–1572. https://doi.org/10.1093/bioinformatics/btr167 (https://pmc.ncbi.nlm.nih.gov/articles/PMC3102221/)
3. Lister R, Pelizzola M, Dowen RH, et al. (2009). Human DNA methylomes at base resolution show widespread epigenomic differences. Nature 462(7271):315–322. https://doi.org/10.1038/nature08514 (https://pubmed.ncbi.nlm.nih.gov/19829295/)
4. Schultz MD, Schmitz RJ, Ecker JR. (2012). 'Leveling' the playing field for analyses of single-base resolution DNA methylomes. Trends in Genetics 28(12):583–585. https://doi.org/10.1016/j.tig.2012.10.012 (https://www.cell.com/trends/genetics/abstract/S0168-9525(12)00171-0)

---

## Change History

- **2026-06-13**: Initial documentation.
