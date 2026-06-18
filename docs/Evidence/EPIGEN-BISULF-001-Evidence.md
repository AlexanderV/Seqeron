# Evidence Artifact: EPIGEN-BISULF-001

**Test Unit ID:** EPIGEN-BISULF-001
**Algorithm:** Bisulfite Sequencing Analysis (in-silico conversion, methylation calling, profile aggregation)
**Date Collected:** 2026-06-13

---

## Online Sources

### Frommer et al. (1992) — Original bisulfite genomic sequencing protocol

**URL:** https://pubmed.ncbi.nlm.nih.gov/1542678/ (DOI https://doi.org/10.1073/pnas.89.5.1827)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed primary paper, PNAS)
**Retrieved how:** WebSearch `Frommer 1992 PNAS "genomic sequencing protocol that yields a positive display of 5-methylcytosine"`, then WebFetch of the PubMed abstract page.

**Key Extracted Points:**

1. **Conversion chemistry:** The method uses "bisulfite-induced modification of genomic DNA, under conditions whereby cytosine is converted to uracil, but 5-methylcytosine remains nonreactive" (abstract, verbatim).
2. **Read-out base:** Uracil reads/amplifies as thymine; 5-methylcytosine reads as cytosine. The protocol "yields a positive display of 5-methylcytosine residues" — i.e. a remaining C marks methylation.
3. **Single-strand / strand-specific:** the protocol yields "strand-specific sequences of individual molecules"; conversion is applied to one strand at a time (top or bottom), each treated independently.

### Krueger & Andrews (2011) — Bismark methylation caller (peer-reviewed)

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC3102221/ (DOI https://doi.org/10.1093/bioinformatics/btr167)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed) / 3 (reference implementation)
**Retrieved how:** WebFetch of the PMC full-text page.

**Key Extracted Points:**

1. **Methylation call rule:** "The methylation state of positions involving cytosines is determined by comparing the read sequence with the corresponding genomic sequence." A cytosine in the read at a reference-C position indicates methylation (protected from conversion); a thymine at that position indicates an unmethylated cytosine (converted).
2. **Context discrimination:** Bismark discriminates between cytosines in CpG, CHG and CHH context.

### Bismark User Guide (v0.15.0) — methylation extraction percentage formula and call symbols

**URL:** https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.15.0.pdf
**Accessed:** 2026-06-13
**Authority rank:** 3 (reference-implementation documentation)
**Retrieved how:** WebSearch `Bismark methylation extraction percentage methylation 100 * count_methylated / (count_methylated + count_unmethylated)`, then WebFetch of the user-guide PDF.

**Key Extracted Points:**

1. **Percentage methylation formula (verbatim):** "methylation percentage = 100 * count methylated / (count methylated + count unmethylated)". Expressed as a fraction in [0,1] this is `methylated / (methylated + unmethylated)`.
2. **Call symbols:** `z/Z` = CpG (unmethylated/methylated), `x/X` = CHG, `h/H` = CHH; lowercase = unmethylated, uppercase = methylated.
3. **Per-context calculation:** "The percentage methylation is calculated individually for each context following the equation: % methylation (context) = 100 * methylated Cs (context) / (methylated Cs (context) + unmethylated Cs (context))" (Bismark methylation-extractor documentation, retrieved via WebSearch of the babraham bismark page).

### Schultz et al. (2012) — Weighted methylation level

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC6686912/ (cites Schultz et al. 2012, Trends Genet. 28:583–585)
**Accessed:** 2026-06-13
**Authority rank:** 1 (peer-reviewed methods note; verified through a PMC article quoting it verbatim)
**Retrieved how:** WebSearch `Schultz 2012 "weighted methylation level"`, then WebFetch of PMC6686912.

**Key Extracted Points:**

1. **Weighted methylation level (verbatim quote of the definition):** "Weighted DNA methylation was calculated for CG sites by dividing the total number of aligned methylated reads by the total number of methylated plus unmethylated reads" (Schultz et al. 2012). For a per-site fraction f = methylated/coverage, this equals Σ(f·coverage) / Σ(coverage) over the context's sites.

---

## Documented Corner Cases and Failure Modes

### From Frommer et al. (1992)

1. **Non-cytosine bases:** A, G, T are unaffected by bisulfite — only cytosines (5mC stays C, unmethylated C → U→T). Implementation must leave non-C bases unchanged.
2. **Strand specificity:** Conversion operates on a single strand; the complementary strand is a separate molecule. Converting the supplied sequence only (no reverse-complement merge) is the correct top-strand behavior.

### From Bismark User Guide / Krueger & Andrews (2011)

1. **Site with zero coverage:** a cytosine seen by no read has no defined methylation percentage (denominator 0) — such sites are not reported / excluded from percentage.
2. **Reads outside reference / past end:** read bases that fall outside the reference are ignored; a CpG needs both C and G present (last reference base cannot start a CpG).
3. **Reference C that is neither C nor T in a read** (e.g. A/G mismatch): not a valid bisulfite call — neither methylated nor unmethylated; excluded from the count.

---

## Test Datasets

### Dataset: Worked bisulfite conversion (derived from Frommer 1992 chemistry)

**Source:** Frommer et al. (1992) PNAS 89:1827–1831 — conversion rule applied by hand.

| Parameter | Value |
|-----------|-------|
| Input sequence | `ACGTCGAA` (cytosines at index 1 and 4) |
| Methylated positions | `{1}` (the CpG C at index 1 protected) |
| Expected output | `ACGTTGAA` — C@1 protected stays `C`; C@4 unmethylated → `T`; non-C bases unchanged |

### Dataset: Methylation calling from reads (Bismark rule)

**Source:** Bismark User Guide v0.15.0 — `% = 100 * meth / (meth+unmeth)`.

| Parameter | Value |
|-----------|-------|
| Reference | `ACGTACGT` (CpG sites at index 1 and 5) |
| Reads at CpG index 1 | one read with `C` (methylated), one with `T` (unmethylated) |
| Expected level @1 | 1/2 = 0.5, coverage 2 |
| Read at CpG index 5 | one read with `T` (unmethylated) |
| Expected level @5 | 0/1 = 0.0, coverage 1 |

### Dataset: Weighted methylation profile (Schultz 2012)

**Source:** Schultz et al. (2012) Trends Genet. 28:583–585.

| Parameter | Value |
|-----------|-------|
| Site A | CpG, level 1.0, coverage 10 → 10 meth / 10 total |
| Site B | CpG, level 0.0, coverage 30 → 0 meth / 30 total |
| Weighted CpG methylation | (1.0·10 + 0.0·30) / (10+30) = 10/40 = 0.25 |
| Unweighted mean (for contrast) | (1.0 + 0.0)/2 = 0.5 |

---

## Assumptions

1. **ASSUMPTION: Coordinate base** — Positions are 0-based offsets into the supplied sequence (matches sibling `FindCpGSites`/`FindMethylationSites`). Frommer/Bismark are silent on the internal index base; this is an API convention, not correctness-affecting for the values.
2. **ASSUMPTION: Single-strand conversion** — `SimulateBisulfiteConversion` converts only the supplied strand and does not synthesize/merge the complementary strand. Justified by Frommer's strand-specific protocol (each strand is a separate molecule); a two-strand merge is a different operation and out of scope.

---

## Recommendations for Test Coverage

1. **MUST Test:** Unmethylated C → T, methylated C protected, non-C unchanged on a known sequence — Evidence: Frommer 1992 (cytosine→uracil→thymine; 5mC nonreactive).
2. **MUST Test:** Methylation level from reads = meth/(meth+unmeth) at each CpG; coverage = total valid calls; zero-coverage sites excluded — Evidence: Bismark formula `100*meth/(meth+unmeth)`.
3. **MUST Test:** Profile per-context weighted level = Σ(level·coverage)/Σ(coverage) (Schultz) with a worked example giving 0.25 ≠ unweighted 0.5 — Evidence: Schultz 2012.
4. **SHOULD Test:** Empty/null sequence → empty result; lowercase handling; T at CpG counts as unmethylated — Rationale: documented failure modes.
5. **COULD Test:** Reads extending past reference end are ignored; reference A/G at the read position not counted — Rationale: Bismark read-boundary handling.

---

## References

1. Frommer M, McDonald LE, Millar DS, Collis CM, Watt F, Grigg GW, Molloy PL, Paul CL (1992). A genomic sequencing protocol that yields a positive display of 5-methylcytosine residues in individual DNA strands. Proc Natl Acad Sci USA 89(5):1827–1831. https://doi.org/10.1073/pnas.89.5.1827
2. Krueger F, Andrews SR (2011). Bismark: a flexible aligner and methylation caller for Bisulfite-Seq applications. Bioinformatics 27(11):1571–1572. https://doi.org/10.1093/bioinformatics/btr167
3. Babraham Bioinformatics. Bismark Bisulfite Mapper – User Guide v0.15.0 (2016). https://www.bioinformatics.babraham.ac.uk/projects/bismark/Bismark_User_Guide_v0.15.0.pdf
4. Schultz MD, Schmitz RJ, Ecker JR (2012). 'Leveling' the playing field for analyses of single-base resolution DNA methylomes. Trends Genet. 28(12):583–585. https://doi.org/10.1016/j.tig.2012.10.012

---

## Change History

- **2026-06-13**: Initial documentation.
