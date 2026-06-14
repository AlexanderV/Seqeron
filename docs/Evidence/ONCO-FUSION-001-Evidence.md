# Evidence Artifact: ONCO-FUSION-001

**Test Unit ID:** ONCO-FUSION-001
**Algorithm:** Fusion Gene Detection (candidate fusion calling from breakpoint-supporting reads)
**Date Collected:** 2026-06-14

---

## Online Sources

### STAR-Fusion (Haas et al.) — main program source

**URL:** https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
**Accessed:** 2026-06-14 (fetched with WebFetch; the file is the executable `STAR-Fusion` Perl driver)
**Authority rank:** 3 (reference implementation in an established bioinformatics tool)

**Key Extracted Points:**

1. **min_junction_reads default:** the source assigns `my $MIN_JUNCTION_READS = 1;` with help text
   "minimum number of junction-spanning reads required. Default: $MIN_JUNCTION_READS". So the default
   minimum number of junction (split) reads is **1**.
2. **min_sum_frags default:** the source assigns `my $MIN_SUM_FRAGS = 2; # requires at least one junction
   read, else see min_spanning_frags_only below.` Help text: "minimum fusion support = ( # junction_reads
   + # spanning_frags ) Default: $MIN_SUM_FRAGS". So the default minimum total support
   (junction reads + spanning frags) is **2**, and at least one junction read is required to use it.
3. **min_spanning_frags_only default:** `my $MIN_SPANNING_FRAGS_ONLY = 5;` Help text: "minimum number of
   rna-seq fragments required as fusion evidence if there are no junction reads". So when there are
   **zero** junction (split) reads, **5** spanning (discordant) fragments are required.
4. **min_FFPM default:** `my $MIN_FFPM = 0.1;` (fusion fragments per million RNA-seq fragments) — abundance
   normalization filter; recorded for completeness, not used by the rule-based count threshold.
5. **Support is counted as a sum of two evidence classes:** "fusion support = ( # junction_reads +
   # spanning_frags )" — junction reads (split reads crossing the breakpoint) plus spanning fragments
   (discordant pairs that bracket the breakpoint).

### STAR-Fusion preprint (Haas et al. 2017, bioRxiv) — method description

**URL:** https://www.biorxiv.org/content/10.1101/120295.full.pdf
**Accessed:** 2026-06-14 (WebFetch returned 403 on the PDF; method summary obtained from the Genome
Biology benchmark below and from the WebSearch result snippet quoting the abstract)
**Authority rank:** 1 (peer-reviewed companion); recorded as partially retrievable.

**Key Extracted Points:**

1. **Evidence classes:** STAR-Fusion "identifies discordant and split-read alignments, maps them to
   reference transcript structure annotations, filters to remove likely artifacts, and scores them
   according to the abundance of fusion-supporting reads" (search-result abstract snippet). Confirms the
   two evidence classes (split reads + discordant reads) and count-based support scoring.

### Haas et al. 2019, Genome Biology — fusion-detection benchmark

**URL:** https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1842-9
**Accessed:** 2026-06-14 (WebSearch result page describing the paper)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Standard evidence model:** the benchmark (23 methods incl. STAR-Fusion, Arriba) confirms fusion
   callers detect fusions from "discordant and split-read alignments" and score by "the abundance of
   fusion-supporting reads", establishing the split-read + discordant-pair + minimum-support paradigm
   used here.

### Arriba — output-file specification (Uhrig et al.)

**URL:** https://github.com/suhrig/arriba/wiki/05-Output-files
**Accessed:** 2026-06-14 (WebFetch of the GitHub wiki page)
**Authority rank:** 3 (reference-implementation documentation)

**Key Extracted Points:**

1. **split_reads1 / split_reads2 (anchor):** "The number of supporting split fragments with an anchor in
   `gene1` or `gene2`, respectively... The gene to which the longer segment of the split read aligns is
   defined as the anchor." So a split read's anchor is the gene the *longer* aligned segment maps to.
2. **discordant_mates:** "This column contains the number of pairs (fragments) of discordant mates (a.k.a.
   spanning reads or bridge reads) supporting the fusion."
3. **Total supporting reads:** "The total number of supporting reads can be obtained by summing up the
   reads given in the columns `split_reads1`, `split_reads2`, `discordant_mates`..." So total support =
   split_reads1 + split_reads2 + discordant_mates.
4. **reading_frame column:** "in-frame" = gene at 3' end fused in correct frame; "out-of-frame" = fused
   out of correct frame; "stop-codon" = a stop codon exists prior to the junction; "." = peptide sequence
   cannot be predicted.

### Arriba paper (Uhrig et al. 2021, Genome Research) — read-evidence definitions

**URL:** https://pmc.ncbi.nlm.nih.gov/articles/PMC7919457/
**Accessed:** 2026-06-14 (WebFetch of the PMC full-text)
**Authority rank:** 1 (peer-reviewed paper)

**Key Extracted Points:**

1. **Split read:** "split reads, that is, reads with two segments aligning in a noncontiguous fashion".
2. **Discordant mates:** "discordant mates (also referred to as spanning reads or bridge reads), which are
   paired-end reads originating from the same fragment but with the mates aligning in a nonlinear way."

### Genomics England — gene-fusion reading-frame definition

**URL:** https://www.genomicsengland.co.uk/blog/gene-fusion-reporting
**Accessed:** 2026-06-14 (WebFetch of the blog post)
**Authority rank:** 4 (curated explainer; the codon-phase rule it states is the standard exon-phase
definition, cross-checked against the Wikipedia "Reading frame" primary citations below)

**Key Extracted Points:**

1. **In-frame rule (exon phase):** "if one exon finishes after the second letter of a triplet (having end
   phase 2), the next one should start with the third letter (or have start phase 2)" to preserve the
   protein reading frame.
2. **Out-of-frame rule:** when "the reading frame of the downstream partner is changed after the fusion
   point, the fusion is called 'out-of-frame', and it is unlikely to generate a protein."
3. **Goal of in-frame:** "reading through the fusion sequence would still generate a protein consisting of
   parts of either protein encoded by both fusion partners."

### Wikipedia — "Reading frame" (citing primary references)

**URL:** https://en.wikipedia.org/wiki/Reading_frame
**Accessed:** 2026-06-14 (WebFetch)
**Authority rank:** 4 (used only for the primary citations it carries: Badger & Olsen 1999, Mol Biol
Evol 16(4):512–524; Lodish, Molecular Cell Biology 6th ed. p.121)

**Key Extracted Points:**

1. **Codon triplets:** a reading frame is "a specific choice out of the possible ways to read the sequence
   of nucleotides... as a sequence of triplets."
2. **Three frames / modulo 3:** "There are three reading frames that can be read in this 5′→3′ direction,
   each beginning from a different nucleotide in a triplet" — i.e. frame is determined by the start
   position modulo 3. This grounds the in-frame test: the junction is in-frame iff
   (coding bases contributed by the 5' partner up to the breakpoint) ≡ 0 (mod 3) relative to the 3'
   partner's coding start, so downstream codons stay in phase.

---

## Documented Corner Cases and Failure Modes

### From STAR-Fusion source

1. **No junction reads:** if a candidate has zero junction (split) reads, the ordinary `min_sum_frags = 2`
   rule does not apply; instead `min_spanning_frags_only = 5` discordant fragments are required. A
   candidate with only 1–4 discordant fragments and no split reads is filtered out (likely artifact).
2. **Low supporting read count:** candidates below the support thresholds are not reported (controls the
   false-positive rate from low-evidence chimeras).

### From Arriba paper / output spec

1. **Read-through transcripts:** adjacent same-strand neighboring genes produce read-through chimeras that
   are common false-positive "fusions"; callers flag/remove them. (Edge case to guard via support and
   distinct-gene rules.)
2. **Stop codon before junction:** even when the frame is numerically preserved, a stop codon upstream of
   the junction means the 3' partner is not translated ("stop-codon" reading_frame value).

### From general nomenclature

1. **Same gene on both sides:** a "fusion" of a gene with itself is not a gene fusion (Registry invariant
   gene5p ≠ gene3p); such candidates must be rejected.

---

## Test Datasets

### Dataset: Synthetic breakpoint-evidence candidates (derived from cited thresholds)

**Source:** STAR-Fusion defaults (MIN_JUNCTION_READS=1, MIN_SUM_FRAGS=2, MIN_SPANNING_FRAGS_ONLY=5) +
Arriba total-support definition (split_reads1+split_reads2+discordant_mates). No published per-read raw
table is openly downloadable without an alignment run, so candidate-level support counts are constructed
directly from the cited rules; each expected call below is derived from the threshold rule, not from code.

| Candidate | gene5p | gene3p | split1 | split2 | discordant | Total | Junction(=split1+split2) | Expected |
|-----------|--------|--------|--------|--------|-----------|-------|--------------------------|----------|
| C1 EML4-ALK | EML4 | ALK | 3 | 2 | 4 | 9 | 5 | DETECTED (junc≥1, sum≥2) |
| C2 low-junc | BCR | ABL1 | 1 | 0 | 0 | 1 | 1 | DETECTED (junc≥1, sum=... see note) |
| C3 junc=1,sum=2 | TMPRSS2 | ERG | 1 | 0 | 1 | 2 | 1 | DETECTED (junc≥1, sum=2) |
| C4 spanning-only pass | CD74 | ROS1 | 0 | 0 | 5 | 5 | 0 | DETECTED (no junc, span=5≥5) |
| C5 spanning-only fail | NCOA4 | RET | 0 | 0 | 4 | 4 | 0 | REJECTED (no junc, span=4<5) |
| C6 sum<2 | KIF5B | RET | 1 | 0 | 0 | 1 | 1 | REJECTED (junc≥1 but sum=1<2) |
| C7 same gene | ALK | ALK | 5 | 5 | 5 | 15 | 10 | REJECTED (gene5p==gene3p) |

Note on C2/C6: STAR-Fusion requires BOTH min_junction_reads (≥1) AND min_sum_frags (≥2) when junction
reads are present. C2/C6 with junction=1 and total=1 fail min_sum_frags (1 < 2) → REJECTED. C3 with
junction=1 and total=2 passes both → DETECTED.

### Dataset: Reading-frame phase examples (derived from codon-phase rule)

**Source:** Genomics England exon-phase rule + Wikipedia "Reading frame" (modulo-3). In-frame iff the
number of coding bases the 5' partner contributes before the breakpoint, taken modulo 3, equals the 3'
partner's coding-start phase so codons remain in phase.

| Example | 5' coding bases before breakpoint | 3' coding-start phase | (5pBases - 3pPhase) mod 3 | Reading frame |
|---------|-----------------------------------|-----------------------|---------------------------|---------------|
| F1 | 300 | 0 | 0 | In-frame |
| F2 | 301 | 0 | 1 | Out-of-frame |
| F3 | 302 | 0 | 2 | Out-of-frame |
| F4 | 301 | 1 | 0 | In-frame |

---

## Assumptions

1. **ASSUMPTION: Candidate-level input granularity** — The unit consumes already-grouped breakpoint
   candidates with per-class supporting-read counts (split_reads1, split_reads2, discordant_mates),
   not raw BAM records. This mirrors the Arriba output schema (split_reads1/split_reads2/discordant_mates
   columns) and lets the rule-based threshold be tested deterministically. Justification: extracting
   chimeric reads from a BAM is `FindChimericReads` (a separate read-extraction method, out of the
   canonical-threshold scope); the formally defined, source-backed decision rule operates on the counts.

2. **ASSUMPTION: In-frame uses coding-base phase, not stop-codon scanning** — The unit computes in-frame
   status from codon phase ((5pCodingBases - 3pStartPhase) mod 3 == 0) per the exon-phase rule, but does
   not scan the spliced transcript for premature stop codons (Arriba's "stop-codon" value). Stop-codon
   detection requires transcript reconstruction (ONCO-FUSION-003 scope). Out of scope here; documented.

---

## Recommendations for Test Coverage

1. **MUST Test:** A candidate with ≥1 junction read AND total support ≥2 is DETECTED — Evidence:
   STAR-Fusion `MIN_JUNCTION_READS=1`, `MIN_SUM_FRAGS=2`.
2. **MUST Test:** A candidate with 0 junction reads is DETECTED only if discordant fragments ≥5; 4 is
   rejected — Evidence: STAR-Fusion `MIN_SPANNING_FRAGS_ONLY=5`.
3. **MUST Test:** A candidate with junction=1 but total support=1 (<2) is REJECTED — Evidence:
   `MIN_SUM_FRAGS=2`.
4. **MUST Test:** total support = split_reads1 + split_reads2 + discordant_mates — Evidence: Arriba output
   spec.
5. **MUST Test:** gene5p == gene3p candidate is REJECTED (invariant gene5p ≠ gene3p) — Evidence: Registry
   invariant + fusion nomenclature.
6. **MUST Test:** in-frame iff (5' coding bases − 3' start phase) mod 3 == 0 — Evidence: exon-phase rule /
   modulo-3 reading frame.
7. **SHOULD Test:** results ordered by descending total support (most-supported fusions first) — Rationale:
   STAR-Fusion scores "according to the abundance of fusion-supporting reads".
8. **SHOULD Test:** null input → ArgumentNullException; negative counts → ArgumentException — Rationale:
   input validation consistent with sibling OncologyAnalyzer methods.
9. **COULD Test:** custom thresholds override the defaults — Rationale: STAR-Fusion exposes these as CLI
   parameters (`--min_junction_reads`, `--min_sum_frags`, `--min_spanning_frags_only`).

---

## References

1. Haas BJ, Dobin A, Li B, Stransky N, Pochet N, Regev A. (2019). Accuracy assessment of fusion transcript
   detection via read-mapping and de novo fusion transcript assembly-based methods. *Genome Biology*
   20:213. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1842-9
2. Haas BJ, Dobin A, Stransky N, et al. (2017). STAR-Fusion: Fast and Accurate Fusion Transcript Detection
   from RNA-Seq. *bioRxiv* 120295. https://www.biorxiv.org/content/10.1101/120295
3. STAR-Fusion source (defaults MIN_JUNCTION_READS=1, MIN_SUM_FRAGS=2, MIN_SPANNING_FRAGS_ONLY=5,
   MIN_FFPM=0.1). https://raw.githubusercontent.com/STAR-Fusion/STAR-Fusion/master/STAR-Fusion
4. Uhrig S, Ellermann J, Walther T, et al. (2021). Accurate and efficient detection of gene fusions from
   RNA sequencing data. *Genome Research* 31(3):448–460. https://pmc.ncbi.nlm.nih.gov/articles/PMC7919457/
5. Arriba output-file documentation (split_reads1/split_reads2/discordant_mates; reading_frame).
   https://github.com/suhrig/arriba/wiki/05-Output-files
6. Genomics England (2021). Improving how we report gene fusion productivity (in-frame / out-of-frame
   exon-phase rule). https://www.genomicsengland.co.uk/blog/gene-fusion-reporting
7. Wikipedia. Reading frame (citing Badger JH, Olsen GJ. 1999. *Mol Biol Evol* 16(4):512–524; Lodish,
   *Molecular Cell Biology* 6th ed., p.121). https://en.wikipedia.org/wiki/Reading_frame

---

## Change History

- **2026-06-14**: Initial documentation.
