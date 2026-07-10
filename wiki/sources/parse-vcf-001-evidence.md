---
type: source
title: "Evidence: PARSE-VCF-001 (VCF parsing — 8 fixed columns + FORMAT/genotype samples, 1-based)"
tags: [validation, file-io]
doc_path: docs/Evidence/PARSE-VCF-001-Evidence.md
sources:
  - docs/Evidence/PARSE-VCF-001-Evidence.md
source_commit: 7dbb2946c73ee25f13cab78aa33fee223dfa48e1
ingested: 2026-07-10
created: 2026-07-10
updated: 2026-07-10
---

# Evidence: PARSE-VCF-001

The validation-evidence artifact for test unit **PARSE-VCF-001** — **VCF (Variant Call Format)
parsing**: read a variant file — `##`-prefixed metadata headers plus tab-delimited records of
**eight fixed columns** (`CHROM` / `POS` / `ID` / `REF` / `ALT` / `QUAL` / `FILTER` / `INFO`),
optionally followed by a `FORMAT` column and one column per genotyped sample — into typed
variant records, classify each variant, and compute site/genotype statistics. This is a
**FileIO** (file-parsing `PARSE-*`) family Evidence file and one instance of the templated
per-algorithm [[algorithm-validation-evidence|evidence artifact]] pattern. Like
[[parse-gff-001-evidence]] it is a **tab-delimited, 1-based-coordinate** format (not an INSDC
flat file such as [[parse-genbank-001-evidence]] / [[parse-embl-001-evidence]]), so it shares
the family anchor [[bed-format-parsing]]; the **coordinate contrast** that matters here — VCF
`POS` is **1-based** whereas BED is 0-based half-open — is already spelled out on that anchor
page (BED explicitly names GFF/GTF/VCF as its 1-based counterparts). See [[test-unit-registry]]
for how units are tracked and [[fuzzing]] for why parsers are the family's hottest
malformed-input target. No dedicated concept page is warranted: VCF's genotype/INFO model and
classification rules are captured concisely below (economical per the ingest directive), and
the coordinate model lives on the shared anchor.

## What this file records

- **Authoritative / online sources:**
  - **Wikipedia — Variant Call Format** (accessed 2026-02-05) — format overview (a text format
    for gene-sequence variations, developed 2010 for the 1000 Genomes Project, current
    **VCFv4.5** Oct 2024); the eight-column specification; header-line grammar; the common
    INFO / FORMAT fields; variant-classification and genotype-notation rules.
  - **Danecek et al. (2011), "The variant call format and VCFtools"**, *Bioinformatics*
    27(15):2156–2158 (DOI 10.1093/bioinformatics/btr330, PMID 21653522) — the authoritative
    paper defining the format and its parsing/validation conventions, including missing-value
    (`.`) semantics.
  - **SAMtools/GA4GH HTS-Specs — VCFv4.3 specification** (`samtools.github.io/hts-specs`) — the
    strict parsing rules, symbolic-allele / breakend (§5.4) and spanning-deletion (§5.5) grammar,
    and error handling for malformed records.
  - **EMBL-EBI VCF training** — usage context: VCF as the standard input to variant-analysis
    tools (VEP) and standard output from callers (GATK).

- **Columns (1-based coordinates):**
  - **8 mandatory:** `CHROM` (sequence name), `POS` (**1-based** position), `ID` (variant id /
    rsID or `.`), `REF` (reference base(s)), `ALT` (comma-separated alternate alleles), `QUAL`
    (Phred-scaled quality or `.`), `FILTER` (`PASS`, or semicolon-separated failed filters, or
    `.`), `INFO` (semicolon-separated key–value pairs).
  - **Optional per-sample block:** column 9 `FORMAT` (colon-separated field keys) plus one
    column per sample carrying colon-separated values in `FORMAT` order.
  - **Common INFO fields:** `DP` (total depth), `AF` (allele frequency), `NS` (samples with
    data), `DB` (dbSNP membership — a flag), `AA` (ancestral allele). Flags with no value are
    stored as `key="true"`.
  - **Common FORMAT fields:** `GT` (genotype), `DP` (read depth), `AD` (allele depth), `GQ`
    (genotype quality), `PL` (Phred-scaled genotype likelihoods).

- **Header lines:** `##`-prefixed metadata — `##fileformat=VCFvX.X`, `##INFO=<…>`,
  `##FORMAT=<…>`, `##FILTER=<…>` structured definitions — and the single `#CHROM …` column-header
  line that names the fixed columns and any samples.

- **Variant classification (Wikipedia / spec):** **SNP** = single-base REF and ALT that differ;
  **MNP** = equal-length (>1) REF/ALT that differ; **Insertion** = ALT longer than REF; **Deletion**
  = REF longer than ALT; **Symbolic** = ALT begins with `<` (e.g. `<DEL>`, `<INS>`, `<DUP>`,
  `<INV>`, `<CNV>`), breakend notation (`]p]t` / `[p[t` / `t[p[` / `t]p]`), or the spanning-deletion
  placeholder `*`.

- **Genotype notation:** `0/0` hom-ref, `0/1`/`1/0` het, `1/1` hom-alt, `1/2` both-alleles-alt
  (multi-allelic); `|` phased vs `/` unphased; `.` = missing allele.

## Spec-compliance points (the theory this unit pins)

The artifact's 2026-03-11 audit corrected five deviations to strict spec compliance — these are
the load-bearing correctness facts:

1. **`FILTER "."` ≠ `PASS`.** `PASS` = the record was evaluated and passed all filters; `.` = no
   filtering has been applied (record unevaluated). Distinct states per VCF 4.3 + GATK.
2. **Missing allele `.` in `GT` → zygosity indeterminate** (Danecek 2011: `.` = missing data);
   `IsHet` returns false when any allele is missing (neither het nor hom).
3. **Ti/Tv iterates all ALT alleles** at multi-allelic sites, not just the first (transitions
   A↔G / C↔T; everything else a transversion). Whole-genome Ti/Tv ≈ 2.0–2.1, exome ≈ 3.0.
4. **`*` spanning deletion (VCF 4.3 §5.5)** = an upstream deletion spanning this position;
   classified **Symbolic**, not SNP/indel.
5. **All four breakend forms (§5.4)** detected via the `[` / `]` characters.

Plus: `QUAL "."` → null; `ID "."` → no identifier; `INFO "."` → no fields; allele frequency =
`alt_count / total_alleles`.

## Documented edge cases and malformed-record handling

- **Missing values** `.` across `QUAL` / `ID` / `FILTER` / `INFO` handled per the rules above.
- **Multi-allelic** `ALT` = `A,C,G` (each a distinct alternate at the same position);
  **complex filters** `LowQual;LowCov`; **missing-allele genotypes** `./0`, `0/.`, `./.`.
- **Malformed records skipped gracefully:** a non-integer `POS`, or fewer than 8 columns → the
  record is dropped (parser returns null for the line) rather than throwing.

## Testing methodology

Required categories mirror the format and its downstream use: parsing (empty/null → empty;
header/metadata extraction; per-sample genotype parsing); variant classification
(SNP/MNP/Ins/Del/Symbolic); filtering (by chromosome, region, quality threshold, `PASS`-only,
variant type); genotype analysis (hom-ref/hom-alt/het, phased vs unphased); statistics (counts
by type, chromosome distribution, **Ti/Tv ratio**, passing-variant counts); and I/O
(file parse + round-trip write/re-read). **Reference-data tests** validate against real
**1000 Genomes** (common SNPs, multi-allelic sites, AF/DP/NS INFO fields) and **ClinVar**
(pathogenic variants, `CLNSIG`/`CLNDN`, symbolic structural variants) records.

## Invariants (source-derived)

`POS` is a 1-based positive integer; `QUAL` non-negative when present else null;
`|len(ALT) − len(REF)|` gives indel length; Ti/Tv accounts for all ALT alleles; missing-allele
genotypes have indeterminate zygosity; `*` = spanning deletion (Symbolic).

## Deviations and assumptions

The artifact records **no residual deviations** after the 2026-03-11 audit ("No deviations. No
tests tailored to implementation quirks."). The implementation follows the **VCFv4.3**
specification; the five previously-flagged points above were corrected *to* the spec rather than
away from it. The only implementation-shaped behaviours are API-contract choices outside the
format spec: **malformed lines (bad `POS` or < 8 columns) are skipped, not raised**, and valueless
**INFO flags are stored as `key="true"`**. No source contradictions.
