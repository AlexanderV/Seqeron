---
type: concept
title: "VCF parsing (FileIO VcfParser — 8 fixed columns + FORMAT/genotype samples, variant classification + Ti/Tv)"
tags: [file-io, algorithm]
sources:
  - docs/algorithms/FileIO/VCF_Parsing.md
  - docs/Evidence/PARSE-VCF-001-Evidence.md
source_commit: a8fa2abd35e44fbdff855a9649086c90dcbb8179
created: 2026-07-13
updated: 2026-07-13
graph:
  relationships:
    - predicate: relates_to
      object: concept:test-unit-registry
      source: parse-vcf-001-evidence
      evidence: "Test Unit ID: PARSE-VCF-001 | Algorithm Group: FileIO (VcfParser)"
      confidence: high
      status: current
---

# VCF parsing (FileIO VcfParser — 8 fixed columns + FORMAT/genotype samples, variant classification + Ti/Tv)

**VCF (Variant Call Format)** is the tab-delimited text format for genomic variation exchange: a
`##`-prefixed metadata **header** section followed by a data section carrying **one record per
variant locus**. Each data row has **eight fixed columns** (`CHROM` / `POS` / `ID` / `REF` / `ALT` /
`QUAL` / `FILTER` / `INFO`), optionally followed by a `FORMAT` column and one column per genotyped
sample. Seqeron's FileIO **`VcfParser`** (test unit **PARSE-VCF-001**, status *Simplified*) parses
this text from files and strings into typed `VcfRecord` values, optionally returns a structured
`VcfHeader`, then layers variant classification, filtering, genotype/sample inspection, summary
statistics, transition/transversion ratios, INFO helpers, and record writing on top.

This concept is the **parser-surface** synthesis of the primary algorithm spec
`docs/algorithms/FileIO/VCF_Parsing.md` — it owns the `VcfParser` contract, the line-oriented parse,
and the classification/statistics helpers. It deliberately **does not re-derive** the
literature-traced format facts (the column specification, the common INFO/FORMAT field tables, the
symbolic/breakend/spanning-deletion grammar, the 1000 Genomes / ClinVar reference-data oracle, and
the five 2026-03-11 spec-compliance corrections) — those live on the Evidence source page
[[parse-vcf-001-evidence]] and are **cross-linked, not duplicated** here. It is a member of the
FileIO file-parsing **`PARSE-*`** family anchored by [[bed-format-parsing]] (where the **coordinate
contrast** — VCF `POS` is 1-based, BED is 0-based half-open — is spelled out); its tab-delimited,
1-based sibling is [[gff-parsing]], and its format siblings are [[fasta-parsing]] and
[[fastq-parsing]]. [[test-unit-registry]] tracks the unit, [[algorithm-validation-evidence]]
describes the evidence-artifact pattern, and [[fuzzing]] explains why parsers are the campaign's
highest-priority malformed-input target.

## The format (what the parser consumes)

```text
##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=q10,Description="Quality below 10">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1
```

Coordinates are **1-based**. The 8 mandatory columns are fixed; column 9 (`FORMAT`, colon-separated
field keys) and the following per-sample columns are optional and parsed **only when the `#CHROM`
header line defines sample names**. The full column / common-INFO / common-FORMAT / genotype-notation
tables live on [[parse-vcf-001-evidence]].

### Invariants (from the spec)

| ID | Invariant |
|----|-----------|
| INV-01 | `POS` is **1-based** (VCF coordinates are 1-based by specification). |
| INV-02 | VCF data lines have **at least 8** mandatory tab-separated columns. |
| INV-03 | `PASS` and `.` are **different** FILTER states — `PASS` = evaluated and passed; `.` = no filtering applied. |

## The parser (line-oriented, four entry points)

`VcfParser` exposes four read entry points that split into *records-only* and *records + structured
header* pairs:

| Entry point | Returns | Null / empty / missing behavior |
|-------------|---------|----------------------------------|
| `Parse(content)` | records | null or empty content → **no records** (early-return guard) |
| `ParseWithHeader(content)` | records + `VcfHeader` | **null throws** (`StringReader` needs non-null); empty content → default `VCFv4.3` header + no records |
| `ParseFile(filePath)` | records | null / empty / missing path → **no records** |
| `ParseFileWithHeader(filePath)` | records + `VcfHeader` | missing path → default `VCFv4.3` header, empty collections, no records |

Parse steps:

1. Read `##` header-metadata lines; promote structured **INFO / FORMAT / FILTER** definitions into
   specialized `VcfHeader` structures and capture every other `##key=value` line as raw
   `OtherMetadata`.
2. Parse the `#CHROM` line to determine **sample names** (this is what enables per-sample parsing).
3. **Split each data line on tabs** and parse the first 8 mandatory columns.
4. Parse optional `FORMAT`/sample columns **only when sample names are available**.
5. Expose the classification / filtering / genotype / statistics / Ti-Tv / writing / INFO helpers.

### Sentinel normalization (`.` handling)

| Field | `.` becomes | Rationale |
|-------|-------------|-----------|
| `QUAL` | **`null`** (`double?`) | missing quality score |
| `FILTER` | **empty array** | distinct state from `PASS` (INV-03) |
| `ID` | no identifier | missing variant id |
| INFO flag with **no `=`** | stored as value **`"true"`** | flag-style INFO handling |

## Contract & surface (`VcfParser`, `VcfRecord`, `VcfHeader`)

Implementation: `src/Seqeron/Algorithms/Seqeron.Genomics.IO/VcfParser.cs`.

| Group | Entry points | Notes |
|-------|--------------|-------|
| **Parse** | `Parse`, `ParseWithHeader`, `ParseFile`, `ParseFileWithHeader` | string / file × records / records+header. |
| **Classify** | `ClassifyVariant`, `IsSNP`, `IsIndel`, `GetVariantLength` | per-allele, `O(1)`. |
| **Filter** | `FilterByChrom`, `FilterByRegion`, `FilterByQuality`, `FilterPassing`, `FilterByType`, `FilterSNPs`, `FilterIndels`, `FilterByInfo` | `O(r)` scan, `O(m)` matches. |
| **Genotype / sample** | `GetGenotype`, `IsHomRef`, `IsHomAlt`, `IsHet`, `GetAlleleDepth`, `GetReadDepth`, `CalculateAlleleFrequency` | `GT` split on `/` or `|`; out-of-range `sampleIndex` → null/false. |
| **Statistics / write** | `CalculateStatistics`, `CalculateTiTvRatio`, `WriteToStream`, `WriteToFile` | one pass over records; writer normalizes formatting. |
| **INFO access** | `GetInfoValue`, `GetInfoInt`, `GetInfoDouble`, `HasInfoFlag` | typed INFO lookups over stringly-typed values. |

`VcfRecord` exposes `Chrom` (`string`), `Pos` (`int`, 1-based), `Id` (`string`), `Ref` (`string`),
`Alt` (`string[]`, multi-allelic), `Qual` (`double?`), `Filter` (`string[]`), `Info`
(`IReadOnlyDictionary<string,string>`, flags as `"true"`), `Format` (`string[]?`), and `Samples`
(`IReadOnlyList<IReadOnlyDictionary<string,string>>?`, dictionaries keyed by FORMAT entries).
`VcfHeader` carries `fileformat`, structured INFO/FORMAT/FILTER, sample names, and `OtherMetadata`.

### Complexity

| Operation | Time | Space |
|-----------|------|-------|
| `Parse` / `ParseWithHeader` | `O(n)` | `O(n)` (n = text size) |
| `ClassifyVariant` / `GetVariantLength` | `O(1)` | `O(1)` |
| `FilterBy*` / `FilterPassing` | `O(r)` | `O(m)` (m = matches) |
| `CalculateStatistics` / `CalculateTiTvRatio` | `O(r)` | `O(1)` auxiliary |

## Variant classification (`ClassifyVariant`)

Classification is by REF/ALT length with a symbolic override, aligned to repository behavior:

| Type | Condition | Example |
|------|-----------|---------|
| `SNP` | `len(REF)=1` and `len(ALT)=1` | `A → G` |
| `MNP` | `len(REF)=len(ALT)>1` | `AT → GC` |
| `Insertion` | `len(REF)=1` and `len(ALT)>1` | `A → ATG` |
| `Deletion` | `len(REF)>1` and `len(ALT)=1` | `ATG → A` |
| `Symbolic` | `ALT` is `*`, starts with `<`, or contains `[` or `]` | `A → <DEL>` |
| `Complex` | any other unequal-length non-symbolic case | `ATG → CT` |

The **`*` spanning deletion** (VCF 4.3 §5.5) and all four **breakend** forms (§5.4, detected via
`[` / `]`) are classified **Symbolic**, not SNP/indel — one of the five spec-compliance corrections
recorded on [[parse-vcf-001-evidence]].

### Ti/Tv ratio (`CalculateTiTvRatio`)

Computed over SNP alleles as transitions ÷ transversions:

$$
Ti/Tv = \frac{count(AG, GA, CT, TC)}{count(AC, CA, AT, TA, GC, CG, GT, TG)}
$$

Ti/Tv **iterates all ALT alleles** at multi-allelic sites (not just the first). Expected values:
whole-genome ≈ 2.0–2.1, exome ≈ 3.0.

### Zygosity (`IsHomRef` / `IsHomAlt` / `IsHet`)

Helpers read the `GT` sample field and split on `/` or `|` to distinguish homozygous reference,
homozygous alternate, and heterozygous calls. A **missing allele `.`** in `GT` makes zygosity
**indeterminate** — `IsHet` returns false when any allele is missing (neither het nor hom). Allele
frequency = `alt_count / total_alleles`.

## Edge cases & intentional simplifications

- **Null / empty content** in `Parse` / `ParseFile` → **no records** (early-return guards).
- **Null content** in `ParseWithHeader` → **throws** (`StringReader` requires non-null); **empty**
  content → default `VCFv4.3` header + no records.
- **Missing path** in `ParseFileWithHeader` → default `VCFv4.3` header, empty collections, no records.
- **Data line with fewer than 8 columns** → **skipped** (minimum structural validation).
- **Non-integer `POS`** → line **skipped** (mandatory numeric coordinate). *Malformed lines are
  skipped, not raised* — an API-contract choice, not a format rule.
- **`QUAL "."`** → `null`; **`FILTER "."`** → empty array; **INFO flag without `=`** → `"true"`.
- **`FilterPassing`** returns only records whose FILTER is exactly `PASS`; records with `.` are
  **excluded** (different VCF state — INV-03).
- **Stringly-typed values** — INFO and sample fields are stored as strings **without** enforcing
  header-declared `Number`/`Type`; semantic inconsistencies remain caller-visible.
- **Round-trip is lossy** — `WriteToStream` always emits a `##fileformat=` line and formats `QUAL`
  with two decimals; it **normalizes** output rather than preserving exact original text. Only INFO,
  FORMAT, and FILTER metadata are promoted into header structures; other metadata stays as raw
  key/value strings.
- **Not implemented** — full allele-content validation and compressed **BCF / BGZF-indexed**
  workflows (no in-repo alternative).

The format-facts, INFO/FORMAT field tables, symbolic/breakend grammar, spec-compliance corrections,
and testing-methodology categories for this unit live on [[parse-vcf-001-evidence]].
