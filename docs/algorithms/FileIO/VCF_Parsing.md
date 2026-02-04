# VCF Parsing

## Overview

VCF (Variant Call Format) is a standardized text file format for storing genetic variation data. It was developed in 2010 for the 1000 Genomes Project and has become the de facto standard for variant data exchange in bioinformatics.

**Specification Version:** VCFv4.3  
**Complexity:** O(n) where n is the file size  
**Reference:** Danecek et al. (2011), Bioinformatics 27(15):2156-2158

---

## File Structure

### Header Section

Lines starting with `##` contain metadata:

```
##fileformat=VCFv4.3
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FILTER=<ID=q10,Description="Quality below 10">
#CHROM  POS  ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  Sample1
```

### Data Section

Tab-delimited records with 8 mandatory columns:

| Column | Name | Type | Description |
|--------|------|------|-------------|
| 1 | CHROM | String | Chromosome identifier |
| 2 | POS | Integer | 1-based position |
| 3 | ID | String | Variant ID or "." |
| 4 | REF | String | Reference allele |
| 5 | ALT | String | Alternate allele(s) |
| 6 | QUAL | Float/. | Quality score |
| 7 | FILTER | String | Filter status |
| 8 | INFO | String | Additional info |
| 9 | FORMAT | String | Sample format (optional) |
| 10+ | SAMPLE | String | Sample data (optional) |

---

## Variant Classification

### Type Determination Algorithm

```
function ClassifyVariant(ref, alt):
    if alt starts with '<' or '[' or ']':
        return Symbolic
    
    refLen = length(ref)
    altLen = length(alt)
    
    if refLen == altLen:
        if refLen == 1:
            return SNP
        else:
            return MNP
    
    if refLen == 1 and altLen > 1:
        return Insertion
    
    if refLen > 1 and altLen == 1:
        return Deletion
    
    return Complex
```

### Variant Types

| Type | Condition | Example |
|------|-----------|---------|
| SNP | ref=1, alt=1, ref≠alt | A→G |
| MNP | ref=n, alt=n, n>1 | AT→GC |
| Insertion | ref<alt | A→ATG |
| Deletion | ref>alt | ATG→A |
| Symbolic | alt starts with < | A→<DEL> |
| Complex | Other | ATG→CT |

---

## Genotype Analysis

### Genotype Notation

| Notation | Meaning |
|----------|---------|
| 0/0 | Homozygous reference |
| 0/1 | Heterozygous |
| 1/1 | Homozygous alternate |
| 1/2 | Heterozygous (both alternate) |
| ./. | Missing genotype |
| 0\|1 | Phased heterozygous |

### Zygosity Classification

```
function ClassifyZygosity(genotype):
    alleles = split(genotype, '/|')
    
    if all alleles == '0':
        return HomozygousReference
    
    if alleles[0] == alleles[1] and alleles[0] != '0':
        return HomozygousAlternate
    
    if alleles[0] != alleles[1]:
        return Heterozygous
```

---

## Transition/Transversion Ratio

### Definition

**Transitions:** Purine↔Purine (A↔G) or Pyrimidine↔Pyrimidine (C↔T)  
**Transversions:** Purine↔Pyrimidine (A↔C, A↔T, G↔C, G↔T)

### Calculation

```
transitions = {AG, GA, CT, TC}
transversions = {AC, CA, AT, TA, GC, CG, GT, TG}

Ti/Tv ratio = count(transitions) / count(transversions)
```

### Expected Values (Source: Danecek et al. 2011)
- Whole genome: ~2.0-2.1
- Exome/coding: ~3.0

---

## Implementation Notes

### Seqeron.Genomics.IO.VcfParser

The implementation provides:

**Parsing Methods:**
- `Parse(content)` - Parse VCF content string
- `ParseFile(filePath)` - Parse VCF file
- `ParseWithHeader(content)` - Parse with header extraction

**Variant Classification:**
- `ClassifyVariant(record)` - Determine variant type
- `IsSNP(record)` - Check if SNP
- `IsIndel(record)` - Check if insertion/deletion
- `GetVariantLength(record)` - Calculate indel length

**Filtering:**
- `FilterByChrom(records, chrom)` - Filter by chromosome
- `FilterByRegion(records, chrom, start, end)` - Filter by position
- `FilterByQuality(records, minQuality)` - Filter by quality score
- `FilterPassing(records)` - Filter PASS only
- `FilterSNPs(records)` - Filter SNPs only
- `FilterIndels(records)` - Filter indels only

**Genotype Analysis:**
- `GetGenotype(record, sampleIndex)` - Get sample genotype
- `IsHomRef(record, sampleIndex)` - Check homozygous reference
- `IsHomAlt(record, sampleIndex)` - Check homozygous alternate
- `IsHet(record, sampleIndex)` - Check heterozygous

**Statistics:**
- `CalculateStatistics(records)` - Compute summary statistics
- `CalculateTiTvRatio(records)` - Calculate Ti/Tv ratio

---

## Edge Cases

| Case | Behavior |
|------|----------|
| Empty content | Return empty enumerable |
| Null content | Return empty enumerable |
| Missing QUAL (.) | Set Qual to null |
| Missing ID (.) | Store "." as ID |
| Missing FILTER (.) | Empty filter array |
| Malformed line | Skip (return null) |
| Non-integer POS | Skip record |
| < 8 columns | Skip record |

---

## References

1. Danecek P, et al. "The variant call format and VCFtools." Bioinformatics. 2011;27(15):2156-2158.
2. SAMtools HTS-specs VCF Specification v4.3. https://samtools.github.io/hts-specs/
3. Wikipedia. "Variant Call Format." https://en.wikipedia.org/wiki/Variant_Call_Format
