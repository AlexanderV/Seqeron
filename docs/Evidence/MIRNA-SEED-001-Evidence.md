# MIRNA-SEED-001 Evidence: Seed Sequence Analysis

## Sources

| Source | Key Finding | URL |
|--------|-------------|-----|
| Wikipedia: MicroRNA | Seed region = positions 2–7 of miRNA; must be perfectly complementary for target recognition | https://en.wikipedia.org/wiki/MicroRNA |
| TargetScan FAQ | Defines canonical seed site types: 8mer, 7mer-m8, 7mer-A1, 6mer based on positions 2-8 | https://www.targetscan.org/faqs.Release_7.html |
| TargetScan Docs (7mer) | 7mer-m8: exact match to positions 2-8; 7mer-A1: positions 2-7 + 'A' | https://www.targetscan.org/docs/7mer.html |
| Lewis et al. (2005) | "Conserved seed pairing, often flanked by adenosines" — defines seed-based target prediction | PMID: 15652477 |
| Bartel (2009) | "MicroRNAs: Target recognition and regulatory functions" — defines seed region and site hierarchy | PMID: 19167326, doi:10.1016/j.cell.2009.01.002 |
| Grimson et al. (2007) | "MicroRNA targeting specificity in mammals: determinants beyond seed pairing" — context scores | PMID: 17612493 |
| Agarwal et al. (2015) | "Predicting effective microRNA target sites in mammalian mRNAs" — TargetScan 7, context++ model | PMID: 26267216, doi:10.7554/eLife.05005 |
| miRBase | Authoritative database of miRNA sequences and nomenclature | https://mirbase.org/ |
| Friedman et al. (2009) | Mammalian miRNAs conserved targets — PCT scoring | PMID: 18955434 |

## Key Definitions (from sources)

### Seed Region
- **Canonical seed**: positions 2–7 of the mature miRNA (6 nucleotides)
- **Extended seed**: positions 2–8 (7 nucleotides, includes m8 position)
- Positions are 1-indexed from 5' end of the mature miRNA

### miRNA Seed Site Types (TargetScan, Bartel 2009)
Hierarchical effectiveness (highest to lowest):
1. **8mer**: match to positions 2-8 (seed + pos 8) + 'A' opposite position 1
2. **7mer-m8**: match to positions 2-8 (seed + position 8)
3. **7mer-A1**: match to positions 2-7 (seed) + 'A' opposite position 1
4. **6mer**: match to positions 2-7 (seed only)

### miRNA Family Concept
- miRNAs sharing the same seed sequence (positions 2-7 or 2-8) belong to the same family (Griffiths-Jones et al., 2006; TargetScan)
- Family members are expected to target overlapping sets of mRNAs

## Reference miRNA Sequences (from miRBase)

| miRNA | Accession | Mature Sequence | Seed (pos 2-8) |
|-------|-----------|-----------------|-----------------|
| hsa-let-7a-5p | MIMAT0000062 | UGAGGUAGUAGGUUGUAUAGUU | GAGGUAG |
| hsa-let-7b-5p | MIMAT0000063 | UGAGGUAGUAGGUUGUGUGGUU | GAGGUAG |
| hsa-let-7c-5p | MIMAT0000064 | UGAGGUAGUAGGUUGUAUGGUU | GAGGUAG |
| hsa-miR-21-5p | MIMAT0000076 | UAGCUUAUCAGACUGAUGUUGA | AGCUUAU |
| hsa-miR-155-5p | MIMAT0000646 | UUAAUGCUAAUUGUGAUAGGGGU | UAAUGCU |
| hsa-miR-1-3p | MIMAT0000416 | UGGAAUGUAAAGAAGUAUGUAU | GGAAUGU |

## Edge Cases (documented in sources)

1. **miRNA shorter than 8 nt**: Cannot extract seed — return empty (ASSUMPTION: graceful handling)
2. **Null/empty input**: Return empty string (defensive programming)
3. **DNA vs RNA input**: T should be handled equivalently to U in CreateMiRna (converts T→U)
4. **Case insensitivity**: Seed extraction normalizes to uppercase
5. **let-7 family**: All members share identical seed — validates family grouping
6. **Same miRNA**: CompareSeedRegions with self should yield 0 mismatches
7. **Completely different seeds**: Maximum mismatch count = seed length

## Known Constraints
- miRNA sequences are typically 18-25 nt in length
- Seed region is the primary determinant of target specificity in animals
- Seed comparison is a string-level operation (Hamming distance over fixed-length region)
