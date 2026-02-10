# MIRNA-TARGET-001 Evidence: Target Site Prediction

## Sources

| Source | Key Finding | URL / Reference |
|--------|-------------|-----------------|
| Bartel (2009) | Defines seed site hierarchy: 8mer > 7mer-m8 > 7mer-A1 > 6mer; seed positions 2–8; antiparallel binding | PMID: 19167326, doi:10.1016/j.cell.2009.01.002 |
| Lewis et al. (2005) | "Conserved seed pairing, often flanked by adenosines" — adenosine at position 1 opposite correlates with efficacy | PMID: 15652477, doi:10.1016/j.cell.2004.12.035 |
| Grimson et al. (2007) | Context scores: site-type efficacy weights (8mer=0.31, 7mer-m8=0.161, 7mer-A1=0.099); AU-rich context favors targeting | PMID: 17612493, doi:10.1016/j.molcel.2007.06.017 |
| Agarwal et al. (2015) | TargetScan 7: context++ model; 6mer has minimal but detectable efficacy | PMID: 26267216, doi:10.7554/eLife.05005 |
| TargetScan 8.0 FAQ | Centered sites removed (not reliably functional); 3 canonical + offset 6mer | https://www.targetscan.org/vert_80/docs/help.html |
| Wikipedia: MicroRNA | Seed region = positions 2–7/8; target recognition in 3' UTR; antiparallel Watson-Crick+ G:U wobble | https://en.wikipedia.org/wiki/MicroRNA |
| miRBase | Authoritative miRNA sequences for test data | https://mirbase.org/ |
| Friedman et al. (2009) | Probability of conserved targeting (PCT); most mammalian mRNAs are miRNA targets | PMID: 18955434 |

## Seed Site Type Definitions (Bartel 2009, TargetScan)

Target sites are on the mRNA and pair **antiparallel** to the miRNA. The seed region is positions 2–8 (1-indexed from 5' end of miRNA). The reverse complement of the seed is sought in the mRNA (5'→3').

### Site Types (decreasing efficacy)

| Type | Definition | miRNA Positions | Additional Requirement |
|------|-----------|-----------------|----------------------|
| **8mer** | Match to miRNA positions 2-8 + A opposite position 1 | 2–8 + pos 1 | 'A' on mRNA at position opposite miRNA pos 1 |
| **7mer-m8** | Match to miRNA positions 2-8 | 2–8 | None |
| **7mer-A1** | Match to miRNA positions 2-7 + A opposite position 1 | 2–7 + pos 1 | 'A' on mRNA at position opposite miRNA pos 1 |
| **6mer** | Match to miRNA positions 2-7 | 2–7 | None |
| **Offset 6mer** | Match to miRNA positions 3-8 | 3–8 | None (marginal efficacy) |

### Orientation Details
- miRNA binds 3'→5' to mRNA 5'→3' (antiparallel)
- The "reverse complement" of the miRNA seed is what appears on the mRNA target
- Position 1 of miRNA is at the 3' end of the target site on mRNA

### Scoring Hierarchy (Grimson 2007)
From experimental data in mammalian cells:
- 8mer most effective (largest median fold repression)
- 7mer-m8 intermediate-high
- 7mer-A1 intermediate
- 6mer weak but detectable
- Offset 6mer marginal

## Worked Example: hsa-let-7a-5p

**miRNA**: `UGAGGUAGUAGGUUGUAUAGUU`
**Seed (pos 2-8)**: `GAGGUAG`
**Reverse complement of seed**: `CUACCUC`

For an 8mer site on mRNA: `...CUACCUCA...` (the 'A' at the 3' end corresponds to position 1 of miRNA)
For a 7mer-m8 site: `...CUACCUC...` (match to pos 2-8 only)
For a 7mer-A1 site: `...UACCUCA...` (match to pos 2-7 + 'A' at pos 1 position)
For a 6mer site: `...UACCUC...` (match to pos 2-7 only)

## Alignment Conventions
- miRNA aligns from 5'→3' (left to right)
- Target aligns from 3'→5' (left to right, i.e., reverse of mRNA direction)
- Watson-Crick pairs: A:U, G:C → '|'
- Wobble pairs: G:U → ':'
- Mismatches: ' '

## Edge Cases & Boundary Conditions

1. **Empty/null inputs**: FindTargetSites returns empty collection
2. **mRNA shorter than seed**: No possible matches, return empty
3. **Multiple sites in same mRNA**: All should be found independently
4. **Overlapping sites**: Same seed RC at positions that overlap — each should be reported
5. **DNA input (T instead of U)**: Implementation converts T→U before matching
6. **Case insensitivity**: Normalized to uppercase
7. **minScore threshold**: Sites below threshold are filtered out
8. **Score range**: All scores should be in [0.0, 1.0]
9. **8mer vs 7mer-m8 at same position**: 8mer should be reported (higher rank)

## Known Numerical Properties
- Score for 8mer ≥ Score for 7mer-m8 ≥ Score for 7mer-A1 ≥ Score for 6mer ≥ Score for offset 6mer
- FreeEnergy for well-paired duplex should be negative
- SeedMatchLength: 8 for 8mer, 7 for 7mer types, 6 for 6mer types
