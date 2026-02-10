# Target Site Prediction

## Algorithm Definition

miRNA target site prediction identifies regions on an mRNA (typically in the 3' UTR) that are complementary to a miRNA's seed region and thus are candidate binding sites for post-transcriptional regulation.

**Sources**: Bartel (2009) "MicroRNAs: Target recognition and regulatory functions" (Cell 136:215–233); Lewis et al. (2005) "Conserved seed pairing, often flanked by adenosines" (Cell 120:15–20); Grimson et al. (2007) "MicroRNA targeting specificity in mammals: determinants beyond seed pairing" (Mol Cell 27:91–105).

### Core Mechanism

The miRNA seed region (positions 2–8 from the 5' end) pairs antiparallel to the mRNA target. The reverse complement of the seed sequence is sought in the mRNA. Different degrees of complementarity define site types with a well-established efficacy hierarchy (Bartel 2009):

| Site Type | Positions Matched | Additional Req. | Relative Efficacy |
|-----------|-------------------|-----------------|-------------------|
| 8mer | 2–8 + A@pos1 | A opposite miRNA pos 1 | Highest |
| 7mer-m8 | 2–8 | None | High |
| 7mer-A1 | 2–7 + A@pos1 | A opposite miRNA pos 1 | Intermediate |
| 6mer | 2–7 | None | Low |
| Offset 6mer | 3–8 | None | Marginal |

## Formal Properties and Invariants

1. **Score monotonicity by site type**: score(8mer) ≥ score(7mer-m8) ≥ score(7mer-A1) ≥ score(6mer) ≥ score(offset 6mer), all else being equal.
2. **Score range**: All scores ∈ [0.0, 1.0].
3. **Antiparallel binding**: The seed reverse complement is the literal string that appears on the mRNA target.
4. **Completeness**: Every occurrence of the seed RC (or subset) in the mRNA should yield a candidate site (subject to minScore filtering).
5. **Site type priority**: When multiple site types match at the same position, the highest-efficacy type is reported (8mer > 7mer-m8 > 7mer-A1 > 6mer > offset 6mer).
6. **Free energy**: Well-paired duplexes yield negative free energy values.
7. **Alignment string**: Character-level alignment (|, :, space) for each paired position.

## Input Domain

- **mRNA sequence**: Non-empty string of {A, U, G, C, T} characters (T converted to U). Typically 3' UTR.
- **miRNA**: `MiRna` record with valid Name, Sequence (≥8 nt), and SeedSequence.
- **minScore**: Double in [0.0, 1.0], default 0.5. Threshold for reporting.

## Failure Modes

- Empty/null mRNA or miRNA → empty result (no exception)
- miRNA shorter than 8 nt → no meaningful seed extraction; may yield no matches
- mRNA shorter than 6 → no possible seed match

## Implementation Notes

- `FindTargetSites`: Scans mRNA linearly; at each position checks for seed RC matches in order of decreasing site type priority.
- `CheckSeedMatch`: Internal method implementing the site type classification.
- `CalculateTargetScore`: Assigns base score by site type, adjusts for 3' supplementary pairing and mismatches. Clamped to [0, 1].
- `AlignMiRnaToTarget`: Full-length miRNA-target alignment; counts matches, mismatches, G:U wobbles.

## Deviations / Notes

- TargetScan 8.0 removed centered sites as unreliable (McGeary et al. 2019). The implementation includes `TargetSiteType.Centered` as an enum value but `CheckSeedMatch` does not generate centered sites. This is consistent with current best practice.
- The energy calculation is simplified (not full nearest-neighbor thermodynamics). Scores are relative, not absolute kcal/mol values.
- Offset 6mer is included despite marginal efficacy, consistent with TargetScan 4–7 inclusion.
