# MIRNA-TARGET-001: Target Site Prediction — Test Specification

## Test Unit
- **ID**: MIRNA-TARGET-001
- **Area**: MiRNA
- **Canonical methods**: `FindTargetSites(mRna, miRna, minScore)`, scoring via `CalculateTargetScore` (internal)
- **Supporting**: `AlignMiRnaToTarget`, `CheckSeedMatch` (internal)

## Canonical Test File
`MiRnaAnalyzer_TargetPrediction_Tests.cs`

## Reference Data
- hsa-let-7a-5p: `UGAGGUAGUAGGUUGUAUAGUU`, seed `GAGGUAG`, seed RC `CUACCUC`
- hsa-miR-21-5p: `UAGCUUAUCAGACUGAUGUUGA`, seed `AGCUUAU`, seed RC `AUAAGCU`

---

## Must Tests (evidence-backed)

### M-001: FindTargetSites — 8mer site detected and scored highest
**Evidence**: Bartel (2009) — 8mer = seed match positions 2-8 + A opposite position 1
- Construct mRNA with exact 8mer site (`CUACCUCA` for let-7a)
- Assert: site found, Type = Seed8mer, SeedMatchLength = 8, score ≥ 0.9

### M-002: FindTargetSites — 7mer-m8 site detected
**Evidence**: Bartel (2009) — 7mer-m8 = match to positions 2-8
- Construct mRNA with seed RC but no trailing A (`CUACCUCG`)
- Assert: site found, Type = Seed7merM8, SeedMatchLength = 7

### M-003: FindTargetSites — 7mer-A1 site detected
**Evidence**: Bartel (2009) — 7mer-A1 = positions 2-7 + A opposite position 1
- Construct mRNA with 6mer seed RC + preceding A (`AUACCUC`)
- Assert: site found, Type = Seed7merA1, SeedMatchLength = 7

### M-004: FindTargetSites — 6mer site detected
**Evidence**: Bartel (2009) — 6mer = match to positions 2-7 only
- Construct mRNA where only 6-mer (pos 2-7 RC) matches, no A context
- Assert: site found, Type = Seed6mer, SeedMatchLength = 6

### M-005: Score monotonicity — 8mer > 7mer-m8 > 7mer-A1 > 6mer
**Evidence**: Grimson (2007) — efficacy hierarchy is well-established
- Create four mRNAs, each with exactly one site type
- Assert: score(8mer) > score(7mer-m8) > score(7mer-A1) > score(6mer)

### M-006: FindTargetSites — empty/null inputs return empty
**Evidence**: Defensive API contract
- Empty mRNA → empty, null miRNA sequence → empty

### M-007: FindTargetSites — no match returns empty
**Evidence**: Trivial correctness
- mRNA with no seed RC occurrence → empty result

### M-008: FindTargetSites — multiple sites found independently
**Evidence**: Bartel (2009) — each seed match site functions independently
- mRNA with 3 copies of seed RC separated by spacers
- Assert: Count ≥ 3 sites

### M-009: Score range is [0.0, 1.0]
**Evidence**: Implementation contract; Grimson (2007) context scores are bounded
- For all site types, assert 0 ≤ score ≤ 1

### M-010: FindTargetSites — minScore filtering works
**Evidence**: API contract
- Use high minScore threshold (0.99); verify low-scoring sites excluded

### M-011: AlignMiRnaToTarget — perfect complementary yields all matches
**Evidence**: Watson-Crick pairing definition
- miRNA = "AAAA", target = "UUUU" → 4 matches, 0 mismatches

### M-012: AlignMiRnaToTarget — G:U wobble detected
**Evidence**: Wobble pair definition (Crick 1966)
- GGGG vs UUUU → 4 wobbles

### M-013: AlignMiRnaToTarget — mismatches counted
**Evidence**: Trivial correctness
- AAAA vs AAAA → 4 mismatches (same bases don't pair)

### M-014: AlignMiRnaToTarget — empty input returns empty duplex
**Evidence**: Defensive API contract

### M-015: Free energy is negative for well-paired duplex
**Evidence**: Thermodynamic principle; Bartel (2009)
- Perfect complement pairing → negative free energy

### M-016: Target site includes alignment string
**Evidence**: API contract
- Found site has non-empty Alignment field

### M-017: DNA input (T) handled equivalently to RNA (U)
**Evidence**: Implementation design — T→U conversion
- mRNA with T bases, same result as U

---

## Should Tests

### S-001: Offset 6mer detected
- mRNA with match to positions 3-8 RC → Type = Offset6mer

### S-002: TargetSite record fields populated correctly
- All fields (Start, End, TargetSequence, MiRnaName, Type, Score, FreeEnergy, Alignment) populated

### S-003: Real miRNA (let-7a) against constructed target
- Integration-style test with biologically plausible sequences

---

## Could Tests

### C-001: Context (AU-rich) may influence scoring indirectly via supplementary pairing bonuses

---

## Audit Notes

### Existing tests in MiRnaAnalyzerTests.cs — Target Site section
| Test | Classification | Decision |
|------|---------------|----------|
| FindTargetSites_PerfectSeedMatch_FindsSite | Weak (vague assertion) | Replace with M-001 |
| FindTargetSites_8mer_HighestScore | Weak (conditional) | Replace with M-005 |
| FindTargetSites_NoMatch_ReturnsEmpty | Covered | Migrate as M-007 |
| FindTargetSites_MultipleSites_FindsAll | Covered | Migrate as M-008 |
| FindTargetSites_EmptySequences_ReturnsEmpty | Covered | Migrate as M-006 |
| FindTargetSites_IncludesAlignment | Weak (conditional) | Replace with M-016 |
| AlignMiRnaToTarget_PerfectMatch_AllMatches | Covered | Migrate as M-011 |
| AlignMiRnaToTarget_AllMismatches_NoMatches | Covered | Migrate as M-013 |
| AlignMiRnaToTarget_WobblePairs_Detected | Covered | Migrate as M-012 |
| AlignMiRnaToTarget_EmptySequences_ReturnsEmptyDuplex | Covered | Migrate as M-014 |
| AlignMiRnaToTarget_CalculatesFreeEnergy | Weak (non-zero) | Replace with M-015 |
| FullWorkflow_PredictTargets | Integration | Replace with S-003 |
| TargetSite_HasAllFields | Covered | Replace with S-002 |
