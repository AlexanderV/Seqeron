# Test Specification: PROBE-VALID-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PROBE-VALID-001 |
| **Area** | MolTools |
| **Title** | Probe Validation |
| **Canonical Class** | `ProbeDesigner` |
| **Canonical Methods** | `ValidateProbe`, `CheckSpecificity`, `ScanOffTargetsGapped`, `ComputeLambdaNucleotide`, `ComputeKarlinAltschul` |
| **Complexity** | O(n × g) ungapped; O(g × n·m) gapped scan |
| **Status** | ☐ Not Started (re-validation pending after limitation fix) |
| **Last Updated** | 2026-06-24 |

---

## Evidence Sources

| Source | Type | Key Information |
|--------|------|-----------------|
| Wikipedia: Hybridization probe | Academic | Probe hybridization, stringency-dependent cross-hybridization |
| Wikipedia: DNA microarray | Academic | Probe specificity, cross-hybridization in high-density arrays |
| Wikipedia: Off-target genome editing | Academic | CRISPR/Cas9 tolerates 3-5 bp mismatches per 20nt guide |
| Wikipedia: Nucleic acid thermodynamics | Academic | Mismatch destabilization, complementarity metrics |
| Wikipedia: Off-target activity | Academic | Off-target detection methods, mismatch tolerance mechanisms |
| Smith & Waterman (1981), J Mol Biol 147:195 | Primary paper | Local-alignment recurrence with zero floor; indel-aware (reused via `SequenceAligner.LocalAlign`) |
| Altschul et al. (1990), J Mol Biol 215:403 | Primary paper | Gapped local alignment finds indels the ungapped scan misses ("BLAST-grade") |
| Kane et al. (2000), Nucleic Acids Res 28(22):4552 | Primary paper | >75% identity over the probe → cross-hybridization (off-target identity threshold) |
| Karlin & Altschul (1990), PNAS 87:2264 | Primary paper | λ = unique positive root of Σ p_i p_j e^{λ s_ij} = 1; E = K·m·n·e^{−λS}; negative-expected-score precondition |
| Altschul et al. (1990), J Mol Biol 215:403 | Primary paper | Bit score S' = (λS − ln K)/ln 2; E = m·n·2^{−S'}; λ ≈ 1.37 / K ≈ 0.711 for +1/−3 (NCBI blastn) |

---

## Invariants

1. **Specificity Range**: 0.0 ≤ specificityScore ≤ 1.0 (Source: Implementation)
2. **Self-Complementarity Range**: 0.0 ≤ selfComplementarity ≤ 1.0 (Source: Mathematical)
3. **Off-Target Non-Negative**: offTargetHits ≥ 0 (Source: Implementation)
4. **Unique Match Specificity**: offTargetHits == 1 → specificityScore == 1.0 (Source: Implementation)
5. **No Match Zero Specificity**: offTargetHits == 0 → specificityScore == 0.0 (Source: Implementation)
6. **Multiple Hits**: offTargetHits > 1 → specificityScore == 1.0 / offTargetHits (Source: Implementation)
7. **Gapped Identity Range**: each gapped hit has 0.0 ≤ Identity ≤ 1.0 and 0.0 ≤ Coverage ≤ 1.0 (Source: Mathematical — identical/ungapped columns ÷ probe length)
8. **On/Off Separation**: the first perfect ungapped full-coverage exact match is on-target and is excluded from `OffTargetHits`; all imperfect/indel hits ≥ minIdentity are off-target (Source: Kane et al. 2000 intended-vs-non-intended signal; on/off labelling assumption)
9. **Indel Detection**: a hit reachable only via an indel has `HasGaps == true` and is found by `ScanOffTargetsGapped` but not by the ungapped `ValidateProbe` Hamming scan (Source: Altschul et al. 1990; Smith & Waterman 1981)
10. **Identity Threshold**: a hit is reported iff Identity ≥ `minIdentity` (default 0.75) (Source: Kane et al. 2000)
11. **Karlin–Altschul λ**: `ComputeLambdaNucleotide` returns the unique positive root of Σ p_i p_j e^{λ s_ij} = 1; for +1/−3, p=0.25 it equals 1.3740631 ≈ published 1.374 (Source: Karlin & Altschul 1990; NCBI blastn cross-check). Requires a positive score and negative expected score.
12. **Karlin–Altschul E-value/bit-score**: E = K·m·n·e^{−λS} = m·n·2^{−S'} with S' = (λS − ln K)/ln 2; E strictly decreases in S and is linear in m·n (Source: Karlin & Altschul 1990; Altschul et al. 1990)

---

## Test Cases

### Must (Required - Evidence-Based)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| M1 | Empty probe returns validation with zero specificity | Boundary condition | Implementation spec |
| M2 | Null/empty references returns validation without off-target hits | Boundary condition | Implementation spec |
| M3 | Unique probe (1 hit) has specificity score 1.0 | Invariant #4 | Implementation |
| M4 | Multiple hits reduce specificity to 1.0/hitCount | Invariant #6 | Implementation |
| M5 | Specificity score always in [0.0, 1.0] | Invariant #1 | Implementation |
| M6 | Self-complementarity in [0.0, 1.0] | Invariant #2 | Mathematical |
| M7 | OffTargetHits is non-negative | Invariant #3 | Implementation |
| M8 | High self-complementarity (>30%) reported in issues | Quality criterion | Implementation spec |
| M9 | CheckSpecificity with suffix tree returns correct score | API contract | Implementation |
| M10 | CheckSpecificity unique match returns 1.0 | Invariant #4 | Implementation |
| M11 | CheckSpecificity multiple matches returns 1.0/count | Invariant #6 | Implementation |
| M12 | Case-insensitive probe handling | Usability | Implementation |
| MG1 | Indel off-target found by gapped scan, missed by ungapped Hamming scan | Invariant #9 | Altschul et al. 1990; hand-derived |
| MG2 | On-target exact match reported separately, excluded from off-target count | Invariant #8 | Kane et al. 2000; on/off separation |
| MG3 | Indel off-target identity/coverage = 1.0, HasGaps, exact aligned strings | Invariant #7 | Smith & Waterman 1981 (hand-derived) |
| MG4 | Indel+mismatch off-target identity = 10/12 = 0.8333 (SW trims mismatched tail) | Invariant #7 | Smith & Waterman 1981 (hand-derived) |
| KA1 | λ for +1/−3, uniform 0.25 = 1.3740631 (≈ published 1.374) within 1e-6 | Invariant #11 | Karlin & Altschul 1990; NCBI blastn |
| KA2 | The solved λ satisfies its defining equation 0.25·e^λ + 0.75·e^{−3λ} = 1 | Invariant #11 | Karlin & Altschul 1990 |
| KA3 | Hand-derived (S=30,m=20,n=1000,K=0.711): S'=59.9627, E=1.7802e−14 | Invariant #12 | Altschul et al. 1990 (hand-derived) |
| KA4 | The two E-value forms agree: K·m·n·e^{−λS} == m·n·2^{−S'} | Invariant #12 | Altschul et al. 1990 |
| KA5 | E strictly decreases as raw score S increases | Invariant #12 | Karlin & Altschul 1990 |
| KA6 | E scales linearly with the search space m·n (double n → double E) | Invariant #12 | Karlin & Altschul 1990 |
| KA7 | λ guards: non-positive match throws; non-negative expected score throws; non-positive m/n/K throws | Invariant #11 | Karlin & Altschul 1990 (preconditions) |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | Secondary structure potential detected for hairpin sequences | Quality metric | Implementation |
| S2 | Issues list populated for problematic probes | User feedback | Implementation |
| S3 | IsValid false when multiple issues exist | Validation logic | Implementation |
| S4 | Approximate matching with maxMismatches works correctly | Off-target detection | Wikipedia (Off-target) |
| SG1 | minIdentity threshold gates hits (0.8333 admitted at 0.75, rejected at 0.90) | Invariant #10 | Kane et al. 2000 |
| SG2 | Specific probe (1 on-target, 0 off-target) → IsSpecific true | API contract | Implementation |

### Could (Optional)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| C1 | Long reference sequences handled efficiently | Performance | Implementation |
| C2 | Multiple references all searched | Completeness | Implementation |

---

## Coverage Classification

| ID | Test Method | Classification | Notes |
|----|------------|----------------|-------|
| M1 | `ValidateProbe_EmptyProbe_ReturnsValidationResult` | ✅ Covered | Exact: specificity=0.0, offTargetHits=0, IsValid=false |
| M2 | `ValidateProbe_EmptyReferences_ReturnsValidationWithNoOffTargetHits` | ✅ Covered | Exact: specificity=0.0, offTargetHits=0 (Inv. #5) |
| — | `ValidateProbe_NullReferences_ThrowsArgumentNullException` | ✅ Covered | ArgumentNullException |
| M3 | `ValidateProbe_UniqueProbe_HasSpecificityScoreOne` | ✅ Covered | Exact: offTargetHits=1, specificity=1.0 |
| M4 | `ValidateProbe_MultipleHits_ReducesSpecificityByHitCount` | ✅ Covered | Exact: offTargetHits=25, specificity=0.04 |
| M5 | ~~`ValidateProbe_AnyInput_SpecificityScoreInValidRange`~~ | 🔁 Deleted | Subsumed by AllInvariants + every exact-value test |
| M6 | ~~`ValidateProbe_AnyInput_SelfComplementarityInValidRange`~~ | 🔁 Deleted | Subsumed by AllInvariants + M8 |
| M7 | ~~`ValidateProbe_AnyInput_OffTargetHitsNonNegative`~~ | 🔁 Deleted | Subsumed by AllInvariants + all specific tests |
| M8 | `ValidateProbe_HighSelfComplementarity_ReportsInIssues` | ✅ Covered | Exact: selfComp=1.0, issues contain "Self-complementarity" |
| M9 | ~~`CheckSpecificity_ResultInValidRange`~~ | 🔁 Deleted | Subsumed by M10 + M11 + NoMatch exact tests |
| M9 | `CheckSpecificity_NoMatch_ReturnsZero` | ✅ Covered | Exact: specificity=0.0 |
| M10 | `CheckSpecificity_UniqueSequence_ReturnsOne` | ✅ Covered | Exact: specificity=1.0 |
| M11 | `CheckSpecificity_MultipleOccurrences_ReturnsOneOverCount` | ✅ Covered | Exact: 1.0/count |
| M12 | `ValidateProbe_MixedCaseProbe_HandledCaseInsensitively` | ✅ Covered | Upper/lower/mixed → same results |
| S1 | `ValidateProbe_PotentialHairpin_DetectsSecondaryStructure` | ✅ Covered | Exact: HasSecondaryStructure=true, issues contain text |
| S2 | `ValidateProbe_ProblematicProbe_PopulatesIssuesList` | ✅ Covered | Exact: offTargetHits=16, issues contain "16 potential off-target sites" |
| S3 | `ValidateProbe_MultipleProblems_IsValidFalse` | ✅ Covered | Exact: hits=12, selfComp=1.0, issues=2, IsValid=false |
| S4 | `ValidateProbe_ApproximateMatching_FindsNearMatches` | ✅ Covered | Exact: 0 hits (strict) vs 1 hit (approx) |
| C1 | `ValidateProbe_LongReference_FindsProbeCorrectly` | ✅ Covered | 20k-char ref, exact: hits=1, specificity=1.0 |
| C2 | `ValidateProbe_MultipleReferences_AccumulatesHits` | ✅ Covered | 3 refs, exact: hits=3, specificity=1/3 |
| — | `ValidateProbe_AllInvariants_HoldForTypicalProbe` | ✅ Covered | All 6 invariants verified for typical input |
| — | `ValidateProbe_ZeroHits_ReturnsZeroSpecificity` | ✅ Covered | Explicit Invariant #5 |
| — | `DesignProbes_WithGenomeIndex_UsesCheckSpecificity` | ✅ Covered | Integration: suffix tree + DesignProbes |
| MG1 | `ScanOffTargetsGapped_IndelOffTarget_FoundByGappedScanMissedByHammingScan` | ✅ Covered | Hamming offTargetHits=1; gapped on=1, off=1, HasGaps |
| MG2 | `ScanOffTargetsGapped_OnTargetExactMatch_NotCountedAsOffTarget` | ✅ Covered | On-target start=5,end=16,identity=1.0,cov=1.0,no gap; not in off-targets |
| MG3 | `ScanOffTargetsGapped_IndelOffTarget_HasExactHandDerivedIdentity` | ✅ Covered | start=27,identity=1.0,cov=1.0,HasGaps; aligned strings exact |
| MG4 | `ScanOffTargetsGapped_IndelPlusMismatch_IdentityIsHandDerivedFraction` | ✅ Covered | identity=10/12=0.8333, HasGaps, off=1 |
| SG1 | `ScanOffTargetsGapped_IdentityThreshold_GatesHits` | ✅ Covered | 0.75→1 hit, 0.90→0 hits |
| SG2 | `ScanOffTargetsGapped_SpecificProbe_IsSpecificTrue` | ✅ Covered | on=1, off=0, IsSpecific |
| — | `ScanOffTargetsGapped_NullProbe_ThrowsArgumentNullException` | ✅ Covered | ArgumentNullException |
| — | `ScanOffTargetsGapped_NullReferences_ThrowsArgumentNullException` | ✅ Covered | ArgumentNullException |
| — | `ScanOffTargetsGapped_EmptyProbe_ReturnsNoHits` | ✅ Covered | Empty probe → no hits |
| KA1 | `ComputeLambdaNucleotide_Plus1Minus3_UniformFrequencies_MatchesPublishedValue` | ✅ Covered | λ = 1.3740631 within 1e-6 |
| KA2 | `ComputeLambdaNucleotide_SolvedRoot_SatisfiesDefiningEquation` | ✅ Covered | 0.25·e^λ+0.75·e^{−3λ}=1 within 1e-9 |
| KA3 | `ComputeKarlinAltschul_HandDerivedExample_MatchesBitScoreAndEValue` | ✅ Covered | S'=59.9627, E=1.7802e−14 |
| KA4 | `ComputeKarlinAltschul_EValue_EqualsSearchSpaceTimesTwoToMinusBitScore` | ✅ Covered | K·m·n·e^{−λS}==m·n·2^{−S'} |
| KA5 | `ComputeKarlinAltschul_EValue_DecreasesAsScoreIncreases` | ✅ Covered | E(S+1) < E(S) |
| KA6 | `ComputeKarlinAltschul_EValue_ScalesLinearlyWithSearchSpace` | ✅ Covered | E(2n) = 2·E(n) |
| KA7 | `ComputeLambdaNucleotide_NonPositiveMatch_Throws` / `..._NonNegativeExpectedScore_Throws` / `ComputeKarlinAltschul_NonPositiveLength_Throws` | ✅ Covered | Preconditions + arg guards throw |

---

## Evidence-Backed Parameters

All configurable parameters have external evidence justification. No assumptions remain.

| Parameter | Default | Evidence | Source |
|-----------|---------|----------|--------|
| `maxMismatches` | 3 | Lower bound of CRISPR/Cas9 off-target tolerance: 3-5 bp mismatches per 20nt guide (Hsu et al. 2013, Fu et al. 2013). Configurable per caller. | Wikipedia: Off-target genome editing |
| `selfComplementarityThreshold` | 0.3 | For random DNA (uniform base distribution), expected self-complementarity ≈ 0.25. Threshold of 0.3 (20% above baseline) detects statistically elevated palindromic character. Per-application defaults: qPCR=0.25, Microarray=0.3, NorthernBlot=0.35, FISH/SouthernBlot=0.4. | Wikipedia: Nucleic acid thermodynamics; `ProbeParameters.Defaults` |
| `minIdentity` (gapped scan) | 0.75 | Kane et al. (2000): non-target transcripts >75% similar over the probe may cross-hybridize. Caller-configurable. | Kane et al. (2000), Nucleic Acids Res 28(22):4552 |
| `scoring` (gapped scan) | `SequenceAligner.BlastDna` (+2/−3, gap −2) | Reuses the BLAST-style DNA scoring already used for the library's gapped ANI alignment (COMPGEN-ANI). | Altschul et al. (1990); reused infrastructure |
| `k` (Karlin–Altschul) | 0.711 | Published nucleotide K for the +1/−3 scheme (NCBI blastn). K's full closed form needs the Karlin–Altschul score-lattice machinery, so it is a caller parameter; λ is computed (not assumed). | Karlin & Altschul (1990); NCBI blastn |
| `baseFrequency` (λ) | 0.25 | Uniform four-base background — the standard assumption used to derive the +1/−3 λ ≈ 1.374. | Karlin & Altschul (1990) |


