# Test Specification: PROBE-VALID-001

## Test Unit Information

| Field | Value |
|-------|-------|
| **Test Unit ID** | PROBE-VALID-001 |
| **Area** | MolTools |
| **Title** | Probe Validation |
| **Canonical Class** | `ProbeDesigner` |
| **Canonical Methods** | `ValidateProbe`, `CheckSpecificity`, `AssessCrossHybridization`, `ScanOffTargetsGapped`, `ComputeLambdaNucleotide`, `ComputeKarlinAltschul` |
| **Complexity** | O(n × g) ungapped; O(g × n·m) gapped scan |
| **Status** | ☑ Reviewed (B07 campaign 2026-09, PROBE-VALID-001) |
| **Last Updated** | 2026-10-08 |

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
| Kane et al. (2000), Nucleic Acids Res 28(22):4552 | Primary paper | >75% identity over the probe or a >15-nt contiguous identical stretch → cross-hybridization |
| Satya et al. (2008), BMC Bioinformatics 9:185 | Secondary | Kane rule as "identity > 75% or contiguous match > 15 bp" |
| Primer3 `libprimer3.cc` / primer3-py 2.3.1 | Reference implementation | Internal-oligo ntthal self-dimer / 3′ self-dimer / hairpin Tm ≤ 47 °C; `calc_homodimer` / `calc_end_stability` / `calc_hairpin` / `calc_heterodimer` |
| Rouillard et al. (2003), NAR 31:3057 (OligoArray 2.0) | Primary paper | Off-target duplex Tm vs a user specificity threshold |
| Biopython 1.88 `PairwiseAligner` (local) | Reference implementation | Local-alignment score / identity cross-check |
| Karlin & Altschul (1990), PNAS 87:2264 | Primary paper | λ = unique positive root of Σ p_i p_j e^{λ s_ij} = 1; E = K·m·n·e^{−λS}; negative-expected-score precondition |
| Altschul et al. (1990), J Mol Biol 215:403 | Primary paper | Bit score S' = (λS − ln K)/ln 2; E = m·n·2^{−S'}; λ ≈ 1.37 / K ≈ 0.711 for +1/−3 (NCBI blastn) |

---

## Invariants

1. **Specificity Range**: 0.0 ≤ specificityScore ≤ 1.0 (Source: library-defined uniqueness score 1/N — not a published metric)
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
13. **Thermodynamic self-structure screen**: for ACGT probes of ≤ `ThermodynamicScreenMaxLength` nt (THAL_MAX_ALIGN, default 60; larger = opt-in, thal.c compiled with a larger THAL_MAX_ALIGN parity) (Thermodynamic screen) SelfDimerTm / SelfEndDimerTm / HairpinTm equal primer3-py calc_homodimer / calc_end_stability / calc_hairpin at the stated conditions; a self-complementarity issue iff max(self-dimer, 3′ self-dimer) Tm > MaxStructureTm (47 °C), HasSecondaryStructure iff hairpin Tm > MaxStructureTm; otherwise the fallback screens (Source: Primer3 internal-oligo screen)
14. **IsValid**: IsValid ⇔ Issues is empty (Source: Primer3 rejects an oligo violating any limit; Kane decision); SpecificityScore never enters it (library convention)
17. **Reference sites (Kane)**: CrossHybridizingHits = #{hits with (L − d)/L > maxNonTargetIdentity ∨ longest identical run > maxContiguousMatch}; off-target issue ⇔ CrossHybridizingHits > 1; for L ≥ 13 and maxMismatches ≤ 3, CrossHybridizingHits = OffTargetHits (Source: Kane et al. 2000; audit round 3, A3-12)
18. **Both strands (opt-in `bothStrands`)**: sites = ∪ over references of the hit positions of the probe and of its reverse complement; a position hit in both orientations is one site, a Kane site when either orientation meets a criterion; a reverse-palindromic probe or site is counted once (= `CheckSpecificity(bothStrands)` for exact matching); default false = given strand only, bit-identical to before; OffTargetHits(both) ≥ OffTargetHits(single), CrossHybridizingHits likewise (Source: blastn 2.12.0+ `-strand both` default; audit round 4, A4-2)
19. **Fallback hairpin screen** (> `ThermodynamicScreenMaxLength` nt, non-ACGT, `Heuristic`): HasSecondaryStructure = `PrimerDesigner.HasHairpinPotential(probe)` — an exactly complementary A·T/G·C stem ≥ 4 bp (library convention, F49) closing a loop ≥ 3 nt (thal.c `min_hrpn_loop` = 3); no mismatch tolerance (audit round 5, A5-2). The flag is reported (HasSecondaryStructure + the `Warnings` entry "Potential secondary structure formation") but is **not an issue** and never makes IsValid false (library convention; audit round 5, A5-5); the ntthal hairpin Tm > MaxStructureTm (PRIMER_INTERNAL_MAX_HAIRPIN_TH) and the fallback self_any / self_end > 12.00 (PRIMER_INTERNAL_MAX_SELF_ANY / _END) remain issues
15. **Kane criteria**: CrossHybridizes ⇔ Identity > maxIdentity (0.75) ∨ LongestContiguousMatch > maxContiguousMatch (15) ∨ (maxDuplexTm given ∧ DuplexTm > maxDuplexTm); AlignmentScore = Biopython local score; LongestContiguousMatch = LCS length; both strands by default (Source: Kane et al. 2000; OligoArray 2.0)
16. **Site duplex Tm**: DuplexTm = primer3-py calc_heterodimer(probe, revcomp(site)).tm (0 if no duplex); null when the probe and the site are both longer than the conditions' `ThermodynamicScreenMaxLength` (THAL_MAX_ALIGN, default 60; thal.c `thal_check_errors` needs only one strand ≤ it — primer3-py raises otherwise), for non-ACGT probes/sites or no site (Source: thal.c, OligoArray 2.0; A3-27)

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
| TH1 | ntthal Tm of GCGC…(20) = 78.85652531616256 / 78.85652531616256 / 87.30265612393043; issues + IsValid false | Invariant #13 | primer3-py 2.3.1 |
| TH2 | Stated conditions (mv 100, dv 2, dntp 0.2, dna 250) → 69.17069845823409 / 69.17069845823409 / 74.99462150250321 | Invariant #13 | primer3-py 2.3.1 |
| TH3 | Fold-back fraction 0.64 but no stable ntthal structure → no self-structure issue; Heuristic (fallback) screen: Primer3 alignment-mode self_any 9.00 / self_end 7.00 ≤ 12.00 → no issue either (audit round 2, A6) | Invariant #13 | primer3-py 2.3.1; dpal.c |
| TH4 | > 60-nt probe → fallback screens, ntthal fields null; (ACGT)16 self_any = self_end = 64.00 > PRIMER_INTERNAL_MAX_SELF_ANY 12.00 → issue | Invariant #13 | thal.c THAL_MAX_ALIGN; dpal.c |
| TH5 | Opt-in `ThermodynamicScreenMaxLength` = 64 / 100 (THAL_MAX_ALIGN override): (ACGT)16 hairpin 77.098245155727511, self-dimer = 3′ self-dimer 72.769499884640766 → both issues; 82-nt stem-loop hairpin 68.102992186803021 (issue), self-dimer 45.91673334496744 (none); 59 → `ArgumentOutOfRangeException` (A3-9, F55; `PrimerDesigner_NtthalMaxAlign_Tests.ValidateProbe_ThermodynamicScreenMaxLength_ScreensLongProbeWithNtthal`) | Invariant #13 | thal.c compiled with -DTHAL_MAX_ALIGN=10000 |
| KN7 | Opt-in THAL_MAX_ALIGN for the site duplex (A3-27, F56): 75-nt probe vs 75-nt site → null by default, 66.457703046655695 with `ThermodynamicScreenMaxLength` = 120; ≤ 60-nt sites (7 / 40 / 13 nt) computed by default = primer3-py calc_heterodimer 19.05924515571178 / 63.99555300270714 / −33.229679404433625; 62-nt probe / 49-nt site 36.11423712379826; `ValidateProbe` passes the conditions' value (hairpin 36.318725375690178, ANY 23.325175622075108, END1 17.579447109879538) | Invariant #16 | thal.c -DTHAL_MAX_ALIGN=10000; primer3-py 2.3.1 |
| HP1 | 74-nt probe, 10-bp perfect stem + 4-nt loop: fallback flags it (was missed); opt-in ntthal hairpin 77.95325865166825 °C (50 mM/50 nM) / 92.272558368972682 °C (Microarray 1 M/1 µM) | Invariant #19 | thal.c `-DTHAL_MAX_ALIGN=10000` |
| HP2 | Loops 3 / 4 / 8 nt all flagged | Invariant #19 | thal.c min_hrpn_loop |
| HP3 | 4-of-5 matched stem, no exact 4-bp stem → not flagged (80 % tolerance dropped) | Invariant #19 | F49 convention |
| HP4 | Fallback = HasHairpinPotential for 34 long / non-ACGT / Heuristic probes (ValidateProbe and DesignProbes warning) | Invariant #19 | canonical screen |
| HP5 | Fallback stem flag does not decide validity (A5-5): random 61-mer (F55 `Random61`, flagged by the stem screen, self_any / self_end ≤ 12) → HasSecondaryStructure true, `Warnings` = ["Potential secondary structure formation"], no issue, IsValid true; with the opt-in (`ThermodynamicScreenMaxLength` = 61) ntthal hairpin 33.980529935122263 °C ≤ 47 → valid. HP1 74-mer → HasSecondaryStructure true + warning; its only default issue is the sourced Primer3 self_any 16.00 > 12.00 (independent Smith–Waterman oracle 16); opt-in (100) → invalid by "ntthal hairpin Tm 78.0°C exceeds 47°C" | Invariant #14/#19 | thal.c `-DTHAL_MAX_ALIGN=10000`; F52 validity rule |
| KN1 | Kane fixtures A–E: score / identical / LCS = Biopython | Invariant #15 | Biopython 1.88 |
| KN2 | Kane thresholds strict (0.80 ↛ > 0.80; 15 nt ↛ > 15) | Invariant #15 | Kane et al. 2000 |
| KN3 | Chunked long non-target score = canonical whole-strand LocalAlignAffine | Invariant #15 | SequenceAligner |
| KN4 | Site duplex Tm = primer3-py calc_heterodimer (36.11423712379826, 66.04038852959525); maxDuplexTm flag | Invariant #16 | primer3-py; OligoArray 2.0 |
| KN5 | ValidateProbe with non-targets records Kane issues, IsValid false | Invariant #14/#15 | Kane et al. 2000 |
| KN6 | Argument guards (null/empty probe, identity ∉ [0,1], negative contiguous) | API | — |
| CS1 | CheckSpecificity bothStrands counts reverse-complement sites, palindromes once | API | blastn strand = both |
| KS1 | 12-mer: 3-mismatch site (9/12 = 0.75) not a Kane site → no off-target issue; + 2-mismatch site (10/12) → 2 Kane sites → issue | Invariant #17 | Kane et al. 2000; Python brute-force oracle |
| KS2 | 40-mer, radius 12: clustered 12-mismatch site (0.70, run 28) counts, spread one (0.70, run 5) not; thresholds 0.69 / 28 move the count 2 → 3 / 1 | Invariant #17 | Kane et al. 2000; Python brute-force oracle |
| KS3 | Guards: identity ∉ [0,1] / NaN, contiguous < 0, illegal stated conditions (Primer3 `_pr_data_control`) throw | API | Primer3 `libprimer3.cc` |
| BS1 | 20-mer + reference with exact site, rc site (2 subst., 0.90) and rc site (5 subst., 0.75, run 3): given strand (1,1); both strands radius 3 → (2,2) + issue, specificity 0.5; radius 5 → (3,2) | Invariant #18 | Python brute force over both strands; blastn 2.12.0+ `-strand plus` vs `both` |
| BS2 | Reverse-palindromic probe `GAATTCCGGAATTC` and reverse-palindromic site `GACGTCAGCTGACGTC` (probe 1 subst. from it, rc too) counted once on both strands | Invariant #18 | Python brute force; CheckSpecificity rule |
| HP1 | `ValidateProbe_LongProbeFallback_FlagsStemLoopWithFourNtLoop` | ✅ Covered | repro + thal.c Tm within 1e-6 |
| HP2 | `ValidateProbe_LongProbeFallback_AnyLoopOfAtLeastThreeIsFlagged` | ✅ Covered | 3 loops |
| HP3 | `ValidateProbe_LongProbeFallback_NoMismatchTolerance` | ✅ Covered | not flagged |
| HP4 | `FallbackHairpinScreen_EqualsHasHairpinPotential_ForLongNonAcgtAndHeuristicProbes`; `MolToolsCombinatorialTests.ProbeValid_SpecificityAndIssues_FollowValidationRules` (updated: GCGC·GCGC / ACGT·ACGT probes now flagged; A5-5: the flag is a warning, not an issue) | ✅ Covered | equality |
| HP5 | `ValidateProbe_LongProbeFallback_StemFlagIsWarning_NotValidityCriterion`, `ValidateProbe_LongProbeFallback_StemLoopRepro_StemFlagIsNotAnIssue_NtthalHairpinDecidesWithOptIn`; MCP `ValidateProbe_LongProbe_FallbackStemFlagIsWarning` | ✅ Covered | A5-5 |
| BS3 | 11 random cases (Python random.Random(7)) where the modes differ: (hits, Kane sites) for both modes = oracle | Invariant #18 | Python brute force over both strands |

### Should (Important)

| ID | Test Case | Rationale | Source |
|----|-----------|-----------|--------|
| S1 | Secondary structure potential detected for hairpin sequences | Quality metric | Implementation |
| S2 | Issues list populated for problematic probes | User feedback | Implementation |
| S3 | IsValid false when multiple issues exist | Validation logic | Implementation |
| S4 | Approximate matching with maxMismatches works correctly | Off-target detection | Wikipedia (Off-target) |
| SG1 | minIdentity threshold gates hits (0.8333 admitted at 0.75, rejected at 0.90) | Invariant #10 | Kane et al. 2000 |
| SG2 | Specific probe (1 on-target, 0 off-target) → IsSpecific true | API contract | Implementation |
| SG3 | Two identical perfect copies: one intended site (first in reference order), the other on-target-class and off-target; counts / IsSpecific independent of reference order (A5-3) | API contract | Implementation |

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
| SG3 | `ScanOffTargetsGapped_TwoPerfectCopies_CountsIndependentOfReferenceOrder` | ✅ Covered | both orders: on=2, off=1 (ref 1, start 5), IsSpecific false |
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
| TH1 | `ValidateProbe_ThermodynamicScreen_ReportsPrimer3NtthalTm` | ✅ Covered | primer3-py values within 1e-9 |
| TH2 | `ValidateProbe_ThermodynamicScreen_UsesStatedConditions` | ✅ Covered | primer3-py values within 1e-9 |
| TH3 | `ValidateProbe_HighFoldBackFractionWithoutStableStructure_PassesThermodynamicScreen` | ✅ Covered | thermo vs heuristic |
| TH4 | `ValidateProbe_ProbeLongerThan60nt_UsesPrimer3AlignmentSelfDimerFallback` | ✅ Covered | > 60 nt fallback (Primer3 alignment self-dimer) |
| KN1 | `AssessCrossHybridization_MatchesBiopythonLocalAlignmentAndLcs` | ✅ Covered | 10 strands |
| KN7 | `PrimerDesigner_NtthalMaxAlign_Tests.AssessCrossHybridization_MaxAlignOptIn_SiteDuplexTm_MatchesThal`, `ValidateProbe_ThermodynamicScreenMaxLength_AppliesToNonTargetDuplexTm`, `AssessCrossHybridization_SiteDuplexTm_MatchesPrimer3CalcHeterodimer` (62-nt probe) | ✅ Covered | thal.c / primer3-py within 1e-9 |
| KN2 | `AssessCrossHybridization_AppliesKaneThresholds` | ✅ Covered | strict thresholds |
| KN3 | `AssessCrossHybridization_LongNonTarget_ChunkedScoreEqualsCanonicalWholeStrandAlignment` | ✅ Covered | 9-kb strand |
| KN4 | `AssessCrossHybridization_SiteDuplexTm_MatchesPrimer3CalcHeterodimer` | ✅ Covered | calc_heterodimer within 1e-9 |
| KN5 | `ValidateProbe_WithNonTargets_RecordsKaneCrossHybridizationIssues` | ✅ Covered | issue text exact |
| KN6 | `AssessCrossHybridization_InvalidArguments_Throw` | ✅ Covered | guards |
| CS1 | `CheckSpecificity_BothStrands_CountsReverseComplementSites` | ✅ Covered | 0 / 1 / 0.5 |
| S3 | `ValidateProbe_MultipleProblems_IsValidFalse` (updated) | ✅ Covered | 3 issues: off-target + ntthal self-dimer 52.76 °C + hairpin 55.85 °C |
| KS1 | `ValidateProbe_ShortProbe_ThreeMismatchSiteAtSeventyFivePercent_IsNotAnOffTargetIssue` | ✅ Covered | hits 2/3, Kane sites 1/2, issue text exact |
| KS2 | `ValidateProbe_ReferenceSites_FollowKaneIdentityAndContiguityCriteria` | ✅ Covered | (1,1), (3,2), (3,3), (3,1) |
| KS3 | `ValidateProbe_InvalidKaneThresholdsOrConditions_Throw` | ✅ Covered | 7 guards |
| BS1 | `ValidateProbe_BothStrands_CountsReverseComplementSitesByKaneCriteria`; MCP `ValidateProbeTests.ValidateProbe_BothStrands_DelegatesReverseComplementScan` | ✅ Covered | (1,1) / (2,2) / (1,1) / (3,2) |
| BS2 | `ValidateProbe_BothStrands_PalindromicProbeOrSiteCountedOnce` | ✅ Covered | (1,1) in both modes |
| BS3 | `ValidateProbe_BothStrands_MatchesPythonBruteForce` (11 cases) | ✅ Covered | oracle tuples |

---

## Evidence-Backed Parameters

Every parameter is either externally sourced or marked as a library convention below.

| Parameter | Default | Evidence | Source |
|-----------|---------|----------|--------|
| `maxMismatches` | 3 | **Library convention** — search radius of the ungapped site scan (the figure comes from CRISPR/Cas9 guides — 3-5 bp mismatches per 20-nt guide, Hsu et al. 2013 — not from hybridization literature; kept for compatibility). Whether a found site is an issue is decided by the Kane criteria (identity > 0.75 or > 15-nt run); ⌈L/4⌉ − 1 reaches every site above 75 % identity. | library convention (Wikipedia: Off-target genome editing for the figure) |
| `SpecificityScore` | 1/N | **Library convention** — no published 1/N specificity score found (Kane 2000, OligoArray 2.0, Li & Stormo 2001 abstract, Primer3 PRIMER_INTERNAL_MAX_LIBRARY_MISHYB judge by identity / ΔG / alignment score); reported only, not used by `IsValid`. | library convention |
| `selfComplementarityThreshold` | 0.3 | Legacy, unused by the screen since audit round 2 (A6); the fold-back fraction is only reported. | library convention (compatibility) |
| `MaxSelfAny` / `MaxSelfEnd` (ProbeParameters) | 12.00 | Fallback self-dimer criterion (> 60 nt, non-ACGT, `Heuristic`): Primer3 alignment-mode PRIMER_INTERNAL_MAX_SELF_ANY / _SELF_END (dpal `self_any` / `self_end`). | Primer3 `libprimer3.cc` `pr_set_default_global_args_1`, `oligo_compl` |
| `MaxStructureTm` | 47 °C | PRIMER_INTERNAL_MAX_SELF_ANY_TH / _SELF_END_TH / _HAIRPIN_TH | Primer3 `libprimer3.cc` |
| `maxNonTargetIdentity` / `maxContiguousMatch` | 0.75 / 15 (strict >) | Kane et al. (2000) | Kane et al. 2000; Satya et al. 2008 |
| `maxDuplexTm` | none | OligoArray: user-set specificity threshold | Rouillard et al. 2003 |
| `bothStrands` (`ValidateProbe`) | false | Opt-in reverse-complement scan of the references (backward-compatible default; blastn defaults to `-strand both`) | blastn 2.12.0+ `-help`; audit round 4, A4-2 |
| `minIdentity` (gapped scan) | 0.75 | Kane et al. (2000): non-target transcripts >75% similar over the probe may cross-hybridize. Caller-configurable. | Kane et al. (2000), Nucleic Acids Res 28(22):4552 |
| `scoring` (gapped scan) | `SequenceAligner.BlastDna` (+2/−3, gap −2) | Reuses the BLAST-style DNA scoring already used for the library's gapped ANI alignment (COMPGEN-ANI). | Altschul et al. (1990); reused infrastructure |
| `k` (Karlin–Altschul) | 0.711 | Published nucleotide K for the +1/−3 scheme (NCBI blastn). K's full closed form needs the Karlin–Altschul score-lattice machinery, so it is a caller parameter; λ is computed (not assumed). | Karlin & Altschul (1990); NCBI blastn |
| `baseFrequency` (λ) | 0.25 | Uniform four-base background — the standard assumption used to derive the +1/−3 λ ≈ 1.374. | Karlin & Altschul (1990) |


