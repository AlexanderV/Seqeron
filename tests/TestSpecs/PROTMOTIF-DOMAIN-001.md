# Test Specification: PROTMOTIF-DOMAIN-001

**Test Unit ID:** PROTMOTIF-DOMAIN-001
**Area:** ProteinMotif
**Algorithm:** Domain Prediction & Signal Peptide Prediction
**Status:** ☑ Complete
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-02-13

---

## 1. Evidence Summary

### 1.1 Authoritative Sources

| # | Source | Authority Rank | DOI or URL | Accessed |
|---|--------|---------------|------------|----------|
| 1 | von Heijne G (1986). Signal sequences. J Mol Biol | 1 | https://doi.org/10.1016/0022-2836(85)90046-4 | 2026-02-13 |
| 2 | von Heijne G (1983). Patterns near cleavage sites. Eur J Biochem | 1 | https://doi.org/10.1111/j.1432-1033.1983.tb07624.x | 2026-02-13 |
| 3 | Walker JE et al. (1982). ATP synthase subunits. EMBO J | 1 | https://doi.org/10.1002/j.1460-2075.1982.tb01276.x | 2026-02-13 |
| 4 | PROSITE PS00028 — Zinc finger C2H2 | 2 | https://prosite.expasy.org/PS00028 | 2026-02-13 |
| 5 | PROSITE PS00017 — P-loop (Walker A) | 2 | https://prosite.expasy.org/PS00017 | 2026-02-13 |
| 6 | Krishna et al. (2003). Zinc finger classification. NAR | 1 | https://doi.org/10.1093/nar/gkg161 | 2026-02-13 |
| 7 | Pfam PF00096, PF00400, PF00018, PF00595, PF00069 | 5 | https://www.ebi.ac.uk/interpro/ | 2026-02-13 |
| 8 | Owji et al. (2018). Signal peptides review. Eur J Cell Biol | 1 | https://doi.org/10.1016/j.ejcb.2018.06.003 | 2026-02-13 |

### 1.2 Key Evidence Points

1. C2H2 zinc finger consensus: `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` — PROSITE PS00028
2. Walker A motif: `[AG]-x(4)-G-K-[ST]` — PROSITE PS00017, Walker et al. (1982)
3. Signal peptide tripartite: n-region (charged), h-region (hydrophobic, 7–15 aa), c-region (cleavage) — von Heijne (1985, 1986)
4. (-1,-3) rule: Small amino acids {A, G, S} at positions -1 and -3 from cleavage — von Heijne (1983)
5. Signal peptide length: 16–30 aa — Owji et al. (2018)
6. H-region weighting: "both necessary and sufficient for membrane targeting" — von Heijne (1985)
7. N-region mean charge: ≈ +2.0 (eukaryotic +1.7, prokaryotic +2.3) — von Heijne (1986)

### 1.3 Documented Corner Cases

1. **Empty/null sequence:** Both methods must handle gracefully (return empty/null).
2. **Short sequence:** Signal peptide requires minimum ~15 aa; domain patterns require minimum pattern length.
3. **No matching pattern:** Charged-only sequences should have no signal peptide; random short peptides should have no domains.
4. **Case insensitivity:** Protein sequences may be upper or lower case.

### 1.4 Known Failure Modes / Pitfalls

1. **False positives from simplified patterns:** Regex-based domain detection may match spurious sequences that do not fold into the expected domain structure. — Inherent limitation of pattern vs. HMM approach.
2. **Signal peptide scoring is evidence-based but simplified:** The 1:2:1 weighting (von Heijne 1985) and region scoring formulas are derived from literature statistics, not trained on specific datasets. Boundary cases may not match empirical Signal P or similar tools.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindDomains(sequence)` | ProteinMotifFinder | Canonical | Domain detection via regex patterns. Mapped from checklist `PredictDomains`. |
| `PredictSignalPeptide(sequence, maxLength)` | ProteinMotifFinder | Canonical | Signal peptide prediction with tripartite scoring. |

---

## 3. Invariants

| ID | Invariant | Verifiable | Evidence |
|----|-----------|------------|----------|
| INV-1 | FindDomains on null/empty → yields nothing | Yes | Trivial |
| INV-2 | FindDomains returns ProteinDomain with non-empty Name, Accession | Yes | Pfam |
| INV-3 | FindDomains Start ≤ End for every result | Yes | Trivial |
| INV-4 | PredictSignalPeptide on null/empty → null | Yes | Trivial |
| INV-5 | PredictSignalPeptide on short (<15 aa) → null | Yes | von Heijne (1986) |
| INV-6 | PredictSignalPeptide cleavage position ∈ [15, 35] | Yes | Implementation constraint |
| INV-7 | PredictSignalPeptide Score ∈ (0, 1] | Yes | Scoring formula |
| INV-8 | PredictSignalPeptide Probability = Score | Yes | Direct quality measure |
| INV-9 | PredictSignalPeptide returns non-empty N, H, C regions when not null | Yes | von Heijne (1986) |
| INV-10 | Case insensitivity: upper and lower input yield same result | Yes | Standard convention |
| INV-11 | H-region length ≥ 7 | Yes | von Heijne (1985) |
| INV-12 | Positions -1 and -3 from cleavage are in {A, G, S} | Yes | von Heijne (1983) |

---

## 4. Test Cases

### 4.1 MUST Tests (Required — every row needs Evidence)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| M1 | FindDomains_ZincFinger_MatchesPrositeConsensus | Sequence matching PS00028 pattern yields Zinc Finger C2H2 domain | Domain with Name="Zinc Finger C2H2", Accession="PF00096" | PROSITE PS00028 |
| M2 | FindDomains_WalkerA_MatchesKinaseDomain | Sequence matching PS00017 pattern yields Protein Kinase ATP-binding | Domain with Accession="PF00069" | PROSITE PS00017, Walker (1982) |
| M3 | FindDomains_EmptySequence_ReturnsEmpty | Empty string input | Empty enumerable | Trivial |
| M4 | FindDomains_NullSequence_ReturnsEmpty | Null input | Empty enumerable | Trivial |
| M5 | FindDomains_DomainMetadata_HasCorrectFields | All returned domains have non-empty Name, Accession, Description | Assert.Multiple on fields | Pfam |
| M6 | FindDomains_StartLessOrEqualEnd | Every domain has Start ≤ End | INV-3 | Trivial |
| M7 | PredictSignalPeptide_TripartiteStructure | Signal peptide with classic n/h/c regions detected | Non-null result with NRegion, HRegion, CRegion | von Heijne (1986) |
| M8 | PredictSignalPeptide_MinusOneMinusThreeRule | Cleavage at position where -1 and -3 are small amino acids | Detected; chars at -1,-3 ∈ {A,G,S} | von Heijne (1983) |
| M9 | PredictSignalPeptide_NoSignal_AllCharged | Fully charged sequence (no hydrophobic region) | Returns null | von Heijne (1986) |
| M10 | PredictSignalPeptide_ShortSequence | Sequence < 15 aa | Returns null | von Heijne (1986) |
| M11 | PredictSignalPeptide_NullInput | Null input | Returns null | Trivial |
| M12 | PredictSignalPeptide_EmptyInput | Empty string | Returns null | Trivial |
| M13 | PredictSignalPeptide_CaseInsensitive | Upper and lower case same sequence produce identical results | signalLower == signalUpper | Convention |
| M14 | PredictSignalPeptide_ScoreAndProbabilityInRange | Score ∈ (0,1], Probability = Score | INV-7, INV-8 | Scoring formula |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S1 | FindDomains_WD40Repeat | WD40 pattern match | Domain detected | Pfam PF00400 |
| S2 | FindDomains_SH3Domain | SH3 signature match | Domain detected | Pfam PF00018 |
| S3 | FindDomains_PDZDomain | PDZ pattern match | Domain detected | Pfam PF00595 |
| S4 | FindDomains_NoMatchingDomains | Random/short sequence | Empty result | Design |
| S5 | PredictSignalPeptide_CleavagePositionRange | Cleavage ∈ [15, 35] | Within bounds | Implementation |
| S6 | PredictSignalPeptide_HRegionMinLength | H-region must be ≥ 7 aa | Verified in result | von Heijne (1985) |
| S7 | FindDomains_CaseInsensitive | Upper and lower case produce identical domains | domainsLower == domainsUpper | INV-10, Convention |

### 4.3 COULD Tests (Nice to have)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| C1 | FindDomains_MultipleDomainTypes | Sequence containing both zinc finger and kinase motifs | Both detected | Architecture |
| C2 | PredictSignalPeptide_MaxLengthParameter | Custom maxLength limits search range | Affects result | API behavior |

---

## 5. Audit of Existing Tests

### 5.1 Discovery Summary

- Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_DomainPrediction_Tests.cs`
- Integration tests remain in `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinderTests.cs`

### 5.2 Coverage Classification

| Area / Test Case ID | Status | Resolution |
|---------------------|--------|------------|
| M1: FindDomains_ZincFinger | ✅ Covered | Exact values, assertion messages, evidence citations |
| M2: FindDomains_WalkerA (P-loop) | ✅ Covered | Exact values, assertion messages, evidence citations |
| M3: FindDomains_EmptySequence | ✅ Covered | Assertion message |
| M4: FindDomains_NullSequence | ✅ Covered | Implemented |
| M5: FindDomains_DomainMetadata | ✅ Covered | Verifies Name, Accession, Description, Score |
| M6: FindDomains_StartLessOrEqualEnd | ✅ Covered | Invariant over all 5 domain pattern types |
| M7: PredictSignalPeptide_TripartiteStructure | ✅ Covered | Exact score ±1e-10, explicit region values |
| M8: PredictSignalPeptide_MinusOneMinusThreeRule | ✅ Covered | Verifies -1,-3 positions ∈ {A,G,S} |
| M9: PredictSignalPeptide_NoSignal | ✅ Covered | All-charged sequence |
| M10: PredictSignalPeptide_ShortSequence | ✅ Covered | Below 15 aa minimum |
| M11: PredictSignalPeptide_NullInput | ✅ Covered | Implemented |
| M12: PredictSignalPeptide_EmptyInput | ✅ Covered | Implemented |
| M13: PredictSignalPeptide_CaseInsensitive | ✅ Covered | Tolerance ±1e-10 for identical computations |
| M14: PredictSignalPeptide_ScoreRange | ✅ Covered | Range (0,1] + Probability == Score (INV-8) |
| S1: FindDomains_WD40 | ✅ Covered | Exact position and name |
| S2: FindDomains_SH3 | ✅ Covered | Exact position and name |
| S3: FindDomains_PDZ | ✅ Covered | Exact position and name |
| S4: FindDomains_NoMatch | ✅ Covered | Short random peptide |
| S5: PredictSignalPeptide_CleavageRange | ✅ Covered | Invariant [15, 35] |
| S6: PredictSignalPeptide_HRegionMinLength | ✅ Covered | ≥ 7 aa per von Heijne (1985) |
| S7: FindDomains_CaseInsensitive | ✅ Covered | Upper/lower identical (INV-10) |
| C1: FindDomains_MultipleDomainTypes | ✅ Covered | Both zinc finger and kinase |
| C2: PredictSignalPeptide_MaxLengthParameter | ✅ Covered | Restricts search range |

**Total: 23 tests — all passing (0 failures, 0 warnings)**
