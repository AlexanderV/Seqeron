# Test Specification: PROTMOTIF-DOMAIN-001

**Test Unit ID:** PROTMOTIF-DOMAIN-001
**Area:** ProteinMotif
**Algorithm:** Protein Domain Identification (exact PROSITE PATTERN based)
**Status:** ☐ Not Started (reset for re-validation 2026-06-24 — exact-pattern fix)
**Owner:** Algorithm QA Architect
**Last Updated:** 2026-06-24

> **Exact-pattern fix (2026-06-24):** `FindDomains` now detects ONLY domains that have an EXACT,
> citable PROSITE **PATTERN**: zinc finger C2H2 (PS00028), WD-repeats (**PS00678**, newly made
> exact — replaced the prior ad-hoc `[LIVMFYWC]-x(5,12)-[WF]-D` regex), and Walker A / P-loop
> ATP-binding (PS00017). **SH3 (PS50002) and PDZ (PS50106) are weight-matrix PROFILES — no
> deterministic PROSITE pattern exists — so their previously shipped, unsourced ad-hoc regexes
> were removed (honest residual).** Pfam PF00018/PF00595 and the WD40 PF00400 full family are HMM
> profiles (trained models) and are not bundled.

> **Scope note (2026-06-24 validation):** This unit covers **`FindDomains`** only — domain
> identification via PROSITE-style consensus regexes. **Signal-peptide prediction was split out
> into its own unit `PROTMOTIF-SP-001`** (`tests/TestSpecs/PROTMOTIF-SP-001.md`,
> `ProteinMotifFinder_PredictSignalPeptide_Tests.cs`) and the implementation was rewritten to the
> von Heijne (1986) position-specific weight-matrix / EMBOSS `sigcleave` method
> (`PredictSignalPeptide(sequence, prokaryote, minWeight)`). The signal-peptide rows below
> (sources 1,2,8; evidence points 3–7; methods/invariants/test-cases referring to a tripartite
> n/h/c heuristic with `Score ∈ (0,1]`/`Probability`) describe the **superseded** combined model
> and are retained here only as historical context — they are validated under PROTMOTIF-SP-001,
> not by this unit's tests.

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
| 9 | PROSITE PS00678 — WD-repeats signature (WD_REPEATS_1) | 2 | https://prosite.expasy.org/PS00678 | 2026-06-24 |
| 10 | PROSITE pattern syntax (ScanProsite doc) | 2 | https://prosite.expasy.org/scanprosite/scanprosite_doc.html | 2026-06-24 |
| 11 | PROSITE PS50002 (SH3) / PS50106 (PDZ) — PROFILES, no pattern | 2 | https://prosite.expasy.org/PDOC50002 ; PDOC50106 | 2026-06-24 |
| 12 | UniProt P62873 (GBB1_HUMAN) — real WD40 β-propeller | 5 | https://rest.uniprot.org/uniprotkb/P62873.fasta | 2026-06-24 |

### 1.2 Key Evidence Points

1. C2H2 zinc finger consensus: `C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H` — PROSITE PS00028
2. Walker A motif: `[AG]-x(4)-G-K-[ST]` — PROSITE PS00017, Walker et al. (1982)
2a. WD-repeats signature (PROSITE PS00678, verbatim):
   `[LIVMSTAC]-[LIVMFYWSTAGC]-[LIMSTAG]-[LIVMSTAGC]-x(2)-[DN]-x-{P}-[LIVMWSTAC]-{DP}-[LIVMFSTAG]-W-[DEN]-[LIVMFSTAGCN]`
   → regex `[LIVMSTAC][LIVMFYWSTAGC][LIMSTAG][LIVMSTAGC].{2}[DN].[^P][LIVMWSTAC][^DP][LIVMFSTAG]W[DEN][LIVMFSTAGCN]`
   (14 elements, fixed **15** residues). Translation per ScanProsite rules: `-` dropped, `x`→`.`,
   `x(2)`→`.{2}`, `{P}`→`[^P]`, `{DP}`→`[^DP]`, `[..]` kept.
2b. SH3 (PS50002) and PDZ (PS50106) are PROSITE **PROFILES** (weight matrices) — no deterministic
   pattern → not reproducible as an exact regex → NOT detected (honest residual).
2c. Real WD40 positive: GBB1_HUMAN (P62873) matches PS00678 at 0-based starts 69, 156, 284.
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

1. **False positives from short patterns:** A deterministic PROSITE pattern is shorter and more permissive than the full Pfam HMM profile, so it may match sequences that do not fold into the expected domain. This is inherent to pattern-vs-HMM detection and is the honest residual: the exact PROSITE pattern is reproduced faithfully, but it is not equivalent to the trained HMM profile.
1a. **Profile-only domains absent:** SH3 (PS50002) and PDZ (PS50106) have no deterministic pattern; `FindDomains` does not detect them. This is intentional — fabricating a pattern would be a defect.
2. **Signal peptide scoring is evidence-based but simplified:** The 1:2:1 weighting (von Heijne 1985) and region scoring formulas are derived from literature statistics, not trained on specific datasets. Boundary cases may not match empirical Signal P or similar tools.

---

## 2. Canonical Methods Under Test

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `FindDomains(sequence)` | ProteinMotifFinder | Canonical | Domain detection via EXACT PROSITE PATTERN regexes (PS00028, PS00678, PS00017). Mapped from checklist `PredictDomains`. |
| `FindMotifByProsite(sequence, prositePattern, name)` | ProteinMotifFinder | Canonical (translation) | Used by M10 to verify the PS00678 PROSITE→regex translation end-to-end. |
| ~~`PredictSignalPeptide(...)`~~ | ProteinMotifFinder | — | **Moved to PROTMOTIF-SP-001** (von Heijne weight matrix). Not under test in this unit. |

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
| M7 | FindDomains_WD40Repeat_MatchesPrositePS00678 | Real GBB1 WD repeat (`LVSASQDGKLIIWDS`) padded → exact PS00678 match | Name="WD40 Repeat", Accession="PF00400", Start=4, End=18 (15 residues) | PROSITE PS00678 |
| M8 | FindDomains_WD40NearMiss_NoConservedTrp_ReturnsNoWD40 | Same segment with the invariant Trp (W) → A | No WD40 (empty) | PS00678 mandates literal W |
| M9 | FindDomains_GBB1Human_DetectsMultipleWD40Repeats | Full GBB1_HUMAN (P62873) sequence | WD40 hits at starts {69,156,284}, each 15 aa | PS00678; UniProt P62873 |
| M10 | FindMotifByProsite_PS00678_Translation_MatchesHandTracedSegment | Verbatim PS00678 string fed through PROSITE→regex translator vs hand-traced 15-mer | 1 match, Start=0, End=14, Sequence="LVSASQDGKLIIWDS" | PROSITE syntax rules + PS00678 |
| ~~M7s–M15s~~ | PredictSignalPeptide_* | **Superseded — moved to PROTMOTIF-SP-001** | — | von Heijne (1986) |
| M8 | PredictSignalPeptide_MinusOneMinusThreeRule | Cleavage at position where -1 and -3 are small amino acids | Detected; chars at -1,-3 ∈ {A,G,S} | von Heijne (1983) |
| M9 | PredictSignalPeptide_NoSignal_AllCharged | Fully charged sequence (no hydrophobic region) | Returns null | von Heijne (1986) |
| M10 | PredictSignalPeptide_ShortSequence | Sequence < 15 aa | Returns null | von Heijne (1986) |
| M11 | PredictSignalPeptide_NullInput | Null input | Returns null | Trivial |
| M12 | PredictSignalPeptide_EmptyInput | Empty string | Returns null | Trivial |
| M13 | PredictSignalPeptide_CaseInsensitive | Upper and lower case same sequence produce identical results | signalLower == signalUpper | Convention |
| M14 | PredictSignalPeptide_ScoreAndProbabilityInRange | Score ∈ (0,1], Probability = Score | INV-7, INV-8 | Scoring formula |
| M15 | PredictSignalPeptide_RejectsThreonineAtMinusThree | Sequence with T at -3 for only cleavage candidate | Returns null | von Heijne (1983) |

### 4.2 SHOULD Tests (Important edge cases)

| ID | Test Case | Description | Expected Outcome | Notes |
|----|-----------|-------------|------------------|-------|
| S4 | FindDomains_NoMatchingDomains | Random/short sequence | Empty result | Design |
| S5 | FindDomains_NoFabricatedSH3OrPDZ_ProfileOnlyDomains | Canonical Src SH3 core sequence | No SH3/PDZ domain reported | PS50002/PS50106 are profiles |
| S5 | PredictSignalPeptide_CleavagePositionRange | Cleavage ∈ [15, 35] on AlternateSignalPeptide | Exact cleavage=18, within bounds | Implementation |
| S6 | PredictSignalPeptide_HRegionMinLength | H-region must be ≥ 7 aa on AlternateSignalPeptide | Exact length=8, ≥ 7 | von Heijne (1985) |
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
| M5: FindDomains_DomainMetadata | ✅ Covered | Verifies Name, Accession, Description, Score across the 3 exact-pattern domain types |
| M6: FindDomains_StartLessOrEqualEnd | ✅ Covered | Invariant over all exact-pattern hits incl. GBB1 |
| M7: FindDomains_WD40Repeat_MatchesPrositePS00678 | ✅ Covered | Exact Start=4/End=18, 15-residue PS00678 window |
| M8: FindDomains_WD40NearMiss_NoConservedTrp | ✅ Covered | Removing invariant Trp abolishes match (empty) |
| M9: FindDomains_GBB1Human_DetectsMultipleWD40Repeats | ✅ Covered | Real β-propeller, hits {69,156,284}, each 15 aa |
| M10: FindMotifByProsite_PS00678_Translation | ✅ Covered | Verbatim PROSITE string → regex; hand-traced 15-mer |
| S4: FindDomains_NoMatch | ✅ Covered | Short random peptide |
| S5: FindDomains_NoFabricatedSH3OrPDZ | ✅ Covered | Src SH3 core not reported (profile-only) |
| S7: FindDomains_CaseInsensitive | ✅ Covered | Upper/lower identical (INV-10) |
| C1: FindDomains_MultipleDomainTypes | ✅ Covered | Both zinc finger and kinase |
| ~~PredictSignalPeptide_*~~ | — | **Moved to PROTMOTIF-SP-001** (not in this fixture) |

**Total: 14 tests in `ProteinMotifFinder_DomainPrediction_Tests.cs` — all passing (0 failures, 0 warnings).**
