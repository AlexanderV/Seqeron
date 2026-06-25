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

---

## Addendum 2026-06-25 — Plan7 profile-HMM engine + bundled Pfam SH3/PDZ/WD40 (opt-in)

Limitation fix: SH3 (PF00018), PDZ (PF00595) and WD40 (PF00400) — which have **no deterministic
PROSITE pattern** — are now detectable via an **opt-in** Plan7 profile-HMM scorer. The exact
PROSITE-pattern `FindDomains` path and its defaults are unchanged. Evidence: Evidence-doc addendum
"Plan7 profile-HMM engine + bundled Pfam SH3/PDZ/WD40".

### A.1 Canonical Methods Under Test (HMM)

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `Plan7ProfileHmm.Parse` | Plan7ProfileHmm | Canonical | HMMER3/f parser |
| `Plan7ProfileHmm.ViterbiScore` | Plan7ProfileHmm | Canonical | glocal log-odds Viterbi (nats) |
| `Plan7ProfileHmm.ForwardScore` | Plan7ProfileHmm | Canonical | glocal log-odds Forward (nats) |
| `ProteinMotifFinder.FindDomainsByHmm` | ProteinMotifFinder | Canonical | SH3/PDZ/WD40 detection |
| `ProteinMotifFinder.ScoreDomainHmm` | ProteinMotifFinder | Canonical | one-profile bit score |

### A.2 Invariants (HMM) — see algorithm doc INV-HMM-01..04

### A.3 MUST Tests (HMM)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| H1 | Exact DP on hand-built HMM | tiny 2-match-state HMM, seq "AC", path B→M1→M2→E | Viterbi = 0.5187937934151676 nats (Within 1e-9) | Durbin §5.4 (Evidence addendum hand-HMM dataset) |
| H2 | Forward ≥ Viterbi | same hand HMM | ForwardScore ≥ ViterbiScore | INV-HMM-01 |
| H3 | SH3 true positive | SRC_HUMAN SH3 core vs PF00018 | bit score > 10 (≈ +60) | UniProt P12931; Pfam PF00018 |
| H4 | PDZ true positive | DLG4_HUMAN PDZ1 (res 61–151) vs PF00595 | bit score > 10 (≈ +83) | UniProt P78352; Pfam PF00595 |
| H5 | WD40 true positive | GBB1_HUMAN vs PF00400 | bit score > 10 (≈ +36) | UniProt P62873; Pfam PF00400 |
| H6 | True negative rejected | low-complexity A14E14K12 vs all 3 | bit score < 0; not reported by FindDomainsByHmm | INV-HMM-02 |
| H7 | Cross-domain specificity | SH3 core vs PF00400 < its PF00018 score | strongly negative for the wrong profile | INV-HMM-02 |
| H8 | `.hmm` parser round-trip | parse embedded PF00018 | Name=SH3_1, Acc=PF00018.35, Length=48, GA=22.9 | HMMER3/f header; PF00018.35 |
| H9 | Parser handles '*' | tiny HMM with a `*` on the only path | path forbidden (−∞ score) | INV-HMM-04; HMMER guide |
| H10 | FindDomainsByHmm assigns correct family | SH3 core / PDZ1 / GBB1 each detected as own family | one matching ProteinDomain each | Pfam |
| H11 | Determinism | re-score same input twice | identical scores | INV-HMM-03 |
| H12 | Null/empty + unknown accession | guards | empty / ArgumentNullException / ArgumentException / FormatException | contract §3.3 |

### A.4 Honest residual

Exact `hmmsearch` bit-score / E-value parity is NOT reproduced (filter pipeline, null2, Gumbel
calibration out of scope); the DP is verified exactly on the hand HMM (1e-9) and by correct
ranking on real true/false positives. Only 3 Pfam domains bundled. This keeps the unit's status
as an honest, partial-but-verified fix (Status remains ☐ in the root registry).

### A.5 Post-implementation coverage (HMM)

Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs`
(17 tests). All H1–H12 implemented and ✅ green:

| Case | Test method | Status |
|------|-------------|--------|
| H1 | ViterbiScore_HandBuiltHmm_MatchesExactDurbinDerivation | ✅ |
| H2 | ForwardScore_HandBuiltHmm_IsAtLeastViterbiScore | ✅ |
| H3 | ScoreDomainHmm_Sh3TruePositive_ScoresAboveThreshold | ✅ |
| H4 | ScoreDomainHmm_PdzTruePositive_ScoresAboveThreshold | ✅ |
| H5 | ScoreDomainHmm_Wd40TruePositive_ScoresAboveThreshold | ✅ |
| H6 | ScoreDomainHmm_TrueNegative_ScoresBelowZeroForAllProfiles; FindDomainsByHmm_TrueNegative_ReportsNoDomain | ✅ |
| H7 | ScoreDomainHmm_CrossDomain_TruePositiveScoresHigherOnItsOwnFamily | ✅ |
| H8 | Parse_EmbeddedSh3Profile_ReadsHeaderFieldsExactly | ✅ |
| H9 | ViterbiScore_StarOnOnlyPath_ForbidsPath_ReturnsNegativeInfinity | ✅ |
| H10 | FindDomainsByHmm_Sh3TruePositive_DetectsOnlySh3; FindDomainsByHmm_Wd40TruePositive_DetectsWd40 | ✅ |
| H11 | ScoreDomainHmm_RepeatedCalls_AreDeterministic | ✅ |
| H12 | FindDomainsByHmm_NullOrEmpty_ReturnsEmpty; ScoreDomainHmm_NullSequence_Throws; ScoreDomainHmm_UnknownAccession_Throws; Parse_NonHmmer3Text_ThrowsFormatException | ✅ |

Full unfiltered suite green (all projects Failed: 0). Branch coverage on `Plan7ProfileHmm` /
`FindDomainsByHmm` / `ScoreDomainHmm` ≥ 80% (parser header/COMPO/`*`/node lines, Viterbi, Forward,
both detection methods, and all guard paths exercised).

---

## Addendum 2026-06-25 — HMMER E-value / P-value statistics from profile STATS (opt-in)

Adds the Gumbel (MSV/Viterbi) and exponential-tail (Forward) P-value and `E = P·Z` E-value layer.
Detection and defaults unchanged. Evidence: Evidence-doc addendum "HMMER E-value / P-value
statistics from profile STATS".

### B.1 Canonical Methods Under Test (statistics)

| Method | Class | Type | Notes |
|--------|-------|------|-------|
| `Plan7ProfileHmm.Statistics` (STATS parse) | Plan7ProfileHmm | Canonical | parses `STATS LOCAL MSV/VITERBI/FORWARD` µ/λ/τ; null if uncalibrated |
| `Plan7ProfileHmm.GumbelSurvival` / `ExponentialSurvival` | Plan7ProfileHmm | Canonical | Easel `esl_gumbel_surv` / `esl_exp_surv` |
| `Plan7ProfileHmm.ViterbiPValue/MsvPValue/ForwardPValue` | Plan7ProfileHmm | Canonical | survival at the profile's µ/λ/τ |
| `Plan7ProfileHmm.ViterbiEValue/ForwardEValue` / `EValue` | Plan7ProfileHmm | Canonical | `E = P·Z` |
| `ProteinMotifFinder.FindDomainHitsByHmm` | ProteinMotifFinder | Canonical | detection + E-value hit |
| `ProteinMotifFinder.ScoreDomainHmmEValue` | ProteinMotifFinder | Canonical | bit score + E-value |

### B.2 Invariants — see algorithm doc INV-HMM-05 (E monotone in S, linear in Z)

### B.3 MUST Tests (statistics)

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| H13 | STATS parsed | parse PF00018 | Statistics = (MSV −8.1284/0.71923, VIT −8.2932/0.71923, FWD −4.5735/0.71923); hand-HMM → null + ViterbiPValue throws | HMMER guide STATS line; `Resources/PF00018_SH3_1.hmm` |
| H14 | Exact Gumbel/exponential P + E | SH3 STATS, S=40, Z=1000 | VIT P=8.227179545686635e-16, E=8.227…e-13; FWD P=1.1943390031599535e-14, E=1.1943…e-11 (Within 1e-9 rel); FWD below τ → P=1 | Easel esl_gumbel_surv/esl_exp_surv; Eddy 2008; hand-derived |
| H15 | Monotonicity / Z-scaling | E(40)<E(20); E(Z=1000)=1000·E(Z=1); Z<0 throws | INV-HMM-05 | HMMER guide `E=P·Z` |
| H16 | End-to-end E-value | SH3 true positive E≪1e-3; negative E>1 (≈Z); FindDomainHitsByHmm reports hit with E + bit Score; negative/empty/null → no hit | UniProt P12931; low-complexity negative; `E=P·Z` |

### B.4 Honest residual (narrowed)

The Gumbel/exponential P-value and `E = P·Z` are now implemented exactly (and read the profile's own
`STATS LOCAL` calibration). What remains out of scope is exact `hmmsearch`-**reported** E-value
*pipeline* parity: HMMER applies these formulas to its local-multihit sequence bit score after the
MSV/bias prefilters and the **null2 biased-composition correction**, which this glocal scorer does
not compute, so absolute bit scores (and hence absolute reported E-values) differ from `hmmsearch`.
Pfam coverage beyond the three bundled (caller-supplied `.hmm`) profiles is likewise out of scope.
Status remains ☐ in the root registry (independent re-validation).

### B.5 Post-implementation coverage (statistics)

Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs`
(now 33 tests). All H13–H16 implemented and ✅ green:

| Case | Test method(s) | Status |
|------|----------------|--------|
| H13 | Parse_EmbeddedSh3Profile_ReadsStatsLinesExactly; Statistics_UncalibratedProfile_IsNull_AndPValueThrows | ✅ |
| H14 | ViterbiPValue_Sh3At40Bits_MatchesHandDerivedGumbel; ViterbiEValue_Sh3At40Bits_Z1000_MatchesHandDerived; ForwardPValue_Sh3At40Bits_MatchesHandDerivedExponential; ForwardEValue_Sh3At40Bits_Z1000_MatchesHandDerived; ForwardPValue_BelowTau_IsClampedToOne; GumbelSurvival_PureFormula_MatchesHandDerived | ✅ |
| H15 | ViterbiEValue_DecreasesAsScoreIncreases; ViterbiEValue_ScalesLinearlyWithDatabaseSize; EValue_NegativeDatabaseSize_Throws | ✅ |
| H16 | ScoreDomainHmmEValue_Sh3TruePositive_HasTinyEValue; ScoreDomainHmmEValue_TrueNegative_HasLargeEValue; FindDomainHitsByHmm_Sh3TruePositive_ReportsHitWithEValue; FindDomainHitsByHmm_TrueNegative_ReportsNoHit; FindDomainHitsByHmm_NullOrEmpty_ReturnsEmpty | ✅ |

Full unfiltered suite green (all projects Failed: 0). Work Queue Remaining = 0.

---

## Addendum 2026-06-25 — HMMER local-multihit Forward + null2 (hmmsearch parity, opt-in)

This addendum records the test spec for the opt-in `hmmsearch`-parity layer. The glocal path and all
existing defaults/tests are unchanged. A real HMMER reference was obtained via **pyhmmer 0.12.1**.

### C.1 Canonical Methods Under Test

- `Plan7ProfileHmm.LocalForwardScore(seq)` — local-multihit Forward score (nats).
- `Plan7ProfileHmm.LocalForwardBitScore(seq)` — `(fwd − nullsc)/ln2` = HMMER `pre_score`.
- `Plan7ProfileHmm.Null2BiasBits(seq)` — null2 biased-composition correction (bits) = HMMER `bias`.
- `Plan7ProfileHmm.HmmSearchBitScore(seq)` — `pre_score − bias` = HMMER reported per-seq score.

### C.2 Invariants

- INV-HMM-06: `HmmSearchBitScore == LocalForwardBitScore − Null2BiasBits` (pipeline identity).
- INV-HMM-07: `Null2BiasBits ≥ 0` (seqbias = logsumexp(0, …) never raises the score); low-complexity
  sequences incur at least as large a correction as a real domain.
- INV-HMM-08: the opt-in local path does not change the glocal `ForwardScore` (regression).

### C.3 MUST Tests

| ID | Test Case | Description | Expected Outcome | Evidence |
|----|-----------|-------------|------------------|----------|
| H18a | Hand-built local-mode pin | 1-node HMM emitting A, B→M1=1, seq "A" | LocalForwardScore = 1.272400756045032 nats; pre = 3.835686260769536 bits (Within 1e-12) | hand derivation from modelconfig.c/generic_fwdback.c |
| H18b | SH3 pre_score parity | SRC_HUMAN core vs PF00018 | LocalForwardBitScore = 68.709740 (Within 1e-4) | pyhmmer 0.12.1 hmmsearch `pre_score` |
| H18c | PDZ pre_score parity | PSD-95 PDZ1 vs PF00595 | = 84.862930 (Within 1e-4) | pyhmmer 0.12.1 |
| H18d | WD40 pre_score parity | GBB1 β-propeller vs PF00400 | = 213.411926 (Within 1e-4) | pyhmmer 0.12.1 |
| H18e | null2 bias parity | SH3 domain envelope (pos 3–50) vs PF00018 | Null2BiasBits = 0.025574 (Within 1e-4) | pyhmmer 0.12.1 reported `bias` |
| H18f | Pipeline identity | SH3 | HmmSearchBitScore = pre − bias (Within 1e-9) | p7_pipeline.c |
| H18g | Bias non-negative | SH3 vs low-complexity | both ≥ 0; low-complexity ≥ real-domain bias | seqbias = logsumexp(0,·) |
| H18h | Glocal unchanged | SH3 | ForwardScore = 41.78685952002655; LocalForwardScore = 42.609594871580114 (≠ glocal) | regression |
| H18i | Guards | null/empty | LocalForwardScore(null) throws; HmmSearchBitScore("") = −∞; Null2BiasBits("") = 0 | contract |

### C.4 Honest residual

`pre_score` (local-multihit Forward) and the null2 correction now reproduce HMMER (pyhmmer 0.12.1) to
single-precision rounding. The remaining residual is HMMER's automatic **multi-domain envelope
decomposition** (region detection + stochastic-traceback clustering): null2 is applied over the
caller-supplied sequence/envelope, so a single well-resolved domain matches `hmmsearch`'s corrected
score, but a multi-domain target must be scored per-envelope. The MSV/bias prefilters (which only
gate which sequences reach Forward, not a hit's score) and the full Pfam library remain out of scope.
Status remains ☐ in the root registry (independent re-validation).

### C.5 Post-implementation coverage

Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs`
(now 43 tests). All H18 implemented and ✅ green:

| Case | Test method(s) | Status |
|------|----------------|--------|
| H18a | LocalForwardScore_HandBuiltOneNode_MatchesExactDerivation | ✅ |
| H18b–d | LocalForwardBitScore_Sh3/Pdz/Wd40TruePositive_MatchesHmmsearchReference | ✅ |
| H18e | Null2BiasBits_Sh3DomainEnvelope_MatchesHmmsearchReportedBias | ✅ |
| H18f | HmmSearchBitScore_EqualsPreScoreMinusNull2Bias | ✅ |
| H18g | Null2BiasBits_IsNonNegative_AndPositiveForBiasedComposition | ✅ |
| H18h | LocalForwardScore_DoesNotChangeGlocalForwardScore | ✅ |
| H18i | LocalForwardScore_NullOrEmpty_GuardsCorrectly; HmmSearchBitScore_Sh3TruePositive_IsHighlyPositive | ✅ |

Full unfiltered suite green (all projects Failed: 0). Work Queue Remaining = 0.

---

## Addendum 2026-06-25 — HMMER multi-domain envelope decomposition (p7_domaindef, opt-in)

This addendum adds HMMER's automatic per-target **domain/envelope decomposition** as an opt-in
extension. Existing methods + defaults unchanged. Evidence: the same Evidence file's
"Addendum 2026-06-25 — HMMER multi-domain envelope decomposition (p7_domaindef)".

### D.1 Evidence summary

`p7_domaindef.c` region identification (`rt1=0.25`, `rt2=0.10`), `is_multidomain_region` (`rt3=0.20`),
`rescore_isolated_domain` (unihit Forward over the envelope at full length n + null2 by expectation);
`generic_decoding.c` `p7_GDomainDecoding` (`btot`/`etot`/`mocc`); `p7_pipeline.c` per-domain bit
score `(envsc + (n−Ld)·ln(n/(n+3)) − (nullsc + dombias))/ln 2`,
`dombias = logsumexp(0, ln(1/256)+domcorrection)`, i-Evalue `= Z·exp(−λ(score−τ))`. GROUND TRUTH:
pyhmmer 0.12.1 `hmmsearch` — GBB1_HUMAN (P62873, L=340) vs PF00400 → **7** domains; SRC_HUMAN SH3
(P12931, L=55) vs PF00018 → **1** domain.

### D.2 Canonical methods under test

- `Plan7ProfileHmm.FindDomains(string)` → `IReadOnlyList<DomainEnvelope>` (envelope coords + score + bias + i-Evalue).
- `ProteinMotifFinder.FindDomainEnvelopes(string[, minBitScore])` / `FindDomainEnvelopes(string, accession)`.

### D.3 Invariants

- INV-HMM-09: domain count and per-envelope `EnvelopeStart`/`EnvelopeEnd` reproduce hmmsearch
  (well-separated domains).
- INV-HMM-10: per-domain `BitScore` reproduces hmmsearch's per-domain score to single precision
  (~1e-2 bits, float32 vs float64); `IndependentEValue` to a 5% band.
- INV-HMM-11: a single well-resolved domain → exactly one envelope; its score ≈ the per-sequence
  `HmmSearchBitScore`.

### D.4 MUST cases

| ID | Scenario | Input | Expected (Evidence) | Source |
|----|----------|-------|---------------------|--------|
| H19a | Domain count (multi) | GBB1 vs PF00400 | `FindDomains` count = 7 | pyhmmer 0.12.1 hmmsearch |
| H19b | Envelope coords | GBB1 vs PF00400 | env = {45-83,87-125,133-170,174-212,216-254,259-298,303-340} exactly | pyhmmer 0.12.1 |
| H19c | Per-domain scores | GBB1 vs PF00400 | {31.139467,19.004278,25.053679,35.552242,40.454269,23.443121,27.824228} (Within 1e-2) | pyhmmer 0.12.1 |
| H19d | Per-domain i-Evalues | GBB1 vs PF00400 | {1.21e-11,8.41e-08,1.02e-09,4.85e-13,1.36e-14,3.31e-09,1.36e-10} (Within 5%) | pyhmmer 0.12.1 |
| H19e | Single domain | SH3 core vs PF00018 | 1 envelope, env 3-50, score 68.540695 (Within 1e-2), i-Evalue 1.4529e-23 (Within 5%) | pyhmmer 0.12.1 |
| H19f | Wrapper ordering | GBB1 via FindDomainEnvelopes | 7 WD40 hits, ascending start, family "WD40" | derived |
| H19g | Empty | "" | no envelopes (both APIs) | contract |
| H19h | Null | null | ArgumentNullException | contract |
| H19i | Unknown accession | (seq, "PF99999") | ArgumentException | only 3 bundled |
| H19j | Envelope vs per-seq | SH3 core | single-envelope score ≈ HmmSearchBitScore (Within 0.5) | INV-HMM-11 |

### D.5 Post-implementation coverage

Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs`
(now 53 tests). All H19 implemented and ✅ green:

| Case | Test method | Status |
|------|-------------|--------|
| H19a | FindDomains_Gbb1Wd40Propeller_RecoversSevenEnvelopes | ✅ |
| H19b | FindDomains_Gbb1Wd40Propeller_EnvelopeCoordinatesMatchHmmsearch | ✅ |
| H19c | FindDomains_Gbb1Wd40Propeller_PerDomainBitScoresMatchHmmsearch | ✅ |
| H19d | FindDomains_Gbb1Wd40Propeller_PerDomainIEvaluesMatchHmmsearch | ✅ |
| H19e | FindDomains_SingleSh3Domain_YieldsOneEnvelopeMatchingHmmsearch | ✅ |
| H19f | FindDomainEnvelopes_Gbb1_ReportsSevenWd40HitsInOrder | ✅ |
| H19g | FindDomains_EmptySequence_ReturnsNoEnvelopes | ✅ |
| H19h | FindDomains_NullSequence_Throws | ✅ |
| H19i | FindDomainEnvelopes_UnknownAccession_Throws | ✅ |
| H19j | FindDomains_SingleDomainEnvelope_ScoreConsistentWithPerSequence | ✅ |

### D.6 Residual

Stochastic-traceback clustering (`region_trace_ensemble`) for regions the `rt3` test flags as
closely-overlapping multi-domain is NOT implemented; such a region is emitted as a single envelope.
Verified path = well-separated domains (tandem repeats, multi-domain propellers). Full unfiltered
suite green. Work Queue Remaining = 0.

---

## Addendum 2026-06-25 — HMMER stochastic-traceback clustering of overlapping domains (`region_trace_ensemble`, opt-in)

This addendum closes the prior D.6 residual: the stochastic-traceback clustering that splits a
*closely-overlapping* multi-domain region is now implemented as the default path of
`Plan7ProfileHmm.FindDomains(seq, clusterOverlapping=true)`. Defaults of all other methods unchanged.

### E.1 Evidence summary

HMMER `region_trace_ensemble` (`p7_domaindef.c`) + `p7_spensemble_Cluster` (`p7_spensemble.c`) +
`p7_GStochasticTrace` (`generic_stotrace.c`) + `p7_trace_Index` (`p7_trace.c`) +
`esl_cluster_SingleLinkage` (`esl_cluster.c`) + the Easel LCG RNG (`esl_random.c` /`esl_mix3` in
`easel.c`) — all retrieved verbatim 2026-06-25 (see Evidence Addendum, refs 31/34–39). Defaults
(verbatim): `nsamples=200`, `min_overlap=0.8`, `of_smaller=TRUE`, `max_diagdiff=4`,
`min_posterior=0.25`, `min_endpointp=0.02`; pipeline RNG `esl_randomness_CreateFast(42)`,
`do_reseeding=TRUE`. Ground truth: pyhmmer 0.12.1 `hmmsearch` (Z=1, domZ=1, seed=42) on closely-
overlapping tandem-SH3 constructs (a truncated SH3 core + a full SH3 core).

### E.2 Canonical method

`Plan7ProfileHmm.FindDomains(string sequence, bool clusterOverlapping = true)` — when a region is
flagged multi-domain by `is_multidomain_region` (`rt3=0.20`), it is resolved into one or more
envelopes by the stochastic-traceback ensemble; `clusterOverlapping:false` keeps the prior single-
envelope behaviour.

### E.3 Invariant

`INV-HMM-12` — a region the `rt3` test flags multi-domain is split into the consensus envelopes of
the trace ensemble; well-separated regions and single domains are unaffected; the ensemble is
deterministic (fixed-seed LCG, reproducible across runs).

### E.4 MUST tests

| ID | Scenario | Input | Expected (Evidence) | Source |
|----|----------|-------|---------------------|--------|
| H20a | Overlapping split (primary) | SH3 core trim=12 + full core (L=84) | 2 envelopes, env 1-37 & 37-84 (overlap at 37); scores 48.047169 & 66.678467 (Within 0.1) | pyhmmer 0.12.1 hmmsearch |
| H20b | Opt-in / opt-out | same (L=84) | clusterOverlapping=false → 1 envelope; =true → 2 | region_trace_ensemble branch |
| H20c | Coords across trims | trims 4/12/16 | {1-46,45-92},{1-37,37-84},{1-33,33-80} exactly | pyhmmer 0.12.1 |
| H20d | Determinism | L=84, two calls | identical envelopes | fixed-seed LCG reseed |
| H20e | Well-separated regression | two full cores (L=96) | 2 envelopes 1-48 & 49-96 (no overlap, flank-bound split) | pyhmmer 0.12.1 |
| H20f | Single-domain regression | SH3 core (L=55) | 1 envelope, env 3-50 | pyhmmer 0.12.1 |

### E.5 Post-implementation coverage

Canonical file: `tests/Seqeron/Seqeron.Genomics.Tests/ProteinMotifFinder_FindDomainsByHmm_Tests.cs`
(now 59 tests). All H20 implemented and ✅ green:

| Case | Test method | Status |
|------|-------------|--------|
| H20a | FindDomains_CloselyOverlappingTandemSh3_SplitsIntoTwoEnvelopesViaEnsemble | ✅ |
| H20b | FindDomains_CloselyOverlappingTandemSh3_OptOutEmitsSingleFusedEnvelope | ✅ |
| H20c | FindDomains_OverlappingTandemSh3_EnvelopeCoordinatesMatchHmmsearchAcrossTrims | ✅ |
| H20d | FindDomains_Ensemble_IsDeterministicAcrossRuns | ✅ |
| H20e | FindDomains_WellSeparatedTandemSh3_StillSplitsByFlankBound | ✅ |
| H20f | FindDomains_SingleDomain_StillYieldsOneEnvelope_UnderEnsembleDefault | ✅ |

### E.6 Residual

Envelope **count** and **coordinates** reproduce `hmmsearch` exactly (scores within ~0.06 bits). The
only residual is **exact RNG-bit parity** of the per-sample trace ensemble: HMMER samples in float32
with the Easel LCG (ported verbatim here), but this engine computes the Forward matrix in float64; the
200-sample consensus is robust to that, so coordinates match but per-sample trace-by-trace bit parity
is not asserted. Full unfiltered suite green. Work Queue Remaining = 0.
