using System;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// PROTMOTIF-FIND-001: Protein Motif Search — canonical test file.
/// Validates FindMotifByPattern and FindCommonMotifs against PROSITE patterns.
/// Evidence: PROSITE Database (https://prosite.expasy.org/), Hulo et al. (2007).
/// </summary>
[TestFixture]
public class ProteinMotifFinder_MotifSearch_Tests
{
    #region PROSITE Pattern Constants — verified against official PROSITE entries

    // PROSITE PS00001: N-glycosylation site — https://prosite.expasy.org/PS00001
    private const string Ps00001Pattern = "N-{P}-[ST]-{P}";
    private const string Ps00001Regex = @"N[^P][ST][^P]";

    // PROSITE PS00005: Protein kinase C phosphorylation site — https://prosite.expasy.org/PS00005
    private const string Ps00005Pattern = "[ST]-x-[RK]";
    private const string Ps00005Regex = @"[ST].[RK]";

    // PROSITE PS00006: Casein kinase II phosphorylation site — https://prosite.expasy.org/PS00006
    private const string Ps00006Pattern = "[ST]-x(2)-[DE]";
    private const string Ps00006Regex = @"[ST].{2}[DE]";

    // PROSITE PS00004: cAMP/cGMP-dependent phosphorylation site — https://prosite.expasy.org/PS00004
    private const string Ps00004Pattern = "[RK](2)-x-[ST]";
    private const string Ps00004Regex = @"[RK]{2}.[ST]";

    // PROSITE PS00007: Tyrosine kinase phosphorylation site — https://prosite.expasy.org/PS00007
    private const string Ps00007Pattern = "[RK]-x(2)-[DE]-x(3)-Y";
    private const string Ps00007Regex = @"[RK].{2}[DE].{3}Y";

    // PROSITE PS00008: N-myristoylation site — https://prosite.expasy.org/PS00008
    private const string Ps00008Pattern = "G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}";
    private const string Ps00008Regex = @"G[^EDRKHPFYW].{2}[STAGCN][^P]";

    // PROSITE PS00009: Amidation site — https://prosite.expasy.org/PS00009
    private const string Ps00009Pattern = "x-G-[RK]-[RK]";
    private const string Ps00009Regex = @".G[RK][RK]";

    // PROSITE PS00016: Cell attachment sequence (RGD) — https://prosite.expasy.org/PS00016
    private const string Ps00016Pattern = "R-G-D";
    private const string Ps00016Regex = @"RGD";

    // PROSITE PS00017: ATP/GTP-binding site motif A (P-loop) — https://prosite.expasy.org/PS00017
    private const string Ps00017Pattern = "[AG]-x(4)-G-K-[ST]";
    private const string Ps00017Regex = @"[AG].{4}GK[ST]";

    // PROSITE PS00018: EF-hand calcium-binding domain — https://prosite.expasy.org/PS00018
    private const string Ps00018Pattern = "D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]";
    private const string Ps00018Regex = @"D[^W][DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE][LIVMFYW]";

    // PROSITE PS00028: Zinc finger C2H2 type — https://prosite.expasy.org/PS00028
    private const string Ps00028Pattern = "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H";
    private const string Ps00028Regex = @"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H";

    // PROSITE PS00029: Leucine zipper pattern — https://prosite.expasy.org/PS00029
    private const string Ps00029Pattern = "L-x(6)-L-x(6)-L-x(6)-L";
    private const string Ps00029Regex = @"L.{6}L.{6}L.{6}L";

    #endregion

    #region M1: RGD Simple Pattern — PROSITE PS00016

    [Test]
    public void FindMotifByPattern_RGD_FindsThreeMatches()
    {
        // R-G-D (PS00016 cell attachment sequence) appears 3 times
        const string protein = "MRGDKLARGDPMRGD";
        var matches = FindMotifByPattern(protein, Ps00016Regex, "RGD", "PS00016").ToList();

        Assert.That(matches, Has.Count.EqualTo(3),
            "PS00016 RGD pattern should find exactly 3 occurrences in MRGDKLARGDPMRGD");
        Assert.That(matches.All(m => m.Sequence == "RGD"), Is.True,
            "All matches should have sequence 'RGD'");
    }

    #endregion

    #region M2: RGD Position Verification — Regex match semantics

    [Test]
    public void FindMotifByPattern_RGD_ReturnsCorrectPositions()
    {
        const string protein = "MRGDKLARGDPMRGD";
        var matches = FindMotifByPattern(protein, Ps00016Regex, "RGD", "PS00016").ToList();

        Assert.Multiple(() =>
        {
            // Position 0: M, Position 1-3: RGD
            Assert.That(matches[0].Start, Is.EqualTo(1),
                "First RGD starts at index 1");
            Assert.That(matches[0].End, Is.EqualTo(3),
                "First RGD ends at index 3");

            // Position 7-9: RGD (after KLAL)
            Assert.That(matches[1].Start, Is.EqualTo(7),
                "Second RGD starts at index 7");
            Assert.That(matches[1].End, Is.EqualTo(9),
                "Second RGD ends at index 9");

            // Position 12-14: RGD
            Assert.That(matches[2].Start, Is.EqualTo(12),
                "Third RGD starts at index 12");
            Assert.That(matches[2].End, Is.EqualTo(14),
                "Third RGD ends at index 14");
        });
    }

    #endregion

    #region M3: N-glycosylation Valid Site — PROSITE PS00001

    [Test]
    public void FindCommonMotifs_NGlycosylation_FindsValidSite()
    {
        // PS00001: N-{P}-[ST]-{P} — N followed by non-P, then S/T, then non-P
        const string protein = "AAAANFTAAAA";
        var motifs = FindCommonMotifs(protein).ToList();

        var glycoMatches = motifs.Where(m =>
            m.MotifName == "ASN_GLYCOSYLATION").ToList();

        Assert.That(glycoMatches, Has.Count.EqualTo(1),
            "Should find exactly one N-glycosylation site (NFTA) per PS00001");
        Assert.That(glycoMatches[0].Start, Is.EqualTo(4),
            "N-glycosylation match should start at position 4 (the N)");
        Assert.That(glycoMatches[0].Sequence, Does.StartWith("N"),
            "Matched sequence should start with N (asparagine)");
    }

    #endregion

    #region M4: N-glycosylation Proline Exclusion — PROSITE PS00001

    [Test]
    public void FindCommonMotifs_NGlycosylation_ExcludesProline()
    {
        // PS00001 uses {P} — proline is excluded at positions 2 and 4
        const string protein = "AAAANPSAAAAANPTAAA";
        var motifs = FindCommonMotifs(protein).ToList();

        var glycoMatches = motifs.Where(m =>
            m.MotifName == "ASN_GLYCOSYLATION").ToList();

        Assert.That(glycoMatches, Has.Count.EqualTo(0),
            "N-P-S and N-P-T should NOT match PS00001 because proline is at excluded position 2");
    }

    #endregion

    #region M5: PKC Phosphorylation — PROSITE PS00005

    [Test]
    public void FindCommonMotifs_PKC_FindsPhosphorylation()
    {
        // PS00005: [ST]-x-[RK] — S or T, any, then R or K
        const string protein = "AAAAASARKAAA";
        var motifs = FindCommonMotifs(protein).ToList();

        var pkcMatches = motifs.Where(m =>
            m.MotifName == "PKC_PHOSPHO_SITE").ToList();

        Assert.That(pkcMatches, Has.Count.GreaterThanOrEqualTo(1),
            "Should find at least one PKC phosphorylation site (SAR) per PS00005");
        Assert.That(pkcMatches.Any(m => m.Sequence.Length == 3),
            Is.True, "PKC match should be 3 residues long: [ST]-x-[RK]");
    }

    #endregion

    #region M6: CK2 Phosphorylation — PROSITE PS00006

    [Test]
    public void FindCommonMotifs_CK2_FindsPhosphorylation()
    {
        // PS00006: [ST]-x(2)-[DE] — S or T, two any, then D or E
        const string protein = "AAAASAAEASDEDAAA";
        var motifs = FindCommonMotifs(protein).ToList();

        var ck2Matches = motifs.Where(m =>
            m.MotifName == "CK2_PHOSPHO_SITE").ToList();

        Assert.That(ck2Matches, Has.Count.GreaterThanOrEqualTo(1),
            "Should find at least one CK2 phosphorylation site per PS00006");
        Assert.That(ck2Matches.Any(m => m.Sequence.Length == 4),
            Is.True, "CK2 match should be 4 residues long: [ST]-x(2)-[DE]");
    }

    #endregion

    #region M7: P-loop Detection — PROSITE PS00017

    [Test]
    public void FindCommonMotifs_PLoop_FindsSite()
    {
        // PS00017: [AG]-x(4)-G-K-[ST] — A or G, four any, then GK[ST]
        const string protein = "AAAAAGXXXXGKSAAAA";
        var motifs = FindCommonMotifs(protein).ToList();

        var pLoopMatches = motifs.Where(m =>
            m.MotifName == "ATP_GTP_A").ToList();

        Assert.That(pLoopMatches, Has.Count.GreaterThanOrEqualTo(1),
            "Should find at least one P-loop (ATP/GTP-binding) site per PS00017");
        Assert.That(pLoopMatches[0].Sequence.Length, Is.EqualTo(8),
            "P-loop match should be 8 residues long: [AG]-x(4)-G-K-[ST]");
    }

    #endregion

    #region M8: Empty Sequence — Trivial

    [Test]
    public void FindMotifByPattern_EmptySequence_ReturnsEmpty()
    {
        var matches = FindMotifByPattern("", Ps00016Regex, "RGD").ToList();

        Assert.That(matches, Is.Empty,
            "Empty sequence should produce no matches");
    }

    #endregion

    #region M9: Null Sequence — Trivial

    [Test]
    public void FindMotifByPattern_NullSequence_ReturnsEmpty()
    {
        var matches = FindMotifByPattern(null!, Ps00016Regex, "RGD").ToList();

        Assert.That(matches, Is.Empty,
            "Null sequence should produce no matches");
    }

    [Test]
    public void FindCommonMotifs_NullSequence_ReturnsEmpty()
    {
        var matches = FindCommonMotifs(null!).ToList();

        Assert.That(matches, Is.Empty,
            "FindCommonMotifs with null should produce no matches");
    }

    #endregion

    #region M10: Invalid Regex — Robustness

    [Test]
    public void FindMotifByPattern_InvalidRegex_ReturnsEmpty()
    {
        var matches = FindMotifByPattern("AAAA", "[[[invalid", "Bad").ToList();

        Assert.That(matches, Is.Empty,
            "Invalid regex pattern should return empty without throwing");
    }

    #endregion

    #region M11: Case Insensitivity — PROSITE convention

    [Test]
    public void FindMotifByPattern_CaseInsensitive_SameResults()
    {
        const string proteinUpper = "AAARGDAAA";
        const string proteinLower = "aaargdaaa";
        const string proteinMixed = "AaaRgDaAa";

        var matchesUpper = FindMotifByPattern(proteinUpper, Ps00016Regex, "RGD").ToList();
        var matchesLower = FindMotifByPattern(proteinLower, Ps00016Regex, "RGD").ToList();
        var matchesMixed = FindMotifByPattern(proteinMixed, Ps00016Regex, "RGD").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matchesUpper, Has.Count.EqualTo(1),
                "Uppercase input should find one RGD match");
            Assert.That(matchesLower, Has.Count.EqualTo(1),
                "Lowercase input should find one RGD match");
            Assert.That(matchesMixed, Has.Count.EqualTo(1),
                "Mixed-case input should find one RGD match");
            Assert.That(matchesUpper[0].Start, Is.EqualTo(matchesLower[0].Start),
                "Match positions should be identical regardless of case");
        });
    }

    #endregion

    #region M12: PS00007 Correct Pattern — PROSITE PS00007

    [Test]
    public void CommonMotifs_PS00007_MatchesOfficialPrositeDefinition()
    {
        // Official PROSITE PS00007: [RK]-x(2)-[DE]-x(3)-Y (fixed repeats, not ranges)
        var ps00007 = CommonMotifs["PS00007"];

        Assert.Multiple(() =>
        {
            Assert.That(ps00007.Pattern, Is.EqualTo(Ps00007Pattern),
                "PS00007 PROSITE pattern must match official definition: [RK]-x(2)-[DE]-x(3)-Y");
            Assert.That(ps00007.RegexPattern, Is.EqualTo(Ps00007Regex),
                "PS00007 regex must be [RK].{2}[DE].{3}Y (exact repeats, not ranges)");
        });
    }

    [Test]
    public void FindMotifByPattern_PS00007_MatchesCorrectLength()
    {
        // [RK]-x(2)-[DE]-x(3)-Y = 1+2+1+3+1 = 8 residues
        const int expectedMatchLength = 8;
        const string protein = "AAARAAEDDDYAAAA";
        var matches = FindMotifByPattern(protein, Ps00007Regex, "TYR_PHOSPHO_SITE").ToList();

        Assert.That(matches, Has.Count.EqualTo(1),
            "Should find exactly one PS00007 match in RAAEDDDY");
        Assert.That(matches[0].Sequence.Length, Is.EqualTo(expectedMatchLength),
            "PS00007 match must be exactly 8 residues: [RK](1)+x(2)+[DE](1)+x(3)+Y(1)");
    }

    #endregion

    #region M13: PS00018 Correct Pattern — PROSITE PS00018

    [Test]
    public void CommonMotifs_PS00018_MatchesOfficialPrositeDefinition()
    {
        // Official PROSITE PS00018: D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]
        var ps00018 = CommonMotifs["PS00018"];

        Assert.Multiple(() =>
        {
            Assert.That(ps00018.Pattern, Is.EqualTo(Ps00018Pattern),
                "PS00018 PROSITE pattern must match official definition with {W} at position 2 and trailing [LIVMFYW]");
            Assert.That(ps00018.RegexPattern, Is.EqualTo(Ps00018Regex),
                "PS00018 regex must include [^W] at position 2 and [LIVMFYW] at end");
        });
    }

    #endregion

    #region M14: All PROSITE Patterns Correct — PROSITE entries

    [Test]
    public void CommonMotifs_AllPrositePatterns_MatchOfficialDefinitions()
    {
        // Verify every PROSITE-sourced entry against its official pattern
        var expectedPatterns = new Dictionary<string, (string PrositePattern, string Regex)>
        {
            ["PS00001"] = (Ps00001Pattern, Ps00001Regex),
            ["PS00005"] = (Ps00005Pattern, Ps00005Regex),
            ["PS00006"] = (Ps00006Pattern, Ps00006Regex),
            ["PS00004"] = (Ps00004Pattern, Ps00004Regex),
            ["PS00007"] = (Ps00007Pattern, Ps00007Regex),
            ["PS00008"] = (Ps00008Pattern, Ps00008Regex),
            ["PS00009"] = (Ps00009Pattern, Ps00009Regex),
            ["PS00016"] = (Ps00016Pattern, Ps00016Regex),
            ["PS00017"] = (Ps00017Pattern, Ps00017Regex),
            ["PS00018"] = (Ps00018Pattern, Ps00018Regex),
            ["PS00028"] = (Ps00028Pattern, Ps00028Regex),
            ["PS00029"] = (Ps00029Pattern, Ps00029Regex),
        };

        Assert.Multiple(() =>
        {
            foreach (var (accession, expected) in expectedPatterns)
            {
                Assert.That(CommonMotifs.ContainsKey(accession), Is.True,
                    $"CommonMotifs should contain PROSITE entry {accession}");

                var motif = CommonMotifs[accession];
                Assert.That(motif.Pattern, Is.EqualTo(expected.PrositePattern),
                    $"{accession} PROSITE pattern mismatch");
                Assert.That(motif.RegexPattern, Is.EqualTo(expected.Regex),
                    $"{accession} regex pattern mismatch");
            }
        });
    }

    #endregion

    #region M15: All Regex Patterns Valid — Regex correctness

    [Test]
    public void CommonMotifs_AllRegexPatterns_CompileSuccessfully()
    {
        Assert.Multiple(() =>
        {
            foreach (var motif in CommonMotifs.Values)
            {
                Assert.DoesNotThrow(
                    () => new System.Text.RegularExpressions.Regex(motif.RegexPattern),
                    $"Regex for {motif.Accession} ({motif.Name}) should compile without error");
            }
        });
    }

    #endregion

    #region M16: Substring Consistency Invariant — INV-3

    [Test]
    public void FindMotifByPattern_MatchSequence_EqualsSubstring()
    {
        const string protein = "AAARGDAAARGDAAA";
        var matches = FindMotifByPattern(protein, Ps00016Regex, "RGD").ToList();

        Assert.Multiple(() =>
        {
            foreach (var match in matches)
            {
                string expected = protein.ToUpperInvariant()
                    .Substring(match.Start, match.End - match.Start + 1);
                Assert.That(match.Sequence, Is.EqualTo(expected),
                    $"Match at [{match.Start}..{match.End}]: Sequence should equal input substring");
            }
        });
    }

    [Test]
    public void FindCommonMotifs_AllMatches_SatisfySubstringInvariant()
    {
        const string protein = "ANFTASARKGXXXXGKSRGD";
        string upper = protein.ToUpperInvariant();
        var matches = FindCommonMotifs(protein).ToList();

        Assert.Multiple(() =>
        {
            foreach (var match in matches)
            {
                Assert.That(match.Start, Is.GreaterThanOrEqualTo(0),
                    $"Match Start should be >= 0 for {match.MotifName}");
                Assert.That(match.End, Is.LessThan(upper.Length),
                    $"Match End should be < sequence length for {match.MotifName}");

                string expected = upper.Substring(match.Start, match.End - match.Start + 1);
                Assert.That(match.Sequence, Is.EqualTo(expected),
                    $"Match Sequence should equal substring for {match.MotifName} at [{match.Start}..{match.End}]");
            }
        });
    }

    #endregion

    #region M17: Score Information Content — Schneider & Stephens 1990

    [Test]
    public void FindMotifByPattern_Score_IsInformationContent()
    {
        // RGD pattern: 3 fixed positions → IC = 3 × log2(20) ≈ 12.97 bits
        const string protein = "AAARGDAAA";
        var matches = FindMotifByPattern(protein, Ps00016Regex, "RGD").ToList();

        Assert.That(matches, Has.Count.EqualTo(1),
            "Should find exactly one RGD match");

        double expectedIC = 3 * Math.Log2(20); // ≈ 12.97
        Assert.Multiple(() =>
        {
            Assert.That(matches[0].Score, Is.EqualTo(expectedIC).Within(0.01),
                "Score should equal total information content (3 fixed positions × log2(20))");
            Assert.That(matches[0].EValue, Is.GreaterThan(0),
                "E-value should be positive for any valid match");
        });
    }

    [Test]
    public void FindMotifByPattern_Score_AccountsForCharacterClasses()
    {
        // PKC pattern [ST].[RK]: IC = log2(20/2) + 0 + log2(20/2) ≈ 6.64 bits
        const string protein = "AAAASAR";
        var matches = FindMotifByPattern(protein, Ps00005Regex, "PKC").ToList();

        Assert.That(matches, Has.Count.GreaterThanOrEqualTo(1));

        double expectedIC = 2 * Math.Log2(20.0 / 2); // ≈ 6.64
        Assert.That(matches[0].Score, Is.EqualTo(expectedIC).Within(0.01),
            "Score for [ST].[RK]: two 2-choice positions (IC=3.32 each) + one any (IC=0)");
    }

    [Test]
    public void FindMotifByPattern_EValue_UsesProperProbability()
    {
        // For RGD in "AAARGDAAA" (length 9, motif length 3):
        // E = (9 - 3 + 1) × 2^(-IC) = 7 × 2^(-12.97) ≈ 7/8000
        const string protein = "AAARGDAAA";
        var matches = FindMotifByPattern(protein, Ps00016Regex, "RGD").ToList();

        Assert.That(matches, Has.Count.EqualTo(1));

        double ic = 3 * Math.Log2(20);
        double expectedEValue = 7 * Math.Pow(2.0, -ic); // 7/8000 ≈ 0.000875
        Assert.That(matches[0].EValue, Is.EqualTo(expectedEValue).Within(1e-6),
            "E-value = (N-L+1) × 2^(-IC) per Altschul et al. 1990");
    }

    #endregion

    #region S1: Multiple Patterns Same Sequence — Integration

    [Test]
    public void FindCommonMotifs_MultiplePatterns_ReturnsMultipleMotifTypes()
    {
        // Sequence containing N-glycosylation (NFTA), PKC (SAR[K]), RGD, and P-loop
        const string protein = "NFTASARKRGDAGXXXXGKS";
        var motifs = FindCommonMotifs(protein).ToList();

        var motifNames = motifs.Select(m => m.MotifName).Distinct().ToList();

        Assert.That(motifNames.Count, Is.GreaterThanOrEqualTo(2),
            "Multiple motif types should be detected in a sequence with diverse sites");
    }

    #endregion

    #region S2: No Match Returns Empty — Edge case

    [Test]
    public void FindMotifByPattern_NoMatch_ReturnsEmpty()
    {
        // Zinc finger pattern requires C...C...H...H — not present in all-A sequence
        const string protein = "AAAAAAAAAAAAAAAA";
        var matches = FindMotifByPattern(protein, Ps00028Regex, "Zinc Finger").ToList();

        Assert.That(matches, Is.Empty,
            "Sequence without pattern elements should return no matches");
    }

    [Test]
    public void FindMotifByPattern_EmptyPattern_ReturnsEmpty()
    {
        var matches = FindMotifByPattern("AAAA", "", "Empty").ToList();

        Assert.That(matches, Is.Empty,
            "Empty pattern should return no matches");
    }

    #endregion

    #region S3: MotifMatch Fields — Consistency

    [Test]
    public void FindMotifByPattern_MatchFields_ArePopulated()
    {
        const string protein = "AAARGDAAA";
        const string motifName = "TestMotif";
        const string patternId = "TEST001";
        var matches = FindMotifByPattern(protein, Ps00016Regex, motifName, patternId).ToList();

        Assert.That(matches, Has.Count.EqualTo(1));

        Assert.Multiple(() =>
        {
            Assert.That(matches[0].MotifName, Is.EqualTo(motifName),
                "MotifName should match the provided name");
            Assert.That(matches[0].Pattern, Is.EqualTo(patternId),
                "Pattern field should match the provided pattern ID");
            Assert.That(matches[0].Sequence, Is.EqualTo("RGD"),
                "Sequence should be the matched text");
        });
    }

    #endregion

    #region S4: Overlapping Match Discovery — ScanProsite consistency

    [Test]
    public void FindMotifByPattern_OverlappingMatches_AllDiscovered()
    {
        // Pattern P-x-x-P applied to PPPPP should find overlapping matches
        // at positions 0 (PPPP) and 1 (PPPP)
        const string protein = "PPPPP";
        var matches = FindMotifByPattern(protein, @"P..P", "Test").ToList();

        Assert.That(matches, Has.Count.EqualTo(2),
            "Overlapping pattern P..P in PPPPP should find 2 matches: pos 0 and pos 1");
        Assert.Multiple(() =>
        {
            Assert.That(matches[0].Start, Is.EqualTo(0), "First overlapping match at position 0");
            Assert.That(matches[1].Start, Is.EqualTo(1), "Second overlapping match at position 1");
        });
    }

    [Test]
    public void FindMotifByPattern_NonOverlapping_SameAsOverlapping()
    {
        // For non-overlapping patterns, results should be identical
        const string protein = "MRGDKLARGDPMRGD";
        var matches = FindMotifByPattern(protein, Ps00016Regex, "RGD").ToList();

        Assert.That(matches, Has.Count.EqualTo(3),
            "Non-overlapping RGD matches should still find all 3 occurrences");
    }

    #endregion

    #region S5: Non-PROSITE Patterns — Literature verification

    [Test]
    public void CommonMotifs_NLS1_MatchesChelskysConsensus()
    {
        // NLS monopartite consensus: K-K/R-X-K/R (Dingwall & Laskey 1991)
        var nls = CommonMotifs["NLS1"];

        Assert.Multiple(() =>
        {
            Assert.That(nls.Pattern, Is.EqualTo("[KR]-[KR]-x-[KR]"),
                "NLS1 pattern should match Chelsky consensus K-K/R-X-K/R");
            Assert.That(nls.RegexPattern, Is.EqualTo(@"[KR][KR].[KR]"),
                "NLS1 regex should match [KR][KR].[KR]");
        });

        // Verify it matches a canonical NLS: KKKR
        var matches = FindMotifByPattern("AAAKKKRAAA", nls.RegexPattern, nls.Name).ToList();
        Assert.That(matches, Has.Count.EqualTo(1), "Should find KKKR as NLS");
        Assert.That(matches[0].Sequence, Is.EqualTo("KKKR"));
    }

    [Test]
    public void CommonMotifs_NES1_MatchesLaCourConsensus()
    {
        // NES consensus: Φ1-x(2,3)-Φ2-x(2,3)-Φ3-x-Φ4 (la Cour et al. 2004)
        var nes = CommonMotifs["NES1"];

        Assert.That(nes.Pattern, Is.EqualTo("L-x(2,3)-[LIVFM]-x(2,3)-L-x-[LI]"),
            "NES1 pattern should follow la Cour et al. 2004 NES consensus");

        // Verify with canonical NES: LAAILAALALI
        var matches = FindMotifByPattern("AAALAALAALALI", nes.RegexPattern, nes.Name).ToList();
        Assert.That(matches, Has.Count.GreaterThanOrEqualTo(1), "Should find NES in leucine-rich sequence");
    }

    [Test]
    public void CommonMotifs_SIM1_MatchesHeckerConsensus()
    {
        // SIM type b: Ψ-x-Ψ-Ψ where Ψ = V/I/L (Hecker et al. 2006, JBC 281:16117-27)
        var sim = CommonMotifs["SIM1"];

        Assert.Multiple(() =>
        {
            Assert.That(sim.Pattern, Is.EqualTo("[VIL]-x-[VIL]-[VIL]"),
                "SIM1 pattern should match Hecker type b consensus");
            Assert.That(sim.RegexPattern, Is.EqualTo(@"[VIL].[VIL][VIL]"),
                "SIM1 regex should be [VIL].[VIL][VIL]");
        });

        // Verify: VAVV is a valid SIM
        var matches = FindMotifByPattern("AAAVAVVAAA", sim.RegexPattern, sim.Name).ToList();
        Assert.That(matches, Has.Count.EqualTo(1), "Should find VAVV as SIM");
        Assert.That(matches[0].Sequence, Is.EqualTo("VAVV"));
    }

    [Test]
    public void CommonMotifs_WW1_MatchesChenSudolConsensus()
    {
        // PPxY motif (Chen & Sudol 1995, PNAS 92:7819-23)
        var ww = CommonMotifs["WW1"];

        Assert.Multiple(() =>
        {
            Assert.That(ww.Pattern, Is.EqualTo("P-P-x-Y"),
                "WW1 pattern should match Chen & Sudol PPxY consensus");
            Assert.That(ww.RegexPattern, Is.EqualTo(@"PP.Y"),
                "WW1 regex should be PP.Y");
        });

        // Verify: PPAY is a valid PY motif
        var matches = FindMotifByPattern("AAAPPAYAAA", ww.RegexPattern, ww.Name).ToList();
        Assert.That(matches, Has.Count.EqualTo(1), "Should find PPAY as PY motif");
        Assert.That(matches[0].Sequence, Is.EqualTo("PPAY"));
    }

    [Test]
    public void CommonMotifs_SH3_1_MatchesMayerClassIConsensus()
    {
        // SH3 Class I: +xxPxxP where + = R/K (Mayer 2001, J Cell Sci 114:1253-63)
        var sh3 = CommonMotifs["SH3_1"];

        Assert.Multiple(() =>
        {
            Assert.That(sh3.Pattern, Is.EqualTo("[RK]-x(2)-P-x(2)-P"),
                "SH3_1 pattern should match Mayer Class I consensus +xxPxxP");
            Assert.That(sh3.RegexPattern, Is.EqualTo(@"[RK].{2}P.{2}P"),
                "SH3_1 regex should be [RK].{2}P.{2}P");
        });

        // Verify: RAAPAAP is a valid Class I SH3 binding motif
        var matches = FindMotifByPattern("AARAAPAAPAAA", sh3.RegexPattern, sh3.Name).ToList();
        Assert.That(matches, Has.Count.EqualTo(1), "Should find RAAPAA P as SH3 class I");
        Assert.That(matches[0].Sequence, Is.EqualTo("RAAPAA P".Replace(" ", "")));
    }

    [Test]
    public void CommonMotifs_AllNonProsite_HaveCorrectPatterns()
    {
        // Comprehensive verification: all 5 non-PROSITE entries present with source-verified patterns
        var expected = new Dictionary<string, (string Pattern, string Regex)>
        {
            ["NLS1"] = ("[KR]-[KR]-x-[KR]", @"[KR][KR].[KR]"),
            ["NES1"] = ("L-x(2,3)-[LIVFM]-x(2,3)-L-x-[LI]", @"L.{2,3}[LIVFM].{2,3}L.[LI]"),
            ["SIM1"] = ("[VIL]-x-[VIL]-[VIL]", @"[VIL].[VIL][VIL]"),
            ["WW1"] = ("P-P-x-Y", @"PP.Y"),
            ["SH3_1"] = ("[RK]-x(2)-P-x(2)-P", @"[RK].{2}P.{2}P"),
        };

        Assert.Multiple(() =>
        {
            foreach (var (key, (pattern, regex)) in expected)
            {
                Assert.That(CommonMotifs.ContainsKey(key), Is.True,
                    $"CommonMotifs should contain non-PROSITE entry {key}");
                Assert.That(CommonMotifs[key].Pattern, Is.EqualTo(pattern),
                    $"{key} PROSITE-style pattern mismatch");
                Assert.That(CommonMotifs[key].RegexPattern, Is.EqualTo(regex),
                    $"{key} regex pattern mismatch");
            }
        });
    }

    #endregion
}
