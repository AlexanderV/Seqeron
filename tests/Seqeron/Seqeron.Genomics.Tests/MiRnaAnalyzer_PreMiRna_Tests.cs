using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// MIRNA-PRECURSOR-001: Pre-miRNA Hairpin Detection
/// Canonical test file for FindPreMiRnaHairpins (and indirectly AnalyzeHairpin).
/// Evidence: Bartel (2004), Ambros (2003), Krol (2004), Bartel (2009), miRBase.
/// </summary>
[TestFixture]
public class MiRnaAnalyzer_PreMiRna_Tests
{
    #region Test Data

    // ── Biologically plausible mixed-base hairpins ─────────────────────────
    // All test hairpins use mixed nucleotide composition (AUGC), matching
    // the diversity found in real pre-miRNAs. Mono-nucleotide sequences
    // (all-G stems, all-C stems) are avoided as they are biologically
    // implausible and can mask implementation bugs.

    // 25 bp Watson-Crick stem + 7 nt loop = 57 nt
    // maxStem = 57/2 - 5 = 23 → effective stem=23, effective loop=11
    private const string ValidHairpin57 =
        "GCAUAGCUAGCUAGCUAGCUAGCUA" +   // 25 nt 5' stem (mixed)
        "GAAAUUU" +                       // 7 nt loop
        "UAGCUAGCUAGCUAGCUAGCUAUGC";     // 25 nt 3' stem (reverse complement)

    // 20 bp Watson-Crick stem + 7 nt loop = 47 nt (below default min=55)
    private const string ShortHairpin47 =
        "GCAUAGCUAGCUAGCUAGCU" +          // 20 nt 5' stem
        "GAAAUUU" +                        // 7 nt loop
        "AGCUAGCUAGCUAGCUAUGC";           // 20 nt 3' stem (reverse complement)

    // 25 bp stem with 4 G:U wobble pairs (at pairing positions 0, 5, 9, 17) + 7 nt loop = 57 nt
    // Wobble pairs are common in real pre-miRNAs — Krol (2004)
    private const string WobbleHairpin57 =
        "GCAUAGCUAGCUAGCUAGCUAGCUA" +     // 25 nt 5' stem (same as ValidHairpin57)
        "GAAAUUU" +                        // 7 nt loop
        "UAGCUAGUUAGCUAGUUAGUUAUGU";      // 25 nt 3' stem (4× C→U = G:U wobble)

    // 15 bp stem + 25 nt loop = 55 nt — tests stem rejection (not n<55 hardcode)
    // AnalyzeHairpin rejects stems < 18 bp; this sequence is ≥55 nt total
    // so the rejection must come from stem length, not the n<55 early-exit.
    private const string ShortStemHairpin55 =
        "GCAUAGCUAGCUAGC" +               // 15 nt 5' stem
        "AUGCAUGCAUGCAUGCAUGCAUGCA" +      // 25 nt loop (mixed, non-self-complementary)
        "GCUAGCUAGCUAUGC";                 // 15 nt 3' stem (reverse complement)

    // No complementarity: alternating purines (A-G never pair with each other)
    private const string NoComplementarity =
        "AGAGAGAGAGAGAGAGAGAGAGAGA" +      // 25 nt
        "GAGAGAGAGAGAGAGAGAGAGAGAG" +      // 25 nt
        "AGAGAGAGAGAGAGAGAGAG";            // 20 nt = 70 total

    // DNA input with T instead of U (same stem structure as ValidHairpin57)
    private const string DnaHairpin57 =
        "GCATAGCTAGCTAGCTAGCTAGCTA" +     // 25 nt 5' stem (DNA: T for U)
        "GAAATTT" +                        // 7 nt loop
        "TAGCTAGCTAGCTAGCTAGCTATGC";       // 25 nt 3' stem

    // Large hairpin: 50 bp mixed stem + 10 nt loop = 110 nt
    private static readonly string LargeHairpin110 = BuildPerfectHairpin(50, 10);

    // Two hairpins separated by non-complementary spacer
    private static readonly string TwoHairpinSequence =
        ValidHairpin57 + "AGAGAGAGAGAGAGAGAGAG" + ValidHairpin57;

    // ── Real miRBase pre-miRNA sequences (known-limitation reference) ──────
    // These sequences are from miRBase v22 and represent real pre-miRNAs.
    // The simplified consecutive-pairing model CANNOT detect them because
    // real pre-miRNAs have internal mismatches and bulges that break
    // consecutive pairing from the ends.

    /// <summary>hsa-mir-21 (MI0000077) — 72 nt, miRBase v22</summary>
    private const string HsaMir21_MI0000077 =
        "UGUCGGGUAGCUUAUCAGACUGAUGUUGA" +
        "CUGUUGAAUCUCAUGGCAACACCAGUCGA" +
        "UGGGCUGUCUGACA";

    /// <summary>hsa-let-7a-1 (MI0000060) — 80 nt, miRBase v22</summary>
    private const string HsaLet7a1_MI0000060 =
        "UGGGAUGAGGUAGUAGGUUGUAUAGUUUU" +
        "AGGGUCACACCCACCACUGGGAGAUAACU" +
        "AUACAAUCUACUGUCUUUCCUA";

    #endregion

    #region M1: Null Input — Defensive coding

    [Test]
    public void FindPreMiRnaHairpins_NullInput_ReturnsEmpty()
    {
        var results = FindPreMiRnaHairpins(null!).ToList();

        Assert.That(results, Is.Empty,
            "Null input must return empty enumerable, not throw");
    }

    #endregion

    #region M2: Empty Input — Defensive coding

    [Test]
    public void FindPreMiRnaHairpins_EmptyInput_ReturnsEmpty()
    {
        var results = FindPreMiRnaHairpins("").ToList();

        Assert.That(results, Is.Empty,
            "Empty string must return empty enumerable");
    }

    #endregion

    #region M3: Short Sequence — Bartel (2004)

    [Test]
    public void FindPreMiRnaHairpins_ShortSequence_ReturnsEmpty()
    {
        // 47 nt hairpin is below default minHairpinLength=55
        var results = FindPreMiRnaHairpins(ShortHairpin47, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Sequence shorter than minHairpinLength (55) must return empty — " +
            "pre-miRNAs are ≥55 nt (Bartel 2004, miRBase)");
    }

    #endregion

    #region M4: Valid Hairpin Detection — Bartel (2004), Krol (2004)

    [Test]
    public void FindPreMiRnaHairpins_ValidHairpin_DetectsCandidate()
    {
        // 57 nt hairpin: 25 bp G-C stem + 7 nt loop
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1),
            "A 57 nt sequence with 25 bp perfect stem and 7 nt loop must be detected — " +
            "meets all structural criteria (stem ≥ 18 bp, loop 3–25 nt, length ≥ 55)");
    }

    #endregion

    #region M5: Position Correct — Array bounds

    [Test]
    public void FindPreMiRnaHairpins_Position_WithinInputBounds()
    {
        // ValidHairpin57 is exactly 57 nt. Scanning with minHairpinLength=55:
        // - i=0, len=57 → stem=23 (full hairpin) → accepted
        // - i=1, len=55 → stem=22 (sub-window) → accepted
        // Other windows rejected (pairing breaks at position 0/1).
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Has.Count.EqualTo(2),
            "57 nt input with min=55 must yield exactly 2 candidates: " +
            "full 57-nt window (stem=23) and 55-nt sub-window at offset 1 (stem=22)");

        // Primary result: full hairpin window
        Assert.Multiple(() =>
        {
            Assert.That(results[0].Start, Is.EqualTo(0),
                "Primary hairpin starts at position 0");
            Assert.That(results[0].End, Is.EqualTo(56),
                "Primary hairpin ends at position 56 (0-indexed, inclusive)");
            Assert.That(results[0].Sequence.Length, Is.EqualTo(57),
                "Primary hairpin is the full 57-nt input");
        });

        // Secondary result: sub-window starting at offset 1
        Assert.Multiple(() =>
        {
            Assert.That(results[1].Start, Is.EqualTo(1),
                "Sub-window hairpin starts at position 1");
            Assert.That(results[1].End, Is.EqualTo(55),
                "Sub-window hairpin ends at position 55");
            Assert.That(results[1].Sequence.Length, Is.EqualTo(55),
                "Sub-window is 55 nt");
        });
    }

    #endregion

    #region M6: Mature From 5' Arm — Bartel (2009)

    [Test]
    public void FindPreMiRnaHairpins_MatureSequence_ExtractedFrom5PrimeArm()
    {
        // For 57-nt hairpin with stem=23: mature = first min(22, 23) = 22 nt
        // Hand-derived: first 22 chars of ValidHairpin57
        const string expectedMature = "GCAUAGCUAGCUAGCUAGCUAG";

        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];

        Assert.Multiple(() =>
        {
            Assert.That(pre.MatureSequence, Is.EqualTo(expectedMature),
                "Mature must be the first 22 nt of 5' arm — " +
                "matureEnd = min(matureLength=22, stemLength=23) = 22 (Bartel 2009)");
            Assert.That(pre.Sequence.StartsWith(pre.MatureSequence, StringComparison.Ordinal),
                "Mature must come from 5' arm (beginning of hairpin sequence)");
        });
    }

    #endregion

    #region M7: Star From 3' Arm — Bartel (2009)

    [Test]
    public void FindPreMiRnaHairpins_StarSequence_ExtractedFrom3PrimeArm()
    {
        // For 57-nt hairpin: star = last 22 chars (positions 35-56)
        // Hand-derived from 3' stem of ValidHairpin57
        const string expectedStar = "CUAGCUAGCUAGCUAGCUAUGC";

        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];

        Assert.Multiple(() =>
        {
            Assert.That(pre.StarSequence, Is.EqualTo(expectedStar),
                "Star must be the last 22 nt from 3' arm — " +
                "starStart = n(57) - matureEnd(22) = 35 (Bartel 2009)");
            Assert.That(pre.StarSequence.Length, Is.EqualTo(pre.MatureSequence.Length),
                "Star and mature must have equal length (duplex symmetry)");
            Assert.That(pre.Sequence.EndsWith(pre.StarSequence, StringComparison.Ordinal),
                "Star must come from 3' arm (end of hairpin sequence)");
        });
    }

    #endregion

    #region M8: Dot-Bracket Structure — Standard RNA notation

    [Test]
    public void FindPreMiRnaHairpins_Structure_MatchesStemLoopPattern()
    {
        // For 57-nt hairpin with stem=23, loop=11:
        // Structure = 23×'(' + 11×'.' + 23×')'
        const string expectedStructure =
            "(((((((((((((((((((((((...........)))))))))))))))))))))))";

        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];

        Assert.Multiple(() =>
        {
            Assert.That(pre.Structure, Is.EqualTo(expectedStructure),
                "Structure must be exactly 23×'(' + 11×'.' + 23×')' — " +
                "stem=23 (maxStem=57/2-5=23), loop=57-2×23=11");
            Assert.That(pre.Structure.Length, Is.EqualTo(57),
                "Structure length must equal sequence length (INV-5)");
        });
    }

    #endregion

    #region M9: Dot-Bracket Balanced — Notation definition

    [Test]
    public void FindPreMiRnaHairpins_Structure_BalancedParentheses()
    {
        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        foreach (var pre in results)
        {
            int openCount = pre.Structure.Count(c => c == '(');
            int closeCount = pre.Structure.Count(c => c == ')');
            int dotCount = pre.Structure.Count(c => c == '.');

            Assert.Multiple(() =>
            {
                Assert.That(openCount, Is.EqualTo(closeCount),
                    "Count of '(' must equal count of ')' — balanced base pairs (INV-7)");
                Assert.That(dotCount, Is.GreaterThanOrEqualTo(3),
                    "Loop must have ≥ 3 unpaired positions (Bartel 2004)");
                Assert.That(openCount + closeCount + dotCount, Is.EqualTo(pre.Structure.Length),
                    "Structure must contain only '(', ')', and '.' characters (INV-6)");
            });
        }
    }

    #endregion

    #region M10: Free Energy Negative — Thermodynamic stability

    [Test]
    public void FindPreMiRnaHairpins_FreeEnergy_IsNegative()
    {
        // Hand-calculated from Turner 2004 NNDB parameters for ValidHairpin57:
        //   stem=23, loop=11
        //   Stacking (22 pairs): 6×GC/CG(-3.42) + 1×CA/GU(-2.11) + 1×AU/UA(-1.10)
        //                      + 5×UA/AU(-1.33) + 5×AG/UC(-2.08) + 4×CU/GA(-2.08) = -49.10
        //   Loop init(11):    +6.60 (NNDB hairpin loop table)
        //   Terminal mismatch(CUAG, closing C-G): -1.00 (NNDB tm-parameters)
        //   AU/GU penalties:  0.00 (both outer G-C and closing C-G are WC, not AU/GU)
        //   Total: -43.50 kcal/mol
        const double expectedEnergy = -43.50;

        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        Assert.That(results[0].FreeEnergy, Is.EqualTo(expectedEnergy).Within(0.01),
            "Free energy must match hand-calculated Turner 2004 NNDB value (-43.50 kcal/mol) — " +
            "sum of 22 stacking pairs + loop init(11) + terminal mismatch(CUAG)");

        foreach (var pre in results)
        {
            Assert.That(pre.FreeEnergy, Is.LessThan(0),
                "All candidates must have negative free energy (thermodynamic stability)");
        }
    }

    #endregion

    #region M11: Free Energy Ordering — Turner (2004) principles

    [Test]
    public void FindPreMiRnaHairpins_FreeEnergy_LongerStemMoreNegative()
    {
        // Hairpin A: 57 nt, stem=23, loop=11
        //   Hand-calculated Turner energy: -43.50 kcal/mol (see M10 for breakdown)
        //
        // Hairpin B: 55 nt, stem=20, loop=15
        //   Stacking (19 pairs): 5×GC/CG(-3.42) + 1×CA/GU(-2.11) + 1×AU/UA(-1.10)
        //                       + 4×UA/AU(-1.33) + 4×AG/UC(-2.08) + 4×CU/GA(-2.08) = -42.27
        //   Loop init(15):    +6.90
        //   Terminal mismatch(UAGA, closing U-A): -1.10
        //   AU/GU penalty (closing U-A): +0.45
        //   Total: -36.02 kcal/mol
        const double expectedEnergyLong = -43.50;
        const double expectedEnergyShort = -36.02;

        string stem20Hairpin55 =
            "GCAUAGCUAGCUAGCUAGCU" +     // 20 nt 5' stem
            "AUGCAUGCAUGCAUG" +           // 15 nt loop
            "AGCUAGCUAGCUAGCUAUGC";       // 20 nt 3' stem (reverse complement)

        var resultsLong = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();
        var resultsShort = FindPreMiRnaHairpins(stem20Hairpin55, minHairpinLength: 55).ToList();

        Assert.That(resultsLong, Is.Not.Empty, "23 bp effective stem hairpin (57 nt) must be detected");
        Assert.That(resultsShort, Is.Not.Empty, "20 bp stem hairpin (55 nt) must be detected");

        Assert.Multiple(() =>
        {
            Assert.That(resultsLong[0].FreeEnergy, Is.EqualTo(expectedEnergyLong).Within(0.01),
                "Stem-23 energy must match hand-calculated Turner 2004 value");
            Assert.That(resultsShort[0].FreeEnergy, Is.EqualTo(expectedEnergyShort).Within(0.01),
                "Stem-20 energy must match hand-calculated Turner 2004 value");
            Assert.That(resultsLong[0].FreeEnergy, Is.LessThan(resultsShort[0].FreeEnergy),
                "Longer effective stem (23 bp, -43.50) must be more negative than shorter stem (20 bp, -36.02) — " +
                "reflects greater thermodynamic stability (Turner 2004)");
        });
    }

    #endregion

    #region M12: Stem Too Short Rejected — Krol (2004)

    [Test]
    public void FindPreMiRnaHairpins_StemTooShort_ReturnsEmpty()
    {
        // 15 bp stem + 25 nt loop = 55 nt total (passes n<55 hardcode)
        // Rejection must come from stemLength (15) < 18, not from length filter
        var results = FindPreMiRnaHairpins(ShortStemHairpin55, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Sequence with only 15 bp stem must be rejected — " +
            "pre-miRNA requires ≥ 18 bp stem (Krol 2004). " +
            "Note: sequence is 55 nt, so this tests stem rejection, not the n<55 early-exit");
    }

    #endregion

    #region M13: Loop Size Boundary Rejected — Bartel (2004)

    [Test]
    public void FindPreMiRnaHairpins_LoopTooLarge_Rejected()
    {
        // Note: Loop-too-small (< 3) is structurally unreachable because
        // maxStem = n/2 - 5 ensures minimum loop ≥ 10 nt for any valid input.
        // Instead test the reachable arm: loop > 25 → rejected.
        // 18 bp mixed stem + 30 nt loop = 66 nt. Stem=18, loop=30 > 25
        // Loop uses all-A to prevent self-complementary pairing within the loop
        string largeLoop =
            "GCAUAGCUAGCUAGCUAG" +             // 18 nt 5' stem (mixed)
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + // 30 nt loop (A-A never pairs)
            "CUAGCUAGCUAGCUAUGC";              // 18 nt 3' stem (reverse complement)
        var results = FindPreMiRnaHairpins(largeLoop, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Candidate with 30 nt loop must be rejected — " +
            "loop must be ≤ 25 nt (Bartel 2004, implementation threshold)");
    }

    #endregion

    #region M14: DNA T→U Conversion — RNA biology standard

    [Test]
    public void FindPreMiRnaHairpins_DnaInput_ConvertsToRna()
    {
        // DNA input: G's + T's (loop) + C's
        var results = FindPreMiRnaHairpins(DnaHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty,
            "DNA input with T must be handled (T→U conversion)");

        var pre = results[0];
        Assert.That(pre.Sequence, Does.Not.Contain("T"),
            "Output sequence must use U not T (RNA convention) — INV-10");
    }

    #endregion

    #region M15: G:U Wobble Pairs in Stem — Krol (2004)

    [Test]
    public void FindPreMiRnaHairpins_GU_WobblePairs_AcceptedInStem()
    {
        // 25 G's + 7 A's loop + 25 U's = 57 nt
        // G-U wobble pairs should be accepted by CanPair
        var results = FindPreMiRnaHairpins(WobbleHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(1),
            "G-U wobble pairs must be accepted as valid base pairs in stem — " +
            "common in real pre-miRNAs (Krol 2004)");
    }

    #endregion

    #region M16: Sequence Length In Range — Scanning invariant

    [Test]
    public void FindPreMiRnaHairpins_SequenceLength_WithinMinMaxRange()
    {
        const int min = 55;
        const int max = 120;

        var results = FindPreMiRnaHairpins(LargeHairpin110, minHairpinLength: min, maxHairpinLength: max).ToList();

        foreach (var pre in results)
        {
            Assert.That(pre.Sequence.Length, Is.InRange(min, max),
                $"Returned hairpin length {pre.Sequence.Length} must be within [{min}, {max}] — INV-1");
        }
    }

    #endregion

    #region M17: All Invariants Hold — Multiple sources

    [Test]
    public void FindPreMiRnaHairpins_AllInvariants_Verified()
    {
        const int min = 55;
        const int max = 120;
        const int matureLen = 22;

        var results = FindPreMiRnaHairpins(ValidHairpin57, min, max, matureLen).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        foreach (var pre in results)
        {
            Assert.Multiple(() =>
            {
                // INV-1: Length in range
                Assert.That(pre.Sequence.Length, Is.InRange(min, max),
                    "INV-1: Sequence length must be within [minHairpinLength, maxHairpinLength]");

                // INV-2: Start/End in bounds
                Assert.That(pre.Start, Is.GreaterThanOrEqualTo(0),
                    "INV-2: Start ≥ 0");
                Assert.That(pre.End, Is.LessThan(ValidHairpin57.Length),
                    "INV-2: End < input length");

                // INV-3: MatureSequence length
                Assert.That(pre.MatureSequence.Length, Is.GreaterThan(0)
                    .And.LessThanOrEqualTo(matureLen),
                    "INV-3: Mature length in (0, matureLength]");

                // INV-4: Star length matches mature
                Assert.That(pre.StarSequence.Length, Is.EqualTo(pre.MatureSequence.Length),
                    "INV-4: Star length must equal mature length");

                // INV-5: Structure length matches sequence
                Assert.That(pre.Structure.Length, Is.EqualTo(pre.Sequence.Length),
                    "INV-5: Structure length must equal sequence length");

                // INV-6: Only valid characters
                Assert.That(pre.Structure, Does.Match(@"^[(.)]+$"),
                    "INV-6: Structure must contain only '(', ')', '.' characters");

                // INV-7: Balanced parentheses
                int opens = pre.Structure.Count(c => c == '(');
                int closes = pre.Structure.Count(c => c == ')');
                Assert.That(opens, Is.EqualTo(closes),
                    "INV-7: Count of '(' must equal count of ')'");

                // INV-8: Negative free energy
                Assert.That(pre.FreeEnergy, Is.LessThan(0),
                    "INV-8: Free energy must be negative for stable hairpin");

                // INV-10: RNA uppercase
                Assert.That(pre.Sequence, Does.Match(@"^[AUGC]+$"),
                    "INV-10: Sequence must be uppercase RNA (A, U, G, C only)");
            });
        }
    }

    #endregion

    #region S1: Multiple Hairpins in Long Sequence — Biological miRNA clusters

    [Test]
    public void FindPreMiRnaHairpins_LongSequence_FindsMultipleCandidates()
    {
        // Two hairpins separated by 20 A's spacer
        var results = FindPreMiRnaHairpins(TwoHairpinSequence, minHairpinLength: 55).ToList();

        Assert.That(results, Has.Count.GreaterThanOrEqualTo(2),
            "Sequence with two distinct hairpin regions should yield ≥ 2 candidates — " +
            "miRNA clusters can contain multiple precursors (Bartel 2004)");
    }

    #endregion

    #region S2: MaxHairpinLength Respected — Parameter contract

    [Test]
    public void FindPreMiRnaHairpins_MaxHairpinLength_EnforcedAsUpperBound()
    {
        // 110 nt hairpin with max=80 should only return sub-windows ≤ 80
        var results = FindPreMiRnaHairpins(LargeHairpin110, minHairpinLength: 55, maxHairpinLength: 80).ToList();

        foreach (var pre in results)
        {
            Assert.That(pre.Sequence.Length, Is.LessThanOrEqualTo(80),
                "No returned hairpin may exceed maxHairpinLength parameter");
        }
    }

    #endregion

    #region S3: Custom MinHairpinLength — Parameter contract

    [Test]
    public void FindPreMiRnaHairpins_CustomMinHairpinLength_Applied()
    {
        // AnalyzeHairpin hardcodes n < 55 rejection, so minHairpinLength < 55
        // has no practical effect. Test above the 55 threshold instead:
        // ValidHairpin57 is exactly 57 nt.
        // With minHairpinLength=57, windows of 57 nt are scanned → detected
        // With minHairpinLength=58, 57 nt input is too short → empty
        var resultsAt57 = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 57).ToList();
        var resultsAt58 = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 58).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(resultsAt57, Is.Not.Empty,
                "57 nt hairpin must be detected with minHairpinLength=57");
            Assert.That(resultsAt58, Is.Empty,
                "57 nt hairpin must be rejected when minHairpinLength=58 (sequence too short)");
        });
    }

    #endregion

    #region S4: No Complementarity — No stem possible

    [Test]
    public void FindPreMiRnaHairpins_NoComplementarity_ReturnsEmpty()
    {
        var results = FindPreMiRnaHairpins(NoComplementarity, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "All-purine (AG) sequence has no self-complementarity → no stem → no hairpin");
    }

    #endregion

    #region Known Limitation: Real miRBase pre-miRNAs — consecutive-pairing model

    [Test]
    public void FindPreMiRnaHairpins_RealMiRBase_HsaMir21_NotDetected_KnownLimitation()
    {
        // hsa-mir-21 (MI0000077, miRBase v22) — 72 nt
        // Real pre-miRNA with internal mismatches and bulges.
        // Consecutive-pairing model gets only 8 pairs from ends (< 18 required).
        // This is a known limitation of the simplified consecutive-pairing model.
        var results = FindPreMiRnaHairpins(HsaMir21_MI0000077, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Known limitation: hsa-mir-21 (MI0000077) is a real pre-miRNA from miRBase " +
            "but is not detected because the consecutive-pairing model " +
            "cannot handle internal mismatches and bulges");
    }

    [Test]
    public void FindPreMiRnaHairpins_RealMiRBase_HsaLet7a1_NotDetected_KnownLimitation()
    {
        // hsa-let-7a-1 (MI0000060, miRBase v22) — 80 nt
        // Real pre-miRNA; consecutive pairing from ends yields only 5 pairs.
        var results = FindPreMiRnaHairpins(HsaLet7a1_MI0000060, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Known limitation: hsa-let-7a-1 (MI0000060) is a real pre-miRNA from miRBase " +
            "but is not detected because the consecutive-pairing model " +
            "cannot handle internal mismatches and bulges");
    }

    #endregion

    #region C1: Case-Insensitive Input — Robustness

    [Test]
    public void FindPreMiRnaHairpins_CaseInsensitive_SameResultAsUppercase()
    {
        // Mixed-case version of ValidHairpin57
        string mixedCase = "gcAuAgCuAgCuAgCuAgCuAgCuA" +
                           "gAaAuUu" +
                           "uAgCuAgCuAgCuAgCuAgCuAuGc";

        var upperResults = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55).ToList();
        var mixedResults = FindPreMiRnaHairpins(mixedCase, minHairpinLength: 55).ToList();

        Assert.That(mixedResults, Has.Count.EqualTo(upperResults.Count),
            "Mixed-case input must produce same number of candidates as uppercase");

        for (int i = 0; i < upperResults.Count; i++)
        {
            Assert.Multiple(() =>
            {
                Assert.That(mixedResults[i].Sequence, Is.EqualTo(upperResults[i].Sequence),
                    $"Candidate {i}: sequences must match (both converted to uppercase RNA)");
                Assert.That(mixedResults[i].FreeEnergy, Is.EqualTo(upperResults[i].FreeEnergy).Within(0.001),
                    $"Candidate {i}: energies must match");
                Assert.That(mixedResults[i].Structure, Is.EqualTo(upperResults[i].Structure),
                    $"Candidate {i}: structures must match");
            });
        }
    }

    #endregion

    #region C2: MatureLength Parameter — Parameterization

    [Test]
    public void FindPreMiRnaHairpins_CustomMatureLength_AffectsMatureAndStarSequence()
    {
        // With matureLength=18 and stem=23: matureEnd = min(18, 23) = 18
        // Mature = first 18 chars: "GCAUAGCUAGCUAGCUAG"
        // Star = last 18 chars (pos 39-56): "CUAGCUAGCUAGCUAUGC"
        const int customMatureLength = 18;
        const string expectedMature18 = "GCAUAGCUAGCUAGCUAG";
        const string expectedStar18 = "CUAGCUAGCUAGCUAUGC";

        var results = FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 55,
            matureLength: customMatureLength).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];

        Assert.Multiple(() =>
        {
            Assert.That(pre.MatureSequence, Is.EqualTo(expectedMature18),
                "Mature must be first 18 nt when matureLength=18");
            Assert.That(pre.MatureSequence.Length, Is.EqualTo(customMatureLength),
                "Mature length must equal custom matureLength parameter");
            Assert.That(pre.StarSequence, Is.EqualTo(expectedStar18),
                "Star must be last 18 nt when matureLength=18");
            Assert.That(pre.StarSequence.Length, Is.EqualTo(pre.MatureSequence.Length),
                "Star and mature must have equal length (duplex symmetry)");
        });
    }

    #endregion

    #region MFE-structure-based detection (opt-in) — reuses RNA-STRUCT-001 Zuker DP folder

    // Expected values below are the EXACT outputs of the validated RNA-STRUCT-001 engine
    // (RnaSecondaryStructure.CalculateMfeStructure / CalculateMinimumFreeEnergy, Turner 2004),
    // reused verbatim — NOT recomputed by an independent heuristic.
    // MFEI = AMFE/(G+C)% ; AMFE = 100·|ΔG°|/length — Zhang et al. (2006), Cell Mol Life Sci 63:246-254.

    // MF1 — A designed perfect-stem hairpin folds to a single dominant hairpin whose ΔG° equals
    //       the engine's CalculateMinimumFreeEnergy, with the expected stem/loop and AMFE/MFEI.
    [Test]
    public void AssessHairpinByMfe_PerfectHairpin_DerivesFeaturesFromRealMfeStructure()
    {
        // ValidHairpin57 folds (Turner 2004 DP) to 27×'(' 3×'.' 27×')', ΔG° = −48.48 kcal/mol.
        // GC% = 25/57 = 43.8596…; AMFE = 100·48.48/57 = 85.052632; MFEI = 85.052632/43.8596 = 1.939200.
        const double expectedDg = -48.48;
        const int expectedStemBp = 27;
        const int expectedLoopSize = 3;
        const double expectedAmfe = 85.052632;   // 100·48.48/57
        const double expectedMfei = 1.939200;    // AMFE/(G+C)%

        double engineDg = Seqeron.Genomics.Analysis.RnaSecondaryStructure
            .CalculateMinimumFreeEnergy(ValidHairpin57);
        var assessed = AssessHairpinByMfe(ValidHairpin57);

        Assert.That(assessed, Is.Not.Null, "Designed perfect hairpin must be accepted by the MFE-fold method.");
        var a = assessed!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(a.FreeEnergy, Is.EqualTo(expectedDg).Within(1e-9),
                "ΔG° must be the exact Turner-2004 MFE of the folded candidate (−48.48 kcal/mol).");
            Assert.That(a.FreeEnergy, Is.EqualTo(engineDg).Within(1e-9),
                "MFE-fold ΔG° must agree with the engine's own CalculateMinimumFreeEnergy.");
            Assert.That(a.StemBasePairs, Is.EqualTo(expectedStemBp),
                "Stem base-pair count must be read from the real MFE dot-bracket (27 pairs).");
            Assert.That(a.TerminalLoopSize, Is.EqualTo(expectedLoopSize),
                "Single terminal loop size must be 3 nt (the apical run of dots).");
            Assert.That(a.DotBracket, Is.EqualTo(new string('(', 27) + "..." + new string(')', 27)),
                "Dot-bracket must be the engine's MFE structure (27 paired, 3-nt apical loop).");
            Assert.That(a.Amfe, Is.EqualTo(expectedAmfe).Within(1e-6),
                "AMFE = 100·|ΔG°|/length (Zhang 2006).");
            Assert.That(a.Mfei, Is.EqualTo(expectedMfei).Within(1e-6),
                "MFEI = AMFE/(G+C)% (Zhang 2006).");
            Assert.That(a.Mfei, Is.GreaterThanOrEqualTo(0.85),
                "Genuine pre-miRNA hairpins have MFEI ≥ 0.85 (Zhang 2006).");
        });
    }

    // MF2 — The MFE-fold method DETECTS real miRBase pre-miRNAs that the consecutive-pairing
    //        heuristic rejects (the limitation this fix removes). hsa-mir-21, MI0000077.
    [Test]
    public void AssessHairpinByMfe_RealMiRBase_HsaMir21_DetectedByMfeFold()
    {
        // Folded by the engine: single hairpin (internal loops/bulges), apical loop 3 nt,
        // 32 stem base pairs, ΔG° = −35.13 kcal/mol. GC% = 35/72 = 48.6111…
        // AMFE = 100·35.13/72 = 48.791667; MFEI = 48.791667/48.6111 = 1.003714.
        const double expectedDg = -35.13;
        const int expectedStemBp = 32;
        const int expectedLoopSize = 3;
        const double expectedMfei = 1.003714;

        double engineDg = Seqeron.Genomics.Analysis.RnaSecondaryStructure
            .CalculateMinimumFreeEnergy(HsaMir21_MI0000077);
        var assessed = AssessHairpinByMfe(HsaMir21_MI0000077);

        Assert.That(assessed, Is.Not.Null,
            "hsa-mir-21 is a real pre-miRNA: the MFE-fold method MUST detect it (heuristic does not).");
        var a = assessed!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(a.FreeEnergy, Is.EqualTo(expectedDg).Within(1e-9),
                "ΔG° must equal the engine's Turner-2004 MFE for hsa-mir-21 (−35.13 kcal/mol).");
            Assert.That(a.FreeEnergy, Is.EqualTo(engineDg).Within(1e-9),
                "MFE-fold ΔG° must agree with CalculateMinimumFreeEnergy.");
            Assert.That(a.StemBasePairs, Is.EqualTo(expectedStemBp),
                "Stem base pairs counted from the real MFE structure = 32.");
            Assert.That(a.TerminalLoopSize, Is.EqualTo(expectedLoopSize),
                "Single terminal loop = 3 nt.");
            Assert.That(a.Mfei, Is.EqualTo(expectedMfei).Within(1e-6),
                "MFEI = 1.003714 (Zhang 2006), above the 0.85 cutoff.");
        });
    }

    // MF3 — hsa-let-7a-1 (MI0000060, 80 nt) is likewise detected by the MFE fold.
    [Test]
    public void AssessHairpinByMfe_RealMiRBase_HsaLet7a1_DetectedByMfeFold()
    {
        // Folded by the engine: single hairpin, apical loop 4 nt, 32 stem base pairs,
        // ΔG° = −34.31 kcal/mol. GC% = 34/80 = 42.5; AMFE = 100·34.31/80 = 42.887500;
        // MFEI = 42.887500/42.5 = 1.009118.
        const double expectedDg = -34.31;
        const int expectedStemBp = 32;
        const int expectedLoopSize = 4;
        const double expectedMfei = 1.009118;

        double engineDg = Seqeron.Genomics.Analysis.RnaSecondaryStructure
            .CalculateMinimumFreeEnergy(HsaLet7a1_MI0000060);
        var assessed = AssessHairpinByMfe(HsaLet7a1_MI0000060);

        Assert.That(assessed, Is.Not.Null,
            "hsa-let-7a-1 is a real pre-miRNA: the MFE-fold method MUST detect it (heuristic does not).");
        var a = assessed!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(a.FreeEnergy, Is.EqualTo(expectedDg).Within(1e-9),
                "ΔG° must equal the engine's Turner-2004 MFE for hsa-let-7a-1 (−34.31 kcal/mol).");
            Assert.That(a.FreeEnergy, Is.EqualTo(engineDg).Within(1e-9),
                "MFE-fold ΔG° must agree with CalculateMinimumFreeEnergy.");
            Assert.That(a.StemBasePairs, Is.EqualTo(expectedStemBp),
                "Stem base pairs counted from the real MFE structure = 32.");
            Assert.That(a.TerminalLoopSize, Is.EqualTo(expectedLoopSize),
                "Single terminal loop = 4 nt.");
            Assert.That(a.Mfei, Is.EqualTo(expectedMfei).Within(1e-6),
                "MFEI = 1.009118 (Zhang 2006), above the 0.85 cutoff.");
        });
    }

    // MF4 — Heuristic vs MFE divergence is locked: the SAME real pre-miRNA that the
    //        consecutive-pairing heuristic rejects is accepted by the MFE-fold method.
    [Test]
    public void HeuristicRejects_ButMfeFoldDetects_SameRealPreMiRna()
    {
        var heuristic = FindPreMiRnaHairpins(HsaMir21_MI0000077, minHairpinLength: 55).ToList();
        var mfe = AssessHairpinByMfe(HsaMir21_MI0000077);

        Assert.Multiple(() =>
        {
            Assert.That(heuristic, Is.Empty,
                "Default consecutive-pairing heuristic rejects hsa-mir-21 (documented limitation).");
            Assert.That(mfe, Is.Not.Null,
                "Opt-in MFE-fold method detects hsa-mir-21 — the limitation this fix removes.");
        });
    }

    // MF5 — A non-complementary sequence has ΔG° = 0 (no pairs) and is rejected.
    [Test]
    public void AssessHairpinByMfe_NoComplementarity_Rejected()
    {
        var assessed = AssessHairpinByMfe(NoComplementarity);

        Assert.That(assessed, Is.Null,
            "A sequence with no foldable stem (ΔG° = 0, all-dots structure) is not a hairpin.");
        Assert.That(Seqeron.Genomics.Analysis.RnaSecondaryStructure
                .CalculateMinimumFreeEnergy(NoComplementarity), Is.EqualTo(0.0).Within(1e-12),
            "Prerequisite: the engine reports ΔG° = 0 for the non-complementary sequence.");
    }

    // MF6 — Structure-based rejection: a multibranch fold is rejected EVEN with a strongly
    //        negative ΔG°, proving acceptance is on structure (single hairpin), not energy alone.
    [Test]
    public void AssessHairpinByMfe_MultibranchStructure_RejectedDespiteStrongEnergy()
    {
        // 120-nt 5S-rRNA-like sequence folds (engine) to a multibranch structure
        // (a ')' followed later by '(') with ΔG° = −47.04 kcal/mol — MORE negative than the
        // accepted ValidHairpin57 (−48.48 is comparable), yet it is NOT a single hairpin.
        const string fiveSLike =
            "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACG" +
            "CCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU";

        var mfe = Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure(fiveSLike);
        var assessed = AssessHairpinByMfe(fiveSLike);

        Assert.Multiple(() =>
        {
            Assert.That(mfe.FreeEnergy, Is.EqualTo(-47.04).Within(1e-9),
                "Engine ΔG° for the 5S-like sequence is −47.04 kcal/mol (a strong fold).");
            Assert.That(assessed, Is.Null,
                "Rejected because the MFE structure is multibranched, not a single dominant hairpin — " +
                "acceptance is structural, not energy-only.");
        });
    }

    // MF7 — FindPreMiRnaHairpinsByMfe window scan yields the accepted candidate for a designed hairpin.
    [Test]
    public void FindPreMiRnaHairpinsByMfe_DesignedHairpin_YieldsAcceptedCandidate()
    {
        var results = FindPreMiRnaHairpinsByMfe(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Window scan must yield the designed hairpin.");
        var full = results.First(r => r.Sequence.Length == 57);
        Assert.Multiple(() =>
        {
            Assert.That(full.FreeEnergy, Is.EqualTo(-48.48).Within(1e-9),
                "Yielded candidate carries the engine's exact ΔG° (−48.48 kcal/mol).");
            Assert.That(full.StemBasePairs, Is.EqualTo(27), "27 stem base pairs from the MFE structure.");
            Assert.That(full.Mfei, Is.GreaterThanOrEqualTo(0.85), "Accepted candidates satisfy MFEI ≥ 0.85.");
        });
    }

    // MF8 — minMfei gate: raising the cutoff above the candidate's MFEI rejects it.
    [Test]
    public void AssessHairpinByMfe_MfeiBelowThreshold_Rejected()
    {
        // hsa-let-7a-1 MFEI = 1.009118. A cutoff of 1.5 must reject it; 0.85 accepts it.
        var accepted = AssessHairpinByMfe(HsaLet7a1_MI0000060, minMfei: 0.85);
        var rejected = AssessHairpinByMfe(HsaLet7a1_MI0000060, minMfei: 1.5);

        Assert.Multiple(() =>
        {
            Assert.That(accepted, Is.Not.Null, "MFEI 1.009 ≥ 0.85 ⇒ accepted at default cutoff.");
            Assert.That(rejected, Is.Null, "MFEI 1.009 < 1.5 ⇒ rejected when the cutoff is raised.");
        });
    }

    // MF9 — CalculateMfeIndex formula: exact Zhang (2006) value and degenerate guards.
    [Test]
    public void CalculateMfeIndex_MatchesZhang2006Formula()
    {
        // ΔG° = −48.48, length 57, GC% = 43.8596491… ⇒ AMFE = 85.052631…, MFEI = 1.939200…
        double mfei = CalculateMfeIndex(-48.48, 57, 100.0 * 25 / 57);

        Assert.Multiple(() =>
        {
            Assert.That(mfei, Is.EqualTo(1.939200).Within(1e-6),
                "MFEI = (100·|ΔG°|/length)/(G+C)% — Zhang (2006).");
            Assert.That(CalculateMfeIndex(-48.48, 0, 43.0), Is.EqualTo(0.0),
                "Zero length ⇒ MFEI 0 (guard).");
            Assert.That(CalculateMfeIndex(-48.48, 57, 0.0), Is.EqualTo(0.0),
                "Zero GC% ⇒ MFEI 0 (guard, no division by zero).");
        });
    }

    // MF10 — null/empty input to the MFE-fold methods.
    [Test]
    public void MfeFoldMethods_NullOrEmpty_HandledGracefully()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AssessHairpinByMfe(null!), Is.Null, "Null candidate ⇒ null.");
            Assert.That(AssessHairpinByMfe(""), Is.Null, "Empty candidate ⇒ null.");
            Assert.That(FindPreMiRnaHairpinsByMfe(null!).ToList(), Is.Empty, "Null sequence ⇒ no candidates.");
            Assert.That(FindPreMiRnaHairpinsByMfe("").ToList(), Is.Empty, "Empty sequence ⇒ no candidates.");
            Assert.That(FindPreMiRnaHairpinsByMfe("ACGU").ToList(), Is.Empty,
                "Sequence shorter than minHairpinLength ⇒ no candidates.");
        });
    }

    #endregion

    #region Drosha/Dicer cleavage-site prediction (Han 2006 / Park 2011)

    // DD1 — Han et al. (2006), Cell 125:887: Drosha cleaves ~11 bp from the basal stem-ssRNA junction.
    // With basalJunction = 0, the Drosha 5' cut (5' end of the 5p mature) must be at index 0 + 11 = 11.
    [Test]
    public void PredictDroshaDicerCleavage_DroshaCut_Is11BpFromBasalJunction()
    {
        // Arrange — synthetic pri-miRNA: 11-nt lower stem then the miR-21 stem region (≥ 11+22+2 nt).
        const int basalJunction = 0;
        string pri = "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";

        // Act
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, basalJunction);

        // Assert
        Assert.That(cut, Is.Not.Null, "A well-formed pri-miRNA hairpin must yield a cleavage prediction.");
        Assert.That(cut!.Value.DroshaCut5Prime, Is.EqualTo(11),
            "Han (2006): Drosha cuts ~11 bp from the basal junction; junction 0 + 11 = index 11. " +
            "A wrong distance (e.g. 0 or 22) would land elsewhere.");
    }

    // DD2 — Park et al. (2011), Nature 475:201: Dicer 5' counting rule fixes the mature at ~22 nt.
    [Test]
    public void PredictDroshaDicerCleavage_MatureLength_Is22Nt()
    {
        // Arrange
        string pri = "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";

        // Act
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, 0)!.Value;

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(cut.MatureSequence.Length, Is.EqualTo(22),
                "Park (2011): the 5' counting rule (∼22 nt) fixes the mature length at 22 nt.");
            Assert.That(cut.MatureEnd - cut.MatureStart + 1, Is.EqualTo(22),
                "Mature span [Start,End] inclusive must be 22 nt.");
            Assert.That(cut.MatureStart, Is.EqualTo(cut.DroshaCut5Prime),
                "The 5p mature begins at the Drosha 5' cut.");
        });
    }

    // DD3 — Lee (2003)/Han (2006): each RNase III cut leaves a 2-nt 3' overhang. The 3' (3p) Drosha-
    // generated end sits 2 nt 3' of the position opposite the 5p mature end.
    [Test]
    public void PredictDroshaDicerCleavage_ThreePrimeOverhang_Is2Nt()
    {
        // Arrange
        string pri = "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";

        // Act
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, 0)!.Value;

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(cut.ThreePrimeOverhang, Is.EqualTo(2),
                "RNase III (Drosha/Dicer) leaves a 2-nt 3' overhang.");
            Assert.That(cut.DroshaCut3Prime - cut.MatureEnd, Is.EqualTo(2),
                "The 3p Drosha-generated end is 2 nt 3' of the 5p mature end (the 2-nt overhang).");
            Assert.That(cut.StarEnd - cut.StarStart + 1, Is.EqualTo(22),
                "The 3p (miRNA*) span is also ~22 nt.");
        });
    }

    // DD4 — miRBase cross-check (hsa-mir-21, MI0000077). With an 11-nt lower stem prepended so the
    // Drosha +11 ruler lands at the annotated miR-21-5p start, the predicted 5p mature must equal the
    // miRBase mature hsa-miR-21-5p sequence (MIMAT0000076) exactly. Source: mirbase.org/hairpin/MI0000077.
    [Test]
    public void PredictDroshaDicerCleavage_HsaMir21_MatchesMirBaseMature5p()
    {
        // Arrange — 11-nt lower stem + miR-21 stem (5p starts immediately after the lower stem).
        string pri = "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";
        const string mirBase5p = "UAGCUUAUCAGACUGAUGUUGA"; // MIMAT0000076, 22 nt

        // Act
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, 0)!.Value;

        // Assert
        Assert.That(cut.MatureSequence, Is.EqualTo(mirBase5p),
            "Predicted 5p mature must equal miRBase hsa-miR-21-5p (MIMAT0000076) exactly.");
    }

    // DD5 — Star (3p) span content is read from the correct 3' arm coordinates.
    [Test]
    public void PredictDroshaDicerCleavage_StarSpan_HasExpectedCoordinatesAndSequence()
    {
        // Arrange
        string pri = "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";

        // Act
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, 0)!.Value;

        // Hand-derived: matureStart=11, matureEnd=32, starEnd=34, starStart=13.
        Assert.Multiple(() =>
        {
            Assert.That(cut.MatureStart, Is.EqualTo(11), "matureStart = junction(0)+11.");
            Assert.That(cut.MatureEnd, Is.EqualTo(32), "matureEnd = 11+22-1.");
            Assert.That(cut.StarEnd, Is.EqualTo(34), "starEnd = matureEnd+2 (2-nt overhang).");
            Assert.That(cut.StarStart, Is.EqualTo(13), "starStart = starEnd-22+1.");
            Assert.That(cut.StarSequence, Is.EqualTo(pri.Substring(13, 22)),
                "StarSequence must be the 3p span pri[13..34].");
        });
    }

    // DD6 — T→U normalisation: DNA input yields RNA-alphabet mature/star sequences.
    [Test]
    public void PredictDroshaDicerCleavage_DnaInput_NormalisesToRna()
    {
        // Arrange — DNA spelling of the DD1 pri-miRNA.
        string priDna = "CCCCCCCCCCC" + "TAGCTTATCAGACTGATGTTGACTGTTGAATCTCATGGCAACACCAGTCGATGGGCTGT";

        // Act
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(priDna, 0)!.Value;

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(cut.MatureSequence, Is.EqualTo("UAGCUUAUCAGACUGAUGUUGA"),
                "T must be read as U; predicted mature is the RNA-alphabet hsa-miR-21-5p.");
            Assert.That(cut.MatureSequence, Does.Not.Contain("T"), "No T in normalised output.");
        });
    }

    // DD7 — null / empty input ⇒ null.
    [Test]
    public void PredictDroshaDicerCleavage_NullOrEmpty_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(null!, 0), Is.Null, "null ⇒ null.");
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage("", 0), Is.Null, "empty ⇒ null.");
        });
    }

    // DD8 — basal junction out of range ⇒ null.
    [Test]
    public void PredictDroshaDicerCleavage_JunctionOutOfRange_ReturnsNull()
    {
        string pri = "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";
        Assert.Multiple(() =>
        {
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, -1), Is.Null, "negative junction ⇒ null.");
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, pri.Length), Is.Null, "junction ≥ length ⇒ null.");
        });
    }

    // DD9 — too short to place the predicted cuts (junction+11+22+2 exceeds length) ⇒ null.
    [Test]
    public void PredictDroshaDicerCleavage_TooShortForCuts_ReturnsNull()
    {
        // 11 + 22 + 2 = 35 nt minimum needed from the junction; give 30 nt.
        string tooShort = new string('A', 30);
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(tooShort, 0), Is.Null,
            "When junction+11+22+2 exceeds the sequence, no full mature/star can be placed ⇒ null.");
    }

    // DD10 — CNNC confidence flag (Auyeung 2013): a C-N-N-C placed 16 nt 3' of the Drosha cut is detected.
    [Test]
    public void PredictDroshaDicerCleavage_CnncMotif_DetectedWhenPresentDownstream()
    {
        // Arrange — Drosha cut at index 11; place a CNNC starting at 11+16 = 27.
        // Build: 11 (lower stem) + 22 (mature) so that index 27 falls within, then ensure C..C at 27.
        var sb = new System.Text.StringBuilder(new string('A', 60));
        // index 27 = 'C', index 30 = 'C' (C-N-N-C with the two N's at 28,29)
        sb[27] = 'C'; sb[30] = 'C';
        string withCnnc = sb.ToString();

        var sbNo = new System.Text.StringBuilder(new string('A', 60)); // no C anywhere ⇒ no CNNC
        string noCnnc = sbNo.ToString();

        // Act
        var withFlag = MiRnaAnalyzer.PredictDroshaDicerCleavage(withCnnc, 0)!.Value;
        var withoutFlag = MiRnaAnalyzer.PredictDroshaDicerCleavage(noCnnc, 0)!.Value;

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(withFlag.HasCnncMotif, Is.True,
                "Auyeung (2013): a CNNC 16 nt 3' of the Drosha cut must set the confidence flag.");
            Assert.That(withoutFlag.HasCnncMotif, Is.False,
                "No CNNC downstream ⇒ flag false (optional signal, not required for prediction).");
        });
    }

    #endregion

    #region Helpers

    /// <summary>
    /// Builds a perfect hairpin with mixed-base Watson-Crick stem.
    /// Pattern: GCAU cycling for 5' stem, all-A loop, reverse complement for 3' stem.
    /// </summary>
    private static string BuildPerfectHairpin(int stemLength, int loopSize)
    {
        var bases = new[] { 'G', 'C', 'A', 'U' };
        var stem5 = new char[stemLength];
        for (int i = 0; i < stemLength; i++)
            stem5[i] = bases[i % 4];

        string loop = new string('A', loopSize);

        var stem3 = new char[stemLength];
        for (int i = 0; i < stemLength; i++)
        {
            stem3[i] = stem5[stemLength - 1 - i] switch
            {
                'A' => 'U',
                'U' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => stem5[stemLength - 1 - i]
            };
        }

        return new string(stem5) + loop + new string(stem3);
    }

    #endregion
}
