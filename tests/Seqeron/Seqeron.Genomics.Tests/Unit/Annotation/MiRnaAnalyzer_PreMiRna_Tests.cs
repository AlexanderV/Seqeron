using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Annotation;

/// <summary>
/// MIRNA-PRECURSOR-001: Pre-miRNA Hairpin Detection
/// Canonical test file for FindPreMiRnaHairpins (and indirectly AnalyzeHairpin).
/// Evidence: Bartel (2004), Ambros (2003), Krol (2004), Bartel (2009), miRBase.
/// </summary>
[TestFixture]
public class MiRnaAnalyzer_PreMiRna_Tests
{
    #region Test Data

    // Most tests of this file lock the CONSECUTIVE-PAIRING SCREEN (the pre-WP-7 default). The default of
    // FindPreMiRnaHairpins is now MFE-fold based; the screen is selected with useConsecutiveScreen: true.
    private static IEnumerable<PreMiRna> Screen(
        string sequence, int minHairpinLength = 55, int maxHairpinLength = 120, int matureLength = 22) =>
        FindPreMiRnaHairpins(sequence, minHairpinLength, maxHairpinLength, matureLength, useConsecutiveScreen: true);

    // ── Biologically plausible mixed-base hairpins ─────────────────────────
    // All test hairpins use mixed nucleotide composition (AUGC), matching
    // the diversity found in real pre-miRNAs. Mono-nucleotide sequences
    // (all-G stems, all-C stems) are avoided as they are biologically
    // implausible and can mask implementation bugs.

    // 25 bp Watson-Crick stem + 7 nt loop = 57 nt. The first two loop bases also pair with the last
    // two (G25·U31 wobble, A26-U30), so the uninterrupted consecutive stem is 27 bp closing a 3-nt
    // "AAU" loop — exactly the ViennaRNA 2.8 / RNA-STRUCT-001 MFE structure of this sequence.
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
    private const string TwoHairpinSequence =
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
        var results = Screen(null!).ToList();

        Assert.That(results, Is.Empty,
            "Null input must return empty enumerable, not throw");
    }

    #endregion

    #region M2: Empty Input — Defensive coding

    [Test]
    public void FindPreMiRnaHairpins_EmptyInput_ReturnsEmpty()
    {
        var results = Screen("").ToList();

        Assert.That(results, Is.Empty,
            "Empty string must return empty enumerable");
    }

    #endregion

    #region M3: Short Sequence — Bartel (2004)

    [Test]
    public void FindPreMiRnaHairpins_ShortSequence_ReturnsEmpty()
    {
        // 47 nt hairpin is below default minHairpinLength=55
        var results = Screen(ShortHairpin47, minHairpinLength: 55).ToList();

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
        var results = Screen(ValidHairpin57, minHairpinLength: 55).ToList();

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
        // - i=0, len=57 → stem=27, loop 3 (full hairpin) → accepted
        // - i=1, len=55 → stem=26, loop 3 (sub-window, same register) → accepted
        // Other windows rejected (pairing breaks at position 0/1).
        var results = Screen(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Has.Count.EqualTo(2),
            "57 nt input with min=55 must yield exactly 2 candidates: " +
            "full 57-nt window (stem=27) and 55-nt sub-window at offset 1 (stem=26)");

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
        // For 57-nt hairpin with stem=27: mature = first min(22, 27) = 22 nt
        // Hand-derived: first 22 chars of ValidHairpin57
        const string expectedMature = "GCAUAGCUAGCUAGCUAGCUAG";

        var results = Screen(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];

        Assert.Multiple(() =>
        {
            Assert.That(pre.MatureSequence, Is.EqualTo(expectedMature),
                "Mature must be the first 22 nt of 5' arm — " +
                "matureEnd = min(matureLength=22, stemLength=27) = 22 (Bartel 2009)");
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

        var results = Screen(ValidHairpin57, minHairpinLength: 55).ToList();

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
        // For the 57-nt hairpin the uninterrupted stem is 27 bp (loop "AAU", 3 nt):
        // Structure = 27×'(' + 3×'.' + 27×')' — identical to ViennaRNA 2.8 RNA.fold(ValidHairpin57)
        // "(((((((((((((((((((((((((((...)))))))))))))))))))))))))))" (−48.40, Vienna's 0.1-kcal tables).
        const string expectedStructure =
            "(((((((((((((((((((((((((((...)))))))))))))))))))))))))))";

        var results = Screen(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        var pre = results[0];

        Assert.Multiple(() =>
        {
            Assert.That(pre.Structure, Is.EqualTo(expectedStructure),
                "Structure must be exactly 27×'(' + 3×'.' + 27×')' — the full consecutive stem " +
                "(loop ≥ 3, NNDB minimum), = ViennaRNA MFE structure");
            Assert.That(pre.Structure.Length, Is.EqualTo(57),
                "Structure length must equal sequence length (INV-5)");
        });
    }

    #endregion

    #region M9: Dot-Bracket Balanced — Notation definition

    [Test]
    public void FindPreMiRnaHairpins_Structure_BalancedParentheses()
    {
        var results = Screen(ValidHairpin57, minHairpinLength: 55).ToList();

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
        // Hand-calculated from Turner 2004 NNDB parameters for ValidHairpin57 (stem=27, loop "AAU"):
        //   Stacking (26 stacks): 6×GC/CG(-3.42) + 1×CA/GU(-2.11) + 1×AU/UA(-1.10) + 6×UA/AU(-1.33)
        //                       + 5×AG/UC(-2.08) + 5×CU/GA(-2.08) + 1×AG/UU(-0.55) + 1×GA/UU(-1.27) = -54.33
        //   AU end penalty (helix end A26-U30 closing the loop): +0.45 (NNDB "per AU end", hairpin example 5.3)
        //   Hairpin initiation(3): +5.40, no first-mismatch term for triloops; "AAAUU" is not a special triloop
        //   Outer G-C end: no penalty.   Total: -48.48 kcal/mol
        // = RnaSecondaryStructure.CalculateMinimumFreeEnergy(ValidHairpin57) (this structure is its MFE).
        // ViennaRNA 2.8 eval_structure: -48.40 (its tables are rounded to 0.1 kcal, TerminalAU 0.50).
        const double expectedEnergy = -48.48;

        var results = Screen(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Not.Empty, "Prerequisite: must detect hairpin");

        Assert.That(results[0].FreeEnergy, Is.EqualTo(expectedEnergy).Within(0.01),
            "Free energy must match hand-calculated Turner 2004 NNDB value (-48.48 kcal/mol) — " +
            "26 stacks + AU end + triloop initiation");

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
        // Hairpin A: 57 nt, stem=27, loop=3
        //   Hand-calculated Turner energy: -48.48 kcal/mol (see M10 for breakdown)
        //
        // Hairpin B: 55 nt, stem=20, loop=15
        //   Stacking (19 pairs): 5×GC/CG(-3.42) + 1×CA/GU(-2.11) + 1×AU/UA(-1.10)
        //                       + 4×UA/AU(-1.33) + 4×AG/UC(-2.08) + 4×CU/GA(-2.08) = -42.27
        //   Loop init(15):    +6.90
        //   Terminal mismatch(UAGA, closing U-A): -1.10
        //   AU/GU penalty (closing U-A): +0.45
        //   Total: -36.02 kcal/mol (ViennaRNA 2.8 eval_structure: -35.80 with its rounded tables)
        const double expectedEnergyLong = -48.48;
        const double expectedEnergyShort = -36.02;

        string stem20Hairpin55 =
            "GCAUAGCUAGCUAGCUAGCU" +     // 20 nt 5' stem
            "AUGCAUGCAUGCAUG" +           // 15 nt loop
            "AGCUAGCUAGCUAGCUAUGC";       // 20 nt 3' stem (reverse complement)

        var resultsLong = Screen(ValidHairpin57, minHairpinLength: 55).ToList();
        var resultsShort = Screen(stem20Hairpin55, minHairpinLength: 55).ToList();

        Assert.That(resultsLong, Is.Not.Empty, "27 bp stem hairpin (57 nt) must be detected");
        Assert.That(resultsShort, Is.Not.Empty, "20 bp stem hairpin (55 nt) must be detected");

        Assert.Multiple(() =>
        {
            Assert.That(resultsLong[0].FreeEnergy, Is.EqualTo(expectedEnergyLong).Within(0.01),
                "Stem-27 energy must match hand-calculated Turner 2004 value");
            Assert.That(resultsShort[0].FreeEnergy, Is.EqualTo(expectedEnergyShort).Within(0.01),
                "Stem-20 energy must match hand-calculated Turner 2004 value");
            Assert.That(resultsLong[0].FreeEnergy, Is.LessThan(resultsShort[0].FreeEnergy),
                "Longer stem (27 bp, -48.48) must be more negative than shorter stem (20 bp, -36.02) — " +
                "reflects greater thermodynamic stability (Turner 2004)");
        });
    }

    #endregion

    #region M12: Stem Too Short Rejected — Krol (2004)

    [Test]
    public void FindPreMiRnaHairpins_StemTooShort_ReturnsEmpty()
    {
        // 15 bp stem + 25 nt loop = 55 nt total
        // Rejection must come from stemLength (15) < 18, not from length filter
        var results = Screen(ShortStemHairpin55, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "Sequence with only 15 bp stem must be rejected — " +
            "consecutive screen requires ≥ 18 bp stem (operational screen parameter). " +
            "Note: sequence is 55 nt, so this tests stem rejection, not the window-length filter");
    }

    #endregion

    #region M13: Loop Size Boundary Rejected — Bartel (2004)

    [Test]
    public void FindPreMiRnaHairpins_LoopTooLarge_Rejected()
    {
        // Note: Loop-too-small (< 3) is structurally unreachable because the stem never
        // extends past the 3-nt NNDB minimum loop. Test the reachable arm: loop > 25 → rejected.
        // 18 bp mixed stem + 30 nt loop = 66 nt. Stem=18, loop=30 > 25
        // Loop uses all-A to prevent self-complementary pairing within the loop
        string largeLoop =
            "GCAUAGCUAGCUAGCUAG" +             // 18 nt 5' stem (mixed)
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA" + // 30 nt loop (A-A never pairs)
            "CUAGCUAGCUAGCUAUGC";              // 18 nt 3' stem (reverse complement)
        var results = Screen(largeLoop, minHairpinLength: 55).ToList();

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
        var results = Screen(DnaHairpin57, minHairpinLength: 55).ToList();

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
        var results = Screen(WobbleHairpin57, minHairpinLength: 55).ToList();

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

        var results = Screen(LargeHairpin110, minHairpinLength: min, maxHairpinLength: max).ToList();

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

        var results = Screen(ValidHairpin57, min, max, matureLen).ToList();

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
        var results = Screen(TwoHairpinSequence, minHairpinLength: 55).ToList();

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
        var results = Screen(LargeHairpin110, minHairpinLength: 55, maxHairpinLength: 80).ToList();

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
        // ValidHairpin57 is exactly 57 nt.
        // With minHairpinLength=57, windows of 57 nt are scanned → detected
        // With minHairpinLength=58, 57 nt input is too short → empty
        var resultsAt57 = Screen(ValidHairpin57, minHairpinLength: 57).ToList();
        var resultsAt58 = Screen(ValidHairpin57, minHairpinLength: 58).ToList();

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
        var results = Screen(NoComplementarity, minHairpinLength: 55).ToList();

        Assert.That(results, Is.Empty,
            "All-purine (AG) sequence has no self-complementarity → no stem → no hairpin");
    }

    #endregion

    #region Real miRBase pre-miRNAs — screen limitation and MFE default (WP-7 / F40)

    // The consecutive-pairing SCREEN (useConsecutiveScreen: true) still cannot see bulged pre-miRNAs.
    [Test]
    public void Screen_RealMiRBase_HsaMir21_NotDetected_KnownLimitation()
    {
        // hsa-mir-21 (MI0000077, miRBase v22) — 72 nt; consecutive pairing from the ends gives only 8 pairs (< 18).
        Assert.That(Screen(HsaMir21_MI0000077, minHairpinLength: 55), Is.Empty);
    }

    [Test]
    public void Screen_RealMiRBase_HsaLet7a1_NotDetected_KnownLimitation()
    {
        // hsa-let-7a-1 (MI0000060) — 80 nt; consecutive pairing from the ends yields only 5 pairs.
        Assert.That(Screen(HsaLet7a1_MI0000060, minHairpinLength: 55), Is.Empty);
    }

    // F40 — the DEFAULT is the MFE fold. Oracle: ViennaRNA 2.8 RNA.fold of the same window — identical
    // dot-bracket; ΔG° engine vs Vienna: hsa-mir-21 [1,71] −35.17 / −34.70, hsa-let-7a-1 [0,79] −34.31 / −34.20
    // (all accepted windows of the 5 pri-miRNA hairpins have the Vienna structure, |ΔG − Vienna| ≤ 0.53).
    [Test]
    public void FindPreMiRnaHairpins_Default_DetectsBulgedMiRBasePreMiRnas_WithViennaStructure()
    {
        var mir21 = FindPreMiRnaHairpins(HsaMir21_MI0000077).Single(p => p.Start == 1 && p.End == 71);
        var let7 = FindPreMiRnaHairpins(HsaLet7a1_MI0000060).Single(p => p.Start == 0 && p.End == 79);
        Assert.Multiple(() =>
        {
            Assert.That(mir21.Structure, Is.EqualTo("(((((((((((((((.(((((.(((((.((((.((...))))))))))).)))))))))))))))))))).")); // RNA.fold
            Assert.That(let7.Structure, Is.EqualTo("(((((.(((((((((((((((((((((.....(((...((((....)))).)))))))))))))))))))))))))))))")); // RNA.fold
            Assert.That(mir21.FreeEnergy, Is.EqualTo(-34.70).Within(0.5), "ViennaRNA 2.8 RNA.fold −34.70");
            Assert.That(let7.FreeEnergy, Is.EqualTo(-34.20).Within(0.15), "ViennaRNA 2.8 RNA.fold −34.20");
            Assert.That(mir21.FreeEnergy, Is.EqualTo(
                Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMinimumFreeEnergy(mir21.Sequence)).Within(1e-9));
        });
    }

    // F40 — star = Dicer 2-nt 3' overhang partner (miRDeep2 find_star), not the blunt partner: the star ends
    // 2 nt beyond the partner of the mature start (clipped to the window end), and its start is
    // partner(mature end − 2 − off) − off. hsa-mir-21 window [1,71]: mature starts at the outermost pair.
    [Test]
    public void FindPreMiRnaHairpins_Default_MatureStar_FollowFindStarOverhangRule()
    {
        var mir21 = FindPreMiRnaHairpins(HsaMir21_MI0000077).Single(p => p.Start == 1 && p.End == 71);
        Assert.Multiple(() =>
        {
            Assert.That(mir21.MatureSequence, Is.EqualTo("GUCGGGUAGCUUAUCAGACUGA"), "22 nt from the 5' base of the outermost pair");
            Assert.That(mir21.StarSequence, Is.EqualTo("AGUCGAUGGGCUGUCUGACA"),
                "find_star: ends at partner(mature start)+2 clipped to the window; starts at partner(end−2−off)−off");
            Assert.That(mir21.StarSequence, Is.Not.EqualTo(mir21.Sequence.Substring(mir21.Sequence.Length - 22)),
                "not the blunt last-22-nt partner of the screen");
        });
        // An internal window with 3' overhang inside the window: the star carries the 2-nt overhang.
        var mir155 = FindPreMiRnaHairpins("AUAUA" + "CUGUUAAUGCUAAUCGUGAUAGGGGUUUUUGCCUCCAACUGACUCCUACAUAUUAGCAUUAACAG" + "AUAUAUA")
            .First(p => p.Start == 5 && p.End == 69);
        Assert.That(mir155.MatureSequence, Is.EqualTo("CUGUUAAUGCUAAUCGUGAUAG"));
    }

    // Both modes fill the same record; the screen keeps the blunt labels.
    [Test]
    public void FindPreMiRnaHairpins_ScreenOption_KeepsBluntLabelsAndConsecutiveStem()
    {
        var screen = FindPreMiRnaHairpins(ValidHairpin57, 55, 120, 22, useConsecutiveScreen: true)
            .Single(p => p.Start == 0 && p.End == 56);
        Assert.That(screen.StarSequence, Is.EqualTo(ValidHairpin57.Substring(57 - 22)), "blunt partner (no overhang)");
        Assert.That(screen.MatureSequence, Is.EqualTo(ValidHairpin57.Substring(0, 22)));
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

        var upperResults = Screen(ValidHairpin57, minHairpinLength: 55).ToList();
        var mixedResults = Screen(mixedCase, minHairpinLength: 55).ToList();

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
        // With matureLength=18 and stem=27: matureEnd = min(18, 27) = 18
        // Mature = first 18 chars: "GCAUAGCUAGCUAGCUAG"
        // Star = last 18 chars (pos 39-56): "CUAGCUAGCUAGCUAUGC"
        const int customMatureLength = 18;
        const string expectedMature18 = "GCAUAGCUAGCUAGCUAG";
        const string expectedStar18 = "CUAGCUAGCUAGCUAUGC";

        var results = Screen(ValidHairpin57, minHairpinLength: 55,
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

    #region B13 F8–F10 — canonical stem-loop energy/structure, full consecutive stem, parameter contract

    // F8 — the heuristic's FreeEnergy is the Turner 2004 energy of its reported structure computed by the
    // canonical RnaSecondaryStructure loop model (FindStemLoops → CalculateStemEnergy + CalculateHairpinLoopEnergy),
    // which includes the NNDB special hairpin loops. Fixture: 20-bp stem closed by C-G around the UUCG loop
    // (NNDB turner_2004_tetraloop_dg.txt "CUUCGG 3.7" — the experimental total replaces initiation + mismatch).
    //   NNDB hand calc: 19 stacks = 4×GC/CG(-3.42) + CA/GU(-2.11) + AU/UA(-1.10) + 4×UA/AU(-1.33) + 4×AG/UC(-2.08)
    //                   + 3×CU/GA(-2.08) + GU/CA(-2.24) + UC/AG(-2.35) = -41.36; + special tetraloop 3.70 = -37.66.
    //   ViennaRNA 2.8 eval_structure: -37.60 (hairpin term 3.70 identical; stacks rounded to 0.1 kcal).
    //   The former inline formula (initiation(4) 5.6 + mismatch CUGG -1.0) would give -36.76 (the old scan
    //   never reached it: n < 55 floor and the n/2 − 5 stem cap).
    private const string UucgTetraloopHairpin44 =
        "GCAUAGCUAGCUAGCUAGUC" + "UUCG" + "GACUAGCUAGCUAGCUAUGC";

    [Test]
    public void FindPreMiRnaHairpins_SpecialTetraloop_UsesCanonicalHairpinLoopModel()
    {
        var full = Screen(UucgTetraloopHairpin44, minHairpinLength: 40)
            .Single(p => p.Start == 0 && p.End == 43);

        Assert.Multiple(() =>
        {
            Assert.That(full.Structure, Is.EqualTo(new string('(', 20) + "...." + new string(')', 20)));
            Assert.That(full.FreeEnergy, Is.EqualTo(-37.66).Within(1e-9),
                "NNDB: 19 stacks (-41.36) + special UUCG tetraloop closed by C-G (+3.70).");
            Assert.That(full.FreeEnergy, Is.EqualTo(-37.60).Within(0.25),
                "ViennaRNA 2.8 eval_structure = -37.60 (0.1-kcal rounded stacks).");
            Assert.That(full.FreeEnergy, Is.EqualTo(Seqeron.Genomics.Analysis.RnaSecondaryStructure
                .CalculateMinimumFreeEnergy(UucgTetraloopHairpin44)).Within(1e-9),
                "This perfect hairpin is its own MFE structure ⇒ energy equals the RNA-STRUCT-001 engine.");
        });
    }

    // F8 — GA first-mismatch bonus (NNDB hairpin rule; RnaSecondaryStructure UU_GA_MismatchBonus −0.9), missing
    // from the former inline formula. 23-bp stem closed by C-G around an 11-nt "GAAAAAAAAAA" loop (57 nt):
    //   NNDB: 22 stacks −49.10 + initiation(11) 6.60 + terminal mismatch CGAG −1.40 + GA bonus −0.90 = −44.80.
    //   ViennaRNA 2.8: −44.70 (hairpin term 4.30 identical). Former code output: −43.90 (no GA bonus).
    private const string GaFirstMismatchHairpin57 =
        "GCAUAGCUAGCUAGCUAGCUAGC" + "GAAAAAAAAAA" + "GCUAGCUAGCUAGCUAGCUAUGC";

    [Test]
    public void FindPreMiRnaHairpins_GaFirstMismatch_BonusApplied()
    {
        var full = Screen(GaFirstMismatchHairpin57).Single(p => p.Start == 0 && p.End == 56);

        Assert.Multiple(() =>
        {
            Assert.That(full.Structure, Is.EqualTo(new string('(', 23) + new string('.', 11) + new string(')', 23)));
            Assert.That(full.FreeEnergy, Is.EqualTo(-44.80).Within(1e-9), "NNDB incl. GA first-mismatch bonus");
            Assert.That(full.FreeEnergy, Is.EqualTo(-44.70).Within(0.25),
                "ViennaRNA 2.8 eval_structure −44.70 (its stacks are rounded to 0.1 kcal; ≤ 0.24 on all fixtures)");
        });
    }

    // F8/F9 — for perfect hairpins the reported consecutive stem is the MFE structure, so the reported
    // structure and energy must coincide with the canonical folder. NNDB hand values:
    //   ValidHairpin57: -54.33 stacks + 0.45 (A-U end) + 5.4 (triloop)                 = -48.48 (Vienna -48.40)
    //   WobbleHairpin57: -46.10 stacks + 0.45 (outer G-U end) + 0.45 (A-U end) + 5.4    = -39.80 (Vienna -39.80)
    //   ShortHairpin47 (min 40): -44.54 stacks + 0.45 (A-U end) + 5.4                  = -38.69 (Vienna -38.60)
    [TestCase(ValidHairpin57, 55, -48.48, 27)]
    [TestCase(WobbleHairpin57, 55, -39.80, 27)]
    [TestCase(ShortHairpin47, 40, -38.69, 22)]
    public void FindPreMiRnaHairpins_PerfectHairpin_StructureAndEnergyEqualCanonicalMfe(
        string hairpin, int minLength, double nndbEnergy, int stem)
    {
        var full = Screen(hairpin, minHairpinLength: minLength)
            .Single(p => p.Start == 0 && p.End == hairpin.Length - 1);
        var mfe = Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure(hairpin);

        Assert.Multiple(() =>
        {
            Assert.That(full.Structure, Is.EqualTo(
                new string('(', stem) + new string('.', hairpin.Length - 2 * stem) + new string(')', stem)));
            Assert.That(full.Structure, Is.EqualTo(mfe.DotBracket), "consecutive stem == MFE structure");
            Assert.That(full.FreeEnergy, Is.EqualTo(nndbEnergy).Within(1e-9), "NNDB hand calculation");
            Assert.That(full.FreeEnergy, Is.EqualTo(mfe.FreeEnergy).Within(1e-9), "== engine MFE");
        });
    }

    // F9 — the stem is no longer truncated at min(n/2 − 5, 35): loops of 3–9 nt are reported (Bartel 2004:
    // 3–25 nt), and stems longer than 35 bp are not mis-reported as a > 25-nt "loop" and rejected.
    // LargeHairpin110 = 50-bp perfect stem + 10-nt all-A loop: previously stem capped at 35 ⇒ "loop" 40 ⇒ rejected.
    [Test]
    public void FindPreMiRnaHairpins_LongPerfectStem_NotTruncated()
    {
        var full = Screen(LargeHairpin110).Single(p => p.Start == 0 && p.End == 109);

        Assert.That(full.Structure, Is.EqualTo(new string('(', 50) + new string('.', 10) + new string(')', 50)));
    }

    // F10 — minHairpinLength below 55 is honoured (formerly a hard-coded n < 55 rejection ignored it).
    [Test]
    public void FindPreMiRnaHairpins_MinHairpinLengthBelow55_Honoured()
    {
        var hits = Screen(ShortHairpin47, minHairpinLength: 40).ToList();

        Assert.That(hits.Select(h => (h.Start, h.End)),
            Is.EqualTo(new[] { (0, 46), (1, 45), (2, 44), (3, 43) }),
            "47-nt hairpin (22-bp consecutive stem, loop 3) and its same-register sub-windows down to 41 nt " +
            "(stem ≥ 19 ≥ 18) are accepted when minHairpinLength = 40; 39-nt window 4..42 has stem 18 but " +
            "length 39 < 40.");
    }

    // F10 — argument validation is eager (before enumeration) and explicit.
    [Test]
    public void PreMiRnaScans_InvalidLengths_ThrowEagerly()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(() => Screen(ValidHairpin57, matureLength: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Screen(ValidHairpin57, minHairpinLength: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => Screen(ValidHairpin57, minHairpinLength: -5));
            // default (MFE) mode validates eagerly as well
            Assert.Throws<ArgumentOutOfRangeException>(() => FindPreMiRnaHairpins(ValidHairpin57, matureLength: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => FindPreMiRnaHairpins(ValidHairpin57, minHairpinLength: 0));
            Assert.Throws<ArgumentOutOfRangeException>(() => FindPreMiRnaHairpinsByMfe(ValidHairpin57, minHairpinLength: 0));
            Assert.That(Screen(ValidHairpin57, minHairpinLength: 60, maxHairpinLength: 55), Is.Empty,
                "max < min is an empty scan (documented), not an error");
        });
    }

    // GC% for MFEI is the canonical GC fraction (non-nucleotide symbols excluded, Biopython gc_fraction
    // "remove"); AMFE keeps the full folded length (Zhang 2006: MFE/length·100).
    [Test]
    public void AssessHairpinByMfe_NonNucleotideSymbol_ExcludedFromGcPercent()
    {
        string withN = ValidHairpin57 + "N";
        var a = AssessHairpinByMfe(withN);
        double dg = Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMinimumFreeEnergy(withN);

        Assert.That(a, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(a!.Value.FreeEnergy, Is.EqualTo(dg).Within(1e-12));
            Assert.That(a.Value.Amfe, Is.EqualTo(100.0 * Math.Abs(dg) / 58).Within(1e-12));
            Assert.That(a.Value.Mfei, Is.EqualTo(100.0 * Math.Abs(dg) / 58 / (100.0 * 25 / 57)).Within(1e-12),
                "GC% = 25 G+C of 57 nucleotides (N excluded).");
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
        // Folded by the engine: single hairpin, apical loop 4 nt, 33 stem base pairs,
        // ΔG° = −34.31 kcal/mol. The structure is identical to ViennaRNA 2.7.2 RNA.fold (d2) —
        // (((((.(((((((((((((((((((((.....(((...((((....)))).))))))))))))))))))))))))))))) —
        // and its energy with the NNDB-precision Turner 2004 parameters is −34.31 (RNA-STRUCT-001 Evidence). GC% = 34/80 = 42.5; AMFE = 100·34.31/80 = 42.887500;
        // MFEI = 42.887500/42.5 = 1.009118.
        const double expectedDg = -34.31;
        const int expectedStemBp = 33;
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
                "Stem base pairs counted from the real MFE structure = 33.");
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
        var heuristic = Screen(HsaMir21_MI0000077, minHairpinLength: 55).ToList();
        var mfe = AssessHairpinByMfe(HsaMir21_MI0000077);
        Assert.That(FindPreMiRnaHairpins(HsaMir21_MI0000077), Is.Not.Empty, "the default is the MFE scan since F40");

        Assert.Multiple(() =>
        {
            Assert.That(heuristic, Is.Empty,
                "The consecutive-pairing screen (useConsecutiveScreen: true) rejects hsa-mir-21 (documented limitation).");
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
        // (a ')' followed later by '(') with ΔG° = −50.57 kcal/mol — comparable to the
        // accepted ValidHairpin57 (−48.48), yet it is NOT a single hairpin. −50.57 = ViennaRNA 2.7.2
        // eval_structure (d2, NNDB-precision Turner 2004 parameters) of the engine's structure; RNA.fold's own
        // optimum (−50.79) closes a CCC hairpin that the NNDB all-C loop penalty (+1.5) excludes here.
        const string fiveSLike =
            "UGCCUGGCGGCCGUAGCGCGGUGGUCCCACCUGACCCCAUGCCGAACUCAGAAGUGAAACG" +
            "CCGUAGCGCCGAUGGUAGUGUGGGGUCUCCCCAUGCGAGAGUAGGGAACUGCCAGGCAU";

        var mfe = Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure(fiveSLike);
        var assessed = AssessHairpinByMfe(fiveSLike);

        Assert.Multiple(() =>
        {
            Assert.That(mfe.FreeEnergy, Is.EqualTo(-50.57).Within(1e-9),
                "Engine ΔG° for the 5S-like sequence is −50.57 kcal/mol (a strong fold).");
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
            // R3-A9: public inputs are validated (doc: "length must be > 0"; Zhang 2006 MFEI = (MFE/length×100)/GC%
            // is undefined for length ≤ 0 or GC% = 0 and meaningless for non-finite ΔG° or GC% outside [0,100]).
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMfeIndex(-48.48, 0, 43.0), "length 0");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMfeIndex(-48.48, -5, 43.0), "negative length");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMfeIndex(-48.48, 57, 0.0), "GC% 0 (division by zero)");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMfeIndex(-48.48, 57, double.NaN), "GC% NaN");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMfeIndex(-48.48, 57, 100.01), "GC% > 100");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMfeIndex(-48.48, 57, -1.0), "GC% < 0");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMfeIndex(double.NaN, 57, 43.0), "ΔG° NaN");
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateMfeIndex(double.NegativeInfinity, 57, 43.0), "ΔG° -∞");
            // The folder-derived path keeps the documented GC = 0 value (MFEI 0): an all-A/U hairpin is never accepted.
            Assert.That(ExtractPreMiRnaFeatures(new string('A', 20) + "UUU" + new string('U', 20))!.Value.Mfei, Is.EqualTo(0.0));
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

    // ── MFE-fold structural-gate rejection branches (each Stage-A threshold in isolation) ──
    // These exercise the stem-bp and loop-size acceptance gates DIRECTLY: each rejected
    // candidate folds to a clean single hairpin with MFEI ≫ 0.85, so rejection is provably
    // due to the structural threshold under test, NOT the energy/MFEI gate. dG° values are the
    // engine's Turner-2004 output; AMFE/MFEI are hand-derived from the Zhang (2006) formula.

    // A perfect Watson–Crick hairpin of the given stem length and loop, mixed GCAU composition.
    private static string PerfectHairpin(int stem, int loop) => BuildPerfectHairpin(stem, loop);

    // MF11 — Stem-base-pair gate (Ambros et al. 2003, ≥16 complementary bases). A 14-bp clean
    //        hairpin (MFEI 1.59 ≫ 0.85) is REJECTED solely on stem < 16; a 16-bp one is accepted.
    [Test]
    public void AssessHairpinByMfe_StemBelow16_RejectedOnStemGateNotEnergy()
    {
        string stem14 = PerfectHairpin(14, 5); // folds to 14 bp + 5-nt loop, ΔG°=−25.44, MFEI=1.59
        string stem16 = PerfectHairpin(16, 5); // folds to 16 bp + 5-nt loop, accepted

        var rejected = AssessHairpinByMfe(stem14);
        var accepted = AssessHairpinByMfe(stem16);

        Assert.Multiple(() =>
        {
            // Prove the MFEI gate is NOT the reason: a 14-bp hairpin has MFEI = 1.59 ≥ 0.85.
            double mfei14 = CalculateMfeIndex(-25.44, stem14.Length, 100.0 * 16 / stem14.Length);
            Assert.That(mfei14, Is.EqualTo(1.59).Within(1e-6),
                "14-bp hairpin MFEI = 1.59 (Zhang 2006), comfortably above the 0.85 cutoff.");
            Assert.That(rejected, Is.Null,
                "Stem of 14 bp < 16 ⇒ rejected on the stem-base-pair gate (Ambros 2003), not on energy.");
            Assert.That(accepted, Is.Not.Null, "Stem of 16 bp meets the ≥16 gate ⇒ accepted.");
            Assert.That(accepted!.Value.StemBasePairs, Is.EqualTo(16),
                "Accepted candidate has exactly 16 stem base pairs (the Ambros 2003 boundary).");
        });
    }

    // MF12 — Terminal-loop gate (Bartel 2004, loop ≤ ~25 nt). A 20-bp hairpin with a 26-nt apical
    //        loop (MFEI 1.732 ≫ 0.85) is REJECTED solely on loop > 25.
    [Test]
    public void AssessHairpinByMfe_LoopAbove25_RejectedOnLoopGateNotEnergy()
    {
        string bigLoop = PerfectHairpin(20, 26); // 20 bp stem + 26-nt loop, ΔG°=−34.64, MFEI=1.732

        var rejected = AssessHairpinByMfe(bigLoop);

        Assert.Multiple(() =>
        {
            double mfei = CalculateMfeIndex(-34.64, bigLoop.Length, 100.0 * 20 / bigLoop.Length);
            Assert.That(mfei, Is.EqualTo(1.732).Within(1e-6),
                "Big-loop hairpin MFEI = 1.732 (Zhang 2006), well above the 0.85 cutoff.");
            Assert.That(rejected, Is.Null,
                "Apical loop of 26 nt > 25 ⇒ rejected on the loop-size gate (Bartel 2004), not on energy.");
        });
    }

    // MF13 — MFEI acceptance gate is a TIGHT boundary at the candidate's own MFEI. hsa-let-7a-1
    //        folds to MFEI = 1.0091176…; a cutoff just below accepts, just above rejects.
    [Test]
    public void AssessHairpinByMfe_MfeiGate_IsTightBoundaryAtCandidateMfei()
    {
        // Independently re-derived: ΔG°=−34.31, len 80, GC=34 ⇒ MFEI = 100·34.31/80 / 42.5 = 1.00911764…
        const double let7Mfei = 1.0091176470588237;

        var atOrBelow = AssessHairpinByMfe(HsaLet7a1_MI0000060, minMfei: let7Mfei);          // == ⇒ accept
        var justBelow = AssessHairpinByMfe(HsaLet7a1_MI0000060, minMfei: let7Mfei - 1e-6);   // <  ⇒ accept
        var justAbove = AssessHairpinByMfe(HsaLet7a1_MI0000060, minMfei: let7Mfei + 1e-6);   // >  ⇒ reject

        Assert.Multiple(() =>
        {
            Assert.That(atOrBelow, Is.Not.Null, "MFEI == cutoff ⇒ accepted (gate is ≥, inclusive).");
            Assert.That(justBelow, Is.Not.Null, "Cutoff just below the candidate's MFEI ⇒ accepted.");
            Assert.That(justAbove, Is.Null, "Cutoff just above the candidate's MFEI ⇒ rejected.");
            Assert.That(atOrBelow!.Value.Mfei, Is.EqualTo(let7Mfei).Within(1e-9),
                "The candidate's own MFEI equals the hand-derived Zhang (2006) value.");
        });
    }

    // MF14 — Non-ACGU symbols: the MFE-fold methods upper-case and read T as U; other letters
    //        (e.g. N) are non-pairing, so they do not create spurious hairpins. A DNA-spelled
    //        copy of the designed hairpin must give the SAME assessment as its RNA spelling.
    [Test]
    public void AssessHairpinByMfe_DnaAndNonAcgu_NormalisedConsistently()
    {
        var rna = AssessHairpinByMfe(ValidHairpin57);
        var dna = AssessHairpinByMfe(DnaHairpin57); // same hairpin spelled with T

        Assert.Multiple(() =>
        {
            Assert.That(rna, Is.Not.Null, "RNA spelling accepted.");
            Assert.That(dna, Is.Not.Null, "DNA spelling (T→U) must fold identically and be accepted.");
            Assert.That(dna!.Value.FreeEnergy, Is.EqualTo(rna!.Value.FreeEnergy).Within(1e-9),
                "T is read as U ⇒ identical ΔG° for DNA and RNA spellings.");
            Assert.That(dna.Value.Mfei, Is.EqualTo(rna.Value.Mfei).Within(1e-9),
                "Identical MFEI for DNA and RNA spellings.");
            Assert.That(dna.Value.Sequence, Does.Not.Contain("T"),
                "Folded candidate is normalised to RNA (no residual T).");
            // A run of non-pairing N's cannot form a stem ⇒ no hairpin, no throw.
            Assert.That(AssessHairpinByMfe(new string('N', 60)), Is.Null,
                "All-N sequence has no foldable stem ⇒ not a hairpin (handled, no throw).");
        });
    }

    // MF15 — FindPreMiRnaHairpinsByMfe length-window bounds: a candidate shorter than
    //        minHairpinLength is never assessed; a maxHairpinLength below the only hairpin's
    //        length yields nothing (the window upper bound is enforced).
    [Test]
    public void FindPreMiRnaHairpinsByMfe_LengthWindowBounds_Enforced()
    {
        // ValidHairpin57 is 57 nt. With min=58 the only candidate window is too short ⇒ empty.
        var tooHighMin = FindPreMiRnaHairpinsByMfe(ValidHairpin57, minHairpinLength: 58).ToList();
        // With max=56 (< 57) no 57-nt window is ever formed ⇒ the full hairpin is never assessed.
        var cappedMax = FindPreMiRnaHairpinsByMfe(ValidHairpin57, minHairpinLength: 55, maxHairpinLength: 56).ToList();
        // Default window admits the 57-nt hairpin.
        var ok = FindPreMiRnaHairpinsByMfe(ValidHairpin57, minHairpinLength: 55).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(tooHighMin, Is.Empty, "min length 58 > sequence length 57 ⇒ no candidates.");
            Assert.That(cappedMax.Any(r => r.Sequence.Length == 57), Is.False,
                "maxHairpinLength 56 caps the window below 57 ⇒ the full hairpin is never yielded.");
            Assert.That(ok.Any(r => r.Sequence.Length == 57), Is.True,
                "Default window admits the 57-nt designed hairpin.");
        });
    }

    #endregion

    #region Drosha/Dicer cleavage-site prediction (Han 2006 / Park 2011)

    // Oracle data (retrieved 2026-11-25): hairpin + mature coordinates from mirbase.org/hairpin/<MI>,
    // genomic flanks from the UCSC hg38 REST API (api.genome.ucsc.edu getData/sequence, minus-strand loci
    // reverse-complemented), structures from ViennaRNA 2.8 RNA.fold (37 °C, Turner 2004), expected cuts from
    // an independent Python re-implementation of the rules: 11 bp counted on the fold (Han 2006) + 22-nt
    // Dicer ruler (Park 2011) + the miRDeep2 find_star 2-nt-overhang rule (rajewsky-lab/mirdeep2
    // src/miRDeep2_core_algorithm.pl). The basal junction passed is the pair 11 bp below the miRBase 5p
    // start on that fold. Flank = nt of genomic sequence on each side of the miRBase hairpin.
    // Columns: name, flank, pri, RNAfold structure, junction, 5p[s,e], 3p[s,e] (5' counting), 5p end and
    // 3p start (3' counting), basal UG, CNNC.
    private static readonly object[][] RealPriMiRnas =
    {
        new object[] { "hsa-mir-21 MI0000077", 30,
            "CAUCGUGACAUCUCCAUGGCUGUACCACCUUGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACAUUUUGGUAUCUUUCAUCUGACCAUCCAUAU",
            "..(((((.......)))))..((((((...((((((((((((((((.(((((.(((((.((((.((...))))))))))).)))))))))))))))))))))...)))))).....................",
            23, 37, 58, 76, 96, 59, 75, false, false },
        new object[] { "hsa-let-7a-1 MI0000060", 30,
            "CACCAUUCACCCUGGAUGUUCUCUUCACUGUGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUAACGUGAUAGAAAAGUCUGCAUCCAGGCGGU",
            ".......(..((((((((((((..((((..(((((.(((((((((((((((((((((.....(((...((((....)))).)))))))))))))))))))))))))))))..)))).))).......)))))))))..).",
            20, 36, 57, 85, 106, 57, 85, false, false },
        new object[] { "hsa-mir-16-1 MI0000070", 30,
            "CUUUUUAUUCAUAGCUCUUAUGAUAGCAAUGUCAGCAGUGCCUUAGCAGCACGUAAAUAUUGGCGUUAAGAUUCUAAAAUUAUCUCCAGUAUUAACUGUGCUGCUGAAGUAAGGUUGACCAUACUCUACAGUUGUGUUUUAAUGUAUAU",
            ".........(((((((.....((.......((((((..(((.((((((((((((.((((((((.(.(((..........))).).)))))))).)).)))))))))).)))..)))))).....))...))))))).............",
            22, 43, 64, 85, 107, 63, 86, true, true },
        new object[] { "hsa-mir-155 MI0000681", 30,
            "GAGUGCUGAAGGCUUGCUGUAGGCUGUAUGCUGUUAAUGCUAAUCGUGAUAGGGGUUUUUGCCUCCAACUGACUCCUACAUAUUAGCAUUAACAGUGUAUGAUGCCUGUUACUAGCAUUCACAUG",
            "((((((((..(((.......((((.(((((((((((((((((((.(((.((((((((.(((....)))..))))))))))))))))))))))))))))))...)))))))..)))))))).....",
            21, 33, 54, 74, 93, 56, 72, false, false },
        new object[] { "hsa-mir-30a MI0000088", 20,
            "UAUAUUGCUGUUGACAGUGAGCGACUGUAAACAUCCUCGACUGGAAGCUGUGAAGCCACAGAUGGGCUUUCAGUCGGAUGUUUGCAGCUGCCUACUGCCUCGGACUUCAAG",
            ".......(((..(.(((((.(((.((((((((((((..((((((((((((((....))).....))))))))))))))))))))))).))).))))))..)))........",
            9, 25, 46, 68, 87, 48, 66, true, true },
    };

    // DD1 — 5' counting (default) on real pri-miRNAs with the RNAfold structure. Offsets of the
    // prediction vs miRBase (hairpin 1-based, 5p start/end | 3p start/end): mir-21 0/−2 | +1/−1;
    // let-7a-1 +1/0 | −1/0 (miRBase 5p start is unpaired; the Drosha ruler lands on the next paired base);
    // mir-16-1 0/0 | 0/+1; mir-30a 0/0 | +2/0; mir-155 0/−2 | +2/0. miRBase boundaries are the dominant
    // read ends (current miRBase lists 23–24-nt 5p isoforms for mir-21, let-7a, mir-155).
    [TestCaseSource(nameof(RealPriMiRnas))]
    public void PredictDroshaDicerCleavage_RealPriMiRna_FivePrimeCounting_MatchesOracle(
        string name, int flank, string pri, string fold, int junction,
        int m5s, int m5e, int m3s, int m3e, int m5eThree, int m3sThree, bool ug, bool cnnc)
    {
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, junction, fold)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That((cut.MatureStart, cut.MatureEnd), Is.EqualTo((m5s, m5e)), name + " 5p");
            Assert.That((cut.StarStart, cut.StarEnd), Is.EqualTo((m3s, m3e)), name + " 3p");
            Assert.That(cut.DroshaCut5Prime, Is.EqualTo(m5s));
            Assert.That(cut.DroshaCut3Prime, Is.EqualTo(m3e));
            Assert.That(cut.MatureSequence, Is.EqualTo(pri.Substring(m5s, 22)), "Park 2011: 22-nt 5p");
            Assert.That(cut.HasBasalUgMotif, Is.EqualTo(ug), name + " basal UG (−14/−13)");
            Assert.That(cut.HasCnncMotif, Is.EqualTo(cnnc), name + " CNNC (+16..+18)");
            Assert.That(cut.Structure, Is.EqualTo(fold));
            Assert.That(cut.DicerRule, Is.EqualTo(DicerCountingRule.FivePrimeCounting));
        });
    }

    // DD2b — apical Drosha ruler (Zeng, Yi & Cullen 2005, EMBO J 24:138: ~22 bp from the loop/stem junction),
    // an opt-in alternative to the 11-bp basal ruler (Han 2006). On the five RNAfold-folded pri-miRNAs the
    // 5p cut is the 5' base of the 22nd pair counted from the pair closing the terminal loop; offsets vs the
    // miRBase 5p start (= the basal-ruler cut of DD1): mir-21 +3, let-7a-1 +6, mir-16-1 +1, mir-155 +3,
    // mir-30a +3 nt (independent Python walk of the same dot-bracket, scratch a10.py). The apical ruler lands
    // 1-6 pairs apical of the annotated Dicer-substrate end for every one of them, i.e. the basal ruler is the
    // better predictor here; the apical ruler is offered because the two publications disagree.
    [TestCase(0, 40, 93, 3)]
    [TestCase(1, 42, 100, 6)]
    [TestCase(2, 44, 106, 1)]
    [TestCase(3, 36, 90, 3)]
    [TestCase(4, 28, 84, 3)]
    public void PredictDroshaDicerCleavage_ApicalRuler_PlacesCut22PairsFromApicalJunction(
        int row, int cut5, int cut3, int offsetVsMirBase)
    {
        var r = RealPriMiRnas[row];
        string pri = (string)r[2], fold = (string)r[3];
        int junction = (int)r[4], m5s = (int)r[5];
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(
            pri, junction, fold, DicerCountingRule.FivePrimeCounting, DroshaCountingRule.ApicalJunction)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(cut.DroshaCut5Prime, Is.EqualTo(cut5));
            Assert.That(cut.DroshaCut3Prime, Is.EqualTo(cut3));
            Assert.That(cut.DroshaCut5Prime - m5s, Is.EqualTo(offsetVsMirBase));
            Assert.That(cut.DroshaRule, Is.EqualTo(DroshaCountingRule.ApicalJunction));
            // exactly 22 pairs between the apical junction pair and the cut (inclusive of the cut pair)
            var pairs = RnaSecondaryStructure.ParseDotBracket(fold).ToDictionary(x => x.Item1, x => x.Item2);
            Assert.That(pairs[cut.ApicalJunction], Is.GreaterThan(cut.ApicalJunction));
            Assert.That(fold.Substring(cut.ApicalJunction + 1, pairs[cut.ApicalJunction] - cut.ApicalJunction - 1),
                Does.Not.Contain("("), "apical junction pair closes a hairpin loop");
            Assert.That(fold.Substring(cut.DroshaCut5Prime, cut.ApicalJunction - cut.DroshaCut5Prime + 1).Count(c => c == '('), Is.EqualTo(22));
            // the basal rule on the same input keeps its cut and reports the same apical junction
            var basal = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, junction, fold)!.Value;
            Assert.That(basal.DroshaCut5Prime, Is.EqualTo(m5s));
            Assert.That(basal.DroshaRule, Is.EqualTo(DroshaCountingRule.BasalJunction));
            Assert.That(basal.ApicalJunction, Is.EqualTo(cut.ApicalJunction));
        });
    }

    // The apical rule needs a 22-pair stem and a defined rule value.
    [Test]
    public void PredictDroshaDicerCleavage_ApicalRuler_ShortStemNullAndBadEnumThrows()
    {
        string stem = new string('G', 15) + "AAAAAAAA" + new string('C', 15); // 15 bp
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(stem, 0, new string('(', 15) + "........" + new string(')', 15),
            DicerCountingRule.FivePrimeCounting, DroshaCountingRule.ApicalJunction), Is.Null);
        var r = RealPriMiRnas[4];
        Assert.Throws<ArgumentOutOfRangeException>(() => MiRnaAnalyzer.PredictDroshaDicerCleavage(
            (string)r[2], (int)r[4], (string)r[3], DicerCountingRule.FivePrimeCounting, (DroshaCountingRule)9));
    }

    // DD2 — 3' counting (Park 2011): the 3p product is 22 nt from the Drosha 3' end; the 5p end is the
    // partner of the 3p start + 2. Offsets vs miRBase (5p end | 3p start): mir-21 −1 | 0; let-7a-1 0 | −1;
    // mir-16-1 −1 | +1; mir-30a +2 | 0; mir-155 0 | 0.
    [TestCaseSource(nameof(RealPriMiRnas))]
    public void PredictDroshaDicerCleavage_RealPriMiRna_ThreePrimeCounting_MatchesOracle(
        string name, int flank, string pri, string fold, int junction,
        int m5s, int m5e, int m3s, int m3e, int m5eThree, int m3sThree, bool ug, bool cnnc)
    {
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, junction, fold, DicerCountingRule.ThreePrimeCounting)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That((cut.MatureStart, cut.MatureEnd), Is.EqualTo((m5s, m5eThree)), name + " 5p");
            Assert.That((cut.StarStart, cut.StarEnd), Is.EqualTo((m3sThree, m3e)), name + " 3p");
            Assert.That(cut.StarSequence.Length, Is.EqualTo(22), "3' counting: 22-nt 3p");
        });
    }

    // DD3 — the 2-argument overload folds with the canonical MFE (RNA-STRUCT-001). For these four pri-miRNAs
    // the canonical fold equals RNAfold 2.8 and the cuts equal DD1. (pri-mir-16-1 is excluded: the canonical
    // folder returns a different structure whose energy it scores −42.23 kcal/mol where RNAfold's eval gives
    // −40.80 — an energy-model discrepancy reported to the RNA-structure owner, not a cleavage-rule issue.)
    private static IEnumerable<object[]> RealPriMiRnasCanonicalFoldAgrees() =>
        RealPriMiRnas.Where(r => !((string)r[0]).StartsWith("hsa-mir-16-1", StringComparison.Ordinal));

    [TestCaseSource(nameof(RealPriMiRnasCanonicalFoldAgrees))]
    public void PredictDroshaDicerCleavage_CanonicalFold_AgreesWithRnaFoldCuts(
        string name, int flank, string pri, string fold, int junction,
        int m5s, int m5e, int m3s, int m3e, int m5eThree, int m3sThree, bool ug, bool cnnc)
    {
        string canonical = Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure(pri).DotBracket;
        Assert.That(canonical, Is.EqualTo(fold), name + ": canonical MFE equals RNAfold");
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, junction)!.Value;
        Assert.That((cut.MatureStart, cut.MatureEnd, cut.StarStart, cut.StarEnd), Is.EqualTo((m5s, m5e, m3s, m3e)));
    }

    // Synthetic perfect hairpin: 10-nt 5' flank, 35-bp stem, 8-nt loop, 35-bp stem, 25-nt 3' flank (113 nt).
    // Junction = 10. Pair k (0-based) = (10+k, 87−k).
    private static string SyntheticPri(char[]? edit = null) => edit is null
        ? new string('A', 10) + new string('G', 35) + new string('A', 8) + new string('C', 35) + new string('A', 25)
        : new string(edit);
    private const string SyntheticFold = ".........." + "(((((((((((((((((((((((((((((((((((" + "........"
        + ")))))))))))))))))))))))))))))))))))" + ".........................";

    // DD4 — hand-derived geometry on the perfect stem: 12th pair 5' base = 10+11 = 21 (Han 11 bp), its
    // partner 87−11 = 76, Drosha 3' end 76+2 = 78; 5p [21,42] (Park 22 nt); 3p start = partner(42−2) =
    // 87−30 = 57 (miRDeep2 find_star) ⇒ 3p [57,78], 22 nt, on the opposite arm.
    [Test]
    public void PredictDroshaDicerCleavage_PerfectStem_HandDerivedGeometry()
    {
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 10, SyntheticFold)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(cut.DroshaCut5Prime, Is.EqualTo(21));
            Assert.That(cut.DroshaCut3Prime, Is.EqualTo(78));
            Assert.That((cut.MatureStart, cut.MatureEnd), Is.EqualTo((21, 42)));
            Assert.That((cut.StarStart, cut.StarEnd), Is.EqualTo((57, 78)));
            Assert.That(cut.ThreePrimeOverhang, Is.EqualTo(2));
            Assert.That(cut.StarSequence, Is.EqualTo(new string('C', 22)), "3p lies on the 3' arm");
        });
        var three = MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 10, SyntheticFold, DicerCountingRule.ThreePrimeCounting)!.Value;
        Assert.That((three.MatureStart, three.MatureEnd, three.StarStart, three.StarEnd), Is.EqualTo((21, 42, 57, 78)),
            "on a perfect duplex both counting rules give the same 22/22 products");
    }

    // DD5 — Han's ruler counts base PAIRS: a 1-nt bulge on the 5' arm of the lower stem (after pair 5)
    // moves the 5p start one nucleotide 3' (12th pair 5' base = 22, not junction+11 = 21; partner 77,
    // Drosha 3' end 79); a 1-nt bulge on the 3' arm of the lower stem changes neither (12th pair = (21, 76)).
    [Test]
    public void PredictDroshaDicerCleavage_LowerStemBulge_CountsBasePairsNotNucleotides()
    {
        string fold5 = ".........." + "(((((" + "." + "((((((((((((((((((((((((((((((" + "........"
            + ")))))))))))))))))))))))))))))))))))" + "........................";
        string fold3 = ".........." + "(((((((((((((((((((((((((((((((((((" + "........"
            + "))))))))))))))))))))))))))))))" + "." + ")))))" + "........................";
        var a = MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 10, fold5)!.Value;
        var b = MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 10, fold3)!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(a.DroshaCut5Prime, Is.EqualTo(22), "5'-arm bulge is skipped: 11 bp = 12 nt");
            Assert.That(a.DroshaCut3Prime, Is.EqualTo(79), "partner(22) = 77, + 2");
            Assert.That(b.DroshaCut5Prime, Is.EqualTo(21));
            Assert.That(b.DroshaCut3Prime, Is.EqualTo(78), "partner(21) = 76 + 2");
        });
    }

    // DD6 — miRDeep2 find_star bulge compensation: with pair (40, 57) opened, the 5p anchor matureEnd−2 = 40
    // is unpaired, so the search moves inward (off = 1 → partner(39) = 58) and compensates (58 − 1 = 57).
    [Test]
    public void PredictDroshaDicerCleavage_UnpairedAnchor_FindStarOffsetCompensated()
    {
        char[] f = SyntheticFold.ToCharArray();
        f[40] = '.'; f[57] = '.'; // open pair (40, 57) — 57 = 87−30
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 10, new string(f))!.Value;
        Assert.That(cut.StarStart, Is.EqualTo(57));
    }

    // DD7 — DNA / lower case is normalised (T→U, upper) before folding and indexing.
    [Test]
    public void PredictDroshaDicerCleavage_DnaLowerCase_SameAsRna()
    {
        string pri = (string)RealPriMiRnas[0][2], fold = (string)RealPriMiRnas[0][3];
        int j = (int)RealPriMiRnas[0][4];
        var rna = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, j, fold)!.Value;
        var dna = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri.Replace('U', 'T').ToLowerInvariant(), j, fold)!.Value;
        Assert.That(dna, Is.EqualTo(rna));
        Assert.That(dna.MatureSequence, Does.Not.Contain("T"));
    }

    // DD8 — null / empty / junction out of range ⇒ null.
    [Test]
    public void PredictDroshaDicerCleavage_NullEmptyOrOutOfRange_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(null!, 0), Is.Null);
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage("", 0), Is.Null);
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), -1, SyntheticFold), Is.Null);
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 113, SyntheticFold), Is.Null);
        });
    }

    // DD9 — the junction must be the 5'-arm base of a pair: an unpaired flank base or a 3'-arm base ⇒ null.
    [Test]
    public void PredictDroshaDicerCleavage_JunctionNotFivePrimeOfPair_ReturnsNull()
    {
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 5, SyntheticFold), Is.Null, "unpaired");
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 87, SyntheticFold), Is.Null, "3' arm");
    }

    // DD10 — fewer than 12 stem pairs above the junction (no room for the 11-bp lower stem) ⇒ null; a
    // sequence with no stem at all (poly-A folds to all dots) ⇒ null.
    [Test]
    public void PredictDroshaDicerCleavage_StemTooShortOrAbsent_ReturnsNull()
    {
        // Junction 15 = pair 5: 30 pairs from it ⇒ valid (5p [26,47], 3p [52,73]); junction 35 = pair 25:
        // only 10 pairs from it ⇒ null.
        var ok = MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 15, SyntheticFold);
        Assert.That((ok!.Value.MatureStart, ok.Value.MatureEnd, ok.Value.StarStart, ok.Value.StarEnd), Is.EqualTo((26, 47, 52, 73)),
            "anchor 45 is in the loop: find_star moves to 44 (off 1), partner 53 − 1 = 52");
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 35, SyntheticFold), Is.Null);
        Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(new string('A', 60), 0), Is.Null, "no stem");
    }

    // DD11 — products not on opposite arms ⇒ null: junction 25 (pair 15) puts the 5p at [36,57], through
    // the loop (45–52) into the 3' arm; find_star then returns 42 < 57.
    [Test]
    public void PredictDroshaDicerCleavage_ProductsOverlap_ReturnsNull()
        => Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 25, SyntheticFold), Is.Null);

    // DD12 — invalid / length-mismatched structure ⇒ ArgumentException.
    [Test]
    public void PredictDroshaDicerCleavage_InvalidStructure_Throws()
    {
        Assert.Throws<ArgumentException>(() => MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 10, "((.."));
        Assert.Throws<ArgumentException>(() => MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 10, SyntheticFold[..^1] + "("));
    }

    // DD13 — Auyeung 2013 motif positions on the perfect stem (5p cut 21, 3p cut 78): basal UG = U at
    // 21−14 = 7, G at 8; CNNC starting at +16..+18 downstream of the 3p cut = indices 94..96. A CNNC at
    // +17 (95, 98) and +18 (96, 99) is detected; +15 (93) and +19 (97) are not. A CNNC 16–18 nt 3' of the
    // 5p cut (old reference point, indices 37–39 inside the 5p) is not a CNNC flag.
    [Test]
    public void PredictDroshaDicerCleavage_AuyeungMotifPositions()
    {
        bool Cnnc(int start)
        {
            char[] s = SyntheticPri().ToCharArray();
            s[start] = 'C'; s[start + 3] = 'C';
            return MiRnaAnalyzer.PredictDroshaDicerCleavage(new string(s), 10, SyntheticFold)!.Value.HasCnncMotif;
        }
        char[] ug = SyntheticPri().ToCharArray(); ug[7] = 'U'; ug[8] = 'G';
        char[] ugOff = SyntheticPri().ToCharArray(); ugOff[6] = 'U'; ugOff[7] = 'G';
        Assert.Multiple(() =>
        {
            Assert.That(Cnnc(94) && Cnnc(95) && Cnnc(96), Is.True, "+16, +17, +18");
            Assert.That(Cnnc(93), Is.False, "+15");
            Assert.That(Cnnc(97), Is.False, "+19");
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(SyntheticPri(), 10, SyntheticFold)!.Value.HasCnncMotif, Is.False);
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(new string(ug), 10, SyntheticFold)!.Value.HasBasalUgMotif, Is.True);
            Assert.That(MiRnaAnalyzer.PredictDroshaDicerCleavage(new string(ugOff), 10, SyntheticFold)!.Value.HasBasalUgMotif, Is.False);
        });
    }

    // DD14 — regression for the removed linear model: on hsa-mir-21 the old code put the "star" at
    // matureEnd+2 on the SAME (5') arm; the fold places it on the 3' arm within 1 nt of miRBase
    // hsa-miR-21-3p (MIMAT0004494, hairpin 46–68 ⇒ pri 75–97 at flank 30).
    [Test]
    public void PredictDroshaDicerCleavage_HsaMir21_ThreePrimeArmNotLinearShift()
    {
        string pri = (string)RealPriMiRnas[0][2];
        var cut = MiRnaAnalyzer.PredictDroshaDicerCleavage(pri, 23, (string)RealPriMiRnas[0][3])!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(cut.StarStart, Is.GreaterThan(cut.MatureEnd), "3p on the opposite arm");
            Assert.That(cut.StarSequence, Is.EqualTo("AACACCAGUCGAUGGGCUGUC"), "miRBase 3p CAACACCAGUCGAUGGGCUGUCU, +1/−1");
            Assert.That(cut.StarStart - cut.MatureStart, Is.Not.EqualTo(2), "not the old linear +2 shift");
        });
    }

    #endregion

    #region Trained natural-vs-background classifier (CL1–CL20)

    // ── Real human pre-miRNA precursors from miRBase (public domain) ───────────────────────
    // Retrieved verbatim from mirbase.org/hairpin/<MI>. Same 13 positives + fixed seed used by
    // the bundled training run (Evidence §Training); re-running the pipeline here reproduces the
    // held-out metric deterministically.
    private static readonly string[] MiRBasePositives =
    {
        "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA",                       // hsa-mir-21
        "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA",               // hsa-let-7a-1
        "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGUUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",            // hsa-let-7b
        "CCUUGGCCAUGUAAAAGUGCUUACAGUGCAGGUAGCUUUUUGAGAUCUACUGCAAUGUAAGCACUUCUUACAUUACCAUGG",              // hsa-mir-106a
        "GGAGAGGAGGCAAGAUGCUGGCAUAGCUGUUGAACUGGGAACCUGCUAUGCCAACAUAUUGCCAUCUUUCC",                        // hsa-mir-31
        "UGUGGGAUGAGGUAGUAGAUUGUAUAGUUUUAGGGUCAUACCCCAUCUUGGAGAUAACUAUACAGUCUACUGUCUUUCCCACG",            // hsa-let-7f-2
        "GUCAGCAGUGCCUUAGCAGCACGUAAAUAUUGGCGUUAAGAUUCUAAAAUUAUCUCCAGUAUUAACUGUGCUGCUGAAGUAAGGUUGAC",      // hsa-mir-16-1
        "CUCCGGUGCCUACUGAGCUGAUAUCAGUUCUCAUUUUACACACUGGCUCAGUUCAGCAGGAACAGGAG",                           // hsa-mir-24-1
        "CUGGGGGCUCCAAAGUGCUGUUCGUGCAGGUAGUGUGAUUACCCAACCUACUGCUGAGCUAGCACUUCCCGAGCCCCCGG",              // hsa-mir-93
        "UGGCCGAUUUUGGCACUAGCACAUUUUUGCUUGUGUCUCUCCGCUCUGAGCAAUCAUGUGCAGUGCCAAUAUGGGAAA",                 // hsa-mir-96
        "AGGAUUCUGCUCAUGCCAGGGUGAGGUAGUAAGUUGUAUUGUUGUGGGGUAGGGAUAUUAGGCCCCAAUUAGAAGAUAACUAUACAACUUACUACUUUCCCUGGUGUGUGGCAUAUUCA", // hsa-mir-98
        "AAUCUAAAGACAACAUUUCUGCACACACACCAGACUAUGGAAGCCAGUGUGUGGAAAUGCUUCUGCUAGAUU",                       // hsa-mir-147a
        "AGGAAGCUUCUGGAGAUCCUGCUCCGUCGCCCCAGUGUUCAGACUACCUGUUCAGGACAAUGCCGUUGUACAGUAGUCUGCACAUUGGUUAGACUGGGCAAGGGAGAGCA", // hsa-mir-199a-2
    };

    private const string HsaMir21 = "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA";
    private const string HsaLet7a1 = "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA";

    private const int TrainSeed = 20060101;
    private const int ShufflesPerPositive = 4;

    // CL1 — Feature extractor pins exactly for hsa-mir-21 (Evidence §Exact feature values).
    // Values are the RNA-STRUCT-001 Turner-2004 MFE structure + composition, computed
    // independently of the trained weights. A wrong fold/AMFE/MFEI/%paired would fail here.
    [Test]
    public void ExtractPreMiRnaFeatures_HsaMir21_PinsExactStructureFeatures()
    {
        var f = ExtractPreMiRnaFeatures(HsaMir21);

        Assert.That(f, Is.Not.Null, "hsa-mir-21 is a valid non-empty sequence.");
        Assert.Multiple(() =>
        {
            Assert.That(f!.Value.FreeEnergy, Is.EqualTo(-35.13).Within(1e-9),
                "ΔG° of the MFE structure (RNA-STRUCT-001 Turner 2004) for hsa-mir-21.");
            Assert.That(f.Value.Amfe, Is.EqualTo(48.79166666666667).Within(1e-9),
                "AMFE = 100·|ΔG°|/length (Zhang 2006): 100·35.13/72.");
            Assert.That(f.Value.Mfei, Is.EqualTo(1.0037142857142858).Within(1e-9),
                "MFEI = AMFE/(G+C)% (Zhang 2006); GC%=48.61 ⇒ 48.7917/48.611.");
            Assert.That(f.Value.GcContent, Is.EqualTo(0.48611111111111116).Within(1e-9),
                "G+C fraction: 35 of 72 nt.");
            Assert.That(f.Value.PairedFraction, Is.EqualTo(0.8888888888888888).Within(1e-9),
                "%paired = 2·bp/length = 2·32/72 from the MFE structure.");
            Assert.That(f.Value.StemBasePairs, Is.EqualTo(32),
                "Base-pair count of the MFE structure.");
            Assert.That(f.Value.BasePairingPropensity, Is.EqualTo(32.0 / 72.0).Within(1e-15),
                "dP = |bp|/L (Schultes 1999; miPred, Ng & Mishra 2007) = PairedFraction/2.");
            Assert.That(f.Value.Length, Is.EqualTo(72),
                "hsa-mir-21 precursor length (miRBase MI0000077).");
        });
    }

    // CL2 — null/empty input returns null (documented failure mode).
    [Test]
    public void ExtractPreMiRnaFeatures_NullOrEmpty_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ExtractPreMiRnaFeatures(null!), Is.Null, "null ⇒ null.");
            Assert.That(ExtractPreMiRnaFeatures(""), Is.Null, "empty ⇒ null.");
        });
    }

    // CL3 — hsa-mir-21 (real precursor) is classified NATURAL with high probability.
    [Test]
    public void ClassifyPreMiRna_HsaMir21_IsNatural()
    {
        var c = ClassifyPreMiRna(HsaMir21);
        Assert.That(c, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(c!.Value.IsNatural, Is.True, "hsa-mir-21 is a genuine pre-miRNA.");
            Assert.That(c.Value.NaturalProbability, Is.GreaterThan(0.95),
                "Real precursor scores near 1 under the trained model.");
        });
    }

    // CL4 — hsa-let-7a-1 (real precursor) is classified NATURAL with high probability.
    [Test]
    public void ClassifyPreMiRna_HsaLet7a1_IsNatural()
    {
        var c = ClassifyPreMiRna(HsaLet7a1);
        Assert.That(c, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(c!.Value.IsNatural, Is.True, "hsa-let-7a-1 is a genuine pre-miRNA.");
            Assert.That(c.Value.NaturalProbability, Is.GreaterThan(0.95),
                "Real precursor scores near 1 under the trained model.");
        });
    }

    // CL5 — a dinucleotide-shuffled hsa-mir-21 (background) is classified BACKGROUND.
    // The shuffle has IDENTICAL base/dinucleotide composition, so the call must rest on STRUCTURE.
    [Test]
    public void ClassifyPreMiRna_DiShuffledHsaMir21_IsBackground()
    {
        var rng = new Random(999); // fixed seed ⇒ deterministic shuffle
        string shuffled = DinucleotideShuffle(HsaMir21, rng);

        var c = ClassifyPreMiRna(shuffled);
        Assert.That(c, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(c!.Value.IsNatural, Is.False,
                "A di-shuffled hsa-mir-21 (composition-matched background) is NOT a natural pre-miRNA.");
            Assert.That(c.Value.NaturalProbability, Is.LessThan(0.05),
                "Background scores near 0 under the trained model.");
        });
    }

    // CL6 — DinucleotideShuffle preserves the EXACT dinucleotide counts, first base, last base, and
    // length (Altschul & Erickson 1985). This is the property that makes the negative set valid.
    [Test]
    public void DinucleotideShuffle_PreservesDinucleotideCountsFirstLastLength()
    {
        string orig = HsaMir21; // already RNA, upper
        var rng = new Random(7);
        string shuffled = DinucleotideShuffle(orig, rng);

        Assert.Multiple(() =>
        {
            Assert.That(shuffled.Length, Is.EqualTo(orig.Length), "Length preserved.");
            Assert.That(shuffled[0], Is.EqualTo(orig[0]), "First base preserved (Eulerian start).");
            Assert.That(shuffled[^1], Is.EqualTo(orig[^1]), "Last base preserved (Eulerian end).");
            Assert.That(DinucleotideCounts(shuffled), Is.EqualTo(DinucleotideCounts(orig)),
                "Exact dinucleotide (doublet) counts preserved — Altschul-Erickson 1985.");
        });
    }

    // CL7 — DinucleotideShuffle is deterministic for a fixed seed (reproducibility of the run).
    [Test]
    public void DinucleotideShuffle_FixedSeed_IsDeterministic()
    {
        string a = DinucleotideShuffle(HsaMir21, new Random(999));
        string b = DinucleotideShuffle(HsaMir21, new Random(999));
        Assert.That(a, Is.EqualTo(b), "Same seed ⇒ identical shuffle (reproducible training).");
    }

    // CL8 — DinucleotideShuffle null RNG throws; short input is returned unchanged.
    [Test]
    public void DinucleotideShuffle_Guards()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => DinucleotideShuffle("ACGU", null!),
                "null RNG ⇒ ArgumentNullException.");
            Assert.That(DinucleotideShuffle("A", new Random(1)), Is.EqualTo("A"),
                "length<2 returned unchanged (no dinucleotides to permute).");
            Assert.That(DinucleotideShuffle("", new Random(1)), Is.EqualTo(""),
                "empty ⇒ empty.");
        });
    }

    // CL9 — HELD-OUT discrimination metric, reproduced deterministically from the public API.
    // Rebuilds the EXACT training dataset (13 miRBase positives + 4 di-shuffles each, fixed seed,
    // fixed 70/30 split), scores the held-out test split with the bundled trained model, and asserts
    // accuracy and AUC exceed the Bonnet-2004 expectation (≳0.90). The metric is reproducible because
    // the seed and split are fixed.
    [Test]
    public void ClassifyPreMiRna_HeldOutSplit_DiscriminatesNaturalFromBackground()
    {
        var rng = new Random(TrainSeed);
        var examples = new List<(double prob, int label)>();

        Assert.That(MiRBasePositives, Is.EqualTo(PreMiRnaTrainingPositives),
            "The test's miRBase strings are the bundled training positives.");
        // Fixed order: each positive, then its K shuffles (identical to the training run).
        foreach (string pos in MiRBasePositives)
        {
            examples.Add((ClassifyPreMiRna(pos)!.Value.NaturalProbability, 1));
            for (int s = 0; s < ShufflesPerPositive; s++)
            {
                string shuf = DinucleotideShuffle(pos, rng);
                examples.Add((ClassifyPreMiRna(shuf)!.Value.NaturalProbability, 0));
            }
        }

        int n = examples.Count;                    // 65
        int nTrain = (int)Math.Floor(n * 0.70);    // 45
        var test = examples.Skip(nTrain).ToList();  // last 30% = held-out

        int correct = test.Count(t => (t.prob >= 0.5 ? 1 : 0) == t.label);
        double accuracy = (double)correct / test.Count;
        double auc = MannWhitneyAuc(test);

        Assert.Multiple(() =>
        {
            Assert.That(test.Count, Is.EqualTo(n - nTrain), "Held-out split is the last 30%.");
            Assert.That(test.Any(t => t.label == 1) && test.Any(t => t.label == 0), Is.True,
                "Held-out split contains both classes (interleaved order guarantees it).");
            // scikit-learn 1.9.1 LogisticRegression(C = 1/(1e-3·45)) on the same 65 feature rows: held-out
            // accuracy 0.9 (hsa-mir-98 P = 0.4978; one shuffle of hsa-mir-147a P = 0.5001) and
            // roc_auc_score 0.984375 (63 of 64 positive/negative pairs ordered).
            Assert.That(accuracy, Is.EqualTo(0.90).Within(1e-12),
                "Held-out accuracy of the re-derived model (18/20).");
            Assert.That(auc, Is.EqualTo(63.0 / 64.0).Within(1e-12),
                "Held-out AUC = 0.984375 (scikit-learn roc_auc_score on the same scores).");
        });
    }

    // CL10 — ScorePreMiRnaFeatures monotonic sanity: a feature vector matching a real precursor scores
    // higher than one matching its shuffle. Pins that the bundled weights point the right way.
    [Test]
    public void ScorePreMiRnaFeatures_RealScoresHigherThanShuffle()
    {
        var real = ExtractPreMiRnaFeatures(HsaMir21)!.Value;
        var shuf = ExtractPreMiRnaFeatures(DinucleotideShuffle(HsaMir21, new Random(999)))!.Value;

        double pReal = ScorePreMiRnaFeatures(real);
        double pShuf = ScorePreMiRnaFeatures(shuf);

        Assert.That(pReal, Is.GreaterThan(pShuf),
            "Real precursor scores strictly higher than its composition-matched shuffle.");
    }

    // CL11 — threshold parameter is honoured (a threshold above the real score flips the call).
    [Test]
    public void ClassifyPreMiRna_Threshold_IsHonoured()
    {
        var def = ClassifyPreMiRna(HsaMir21, threshold: 0.5)!.Value;
        var strict = ClassifyPreMiRna(HsaMir21, threshold: 1.0)!.Value;

        Assert.Multiple(() =>
        {
            Assert.That(def.IsNatural, Is.True, "At 0.5 a real precursor (P≈1) is natural.");
            Assert.That(strict.IsNatural, Is.False,
                "At threshold 1.0 even P≈0.99999 is below cutoff ⇒ call flips (threshold honoured).");
            Assert.That(strict.NaturalProbability, Is.EqualTo(def.NaturalProbability).Within(1e-12),
                "Probability is threshold-independent; only the boolean call changes.");
        });
    }

    // CL12 — null/empty classification returns null (documented failure mode).
    [Test]
    public void ClassifyPreMiRna_NullOrEmpty_ReturnsNull()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ClassifyPreMiRna(null!), Is.Null, "null ⇒ null.");
            Assert.That(ClassifyPreMiRna(""), Is.Null, "empty ⇒ null.");
        });
    }

    // CL13 — non-complementary sequence (no foldable stem, ΔG° = 0) is classified BACKGROUND.
    // Unlike AssessHairpinByMfe (which returns null for non-hairpins), ClassifyPreMiRna SCORES any
    // sequence: a structureless candidate must score ≈ 0 (well below threshold), never natural.
    // Hand-derived expectation: with ΔG°=0 ⇒ AMFE=0, MFEI=0, %paired=0, z = bias + Σ w_j·(0−mean_j)/std_j
    // ≈ −31.0 ⇒ P ≈ e^z ≈ 8e-16 (computed from the bundled weights), so P < 0.05 and IsNatural=False.
    [Test]
    public void ClassifyPreMiRna_NonComplementary_IsBackground()
    {
        // All-purine: A and G never Watson–Crick/​wobble pair with each other ⇒ no stem ⇒ ΔG° = 0.
        const string allPurine = "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG";
        var c = ClassifyPreMiRna(allPurine);

        Assert.That(c, Is.Not.Null, "ClassifyPreMiRna scores any non-empty sequence (it is not a hairpin gate).");
        Assert.Multiple(() =>
        {
            Assert.That(c!.Value.Features.FreeEnergy, Is.EqualTo(0.0).Within(1e-12),
                "A non-complementary sequence folds to ΔG° = 0 (no base pairs).");
            Assert.That(c.Value.Features.PairedFraction, Is.EqualTo(0.0).Within(1e-12),
                "%paired = 0 for the all-dots MFE structure.");
            Assert.That(c.Value.NaturalProbability, Is.LessThan(0.05),
                "Structureless candidate scores ≈ 0 ⇒ background, never natural.");
            Assert.That(c.Value.IsNatural, Is.False,
                "A sequence with no hairpin structure is not a natural pre-miRNA.");
        });
    }

    // CL14 — di-shuffle discrimination is NOT seed-specific: an INDEPENDENT fresh di-shuffle
    // (different seed than CL5's 999) of hsa-mir-21 also scores BACKGROUND. Confirms the negative
    // call rests on the (degraded) structure of any composition-matched shuffle, not a lucky draw.
    [Test]
    public void ClassifyPreMiRna_IndependentDiShuffle_IsBackground()
    {
        string shuffled = DinucleotideShuffle(HsaMir21, new Random(12345)); // different seed from CL5
        var c = ClassifyPreMiRna(shuffled);

        Assert.That(c, Is.Not.Null);
        Assert.Multiple(() =>
        {
            // Composition is preserved by the di-shuffle ⇒ the call must come from structure.
            Assert.That(DinucleotideCounts(shuffled), Is.EqualTo(DinucleotideCounts(HsaMir21)),
                "Independent shuffle still preserves exact dinucleotide composition (Altschul–Erickson).");
            Assert.That(c!.Value.NaturalProbability, Is.LessThan(0.05),
                "A composition-matched shuffle (different seed) also scores background.");
            Assert.That(c.Value.IsNatural, Is.False,
                "Discrimination is structural, not seed-specific.");
        });
    }

    // CL15 — the bundled model is reproducible from the repository: re-running the documented training
    // (BuildPreMiRnaTrainingSet → standardise → ridge logistic regression by Newton–Raphson) yields the
    // bundled constants. Reference: scikit-learn 1.9.1 LogisticRegression(C = 1/(λ·45), solver
    // newton-cg, tol 1e-14) on the same 45 standardised rows: b0 −4.278229873240855, b = [−0.9994063824008977,
    // 3.304516383745917, 3.63652043260533, −0.8781055255902053, −0.26063466887427356]; numpy mean/std(ddof 0)
    // identical to the printed digits. (Pre-fix: the constants were a 20 000-epoch gradient-ascent fit on
    // features of the pre-B12 folder; 61 of the 65 MFEs have since changed.)
    [Test]
    public void TrainPreMiRnaClassifier_ReproducesBundledModel_AndScikitLearn()
    {
        var m = TrainPreMiRnaClassifier();
        double[] skW = { -0.9994063824008977, 3.304516383745917, 3.63652043260533, -0.8781055255902053, -0.26063466887427356 };
        Assert.Multiple(() =>
        {
            Assert.That(m.Bias, Is.EqualTo(LrBias).Within(1e-12), "bias = bundled constant");
            Assert.That(m.Weights, Is.EqualTo(LrW).Within(1e-12), "weights = bundled constants");
            Assert.That(m.Mean, Is.EqualTo(LrMean).Within(1e-12), "training mean = bundled constants");
            Assert.That(m.Std, Is.EqualTo(LrStd).Within(1e-12), "training population std = bundled constants");
            Assert.That(m.Bias, Is.EqualTo(-4.278229873240855).Within(1e-8), "bias = scikit-learn");
            Assert.That(m.Weights, Is.EqualTo(skW).Within(1e-8), "weights = scikit-learn");
        });
    }

    // CL16 — the fitter returns the maximum-likelihood estimate. Hand-solvable case: x = −1 with y = {0,0,1} and
    // x = +1 with y = {1,1,0}. The model is saturated, so the MLE reproduces the group frequencies:
    // σ(b0 − b) = 1/3 and σ(b0 + b) = 2/3 ⇒ b0 = 0, b = ln 2 (λ = 0). With λ > 0 the score equation
    // mean((y − p)·x) = λ·b holds at the returned b.
    [Test]
    public void FitRidgeLogisticRegression_SaturatedCase_EqualsMle()
    {
        var x = new[] { new[] { -1.0 }, new[] { -1.0 }, new[] { -1.0 }, new[] { 1.0 }, new[] { 1.0 }, new[] { 1.0 } };
        var y = new[] { 0, 0, 1, 1, 1, 0 };
        var (b0, w) = FitRidgeLogisticRegression(x, y, 0.0);
        var (rb0, rw) = FitRidgeLogisticRegression(x, y, 0.1);
        double score = 0;
        for (int i = 0; i < x.Length; i++) score += (y[i] - 1.0 / (1.0 + Math.Exp(-(rb0 + rw[0] * x[i][0])))) * x[i][0];
        Assert.Multiple(() =>
        {
            Assert.That(b0, Is.EqualTo(0.0).Within(1e-12), "symmetric data ⇒ zero intercept");
            Assert.That(w[0], Is.EqualTo(Math.Log(2.0)).Within(1e-12), "MLE slope = ln 2");
            Assert.That(score / x.Length, Is.EqualTo(0.1 * rw[0]).Within(1e-12), "ridge score equation");
            Assert.That(rw[0], Is.LessThan(w[0]), "the penalty shrinks the slope");
        });
    }

    // CL17 — DinucleotideShuffle consumes random draws exactly as uShuffle 1.1.2 (ushuffle.c shuffle2: Wilson
    // arborescence, swap the last-exit edge to the end, permutei, walk). Expected strings: ushuffle.c
    // compiled with set_randfunc() replaying the listed draws (each transcript fully consumed).
    [TestCase("AAUAAGAUAGUCA", new[] { 2, 0, 1, 3, 1, 2, 2, 0, 1 }, "AUAAGAAUAGUCA")]
    [TestCase("acgtacgtnnACGTTTGCAAAGGGCCCNNNNacgt",
        new[] { 1, 1, 6, 0, 5, 3, 7, 1, 4, 0, 0, 1, 0, 3, 3, 1, 0, 6, 4, 2, 2, 0, 1, 2, 4, 3, 2, 0, 1, 2, 1, 2, 0, 0, 3, 0, 0 },
        "ACAGGUUACGUGCGCCGGUNAAACCNNNNNACGUU")]
    [TestCase("UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA",
        new[] { 12, 13, 15, 8, 2, 7, 0, 14, 6, 12, 11, 12, 6, 11, 10, 0, 5, 7, 5, 4, 2, 2, 0, 1, 13, 16, 7, 14, 14, 11, 7, 3, 5,
                1, 8, 3, 5, 2, 4, 2, 2, 0, 5, 7, 4, 0, 7, 7, 4, 3, 4, 0, 0, 2, 1, 0, 9, 0, 8, 3, 8, 1, 6, 4, 1, 2, 1, 1, 1 },
        "UUCAGUGUCUAUACACUUCAGCAUCUGAUGCGGAAUGUCCUGACGGACUUGUGCAUCUGGGACAAGUGGUGA")]
    public void DinucleotideShuffle_ScriptedDraws_EqualUShuffle(string input, int[] draws, string expected)
    {
        var rng = new ScriptedRandom(draws);
        Assert.Multiple(() =>
        {
            Assert.That(DinucleotideShuffle(input, rng), Is.EqualTo(expected), "uShuffle output for the same draws");
            Assert.That(rng.Remaining, Is.EqualTo(0), "same number of draws as uShuffle");
        });
    }

    // CL18 — no silent identity fallback. "A"×500 + "G" + "A"×500 + "C" has 999 valid shuffles (A^i G A^j C,
    // i + j = 1000, i, j ≥ 1). The former rejection loop needed A's last-exit edge to be the single A→C edge
    // (p = 1/1000 per draw) and returned the INPUT after 1000 failed draws: 69 of 200 seeds gave the input back.
    // Uniform sampling returns it with p = 1/999.
    [Test]
    public void DinucleotideShuffle_LowAcceptanceGraph_NotReturnedUnshuffled()
    {
        string s = new string('A', 500) + "G" + new string('A', 500) + "C";
        int unchanged = 0;
        var distinct = new HashSet<string>();
        for (int seed = 0; seed < 200; seed++)
        {
            string t = DinucleotideShuffle(s, new Random(seed));
            Assert.That(DinucleotideCounts(t), Is.EqualTo(DinucleotideCounts(s)));
            distinct.Add(t);
            if (t == s) unchanged++;
        }
        Assert.Multiple(() =>
        {
            Assert.That(unchanged, Is.LessThanOrEqualTo(3), "P(input) = 1/999 per draw (expected 0.2 of 200)");
            Assert.That(distinct.Count, Is.GreaterThan(170), "≈ 200 distinct draws out of 999");
        });
    }

    // CL19 — uniformity (Altschul & Erickson 1985: every sequence with the same doublet counts and ends is
    // equally likely). "ACGUGUCAGGAUCC" has exactly 612 such sequences (Python Eulerian-trail enumeration);
    // 200 000 draws (seed 42) hit all 612 and χ² = 596.49 < χ²(0.999; 611) = 724.75 (scipy).
    [Test]
    public void DinucleotideShuffle_IsUniformOverValidSequences()
    {
        const string s = "ACGUGUCAGGAUCC";
        const int draws = 200_000;
        var counts = new Dictionary<string, int>();
        var rng = new Random(42);
        var reference = DinucleotideCounts(s);
        for (int i = 0; i < draws; i++)
        {
            string t = DinucleotideShuffle(s, rng);
            counts[t] = counts.GetValueOrDefault(t) + 1;
        }
        Assert.That(counts.Keys.All(k => k[0] == s[0] && k[^1] == s[^1]
            && DinucleotideCounts(k).OrderBy(kv => kv.Key).SequenceEqual(reference.OrderBy(kv => kv.Key))), Is.True,
            "every output preserves ends and doublet counts");
        double expected = (double)draws / 612;
        double chi2 = counts.Values.Sum(c => (c - expected) * (c - expected) / expected);
        Assert.Multiple(() =>
        {
            Assert.That(counts.Count, Is.EqualTo(612), "all 612 valid sequences are reachable");
            Assert.That(chi2, Is.LessThan(724.7478), "χ² below the 0.1 % critical value (df 611)");
        });
    }

    // CL20 — randfold (Bonnet et al. 2004; randfold.c dinucleotide mode): P = (1 + #{ΔG°shuffle ≤ ΔG°})/(N + 1).
    // Reference: the same 199 shuffles (Random(11)) folded by the canonical engine and counted by a Python
    // transcription of randfold.c's rule: 20 ⇒ 21/200 = 0.105. ViennaRNA 2.8 RNA.fold on the same shuffles:
    // ΔG° −14.40 (engine −14.45), 22 ⇒ 0.115 (engine-parameter residual of RNA-STRUCT-001).
    // hsa-mir-21 with 99 shuffles (Random(7)): no shuffle folds as low ⇒ 1/100 (Vienna: also 0 of 99).
    [Test]
    public void CalculateRandfoldPValue_CountsAsRandfold()
    {
        var r = CalculateRandfoldPValue(
            "GAUCAGUCAUGACUAGCUAGAUUGACAAGCUUACGAUCAGGACUUAGCAUCGUAGCAUGACUAG", 199, new Random(11))!.Value;
        var mir21 = CalculateRandfoldPValue(HsaMir21, 99, new Random(7))!.Value;
        Assert.Multiple(() =>
        {
            Assert.That(r.LowerOrEqualCount, Is.EqualTo(20));
            Assert.That(r.PValue, Is.EqualTo(21.0 / 200.0).Within(1e-15));
            Assert.That(r.FreeEnergy, Is.EqualTo(-14.45).Within(1e-9));
            Assert.That(r.Iterations, Is.EqualTo(199));
            Assert.That(mir21.PValue, Is.EqualTo(0.01).Within(1e-15), "pre-miRNA: lower MFE than all shuffles");
            Assert.That(mir21.FreeEnergy, Is.EqualTo(-35.13).Within(1e-9));
        });
    }

    [Test]
    public void CalculateRandfoldPValue_Guards()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CalculateRandfoldPValue(HsaMir21, 0, new Random(1))!.Value.PValue, Is.EqualTo(1.0),
                "randfold with 0 randomisations: (ct_eq = 1)/(0 + 1) = 1");
            Assert.That(CalculateRandfoldPValue("", 10, new Random(1)), Is.Null);
            Assert.That(CalculateRandfoldPValue(null!, 10, new Random(1)), Is.Null);
            Assert.Throws<ArgumentOutOfRangeException>(() => CalculateRandfoldPValue(HsaMir21, -1, new Random(1)));
            Assert.Throws<ArgumentNullException>(() => CalculateRandfoldPValue(HsaMir21, 10, null!));
        });
    }

    private sealed class ScriptedRandom : Random
    {
        private readonly Queue<int> _draws;
        public ScriptedRandom(IEnumerable<int> draws) => _draws = new Queue<int>(draws);
        public int Remaining => _draws.Count;
        public override int Next(int maxValue)
        {
            int v = _draws.Dequeue();
            if (v < 0 || v >= maxValue) throw new InvalidOperationException($"scripted draw {v} not in [0,{maxValue})");
            return v;
        }
    }

    private static Dictionary<string, int> DinucleotideCounts(string s)
    {
        var d = new Dictionary<string, int>();
        for (int i = 0; i < s.Length - 1; i++)
        {
            string k = s.Substring(i, 2);
            d[k] = d.GetValueOrDefault(k) + 1;
        }
        return d;
    }

    private static double MannWhitneyAuc(List<(double prob, int label)> data)
    {
        var pos = data.Where(t => t.label == 1).Select(t => t.prob).ToList();
        var neg = data.Where(t => t.label == 0).Select(t => t.prob).ToList();
        if (pos.Count == 0 || neg.Count == 0) return double.NaN;
        double sum = 0;
        foreach (double p in pos)
            foreach (double q in neg)
                sum += p > q ? 1.0 : (p == q ? 0.5 : 0.0);
        return sum / (pos.Count * neg.Count);
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

    #region ScorePreMiRnaFeatures / ClassifyPreMiRna — exact logistic-model killers (MIRNA-CLASSIFY-001)

    // Bundled logistic-regression model (MiRnaAnalyzer §"Trained logistic-regression model"; re-derived by
    // TrainPreMiRnaClassifier and by scikit-learn 1.9.1, see CL15): P(natural) = sigmoid(b0 + Σ b_j·(x_j − mean_j)/std_j) over the 5 model features
    // [FreeEnergy, Amfe, Mfei, GcContent, PairedFraction]. Re-stated here as an INDEPENDENT oracle so the
    // production standardisation / weighting / bias / sigmoid are pinned exactly (not characterised).
    private const double LrBias = -4.278229873244336;
    private static readonly double[] LrW =
        { -0.9994063824017212, 3.3045163837473233, 3.636520432606299, -0.8781055255888673, -0.260634668873017 };
    private static readonly double[] LrMean =
        { -22.93488888888889, 29.25596144559518, 0.6119707188599345, 0.477134229110014, 0.6329947414221273 };
    private static readonly double[] LrStd =
        { 8.568968205027105, 10.686847007183342, 0.20989563586718338, 0.061249038538481494, 0.10670275339810077 };

    private static PreMiRnaFeatures FeaturesFromModelInputs(double fe, double amfe, double mfei, double gc, double paired)
        => new(FreeEnergy: fe, Amfe: amfe, Mfei: mfei, GcContent: gc, PairedFraction: paired,
               StemBasePairs: 20, LoopSize: 8, Length: 80);

    // K1 — At the training MEAN of every model feature, all standardised z_j = 0, so the score is exactly
    // sigmoid(bias). Pins the bias term and the sigmoid link (independent of the weights).
    [Test]
    public void ScorePreMiRnaFeatures_AtFeatureMean_EqualsSigmoidOfBias()
    {
        var f = FeaturesFromModelInputs(LrMean[0], LrMean[1], LrMean[2], LrMean[3], LrMean[4]);

        double expected = 1.0 / (1.0 + Math.Exp(-LrBias));
        Assert.That(ScorePreMiRnaFeatures(f), Is.EqualTo(expected).Within(1e-12),
            "all features at the training mean ⇒ P = sigmoid(b0).");
    }

    // K2 — One standard deviation ABOVE the mean on every feature makes each standardised z_j = 1, so the
    // score is sigmoid(b0 + Σ b_j). Pins the full weight vector AND the per-feature standardisation (std).
    [Test]
    public void ScorePreMiRnaFeatures_AtMeanPlusOneStd_EqualsSigmoidOfBiasPlusWeightSum()
    {
        var f = FeaturesFromModelInputs(
            LrMean[0] + LrStd[0], LrMean[1] + LrStd[1], LrMean[2] + LrStd[2],
            LrMean[3] + LrStd[3], LrMean[4] + LrStd[4]);

        double z = LrBias;
        for (int j = 0; j < LrW.Length; j++) z += LrW[j]; // each standardised feature contributes weight·1
        double expected = 1.0 / (1.0 + Math.Exp(-z));

        Assert.That(ScorePreMiRnaFeatures(f), Is.EqualTo(expected).Within(1e-12),
            "+1σ on every feature ⇒ P = sigmoid(b0 + Σ b_j).");
    }

    // K3 — A single feature offset isolates ONE weight/std pair: only Mfei (index 2) moved +2σ, so
    // z = bias + b2·2. Pins weight[2] and std[2] specifically (a swapped/duplicated weight would diverge).
    [Test]
    public void ScorePreMiRnaFeatures_OnlyMfeiOffset_IsolatesThatFeaturesWeight()
    {
        var f = FeaturesFromModelInputs(LrMean[0], LrMean[1], LrMean[2] + 2.0 * LrStd[2], LrMean[3], LrMean[4]);

        double expected = 1.0 / (1.0 + Math.Exp(-(LrBias + LrW[2] * 2.0)));
        Assert.That(ScorePreMiRnaFeatures(f), Is.EqualTo(expected).Within(1e-12),
            "only Mfei moved +2σ ⇒ P = sigmoid(b0 + b2·2).");
    }

    // K4 — ClassifyPreMiRna's IsNatural is "P(natural) >= threshold" (inclusive). At threshold == the exact
    // probability the call must be TRUE; this pins the boundary against a strict '>' mutant.
    [Test]
    public void ClassifyPreMiRna_ThresholdEqualToProbability_IsNaturalInclusive()
    {
        var c0 = ClassifyPreMiRna(HsaMir21, threshold: 0.0);
        Assert.That(c0, Is.Not.Null);
        double p = c0!.Value.NaturalProbability;

        var atP = ClassifyPreMiRna(HsaMir21, threshold: p);
        Assert.That(atP!.Value.IsNatural, Is.True,
            "P >= threshold is inclusive: threshold == P ⇒ natural (kills the '>' boundary mutant).");

        var aboveP = ClassifyPreMiRna(HsaMir21, threshold: Math.BitIncrement(p));
        Assert.That(aboveP!.Value.IsNatural, Is.False, "threshold just above P ⇒ not natural.");
    }

    #endregion
}
