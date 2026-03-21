using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for repeat finder algorithms.
/// Verifies invariants for palindromes, inverted repeats, direct repeats, and microsatellites.
///
/// Test Units: REP-STR-001, REP-PALIN-001, REP-TANDEM-001, REP-INV-001, REP-DIRECT-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class RepeatFinderProperties
{
    // Sequences containing known repeats for reliable testing
    private const string PalindromeSequence = "ACGTGAATTCACGTACGTGATATCACGT";
    private const string MicrosatelliteSequence = "ACACACACACACGTGTGTGTGTAAAAAAAAAA";

    private static Arbitrary<string> DnaArbitrary(int minLen = 20) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    #region REP-STR-001: R: positions ≥ 0; M: lower minRepeats → ≥ results; R: repeat count ≥ minRepeats; P: unit len in range

    /// <summary>
    /// INV-1: Microsatellite positions are non-negative.
    /// Evidence: Positions are indices into a 0-based string.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Microsatellite_Positions_NonNegative()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var results = RepeatFinder.FindMicrosatellites(seq, 1, 4, minRepeats: 3).ToList();
            return results.All(r => r.Position >= 0)
                .Label("All positions must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-2: Microsatellite positions + total length do not exceed sequence length.
    /// Evidence: A repeat at position p with length L requires p + L ≤ seqLen.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Microsatellite_Positions_WithinBounds()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var results = RepeatFinder.FindMicrosatellites(seq, 1, 4, minRepeats: 3).ToList();
            return results.All(r => r.Position + r.TotalLength <= seq.Length)
                .Label("Position + TotalLength must not exceed sequence length");
        });
    }

    /// <summary>
    /// INV-3: Repeat count is at least the requested minimum.
    /// Evidence: FindMicrosatellites filters by minRepeats parameter.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Microsatellite_RepeatCount_MeetsMinimum()
    {
        int minRepeats = 3;
        var results = RepeatFinder.FindMicrosatellites(MicrosatelliteSequence, 1, 4, minRepeats).ToList();

        foreach (var r in results)
            Assert.That(r.RepeatCount, Is.GreaterThanOrEqualTo(minRepeats),
                $"Repeat count {r.RepeatCount} must be ≥ {minRepeats}");
    }

    /// <summary>
    /// INV-4: Repeat unit length is within the requested [minUnit, maxUnit] range.
    /// Evidence: FindMicrosatellites scans only for units in the specified range.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Microsatellite_UnitLength_InRange()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            int minUnit = 1, maxUnit = 4;
            var results = RepeatFinder.FindMicrosatellites(seq, minUnit, maxUnit, minRepeats: 3).ToList();
            return results.All(r => r.RepeatUnit.Length >= minUnit && r.RepeatUnit.Length <= maxUnit)
                .Label($"Unit length must be in [{minUnit}, {maxUnit}]");
        });
    }

    /// <summary>
    /// INV-5: Lower minRepeats yields more or equal results (monotonicity).
    /// Evidence: Relaxing the minimum threshold expands the result set.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Microsatellite_LowerMinRepeats_MoreOrEqualResults()
    {
        var resultsStrict = RepeatFinder.FindMicrosatellites(
            MicrosatelliteSequence, 1, 4, minRepeats: 5).ToList();
        var resultsRelaxed = RepeatFinder.FindMicrosatellites(
            MicrosatelliteSequence, 1, 4, minRepeats: 3).ToList();

        Assert.That(resultsRelaxed.Count, Is.GreaterThanOrEqualTo(resultsStrict.Count),
            $"minRepeats=3 → {resultsRelaxed.Count} must be ≥ minRepeats=5 → {resultsStrict.Count}");
    }

    /// <summary>
    /// INV-6: Microsatellite detection is deterministic.
    /// Evidence: FindMicrosatellites is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Microsatellite_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var results1 = RepeatFinder.FindMicrosatellites(seq, 1, 4, minRepeats: 3).ToList();
            var results2 = RepeatFinder.FindMicrosatellites(seq, 1, 4, minRepeats: 3).ToList();
            return results1.SequenceEqual(results2)
                .Label("FindMicrosatellites must be deterministic");
        });
    }

    #endregion

    #region REP-PALIN-001: P: palindrome = revcomp of self; R: len ∈ [minLen, maxLen]; R: positions valid

    /// <summary>
    /// INV-1: Palindrome positions are within sequence bounds.
    /// Evidence: A palindrome at position p with length L requires 0 ≤ p and p + L ≤ seqLen.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Palindrome_Positions_WithinBounds()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var palindromes = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12).ToList();
            return palindromes.All(p => p.Position >= 0 && p.Position + p.Length <= seq.Length)
                .Label("All palindrome positions must be within sequence bounds");
        });
    }

    /// <summary>
    /// INV-2: Each palindrome sequence equals its own reverse complement.
    /// Evidence: DNA palindrome definition — reads the same on both strands 5'→3'.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Palindrome_EqualsReverseComplement()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var palindromes = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12).ToList();
            return palindromes.All(p =>
            {
                string palSeq = seq.Substring(p.Position, p.Length);
                return palSeq == ReverseComplement(palSeq);
            }).Label("Each palindrome must equal its reverse complement");
        });
    }

    /// <summary>
    /// INV-3: Palindrome lengths are within the requested [minLen, maxLen] range.
    /// Evidence: FindPalindromes restricts output to the specified length range.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Palindrome_LengthInRange()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            int minLen = 4, maxLen = 12;
            var palindromes = RepeatFinder.FindPalindromes(seq, minLen, maxLen).ToList();
            return palindromes.All(p => p.Length >= minLen && p.Length <= maxLen)
                .Label($"Palindrome length must be in [{minLen}, {maxLen}]");
        });
    }

    /// <summary>
    /// INV-4: Palindrome lengths are even (DNA palindromes pair symmetrically).
    /// Evidence: Each base on one half pairs with a complement on the other half.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Palindrome_LengthsAreEven()
    {
        var palindromes = RepeatFinder.FindPalindromes(PalindromeSequence, minLength: 4, maxLength: 12).ToList();

        foreach (var p in palindromes)
            Assert.That(p.Length % 2, Is.EqualTo(0),
                $"Palindrome length {p.Length} must be even");
    }

    /// <summary>
    /// INV-5: Palindrome detection is deterministic.
    /// Evidence: FindPalindromes is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Palindrome_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(20), seq =>
        {
            var results1 = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12).ToList();
            var results2 = RepeatFinder.FindPalindromes(seq, minLength: 4, maxLength: 12).ToList();
            return results1.SequenceEqual(results2)
                .Label("FindPalindromes must be deterministic");
        });
    }

    #endregion

    #region REP-TANDEM-001: R: repeat count ≥ minReps; M: wider unit range → ≥ results; R: positions valid; D: deterministic

    /// <summary>
    /// INV-1: Tandem repeat summary total repeats is non-negative.
    /// Evidence: Count of microsatellites found cannot be negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property TandemSummary_TotalRepeats_NonNegative()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var summary = RepeatFinder.GetTandemRepeatSummary(new DnaSequence(seq), minRepeats: 3);
            return (summary.TotalRepeats >= 0)
                .Label($"TotalRepeats={summary.TotalRepeats} must be ≥ 0");
        });
    }

    /// <summary>
    /// INV-2: Tandem repeat bases percentage is in [0, 100].
    /// Evidence: PercentageOfSequence = totalBases / seqLen * 100.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property TandemSummary_Percentage_InRange()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var summary = RepeatFinder.GetTandemRepeatSummary(new DnaSequence(seq), minRepeats: 3);
            return (summary.PercentageOfSequence >= 0 && summary.PercentageOfSequence <= 100)
                .Label($"PercentageOfSequence={summary.PercentageOfSequence} must be in [0, 100]");
        });
    }

    /// <summary>
    /// INV-3: Sub-type counts sum to total repeats.
    /// Evidence: TotalRepeats = Mono + Di + Tri + Tetra + remaining types.
    /// </summary>
    [Test]
    [Category("Property")]
    public void TandemSummary_SubTypeCounts_SumToTotal()
    {
        var summary = RepeatFinder.GetTandemRepeatSummary(
            new DnaSequence(MicrosatelliteSequence), minRepeats: 3);

        int subSum = summary.MononucleotideRepeats + summary.DinucleotideRepeats +
                     summary.TrinucleotideRepeats + summary.TetranucleotideRepeats;

        Assert.That(subSum, Is.LessThanOrEqualTo(summary.TotalRepeats),
            "Mono+Di+Tri+Tetra sub-types must be ≤ TotalRepeats (penta/hexa/complex may exist)");
    }

    /// <summary>
    /// INV-4: Lower minRepeats yields more or equal tandem repeats (monotonicity).
    /// Evidence: Relaxing the repeat threshold cannot remove existing results.
    /// </summary>
    [Test]
    [Category("Property")]
    public void TandemSummary_LowerMinRepeats_Monotonic()
    {
        var dna = new DnaSequence(MicrosatelliteSequence);
        var strict = RepeatFinder.GetTandemRepeatSummary(dna, minRepeats: 5);
        var relaxed = RepeatFinder.GetTandemRepeatSummary(dna, minRepeats: 3);

        Assert.That(relaxed.TotalRepeats, Is.GreaterThanOrEqualTo(strict.TotalRepeats),
            $"minRepeats=3 → {relaxed.TotalRepeats} must be ≥ minRepeats=5 → {strict.TotalRepeats}");
    }

    /// <summary>
    /// INV-5: Tandem repeat summary is deterministic.
    /// Evidence: GetTandemRepeatSummary is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property TandemSummary_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var dna = new DnaSequence(seq);
            var r1 = RepeatFinder.GetTandemRepeatSummary(dna, minRepeats: 3);
            var r2 = RepeatFinder.GetTandemRepeatSummary(dna, minRepeats: 3);
            return (r1.TotalRepeats == r2.TotalRepeats &&
                    r1.TotalRepeatBases == r2.TotalRepeatBases)
                .Label("GetTandemRepeatSummary must be deterministic");
        });
    }

    #endregion

    #region REP-INV-001: P: right arm = revcomp(left arm); R: positions valid; R: arm len ≥ minLen; D: deterministic

    /// <summary>
    /// INV-1: Inverted repeat positions are within sequence bounds.
    /// Evidence: LeftArmStart ≥ 0 and RightArmStart + ArmLength ≤ seqLen.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvertedRepeat_Positions_WithinBounds_Property()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var repeats = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4).ToList();
            return repeats.All(r =>
                r.LeftArmStart >= 0 &&
                r.RightArmStart >= 0 &&
                r.RightArmStart + r.ArmLength <= seq.Length)
                .Label("All inverted repeat positions must be within sequence bounds");
        });
    }

    /// <summary>
    /// INV-2: Right arm equals reverse complement of left arm.
    /// Evidence: Definition of inverted repeat — arms are reverse complementary.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvertedRepeat_RightArm_IsRevCompOfLeftArm()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var repeats = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4).ToList();
            return repeats.All(r => r.RightArm == ReverseComplement(r.LeftArm))
                .Label("Right arm must equal reverse complement of left arm");
        });
    }

    /// <summary>
    /// INV-3: Arm length meets the requested minimum.
    /// Evidence: FindInvertedRepeats filters by minArmLength parameter.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvertedRepeat_ArmLength_MeetsMinimum()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            int minArm = 4;
            var repeats = RepeatFinder.FindInvertedRepeats(seq, minArmLength: minArm).ToList();
            return repeats.All(r => r.ArmLength >= minArm)
                .Label($"Arm length must be ≥ {minArm}");
        });
    }

    /// <summary>
    /// INV-4: Loop length is within the requested [minLoop, maxLoop] range.
    /// Evidence: FindInvertedRepeats constrains loop size.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvertedRepeat_LoopLength_InRange()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            int minLoop = 3, maxLoop = 50;
            var repeats = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4,
                maxLoopLength: maxLoop, minLoopLength: minLoop).ToList();
            return repeats.All(r => r.LoopLength >= minLoop && r.LoopLength <= maxLoop)
                .Label($"Loop length must be in [{minLoop}, {maxLoop}]");
        });
    }

    /// <summary>
    /// INV-5: Inverted repeat detection is deterministic.
    /// Evidence: FindInvertedRepeats is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property InvertedRepeat_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var r1 = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4).ToList();
            var r2 = RepeatFinder.FindInvertedRepeats(seq, minArmLength: 4).ToList();
            return r1.SequenceEqual(r2)
                .Label("FindInvertedRepeats must be deterministic");
        });
    }

    #endregion

    #region REP-DIRECT-001: R: positions valid; M: lower minLen → ≥ results; P: two copies identical; D: deterministic

    /// <summary>
    /// INV-1: Direct repeat positions are within sequence bounds.
    /// Evidence: Both FirstPosition and SecondPosition + Length ≤ seqLen.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DirectRepeat_Positions_WithinBounds_Property()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var repeats = RepeatFinder.FindDirectRepeats(seq, minLength: 5).ToList();
            return repeats.All(r =>
                r.FirstPosition >= 0 &&
                r.SecondPosition >= 0 &&
                r.FirstPosition + r.Length <= seq.Length &&
                r.SecondPosition + r.Length <= seq.Length)
                .Label("All direct repeat positions must be within sequence bounds");
        });
    }

    /// <summary>
    /// INV-2: Both copies of a direct repeat are identical.
    /// Evidence: A direct repeat has two identical subsequences at different positions.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DirectRepeat_BothCopies_AreIdentical()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var repeats = RepeatFinder.FindDirectRepeats(seq, minLength: 5).ToList();
            return repeats.All(r =>
            {
                string copy1 = seq.Substring(r.FirstPosition, r.Length);
                string copy2 = seq.Substring(r.SecondPosition, r.Length);
                return copy1 == copy2;
            }).Label("Both copies of a direct repeat must be identical");
        });
    }

    /// <summary>
    /// INV-3: Repeat length meets the requested minimum.
    /// Evidence: FindDirectRepeats filters by minLength parameter.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DirectRepeat_Length_MeetsMinimum()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            int minLen = 5;
            var repeats = RepeatFinder.FindDirectRepeats(seq, minLength: minLen).ToList();
            return repeats.All(r => r.Length >= minLen)
                .Label($"Repeat length must be ≥ {minLen}");
        });
    }

    /// <summary>
    /// INV-4: Lower minLength yields more or equal direct repeats (monotonicity).
    /// Evidence: Relaxing the length threshold expands the result set.
    /// </summary>
    [Test]
    [Category("Property")]
    public void DirectRepeat_LowerMinLength_Monotonic()
    {
        string seq = "ACGTACGTTTTTTTTTTTACGTACGTTTTTTT";
        var strict = RepeatFinder.FindDirectRepeats(seq, minLength: 8).ToList();
        var relaxed = RepeatFinder.FindDirectRepeats(seq, minLength: 5).ToList();

        Assert.That(relaxed.Count, Is.GreaterThanOrEqualTo(strict.Count),
            $"minLength=5 → {relaxed.Count} must be ≥ minLength=8 → {strict.Count}");
    }

    /// <summary>
    /// INV-5: Direct repeat detection is deterministic.
    /// Evidence: FindDirectRepeats is a pure function.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property DirectRepeat_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(30), seq =>
        {
            var r1 = RepeatFinder.FindDirectRepeats(seq, minLength: 5).ToList();
            var r2 = RepeatFinder.FindDirectRepeats(seq, minLength: 5).ToList();
            return r1.SequenceEqual(r2)
                .Label("FindDirectRepeats must be deterministic");
        });
    }

    #endregion

    private static string ReverseComplement(string dna)
    {
        var comp = dna.Select(c => c switch
        {
            'A' => 'T',
            'T' => 'A',
            'G' => 'C',
            'C' => 'G',
            _ => c
        }).Reverse().ToArray();
        return new string(comp);
    }
}
