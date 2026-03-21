using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for exact/IUPAC/PWM pattern matching.
///
/// Test Units: PAT-EXACT-001, PAT-IUPAC-001, PAT-PWM-001 (Property Extensions)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Matching")]
public class PatternMatchingProperties
{
    private static Arbitrary<string> DnaArbitrary(int minLen = 8) =>
        Gen.Elements('A', 'C', 'G', 'T')
            .ArrayOf()
            .Where(a => a.Length >= minLen)
            .Select(a => new string(a))
            .ToArbitrary();

    #region PAT-EXACT-001: R: positions ∈ [0, len-patLen]; M: substring → count≥1; D: deterministic; P: total ≤ len-patLen+1

    /// <summary>
    /// INV-1: Every occurrence position is within valid bounds [0, seqLen - patLen].
    /// Evidence: A match at position p requires p + patLen ≤ seqLen.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_Positions_WithinBounds_Property()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            string pattern = seq.Substring(0, Math.Min(3, seq.Length));
            var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();

            return positions.All(p => p >= 0 && p + pattern.Length <= seq.Length)
                .Label($"All positions must be in [0, {seq.Length - pattern.Length}]");
        });
    }

    /// <summary>
    /// INV-2: Substring at each found position equals the pattern.
    /// Evidence: Exact match guarantees character-by-character equality.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_FoundSubstrings_EqualPattern_Property()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            string pattern = seq.Substring(0, Math.Min(3, seq.Length));
            var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();

            return positions.All(p => seq.Substring(p, pattern.Length) == pattern)
                .Label("Substring at each position must equal the pattern");
        });
    }

    /// <summary>
    /// INV-3: If the pattern is a known substring, at least one match exists.
    /// Evidence: FindExactMotif must find the pattern at its known position.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_KnownSubstring_AtLeastOneMatch()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            int start = Math.Min(2, seq.Length - 1);
            int len = Math.Min(3, seq.Length - start);
            string pattern = seq.Substring(start, len);
            var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();

            return (positions.Count >= 1)
                .Label($"Pattern '{pattern}' is a substring at {start}, expected ≥1 match, got {positions.Count}");
        });
    }

    /// <summary>
    /// INV-4: Total matches ≤ seqLen - patLen + 1 (upper bound on overlapping matches).
    /// Evidence: Maximum matches occur when every sliding window matches.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_TotalMatches_AtMost_LenMinusPatLenPlus1()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            string pattern = seq.Substring(0, Math.Min(3, seq.Length));
            var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();
            int maxPossible = seq.Length - pattern.Length + 1;

            return (positions.Count <= maxPossible)
                .Label($"Matches={positions.Count} must be ≤ {maxPossible}");
        });
    }

    /// <summary>
    /// INV-5: Exact match is deterministic — same input always yields same positions.
    /// Evidence: FindExactMotif delegates to SuffixTree which is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ExactMatch_IsDeterministic()
    {
        return Prop.ForAll(DnaArbitrary(8), seq =>
        {
            var dna = new DnaSequence(seq);
            string pattern = seq.Substring(0, Math.Min(3, seq.Length));
            var positions1 = MotifFinder.FindExactMotif(dna, pattern).ToList();
            var positions2 = MotifFinder.FindExactMotif(dna, pattern).ToList();

            return positions1.SequenceEqual(positions2)
                .Label("FindExactMotif must be deterministic");
        });
    }

    /// <summary>
    /// A pattern not in the sequence yields no results.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ExactMatch_PatternNotPresent_Empty()
    {
        var dna = new DnaSequence("AAAAAAAAAA");
        var positions = MotifFinder.FindExactMotif(dna, "CCCC").ToList();
        Assert.That(positions, Is.Empty);
    }

    #endregion

    // -- PAT-IUPAC-001 --

    /// <summary>
    /// IUPAC motif 'N' matches any base.
    /// </summary>
    [Test]
    [Category("Property")]
    public void IupacMatch_N_MatchesAnyBase()
    {
        var dna = new DnaSequence("ACGTACGT");
        var matches = MotifFinder.FindDegenerateMotif(dna, "NNNN").ToList();
        Assert.That(matches.Count, Is.EqualTo(5), "'NNNN' on 8-char sequence has exactly 5 positions");
    }

    /// <summary>
    /// IUPAC 'R' matches A or G only.
    /// </summary>
    [Test]
    [Category("Property")]
    public void IupacMatch_R_MatchesPurines()
    {
        string seq = "AAGGCCTT";
        var dna = new DnaSequence(seq);
        var matches = MotifFinder.FindDegenerateMotif(dna, "R").ToList();

        foreach (var m in matches)
        {
            char c = seq[m.Position];
            Assert.That(c, Is.AnyOf('A', 'G'),
                $"IUPAC 'R' matched '{c}' at position {m.Position}");
        }
    }

    /// <summary>
    /// IUPAC match positions are within bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void IupacMatch_Positions_WithinBounds()
    {
        var dna = new DnaSequence("ACGTACGTACGT");
        var matches = MotifFinder.FindDegenerateMotif(dna, "RYN").ToList();

        foreach (var m in matches)
        {
            Assert.That(m.Position, Is.GreaterThanOrEqualTo(0));
            Assert.That(m.Position + m.MatchedSequence.Length, Is.LessThanOrEqualTo(dna.Length));
        }
    }

    // -- PAT-PWM-001 --

    /// <summary>
    /// PWM scores are finite numbers at every position.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PwmScan_Scores_AreFinite()
    {
        var training = new[] { "ACGTAC", "ACGTGC", "ACGAAC" };
        var pwm = MotifFinder.CreatePwm(training);
        var dna = new DnaSequence("ACGTACGTGCACGAACACGTAC");
        var matches = MotifFinder.ScanWithPwm(dna, pwm, threshold: -100).ToList();

        foreach (var m in matches)
            Assert.That(double.IsFinite(m.Score), Is.True,
                $"PWM score at {m.Position} is not finite: {m.Score}");
    }

    /// <summary>
    /// PWM created from identical sequences scores the consensus highest.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PwmScan_ConsensusScoresHighest()
    {
        string motif = "ACGTAC";
        var pwm = MotifFinder.CreatePwm(new[] { motif, motif, motif });
        string seq = "TTTT" + motif + "TTTT";
        var dna = new DnaSequence(seq);
        var matches = MotifFinder.ScanWithPwm(dna, pwm, threshold: -100).ToList();

        if (matches.Count > 0)
        {
            var best = matches.OrderByDescending(m => m.Score).First();
            Assert.That(best.Position, Is.EqualTo(4),
                $"Expected consensus at position 4, got {best.Position}");
        }
    }

    /// <summary>
    /// PWM length matches input sequence length.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Pwm_LengthMatchesInput()
    {
        var training = new[] { "ACGTAC", "ACGTGC", "ACGAAC" };
        var pwm = MotifFinder.CreatePwm(training);
        Assert.That(pwm.Length, Is.EqualTo(6));
    }
}
