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
    // -- PAT-EXACT-001 --

    /// <summary>
    /// Every occurrence position is within valid bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ExactMatch_Positions_WithinBounds()
    {
        var dna = new DnaSequence("ACGTACGTACGT");
        var positions = MotifFinder.FindExactMotif(dna, "ACGT").ToList();

        foreach (int p in positions)
        {
            Assert.That(p, Is.GreaterThanOrEqualTo(0));
            Assert.That(p + 4, Is.LessThanOrEqualTo(dna.Length));
        }
    }

    /// <summary>
    /// Substring at each found position equals the pattern.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ExactMatch_FoundSubstrings_EqualPattern()
    {
        string seq = "ACGTACGTACGT";
        var dna = new DnaSequence(seq);
        string pattern = "ACGT";
        var positions = MotifFinder.FindExactMotif(dna, pattern).ToList();

        foreach (int p in positions)
            Assert.That(seq.Substring(p, pattern.Length), Is.EqualTo(pattern));
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
        Assert.That(matches.Count, Is.GreaterThan(0), "'NNNN' should match everywhere");
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
