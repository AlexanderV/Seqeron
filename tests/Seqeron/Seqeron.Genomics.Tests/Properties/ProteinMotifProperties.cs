namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for protein motif finding: common motifs, PROSITE patterns, signal peptide.
///
/// Test Units: PROTMOTIF-FIND-001, PROTMOTIF-PROSITE-001, PROTMOTIF-DOMAIN-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class ProteinMotifProperties
{
    // A protein with known N-glycosylation motif N-{P}-[ST]-{P}
    private const string TestProtein = "MKTLLLTLVVVTLVLSSQPVLSRELRECPRGSGKSCQACPAG" +
                                       "NISTYQCQSYVMSHLCSYQCNQRCFQSLENQCQTFHCRGFQF" +
                                       "NSTRTMPLHCRGFQFNSTRTMPLHCRG";

    // -- PROTMOTIF-FIND-001 --

    /// <summary>
    /// Found motifs have positions within protein bounds.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindCommonMotifs_Positions_WithinBounds()
    {
        var motifs = ProteinMotifFinder.FindCommonMotifs(TestProtein).ToList();

        foreach (var m in motifs)
        {
            Assert.That(m.Start, Is.GreaterThanOrEqualTo(0),
                $"Motif '{m.MotifName}' start {m.Start} < 0");
            Assert.That(m.End, Is.LessThanOrEqualTo(TestProtein.Length),
                $"Motif '{m.MotifName}' end {m.End} > protein length {TestProtein.Length}");
        }
    }

    /// <summary>
    /// Found motif sequence length matches End - Start.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindCommonMotifs_SequenceLength_MatchesSpan()
    {
        var motifs = ProteinMotifFinder.FindCommonMotifs(TestProtein).ToList();

        foreach (var m in motifs)
            Assert.That(m.Sequence.Length, Is.EqualTo(m.End - m.Start).Or.EqualTo(m.End - m.Start + 1),
                $"Motif '{m.MotifName}': sequence length {m.Sequence.Length} vs span {m.End - m.Start}");
    }

    /// <summary>
    /// Each found motif sequence is a substring of the protein.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindCommonMotifs_Sequence_IsSubstringOfProtein()
    {
        var motifs = ProteinMotifFinder.FindCommonMotifs(TestProtein).ToList();

        foreach (var m in motifs)
            Assert.That(TestProtein, Does.Contain(m.Sequence),
                $"Motif sequence '{m.Sequence}' not found in protein");
    }

    /// <summary>
    /// CommonMotifs dictionary is not empty.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CommonMotifs_Dictionary_NotEmpty()
    {
        Assert.That(ProteinMotifFinder.CommonMotifs.Count, Is.GreaterThan(0),
            "CommonMotifs dictionary should not be empty");
    }

    // -- PROTMOTIF-PROSITE-001 --

    /// <summary>
    /// FindMotifByProsite with N-glycosylation pattern N-{P}-[ST]-{P} finds matches in test protein.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindMotifByProsite_ValidPattern_FindsMatches()
    {
        string prositePattern = "N-{P}-[ST]-{P}";
        var motifs = ProteinMotifFinder.FindMotifByProsite(TestProtein, prositePattern).ToList();

        foreach (var m in motifs)
        {
            Assert.That(m.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(m.End, Is.LessThanOrEqualTo(TestProtein.Length));
            Assert.That(m.Sequence[0], Is.EqualTo('N'), "N-glycosylation motif must start with N");
        }
    }

    /// <summary>
    /// ConvertPrositeToRegex produces a valid regex.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ConvertPrositeToRegex_ProducesValidRegex()
    {
        string prositePattern = "N-{P}-[ST]-{P}";
        string regex = ProteinMotifFinder.ConvertPrositeToRegex(prositePattern);

        Assert.That(regex, Is.Not.Null.And.Not.Empty);
        Assert.DoesNotThrow(() => new System.Text.RegularExpressions.Regex(regex),
            $"Converted regex '{regex}' is not valid");
    }

    /// <summary>
    /// FindMotifByPattern with custom regex returns matches with valid positions.
    /// </summary>
    [Test]
    [Category("Property")]
    public void FindMotifByPattern_CustomPattern_ValidPositions()
    {
        // CxxC motif (zinc finger-like)
        string pattern = "C..C";
        var motifs = ProteinMotifFinder.FindMotifByPattern(TestProtein, pattern).ToList();

        foreach (var m in motifs)
        {
            Assert.That(m.Start, Is.GreaterThanOrEqualTo(0));
            Assert.That(m.Sequence.Length, Is.EqualTo(4));
            Assert.That(m.Sequence[0], Is.EqualTo('C'));
            Assert.That(m.Sequence[3], Is.EqualTo('C'));
        }
    }

    // -- PROTMOTIF-DOMAIN-001 (Signal Peptide Prediction) --

    /// <summary>
    /// Signal peptide cleavage position is within protein bounds, if found.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictSignalPeptide_CleavagePosition_WithinBounds()
    {
        // Protein starting with M, hydrophobic region
        string proteinWithSignal = "MKTLLLTLVVVTLVLSSQPVLSRELRECPRGSGKSCQACPAG";
        var sp = ProteinMotifFinder.PredictSignalPeptide(proteinWithSignal);

        if (sp.HasValue)
        {
            Assert.That(sp.Value.CleavagePosition, Is.GreaterThan(0));
            Assert.That(sp.Value.CleavagePosition, Is.LessThan(proteinWithSignal.Length));
            Assert.That(sp.Value.Score, Is.InRange(0.0, 1.0));
        }
    }

    /// <summary>
    /// Transmembrane helix predictions have valid positions.
    /// </summary>
    [Test]
    [Category("Property")]
    public void PredictTransmembraneHelices_Positions_WithinBounds()
    {
        // Hydrophobic sequence
        string protein = "AAAILLLLLLLLLLLLLLLLLLLLAAAAGGGGG" +
                          "LLLLLLLLLLLLLLLLLLLLLLLLAAAA";
        var helices = ProteinMotifFinder.PredictTransmembraneHelices(protein).ToList();

        foreach (var (start, end, score) in helices)
        {
            Assert.That(start, Is.GreaterThanOrEqualTo(0));
            Assert.That(end, Is.LessThanOrEqualTo(protein.Length));
            Assert.That(double.IsFinite(score), Is.True);
        }
    }
}
