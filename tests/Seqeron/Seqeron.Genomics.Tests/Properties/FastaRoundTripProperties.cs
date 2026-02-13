namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for FASTA round-trip.
/// Verifies: Parse(ToFasta(entries)) preserves all entries.
///
/// Test Unit: PARSE-FASTA-001 (Property Extension)
/// </summary>
[TestFixture]
[Category("Property")]
[Category("IO")]
public class FastaRoundTripProperties
{
    /// <summary>
    /// Round-trip: Parse(ToFasta(entries)) preserves IDs.
    /// </summary>
    [Test]
    [Category("Property")]
    public void RoundTrip_PreservesIds()
    {
        var entries = new[]
        {
            new FastaEntry("seq1", "First sequence", new DnaSequence("ACGTACGTACGT")),
            new FastaEntry("seq2", "Second sequence", new DnaSequence("TTTGGGCCCAAA")),
            new FastaEntry("seq3", null, new DnaSequence("AAACCCGGGTTT"))
        };

        string fasta = FastaParser.ToFasta(entries);
        var parsed = FastaParser.Parse(fasta).ToList();

        Assert.That(parsed.Count, Is.EqualTo(entries.Length));
        for (int i = 0; i < entries.Length; i++)
            Assert.That(parsed[i].Id, Is.EqualTo(entries[i].Id));
    }

    /// <summary>
    /// Round-trip: Parse(ToFasta(entries)) preserves sequences.
    /// </summary>
    [Test]
    [Category("Property")]
    public void RoundTrip_PreservesSequences()
    {
        var entries = new[]
        {
            new FastaEntry("s1", "desc1", new DnaSequence("ACGTACGTACGT")),
            new FastaEntry("s2", "desc2", new DnaSequence("TTTGGGCCCAAA"))
        };

        string fasta = FastaParser.ToFasta(entries);
        var parsed = FastaParser.Parse(fasta).ToList();

        for (int i = 0; i < entries.Length; i++)
            Assert.That(parsed[i].Sequence.ToString(), Is.EqualTo(entries[i].Sequence.ToString()));
    }

    /// <summary>
    /// Parsed FASTA never returns null entries.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Parse_NeverReturnsNull()
    {
        string fasta = ">id1\nACGT\n>id2\nTTTT\n";
        var entries = FastaParser.Parse(fasta).ToList();

        Assert.That(entries, Has.All.Not.Null);
        Assert.That(entries.Select(e => e.Id), Has.All.Not.Null);
        Assert.That(entries.Select(e => e.Sequence), Has.All.Not.Null);
    }

    /// <summary>
    /// ToFasta output always starts with '>' and contains all IDs.
    /// </summary>
    [Test]
    [Category("Property")]
    public void ToFasta_ContainsAllIds()
    {
        var entries = new[]
        {
            new FastaEntry("alpha", "one", new DnaSequence("ACGT")),
            new FastaEntry("beta", "two", new DnaSequence("TTTT"))
        };

        string fasta = FastaParser.ToFasta(entries);
        Assert.That(fasta, Does.StartWith(">"));
        foreach (var entry in entries)
            Assert.That(fasta, Does.Contain(entry.Id));
    }

    /// <summary>
    /// Double round-trip: ToFasta(Parse(ToFasta(entries))) == ToFasta(entries).
    /// </summary>
    [Test]
    [Category("Property")]
    public void DoubleRoundTrip_IsStable()
    {
        var entries = new[]
        {
            new FastaEntry("x1", "some desc", new DnaSequence("ACGTACGT")),
            new FastaEntry("x2", null, new DnaSequence("GGGGCCCC"))
        };

        string fasta1 = FastaParser.ToFasta(entries);
        var parsed = FastaParser.Parse(fasta1).ToList();
        string fasta2 = FastaParser.ToFasta(parsed);

        Assert.That(fasta2, Is.EqualTo(fasta1));
    }
}
