using FsCheck;
using FsCheck.Fluent;

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
    #region Generators

    // FASTA identifiers are the first whitespace-delimited token of the header, so a valid Id
    // contains no space/tab/newline and no '>'.
    private static Gen<string> IdGen() =>
        from chars in Gen.Elements(
                "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789_.|-".ToCharArray())
            .ArrayOf().Where(a => a.Length is >= 1 and <= 12)
        select new string(chars);

    // The description is everything after the first space; it must not contain a newline (which would
    // break the header line) and is recovered verbatim. We trim edge whitespace to keep a canonical form.
    private static Gen<string?> DescriptionGen() =>
        Gen.OneOf(
            Gen.Constant<string?>(null),
            (from chars in Gen.Elements(
                    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789 _.-".ToCharArray())
                .ArrayOf().Where(a => a.Length >= 1)
             let s = new string(chars).Trim()
             where s.Length >= 1
             select (string?)s));

    // A non-empty DNA sequence — the parser drops records whose sequence body is empty.
    private static Gen<DnaSequence> DnaSeqGen() =>
        from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length is >= 1 and <= 200)
        select new DnaSequence(new string(chars));

    private static Arbitrary<List<FastaEntry>> FastaEntriesArbitrary() =>
        (from n in Gen.Choose(1, 6)
         from entries in (from id in IdGen()
                          from desc in DescriptionGen()
                          from seq in DnaSeqGen()
                          select new FastaEntry(id, desc, seq)).ArrayOf(n)
         select entries.ToList()).ToArbitrary();

    #endregion

    #region PARSE-FASTA-001 — Property-based round-trip

    /// <summary>
    /// RT + P (header + sequence preserved): for any list of FASTA entries with valid identifiers,
    /// Parse(ToFasta(entries)) reproduces the entries exactly — count, Id, Description and Sequence.
    /// This is the defining inverse-pair law of a serializer/parser. Source: NCBI/Wikipedia FASTA.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property RoundTrip_PreservesAllEntries()
    {
        return Prop.ForAll(FastaEntriesArbitrary(), entries =>
        {
            string fasta = FastaParser.ToFasta(entries);
            var parsed = FastaParser.Parse(fasta).ToList();

            if (parsed.Count != entries.Count)
                return false.Label($"entry count {entries.Count} → {parsed.Count}");

            for (int i = 0; i < entries.Count; i++)
            {
                if (parsed[i].Id != entries[i].Id)
                    return false.Label($"Id[{i}] '{entries[i].Id}' → '{parsed[i].Id}'");
                if (parsed[i].Description != entries[i].Description)
                    return false.Label($"Description[{i}] '{entries[i].Description}' → '{parsed[i].Description}'");
                if (parsed[i].Sequence.ToString() != entries[i].Sequence.ToString())
                    return false.Label($"Sequence[{i}] mismatch for Id '{entries[i].Id}'");
            }
            return true.ToProperty();
        });
    }

    /// <summary>
    /// D (determinism): ToFasta is a pure function — serializing the same entries twice yields
    /// byte-identical output, and re-parsing is stable.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property ToFasta_IsDeterministic()
    {
        return Prop.ForAll(FastaEntriesArbitrary(), entries =>
        {
            string a = FastaParser.ToFasta(entries);
            string b = FastaParser.ToFasta(entries);
            return (a == b).Label("ToFasta must be deterministic for identical input");
        });
    }

    #endregion

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
