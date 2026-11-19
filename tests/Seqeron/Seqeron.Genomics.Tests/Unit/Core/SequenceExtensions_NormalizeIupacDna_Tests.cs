namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// <see cref="SequenceExtensions.NormalizeIupacDna"/> (B09 A1-23, F79): the canonical validate-and-normalise step of the IUPAC-aware string
/// overloads (<c>Translator.FindOrfs(string, ...)</c>, <c>GenomicAnalyzer.FindOpenReadingFrames(string, ...)</c>): ASCII upper case, U read as T,
/// the first non-IUPAC symbol rejected with its position.
/// </summary>
[TestFixture]
public class SequenceExtensions_NormalizeIupacDna_Tests
{
    [Test]
    public void Normalize_UpperCasesAndReadsUAsT_EveryIupacCodeAccepted()
    {
        Assert.That(SequenceExtensions.NormalizeIupacDna("acgtu"), Is.EqualTo("ACGTT"));
        Assert.That(SequenceExtensions.NormalizeIupacDna("RYSWKMBDHVNryswkmbdhvn"), Is.EqualTo("RYSWKMBDHVNRYSWKMBDHVN"));
        Assert.That(SequenceExtensions.NormalizeIupacDna(""), Is.EqualTo(""));
    }

    [TestCase("ACGX", 'X', 3)]
    [TestCase("AC-GT", '-', 2)]
    [TestCase("ACGT.", '.', 4)]
    [TestCase("ATGJ", 'J', 3)]
    [TestCase("A C", ' ', 1)]
    public void Normalize_RejectsTheFirstNonIupacSymbol_WithItsPosition(string input, char symbol, int position)
    {
        var ex = Assert.Throws<ArgumentException>(() => SequenceExtensions.NormalizeIupacDna(input));
        Assert.That(ex!.ParamName, Is.EqualTo("sequence"));
        Assert.That(ex.Message, Does.Contain($"'{symbol}' at position {position}"));
    }

    [Test]
    public void Normalize_ReportsTheCallersParameterName_AndNullThrows()
    {
        var ex = Assert.Throws<ArgumentException>(() => SequenceExtensions.NormalizeIupacDna("AQ", "dna"));
        Assert.That(ex!.ParamName, Is.EqualTo("dna"));
        var nullEx = Assert.Throws<ArgumentNullException>(() => SequenceExtensions.NormalizeIupacDna(null!, "dna"));
        Assert.That(nullEx!.ParamName, Is.EqualTo("dna"));
    }

    [Test]
    public void Normalize_NonAsciiLookalikes_AreRejected()
    {
        // U+017F (long s) upper-cases to 'S' under ToUpperInvariant; the ASCII-only fold must not accept it.
        var ex = Assert.Throws<ArgumentException>(() => SequenceExtensions.NormalizeIupacDna("ACſT"));
        Assert.That(ex!.Message, Does.Contain("at position 2"));
    }

    [Test]
    public void Callers_KeepTheirContract_TranslatorAndGenomicAnalyzerStringOverloads()
    {
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindOpenReadingFrames("ATGXTAA", 3, new OrfSearchOptions()).ToList());
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames("augaaauaa", 3, new OrfSearchOptions()).Select(o => o.Sequence).ToList(), Does.Contain("ATGAAATAA"));
        Assert.Throws<ArgumentException>(() => Translator.FindOrfs("ATGXTAA").ToList());
    }
}
