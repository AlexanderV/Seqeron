namespace Seqeron.Genomics.Tests.Unit.Core;

/// <summary>
/// Review 2026-09, B02 duplication sweep: the canonical Core codon splitter / counter
/// (<see cref="Translator.SplitInFrameCodons"/>, <see cref="Translator.CountCodons"/>) that
/// CodonUsageAnalyzer delegates to and that Analysis/Annotation can call (B03 request R16).
/// Contract: EMBOSS cusp / ajCodSetTripletsS (ambiguous triplets skipped without shifting the
/// frame, trailing partial triplet dropped), CodonW ident_codon (U read as T), EMBOSS compseq -frame.
/// </summary>
[TestFixture]
public class Translator_CodonSplitting_Tests
{
    [Test]
    public void SplitInFrameCodons_AmbiguousTripletIsNull_FrameKept()
    {
        Assert.That(Translator.SplitInFrameCodons("augNNNgcuRYUgc"),
            Is.EqualTo(new string?[] { "ATG", null, "GCT", null }));
    }

    [TestCase("ATGATGAAA", 1, new[] { "TGA", "TGA" })]
    [TestCase("ATGATGAAA", 2, new[] { "GAT", "GAA" })]
    [TestCase("ATGATGAAA", 3, new[] { "ATG", "AAA" })]
    public void SplitInFrameCodons_FrameOffset(string sequence, int frame, string[] expected)
    {
        Assert.That(Translator.SplitInFrameCodons(sequence, frame), Is.EqualTo(expected));
    }

    [Test]
    public void SplitInFrameCodons_NegativeFrame_Throws()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() => Translator.SplitInFrameCodons("ATG", -1));
        Assert.Throws<ArgumentOutOfRangeException>(() => Translator.CountCodons("ATG", -1));
    }

    [TestCase(null)]
    [TestCase("")]
    [TestCase("AT")]
    public void CountCodons_NoCompleteCodon_Empty(string? sequence)
    {
        Assert.That(Translator.CountCodons(sequence), Is.Empty);
    }

    [Test]
    public void CountCodons_RnaAndDnaSpellingsEqual()
    {
        Assert.That(Translator.CountCodons("AUGAAAUGA"), Is.EqualTo(Translator.CountCodons("atgaaatga")));
        Assert.That(Translator.CountCodons("AUGAAAUGA"),
            Is.EqualTo(new Dictionary<string, int> { ["ATG"] = 1, ["AAA"] = 1, ["TGA"] = 1 }));
    }

    [Test]
    public void CountCodons_EqualsCodonUsageAnalyzerAndSequenceStatistics_OnRandomInput()
    {
        const string alphabet = "ACGTUacgtuNRY-";
        var random = new Random(28092026);
        for (int n = 0; n < 300; n++)
        {
            var chars = new char[random.Next(0, 90)];
            for (int i = 0; i < chars.Length; i++)
                chars[i] = alphabet[random.Next(alphabet.Length)];
            string seq = new(chars);

            var counts = Translator.CountCodons(seq);
            Assert.That(counts, Is.EqualTo(CodonUsageAnalyzer.CountCodons(seq)), seq);

            for (int frame = 0; frame < 3; frame++)
            {
                var framed = Translator.CountCodons(seq, frame);
                int total = framed.Values.Sum();
                var frequencies = SequenceStatistics.CalculateCodonFrequencies(seq, frame);
                Assert.That(frequencies.Keys, Is.EquivalentTo(framed.Keys), seq);
                foreach (var (codon, count) in framed)
                    Assert.That(frequencies[codon], Is.EqualTo((double)count / total).Within(1e-15), seq);
            }
        }
    }
}
