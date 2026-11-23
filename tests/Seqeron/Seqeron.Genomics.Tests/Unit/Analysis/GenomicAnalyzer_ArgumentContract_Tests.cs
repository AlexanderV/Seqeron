using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// GENOMIC-REPEAT-001 / GENOMIC-COMMON-001 / GENOMIC-ORF-001 / GENOMIC-TANDEM-001 (audit R14-22, B09 F295): ONE argument contract for the genomic-analysis methods, documented in
/// docs/algorithms/Sequence_Comparison/Genomic_Argument_Contract.md and locked here. A <b>minimum length of a reported substring / region / ORF</b> (<c>minLength</c>) that is ≤ 0 means "no
/// lower bound": the library accepts it and behaves exactly as <c>minLength = 1</c> (repeats and regions cannot be empty; every ORF is at least 3 nt long); a <b>structural parameter</b>
/// (tandem unit length, repetitions, unit-length bounds) must be ≥ 1 (≥ 2 for repetitions) and throws <see cref="ArgumentOutOfRangeException"/>. The MCP tools are strict about the first kind: they refuse
/// <c>minLength &lt; 1</c> with <see cref="ArgumentOutOfRangeException"/> (the same message shape in <c>find_repeats</c>, <c>find_common_regions</c>, <c>find_open_reading_frames</c>).
/// </summary>
[TestFixture]
[Category("GENOMIC-REPEAT-001")]
public class GenomicAnalyzer_ArgumentContract_Tests
{
    private static readonly int[] NonPositive = [0, -1, -5, int.MinValue];
    private const string Text = "ATGAAACCCGGGTTTAAATAGATGCCCAAATTTGGGTAACCCATGAAACCCGGGTTTAAATAG";
    private const string Other = "CCCGGGTTTAAATAGATGCCCAAATTTGGGTAACC";

    private static string Show(IEnumerable<RepeatInfo> r) => string.Join(";", r.Select(x => x.Sequence + "@" + string.Join(",", x.Positions)));

    [Test]
    public void Repeats_MinLengthAtMostZero_BehavesAsOne_InEveryOverload()
    {
        var dna = new DnaSequence(Text);
        string expected = Show(GenomicAnalyzer.FindRepeats(dna, 1));
        string expectedMaximal = Show(GenomicAnalyzer.FindMaximalRepeats(dna, 1));
        Assert.That(expected, Is.Not.Empty);
        foreach (int v in NonPositive)
        {
            Assert.That(Show(GenomicAnalyzer.FindRepeats(dna, v)), Is.EqualTo(expected), $"FindRepeats {v}");
            Assert.That(Show(GenomicAnalyzer.FindMaximalRepeats(dna, v)), Is.EqualTo(expectedMaximal), $"FindMaximalRepeats {v}");
            Assert.That(Show(GenomicAnalyzer.FindRepeats(Text, v)), Is.EqualTo(expected), $"FindRepeats(string) {v}");
            Assert.That(Show(GenomicAnalyzer.FindMaximalRepeats(Text, v)), Is.EqualTo(expectedMaximal), $"FindMaximalRepeats(string) {v}");
            Assert.That(GenomicAnalyzer.FindRepeatsBothStrands(dna, v).Select(r => r.Sequence), Is.EqualTo(GenomicAnalyzer.FindRepeatsBothStrands(dna, 1).Select(r => r.Sequence)), $"both strands {v}");
        }
    }

    [Test]
    public void CommonRegions_MinLengthAtMostZero_BehavesAsOne_InEveryOverload()
    {
        var a = new DnaSequence(Text);
        var b = new DnaSequence(Other);
        var expected = GenomicAnalyzer.FindCommonRegions(a, b, 1).Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();
        Assert.That(expected, Is.Not.Empty);
        foreach (int v in NonPositive)
        {
            Assert.That(GenomicAnalyzer.FindCommonRegions(a, b, v).Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)), Is.EqualTo(expected), $"FindCommonRegions {v}");
            Assert.That(GenomicAnalyzer.FindCommonRegions(Text, Other, v).Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)), Is.EqualTo(expected), $"string {v}");
            Assert.That(GenomicAnalyzer.FindCommonRegionSpans(a, b, v), Is.EqualTo(GenomicAnalyzer.FindCommonRegionSpans(a, b, 1)), $"spans {v}");
            Assert.That(GenomicAnalyzer.FindCommonRegions(a, b, v, CommonRegionMode.MaximalExactMatches).Select(r => r.Sequence),
                Is.EqualTo(GenomicAnalyzer.FindCommonRegions(a, b, 1, CommonRegionMode.MaximalExactMatches).Select(r => r.Sequence)), $"MEM {v}");
        }
    }

    [Test]
    public void OpenReadingFrames_MinLengthAtMostZero_AdmitsEveryOrf()
    {
        var dna = new DnaSequence(Text);
        var expected = GenomicAnalyzer.FindOpenReadingFrames(dna, 1).Select(o => (o.Sequence, o.Position, o.Frame, o.IsReverseComplement)).ToList();
        Assert.That(expected, Is.Not.Empty);
        foreach (int v in NonPositive)
        {
            Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, v).Select(o => (o.Sequence, o.Position, o.Frame, o.IsReverseComplement)), Is.EqualTo(expected), $"ORF {v}");
            Assert.That(GenomicAnalyzer.FindOpenReadingFrames(Text, v).Select(o => (o.Sequence, o.Position, o.Frame, o.IsReverseComplement)), Is.EqualTo(expected), $"ORF(string) {v}");
        }

        Assert.That(expected.All(o => o.Sequence.Length >= 3), Is.True, "every ORF is at least one codon: a minimum of 1 and a minimum of 3 admit the same set");
        Assert.That(GenomicAnalyzer.FindOpenReadingFrames(dna, 3).Select(o => o.Sequence), Is.EqualTo(expected.Select(o => o.Sequence)));
    }

    [Test]
    public void TandemRepeats_StructuralParameters_MustBePositive()
    {
        foreach (int v in NonPositive)
            Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindTandemRepeats(Text, v, 2, 6).ToList(), $"minUnitLength {v}");
        foreach (int v in new[] { 1, 0, -1, int.MinValue })
            Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindTandemRepeats(Text, 1, v, 6).ToList(), $"minRepetitions {v}");
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindTandemRepeats(Text, 3, 2, 2).ToList(), "maxUnitLength < minUnitLength");
        Assert.DoesNotThrow(() => GenomicAnalyzer.FindTandemRepeats(Text, 1, 2, 6).ToList());
    }
}
