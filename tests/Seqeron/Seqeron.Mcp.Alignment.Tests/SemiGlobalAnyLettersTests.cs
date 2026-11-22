using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>
/// R14-2 (B09 F289): semi_global_align / semi_global_align_affine / general_gap_align accept any letters like global_align / local_align (N, IUPAC codes, protein), instead of the strict ACGT
/// <c>DnaSequence.TryCreate</c> ("Invalid DNA in sequence1."). Scores are Biopython 1.88 PairwiseAligner, mode 'global', match 1, mismatch -1, linear gap -1 (affine: open -2, extend -1,
/// i.e. open_gap_score = -3), left_insertion_score = right_insertion_score = 0 (the fitting policy): align(a, b).score.
/// </summary>
[TestFixture]
public class SemiGlobalAnyLettersTests
{
    [TestCase("ACNGT", "TTACAGTTT", 3)]       // N against A is an ordinary mismatch of the scalar model
    [TestCase("ACNGT", "TTACNGTTT", 5)]       // N against N is a match
    [TestCase("RYACGT", "TTRYACGTT", 6)]      // IUPAC letters
    [TestCase("MKTAY", "GGMKTAYGG", 5)]       // protein letters
    [TestCase("MKTWY", "GGMKTAYGG", 3)]
    public void SemiGlobalAlign_AcceptsAnyLetters_ScoreEqualsBiopython(string query, string reference, int expected)
    {
        var r = AlignmentTools.SemiGlobalAlign(query, reference);
        Assert.That(r.Score, Is.EqualTo(expected));
        Assert.That(r.AlignedSequence1.Replace("-", ""), Is.EqualTo(query));
        Assert.That(r.AlignedSequence2.Replace("-", ""), Is.EqualTo(reference));
    }

    [TestCase("ACNGT", "TTACAGTTT", 3)]
    [TestCase("ACNGGT", "TTACNGTTT", 4)]
    public void SemiGlobalAlignAffine_AcceptsAnyLetters_ScoreEqualsBiopython(string query, string reference, int expected)
    {
        var r = AlignmentTools.SemiGlobalAlignAffine(query, reference, gapOpen: -2, gapExtend: -1);
        Assert.That(r.Score, Is.EqualTo(expected));
    }

    [Test]
    public void GeneralGapAlign_AcceptsAnyLetters_AllModes()
    {
        double[] linear = [.. Enumerable.Range(1, 20).Select(k => -(double)k)];
        Assert.That(AlignmentTools.GeneralGapAlign("ACNGT", "TTACAGTTT", linear, "semi_global").Score, Is.EqualTo(3));
        Assert.That(AlignmentTools.GeneralGapAlign("MKTAY", "MKTAY", linear, "global").Score, Is.EqualTo(5));
        Assert.That(AlignmentTools.GeneralGapAlign("xxMKTAYxx", "MKTAY", linear, "local").Score, Is.EqualTo(5));
    }

    [Test]
    public void SemiGlobalAlign_StringPath_EqualsTheLibraryStringOverload_AndLowerCaseIsUppercased()
    {
        var tool = AlignmentTools.SemiGlobalAlign("acngt", "ttacagttt");
        var lib = Seqeron.Genomics.Alignment.SequenceAligner.SemiGlobalAlign("ACNGT", "TTACAGTTT", new Seqeron.Genomics.Infrastructure.ScoringMatrix(1, -1, -2, -1), Seqeron.Genomics.Infrastructure.EndGaps.Fitting);
        Assert.That((tool.AlignedSequence1, tool.AlignedSequence2, tool.Score), Is.EqualTo((lib.AlignedSequence1, lib.AlignedSequence2, lib.Score)));
    }
}
