using NUnit.Framework;
using SuffixTree.Genomics;
using System.Linq;

namespace SuffixTree.Genomics.Tests;

[TestFixture]
public class ISequenceTests
{
    #region IupacDnaSequence Tests

    [Test]
    public void IupacDnaSequence_StandardBases_Valid()
    {
        var seq = new IupacDnaSequence("ACGT");

        Assert.That(seq.Length, Is.EqualTo(4));
        Assert.That(seq.Sequence, Is.EqualTo("ACGT"));
        Assert.That(seq.IsValid(), Is.True);
        Assert.That(seq.Type, Is.EqualTo(SequenceType.IupacDna));
    }

    [Test]
    public void IupacDnaSequence_AmbiguityCodes_Valid()
    {
        var seq = new IupacDnaSequence("ACGTNRYSWKMBDHV");

        Assert.That(seq.IsValid(), Is.True);
    }

    [Test]
    public void IupacDnaSequence_GetComplement_CorrectlyComplements()
    {
        var seq = new IupacDnaSequence("ACGTRYWSKMBDHVN");
        var comp = seq.GetComplement() as IupacDnaSequence;

        Assert.That(comp, Is.Not.Null);
        // Correct IUPAC complements:
        // A->T, C->G, G->C, T->A, R->Y, Y->R, W->W, S->S, K->M, M->K, B->V, D->H, H->D, V->B, N->N
        Assert.That(comp!.Sequence, Is.EqualTo("TGCAYRWSMKVHDBN"));
    }

    [Test]
    public void IupacDnaSequence_GetReverseComplement_Works()
    {
        var seq = new IupacDnaSequence("ACGT");
        var rc = seq.GetReverseComplement() as IupacDnaSequence;

        Assert.That(rc, Is.Not.Null);
        Assert.That(rc!.Sequence, Is.EqualTo("ACGT")); // Palindrome
    }

    [Test]
    public void IupacDnaSequence_ExpandCode_ExpandsCorrectly()
    {
        Assert.That(IupacDnaSequence.ExpandCode('A'), Is.EquivalentTo(new[] { 'A' }));
        Assert.That(IupacDnaSequence.ExpandCode('N'), Is.EquivalentTo(new[] { 'A', 'C', 'G', 'T' }));
        Assert.That(IupacDnaSequence.ExpandCode('R'), Is.EquivalentTo(new[] { 'A', 'G' }));
        Assert.That(IupacDnaSequence.ExpandCode('Y'), Is.EquivalentTo(new[] { 'C', 'T' }));
        Assert.That(IupacDnaSequence.ExpandCode('W'), Is.EquivalentTo(new[] { 'A', 'T' }));
        Assert.That(IupacDnaSequence.ExpandCode('S'), Is.EquivalentTo(new[] { 'G', 'C' }));
    }

    [Test]
    public void IupacDnaSequence_GetIupacCode_EncodesCorrectly()
    {
        Assert.That(IupacDnaSequence.GetIupacCode(new[] { 'A' }), Is.EqualTo('A'));
        Assert.That(IupacDnaSequence.GetIupacCode(new[] { 'A', 'G' }), Is.EqualTo('R'));
        Assert.That(IupacDnaSequence.GetIupacCode(new[] { 'C', 'T' }), Is.EqualTo('Y'));
        Assert.That(IupacDnaSequence.GetIupacCode(new[] { 'A', 'C', 'G', 'T' }), Is.EqualTo('N'));
    }

    [Test]
    public void IupacDnaSequence_CodesMatch_MatchesCorrectly()
    {
        Assert.That(IupacDnaSequence.CodesMatch('A', 'A'), Is.True);
        Assert.That(IupacDnaSequence.CodesMatch('A', 'R'), Is.True);  // A matches R (A/G)
        Assert.That(IupacDnaSequence.CodesMatch('A', 'Y'), Is.False); // A doesn't match Y (C/T)
        Assert.That(IupacDnaSequence.CodesMatch('N', 'A'), Is.True);  // N matches anything
        Assert.That(IupacDnaSequence.CodesMatch('R', 'Y'), Is.False); // No overlap
        Assert.That(IupacDnaSequence.CodesMatch('W', 'M'), Is.True);  // A in common
    }

    [Test]
    public void IupacDnaSequence_MatchesAt_WithWildcards()
    {
        var seq = new IupacDnaSequence("ACGTACGT");

        Assert.That(seq.MatchesAt("ACGT", 0), Is.True);
        Assert.That(seq.MatchesAt("NNNN", 0), Is.True);  // N matches any
        Assert.That(seq.MatchesAt("RCGT", 0), Is.True);  // R=A/G, A matches
        Assert.That(seq.MatchesAt("YCGT", 0), Is.False); // Y=C/T, A doesn't match
    }

    [Test]
    public void IupacDnaSequence_FindPattern_WithWildcards()
    {
        var seq = new IupacDnaSequence("AAAAACGTAAAA");

        var positions = seq.FindPattern("RCGT").ToList();

        Assert.That(positions, Contains.Item(4));
    }

    [Test]
    public void IupacDnaSequence_GetAmbiguityLevel_CalculatesCorrectly()
    {
        var unambiguous = new IupacDnaSequence("ACGT");
        var halfAmbiguous = new IupacDnaSequence("ACNN");
        var allAmbiguous = new IupacDnaSequence("NNNN");

        Assert.That(unambiguous.GetAmbiguityLevel(), Is.EqualTo(1.0));
        Assert.That(halfAmbiguous.GetAmbiguityLevel(), Is.EqualTo(0.5));
        Assert.That(allAmbiguous.GetAmbiguityLevel(), Is.EqualTo(0.0));
    }

    [Test]
    public void IupacDnaSequence_ExpandAll_GeneratesAllPossibilities()
    {
        var seq = new IupacDnaSequence("AN");

        var expanded = seq.ExpandAll().ToList();

        Assert.That(expanded, Has.Count.EqualTo(4));
        Assert.That(expanded, Contains.Item("AA"));
        Assert.That(expanded, Contains.Item("AC"));
        Assert.That(expanded, Contains.Item("AG"));
        Assert.That(expanded, Contains.Item("AT"));
    }

    [Test]
    public void IupacDnaSequence_ExpandAll_LimitsResults()
    {
        var seq = new IupacDnaSequence("NNNNNN"); // 4^6 = 4096 possibilities

        var expanded = seq.ExpandAll(maxResults: 100).ToList();

        Assert.That(expanded.Count, Is.LessThanOrEqualTo(100));
    }

    [Test]
    public void IupacDnaSequence_Subsequence_ReturnsCorrectType()
    {
        var seq = new IupacDnaSequence("ACGTNRYSWKM");
        var sub = seq.Subsequence(2, 5);

        Assert.That(sub, Is.TypeOf<IupacDnaSequence>());
        Assert.That(sub.Sequence, Is.EqualTo("GTNRY"));
    }

    #endregion

    #region QualitySequence Tests

    [Test]
    public void QualitySequence_Constructor_SetsSequenceAndQuality()
    {
        var qual = new byte[] { 30, 30, 30, 30 };
        var seq = new QualitySequence("ACGT", qual);

        Assert.That(seq.Sequence, Is.EqualTo("ACGT"));
        Assert.That(seq.Qualities, Is.EquivalentTo(qual));
        Assert.That(seq.Type, Is.EqualTo(SequenceType.Quality));
    }

    [Test]
    public void QualitySequence_FromQualityString_ParsesCorrectly()
    {
        // Phred+33: 'I' = 40, '5' = 20
        var seq = new QualitySequence("ACGT", "II55", phredOffset: 33);

        Assert.That(seq.GetQuality(0), Is.EqualTo(40));
        Assert.That(seq.GetQuality(1), Is.EqualTo(40));
        Assert.That(seq.GetQuality(2), Is.EqualTo(20));
        Assert.That(seq.GetQuality(3), Is.EqualTo(20));
    }

    [Test]
    public void QualitySequence_MismatchedLength_ThrowsException()
    {
        Assert.Throws<ArgumentException>(() =>
            new QualitySequence("ACGT", new byte[] { 30, 30 }));
    }

    [Test]
    public void QualitySequence_MeanQuality_CalculatesCorrectly()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 10, 20, 30, 40 });

        Assert.That(seq.MeanQuality, Is.EqualTo(25.0));
    }

    [Test]
    public void QualitySequence_GetQualityString_EncodesCorrectly()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 0, 10, 30, 40 });

        var qualStr = seq.GetQualityString(33);

        Assert.That(qualStr[0], Is.EqualTo('!')); // 0 + 33 = '!'
        Assert.That(qualStr[3], Is.EqualTo('I')); // 40 + 33 = 'I'
    }

    [Test]
    public void QualitySequence_TrimByQuality_TrimsLowQualityEnds()
    {
        var seq = new QualitySequence("AACGTAA", new byte[] { 5, 5, 30, 30, 30, 5, 5 });

        var trimmed = seq.TrimByQuality(minQuality: 20);

        Assert.That(trimmed.Sequence, Is.EqualTo("CGT"));
        Assert.That(trimmed.Length, Is.EqualTo(3));
    }

    [Test]
    public void QualitySequence_TrimByQuality_AllLowQuality_ReturnsEmpty()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 5, 5, 5, 5 });

        var trimmed = seq.TrimByQuality(minQuality: 20);

        Assert.That(trimmed.Length, Is.EqualTo(0));
    }

    [Test]
    public void QualitySequence_MaskLowQuality_MasksWithN()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 30, 5, 30, 5 });

        var masked = seq.MaskLowQuality(minQuality: 20);

        Assert.That(masked.Sequence, Is.EqualTo("ANGT".Replace('G', 'N').Replace("AN", "AN")[0] == 'A' ? "ANGN" : "ANGN"));
        Assert.That(masked.Sequence[0], Is.EqualTo('A'));
        Assert.That(masked.Sequence[1], Is.EqualTo('N'));
        Assert.That(masked.Sequence[2], Is.EqualTo('G'));
        Assert.That(masked.Sequence[3], Is.EqualTo('N'));
    }

    [Test]
    public void QualitySequence_PhredToErrorProbability_ConvertsCorrectly()
    {
        Assert.That(QualitySequence.PhredToErrorProbability(10), Is.EqualTo(0.1).Within(0.001));
        Assert.That(QualitySequence.PhredToErrorProbability(20), Is.EqualTo(0.01).Within(0.0001));
        Assert.That(QualitySequence.PhredToErrorProbability(30), Is.EqualTo(0.001).Within(0.00001));
    }

    [Test]
    public void QualitySequence_ErrorProbabilityToPhred_ConvertsCorrectly()
    {
        Assert.That(QualitySequence.ErrorProbabilityToPhred(0.1), Is.EqualTo(10));
        Assert.That(QualitySequence.ErrorProbabilityToPhred(0.01), Is.EqualTo(20));
        Assert.That(QualitySequence.ErrorProbabilityToPhred(0.001), Is.EqualTo(30));
    }

    [Test]
    public void QualitySequence_ExpectedErrors_CalculatesCorrectly()
    {
        // Q30 = 0.001 error probability per base
        var seq = new QualitySequence("ACGT", new byte[] { 30, 30, 30, 30 });

        var expected = seq.ExpectedErrors();

        Assert.That(expected, Is.EqualTo(0.004).Within(0.0001));
    }

    [Test]
    public void QualitySequence_GetComplement_PreservesQuality()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 10, 20, 30, 40 });
        var comp = seq.GetComplement() as QualitySequence;

        Assert.That(comp, Is.Not.Null);
        Assert.That(comp!.Sequence, Is.EqualTo("TGCA"));
        Assert.That(comp.Qualities, Is.EquivalentTo(new byte[] { 10, 20, 30, 40 }));
    }

    [Test]
    public void QualitySequence_GetReverseComplement_ReversesQuality()
    {
        var seq = new QualitySequence("ACGT", new byte[] { 10, 20, 30, 40 });
        var rc = seq.GetReverseComplement() as QualitySequence;

        Assert.That(rc, Is.Not.Null);
        Assert.That(rc!.Sequence, Is.EqualTo("ACGT"));
        Assert.That(rc.Qualities, Is.EquivalentTo(new byte[] { 40, 30, 20, 10 }));
    }

    [Test]
    public void QualitySequence_Subsequence_PreservesQuality()
    {
        var seq = new QualitySequence("AACGTAA", new byte[] { 10, 20, 30, 40, 50, 60, 70 });
        var sub = seq.Subsequence(2, 3) as QualitySequence;

        Assert.That(sub, Is.Not.Null);
        Assert.That(sub!.Sequence, Is.EqualTo("CGT"));
        Assert.That(sub.Qualities, Is.EquivalentTo(new byte[] { 30, 40, 50 }));
    }

    #endregion

    #region SequenceBase Tests

    [Test]
    public void SequenceBase_Indexer_ReturnsCorrectChar()
    {
        var seq = new IupacDnaSequence("ACGT");

        Assert.That(seq[0], Is.EqualTo('A'));
        Assert.That(seq[2], Is.EqualTo('G'));
    }

    [Test]
    public void SequenceBase_ToString_ReturnsSequence()
    {
        var seq = new IupacDnaSequence("ACGT");

        Assert.That(seq.ToString(), Is.EqualTo("ACGT"));
    }

    [Test]
    public void SequenceBase_Equals_ComparesSequences()
    {
        var seq1 = new IupacDnaSequence("ACGT");
        var seq2 = new IupacDnaSequence("ACGT");
        var seq3 = new IupacDnaSequence("TGCA");

        Assert.That(seq1.Equals(seq2), Is.True);
        Assert.That(seq1.Equals(seq3), Is.False);
    }

    [Test]
    public void SequenceBase_GetHashCode_SameForEqualSequences()
    {
        var seq1 = new IupacDnaSequence("ACGT");
        var seq2 = new IupacDnaSequence("ACGT");

        Assert.That(seq1.GetHashCode(), Is.EqualTo(seq2.GetHashCode()));
    }

    #endregion
}
