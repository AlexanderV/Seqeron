using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceAligner_GlobalAlign_Tests
{
    private static readonly ScoringMatrix WikipediaSimpleScoring = new(
        Match: 1,
        Mismatch: -1,
        GapOpen: -1,
        GapExtend: -1);

    [Test]
    public void GlobalAlign_WikipediaExample_ReconstructsInputsAndScoreMatches()
    {
        var seq1 = new DnaSequence("GCATGCG");
        var seq2 = new DnaSequence("GATTACA");

        var result = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaSimpleScoring);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignmentType, Is.EqualTo(AlignmentType.Global));
            Assert.That(result.AlignedSequence1.Length, Is.EqualTo(result.AlignedSequence2.Length));
            Assert.That(RemoveGaps(result.AlignedSequence1), Is.EqualTo(seq1.Sequence));
            Assert.That(RemoveGaps(result.AlignedSequence2), Is.EqualTo(seq2.Sequence));
            Assert.That(result.Score, Is.EqualTo(CalculateScore(result.AlignedSequence1, result.AlignedSequence2, WikipediaSimpleScoring)));
        });
    }

    [Test]
    public void GlobalAlign_StringOverload_MatchesDnaSequenceResult()
    {
        const string seq1 = "GCATGCG";
        const string seq2 = "GATTACA";

        var dnaResult = SequenceAligner.GlobalAlign(new DnaSequence(seq1), new DnaSequence(seq2), WikipediaSimpleScoring);
        var stringResult = SequenceAligner.GlobalAlign(seq1, seq2, WikipediaSimpleScoring);

        Assert.Multiple(() =>
        {
            Assert.That(stringResult.AlignmentType, Is.EqualTo(AlignmentType.Global));
            Assert.That(stringResult.Score, Is.EqualTo(dnaResult.Score));
            Assert.That(stringResult.AlignedSequence1, Is.EqualTo(dnaResult.AlignedSequence1));
            Assert.That(stringResult.AlignedSequence2, Is.EqualTo(dnaResult.AlignedSequence2));
        });
    }

    [Test]
    public void GlobalAlign_StringOverload_EmptyInput_ReturnsEmpty()
    {
        var result = SequenceAligner.GlobalAlign(string.Empty, "GATTACA", WikipediaSimpleScoring);

        Assert.That(result, Is.EqualTo(AlignmentResult.Empty));
    }

    [Test]
    public void GlobalAlign_NullSequence_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.GlobalAlign((DnaSequence)null!, new DnaSequence("GATTACA"), WikipediaSimpleScoring));
    }

    private static string RemoveGaps(string alignedSequence)
    {
        return new string(alignedSequence.Where(c => c != '-').ToArray());
    }

    private static int CalculateScore(string aligned1, string aligned2, ScoringMatrix scoring)
    {
        int score = 0;

        for (int i = 0; i < aligned1.Length; i++)
        {
            char a = aligned1[i];
            char b = aligned2[i];

            if (a == '-' || b == '-')
            {
                score += scoring.GapExtend;
            }
            else if (a == b)
            {
                score += scoring.Match;
            }
            else
            {
                score += scoring.Mismatch;
            }
        }

        return score;
    }
}
