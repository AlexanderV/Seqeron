namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A6-8d, B09): the scenarios of the former <c>SequenceAligner_MultipleAlign_Benchmark</c> fixture ([Explicit] - no tier ever ran it, and it ended in <c>Assert.Pass</c>) as deterministic tests: the exact center-star
/// alignment (<see cref="SequenceAligner.MultipleAlign(IEnumerable{DnaSequence}, ScoringMatrix?)"/>) and the anchor-based one (<see cref="SequenceAligner.MultipleAlignAnchored(IEnumerable{DnaSequence}, ScoringMatrix?)"/>) on related sequences
/// (95 / 90 / 85 / 80 / 60 % similarity, point mutations and ~1 % indels) and on unrelated random sequences (the worst case: no anchors, the NW fallback) must give a valid MSA: every input aligned, equal row lengths, the gaps removed recover the input
/// byte for byte. The timing table moved to apps/Seqeron.Alignment.Benchmarks (wall-clock numbers are not a test property).
/// </summary>
[TestFixture]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlignAnchored_Scenarios_Tests
{
    private static string Dna(int length, int seed)
    {
        var rng = new Random(seed);
        return new string(Enumerable.Range(0, length).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
    }

    private static List<DnaSequence> Related(int count, int length, double similarity, int seed)
    {
        var rng = new Random(seed);
        string reference = Dna(length, seed);
        var sequences = new List<DnaSequence> { new(reference) };
        for (int s = 1; s < count; s++)
        {
            var chars = reference.ToCharArray();
            for (int m = 0; m < (int)(length * (1.0 - similarity)); m++)
            {
                int pos = rng.Next(length);
                char replacement;
                do
                {
                    replacement = "ACGT"[rng.Next(4)];
                }
                while (replacement == chars[pos]);
                chars[pos] = replacement;
            }

            var sb = new System.Text.StringBuilder(new string(chars));
            for (int d = 0; d < Math.Max(1, length / 100) && sb.Length > 10; d++)
            {
                int pos = rng.Next(sb.Length);
                if (rng.NextDouble() < 0.5)
                    sb.Insert(pos, "ACGT"[rng.Next(4)]);
                else
                    sb.Remove(pos, 1);
            }

            sequences.Add(new DnaSequence(sb.ToString()));
        }

        return sequences;
    }

    private static List<DnaSequence> RandomSet(int count, int length, int seed)
        => [.. Enumerable.Range(0, count).Select(i => new DnaSequence(Dna(length, seed + (i * 1000))))];

    private static void AssertValidMsa(MultipleAlignmentResult result, IReadOnlyList<DnaSequence> input, string who)
    {
        Assert.That(result.AlignedSequences, Has.Length.EqualTo(input.Count), who);
        Assert.That(result.AlignedSequences.Select(r => r.Length).Distinct().Count(), Is.EqualTo(1), $"{who}: equal row lengths");
        for (int i = 0; i < input.Count; i++)
            Assert.That(result.AlignedSequences[i].Replace("-", ""), Is.EqualTo(input[i].Sequence), $"{who}: row {i} without gaps is the input");
    }

    private static void Check(List<DnaSequence> input, bool classic = true)
    {
        AssertValidMsa(SequenceAligner.MultipleAlignAnchored(input, SequenceAligner.SimpleDna), input, "anchored");
        if (classic)
            AssertValidMsa(SequenceAligner.MultipleAlign(input, SequenceAligner.SimpleDna), input, "center star");
    }

    [Test] public void Small_5x200_90pct() => Check(Related(5, 200, 0.90, 42));
    [Test] public void Medium_10x500_85pct() => Check(Related(10, 500, 0.85, 123));
    [Test] public void Large_10x1000_90pct() => Check(Related(10, 1000, 0.90, 456));
    [Test] public void HighSimilarity_10x500_95pct() => Check(Related(10, 500, 0.95, 111));
    [Test] public void MediumSimilarity_10x500_80pct() => Check(Related(10, 500, 0.80, 222));
    [Test] public void LowSimilarity_10x500_60pct() => Check(Related(10, 500, 0.60, 333));
    [Test] public void Random_10x500_NoSimilarity_FallsBackToNeedlemanWunsch() => Check(RandomSet(10, 500, 999));

    [Test]
    public void ExtraLarge_20x2000_90pct_Anchored()
    {
        // the stress size of the old benchmark: only the anchor-based aligner (the center star computes 19 exact 2000 x 2000 alignments)
        Check(Related(20, 2000, 0.90, 789), classic: false);
    }
}
