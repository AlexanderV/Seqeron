namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic relations for the RSAT <c>oligo-analysis</c> statistics (MOTIF-DISCOVER-001 / MOTIF-SHARED-001).
///   • RC:   with both strands (<c>-2str</c>) and a strand-symmetric background (equiprobable, input Bernoulli with
///           complementary residues pooled), reverse-complementing the input leaves every pair's occ, exp_freq and
///           occ_P unchanged (the pair W|W' is the same set of windows read on the other strand).
///   • PERM: permuting the input sequences of a matching-sequence run leaves every word's mseq, exp_ms and ms_P
///           unchanged (only the sequence indices are relabelled).
///   • SHORT: adding a sequence shorter than k adds one trial to the ms binomial (RSAT counts it in
///           sequence_number) without adding positions: S → S + 1, π = nb_pos / S shrinks accordingly.
/// Source: van Helden et al. 1998; RSAT oligo-analysis (CountOligos, CalcExpected, CalcProba).
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Matching")]
public class MotifOligoAnalysisMetamorphicTests
{
    private static string RandomDna(Random rng, int n) => new(Enumerable.Range(0, n).Select(_ => "ACGT"[rng.Next(4)]).ToArray());

    [Test]
    public void BothStrands_ReverseComplementInput_IsInvariant([Values(1, 2, 3, 4, 5)] int k)
    {
        var rng = new Random(100 + k);
        for (int trial = 0; trial < 20; trial++)
        {
            string seq = RandomDna(rng, rng.Next(k, 250));
            string rc = DnaSequence.GetReverseComplementString(seq);
            foreach (var bg in new[] { OligoBackgroundModel.Equiprobable, OligoBackgroundModel.BernoulliFromInput })
            {
                var a = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k, 1, bg, OligoStrandMode.Both).Motifs.ToDictionary(m => m.Sequence);
                var b = MotifFinder.DiscoverMotifs(new DnaSequence(rc), k, 1, bg, OligoStrandMode.Both).Motifs.ToDictionary(m => m.Sequence);
                Assert.That(b.Keys, Is.EquivalentTo(a.Keys));
                foreach (var (word, m) in a)
                {
                    Assert.That(b[word].Count, Is.EqualTo(m.Count), word);
                    Assert.That(b[word].ExpectedFrequency, Is.EqualTo(m.ExpectedFrequency).Within(1e-12).Percent, word);
                    Assert.That(b[word].OccurrenceProbability, Is.EqualTo(m.OccurrenceProbability).Within(1e-10).Percent, word);
                    // Positions are mirrored: window i on the input is window N - k - i on the reverse complement.
                    Assert.That(b[word].Positions, Is.EqualTo(m.Positions.Select(p => seq.Length - k - p).Order()), word);
                }
            }
        }
    }

    [Test]
    public void Shared_PermutingSequences_LeavesStatisticsUnchanged([Values] OligoStrandMode strands)
    {
        var rng = new Random(7);
        for (int trial = 0; trial < 20; trial++)
        {
            var seqs = Enumerable.Range(0, rng.Next(2, 7)).Select(_ => RandomDna(rng, rng.Next(3, 40))).ToList();
            var perm = seqs.OrderBy(_ => rng.Next()).ToList();
            var a = MotifFinder.FindSharedMotifs(seqs.Select(s => new DnaSequence(s)), 3, 2, OligoBackgroundModel.BernoulliFromInput, strands)
                .Motifs.ToDictionary(m => m.Sequence);
            var b = MotifFinder.FindSharedMotifs(perm.Select(s => new DnaSequence(s)), 3, 2, OligoBackgroundModel.BernoulliFromInput, strands)
                .Motifs.ToDictionary(m => m.Sequence);
            Assert.That(b.Keys, Is.EquivalentTo(a.Keys));
            foreach (var (word, m) in a)
            {
                Assert.That(b[word].SequenceIndices.Select(i => perm[i]), Is.EquivalentTo(m.SequenceIndices.Select(i => seqs[i])), word);
                Assert.That(b[word].ExpectedMatchingSequences, Is.EqualTo(m.ExpectedMatchingSequences).Within(1e-12).Percent, word);
                Assert.That(b[word].MatchingSequenceProbability, Is.EqualTo(m.MatchingSequenceProbability).Within(1e-10).Percent, word);
            }
        }
    }

    [Test]
    public void Shared_AddingShortSequence_AddsOneTrialWithoutPositions()
    {
        // RSAT: sequence_number counts every sequence, nb_possible_pos only those with length >= k.
        var seqs = new[] { "ACGTACGTTAGC", "TTACGTAGCAAC", "GGTAGCACGTTT", "CATTTTACG" }.Select(s => new DnaSequence(s)).ToList();
        var baseRun = MotifFinder.FindSharedMotifs(seqs, 4, 2, OligoBackgroundModel.Equiprobable);
        var withShort = MotifFinder.FindSharedMotifs(seqs.Append(new DnaSequence("ACG")), 4, 2, OligoBackgroundModel.Equiprobable);

        Assert.Multiple(() =>
        {
            Assert.That(withShort.SequenceCount, Is.EqualTo(baseRun.SequenceCount + 1));
            Assert.That(withShort.PossiblePositions, Is.EqualTo(baseRun.PossiblePositions));
            Assert.That(withShort.Motifs.Select(m => m.Sequence), Is.EqualTo(baseRun.Motifs.Select(m => m.Sequence)));
            foreach (var (a, b) in baseRun.Motifs.Zip(withShort.Motifs))
            {
                double p = a.ExpectedFrequency;
                double expected = 5 * -StatisticsHelper.ExpM1(33.0 / 5 * StatisticsHelper.Log1P(-p));
                Assert.That(b.ExpectedMatchingSequences, Is.EqualTo(expected).Within(1e-12).Percent, a.Sequence);
                Assert.That(b.MatchingSequenceProbability,
                    Is.EqualTo(StatisticsHelper.BinomialUpperTail(a.SequenceIndices.Count, 5, expected / 5)).Within(1e-10).Percent);
            }
        });
    }
}
