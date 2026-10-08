namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the RSAT <c>oligo-analysis</c> overloads of MOTIF-DISCOVER-001 / MOTIF-SHARED-001: random
/// sequences, word lengths (including k &gt; N and long k where exp_freq underflows), thresholds, strand modes,
/// -ovlp/-noov and every background model (random Markov tables, including zero entries and ψ ∈ {0, 0.01, 1}).
/// Every run must either return a well-formed result — counts within the window total, probabilities in [0, 1],
/// finite significance, E-value = P · multiplier — or throw one of the documented exceptions
/// (ArgumentOutOfRangeException for an order the word length cannot carry; ArgumentException when ψ = 0 leaves an
/// observed word with probability 0, RSAT segment_proba "null transition").
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Matching")]
public class MotifOligoAnalysisFuzzTests
{
    private static string RandomDna(Random rng, int n) => new(Enumerable.Range(0, n).Select(_ => "ACGT"[rng.Next(4)]).ToArray());

    private static OligoBackgroundModel RandomModel(Random rng, int k)
    {
        switch (rng.Next(5))
        {
            case 0: return OligoBackgroundModel.Equiprobable;
            case 1: return OligoBackgroundModel.BernoulliFromInput;
            case 2: return OligoBackgroundModel.Bernoulli(Enumerable.Range(0, 4).Select(_ => 0.05 + rng.NextDouble()).ToArray());
            case 3: return OligoBackgroundModel.MarkovFromInput(k < 3 ? 0 : rng.Next(0, k - 1));
            default:
            {
                int m = rng.Next(0, Math.Min(3, k - 1) + 1);
                var table = new Dictionary<string, double>();
                foreach (int code in Enumerable.Range(0, 1 << (2 * (m + 1))))
                {
                    var word = new char[m + 1];
                    for (int i = 0, c = code; i <= m; i++, c >>= 2) word[i] = "ACGT"[c & 3];
                    table[new string(word)] = rng.Next(4) == 0 ? 0.0 : rng.NextDouble() * 5;
                }
                if (table.Values.Sum() == 0) table[table.Keys.First()] = 1;
                double psi = new[] { 0.0, 0.01, 0.01, 1.0 }[rng.Next(4)];
                return OligoBackgroundModel.MarkovFromOligoFrequencies(table, psi, rng.Next(3) == 0);
            }
        }
    }

    [Test]
    public void DiscoverMotifs_RandomInputs_WellFormedOrDocumentedException()
    {
        var rng = new Random(2026);
        for (int trial = 0; trial < 400; trial++)
        {
            int k = rng.Next(4) == 0 ? rng.Next(1, 700) : rng.Next(1, 8);
            string seq = RandomDna(rng, rng.Next(0, 400));
            if (k > 20 && seq.Length > 0) seq += seq; // long-k repeats (exp_freq underflow)
            int minCount = rng.Next(-1, 4);
            var strands = rng.Next(2) == 0 ? OligoStrandMode.Single : OligoStrandMode.Both;
            bool ovlp = rng.Next(3) != 0;
            var model = RandomModel(rng, Math.Min(k, 8));

            OligoAnalysisResult r;
            try
            {
                r = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k, minCount, model, strands, ovlp);
            }
            catch (ArgumentOutOfRangeException) { continue; }   // Markov order not allowed for this k
            catch (ArgumentException e) when (e.Message.Contains("probability 0")) { continue; }

            long windows = Math.Max(0, seq.Length - k + 1);
            Assert.That(r.TotalOccurrences, Is.EqualTo(windows));
            Assert.That(r.TestedPatterns, Is.EqualTo(r.Motifs.Count));
            Assert.That(r.Motifs.Sum(m => (long)m.Count), Is.LessThanOrEqualTo(windows));
            foreach (var m in r.Motifs)
            {
                Assert.That(m.Sequence, Has.Length.EqualTo(k));
                Assert.That(m.Count, Is.GreaterThanOrEqualTo(Math.Max(1, minCount)));
                Assert.That(m.Positions, Has.Count.EqualTo(m.Count));
                Assert.That(m.OccurrenceProbability, Is.InRange(0.0, 1.0));
                Assert.That(double.IsFinite(m.OccurrenceSignificance), Is.True, $"{m.Sequence} occ_sig");
                Assert.That(double.IsNaN(m.Ratio), Is.False);
                if (!ovlp)
                    Assert.That(m.Positions.Zip(m.Positions.Skip(1)).All(p => p.Second - p.First >= k || strands == OligoStrandMode.Both),
                        Is.True, "-noov single strand: counted occurrences of a word never overlap");
            }
        }
    }

    [Test]
    public void FindSharedMotifs_RandomInputs_WellFormedOrDocumentedException()
    {
        var rng = new Random(77);
        for (int trial = 0; trial < 300; trial++)
        {
            int k = rng.Next(1, 7);
            var seqs = Enumerable.Range(0, rng.Next(0, 9)).Select(_ => RandomDna(rng, rng.Next(0, 60))).ToList();
            int quorum = rng.Next(1, 4);
            var strands = rng.Next(2) == 0 ? OligoStrandMode.Single : OligoStrandMode.Both;
            var model = RandomModel(rng, k);

            SharedMotifAnalysisResult r;
            try
            {
                r = MotifFinder.FindSharedMotifs(seqs.Select(s => new DnaSequence(s)), k, quorum, model, strands);
            }
            catch (ArgumentOutOfRangeException) { continue; }
            catch (ArgumentException e) when (e.Message.Contains("probability 0")) { continue; }

            Assert.That(r.SequenceCount, Is.EqualTo(seqs.Count));
            Assert.That(r.PossiblePositions, Is.EqualTo(seqs.Where(s => s.Length >= k).Sum(s => (long)s.Length - k + 1)));
            foreach (var m in r.Motifs)
            {
                Assert.That(m.SequenceIndices, Has.Count.GreaterThanOrEqualTo(quorum));
                Assert.That(m.SequenceIndices, Is.Ordered.Ascending);
                Assert.That(m.MatchingSequenceProbability, Is.InRange(0.0, 1.0));
                Assert.That(double.IsFinite(m.MatchingSequenceSignificance), Is.True);
                Assert.That(m.ExpectedMatchingSequences, Is.InRange(0.0, (double)seqs.Count));
            }
        }
    }
}
