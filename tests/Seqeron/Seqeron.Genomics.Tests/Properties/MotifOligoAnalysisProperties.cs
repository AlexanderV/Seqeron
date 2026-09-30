using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the RSAT <c>oligo-analysis</c> statistics of MOTIF-DISCOVER-001 / MOTIF-SHARED-001
/// (<see cref="MotifFinder.DiscoverMotifs(DnaSequence, int, int, OligoBackgroundModel, OligoStrandMode, bool)"/>,
/// <see cref="MotifFinder.FindSharedMotifs(IEnumerable{DnaSequence}, int, int, OligoBackgroundModel, OligoStrandMode)"/>).
/// Oracles are the definitions (van Helden et al. 1998; RSAT manual): occ_P = P(X ≥ occ), X ~ Bin(n, exp_freq);
/// occ_E = occ_P · tested; occ_sig = −log₁₀ occ_E; -2str occ(W|W') = occ(W) + occ(W'); Markov-0 from input =
/// single-strand residue composition.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Matching")]
public class MotifOligoAnalysisProperties
{
    private static Arbitrary<string> Dna(int minLen, int maxLen) =>
        (from n in Gen.Choose(minLen, maxLen)
         from chars in Gen.Elements('A', 'C', 'G', 'T').ArrayOf(n)
         select new string(chars)).ToArbitrary();

    /// <summary>R: 0 ≤ occ_P ≤ 1, occ_E = occ_P · tested, occ_sig = −log₁₀ occ_E, ratio = occ / exp_occ.</summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property Significance_IsConsistentWithDefinitions()
    {
        return Prop.ForAll(Dna(10, 200), seq =>
        {
            var r = MotifFinder.DiscoverMotifs(new DnaSequence(seq), 3, 1, OligoBackgroundModel.BernoulliFromInput);
            bool ok = r.Motifs.All(m =>
                m.OccurrenceProbability is >= 0 and <= 1 &&
                Math.Abs(m.OccurrenceEValue - m.OccurrenceProbability * r.TestedPatterns) <= 1e-12 * Math.Max(1, m.OccurrenceEValue) &&
                Math.Abs(m.OccurrenceSignificance + Math.Log10(m.OccurrenceEValue)) <= 1e-9 &&
                Math.Abs(m.Ratio - m.Count / m.ExpectedOccurrences) <= 1e-9 * m.Ratio);
            return ok.Label("occ_P / occ_E / occ_sig / ratio inconsistent with their definitions");
        });
    }

    /// <summary>M: with a common expected frequency (equiprobable), significance is non-decreasing in occ.</summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property Significance_MonotoneInOccurrences_EquiprobableBackground()
    {
        return Prop.ForAll(Dna(10, 300), seq =>
        {
            var ordered = MotifFinder.DiscoverMotifs(new DnaSequence(seq), 2, 1, OligoBackgroundModel.Equiprobable)
                .Motifs.OrderBy(m => m.Count).ToList();
            bool ok = ordered.Zip(ordered.Skip(1)).All(p =>
                p.Second.OccurrenceSignificance >= p.First.OccurrenceSignificance - 1e-12 &&
                p.Second.OccurrenceProbability <= p.First.OccurrenceProbability * (1 + 1e-12));
            return ok.Label("a higher count gave a less significant word");
        });
    }

    /// <summary>E: Markov order 0 estimated from the input equals the input Bernoulli model on one strand.</summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property MarkovOrder0_EqualsInputBernoulli_SingleStrand()
    {
        return Prop.ForAll(Dna(6, 200), seq =>
        {
            var dna = new DnaSequence(seq);
            var markov = MotifFinder.DiscoverMotifs(dna, 3, 1, OligoBackgroundModel.MarkovFromInput(0)).Motifs;
            var bernoulli = MotifFinder.DiscoverMotifs(dna, 3, 1, OligoBackgroundModel.BernoulliFromInput).Motifs;
            bool ok = markov.Zip(bernoulli).All(p =>
                p.First.Sequence == p.Second.Sequence &&
                Math.Abs(p.First.ExpectedFrequency - p.Second.ExpectedFrequency) <= 1e-12 * p.Second.ExpectedFrequency &&
                Math.Abs(p.First.OccurrenceSignificance - p.Second.OccurrenceSignificance) <= 1e-9);
            return ok.Label("Markov-0 differs from the input Bernoulli model");
        });
    }

    /// <summary>E: on a sequence with uniform composition, Markov-0 = input Bernoulli = equiprobable (4^−k).</summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property UniformComposition_Markov0_EqualsEquiprobable()
    {
        var gen = (from n in Gen.Choose(2, 30)
                   from keys in Gen.Choose(0, 1000).ArrayOf(4 * n)
                   select new string(Enumerable.Range(0, 4 * n).Select(i => "ACGT"[i % 4]).Zip(keys)
                       .OrderBy(t => t.Second).Select(t => t.First).ToArray())).ToArbitrary();
        return Prop.ForAll(gen, seq =>
        {
            var markov = MotifFinder.DiscoverMotifs(new DnaSequence(seq), 3, 1, OligoBackgroundModel.MarkovFromInput(0)).Motifs;
            return markov.All(m => Math.Abs(m.ExpectedFrequency - 1.0 / 64) <= 1e-15)
                .Label("uniform composition must give exp_freq = 4^-3");
        });
    }

    /// <summary>
    /// R: -2str count = forward count + reverse-complement count (palindromes once); the pair counts partition the
    /// N − k + 1 windows; positions = merged positions of both members.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 80)]
    public Property BothStrands_CountIsForwardPlusReverseComplement()
    {
        return Prop.ForAll(Dna(4, 150), Gen.Choose(1, 4).ToArbitrary(), (seq, k) =>
        {
            var dna = new DnaSequence(seq);
            var single = MotifFinder.DiscoverMotifs(dna, k, 1, OligoBackgroundModel.Equiprobable).Motifs.ToDictionary(m => m.Sequence);
            var both = MotifFinder.DiscoverMotifs(dna, k, 1, OligoBackgroundModel.Equiprobable, OligoStrandMode.Both);
            bool ok = both.Motifs.All(m =>
            {
                string rc = DnaSequence.GetReverseComplementString(m.Sequence);
                var fwd = single.TryGetValue(m.Sequence, out var f) ? f.Positions : Array.Empty<int>();
                var rev = rc != m.Sequence && single.TryGetValue(rc, out var r) ? r.Positions : Array.Empty<int>();
                return m.ReverseComplement == rc && m.Count == fwd.Count + rev.Count &&
                       m.Positions.SequenceEqual(fwd.Concat(rev).Order());
            });
            long windows = Math.Max(0, seq.Length - k + 1);
            return (ok && both.Motifs.Sum(m => m.Count) == windows).Label("-2str pair count != forward + reverse complement");
        });
    }

    /// <summary>R: ms_P ∈ [0, 1]; ms_E = ms_P · NPO; the quorum words equal the quorum-only overload's.</summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property Shared_Significance_IsConsistent()
    {
        var gen = Dna(3, 30).Generator.ArrayOf().Where(a => a.Length is >= 1 and <= 6).ToArbitrary();
        return Prop.ForAll(gen, seqs =>
        {
            var dnas = seqs.Select(s => new DnaSequence(s)).ToList();
            var r = MotifFinder.FindSharedMotifs(dnas, 3, 2, OligoBackgroundModel.BernoulliFromInput);
            var quorum = MotifFinder.FindSharedMotifs(dnas, 3, 2).Select(m => m.Sequence).ToList();
            bool ok = r.Motifs.Select(m => m.Sequence).SequenceEqual(quorum) && r.Motifs.All(m =>
                m.MatchingSequenceProbability is >= 0 and <= 1 &&
                Math.Abs(m.MatchingSequenceEValue - m.MatchingSequenceProbability * r.PossibleOligos) <= 1e-12 * Math.Max(1, m.MatchingSequenceEValue) &&
                m.ExpectedMatchingSequences <= r.SequenceCount);
            return ok.Label("ms statistics inconsistent");
        });
    }
}
