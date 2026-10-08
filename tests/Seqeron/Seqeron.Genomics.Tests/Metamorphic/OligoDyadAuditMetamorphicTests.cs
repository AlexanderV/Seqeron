namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic / property relations for the review 2026-09 B05 audit additions on MOTIF-DISCOVER-001 /
/// MOTIF-SHARED-001 (docs/Validation/review-2026-09/B05.md F31 RSAT oligo-analysis options + dyad-analysis, F33
/// <c>-seqtype</c>, F35 NPO overflow).
///   • OLIGO-DEF     AnalyzeOligos with default-equivalent options (no pseudo, no degeneracy, no calibration) ==
///                   DiscoverMotifs field for field, over random inputs / backgrounds / strands / overlap modes.
///   • OLIGO-PSI0    PseudoFrequency = 0 is the default (identical report).
///   • SEQTYPE-DNA   AnalyzeOligoStrings(SequenceType = Dna) on pure A/C/G/T (any case) == AnalyzeOligos(DnaSequence).
///   • SEQTYPE-OTHER on pure A/C/G/T with an equiprobable background, -seqtype other learns the alphabet from the
///                   input: with all four letters present the counts and occ_P equal the DNA run (|A| = 4).
///   • NPO-INF       k ≥ 512 (4^k overflows): every statistic of a -pseudo / both-strand run is non-NaN (F35).
///   • DYAD-BRUTE    single strand, overlapping counting: every dyad's Occurrences / Positions == brute-force count
///                   of M₁ at i and M₂ at i + m + s.
///   • DYAD-RANGE    widening the spacing range keeps every dyad of the narrower range with the same occurrences.
///   • DYAD-TYPE     DirectRepeat / InvertedRepeat dyads are the Any dyads with M₂ = M₁ / rc(M₁), same occurrences.
///   • DYAD-MINCOUNT MinCount = c + 1 keeps exactly the MinCount = c dyads with occ ≥ c + 1.
/// Fixed seeds, bounded sizes.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Matching")]
public class OligoDyadAuditMetamorphicTests
{
    private static string RandomDna(Random rng, int n, string alphabet = "ACGT") =>
        new(Enumerable.Range(0, n).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static OligoBackgroundModel RandomBackground(Random rng, int k) => rng.Next(4) switch
    {
        0 => OligoBackgroundModel.Equiprobable,
        1 => OligoBackgroundModel.BernoulliFromInput,
        2 when k > 2 => OligoBackgroundModel.MarkovFromInput(rng.Next(1, k - 1)),
        _ => OligoBackgroundModel.Bernoulli(new[] { 0.3, 0.2, 0.2, 0.3 }),
    };

    private static void AssertSameReport(OligoAnalysisReport a, OligoAnalysisReport b, string msg)
    {
        Assert.That(b.TestedPatterns, Is.EqualTo(a.TestedPatterns), msg);
        Assert.That(b.TotalOccurrences, Is.EqualTo(a.TotalOccurrences), msg);
        Assert.That(b.PossiblePositions, Is.EqualTo(a.PossiblePositions), msg);
        Assert.That(b.PossibleOligos, Is.EqualTo(a.PossibleOligos), msg);
        Assert.That(b.Patterns.Count, Is.EqualTo(a.Patterns.Count), msg);
        for (int i = 0; i < a.Patterns.Count; i++)
        {
            var x = a.Patterns[i];
            var y = b.Patterns[i];
            Assert.That(y.Pattern, Is.EqualTo(x.Pattern), msg);
            Assert.That(y.Occurrences, Is.EqualTo(x.Occurrences), msg);
            Assert.That(y.Positions, Is.EqualTo(x.Positions), msg);
            Assert.That(y.ExpectedFrequency, Is.EqualTo(x.ExpectedFrequency), msg);
            Assert.That(y.ExpectedVariance, Is.EqualTo(x.ExpectedVariance), msg);
            Assert.That(y.ZScore, Is.EqualTo(x.ZScore), msg);
            Assert.That(y.OccurrenceProbability, Is.EqualTo(x.OccurrenceProbability), msg);
            Assert.That(y.OccurrenceSignificance, Is.EqualTo(x.OccurrenceSignificance), msg);
        }
    }

    [Test]
    public void AnalyzeOligos_DefaultOptions_EqualDiscoverMotifs_AndPseudoZeroIsDefault()
    {
        var rng = new Random(310_001);
        for (int trial = 0; trial < 150; trial++)
        {
            int k = rng.Next(1, 6);
            string seq = RandomDna(rng, rng.Next(1, 200), trial % 4 == 0 ? "AACGTT" : "ACGT");
            var bg = RandomBackground(rng, k);
            var strands = rng.Next(2) == 0 ? OligoStrandMode.Single : OligoStrandMode.Both;
            bool overlapping = rng.Next(2) == 0;
            int minCount = rng.Next(0, 4);
            string msg = $"trial {trial} k={k} {strands} ov={overlapping} min={minCount}";

            var a = MotifFinder.DiscoverMotifs(new DnaSequence(seq), k, minCount, bg, strands, overlapping);
            var options = new OligoAnalysisOptions { Background = bg, Strands = strands, CountOverlapping = overlapping, MinCount = minCount };
            var b = MotifFinder.AnalyzeOligos(new[] { new DnaSequence(seq) }, k, options);
            Assert.That(b.TestedPatterns, Is.EqualTo(a.TestedPatterns), msg);
            Assert.That(b.TotalOccurrences, Is.EqualTo(a.TotalOccurrences), msg);
            Assert.That(b.PossibleOligos, Is.EqualTo(a.PossibleOligos), msg);
            Assert.That(b.Patterns.Count, Is.EqualTo(a.Motifs.Count), msg);
            for (int i = 0; i < a.Motifs.Count; i++)
            {
                var x = a.Motifs[i];
                var y = b.Patterns[i];
                Assert.That(y.Pattern, Is.EqualTo(x.Sequence), msg);
                Assert.That(y.ReverseComplement, Is.EqualTo(x.ReverseComplement), msg);
                Assert.That(y.Occurrences, Is.EqualTo(x.Count), msg);
                Assert.That(y.Positions.Select(p => p.Position), Is.EqualTo(x.Positions), msg);
                Assert.That(y.ExpectedFrequency, Is.EqualTo(x.ExpectedFrequency), msg);
                Assert.That(y.ExpectedOccurrences, Is.EqualTo(x.ExpectedOccurrences), msg);
                Assert.That(y.Ratio, Is.EqualTo(x.Ratio), msg);
                Assert.That(y.OccurrenceProbability, Is.EqualTo(x.OccurrenceProbability), msg);
                Assert.That(y.OccurrenceEValue, Is.EqualTo(x.OccurrenceEValue), msg);
                Assert.That(y.OccurrenceSignificance, Is.EqualTo(x.OccurrenceSignificance), msg);
            }

            AssertSameReport(b, MotifFinder.AnalyzeOligos(new[] { new DnaSequence(seq) }, k, options with { PseudoFrequency = 0 }), msg + " psi0");
        }
    }

    [Test]
    public void AnalyzeOligoStrings_DnaOnPureAcgt_EqualsDnaSequencePath()
    {
        var rng = new Random(330_001);
        for (int trial = 0; trial < 120; trial++)
        {
            int k = rng.Next(1, 6);
            var seqs = Enumerable.Range(0, rng.Next(1, 5)).Select(_ => RandomDna(rng, rng.Next(1, 120))).ToList();
            var options = new OligoAnalysisOptions
            {
                Background = RandomBackground(rng, k),
                Strands = rng.Next(2) == 0 ? OligoStrandMode.Single : OligoStrandMode.Both,
                CountOverlapping = rng.Next(2) == 0,
                MinCount = rng.Next(0, 3),
                PseudoFrequency = rng.Next(2) == 0 ? 0 : 0.01,
                Degeneracy = (OligoDegeneracy)rng.Next(3),
            };
            string msg = $"trial {trial} k={k} {options}";
            var dna = MotifFinder.AnalyzeOligos(seqs.Select(s => new DnaSequence(s)), k, options);
            var strings = MotifFinder.AnalyzeOligoStrings(seqs, k, options);
            AssertSameReport(dna, strings, msg);
            Assert.That(strings.AlphabetSize, Is.EqualTo(4));
            AssertSameReport(dna, MotifFinder.AnalyzeOligoStrings(seqs.Select(s => s.ToLowerInvariant()), k, options), msg + " lower");
        }
    }

    [Test]
    public void AnalyzeOligoStrings_OtherOnAcgt_CountsEqualDna()
    {
        var rng = new Random(330_002);
        for (int trial = 0; trial < 80; trial++)
        {
            int k = rng.Next(1, 4);
            var seqs = new List<string> { "ACGT" + RandomDna(rng, rng.Next(0, 80)) };
            var dnaOptions = new OligoAnalysisOptions { Background = OligoBackgroundModel.Equiprobable };
            var otherOptions = dnaOptions with { SequenceType = OligoSequenceType.Other };
            var dna = MotifFinder.AnalyzeOligoStrings(seqs, k, dnaOptions);
            var other = MotifFinder.AnalyzeOligoStrings(seqs, k, otherOptions);
            string msg = $"trial {trial} k={k} {seqs[0]}";
            Assert.That(other.AlphabetSize, Is.EqualTo(4), msg);
            Assert.That(other.PossibleOligos, Is.EqualTo(dna.PossibleOligos), msg);
            Assert.That(other.TotalOccurrences, Is.EqualTo(dna.TotalOccurrences), msg);
            var d = dna.Patterns.ToDictionary(p => p.Pattern.ToUpperInvariant());
            Assert.That(other.Patterns.Select(p => p.Pattern.ToUpperInvariant()), Is.EquivalentTo(d.Keys), msg);
            foreach (var p in other.Patterns)
            {
                var q = d[p.Pattern.ToUpperInvariant()];
                Assert.That(p.Occurrences, Is.EqualTo(q.Occurrences), msg);
                Assert.That(p.ExpectedFrequency, Is.EqualTo(q.ExpectedFrequency).Within(1e-12).Percent, msg);
            }
        }
    }

    [Test]
    public void AnalyzeOligos_PossibleOligosOverflow_NoNaNStatistics()
    {
        var rng = new Random(350_001);
        foreach (int k in new[] { 512, 520 })
        {
            var seqs = Enumerable.Range(0, 3).Select(_ => new DnaSequence(RandomDna(rng, 560))).ToList();
            seqs.Add(seqs[0]); // a repeated word so that some occurrence counts exceed 1
            foreach (var strands in new[] { OligoStrandMode.Single, OligoStrandMode.Both })
            {
                var r = MotifFinder.AnalyzeOligos(seqs, k, new OligoAnalysisOptions
                {
                    Background = OligoBackgroundModel.Equiprobable, Strands = strands, PseudoFrequency = 0.01,
                });
                Assert.That(r.PossibleOligos, Is.EqualTo(double.PositiveInfinity), $"k={k} {strands}");
                Assert.That(r.Patterns, Is.Not.Empty);
                foreach (var p in r.Patterns)
                {
                    Assert.That(double.IsNaN(p.ExpectedFrequency), Is.False, $"k={k} {strands} exp_freq");
                    Assert.That(double.IsNaN(p.ExpectedOccurrences), Is.False, $"k={k} {strands} exp_occ");
                    Assert.That(double.IsNaN(p.OccurrenceProbability), Is.False, $"k={k} {strands} occ_P");
                    Assert.That(double.IsNaN(p.OccurrenceSignificance), Is.False, $"k={k} {strands} occ_sig");
                }
            }
        }
    }

    private static DyadAnalysisReport Dyads(IReadOnlyList<string> seqs, DyadAnalysisOptions o) =>
        MotifFinder.AnalyzeDyads(seqs.Select(s => new DnaSequence(s)), o);

    [Test]
    public void Dyads_SingleStrandOverlapping_EqualBruteForce()
    {
        var rng = new Random(310_101);
        for (int trial = 0; trial < 60; trial++)
        {
            int m = rng.Next(1, 4);
            int minS = rng.Next(0, 4), maxS = minS + rng.Next(0, 5);
            var seqs = Enumerable.Range(0, rng.Next(1, 4)).Select(_ => RandomDna(rng, rng.Next(1, 60), "ACGT")).ToList();
            var r = Dyads(seqs, new DyadAnalysisOptions
            {
                MonadLength = m, MinSpacing = minS, MaxSpacing = maxS, Strands = OligoStrandMode.Single,
                CountOverlapping = true, MinCount = 1,
            });
            var expected = new Dictionary<string, List<OligoOccurrence>>();
            for (int si = 0; si < seqs.Count; si++)
            {
                string s = seqs[si];
                for (int sp = minS; sp <= maxS; sp++)
                    for (int i = 0; i + 2 * m + sp <= s.Length; i++)
                    {
                        string key = $"{s.Substring(i, m)}n{{{sp}}}{s.Substring(i + m + sp, m)}";
                        if (!expected.TryGetValue(key, out var list)) expected[key] = list = new List<OligoOccurrence>();
                        list.Add(new OligoOccurrence(si, i));
                    }
            }
            string msg = $"trial {trial} m={m} sp={minS}-{maxS}";
            Assert.That(r.Dyads.Select(d => d.Pattern), Is.EquivalentTo(expected.Keys), msg);
            foreach (var d in r.Dyads)
            {
                Assert.That(d.Occurrences, Is.EqualTo(expected[d.Pattern].Count), msg + " " + d.Pattern);
                Assert.That(d.Positions, Is.EqualTo(expected[d.Pattern]), msg + " " + d.Pattern);
                Assert.That(d.Pattern, Is.EqualTo($"{d.FirstMonad}n{{{d.Spacing}}}{d.SecondMonad}"), msg);
            }
        }
    }

    [Test]
    public void Dyads_RangeTypeAndMinCount_Monotone()
    {
        var rng = new Random(310_102);
        for (int trial = 0; trial < 60; trial++)
        {
            int m = rng.Next(1, 4);
            var seqs = Enumerable.Range(0, rng.Next(1, 4)).Select(_ => RandomDna(rng, rng.Next(5, 80))).ToList();
            var strands = rng.Next(2) == 0 ? OligoStrandMode.Single : OligoStrandMode.Both;
            bool overlapping = rng.Next(2) == 0;
            int minS = rng.Next(0, 4), maxS = minS + rng.Next(0, 4);
            var narrow = new DyadAnalysisOptions
            {
                MonadLength = m, MinSpacing = minS, MaxSpacing = maxS, Strands = strands, CountOverlapping = overlapping, MinCount = 1,
            };
            var wide = narrow with { MinSpacing = Math.Max(0, minS - rng.Next(0, 3)), MaxSpacing = maxS + rng.Next(0, 4) };
            string msg = $"trial {trial} m={m} {strands} ov={overlapping}";

            var rn = Dyads(seqs, narrow);
            var rw = Dyads(seqs, wide).Dyads.ToDictionary(d => d.Pattern);
            foreach (var d in rn.Dyads)
            {
                Assert.That(rw.ContainsKey(d.Pattern), Is.True, msg + " range " + d.Pattern);
                Assert.That(rw[d.Pattern].Occurrences, Is.EqualTo(d.Occurrences), msg + " range " + d.Pattern);
                Assert.That(rw[d.Pattern].Positions, Is.EqualTo(d.Positions), msg + " range " + d.Pattern);
            }

            var any = rn.Dyads.ToDictionary(d => d.Pattern);
            foreach (var (type, keep) in new (DyadType, Func<DyadStatistics, bool>)[]
                     {
                         (DyadType.DirectRepeat, d => d.SecondMonad == d.FirstMonad),
                         (DyadType.InvertedRepeat, d => d.SecondMonad == DnaSequence.GetReverseComplementString(d.FirstMonad)),
                     })
            {
                var typed = Dyads(seqs, narrow with { Type = type }).Dyads;
                foreach (var d in typed)
                {
                    Assert.That(keep(d), Is.True, $"{msg} {type} {d.Pattern}");
                    Assert.That(any.ContainsKey(d.Pattern), Is.True, $"{msg} {type} {d.Pattern}");
                    Assert.That(d.Occurrences, Is.EqualTo(any[d.Pattern].Occurrences), $"{msg} {type} {d.Pattern}");
                }
                if (strands == OligoStrandMode.Single)
                    Assert.That(typed.Select(d => d.Pattern), Is.EquivalentTo(rn.Dyads.Where(keep).Select(d => d.Pattern)), $"{msg} {type}");
            }

            int c = rng.Next(1, 4);
            var low = Dyads(seqs, narrow with { MinCount = c }).Dyads;
            var high = Dyads(seqs, narrow with { MinCount = c + 1 }).Dyads;
            Assert.That(high.Select(d => d.Pattern), Is.EquivalentTo(low.Where(d => d.Occurrences >= c + 1).Select(d => d.Pattern)), msg + " mincount");
        }
    }
}
