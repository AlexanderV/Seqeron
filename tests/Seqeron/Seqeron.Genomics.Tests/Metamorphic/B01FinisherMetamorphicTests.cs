namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the behaviours added by the B01 finisher
/// (docs/Validation/review-2026-09/B01.md, F14–F31). Fixed-seed random inputs; every relation is exact
/// (all skew sums are sums of identical doubles in identical order, or of integers).
///
/// Test Units: SEQ-GCSKEW-001, SEQ-ATSKEW-001, SEQ-REPLICATION-001, SEQ-GC-001, SEQ-GC-ANALYSIS-001, B01-SWEEP.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
[Category("Composition")]
public class B01FinisherMetamorphicTests
{
    private const int Seed = 20261009;

    private static string Random(Random rng, string alphabet, int length) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static string Map(string s, Func<char, char> f) => new(s.Select(f).ToArray());

    // A→G, T→C, G→A, C→T (both cases): GC counts of the image are the AT counts of the source.
    private static char AtToGc(char c) => c switch
    {
        'A' => 'G', 'T' => 'C', 'G' => 'A', 'C' => 'T',
        'a' => 'g', 't' => 'c', 'g' => 'a', 'c' => 't',
        _ => c,
    };

    private static char SwapGc(char c) => c switch { 'G' => 'C', 'C' => 'G', 'g' => 'c', 'c' => 'g', _ => c };

    /// <summary>
    /// MR (F15): AT skew of s equals GC skew of the A→G, T→C (G→A, C→T) image of s, for windowed (any step,
    /// both partial flags) and cumulative forms; positions and window bounds are unchanged.
    /// </summary>
    [Test]
    public void AtSkew_EqualsGcSkewOfMappedSequence()
    {
        var rng = new Random(Seed);
        for (int trial = 0; trial < 400; trial++)
        {
            string s = Random(rng, "ACGTacgtNU-", rng.Next(0, 200));
            string m = Map(s, AtToGc);
            int w = rng.Next(1, 30), step = rng.Next(1, 30);
            bool partial = rng.Next(2) == 0;

            GcSkewCalculator.CalculateWindowedAtSkew(s, w, step, partial)
                .Select(p => (p.Position, p.AtSkew, p.WindowStart, p.WindowEnd))
                .Should().Equal(GcSkewCalculator.CalculateWindowedGcSkew(m, w, step, partial)
                    .Select(p => (p.Position, p.GcSkew, p.WindowStart, p.WindowEnd)), s);
            GcSkewCalculator.CalculateCumulativeAtSkew(s, w, partial)
                .Select(p => (p.Position, p.AtSkew, p.CumulativeAtSkew))
                .Should().Equal(GcSkewCalculator.CalculateCumulativeGcSkew(m, w, partial)
                    .Select(p => (p.Position, p.GcSkew, p.CumulativeGcSkew)), s);
            GcSkewCalculator.CalculateAtSkew(s).Should().Be(GcSkewCalculator.CalculateGcSkew(m), s);
        }
    }

    /// <summary>
    /// MR (F16/F17): swapping G↔C negates every prefix skew, so minimizer and maximizer sets exchange
    /// (linear and circular), the per-base and windowed origin/terminus exchange with negated skews,
    /// and the SkewIT Skew Index (|Σ−Σ| of flipped signs) is unchanged.
    /// </summary>
    [Test]
    public void GcSwap_ExchangesExtrema_AndPreservesSkewIndex()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = Random(rng, "ACGTacgtN", rng.Next(1, 300));
            string t = Map(s, SwapGc);

            foreach (bool circular in new[] { false, true })
            {
                GcSkewCalculator.FindMinimumSkewPositions(t, circular).Should()
                    .Equal(GcSkewCalculator.FindMaximumSkewPositions(s, circular), s);
                GcSkewCalculator.FindMaximumSkewPositions(t, circular).Should()
                    .Equal(GcSkewCalculator.FindMinimumSkewPositions(s, circular), s);
                var a = GcSkewCalculator.PredictReplicationOrigin(s, circular);
                var b = GcSkewCalculator.PredictReplicationOrigin(t, circular);
                (b.PredictedOrigin, b.PredictedTerminus, b.OriginSkew, b.TerminusSkew)
                    .Should().Be((a.PredictedTerminus, a.PredictedOrigin, -a.TerminusSkew, -a.OriginSkew), s);
            }

            int w = rng.Next(1, 40);
            var wa = GcSkewCalculator.PredictReplicationOrigin(s, w);
            var wb = GcSkewCalculator.PredictReplicationOrigin(t, w);
            (wb.PredictedOrigin, wb.PredictedTerminus, wb.OriginSkew, wb.TerminusSkew, wb.IsSignificant)
                .Should().Be((wa.PredictedTerminus, wa.PredictedOrigin, -wa.TerminusSkew, -wa.OriginSkew, wa.IsSignificant), $"{s} w={w}");

            int k = rng.Next(1, 12);
            GcSkewCalculator.CalculateSkewIndex(t, k).Should().Be(GcSkewCalculator.CalculateSkewIndex(s, k), $"{s} k={k}");
        }
    }

    /// <summary>
    /// MR (F17): on a circular genome with total G−C = 0, rotating by r maps the circular minimizer/maximizer
    /// sets by p ↦ (p − r) mod n and shifts the extreme values by −Skew_r; PredictReplicationOrigin(circular)
    /// reports the first element of the mapped set.
    /// </summary>
    [Test]
    public void CircularRotation_BalancedGenome_IsEquivariant()
    {
        var rng = new Random(Seed + 2);
        for (int trial = 0; trial < 200; trial++)
        {
            int pairs = rng.Next(0, 60);
            var chars = Enumerable.Repeat('G', pairs).Concat(Enumerable.Repeat('C', pairs))
                .Concat(Random(rng, "AT", rng.Next(1, 80))).OrderBy(_ => rng.Next()).ToArray();
            string s = new(chars);
            int n = s.Length, r = rng.Next(0, n);
            string rot = s[r..] + s[..r];
            int skewR = s[..r].Count(c => c == 'G') - s[..r].Count(c => c == 'C');

            List<int> MapSet(IReadOnlyList<int> ps) => ps.Select(p => ((p - r) % n + n) % n).OrderBy(p => p).ToList();
            var minRot = GcSkewCalculator.FindMinimumSkewPositions(rot, circular: true);
            var maxRot = GcSkewCalculator.FindMaximumSkewPositions(rot, circular: true);
            minRot.Should().Equal(MapSet(GcSkewCalculator.FindMinimumSkewPositions(s, circular: true)), $"{s} r={r}");
            maxRot.Should().Equal(MapSet(GcSkewCalculator.FindMaximumSkewPositions(s, circular: true)), $"{s} r={r}");

            var a = GcSkewCalculator.PredictReplicationOrigin(s, circular: true);
            var b = GcSkewCalculator.PredictReplicationOrigin(rot, circular: true);
            b.PredictedOrigin.Should().Be(minRot[0]);
            b.PredictedTerminus.Should().Be(maxRot[0]);
            b.OriginSkew.Should().Be(a.OriginSkew - skewR);
            b.TerminusSkew.Should().Be(a.TerminusSkew - skewR);
        }
    }

    /// <summary>
    /// MR (F14): appending a suffix never changes the windows that were complete before (prefix stability of the
    /// default output), and the partial output of s is the default output of s followed by the windows whose
    /// value equals the GC skew of the corresponding suffix.
    /// </summary>
    [Test]
    public void PartialWindows_PrefixStableUnderAppend()
    {
        var rng = new Random(Seed + 3);
        for (int trial = 0; trial < 300; trial++)
        {
            string s = Random(rng, "ACGTacgtN", rng.Next(0, 150));
            string tail = Random(rng, "ACGTN", rng.Next(0, 60));
            int w = rng.Next(1, 25), step = rng.Next(1, 25);

            var before = GcSkewCalculator.CalculateWindowedGcSkew(s, w, step).ToList();
            var after = GcSkewCalculator.CalculateWindowedGcSkew(s + tail, w, step).ToList();
            after.Take(before.Count).Should().Equal(before, $"{s}|{tail} w={w} step={step}");

            var partial = GcSkewCalculator.CalculateWindowedGcSkew(s, w, step, includePartialWindow: true).ToList();
            partial.Take(before.Count).Should().Equal(before);
            foreach (var p in partial.Skip(before.Count))
                p.GcSkew.Should().Be(GcSkewCalculator.CalculateGcSkew(s[p.WindowStart..]));
        }
    }

    /// <summary>
    /// MR (F18/F19/F20/F21): inserting non-ASCII letters ('ſ' U+017F, 'ı' U+0131, Kelvin U+212A) never adds
    /// G/C or IUPAC content — Remove-mode GC fraction is unchanged, Ignore/Weighted equal the original numerator
    /// over the longer length, IUPAC validation fails at the first inserted char, CountKmersSpan(1) gains exactly
    /// the inserted symbols as new keys, and Hamming distance vs the original with each inserted position
    /// replaced by 'S' rises by exactly the number of insertions.
    /// </summary>
    [Test]
    public void NonAsciiInsertion_NeverCountsAsNucleotide()
    {
        var rng = new Random(Seed + 4);
        const string nonAscii = "\u017F\u0131\u212A"; // ſ, ı, Kelvin sign
        for (int trial = 0; trial < 500; trial++)
        {
            string s = Random(rng, "ACGTSWNRYacgtswnry", rng.Next(0, 60));
            var sb = new System.Text.StringBuilder(s);
            var inserted = new List<int>();
            int count = rng.Next(1, 6);
            for (int j = 0; j < count; j++)
            {
                int at = rng.Next(0, sb.Length + 1);
                sb.Insert(at, nonAscii[rng.Next(nonAscii.Length)]);
            }
            string x = sb.ToString();
            for (int i = 0; i < x.Length; i++) if (x[i] > 127) inserted.Add(i);

            x.CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Remove)
                .Should().Be(s.CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Remove), x);
            int gcLiteral = s.Count(c => "CGScgs".Contains(c));
            x.CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Ignore).Should().Be((double)gcLiteral / x.Length, x);
            x.CalculateGcFraction().Should().Be(s.CalculateGcFraction(), x);

            x.AsSpan().IndexOfInvalidIupacDna().Should().Be(inserted[0], x);
            x.AsSpan().IndexOfInvalidIupacRna().Should().BeLessThanOrEqualTo(inserted[0], x);
            new IupacDnaSequence(x).GetDegenerateFraction().Should()
                .Be(s.Length == 0 ? 0.0 : s.ToUpperInvariant().Count(c => "RYSWKMBDHVN".Contains(c)) / (double)x.Length, x);

            var k1 = x.AsSpan().CountKmersSpan(1);
            var k1s = s.AsSpan().CountKmersSpan(1);
            foreach (char c in nonAscii)
                k1.GetValueOrDefault(c.ToString()).Should().Be(x.Count(ch => ch == c), x);
            k1.Where(kv => kv.Key[0] < 128).ToDictionary(kv => kv.Key, kv => kv.Value).Should().Equal(k1s, x);

            string replaced = new(x.Select(c => c > 127 ? 'S' : c).ToArray());
            x.AsSpan().HammingDistance(replaced).Should().Be(inserted.Count, x);
            x.AsSpan().HammingDistance(replaced.ToLowerInvariant()).Should()
                .Be(inserted.Count, "ASCII case never counts; only the non-ASCII positions differ");
        }
    }

    /// <summary>
    /// MR (F27/F30): AnalyzeGcContent and the windowed GC driver are invariant under ASCII case swapping in every
    /// mode, and on pure A/C/G/T(U) input every ambiguity mode equals the default.
    /// </summary>
    [Test]
    public void GcAnalysis_CaseSwapAndAcgtModeInvariance()
    {
        var rng = new Random(Seed + 5);
        var modes = new[] { SequenceExtensions.GcAmbiguityMode.Remove, SequenceExtensions.GcAmbiguityMode.Ignore, SequenceExtensions.GcAmbiguityMode.Weighted };
        for (int trial = 0; trial < 300; trial++)
        {
            string s = Random(rng, "ACGTUSWNRYBDHVKMX-acgtuswnrybdhvkmx", rng.Next(0, 120));
            string swapped = Map(s, c => c is >= 'a' and <= 'z' ? (char)(c - 32) : c is >= 'A' and <= 'Z' ? (char)(c + 32) : c);
            int w = rng.Next(1, 30), step = rng.Next(1, 30);
            bool fraction = rng.Next(2) == 0;
            foreach (var m in modes)
            {
                var a = GcSkewCalculator.AnalyzeGcContent(s, w, step, fraction, m);
                var b = GcSkewCalculator.AnalyzeGcContent(swapped, w, step, fraction, m);
                b.OverallGcContent.Should().Be(a.OverallGcContent, s);
                b.WindowedGcContent.Should().Equal(a.WindowedGcContent, s);
                b.WindowedGcSkew.Should().Equal(a.WindowedGcSkew, s);
                b.GcContentVariance.Should().Be(a.GcContentVariance, s);
            }

            string acgt = Random(rng, "ACGTacgt", rng.Next(0, 120));
            var d = GcSkewCalculator.AnalyzeGcContent(acgt, w, step, fraction);
            foreach (var m in modes)
            {
                var r = GcSkewCalculator.AnalyzeGcContent(acgt, w, step, fraction, m);
                r.OverallGcContent.Should().Be(d.OverallGcContent, acgt);
                r.WindowedGcContent.Should().Equal(d.WindowedGcContent, acgt);
                r.GcContentVariance.Should().Be(d.GcContentVariance, acgt);
            }
        }
    }
}
