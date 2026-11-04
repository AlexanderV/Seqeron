namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Heavy-tier fuzz tests for the entry points added by the B01 finisher
/// (docs/Validation/review-2026-09/B01.md, F14–F31). Random UTF-16 garbage (incl. surrogate halves, '\0',
/// U+017F 'ſ', U+212A Kelvin) and random sizes (incl. 0, negative, int.MaxValue) are fed to every new
/// surface; each call must return a value satisfying its documented contract or throw the documented
/// exception (ArgumentException family; FormatException for a malformed SkewI table) — never hang and never
/// leak another runtime exception.
///
/// Test Units: SEQ-GCSKEW-001, SEQ-ATSKEW-001, SEQ-REPLICATION-001, SEQ-GC-001, SEQ-GC-ANALYSIS-001,
/// SEQ-VALID-001, B01-SWEEP.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Composition")]
public class B01FinisherFuzzTests
{
    private const int Seed = 20261009;

    private static string Garbage(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
        {
            chars[i] = rng.Next(6) switch
            {
                0 => (char)rng.Next(0, 0x10000),
                1 => (char)rng.Next(0, 128),
                2 => "\u017F\u212A\u0131\u0000\uFFFF\uD800\uDC00"[rng.Next(7)],
                _ => "ACGTUNacgtunRYSWKMBDHVXx-."[rng.Next(26)],
            };
        }
        return new string(chars);
    }

    private static int Size(Random rng) => rng.Next(10) switch
    {
        0 => rng.Next(-3, 1),          // invalid: ≤ 0
        1 => int.MaxValue,
        2 => rng.Next(1, 5000),
        _ => rng.Next(1, 40),
    };

    /// <summary>Runs <paramref name="act"/>; a throw is allowed only when <paramref name="invalid"/> and of the given type.</summary>
    private static void Guarded<TException>(bool invalid, Action act, string context) where TException : Exception
    {
        if (invalid)
            act.Should().Throw<TException>(context);
        else
            act.Should().NotThrow(context);
    }

    /// <summary>
    /// Windowed / cumulative GC and AT skew with includePartialWindow on garbage and random sizes: sizes &lt; 1
    /// throw ArgumentOutOfRangeException eagerly; otherwise values are in [−1, 1], windows lie inside the input,
    /// cumulative values are finite, and a huge step terminates after the first window.
    /// </summary>
    [Test]
    [CancelAfter(120_000)]
    public void SkewWindows_Garbage_ContractHolds()
    {
        var rng = new Random(Seed);
        for (int trial = 0; trial < 1500; trial++)
        {
            string s = Garbage(rng, rng.Next(0, 120));
            int w = Size(rng), step = Size(rng);
            bool partial = rng.Next(2) == 0;
            string ctx = $"len={s.Length} w={w} step={step} partial={partial}";
            bool invalid = w < 1 || step < 1;

            Guarded<ArgumentOutOfRangeException>(invalid, () =>
            {
                var gc = GcSkewCalculator.CalculateWindowedGcSkew(s, w, step, partial).ToList();
                var at = GcSkewCalculator.CalculateWindowedAtSkew(s, w, step, partial).ToList();
                gc.Should().HaveSameCount(at);
                foreach (var p in gc)
                {
                    p.GcSkew.Should().BeInRange(-1, 1);
                    p.WindowStart.Should().BeInRange(0, s.Length - 1);
                    p.WindowEnd.Should().BeInRange(p.WindowStart, s.Length - 1);
                    p.Position.Should().BeInRange(p.WindowStart, p.WindowEnd + 1);
                }
                at.Should().OnlyContain(p => p.AtSkew >= -1 && p.AtSkew <= 1);
            }, ctx);

            Guarded<ArgumentOutOfRangeException>(w < 1, () =>
            {
                var cg = GcSkewCalculator.CalculateCumulativeGcSkew(s, w, partial).ToList();
                var ca = GcSkewCalculator.CalculateCumulativeAtSkew(s, w, partial).ToList();
                cg.Should().OnlyContain(p => double.IsFinite(p.CumulativeGcSkew) && Math.Abs(p.CumulativeGcSkew) <= s.Length);
                ca.Should().OnlyContain(p => double.IsFinite(p.CumulativeAtSkew) && Math.Abs(p.CumulativeAtSkew) <= s.Length);
            }, ctx);
        }
    }

    /// <summary>
    /// Extremum positions, circular and windowed origin on garbage: never throw for any string (incl. null),
    /// positions ascending and inside [0, n] (linear) / [0, n−1] (circular), first element = predicted origin;
    /// windowed origin with w &lt; 1 throws ArgumentOutOfRangeException.
    /// </summary>
    [Test]
    [CancelAfter(120_000)]
    public void ReplicationOrigin_Garbage_ContractHolds()
    {
        var rng = new Random(Seed + 1);
        for (int trial = 0; trial < 1500; trial++)
        {
            string? s = rng.Next(20) == 0 ? null : Garbage(rng, rng.Next(0, 150));
            int n = s?.Length ?? 0;
            foreach (bool circular in new[] { false, true })
            {
                var min = GcSkewCalculator.FindMinimumSkewPositions(s!, circular);
                var max = GcSkewCalculator.FindMaximumSkewPositions(s!, circular);
                int upper = n == 0 ? 0 : circular ? n - 1 : n;
                min.Should().NotBeEmpty().And.BeInAscendingOrder().And.OnlyHaveUniqueItems();
                max.Should().NotBeEmpty().And.BeInAscendingOrder().And.OnlyHaveUniqueItems();
                min.Should().OnlyContain(p => p >= 0 && p <= upper);
                max.Should().OnlyContain(p => p >= 0 && p <= upper);
                var pred = GcSkewCalculator.PredictReplicationOrigin(s!, circular);
                pred.PredictedOrigin.Should().Be(min[0]);
                pred.PredictedTerminus.Should().Be(max[0]);
                pred.OriginSkew.Should().BeLessThanOrEqualTo(0);
                pred.TerminusSkew.Should().BeGreaterThanOrEqualTo(0);
            }

            int w = Size(rng);
            Guarded<ArgumentOutOfRangeException>(w < 1, () =>
            {
                var r = GcSkewCalculator.PredictReplicationOrigin(s!, w);
                r.OriginSkew.Should().BeLessThanOrEqualTo(r.TerminusSkew);
                r.IsSignificant.Should().Be(r.TerminusSkew > r.OriginSkew);
                if (n >= w)
                {
                    r.PredictedOrigin.Should().BeInRange(0, n - 1);
                    r.PredictedTerminus.Should().BeInRange(0, n - 1);
                }
                else
                {
                    r.Should().Be(new ReplicationOriginPrediction(0, 0, 0, 0, false));
                }
            }, $"len={n} w={w}");
        }
    }

    /// <summary>
    /// Skew Index and threshold decision on garbage: SkewI is null or in (0, 1]; k &lt; 1 throws
    /// ArgumentOutOfRangeException; a non-finite threshold throws ArgumentOutOfRangeException; otherwise the
    /// decision is null exactly when SkewI is null.
    /// </summary>
    [Test]
    [CancelAfter(180_000)]
    public void SkewIndex_Garbage_ContractHolds()
    {
        var rng = new Random(Seed + 2);
        double[] specials = { double.NaN, double.PositiveInfinity, double.NegativeInfinity };
        for (int trial = 0; trial < 800; trial++)
        {
            string? s = rng.Next(20) == 0 ? null : Garbage(rng, rng.Next(0, 300));
            int k = rng.Next(10) switch { 0 => rng.Next(-3, 1), 1 => int.MaxValue, _ => rng.Next(1, 15) };
            double th = rng.Next(5) == 0 ? specials[rng.Next(3)] : rng.NextDouble() * 1.4 - 0.3;

            double? skewI = null;
            Guarded<ArgumentOutOfRangeException>(k < 1, () => skewI = GcSkewCalculator.CalculateSkewIndex(s!, k), $"k={k}");
            if (k < 1) continue;
            (skewI is null || (skewI > 0 && skewI <= 1)).Should().BeTrue($"SkewI {skewI}");

            bool? decision = null;
            Guarded<ArgumentOutOfRangeException>(!double.IsFinite(th),
                () => decision = GcSkewCalculator.IsSkewIBelowThreshold(s!, th, k), $"th={th}");
            if (double.IsFinite(th))
                decision.Should().Be(skewI.HasValue ? skewI.Value < th : null);
        }
    }

    /// <summary>
    /// SkewI threshold table parser on garbage text: returns a table (values finite, keys without g__) or throws
    /// FormatException (documented: malformed value, empty/duplicate genus); lookups on the result never throw.
    /// </summary>
    [Test]
    [CancelAfter(60_000)]
    public void SkewIThresholdTable_Garbage_ParsesOrThrowsFormatException()
    {
        var rng = new Random(Seed + 3);
        string[] cells = { "g__Escherichia", "Bacillus", "", " ", "0.711", "-0.2", "1e3", "NaN", "abc", "\u0000", "g__", "934", "0.0688" };
        for (int trial = 0; trial < 2000; trial++)
        {
            var sb = new System.Text.StringBuilder();
            if (rng.Next(2) == 0) sb.Append("Genus\tNum_Genomes\tMean\tSTDEV \tThreshold\r\n");
            int lines = rng.Next(0, 6);
            for (int l = 0; l < lines; l++)
            {
                int cols = rng.Next(0, 7);
                for (int c = 0; c < cols; c++)
                {
                    if (c > 0) sb.Append('\t');
                    sb.Append(rng.Next(4) == 0 ? Garbage(rng, rng.Next(0, 8)) : cells[rng.Next(cells.Length)]);
                }
                sb.Append(rng.Next(2) == 0 ? "\r\n" : "\n");
            }
            string text = sb.ToString();

            IReadOnlyDictionary<string, double>? table = null;
            try
            {
                table = GcSkewCalculator.ParseSkewIGenusThresholds(new StringReader(text));
            }
            catch (FormatException)
            {
                continue;
            }

            table.Values.Should().OnlyContain(v => double.IsFinite(v), text);
            table.Keys.Should().OnlyContain(key => key.Length > 0 && !key.StartsWith("g__", StringComparison.OrdinalIgnoreCase), text);
            foreach (var key in table.Keys)
            {
                GcSkewCalculator.TryGetSkewIThreshold(table, "g__" + new string(key.Select(c => c is >= 'a' and <= 'z' ? (char)(c - 32) : c).ToArray()), out double v).Should().BeTrue(text);
                v.Should().Be(table[key]);
            }
            var act = () => GcSkewCalculator.TryGetSkewIThreshold(table, Garbage(rng, rng.Next(0, 10)), out _);
            act.Should().NotThrow(text);
        }

        var nullTable = () => GcSkewCalculator.TryGetSkewIThreshold(null!, "Escherichia", out _);
        nullTable.Should().Throw<ArgumentNullException>();
        var nullReader = () => GcSkewCalculator.ParseSkewIGenusThresholds(null!);
        nullReader.Should().Throw<ArgumentNullException>();
    }

    /// <summary>
    /// ASCII-only folding surfaces on garbage: CalculateGcFraction(mode) in [0, 1] (string = span); ExpandCode
    /// always non-empty and only A/C/G/T/N; Hamming in [0, n], 0 on self, ArgumentException on unequal lengths;
    /// CountKmersSpan counts sum to n − k + 1 with k ≤ 0 → ArgumentOutOfRangeException; definite + degenerate
    /// fractions in [0, 1] with sum ≤ 1; IndexOfInvalidIupac* in [−1, n) and consistent with IsValidIupac*.
    /// </summary>
    [Test]
    [CancelAfter(120_000)]
    public void AsciiFoldingSurfaces_Garbage_ContractHolds()
    {
        var rng = new Random(Seed + 4);
        var modes = new[] { SequenceExtensions.GcAmbiguityMode.Remove, SequenceExtensions.GcAmbiguityMode.Ignore, SequenceExtensions.GcAmbiguityMode.Weighted };
        for (int trial = 0; trial < 2000; trial++)
        {
            string s = Garbage(rng, rng.Next(0, 80));
            foreach (var m in modes)
            {
                double f = s.CalculateGcFraction(m);
                f.Should().BeInRange(0, 1, s);
                s.AsSpan().CalculateGcFraction(m).Should().Be(f);
            }

            char ch = s.Length > 0 ? s[rng.Next(s.Length)] : (char)rng.Next(0, 0x10000);
            IupacDnaSequence.ExpandCode(ch).Should().NotBeEmpty().And.OnlyContain(b => "ACGTN".Contains(b));

            string other = Garbage(rng, rng.Next(2) == 0 ? s.Length : rng.Next(0, 80));
            if (other.Length == s.Length)
            {
                s.AsSpan().HammingDistance(other).Should().BeInRange(0, s.Length);
                s.AsSpan().HammingDistance(s).Should().Be(0);
            }
            else
            {
                var act = () => s.AsSpan().HammingDistance(other);
                act.Should().Throw<ArgumentException>();
            }

            int k = rng.Next(-2, 12);
            Guarded<ArgumentOutOfRangeException>(k <= 0, () =>
            {
                var counts = s.AsSpan().CountKmersSpan(k);
                counts.Values.Sum().Should().Be(Math.Max(0, s.Length - k + 1));
                counts.Keys.Should().OnlyContain(key => key.Length == k && key.All(c => c < 'a' || c > 'z'));
            }, $"k={k}");

            var iupac = new IupacDnaSequence(s);
            double def = iupac.GetAmbiguityLevel(), deg = iupac.GetDegenerateFraction();
            def.Should().BeInRange(0, 1);
            deg.Should().BeInRange(0, 1);
            if (s.Length > 0)
                (Math.Round(def * s.Length) + Math.Round(deg * s.Length)).Should().BeLessThanOrEqualTo(s.Length);

            int iDna = s.AsSpan().IndexOfInvalidIupacDna(), iRna = s.AsSpan().IndexOfInvalidIupacRna();
            iDna.Should().BeInRange(-1, s.Length - 1);
            iRna.Should().BeInRange(-1, s.Length - 1);
            (iDna < 0).Should().Be(s.AsSpan().IsValidIupacDna());
            (iRna < 0).Should().Be(s.AsSpan().IsValidIupacRna());
        }
    }

    /// <summary>
    /// AnalyzeGcContent(mode) and CalculateWindowedGcContent (all overloads) on garbage with random sizes: sizes
    /// &lt; 1 throw ArgumentOutOfRangeException eagerly; otherwise GC values lie in [0, 1] (fraction) or [0, 100],
    /// windows are complete and inside the input, variances are finite and ≥ 0.
    /// </summary>
    [Test]
    [CancelAfter(120_000)]
    public void GcContentDrivers_Garbage_ContractHolds()
    {
        var rng = new Random(Seed + 5);
        var modes = new[] { SequenceExtensions.GcAmbiguityMode.Remove, SequenceExtensions.GcAmbiguityMode.Ignore, SequenceExtensions.GcAmbiguityMode.Weighted };
        for (int trial = 0; trial < 1500; trial++)
        {
            string? s = rng.Next(20) == 0 ? null : Garbage(rng, rng.Next(0, 120));
            int w = Size(rng), step = Size(rng);
            bool fraction = rng.Next(2) == 0;
            var m = modes[rng.Next(3)];
            double hi = fraction ? 1.0 : 100.0;
            bool invalid = w < 1 || step < 1;
            string ctx = $"len={s?.Length} w={w} step={step} fraction={fraction} mode={m}";

            Guarded<ArgumentOutOfRangeException>(invalid, () =>
            {
                var r = GcSkewCalculator.AnalyzeGcContent(s!, w, step, fraction, m);
                r.OverallGcContent.Should().BeInRange(0, hi);
                r.GcContentVariance.Should().BeGreaterThanOrEqualTo(0);
                r.GcSkewVariance.Should().BeGreaterThanOrEqualTo(0);
                double.IsFinite(r.GcContentVariance).Should().BeTrue();
                r.WindowedGcContent.Should().HaveSameCount(r.WindowedGcSkew);
                foreach (var p in r.WindowedGcContent)
                {
                    p.GcContent.Should().BeInRange(0, hi);
                    p.WindowEnd.Should().Be(p.WindowStart + w - 1);
                    p.WindowEnd.Should().BeLessThan(s!.Length);
                }
            }, ctx);

            Guarded<ArgumentOutOfRangeException>(invalid, () =>
            {
                var a = GcSkewCalculator.CalculateWindowedGcContent(s!, w, step, fraction).ToList();
                var b = GcSkewCalculator.CalculateWindowedGcContent(s!, w, step, fraction, m).ToList();
                a.Should().HaveSameCount(b);
                a.Concat(b).Should().OnlyContain(p => p.GcContent >= 0 && p.GcContent <= hi && p.WindowEnd - p.WindowStart == w - 1);
            }, ctx);
        }

        var nullDna = () => GcSkewCalculator.CalculateWindowedGcContent((DnaSequence)null!, 10, 5);
        nullDna.Should().Throw<ArgumentNullException>();
    }
}
