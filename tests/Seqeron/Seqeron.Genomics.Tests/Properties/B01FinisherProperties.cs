using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the behaviours added by the B01 finisher
/// (docs/Validation/review-2026-09/B01.md, F14–F31): partial-window GC skew (F14), windowed/cumulative
/// AT skew (F15), all skew minimizers/maximizers and Grigoriev windowed origin (F16), circular walk and
/// SkewIT Skew Index (F17), SkewI threshold decision (F31), ASCII-only case folding (F18–F21),
/// definite/degenerate fractions (F22), IndexOfInvalidIupac* (F26), public windowed GC driver (F27) and
/// the AnalyzeGcContent ambiguity mode (F30).
///
/// Every property compares against an independent re-implementation of the cited reference
/// (Biopython 1.88 <c>GC_skew</c>/<c>gc_fraction</c>, Rosalind BA1F brute force, SkewIT <c>skewi.py</c>,
/// scikit-bio <c>definites()/degenerates()</c>, Python <c>Counter</c> over substrings) or asserts an
/// exact algebraic relation. Only exact invariants are asserted (no tolerances).
///
/// Test Units: SEQ-GCSKEW-001, SEQ-ATSKEW-001, SEQ-REPLICATION-001, SEQ-GC-001, SEQ-GC-ANALYSIS-001,
/// SEQ-VALID-001, B01-SWEEP.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Composition")]
public class B01FinisherProperties
{
    // Skew alphabet: both cases, IUPAC codes that are not counted, junk and non-ASCII letters.
    private static readonly char[] SkewAlphabet = "ACGTacgtGGCCNnSsU-ſı".ToCharArray();

    // GC-fraction alphabet: every Biopython gc_fraction class (strong, weak, weighted codes, other) in both cases.
    private static readonly char[] GcAlphabet = "ACGTUSWacgtuswBDHKMNRVXYbdhkmnrvxyZz-.1\u017F\u0131\u212A\u00C5".ToCharArray();

    private static Arbitrary<string> StringOver(char[] alphabet, int minLen = 0, int maxLen = 80) =>
        (from n in Gen.Choose(minLen, maxLen)
         from chars in Gen.Elements(alphabet).ArrayOf(n)
         select new string(chars)).ToArbitrary();

    private static Gen<string> StringGen(char[] alphabet, int minLen, int maxLen) => StringOver(alphabet, minLen, maxLen).Generator;

    private static char UpperAscii(char c) => c is >= 'a' and <= 'z' ? (char)(c - 32) : c;

    private static string UpperAscii(string s) => new(s.Select(UpperAscii).ToArray());

    /// <summary>Biopython 1.88 GC_skew window value: (G−C)/(G+C) case-insensitive, ZeroDivisionError → 0.0.</summary>
    private static double ReferenceSkew(string window, char plus = 'G', char minus = 'C')
    {
        int p = window.Count(c => UpperAscii(c) == plus);
        int m = window.Count(c => UpperAscii(c) == minus);
        return p + m == 0 ? 0.0 : (double)(p - m) / (p + m);
    }

    /// <summary>Window starts 0, step, … (Python <c>range(0, n, step)</c> slicing when partial; complete windows otherwise).</summary>
    private static List<(int Start, int Len)> ReferenceWindows(int n, int w, int step, bool partial)
    {
        var list = new List<(int, int)>();
        for (long i = 0; partial ? i < n : i + w <= n; i += step)
            list.Add(((int)i, (int)Math.Min(w, n - i)));
        return list;
    }

    /// <summary>Rosalind BA1F prefix skew Skew_0..Skew_n (G +1, C −1, case-insensitive).</summary>
    private static int[] PrefixSkew(string s)
    {
        var skew = new int[s.Length + 1];
        for (int i = 0; i < s.Length; i++)
        {
            char c = UpperAscii(s[i]);
            skew[i + 1] = skew[i] + (c == 'G' ? 1 : c == 'C' ? -1 : 0);
        }
        return skew;
    }

    private static Gen<(string S, int W, int Step)> WindowedInput(char[] alphabet, int maxLen = 120) =>
        from s in StringGen(alphabet, 0, maxLen)
        from w in Gen.Choose(1, 25)
        from step in Gen.Choose(1, 25)
        select (s, w, step);

    #region F14 — partial-window GC skew

    /// <summary>
    /// With step == window and includePartialWindow, the window values, starts and ends equal Biopython 1.88
    /// <c>GC_skew(seq, window)</c> (<c>for i in range(0, len(seq), window): s = seq[i:i+window]</c>), and the
    /// position is start + actualLength / 2. The DnaSequence overload agrees on A/C/G/T input.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property PartialWindowedGcSkew_StepEqualsWindow_EqualsBiopythonGcSkew()
    {
        var gen = from s in StringGen(SkewAlphabet, 0, 150)
                  from w in Gen.Choose(1, 30)
                  select (s, w);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var actual = GcSkewCalculator.CalculateWindowedGcSkew(t.s, t.w, t.w, includePartialWindow: true).ToList();
            var expected = ReferenceWindows(t.s.Length, t.w, t.w, partial: true)
                .Select(x => new GcSkewPoint(x.Start + x.Len / 2, ReferenceSkew(t.s.Substring(x.Start, x.Len)), x.Start, x.Start + x.Len - 1))
                .ToList();
            bool ok = actual.SequenceEqual(expected);

            string acgt = new(t.s.Where(c => "ACGT".Contains(c)).ToArray());
            ok &= GcSkewCalculator.CalculateWindowedGcSkew(new DnaSequence(acgt), t.w, t.w, true)
                .SequenceEqual(GcSkewCalculator.CalculateWindowedGcSkew(acgt, t.w, t.w, true));
            return ok.Label($"GC_skew('{t.s}', {t.w}): got [{string.Join(", ", actual.Select(p => p.GcSkew))}], " +
                            $"expected [{string.Join(", ", expected.Select(p => p.GcSkew))}]");
        });
    }

    /// <summary>
    /// Prefix property (any step): the partial result is the default (complete-window) result followed only by
    /// truncated trailing windows, each ending at n − 1, starting on the step lattice and equal to the skew of the
    /// sequence suffix from its start.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property PartialWindowedGcSkew_IsDefaultPlusTrailingTruncatedWindows()
    {
        return Prop.ForAll(WindowedInput(SkewAlphabet).ToArbitrary(), t =>
        {
            var full = GcSkewCalculator.CalculateWindowedGcSkew(t.S, t.W, t.Step).ToList();
            var partial = GcSkewCalculator.CalculateWindowedGcSkew(t.S, t.W, t.Step, includePartialWindow: true).ToList();
            int n = t.S.Length;
            bool ok = partial.Count >= full.Count && partial.Take(full.Count).SequenceEqual(full);
            ok &= partial.Count == (n == 0 ? 0 : (n - 1) / t.Step + 1);
            foreach (var p in partial.Skip(full.Count))
            {
                ok &= p.WindowStart % t.Step == 0
                      && p.WindowStart + t.W > n
                      && p.WindowEnd == n - 1
                      && p.Position == p.WindowStart + (n - p.WindowStart) / 2
                      && p.GcSkew == ReferenceSkew(t.S[p.WindowStart..]);
            }
            return ok.Label($"'{t.S}' w={t.W} step={t.Step}: full {full.Count}, partial {partial.Count}");
        });
    }

    /// <summary>
    /// Cumulative GC skew (both flags): GcSkew equals the windowed values with step == window, and
    /// CumulativeGcSkew is their running sum in order (<c>itertools.accumulate</c>); the default is a prefix of
    /// the partial result (at most one extra, truncated window).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property CumulativeGcSkew_IsRunningSumOfWindowedSkew_BothFlags()
    {
        return Prop.ForAll(WindowedInput(SkewAlphabet).ToArbitrary(), t =>
        {
            bool ok = true;
            foreach (bool partial in new[] { false, true })
            {
                var cum = GcSkewCalculator.CalculateCumulativeGcSkew(t.S, t.W, partial).ToList();
                var win = GcSkewCalculator.CalculateWindowedGcSkew(t.S, t.W, t.W, partial).ToList();
                double sum = 0;
                ok &= cum.Count == win.Count;
                for (int i = 0; ok && i < cum.Count; i++)
                {
                    sum += win[i].GcSkew;
                    ok &= cum[i].Position == win[i].Position && cum[i].GcSkew == win[i].GcSkew && cum[i].CumulativeGcSkew == sum;
                }
            }
            var d = GcSkewCalculator.CalculateCumulativeGcSkew(t.S, t.W).ToList();
            var p2 = GcSkewCalculator.CalculateCumulativeGcSkew(t.S, t.W, includePartialWindow: true).ToList();
            ok &= p2.Take(d.Count).SequenceEqual(d) && p2.Count - d.Count == (t.S.Length % t.W == 0 ? 0 : 1);
            return ok.Label($"cumulative '{t.S}' w={t.W}");
        });
    }

    #endregion

    #region F15 — windowed / cumulative AT skew

    /// <summary>
    /// Windowed and cumulative AT skew equal an independent (A−T)/(A+T) reference over the same windows
    /// (Charneski et al. 2011), with running sum for the cumulative form; the DnaSequence overloads agree.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property AtSkew_WindowedAndCumulative_EqualIndependentReference()
    {
        var gen = from x in WindowedInput(SkewAlphabet)
                  from partial in Gen.Elements(false, true)
                  select (x.S, x.W, x.Step, partial);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var win = GcSkewCalculator.CalculateWindowedAtSkew(t.S, t.W, t.Step, t.partial).ToList();
            var expWin = ReferenceWindows(t.S.Length, t.W, t.Step, t.partial)
                .Select(x => new AtSkewPoint(x.Start + x.Len / 2, ReferenceSkew(t.S.Substring(x.Start, x.Len), 'A', 'T'), x.Start, x.Start + x.Len - 1))
                .ToList();
            bool ok = win.SequenceEqual(expWin);

            var cum = GcSkewCalculator.CalculateCumulativeAtSkew(t.S, t.W, t.partial).ToList();
            double sum = 0;
            var expCum = ReferenceWindows(t.S.Length, t.W, t.W, t.partial).Select(x =>
            {
                double v = ReferenceSkew(t.S.Substring(x.Start, x.Len), 'A', 'T');
                sum += v;
                return new CumulativeAtSkewPoint(x.Start + x.Len / 2, v, sum);
            }).ToList();
            ok &= cum.SequenceEqual(expCum);

            string acgt = new(t.S.Where(c => "ACGT".Contains(c)).ToArray());
            ok &= GcSkewCalculator.CalculateWindowedAtSkew(new DnaSequence(acgt), t.W, t.Step, t.partial)
                      .SequenceEqual(GcSkewCalculator.CalculateWindowedAtSkew(acgt, t.W, t.Step, t.partial))
                  && GcSkewCalculator.CalculateCumulativeAtSkew(new DnaSequence(acgt), t.W, t.partial)
                      .SequenceEqual(GcSkewCalculator.CalculateCumulativeAtSkew(acgt, t.W, t.partial));
            return ok.Label($"AT skew '{t.S}' w={t.W} step={t.Step} partial={t.partial}");
        });
    }

    #endregion

    #region F16/F17 — all extrema, circular walk, windowed origin

    /// <summary>
    /// FindMinimum/MaximumSkewPositions equal the brute-force argmin/argmax sets of Skew_0..Skew_n (Rosalind BA1F,
    /// ascending, never empty); their first elements and extreme values equal PredictReplicationOrigin.
    /// Circular: the same over Skew_0..Skew_{n−1}, matching PredictReplicationOrigin(seq, circular: true).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property SkewExtremaPositions_EqualBruteForce_LinearAndCircular()
    {
        return Prop.ForAll(StringOver(SkewAlphabet, 0, 150), s =>
        {
            int[] skew = PrefixSkew(s);
            bool ok = true;
            foreach (bool circular in new[] { false, true })
            {
                int last = s.Length == 0 ? 0 : circular ? s.Length - 1 : s.Length;
                var values = skew.Take(last + 1).ToArray();
                int min = values.Min(), max = values.Max();
                var expMin = Enumerable.Range(0, values.Length).Where(i => values[i] == min).ToList();
                var expMax = Enumerable.Range(0, values.Length).Where(i => values[i] == max).ToList();

                var actMin = GcSkewCalculator.FindMinimumSkewPositions(s, circular);
                var actMax = GcSkewCalculator.FindMaximumSkewPositions(s, circular);
                var pred = circular ? GcSkewCalculator.PredictReplicationOrigin(s, circular: true)
                                    : GcSkewCalculator.PredictReplicationOrigin(s);
                ok &= actMin.SequenceEqual(expMin) && actMax.SequenceEqual(expMax)
                      && pred.PredictedOrigin == expMin[0] && pred.PredictedTerminus == expMax[0]
                      && pred.OriginSkew == min && pred.TerminusSkew == max
                      && pred.IsSignificant == (max > min);
            }
            return ok.Label($"extrema of '{s}'");
        });
    }

    /// <summary>
    /// Circular walk with total G−C = 0 is rotation-equivariant (Skew'_j = Skew_{(j+r) mod n} − Skew_r):
    /// the minimizer/maximizer sets of the rotation by r are {(p − r) mod n}.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property CircularExtrema_BalancedGc_AreRotationEquivariant()
    {
        var gen = from pairs in Gen.Choose(0, 30)
                  from other in StringGen("ATNat".ToCharArray(), 1, 40)
                  from keys in Gen.Choose(0, int.MaxValue).ArrayOf(2 * pairs + other.Length)
                  from r in Gen.Choose(0, 1000)
                  select (pairs, other, keys, r);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            // Balanced multiset: `pairs` G's and `pairs` C's (mixed case) plus A/T/N filler, shuffled by the keys.
            var chars = Enumerable.Repeat('G', t.pairs).Concat(Enumerable.Repeat('c', t.pairs)).Concat(t.other).ToArray();
            string s = new(chars.Zip(t.keys).OrderBy(z => z.Second).Select(z => z.First).ToArray());
            int n = s.Length, r = t.r % n;
            string rotated = s[r..] + s[..r];

            static List<int> Map(IReadOnlyList<int> ps, int r, int n) => ps.Select(p => ((p - r) % n + n) % n).OrderBy(p => p).ToList();
            bool ok = GcSkewCalculator.FindMinimumSkewPositions(rotated, circular: true)
                          .SequenceEqual(Map(GcSkewCalculator.FindMinimumSkewPositions(s, circular: true), r, n))
                      && GcSkewCalculator.FindMaximumSkewPositions(rotated, circular: true)
                          .SequenceEqual(Map(GcSkewCalculator.FindMaximumSkewPositions(s, circular: true), r, n));
            return ok.Label($"rotation {r} of '{s}'");
        });
    }

    /// <summary>
    /// Windowed PredictReplicationOrigin(seq, w) equals the first argmin/argmax of
    /// <c>numpy.cumsum(GC_skew(seq, w)[:n//w])</c> (Grigoriev 1998; complete windows only), reported at the
    /// window centre start + w/2; no complete window ⇒ (0, 0, 0, 0, false).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property WindowedOrigin_EqualsArgExtremaOfCumsumOfWindowedSkew()
    {
        var gen = from s in StringGen(SkewAlphabet, 0, 200)
                  from w in Gen.Choose(1, 20)
                  select (s, w);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var windows = ReferenceWindows(t.s.Length, t.w, t.w, partial: false);
            var expected = new ReplicationOriginPrediction(0, 0, 0, 0, false);
            if (windows.Count > 0)
            {
                double sum = 0, min = double.PositiveInfinity, max = double.NegativeInfinity;
                int minPos = 0, maxPos = 0;
                foreach (var (start, len) in windows)
                {
                    sum += ReferenceSkew(t.s.Substring(start, len));
                    if (sum < min) { min = sum; minPos = start + t.w / 2; }
                    if (sum > max) { max = sum; maxPos = start + t.w / 2; }
                }
                expected = new ReplicationOriginPrediction(minPos, maxPos, min, max, max > min);
            }
            var actual = GcSkewCalculator.PredictReplicationOrigin(t.s, t.w);
            return (actual == expected).Label($"windowed origin '{t.s}' w={t.w}: {actual} vs {expected}");
        });
    }

    #endregion

    #region F17/F31 — SkewIT Skew Index and threshold decision

    /// <summary>
    /// Literal re-implementation of SkewIT <c>src/skewi.py</c> (master): sign per k-window (partial tail kept),
    /// <c>half_len = round(L/2)</c>, <c>curr_range = round(L·0.04)</c> (Python 3 banker's rounding),
    /// <c>skew += skew[:full_len]</c>, naive <c>sum(skew[i:t])</c> / <c>sum(skew[t:i+full_len])</c> with Python
    /// slice clamping, <c>skewi = max_diff / len(seq) · k</c> capped at 1, nothing reported when max_diff ≤ 0.
    /// Upper-cases first (documented Seqeron deviation; skewi.py counts upper case only).
    /// </summary>
    private static double? SkewIPy(string seq, int k)
    {
        string u = UpperAscii(seq);
        var skew = new List<int>();
        for (int i = 0; i < u.Length; i += k)
        {
            string win = u.Substring(i, Math.Min(k, u.Length - i));
            int g = win.Count(c => c == 'G'), c = win.Count(ch => ch == 'C');
            skew.Add(g - c > 0 ? 1 : g - c < 0 ? -1 : 0);
        }
        int fullLen = skew.Count;
        int halfLen = (int)Math.Round(fullLen / 2.0, MidpointRounding.ToEven);
        skew.AddRange(skew.Take(fullLen).ToList());
        int maxDiff = -1;
        int currRange = (int)Math.Round(fullLen * 0.04, MidpointRounding.ToEven);
        int PySum(int a, int b)
        {
            int lo = Math.Clamp(a, 0, skew.Count), hi = Math.Clamp(b, 0, skew.Count);
            int total = 0;
            for (int j = lo; j < hi; j++) total += skew[j];
            return total;
        }
        for (int i = 0; i < fullLen; i++)
        {
            for (int t = i + halfLen - currRange; t < i + halfLen + currRange; t++)
            {
                int x = PySum(i, t), y = PySum(t, i + fullLen);
                if (Math.Abs(x - y) > maxDiff) maxDiff = Math.Abs(x - y);
            }
        }
        if (maxDiff <= 0) return null;
        double skewi = (double)maxDiff / u.Length * k;
        return skewi > 1 ? 1.0 : skewi;
    }

    /// <summary>CalculateSkewIndex equals the skewi.py re-implementation exactly and is null or in (0, 1].</summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property SkewIndex_EqualsSkewItReimplementation_AndIsBounded()
    {
        var gen = from s in StringGen(SkewAlphabet, 1, 400)
                  from k in Gen.Choose(1, 12)
                  select (s, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double? actual = GcSkewCalculator.CalculateSkewIndex(t.s, t.k);
            double? expected = SkewIPy(t.s, t.k);
            bool ok = actual == expected && (actual is null || (actual > 0 && actual <= 1));
            string acgt = new(t.s.Where(c => "ACGT".Contains(c)).ToArray());
            ok &= acgt.Length == 0 || GcSkewCalculator.CalculateSkewIndex(new DnaSequence(acgt), t.k) == SkewIPy(acgt, t.k);
            return ok.Label($"SkewI('{t.s}', k={t.k}) = {actual?.ToString("R") ?? "null"}, skewi.py {expected?.ToString("R") ?? "null"}");
        });
    }

    /// <summary>
    /// SkewIT decision consistency (F31): IsSkewIBelowThreshold = (SkewI &lt; threshold) or null when SkewI is null;
    /// IsSkewIBelowGenusThreshold = the same with the table's threshold, null for a genus absent from the table;
    /// lookup is case-insensitive and accepts the <c>g__</c> prefix.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property SkewIThresholdDecision_IsConsistentWithSkewIndex()
    {
        var gen = from s in StringGen(SkewAlphabet, 0, 200)
                  from k in Gen.Choose(1, 10)
                  from th in Gen.Choose(-300, 1100)
                  from known in Gen.Elements(true, false)
                  select (s, k, th: th / 1000.0, known);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double? skewI = GcSkewCalculator.CalculateSkewIndex(t.s, t.k);
            bool? expected = skewI.HasValue ? skewI.Value < t.th : null;
            var table = new Dictionary<string, double> { ["Escherichia"] = t.th, ["Streptomyces"] = 0.046 };
            string genus = t.known ? "g__eSCHERICHIA" : "Bacillus";
            bool ok = GcSkewCalculator.IsSkewIBelowThreshold(t.s, t.th, t.k) == expected
                      && GcSkewCalculator.IsSkewIBelowGenusThreshold(t.s, genus, table, t.k) == (t.known ? expected : null)
                      && GcSkewCalculator.TryGetSkewIThreshold(table, genus, out double got) == t.known
                      && (!t.known || got == t.th);
            return ok.Label($"decision '{t.s}' k={t.k} th={t.th} known={t.known}: SkewI {skewI}");
        });
    }

    #endregion

    #region F18/F19/F21 — ASCII-only folding

    /// <summary>
    /// Independent Biopython 1.88 <c>gc_fraction</c> reference: literal counts of "CGScgs" (numerator),
    /// "ATWUatwu" (+ GC → Remove denominator), <c>_gc_values</c> for "BDHKMNRVXY" (both cases, Weighted),
    /// full length for Ignore/Weighted. Weighted contributions are added per character in sequence order.
    /// </summary>
    private static double BiopythonGcFraction(string s, SequenceExtensions.GcAmbiguityMode mode)
    {
        if (s.Length == 0) return 0;
        double gc = 0;
        int strongWeak = 0;
        foreach (char c in s)
        {
            if ("CGScgs".Contains(c)) { gc += 1.0; strongWeak++; }
            else if ("ATWUatwu".Contains(c)) strongWeak++;
            else if (mode == SequenceExtensions.GcAmbiguityMode.Weighted)
            {
                gc += c switch
                {
                    'V' or 'v' or 'B' or 'b' => 2.0 / 3.0,
                    'H' or 'h' or 'D' or 'd' => 1.0 / 3.0,
                    'M' or 'm' or 'R' or 'r' or 'Y' or 'y' or 'K' or 'k' or 'X' or 'x' or 'N' or 'n' => 0.5,
                    _ => 0.0,
                };
            }
        }
        double len = mode == SequenceExtensions.GcAmbiguityMode.Remove ? strongWeak : s.Length;
        return len > 0 ? gc / len : 0;
    }

    private static readonly SequenceExtensions.GcAmbiguityMode[] AllModes =
        { SequenceExtensions.GcAmbiguityMode.Remove, SequenceExtensions.GcAmbiguityMode.Ignore, SequenceExtensions.GcAmbiguityMode.Weighted };

    /// <summary>CalculateGcFraction(mode) (string and span) equals the Biopython reference on input with non-ASCII letters.</summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property GcFractionModes_EqualBiopythonReference_NonAsciiNeverCounted()
    {
        return Prop.ForAll(StringOver(GcAlphabet, 0, 80), s =>
        {
            bool ok = AllModes.All(m =>
                s.CalculateGcFraction(m) == BiopythonGcFraction(s, m)
                && s.AsSpan().CalculateGcFraction(m) == BiopythonGcFraction(s, m));
            return ok.Label($"gc_fraction('{s}')");
        });
    }

    /// <summary>
    /// ExpandCode folds ASCII only: lower-case ASCII equals upper-case, and every non-ASCII char (incl. 'ſ', 'ı',
    /// Kelvin 'K') expands to { 'N' } (not a code; Biopython <c>ambiguous_dna_values</c> has no such key).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property ExpandCode_FoldsAsciiOnly()
    {
        return Prop.ForAll(Gen.Choose(0, 0xFFFF).Select(i => (char)i).ToArbitrary(), c =>
        {
            char[] got = IupacDnaSequence.ExpandCode(c);
            bool ok = c > 127
                ? got.SequenceEqual(new[] { 'N' })
                : got.SequenceEqual(IupacDnaSequence.ExpandCode(UpperAscii(c)));
            return ok.Label($"ExpandCode(U+{(int)c:X4}) = [{new string(got)}]");
        });
    }

    /// <summary>
    /// Span HammingDistance equals a char-wise ASCII-folded reference (scipy <c>hamming(list(a.upper()), …)·n</c>
    /// for ASCII); 'ſ' vs 'S'/'s' and Kelvin vs 'k' are mismatches; symmetric.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property HammingDistance_EqualsAsciiFoldedReference()
    {
        var alphabet = "ACGTNSsacgtn\u017F\u212Ak\u0131".ToCharArray();
        var gen = from n in Gen.Choose(0, 60)
                  from a in Gen.Elements(alphabet).ArrayOf(n)
                  from b in Gen.Elements(alphabet).ArrayOf(n)
                  select (new string(a), new string(b));
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            int expected = t.Item1.Zip(t.Item2).Count(z => UpperAscii(z.First) != UpperAscii(z.Second));
            int d = t.Item1.AsSpan().HammingDistance(t.Item2);
            bool ok = d == expected && t.Item2.AsSpan().HammingDistance(t.Item1) == d;
            return ok.Label($"Hamming('{t.Item1}', '{t.Item2}') = {d}, expected {expected}");
        });
    }

    #endregion

    #region F20 — CountKmersSpan

    /// <summary>
    /// CountKmersSpan equals a substring-based reference (Python <c>Counter(s.upper()[i:i+k])</c> with ASCII-only
    /// upper-casing) in keys, counts and first-occurrence order — across the stackalloc/ArrayPool threshold — and
    /// never modifies its input.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property CountKmersSpan_EqualsSubstringReference_IncludingOrder()
    {
        var gen = from s in Gen.OneOf(StringGen("ACGTNacgtnſ".ToCharArray(), 0, 60), StringGen("ACgt".ToCharArray(), 240, 600))
                  from k in Gen.Choose(1, 10)
                  select (s, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string before = new(t.s.ToCharArray());
            var actual = t.s.AsSpan().CountKmersSpan(t.k).ToList();
            string u = UpperAscii(t.s);
            var expected = new List<KeyValuePair<string, int>>();
            var index = new Dictionary<string, int>();
            for (int i = 0; i + t.k <= u.Length; i++)
            {
                string key = u.Substring(i, t.k);
                if (index.TryGetValue(key, out int at)) expected[at] = new(key, expected[at].Value + 1);
                else { index[key] = expected.Count; expected.Add(new(key, 1)); }
            }
            bool ok = actual.SequenceEqual(expected) && t.s == before;
            return ok.Label($"CountKmersSpan(len {t.s.Length}, k={t.k}): {actual.Count} keys vs {expected.Count}");
        });
    }

    #endregion

    #region F22 — definite / degenerate fractions

    /// <summary>
    /// GetAmbiguityLevel = #ACGT / L (scikit-bio <c>definites().mean()</c>), GetDegenerateFraction =
    /// #{RYSWKMBDHVN} / L (<c>degenerates().mean()</c>) over the ASCII-upper-cased container; their counts never
    /// exceed L, so the sum ≤ 1 with equality iff every symbol is one of the 15 codes. Empty: 1.0 and 0.0.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property DefiniteAndDegenerateFractions_MatchCounts()
    {
        return Prop.ForAll(StringOver("ACGTUacgtuRYSWKMBDHVNrysn-.Xſ".ToCharArray(), 0, 60), s =>
        {
            var seq = new IupacDnaSequence(s);
            string u = UpperAscii(s);
            int def = u.Count(c => "ACGT".Contains(c));
            int deg = u.Count(c => "RYSWKMBDHVN".Contains(c));
            double level = seq.GetAmbiguityLevel(), degenerate = seq.GetDegenerateFraction();
            bool ok = s.Length == 0
                ? level == 1.0 && degenerate == 0.0
                : level == def / (double)s.Length && degenerate == deg / (double)s.Length
                  && def + deg <= s.Length
                  && ((def + deg == s.Length) == u.All(IupacHelper.IsNucleotideCode));
            return ok.Label($"'{s}': definite {level} (#{def}), degenerate {degenerate} (#{deg})");
        });
    }

    #endregion

    #region F26/F27/F30 — IndexOfInvalidIupac*, windowed GC driver, AnalyzeGcContent ambiguity mode

    /// <summary>IndexOfInvalidIupacDna/Rna is the first index outside the 15 ASCII codes (T↔U per alphabet), −1 ⇔ IsValidIupac*.</summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property IndexOfInvalidIupac_IsFirstOffender_AndAgreesWithIsValid()
    {
        return Prop.ForAll(StringOver(GcAlphabet, 0, 40), s =>
        {
            int refDna = Array.FindIndex(s.ToCharArray(), c => !"ACGTRYSWKMBDHVN".Contains(UpperAscii(c)));
            int refRna = Array.FindIndex(s.ToCharArray(), c => !"ACGURYSWKMBDHVN".Contains(UpperAscii(c)));
            int iDna = s.AsSpan().IndexOfInvalidIupacDna(), iRna = s.AsSpan().IndexOfInvalidIupacRna();
            bool ok = iDna == refDna && iRna == refRna
                      && (iDna < 0) == s.AsSpan().IsValidIupacDna() && (iRna < 0) == s.AsSpan().IsValidIupacRna();
            return ok.Label($"'{s}': DNA {iDna} (ref {refDna}), RNA {iRna} (ref {refRna})");
        });
    }

    /// <summary>
    /// CalculateWindowedGcContent (default and every mode, percent and fraction) equals per-window
    /// CalculateGcFraction over complete windows (Biopython <c>gc_fraction(s[i:i+w])</c> for
    /// <c>i in range(0, n−w+1, step)</c>), positioned at start + w/2.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property WindowedGcContent_EqualsPerWindowGcFraction()
    {
        var gen = from x in WindowedInput(GcAlphabet)
                  from fraction in Gen.Elements(false, true)
                  select (x.S, x.W, x.Step, fraction);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var windows = ReferenceWindows(t.S.Length, t.W, t.Step, partial: false);
            double scale = t.fraction ? 1.0 : 100.0;
            GcContentPoint Point(int start, double f) =>
                new(start + t.W / 2, t.fraction ? f : f * scale, start, start + t.W - 1);

            bool ok = GcSkewCalculator.CalculateWindowedGcContent(t.S, t.W, t.Step, t.fraction)
                .SequenceEqual(windows.Select(w => Point(w.Start, t.S.Substring(w.Start, t.W).AsSpan().CalculateGcFraction())));
            foreach (var m in AllModes)
            {
                ok &= GcSkewCalculator.CalculateWindowedGcContent(t.S, t.W, t.Step, t.fraction, m)
                    .SequenceEqual(windows.Select(w => Point(w.Start, BiopythonGcFraction(t.S.Substring(w.Start, t.W), m))));
            }
            return ok.Label($"windowed GC '{t.S}' w={t.W} step={t.Step} fraction={t.fraction}");
        });
    }

    /// <summary>
    /// AnalyzeGcContent(…, mode): OverallGcContent = CalculateGcFraction(mode) (×100 unless fraction), windows =
    /// CalculateWindowedGcContent(…, mode), variance = numpy.var of those windows; the skew fields are identical
    /// to the default overload (the mode affects GC content only).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property AnalyzeGcContent_AmbiguityMode_EqualsGcFractionOverallAndPerWindow()
    {
        var gen = from x in WindowedInput(GcAlphabet)
                  from fraction in Gen.Elements(false, true)
                  from m in Gen.Elements(AllModes)
                  select (x.S, x.W, x.Step, fraction, m);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var r = GcSkewCalculator.AnalyzeGcContent(t.S, t.W, t.Step, t.fraction, t.m);
            var d = GcSkewCalculator.AnalyzeGcContent(t.S, t.W, t.Step, t.fraction);
            double overall = t.S.CalculateGcFraction(t.m);
            var windows = GcSkewCalculator.CalculateWindowedGcContent(t.S, t.W, t.Step, t.fraction, t.m).ToList();
            bool ok = r.OverallGcContent == (t.fraction ? overall : overall * 100.0)
                      && r.WindowedGcContent.SequenceEqual(windows)
                      && r.GcContentVariance == StatisticsHelper.PopulationVariance(windows.Select(w => w.GcContent).ToList())
                      && r.OverallGcSkew == d.OverallGcSkew && r.OverallAtSkew == d.OverallAtSkew
                      && r.GcSkewVariance == d.GcSkewVariance && r.WindowedGcSkew.SequenceEqual(d.WindowedGcSkew)
                      && r.SequenceLength == t.S.Length;
            return ok.Label($"AnalyzeGcContent('{t.S}', {t.W}, {t.Step}, {t.fraction}, {t.m})");
        });
    }

    #endregion
}
