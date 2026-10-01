using System.Numerics;
using FsCheck;
using FsCheck.Fluent;
using Seqeron.Genomics.Tests.Unit.Analysis;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Heavy-tier property tests for the behaviour added by the B04 completeness audit (docs/Validation/review-2026-09/B04.md,
/// fixes F34–F65) in <see cref="SequenceComplexity"/> and <see cref="RepeatFinder"/>:
/// DUST word size / sdust linker + soft mask + per-ACGT-run string path (F34–F36), the dustmasker engine (F53),
/// longdust (F34), fixed-alphabet LC (F37), windowed LC / string LCR / BBDuk entropy (F38, F39, F55), k-mer entropy
/// corrections + normalisation + Digamma (F54), the TRF parameter API, EntropyTrf, mask and flanks (F40–F42),
/// caller-supplied apparent-size tables (F56), the TRF .dat / HTML / alignment-page formatters (F57, F64, F65),
/// degenerate repeats + Vmatch reporting (F47, F58, F59) and the MISA scan mode / unit sizes / string summary
/// (F48, F51, F61). Every property is an exact invariant from the sourced definition or a brute-force oracle;
/// expected values are never derived from the code under test. Sizes are kept small so the fixture stays fast.
///
/// Test Units: SEQ-COMPLEX-DUST-001, SEQ-COMPLEX-001, SEQ-COMPLEX-WINDOW-001, SEQ-COMPLEX-KMER-001, REP-APPROX-001,
/// REP-DIRECT-001, REP-STR-001, REP-TANDEM-001.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Complexity")]
public class B04AuditProperties
{
    #region Generators / helpers

    private static readonly Gen<int> Seeds = Gen.Choose(0, int.MaxValue - 1);

    /// <summary>
    /// Low-complexity-rich DNA built from segments: random ACGT, short tandem units with point mutations,
    /// and (when <paramref name="extra"/> is non-empty) runs of the extra symbols.
    /// </summary>
    internal static string Segmented(Random rng, int length, string extra = "")
    {
        var sb = new System.Text.StringBuilder(length + 64);
        while (sb.Length < length)
        {
            int kind = rng.Next(extra.Length > 0 ? 4 : 3);
            if (kind == 0)
            {
                int len = rng.Next(3, 40);
                for (int i = 0; i < len; i++) sb.Append("ACGT"[rng.Next(4)]);
            }
            else if (kind == 3)
            {
                int len = rng.Next(1, 12);
                for (int i = 0; i < len; i++) sb.Append(extra[rng.Next(extra.Length)]);
            }
            else
            {
                string unit = new(Enumerable.Range(0, rng.Next(1, 7)).Select(_ => "ACGT"[rng.Next(4)]).ToArray());
                int copies = rng.Next(3, 15);
                for (int c = 0; c < copies; c++)
                    foreach (char ch in unit)
                        sb.Append(rng.Next(25) == 0 ? "ACGT"[rng.Next(4)] : ch);
            }
        }
        return sb.ToString(0, length);
    }

    private static string RandomOver(Random rng, string alphabet, int length) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static bool IsAcgt(char c) => c is 'A' or 'C' or 'G' or 'T';

    private static char Complement(char c) => c switch { 'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => '?' };

    private static string ReverseComplement(string s) => new(s.Reverse().Select(Complement).ToArray());

    private static HashSet<int> Covered(IEnumerable<(int Start, int End)> intervals)
    {
        var set = new HashSet<int>();
        foreach (var (s, e) in intervals)
            for (int i = s; i < e; i++) set.Add(i);
        return set;
    }

    private static bool WellFormed(IReadOnlyList<(int Start, int End)> iv, int n)
    {
        for (int i = 0; i < iv.Count; i++)
        {
            if (iv[i].Start < 0 || iv[i].End > n || iv[i].Start >= iv[i].End) return false;
            if (i > 0 && iv[i - 1].End > iv[i].Start) return false;
        }
        return true;
    }

    /// <summary>Unit-cost edit distance where only equal ACGT symbols match (Vmatch: a wildcard never matches).</summary>
    private static int WildcardEditDistance(string a, string b)
    {
        var prev = new int[b.Length + 1];
        var cur = new int[b.Length + 1];
        for (int j = 0; j <= b.Length; j++) prev[j] = j;
        for (int i = 1; i <= a.Length; i++)
        {
            cur[0] = i;
            for (int j = 1; j <= b.Length; j++)
            {
                bool match = IsAcgt(a[i - 1]) && a[i - 1] == b[j - 1];
                cur[j] = Math.Min(prev[j - 1] + (match ? 0 : 1), Math.Min(prev[j], cur[j - 1]) + 1);
            }
            (prev, cur) = (cur, prev);
        }
        return prev[b.Length];
    }

    private static int WildcardHamming(string a, string b) =>
        a.Length != b.Length ? int.MaxValue : a.Zip(b).Count(p => !(IsAcgt(p.First) && p.First == p.Second));

    #endregion

    #region DUST: sdust linker / soft mask / string path (F34–F36) and dustmasker engine (F53)

    /// <summary>
    /// F35/F36/F53: for both engines and random (W, level, linker) the intervals are sorted, disjoint, half-open and in
    /// range; the hard mask (mask char 'X') marks exactly their union; the soft mask lower-cases exactly the same
    /// positions and upper-cases back to the input; the DnaSequence overload equals the string one on ACGT input.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property Dust_Intervals_HardMask_SoftMask_AreConsistent()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(0, 400)
                  from w in Gen.Choose(8, 64)
                  from level in Gen.Choose(2, 64)
                  from linker in Gen.Choose(1, 32)
                  from engine in Gen.Elements(DustEngine.Sdust, DustEngine.Dustmasker)
                  from withN in Gen.Elements(true, false)
                  select (seed, n, w, level, linker, engine, withN);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var rng = new Random(t.seed);
            string s = Segmented(rng, t.n, t.withN ? "NNNRY" : "");
            double threshold = t.level / 10.0;
            var iv = SequenceComplexity.FindLowComplexityIntervals(s, t.w, threshold, t.linker, t.engine);
            var cov = Covered(iv);
            string hard = SequenceComplexity.MaskLowComplexity(s, t.w, threshold, 'X', t.linker, false, t.engine);
            string soft = SequenceComplexity.MaskLowComplexity(s, t.w, threshold, 'X', t.linker, true, t.engine);
            bool ok = WellFormed(iv, s.Length) && hard.Length == s.Length && soft.Length == s.Length
                && Enumerable.Range(0, s.Length).All(i => (hard[i] == 'X') == cov.Contains(i) && (cov.Contains(i) || hard[i] == s[i]))
                && Enumerable.Range(0, s.Length).All(i => char.IsLower(soft[i]) == cov.Contains(i))
                && soft.ToUpperInvariant() == s;
            if (ok && !t.withN && s.Length > 0)
                ok = SequenceComplexity.FindLowComplexityIntervals(new DnaSequence(s), t.w, threshold, t.linker, t.engine).SequenceEqual(iv)
                     && SequenceComplexity.MaskLowComplexity(new DnaSequence(s), t.w, threshold, 'X', t.linker, true, t.engine) == soft;
            return ok.Label($"{t.engine} W={t.w} T={threshold} linker={t.linker} s={s}");
        });
    }

    /// <summary>
    /// F35/F53: dustmasker's linker only joins masked intervals (gap &lt; linker), so the masked set never shrinks as the
    /// linker grows; linker 1 equals the 4-argument (pre-F35) overload.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property Dust_Linker_IsMonotone_AndLinkerOneIsLegacy()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 400)
                  from w in Gen.Choose(8, 64)
                  from level in Gen.Choose(2, 64)
                  from l1 in Gen.Choose(1, 31)
                  from dl in Gen.Choose(1, 32)
                  from engine in Gen.Elements(DustEngine.Sdust, DustEngine.Dustmasker)
                  select (seed, n, w, level, l1, l2: Math.Min(32, l1 + dl), engine);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n);
            double threshold = t.level / 10.0;
            var a = Covered(SequenceComplexity.FindLowComplexityIntervals(s, t.w, threshold, t.l1, t.engine));
            var b = Covered(SequenceComplexity.FindLowComplexityIntervals(s, t.w, threshold, t.l2, t.engine));
            bool legacy = SequenceComplexity.MaskLowComplexity(new DnaSequence(s), t.w, threshold, 'N')
                          == SequenceComplexity.MaskLowComplexity(new DnaSequence(s), t.w, threshold, 'N', 1);
            return (a.IsSubsetOf(b) && legacy).Label($"{t.engine} W={t.w} T={threshold} {t.l1}->{t.l2} s={s}");
        });
    }

    /// <summary>
    /// F36: the sdust engine scans every maximal ACGT run independently (sdust's documented contract): the string result
    /// with N / IUPAC symbols equals the per-run results shifted by the run offsets, after which only the dustmasker
    /// linker rule (merge when the gap is &lt; linker; a no-op for linker 1) can join intervals of neighbouring runs.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property Sdust_NonAcgtInput_EqualsPerAcgtRunResults_JoinedOnlyByLinker()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(0, 400)
                  from w in Gen.Choose(3, 64)
                  from level in Gen.Choose(0, 64)
                  from linker in Gen.Elements(1, 1, 2, 5, 32)
                  select (seed, n, w, level, linker);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n, "NNRYKM");
            double threshold = t.level / 10.0;
            var actual = SequenceComplexity.FindLowComplexityIntervals(s, t.w, threshold, t.linker);
            var expected = new List<(int Start, int End)>();
            int i = 0;
            while (i < s.Length)
            {
                if (!IsAcgt(s[i])) { i++; continue; }
                int j = i;
                while (j < s.Length && IsAcgt(s[j])) j++;
                foreach (var (a, b) in SequenceComplexity.FindLowComplexityIntervals(s[i..j], t.w, threshold, t.linker))
                {
                    if (expected.Count > 0 && i + a - expected[^1].End < t.linker)
                        expected[^1] = (expected[^1].Start, Math.Max(expected[^1].End, i + b));
                    else
                        expected.Add((i + a, i + b));
                }
                i = j;
            }
            return actual.SequenceEqual(expected).Label($"W={t.w} T={threshold} linker={t.linker} s={s}");
        });
    }

    /// <summary>
    /// F53: on ACGT input the dustmasker engine (symdust) and sdust use the same score and the same masking rule; they
    /// differ only by symdust's "single triplet value" shortcut, which masks a homopolymer window without a score
    /// test. A W-base homopolymer window scores (W−2)/2, so for level = 10·T &lt; 5(W−2) the shortcut never changes the
    /// outcome and both engines agree interval for interval.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property Dustmasker_EqualsSdust_OnAcgt_WhenLevelBelowHomopolymerScore()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 500)
                  from w in Gen.Choose(8, 64)
                  from level in Gen.Choose(2, Math.Min(64, 5 * (w - 2) - 1))
                  from linker in Gen.Choose(1, 32)
                  select (seed, n, w, level, linker);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n);
            double threshold = t.level / 10.0;
            var sd = SequenceComplexity.FindLowComplexityIntervals(s, t.w, threshold, t.linker, DustEngine.Sdust);
            var dm = SequenceComplexity.FindLowComplexityIntervals(s, t.w, threshold, t.linker, DustEngine.Dustmasker);
            return sd.SequenceEqual(dm).Label($"W={t.w} level={t.level} linker={t.linker} s={s}");
        });
    }

    /// <summary>
    /// F53: dustmasker (GetDustMasks_SkipNs) reports every N run longer than the window and every leading / trailing N
    /// run as masked.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property Dustmasker_LongAndTerminalNRuns_AreMasked()
    {
        var gen = from seed in Seeds
                  from w in Gen.Choose(8, 40)
                  from level in Gen.Choose(2, 64)
                  select (seed, w, level);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var rng = new Random(t.seed);
            var sb = new System.Text.StringBuilder();
            int parts = rng.Next(1, 6);
            for (int p = 0; p < parts; p++)
            {
                if (rng.Next(2) == 0) sb.Append('N', rng.Next(1, 2 * t.w + 5));
                sb.Append(Segmented(rng, rng.Next(1, 120), "RY"));
            }
            if (rng.Next(2) == 0) sb.Append('N', rng.Next(1, 2 * t.w));
            string s = sb.ToString();
            var cov = Covered(SequenceComplexity.FindLowComplexityIntervals(s, t.w, t.level / 10.0, 1, DustEngine.Dustmasker));
            bool ok = true;
            int i = 0;
            while (i < s.Length)
            {
                if (s[i] != 'N') { i++; continue; }
                int j = i;
                while (j < s.Length && s[j] == 'N') j++;
                if (j - i > t.w || i == 0 || j == s.Length)
                    ok &= Enumerable.Range(i, j - i).All(cov.Contains);
                i = j;
            }
            return ok.Label($"W={t.w} level={t.level} s={s}");
        });
    }

    #endregion

    #region Longdust (F34)

    /// <summary>
    /// F34 (longdust 1.4-r97): the default result is the union of the forward scan and the reverse-complement scan, so
    /// the forward-only result is contained in it; all intervals are sorted, disjoint and in range.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property Longdust_ForwardOnly_IsSubsetOfBothStrands()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(20, 700)
                  from k in Gen.Choose(2, 7)
                  from w in Gen.Choose(20, 400)
                  from withN in Gen.Elements(true, false)
                  select (seed, n, k, w, withN);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n, t.withN ? "N" : "");
            var both = SequenceComplexity.FindLongdustRegions(s, t.k, t.w);
            var fwd = SequenceComplexity.FindLongdustRegions(s, t.k, t.w, forwardOnly: true);
            return (WellFormed(both, s.Length) && WellFormed(fwd, s.Length) && Covered(fwd).IsSubsetOf(Covered(both)))
                .Label($"k={t.k} w={t.w} s={s}");
        });
    }

    /// <summary>
    /// F34: S_L(x) = Σ log c(t)! − f(ℓ/4^k) with uniform base composition depends only on the multiset of k-mer counts,
    /// so a bijective relabelling of A/C/G/T and the reverse complement leave it unchanged; it is 0 when ℓ ≤ 0.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property LongdustScore_InvariantUnderRelabellingAndReverseComplement()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(0, 300)
                  from k in Gen.Choose(1, 8)
                  select (seed, n, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n);
            double v = SequenceComplexity.CalculateLongdustScore(s, t.k);
            string relabelled = new(s.Select(c => c switch { 'A' => 'G', 'C' => 'T', 'G' => 'C', _ => 'A' }).ToArray());
            double r = SequenceComplexity.CalculateLongdustScore(relabelled, t.k);
            double rc = SequenceComplexity.CalculateLongdustScore(ReverseComplement(s), t.k);
            bool ok = Math.Abs(v - r) <= 1e-9 * Math.Max(1, Math.Abs(v)) && Math.Abs(v - rc) <= 1e-9 * Math.Max(1, Math.Abs(v))
                      && (s.Length - t.k + 1 > 0 || v == 0);
            return ok.Label($"k={t.k} S={v} relabel={r} rc={rc} s={s}");
        });
    }

    #endregion

    #region Linguistic complexity, windowed LC, LCR and BBDuk (F37–F39, F55)

    /// <summary>
    /// F37: LC = Σ V_i / Σ min(a^i, N−i+1) with a fixed alphabet a lies in [0, 1] and is non-increasing in a (only the
    /// denominator depends on a); a = 4 equals the default overload on ACGT input; the DnaSequence overload agrees.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property LinguisticComplexity_FixedAlphabet_BoundedAndNonIncreasingInA()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 120)
                  from m in Gen.Choose(1, 14)
                  from alpha in Gen.Elements("A", "AC", "ACG", "ACGT")
                  select (seed, n, m, alpha);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = RandomOver(new Random(t.seed), t.alpha, t.n);
            int distinct = s.Distinct().Count();
            double prev = double.PositiveInfinity;
            bool ok = true;
            for (int a = distinct; a <= distinct + 5; a++)
            {
                double lc = SequenceComplexity.CalculateLinguisticComplexity(s, t.m, a);
                ok &= lc >= 0 && lc <= 1 && lc <= prev + 1e-15
                      && lc == SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(s), t.m, a);
                prev = lc;
            }
            ok &= SequenceComplexity.CalculateLinguisticComplexity(s, t.m, 4) == SequenceComplexity.CalculateLinguisticComplexity(s, t.m)
                  && SequenceComplexity.CalculateLinguisticComplexity(s, t.m) == SequenceComplexity.CalculateLinguisticComplexity(new DnaSequence(s), t.m);
            return ok.Label($"m={t.m} s={s}");
        });
    }

    /// <summary>
    /// F38/F55: every windowed point equals the scalar LC (word cap min(m, w)) and Shannon entropy of its window,
    /// bit for bit, for every m (hash path m ≤ 3, suffix-tree path m ≥ 4); windows holding a non-ACGTU symbol are skipped
    /// and all others are reported; a cap m ≥ w gives the all-length LC.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property WindowedComplexity_EqualsScalarPerWindow_ForEveryWordCap()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(0, 250)
                  from w in Gen.Choose(1, 40)
                  from step in Gen.Choose(1, 12)
                  from m in Gen.Choose(1, 10)
                  select (seed, n, w, step, m);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n, "N");
            var points = SequenceComplexity.CalculateWindowedComplexity(s, t.w, t.step, t.m).ToList();
            var expectedStarts = new List<int>();
            for (int i = 0; i + t.w <= s.Length; i += t.step)
                if (s.Substring(i, t.w).All(IsAcgt)) expectedStarts.Add(i);
            bool ok = points.Select(p => p.WindowStart).SequenceEqual(expectedStarts);
            foreach (var p in points)
            {
                string window = s.Substring(p.WindowStart, t.w);
                ok &= p.WindowEnd == p.WindowStart + t.w - 1 && p.Position == p.WindowStart + t.w / 2
                      && p.LinguisticComplexity == SequenceComplexity.CalculateLinguisticComplexity(window, Math.Min(t.m, t.w))
                      && p.ShannonEntropy == SequenceComplexity.CalculateShannonEntropy(window);
            }
            var capped = SequenceComplexity.CalculateWindowedComplexity(s, t.w, t.step, t.w).ToList();
            var over = SequenceComplexity.CalculateWindowedComplexity(s, t.w, t.step, t.w + 7).ToList();
            ok &= capped.SequenceEqual(over);
            return ok.Label($"w={t.w} step={t.step} m={t.m} s={s}");
        });
    }

    /// <summary>
    /// F39: the per-base Shannon LCR scan with window w and threshold t equals BBDuk's entropy masking with k = 1 and
    /// cutoff t / log2 w (BBDuk normalises by ln(window k-mers); t ≤ log2 w so the cutoff is in [0, 1]); no region contains a non-ACGTU symbol.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property Lcr_EqualsBbdukK1_AndNeverSpansUndefinedSymbols()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(0, 300)
                  from w in Gen.Choose(2, 64)
                  from tenths in Gen.Choose(1, Math.Min(20, (int)Math.Floor(10 * Math.Log2(w))))
                  select (seed, n, w, t: tenths / 10.0);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n, "NRY");
            var lcr = SequenceComplexity.FindLowComplexityRegions(s, t.w, t.t).Select(r => (r.Start, r.End)).ToList();
            var bb = SequenceComplexity.FindLowEntropyRegionsBbduk(s, t.t / Math.Log2(t.w), t.w, 1).Select(r => (r.Start, r.End)).ToList();
            bool clean = SequenceComplexity.FindLowComplexityRegions(s, t.w, t.t).All(r => r.Sequence.All(c => IsAcgt(c) || c == 'U'));
            return (lcr.SequenceEqual(bb) && clean).Label($"w={t.w} t={t.t} s={s}");
        });
    }

    /// <summary>
    /// F39: BBDuk masks a window when its normalised entropy is below the cutoff, so a higher cutoff masks a superset;
    /// the result is invariant under case and U ↔ T (BBTools baseToNumber).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property Bbduk_MonotoneInCutoff_InvariantUnderCaseAndU()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(0, 300)
                  from w in Gen.Choose(5, 60)
                  from k in Gen.Choose(1, 5)
                  from c1 in Gen.Choose(0, 100)
                  from dc in Gen.Choose(0, 100)
                  select (seed, n, w, k, c1: c1 / 100.0, c2: Math.Min(1.0, (c1 + dc) / 100.0));
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n, "N");
            HashSet<int> Mask(string x, double c) =>
                Covered(SequenceComplexity.FindLowEntropyRegionsBbduk(x, c, t.w, t.k).Select(r => (r.Start, r.End + 1)));
            var a = Mask(s, t.c1);
            var b = Mask(s, t.c2);
            string variant = s.Replace('T', 'U').ToLowerInvariant();
            return (a.IsSubsetOf(b) && Mask(variant, t.c1).SetEquals(a)).Label($"w={t.w} k={t.k} {t.c1}<={t.c2} s={s}");
        });
    }

    #endregion

    #region k-mer entropy corrections, normalisation, Digamma (F54)

    /// <summary>
    /// F54: Miller–Madow adds exactly (D − 1)/(2N) nats to the plug-in estimate (D observed k-mers, N = L − k + 1);
    /// None equals the legacy overload; normalisation divides by log2 N (0 for N ≤ 1); DnaSequence and string agree.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property KmerEntropy_MillerMadowOffset_Normalisation_AndLegacyEquality()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 200)
                  from k in Gen.Choose(1, 6)
                  from alpha in Gen.Elements("A", "AC", "ACGT", "ACGT")
                  select (seed, n, k, alpha);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Segmented(new Random(t.seed), t.n);
            if (t.alpha.Length < 4) s = RandomOver(new Random(t.seed), t.alpha, t.n);
            double plug = SequenceComplexity.CalculateKmerEntropy(s, t.k, KmerEntropyCorrection.None);
            double mm = SequenceComplexity.CalculateKmerEntropy(s, t.k, KmerEntropyCorrection.MillerMadow);
            bool ok = plug == SequenceComplexity.CalculateKmerEntropy(s, t.k);
            int windows = s.Length - t.k + 1;
            if (windows >= 1)
            {
                int d = Enumerable.Range(0, windows).Select(i => s.Substring(i, t.k)).Distinct().Count();
                ok &= Math.Abs(mm - plug - (d - 1) / (2.0 * windows * Math.Log(2))) <= 1e-12;
            }
            else
            {
                ok &= plug == 0 && mm == 0;
            }
            foreach (var c in new[] { KmerEntropyCorrection.None, KmerEntropyCorrection.MillerMadow, KmerEntropyCorrection.Grassberger })
            {
                double v = SequenceComplexity.CalculateKmerEntropy(s, t.k, c);
                double norm = SequenceComplexity.CalculateKmerEntropy(s, t.k, c, normalize: true);
                double expected = windows > 1 ? v / Math.Log2(windows) : 0;
                ok &= Math.Abs(norm - expected) <= 1e-12
                      && v == SequenceComplexity.CalculateKmerEntropy(new DnaSequence(s), t.k, c)
                      && norm == SequenceComplexity.CalculateKmerEntropy(new DnaSequence(s), t.k, c, true);
            }
            return ok.Label($"k={t.k} s={s}");
        });
    }

    /// <summary>
    /// F54: ψ(x + 1) = ψ(x) + 1/x (Abramowitz &amp; Stegun 6.3.5); Grassberger's G (2003, eq. 35) satisfies
    /// G(2n+1) = G(2n), G(2n+2) = G(2n) + 2/(2n+1) and G(1) = −γ − ln 2.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Digamma_Recurrence_AndGrassbergerGRecurrence()
    {
        var gen = from num in Gen.Choose(1, 2_000_000)
                  from n in Gen.Choose(1, 5000)
                  select (x: num / 1000.0, n);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            double lhs = StatisticsHelper.Digamma(t.x + 1);
            double rhs = StatisticsHelper.Digamma(t.x) + 1 / t.x;
            bool digamma = Math.Abs(lhs - rhs) <= 1e-10 * Math.Max(1, Math.Abs(lhs));
            double g2n = SequenceComplexity.GrassbergerG(2 * t.n);
            bool odd = Math.Abs(SequenceComplexity.GrassbergerG(2 * t.n + 1) - g2n) <= 1e-10;
            bool even = Math.Abs(SequenceComplexity.GrassbergerG(2 * t.n + 2) - g2n - 2.0 / (2 * t.n + 1)) <= 1e-10;
            bool g1 = Math.Abs(SequenceComplexity.GrassbergerG(1) - (-0.57721566490153286 - Math.Log(2))) <= 1e-12;
            return (digamma && odd && even && g1).Label($"x={t.x} n={t.n}");
        });
    }

    #endregion

    #region TRF parameter API, statistics, mask, flanks, apparent-size tables (F40–F43, F56)

    /// <summary>Tandem-repeat-rich input: planted (possibly mutated) tandem arrays in random flanks, optionally with N.</summary>
    internal static string TandemRich(Random rng, int length, bool withN)
    {
        var sb = new System.Text.StringBuilder(length + 64);
        while (sb.Length < length)
        {
            sb.Append(RandomOver(rng, "ACGT", rng.Next(0, 30)));
            string unit = RandomOver(rng, "ACGT", rng.Next(1, rng.Next(2) == 0 ? 8 : 40));
            int copies = rng.Next(2, 12);
            for (int c = 0; c < copies; c++)
            {
                foreach (char ch in unit)
                {
                    int r = rng.Next(30);
                    if (r == 0) continue;
                    sb.Append(r == 1 ? "ACGT"[rng.Next(4)] : r == 2 && withN ? 'N' : ch);
                    if (r == 3) sb.Append("ACGT"[rng.Next(4)]);
                }
            }
        }
        return sb.ToString(0, length);
    }

    internal static Gen<TandemRepeatsFinderParameters> TrfParameterSets =>
        from match in Gen.Choose(1, 3)
        from mismatch in Gen.Choose(3, 7)
        from delta in Gen.Choose(3, 7)
        from pm in Gen.Elements(75, 80)
        from pi in Gen.Elements(10, 20)
        from minScore in Gen.Choose(20, 80)
        from maxPeriod in Gen.Elements(3, 10, 50, 100, 500)
        from r in Gen.Elements(true, true, false)
        from flank in Gen.Elements(0, 0, 5, 50)
        select new TandemRepeatsFinderParameters
        {
            MatchWeight = match, MismatchPenalty = mismatch, IndelPenalty = delta, MatchProbability = pm,
            IndelProbability = pi, MinScore = minScore, MaxPeriod = maxPeriod, EliminateRedundancy = r, FlankLength = flank,
        };

    internal static string CheckTrfResult(string s, TandemRepeatsFinderParameters p, ApproximateTandemRepeatResult r, int outputCount)
    {
        string upper = s.ToUpperInvariant();
        int n = s.Length;
        if (r.Start < 0 || r.SpanLength < 1 || r.Start + r.SpanLength > n) return "span";
        if (r.Period < 1 || r.Period > p.MaxPeriod) return "period";
        if (r.ConsensusSize < 1 || r.Consensus.Length != r.ConsensusSize) return "consensus";
        if (r.AlignmentScore < p.MinScore) return "score";
        string region = upper.Substring(r.Start, r.SpanLength);
        if (r.AlignedSequence is null || r.AlignedConsensus is null || r.AlignedSequence.Length != r.AlignedConsensus.Length)
            return "rows";
        if (r.AlignedSequence.Replace("-", "") != region) return "aligned sequence";
        int trials = r.CopyMatches + r.CopyMismatches + r.CopyIndels;
        if (trials > 0 && (r.PercentMatches != 100.0 * r.CopyMatches / trials || r.PercentIndels != 100.0 * r.CopyIndels / trials))
            return "percentages";
        int a = region.Count(c => c == 'A'), c = region.Count(ch => ch == 'C'), g = region.Count(ch => ch == 'G'), t = region.Count(ch => ch == 'T');
        int span = r.SpanLength;
        if (r.PercentA != 100.0 * a / span || r.PercentC != 100.0 * c / span || r.PercentG != 100.0 * g / span || r.PercentT != 100.0 * t / span)
            return "base percentages";
        // EntropyTrf uses p_b = count_b / span (TRF get_statistics), Entropy p_b = count_b / ACGT count:
        // EntropyTrf = f·H − f·log2 f with f = ACGT fraction; equal when the region is pure ACGT.
        double f = (double)(a + c + g + t) / span;
        double expectedTrf = f == 0 ? 0 : f * r.Entropy - f * Math.Log2(f);
        if (Math.Abs(r.EntropyTrf - expectedTrf) > 1e-9) return $"EntropyTrf {r.EntropyTrf} vs {expectedTrf}";
        if (r.OutputIndex < 1 || r.OutputIndex > r.OutputCount || r.OutputCount != outputCount) return "output index";
        if (p.FlankLength > 0)
        {
            int ls = Math.Max(0, r.Start - p.FlankLength), re = Math.Min(n, r.Start + r.SpanLength + p.FlankLength);
            if (r.LeftFlank != upper[ls..r.Start] || r.RightFlank != upper[(r.Start + r.SpanLength)..re]) return "flanks";
        }
        else if (r.LeftFlank is not null || r.RightFlank is not null)
        {
            return "flanks without -f";
        }
        return "";
    }

    /// <summary>
    /// F40–F43/F57/F64: for random TRF parameter sets (PM 75/80, PI 10/20, -r, -f) every reported repeat lies in the
    /// sequence, respects MaxPeriod and Minscore, its ungapped alignment row is the region, its percentages equal the
    /// adjacent-copy counts, its base percentages and EntropyTrf follow TRF get_statistics, its flanks are the -f
    /// windows and 1 ≤ OutputIndex ≤ OutputCount = the out count; output indices are distinct.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 60)]
    public Property Trf_RandomParameterSets_ResultInvariants()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 400)
                  from withN in Gen.Elements(true, false)
                  from lower in Gen.Elements(true, false, false)
                  from p in TrfParameterSets
                  select (seed, n, withN, lower, p);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = TandemRich(new Random(t.seed), t.n, t.withN);
            if (t.lower) s = s.ToLowerInvariant();
            var rows = RepeatFinder.FindApproximateTandemRepeats(s, t.p, out int outputCount);
            string errors = string.Join(",", rows.Select(r => CheckTrfResult(s, t.p, r, outputCount)).Where(e => e.Length > 0));
            bool distinct = rows.Select(r => r.OutputIndex).Distinct().Count() == rows.Count;
            bool same = rows.SequenceEqual(RepeatFinder.FindApproximateTandemRepeats(s, t.p));
            return (errors.Length == 0 && distinct && same).Label($"{errors} p={t.p} s={s}");
        });
    }

    private static (int, int, int, string, int) Key(ApproximateTandemRepeatResult r) =>
        (r.Start, r.SpanLength, r.Period, r.Consensus, r.AlignmentScore);

    /// <summary>
    /// F41: <c>-r</c> only disables redundancy elimination, so its rows contain the default rows; a maximum TR length
    /// (<c>-l</c>) of at least n never truncates and equals the default; supplying the exact apparent-size table equals
    /// null (F56).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 50)]
    public Property Trf_NoRedundancyElimination_Superset_MaxLengthAndExactTableNeutral()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 350)
                  from p in TrfParameterSets
                  select (seed, n, p: p with { EliminateRedundancy = true });
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = TandemRich(new Random(t.seed), t.n, withN: false);
            var def = RepeatFinder.FindApproximateTandemRepeats(s, t.p).ToList();
            var all = RepeatFinder.FindApproximateTandemRepeats(s, t.p with { EliminateRedundancy = false }).Select(Key).ToHashSet();
            bool superset = def.Select(Key).All(all.Contains);
            bool maxLength = RepeatFinder.FindApproximateTandemRepeats(s, t.p with { MaxRepeatLength = s.Length + 1 }).SequenceEqual(def);
            var exact = TandemRepeatsFinderParameters.ExactApparentSizeTable(t.p.MatchProbability);
            bool table = RepeatFinder.FindApproximateTandemRepeats(s, t.p with { ApparentSizeTable = exact }).SequenceEqual(def);
            return (superset && maxLength && table).Label($"superset={superset} maxLength={maxLength} table={table} p={t.p} s={s}");
        });
    }

    /// <summary>
    /// F42: TRF -m masks exactly the union of the reported spans (hard: N, soft: lower case), keeps every other symbol,
    /// and the parameter overload equals masking its own results.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 50)]
    public Property TrfMask_IsUnionOfReportedSpans()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(0, 350)
                  from p in TrfParameterSets
                  select (seed, n, p);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = TandemRich(new Random(t.seed), t.n, withN: false);
            var rows = RepeatFinder.FindApproximateTandemRepeats(s, t.p).ToList();
            var cov = Covered(rows.Select(r => (r.Start, r.Start + r.SpanLength)));
            string hard = RepeatFinder.MaskApproximateTandemRepeats(s, t.p);
            string soft = RepeatFinder.MaskApproximateTandemRepeats(s, t.p, softMask: true);
            bool ok = hard.Length == s.Length && soft.ToUpperInvariant() == s
                && Enumerable.Range(0, s.Length).All(i => (hard[i] == 'N') == cov.Contains(i) && char.IsLower(soft[i]) == cov.Contains(i))
                && hard == RepeatFinder.MaskApproximateTandemRepeats(s, rows);
            return ok.Label($"p={t.p} s={s}");
        });
    }

    /// <summary>
    /// F56: the exact tables lie in TRF's admissible range 0..max(d,20)−1; a table built from waiting times w is
    /// max(d,20) − w − 1 entry by entry; any valid random table is accepted and keeps the result invariants.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 30)]
    public Property ApparentSizeTables_RangeConversionAndRandomTablesAccepted()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 300)
                  from p in TrfParameterSets
                  select (seed, n, p);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var rng = new Random(t.seed);
            int len = TandemRepeatsFinderParameters.ApparentSizeTableLength;
            bool ok = true;
            foreach (int pm in new[] { 75, 80 })
            {
                var exact = TandemRepeatsFinderParameters.ExactApparentSizeTable(pm);
                ok &= exact.Length == len && Enumerable.Range(1, len - 1).All(d => exact[d] >= 0 && exact[d] <= Math.Max(d, 20) - 1);
            }
            var waits = new int[len];
            for (int d = 1; d < len; d++) waits[d] = rng.Next(0, Math.Max(d, 20));
            var fromWaits = TandemRepeatsFinderParameters.ApparentSizeTableFromWaitingTimes(waits);
            ok &= Enumerable.Range(1, len - 1).All(d => fromWaits[d] == Math.Max(d, 20) - waits[d] - 1);

            string s = TandemRich(rng, t.n, withN: rng.Next(2) == 0);
            var p = t.p with { ApparentSizeTable = fromWaits };
            var rows = RepeatFinder.FindApproximateTandemRepeats(s, p, out int count);
            ok &= rows.All(r => CheckTrfResult(s, p, r, count).Length == 0);
            return ok.Label($"p={t.p} s={s}");
        });
    }

    #endregion

    #region TRF output formatters (F57, F64, F65)

    private static List<ApproximateTandemRepeatResult> Replicated(IReadOnlyList<ApproximateTandemRepeatResult> rows, int times)
    {
        var list = new List<ApproximateTandemRepeatResult>();
        for (int k = 0; k < times; k++) list.AddRange(rows);
        int total = list.Count;
        return list.Select((r, i) => r with { OutputIndex = i + 1, OutputCount = total }).ToList();
    }

    private static int CountOf(string text, string part)
    {
        int count = 0;
        for (int i = text.IndexOf(part, StringComparison.Ordinal); i >= 0; i = text.IndexOf(part, i + part.Length, StringComparison.Ordinal))
            count++;
        return count;
    }

    /// <summary>
    /// F57/F64 (TestSpec INV-13–INV-15, input without D/F letters): .dat rows follow OutputIndex order with 15 fields
    /// (-ngs: 17); HTML table and alignment pages number max(1, ⌈rows/120⌉); each repeat gets exactly one "Found at"
    /// on its alignment page, in OutputIndex order, and every table anchor resolves to a NAME on the matching page.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 25)]
    public Property TrfFormatters_RowOrder_FieldCounts_PageCounts_AndAnchors()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(40, 400)
                  from times in Gen.Elements(1, 1, 3, 40, 130)
                  from p in TrfParameterSets
                  select (seed, n, times, p);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = TandemRich(new Random(t.seed), t.n, withN: true);
            var found = RepeatFinder.FindApproximateTandemRepeats(s, t.p, out int outputCount);
            var rows = t.times == 1 ? found.ToList() : Replicated(found, t.times);
            int count = t.times == 1 ? outputCount : rows.Count;
            var ordered = rows.OrderBy(r => r.OutputIndex).ToList();

            var datLines = RepeatFinder.FormatTrfDatLines(s, rows, "seq", t.p).Split('\n')
                .Where(l => l.Length > 0 && char.IsDigit(l[0])).ToList();
            var ngsLines = RepeatFinder.FormatTrfDatLines(s, rows, "seq", t.p, TrfDatLayout.Ngs).Split('\n')
                .Where(l => l.Length > 0 && char.IsDigit(l[0])).ToList();
            bool dat = datLines.Count == rows.Count && ngsLines.Count == rows.Count
                && datLines.All(l => l.Split(' ').Length == 15) && ngsLines.All(l => l.Split(' ').Length == 17)
                && datLines.Select(l => l.Split(' ')[0] + "-" + l.Split(' ')[1])
                    .SequenceEqual(ordered.Select(r => $"{r.Start + 1}-{r.Start + r.SpanLength}"));

            int expectedPages = Math.Max(1, (rows.Count + 119) / 120);
            var tables = RepeatFinder.FormatTrfHtmlTables(s, rows, "seq", t.p, "f");
            var pages = RepeatFinder.FormatTrfAlignmentPages(s, rows, "seq", t.p, "f", count);
            bool pageCounts = tables.Count == expectedPages && pages.Count == expectedPages;

            bool foundAt = pages.Sum(pg => CountOf(pg.Html, "\nFound at i:")) == rows.Count;
            var names = pages.SelectMany(pg => System.Text.RegularExpressions.Regex.Matches(pg.Html, "<A NAME=\"([^\"]+)\"></A>")
                .Select(m => m.Groups[1].Value)).ToList();
            bool order = names.Select(x => int.Parse(x[(x.LastIndexOf(',') + 1)..])).SequenceEqual(ordered.Select(r => r.OutputIndex));

            bool anchors = true;
            for (int k = 0; k < tables.Count; k++)
            {
                string pageName = pages[k].FileName;
                foreach (System.Text.RegularExpressions.Match m in System.Text.RegularExpressions.Regex.Matches(tables[k].Html, "<A HREF=\"([^\"#]+)#([^\"]+)\">"))
                    anchors &= m.Groups[1].Value == pageName && pages[k].Html.Contains($"<A NAME=\"{m.Groups[2].Value}\"></A>", StringComparison.Ordinal);
            }
            return (dat && pageCounts && foundAt && order && anchors)
                .Label($"dat={dat} pages={pageCounts} foundAt={foundAt} order={order} anchors={anchors} rows={rows.Count} p={t.p} s={s}");
        });
    }

    /// <summary>
    /// F57: FormatCFixed reproduces C printf("%.Nf"): the exact binary value of the double, rounded to N places with
    /// ties to even. Oracle: exact rational arithmetic on the IEEE-754 decomposition.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 500)]
    public Property FormatCFixed_EqualsExactRationalRounding()
    {
        var gen = from mant in Gen.Choose(0, int.MaxValue - 1)
                  from scale in Gen.Choose(0, 9)
                  from decimals in Gen.Choose(0, 4)
                  from half in Gen.Elements(false, true)
                  select (value: half ? (mant % 100000) / 8.0 + 0.125 * (mant % 2) : mant / Math.Pow(10, scale), decimals);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string actual = RepeatFinder.FormatCFixed(t.value, t.decimals);
            string expected = ExactFixed(t.value, t.decimals);
            return (actual == expected).Label($"{t.value:R} %.{t.decimals}f: {actual} vs {expected}");
        });
    }

    private static string ExactFixed(double value, int decimals)
    {
        long bits = BitConverter.DoubleToInt64Bits(value);
        int exp = (int)((bits >> 52) & 0x7FF);
        long man = bits & 0xF_FFFF_FFFF_FFFFL;
        if (exp == 0) exp = 1; else man |= 1L << 52;
        exp -= 1075;
        BigInteger num = man * BigInteger.Pow(10, decimals), den = BigInteger.One;
        if (exp >= 0) num <<= exp; else den <<= -exp;
        var q = BigInteger.DivRem(num, den, out var rem);
        int cmp = (rem * 2).CompareTo(den);
        if (cmp > 0 || (cmp == 0 && !q.IsEven)) q += 1;
        string digits = q.ToString().PadLeft(decimals + 1, '0');
        return decimals == 0 ? digits : digits[..^decimals] + "." + digits[^decimals..];
    }

    #endregion

    #region Degenerate repeats and Vmatch reporting (F47, F58, F59)

    private static Gen<(string s, int k, int minLength, bool edit, bool pal)> SmallDegenerateCases =>
        from seed in Seeds
        from n in Gen.Choose(6, 22)
        from k in Gen.Choose(1, 2)
        from extra in Gen.Choose(1, 5)
        from edit in Gen.Elements(true, false)
        from pal in Gen.Elements(true, false)
        from alpha in Gen.Elements("ACGT", "AC", "ACGTN", "AT")
        select (RandomOver(new Random(seed), alpha, n), k, k + extra, edit, pal);

    /// <summary>
    /// F47: FindDegenerateRepeats equals an independent brute force of Vmatch's App. A definitions (all four modes:
    /// edit / Hamming × direct / palindromic, with wildcards) on small random inputs.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property DegenerateRepeats_EqualBruteForce_AllModes()
    {
        return Prop.ForAll(SmallDegenerateCases.ToArbitrary(), t =>
        {
            var expected = RepeatFinder_DegenerateRepeats_Tests.BruteForce(t.s, t.k, t.minLength, t.edit, t.pal).OrderBy(x => x).ToList();
            var actual = RepeatFinder.FindDegenerateRepeats(t.s, t.minLength, t.k,
                    t.edit ? ApproximateRepeatDistance.Edit : ApproximateRepeatDistance.Hamming, t.pal, int.MaxValue, int.MinValue)
                .Select(x => (x.FirstPosition, x.SecondPosition, x.FirstLength, x.SecondLength, x.Distance)).OrderBy(x => x).ToList();
            return actual.SequenceEqual(expected).Label($"{t.s} k={t.k} l={t.minLength} edit={t.edit} pal={t.pal}");
        });
    }

    private static Gen<(string s, int k, int minLength, ApproximateRepeatDistance d, bool pal)> MediumDegenerateCases =>
        from seed in Seeds
        from n in Gen.Choose(20, 120)
        from k in Gen.Choose(1, 3)
        from extra in Gen.Choose(2, 10)
        from d in Gen.Elements(ApproximateRepeatDistance.Edit, ApproximateRepeatDistance.Hamming)
        from pal in Gen.Elements(true, false)
        select (Segmented(new Random(seed), n, "N"), k, k + extra, d, pal);

    private static int TrueDistance(DegenerateRepeatResult x, ApproximateRepeatDistance d)
    {
        string second = x.IsReverseComplement ? ReverseComplement(x.SecondCopy) : x.SecondCopy;
        return d == ApproximateRepeatDistance.Edit ? WildcardEditDistance(x.FirstCopy, second) : WildcardHamming(x.FirstCopy, second);
    }

    /// <summary>
    /// F47/F58/F59: in every mode and reporting option, copies are the stated substrings, both instances are at least
    /// minLength long, k ≥ Distance ≥ the true (wildcard) distance of the copies; the definition-correct AllMaximal mode
    /// reports exactly the true distance; Hamming rows have equal lengths; the result does not depend on input case.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 80)]
    public Property DegenerateRepeats_RowsAreValidMatches_InEveryReportingMode()
    {
        return Prop.ForAll(MediumDegenerateCases.ToArbitrary(), t =>
        {
            string upper = t.s;
            var problems = new List<string>();
            foreach (var reporting in new[] { DegenerateRepeatReporting.AllMaximal, DegenerateRepeatReporting.BestPerSeed })
            {
                foreach (bool compat in new[] { false, true })
                {
                    var rows = RepeatFinder.FindDegenerateRepeats(upper, t.minLength, t.k, t.d, t.pal, int.MaxValue, int.MinValue, reporting, compat).ToList();
                    foreach (var x in rows)
                    {
                        int truth = TrueDistance(x, t.d);
                        bool ok = x.FirstCopy == upper.Substring(x.FirstPosition, x.FirstLength)
                            && x.SecondCopy == upper.Substring(x.SecondPosition, x.SecondLength)
                            && x.FirstLength >= t.minLength && x.SecondLength >= t.minLength
                            && x.Distance <= t.k && truth <= x.Distance && x.IsReverseComplement == t.pal
                            && x.Spacing == x.SecondPosition - x.FirstPosition - x.FirstLength
                            && (t.d == ApproximateRepeatDistance.Edit || x.FirstLength == x.SecondLength)
                            && (reporting != DegenerateRepeatReporting.AllMaximal || compat || truth == x.Distance);
                        if (!ok) problems.Add($"{reporting}/{compat}: {x} truth={truth}");
                    }
                    var lowerRows = RepeatFinder.FindDegenerateRepeats(upper.ToLowerInvariant(), t.minLength, t.k, t.d, t.pal,
                        int.MaxValue, int.MinValue, reporting, compat).ToList();
                    if (!lowerRows.SequenceEqual(rows)) problems.Add($"{reporting}/{compat}: case");
                }
            }
            return (problems.Count == 0).Label($"{string.Join("; ", problems.Take(3))} s={t.s} k={t.k} l={t.minLength} {t.d} pal={t.pal}");
        });
    }

    /// <summary>
    /// F59: the stock-Vmatch seed shortcut does not change the Hamming maximal set (extendHD.c), so compat AllMaximal
    /// equals the default for Hamming distance.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 80)]
    public Property DegenerateRepeats_HammingCompatAllMaximal_EqualsDefault()
    {
        return Prop.ForAll(MediumDegenerateCases.ToArbitrary(), t =>
        {
            var def = RepeatFinder.FindDegenerateRepeats(t.s, t.minLength, t.k, ApproximateRepeatDistance.Hamming, t.pal, int.MaxValue, int.MinValue,
                DegenerateRepeatReporting.AllMaximal).ToList();
            var compat = RepeatFinder.FindDegenerateRepeats(t.s, t.minLength, t.k, ApproximateRepeatDistance.Hamming, t.pal, int.MaxValue, int.MinValue,
                DegenerateRepeatReporting.AllMaximal, vmatchCompatible: true).ToList();
            return def.SequenceEqual(compat).Label($"s={t.s} k={t.k} l={t.minLength} pal={t.pal}");
        });
    }

    /// <summary>
    /// F47: a k-mismatch match is also a k-differences match, so every non-overlapping direct Hamming repeat is contained
    /// (both instances) in some maximal direct edit repeat.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 80)]
    public Property DegenerateRepeats_HammingRepeat_LiesInsideAnEditRepeat()
    {
        return Prop.ForAll(MediumDegenerateCases.ToArbitrary(), t =>
        {
            var edit = RepeatFinder.FindDegenerateRepeats(t.s, t.minLength, t.k, ApproximateRepeatDistance.Edit, false, int.MaxValue, int.MinValue).ToList();
            var missing = RepeatFinder.FindDegenerateRepeats(t.s, t.minLength, t.k, ApproximateRepeatDistance.Hamming, false, int.MaxValue, int.MinValue)
                .Where(h => h.FirstPosition + h.FirstLength <= h.SecondPosition)
                .Where(h => !edit.Any(e => e.FirstPosition <= h.FirstPosition && e.FirstPosition + e.FirstLength >= h.FirstPosition + h.FirstLength
                                        && e.SecondPosition <= h.SecondPosition && e.SecondPosition + e.SecondLength >= h.SecondPosition + h.SecondLength))
                .ToList();
            return (missing.Count == 0).Label($"{string.Join(";", missing.Take(3))} s={t.s} k={t.k} l={t.minLength}");
        });
    }

    #endregion

    #region MISA scan mode, unit sizes, summaries (F48, F51, F61)

    private static Gen<Dictionary<int, int>> MisaMaps =>
        from count in Gen.Choose(1, 6)
        from sizes in Gen.Choose(1, 10).ArrayOf(count)
        from mins in Gen.Choose(2, 6).ArrayOf(count)
        select sizes.Zip(mins).GroupBy(p => p.First).ToDictionary(g => g.Key, g => g.First().Second);

    /// <summary>SSR-rich input with unit sizes 1–10, optional N / IUPAC / lower case.</summary>
    internal static string SsrRich(Random rng, int length, bool dirty)
    {
        var sb = new System.Text.StringBuilder();
        while (sb.Length < length)
        {
            sb.Append(RandomOver(rng, dirty ? "ACGTACGTNR" : "ACGT", rng.Next(0, 12)));
            string unit = RandomOver(rng, "ACGT", rng.Next(1, 11));
            sb.Append(string.Concat(Enumerable.Repeat(unit, rng.Next(2, 9))));
        }
        string s = sb.ToString(0, length);
        return dirty && rng.Next(2) == 0 ? s.ToLowerInvariant() : s;
    }

    /// <summary>
    /// F48/F61: the MisaRegex scan equals the line-by-line transcription of misa.pl's scan loop (lines 101–125) for
    /// random definitions including unit sizes 7–10.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property MisaRegex_EqualsMisaPlTranscription_ForRandomDefinitions()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(0, 300)
                  from map in MisaMaps
                  from dirty in Gen.Elements(true, false)
                  select (seed, n, map, dirty);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = SsrRich(new Random(t.seed), t.n, t.dirty);
            string actual = RepeatFinder_MisaScan_Tests.Ssrs(RepeatFinder.FindMicrosatellites(s, t.map, MicrosatelliteScanMode.MisaRegex));
            string expected = RepeatFinder_MisaScan_Tests.Ssrs(RepeatFinder_MisaScan_Tests.MisaPlScan(s, t.map));
            return (actual == expected).Label($"map={string.Join(' ', t.map.Select(kv => $"{kv.Key}-{kv.Value}"))} s={s}");
        });
    }

    /// <summary>
    /// F61: CountsByUnitLength has exactly the definition's sizes as keys and sums to TotalRepeats; sizes 1–6 equal the
    /// named class counts (0 when not searched); PercentageOfSequence is in [0, 100]; equality and hash agree between
    /// the string and DnaSequence paths and with a re-ordered copy of the counts.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property TandemSummary_CountsByUnitLength_Consistent_EqualityAndHash()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 300)
                  from map in MisaMaps
                  from mode in Gen.Elements(MicrosatelliteScanMode.MaximalRuns, MicrosatelliteScanMode.MisaRegex)
                  select (seed, n, map, mode);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = SsrRich(new Random(t.seed), t.n, dirty: false);
            var sum = RepeatFinder.GetTandemRepeatSummary(s, t.map, t.mode);
            var counts = sum.CountsByUnitLength;
            int[] named = [sum.MononucleotideRepeats, sum.DinucleotideRepeats, sum.TrinucleotideRepeats,
                sum.TetranucleotideRepeats, sum.PentanucleotideRepeats, sum.HexanucleotideRepeats];
            bool ok = counts.Keys.OrderBy(x => x).SequenceEqual(t.map.Keys.OrderBy(x => x))
                && counts.Values.Sum() == sum.TotalRepeats
                && Enumerable.Range(1, 6).All(u => named[u - 1] == counts.GetValueOrDefault(u))
                && sum.PercentageOfSequence >= 0 && sum.PercentageOfSequence <= 100;
            var dna = RepeatFinder.GetTandemRepeatSummary(new DnaSequence(s), t.map, t.mode);
            var reordered = sum with { CountsByUnitLength = counts.Reverse().ToDictionary(kv => kv.Key, kv => kv.Value) };
            ok &= dna.Equals(sum) && dna.GetHashCode() == sum.GetHashCode()
                  && reordered.Equals(sum) && reordered.GetHashCode() == sum.GetHashCode();
            if (sum.TotalRepeats > 0)
                ok &= !(sum with { CountsByUnitLength = counts.ToDictionary(kv => kv.Key, kv => kv.Value + 1) }).Equals(sum);
            return ok.Label($"{t.mode} s={s}");
        });
    }

    /// <summary>
    /// F51: N / IUPAC never form or extend an SSR, so the string summary of a sequence equals the sum of the summaries of
    /// its maximal ACGT pieces (counts per unit size, total repeats, total bases, covered bases), with the percentage
    /// taken over the full length (misa.pl <c>length $seq</c>).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 120)]
    public Property StringSummary_IsSumOverAcgtPieces()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 300)
                  from map in MisaMaps
                  from mode in Gen.Elements(MicrosatelliteScanMode.MaximalRuns, MicrosatelliteScanMode.MisaRegex)
                  select (seed, n, map, mode);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = SsrRich(new Random(t.seed), t.n, dirty: true);
            var whole = RepeatFinder.GetTandemRepeatSummary(s, t.map, t.mode);
            string upper = s.ToUpperInvariant();
            var pieces = new List<(TandemRepeatSummary Summary, int Length)>();
            int i = 0;
            while (i < upper.Length)
            {
                if (!IsAcgt(upper[i])) { i++; continue; }
                int j = i;
                while (j < upper.Length && IsAcgt(upper[j])) j++;
                pieces.Add((RepeatFinder.GetTandemRepeatSummary(upper[i..j], t.map, t.mode), j - i));
                i = j;
            }
            double covered = pieces.Sum(p => p.Summary.PercentageOfSequence * p.Length / 100.0);
            bool ok = whole.TotalRepeats == pieces.Sum(p => p.Summary.TotalRepeats)
                && whole.TotalRepeatBases == pieces.Sum(p => p.Summary.TotalRepeatBases)
                && t.map.Keys.All(u => whole.CountsByUnitLength[u] == pieces.Sum(p => p.Summary.CountsByUnitLength[u]))
                && Math.Abs(whole.PercentageOfSequence - 100.0 * covered / s.Length) <= 1e-9
                && whole.PercentageOfSequence >= 0 && whole.PercentageOfSequence <= 100;
            return ok.Label($"{t.mode} s={s}");
        });
    }

    /// <summary>
    /// F51: on ACGT input the string summary overloads equal the DnaSequence ones (uniform minimum, map, map + both scan
    /// modes).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property StringSummaryOverloads_EqualDnaSequenceOverloads_OnAcgt()
    {
        var gen = from seed in Seeds
                  from n in Gen.Choose(1, 300)
                  from min in Gen.Choose(2, 6)
                  from map in MisaMaps
                  select (seed, n, min, map);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = SsrRich(new Random(t.seed), t.n, dirty: false);
            var dna = new DnaSequence(s);
            bool ok = RepeatFinder.GetTandemRepeatSummary(s, t.min).Equals(RepeatFinder.GetTandemRepeatSummary(dna, t.min))
                && RepeatFinder.GetTandemRepeatSummary(s, t.map).Equals(RepeatFinder.GetTandemRepeatSummary(dna, t.map))
                && RepeatFinder.GetTandemRepeatSummary(s, t.map, MicrosatelliteScanMode.MaximalRuns)
                    .Equals(RepeatFinder.GetTandemRepeatSummary(dna, t.map, MicrosatelliteScanMode.MaximalRuns))
                && RepeatFinder.GetTandemRepeatSummary(s, t.map, MicrosatelliteScanMode.MisaRegex)
                    .Equals(RepeatFinder.GetTandemRepeatSummary(dna, t.map, MicrosatelliteScanMode.MisaRegex));
            return ok.Label($"s={s}");
        });
    }

    /// <summary>F61: formatting a definition as misa.ini <c>size-min</c> pairs (optionally with the <c>def</c> token, commas) parses back to it.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property ParseMisaDefinition_RoundTrip()
    {
        var gen = from map in MisaMaps
                  from prefix in Gen.Elements("", "definition(unit_size,min_repeats): ", "def ")
                  from sep in Gen.Elements(" ", ", ", "\t", ",")
                  select (map, prefix, sep);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string text = t.prefix + string.Join(t.sep, t.map.Select(kv => $"{kv.Key}-{kv.Value}"));
            var parsed = RepeatFinder.ParseMisaDefinition(text);
            bool ok = parsed.Count == t.map.Count && t.map.All(kv => parsed[kv.Key] == kv.Value)
                      && parsed.Keys.SequenceEqual(parsed.Keys.OrderBy(x => x));
            return ok.Label(text);
        });
    }

    #endregion
}
