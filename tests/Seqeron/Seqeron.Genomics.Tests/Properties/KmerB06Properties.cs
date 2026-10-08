using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// B06 heavy tier — property-based tests for the behaviour added to <see cref="KmerAnalyzer"/> by batch B06
/// (docs/Validation/review-2026-09/B06.md, F1–F38): option-aware counting (Jellyfish ACGT-only window rule and
/// canonical <c>-C</c>, F10), parallel counting (F14), Jellyfish <c>histo</c> (F12), clump windows (F12), the
/// <c>dump -L/-U</c> range filter and its order contract (F9), exact Jaccard / containment / Mash distance (F11, F24),
/// Mash MinHash sketches (F23), sourmash FracMinHash sketches, downsampling and abundances (F30, F32, F33), the KMP
/// position scan (F5) and the O(k) generator (F4).
///
/// Every property states an identity that follows from the definition cited in the method's XML documentation
/// (the oracle is an independent re-computation, never the observed output).
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class KmerB06Properties
{
    #region Generators / helpers

    private static readonly Gen<int> Seeds = Gen.Choose(0, int.MaxValue - 1);

    private static readonly Gen<KmerCountingOptions> AnyOptions =
        Gen.Elements(
            KmerCountingOptions.Default,
            new KmerCountingOptions(AcgtOnly: true),
            new KmerCountingOptions(Canonical: true),
            new KmerCountingOptions(Canonical: true, AcgtOnly: true));

    /// <summary>
    /// Repeat-rich sequence: random segments over <paramref name="alphabet"/> interleaved with tandem copies of short
    /// units (so k-mer counts above 1, clumps and shared k-mers actually occur).
    /// </summary>
    internal static string Mixed(Random rng, int length, string alphabet)
    {
        var sb = new System.Text.StringBuilder(length + 64);
        while (sb.Length < length)
        {
            if (rng.Next(2) == 0)
            {
                int len = rng.Next(1, 30);
                for (int i = 0; i < len; i++) sb.Append(alphabet[rng.Next(alphabet.Length)]);
            }
            else
            {
                string unit = new(Enumerable.Range(0, rng.Next(1, 6)).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());
                int copies = rng.Next(2, 10);
                for (int c = 0; c < copies; c++) sb.Append(unit);
            }
        }
        return sb.ToString(0, length);
    }

    /// <summary>A point-mutated copy of <paramref name="s"/> (substitution rate 1/<paramref name="every"/>).</summary>
    internal static string Mutate(Random rng, string s, int every, string alphabet)
    {
        var chars = s.ToCharArray();
        for (int i = 0; i < chars.Length; i++)
            if (rng.Next(every) == 0)
                chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    private static bool SameTable(IReadOnlyDictionary<string, int> a, IReadOnlyDictionary<string, int> b) =>
        a.Count == b.Count && a.All(kv => b.TryGetValue(kv.Key, out var v) && v == kv.Value);

    private static bool IsAcgtUpper(char c) => char.ToUpperInvariant(c) is 'A' or 'C' or 'G' or 'T';

    #endregion

    #region Counting options and parallel counting (F10, F13, F14)

    /// <summary>
    /// F14: <c>CountKmersParallel</c> is documented to return exactly the serial <c>CountKmers</c> table for every
    /// input, option and degree (counting is a sum over a partition of the windows). Inputs above
    /// 2·<see cref="KmerAnalyzer.ParallelMinWindowsPerRange"/> windows exercise the multi-range path.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 6)]
    public Property CountKmersParallel_EqualsSerial_MultiRange()
    {
        var gen = from seed in Seeds
                  from extra in Gen.Choose(20, 60000)
                  from k in Gen.Choose(1, 14)
                  from options in AnyOptions
                  from degree in Gen.Elements(-1, 2, 3, 4)
                  select (seed, len: 2 * KmerAnalyzer.ParallelMinWindowsPerRange + extra, k, options, degree);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTACGTACGTNacgtR");
            var serial = KmerAnalyzer.CountKmers(s, t.k, t.options);
            var parallel = KmerAnalyzer.CountKmersParallel(s, t.k, t.options, t.degree);
            return SameTable(serial, parallel).Label($"len={t.len} k={t.k} {t.options} degree={t.degree}");
        });
    }

    /// <summary>F14: below the range threshold the parallel overload falls back to the serial count — same table.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property CountKmersParallel_EqualsSerial_SmallInputs()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 400)
                  from k in Gen.Choose(1, 9)
                  from options in AnyOptions
                  from degree in Gen.Elements(-1, 1, 2, 8)
                  select (seed, len, k, options, degree);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTNacgtRY-");
            return SameTable(KmerAnalyzer.CountKmers(s, t.k, t.options),
                    KmerAnalyzer.CountKmersParallel(s, t.k, t.options, t.degree))
                .Label($"len={t.len} k={t.k} {t.options}");
        });
    }

    /// <summary>
    /// F10: canonical counting keys each k-mer by min(w, RC(w)) over the ACGT windows (Jellyfish <c>-C</c>), so the
    /// canonical table of RC(s) equals that of s — every ACGT window of s is an ACGT window of RC(s) with reverse-complemented
    /// content, and the IUPAC complement keeps non-ACGT symbols non-ACGT. Holds with lower case, N, IUPAC and gaps.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property CanonicalCounts_InvariantUnderReverseComplement()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 300)
                  from k in Gen.Choose(1, 8)
                  select (seed, len, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTACGTacgtNRYKM-");
            string rc = DnaSequence.GetReverseComplementString(s) ?? string.Empty;
            var canonical = new KmerCountingOptions(Canonical: true);
            return SameTable(KmerAnalyzer.CountKmers(s, t.k, canonical), KmerAnalyzer.CountKmers(rc, t.k, canonical))
                .Label($"k={t.k} s={s}");
        });
    }

    /// <summary>
    /// F10: the ACGT-only rule only removes windows that contain a non-ACGT symbol, so on pure (mixed-case) ACGT input
    /// the ACGT-only table equals the literal table.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property AcgtOnlyCounts_EqualLiteral_OnPureAcgt()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 300)
                  from k in Gen.Choose(1, 10)
                  select (seed, len, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTacgt");
            return SameTable(KmerAnalyzer.CountKmers(s, t.k), KmerAnalyzer.CountKmers(s, t.k, new KmerCountingOptions(AcgtOnly: true)))
                .Label($"k={t.k} s={s}");
        });
    }

    /// <summary>
    /// F10: Σ canonical counts = Σ ACGT-only counts = the number of windows whose k symbols are all A/C/G/T after case
    /// folding (independent scan), and every canonical key is upper case, ACGT-only and ≤ its reverse complement (ordinal).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property CanonicalTotal_EqualsAcgtWindowCount_KeysAreCanonical()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 300)
                  from k in Gen.Choose(1, 8)
                  select (seed, len, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTACGTacgtNRW");
            int acgtWindows = Enumerable.Range(0, Math.Max(0, s.Length - t.k + 1))
                .Count(i => s.Substring(i, t.k).All(IsAcgtUpper));
            var canonical = KmerAnalyzer.CountKmers(s, t.k, new KmerCountingOptions(Canonical: true));
            var acgtOnly = KmerAnalyzer.CountKmers(s, t.k, new KmerCountingOptions(AcgtOnly: true));
            bool keysOk = canonical.Keys.All(w => w.All(c => c is 'A' or 'C' or 'G' or 'T')
                && string.CompareOrdinal(w, DnaSequence.GetReverseComplementString(w)) <= 0);
            return (canonical.Values.Sum() == acgtWindows && acgtOnly.Values.Sum() == acgtWindows && keysOk)
                .Label($"k={t.k} windows={acgtWindows} canon={canonical.Values.Sum()} acgt={acgtOnly.Values.Sum()}");
        });
    }

    /// <summary>F10: <c>DistinctKmers(s, k, options)</c> is documented as the key set of <c>CountKmers(s, k, options)</c>.</summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property DistinctKmers_EqualsCountKmersKeys()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 300)
                  from k in Gen.Choose(1, 8)
                  from options in AnyOptions
                  select (seed, len, k, options);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTacgtN");
            return KmerAnalyzer.DistinctKmers(s, t.k, t.options).SetEquals(KmerAnalyzer.CountKmers(s, t.k, t.options).Keys)
                .Label($"k={t.k} {t.options}");
        });
    }

    #endregion

    #region Jellyfish histo, clump windows, dump -L/-U (F9, F12)

    /// <summary>
    /// F12: each distinct k-mer is tallied in exactly one bucket (Jellyfish <c>histo_main.cc</c>), so the frequencies sum to
    /// the number of distinct k-mers (with and without <c>--full</c>); with low = 1 and high ≥ the maximum count every
    /// bucket label is the multiplicity itself, so Σ bin·frequency = Σ counts, and the default call matches the spectrum.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property Histogram_BinsSumToDistinct_AndWeightedSumToTotal()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 400)
                  from k in Gen.Choose(1, 6)
                  from options in AnyOptions
                  select (seed, len, k, options);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTacgtN");
            var counts = KmerAnalyzer.CountKmers(s, t.k, t.options);
            int distinct = counts.Count;
            long total = counts.Values.Sum();
            var hist = KmerAnalyzer.GetKmerHistogram(s, t.k, t.options, low: 1, high: 1_000_000);
            long maxCount = counts.Count == 0 ? 1 : counts.Values.Max();
            var full = KmerAnalyzer.GetKmerHistogram(s, t.k, t.options, low: 1, high: maxCount, full: true);
            var spectrum = KmerAnalyzer.GetKmerSpectrum(s, t.k, t.options);
            var dflt = KmerAnalyzer.GetKmerHistogram(s, t.k, t.options);
            bool spectrumOk = dflt.Select(b => (b.Bin, b.Frequency))
                .SequenceEqual(spectrum.OrderBy(kv => kv.Key).Select(kv => ((long)kv.Key, (long)kv.Value)));
            return (hist.Sum(b => b.Frequency) == distinct
                    && hist.Sum(b => b.Bin * b.Frequency) == total
                    && full.Sum(b => b.Frequency) == distinct
                    && hist.Select(b => b.Bin).SequenceEqual(hist.Select(b => b.Bin).Order())
                    && spectrumOk)
                .Label($"k={t.k} {t.options} distinct={distinct} total={total}");
        });
    }

    /// <summary>
    /// F12: <c>FindClumpWindows</c> reports the same k-mer set as <c>FindClumps</c>, and its runs are exactly the maximal runs
    /// of window starts i ∈ [0, n − L] in which the k-mer has ≥ t occurrences at starts p ∈ [i, i + L − k] (brute force over
    /// every window, Compeau &amp; Pevzner ch. 1 / Rosalind BA1E). Order: FirstWindowStart, then ordinal k-mer; runs ascending,
    /// disjoint and non-adjacent.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property ClumpWindows_MatchBruteForce_AndFindClumpsSet()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 120)
                  from k in Gen.Choose(1, 4)
                  from extraL in Gen.Choose(0, 40)
                  from t in Gen.Choose(1, 4)
                  select (seed, len, k, L: k + extraL, t);
        return Prop.ForAll(gen.ToArbitrary(), p =>
        {
            string s = Mixed(new Random(p.seed), p.len, "ACGTacgt");
            string u = s.ToUpperInvariant();
            var windows = KmerAnalyzer.FindClumpWindows(s, p.k, p.L, p.t);
            var clumps = KmerAnalyzer.FindClumps(s, p.k, p.L, p.t).ToHashSet(StringComparer.Ordinal);

            var expected = new SortedSet<(string, int)>();
            for (int i = 0; i + p.L <= u.Length; i++)
            {
                var c = new Dictionary<string, int>(StringComparer.Ordinal);
                for (int q = i; q <= i + p.L - p.k; q++)
                {
                    string w = u.Substring(q, p.k);
                    c[w] = c.GetValueOrDefault(w) + 1;
                }
                foreach (var (w, n) in c)
                    if (n >= p.t)
                        expected.Add((w, i));
            }
            var actual = new SortedSet<(string, int)>();
            bool runsOk = true;
            foreach (var clump in windows)
            {
                for (int r = 0; r < clump.WindowRuns.Count; r++)
                {
                    var run = clump.WindowRuns[r];
                    runsOk &= run.FirstWindowStart <= run.LastWindowStart;
                    if (r > 0) runsOk &= run.FirstWindowStart > clump.WindowRuns[r - 1].LastWindowStart + 1;
                    for (int i = run.FirstWindowStart; i <= run.LastWindowStart; i++)
                        actual.Add((clump.Kmer, i));
                }
            }
            bool orderOk = windows.Zip(windows.Skip(1)).All(z =>
                z.First.FirstWindowStart < z.Second.FirstWindowStart
                || (z.First.FirstWindowStart == z.Second.FirstWindowStart && string.CompareOrdinal(z.First.Kmer, z.Second.Kmer) < 0));
            return (actual.SetEquals(expected) && runsOk && orderOk
                    && clumps.SetEquals(windows.Select(w => w.Kmer)))
                .Label($"k={p.k} L={p.L} t={p.t} s={s}");
        });
    }

    /// <summary>
    /// F9: <c>FindKmersWithMinCount(min, max, options)</c> = Jellyfish <c>dump -L min -U max</c>: exactly the k-mers of the
    /// option-aware table whose count lies in [min, max], with their counts, ordered by count descending then ordinal k-mer.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property MinCountRange_IsExactSubset_InDocumentedOrder()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 300)
                  from k in Gen.Choose(1, 5)
                  from min in Gen.Choose(-2, 6)
                  from max in Gen.Choose(0, 8)
                  from options in AnyOptions
                  select (seed, len, k, min, max, options);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTacgtN");
            var counts = KmerAnalyzer.CountKmers(s, t.k, t.options);
            var result = KmerAnalyzer.FindKmersWithMinCount(s, t.k, t.min, t.max, t.options).ToList();
            int expectedCount = counts.Count(kv => kv.Value >= t.min && kv.Value <= t.max);
            bool exact = result.All(r => r.Count >= t.min && r.Count <= t.max && counts[r.Kmer] == r.Count);
            bool ordered = result.Zip(result.Skip(1)).All(z =>
                z.First.Count > z.Second.Count
                || (z.First.Count == z.Second.Count && string.CompareOrdinal(z.First.Kmer, z.Second.Kmer) < 0));
            var unique = KmerAnalyzer.FindUniqueKmers(s, t.k, t.options).ToList();
            bool uniqueOk = unique.SequenceEqual(counts.Where(kv => kv.Value == 1).Select(kv => kv.Key).Order(StringComparer.Ordinal));
            return (result.Count == expectedCount && exact && ordered && uniqueOk)
                .Label($"k={t.k} [{t.min},{t.max}] {t.options}");
        });
    }

    #endregion

    #region Exact Jaccard / containment / Mash distance (F11, F24)

    /// <summary>
    /// F11/F24: J(A, B) = |A∩B| / |A∪B| is symmetric, in [0, 1], 1 for identical non-empty sets (0 by convention for
    /// two empty sets); C(A, A) = 1 for a non-empty set; C(A, B)·|A| = C(B, A)·|B| = |A∩B|; the Mash distance is
    /// symmetric, in [0, 1], and 0 for identical sets.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 300)]
    public Property Jaccard_Containment_Mash_SetIdentities()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 200)
                  from k in Gen.Choose(1, 8)
                  from options in AnyOptions
                  select (seed, len, k, options);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var rng = new Random(t.seed);
            string a = Mixed(rng, t.len, "ACGTacgtN");
            string b = Mutate(rng, a, 6, "ACGTN");
            var setA = KmerAnalyzer.DistinctKmers(a, t.k, t.options);
            var setB = KmerAnalyzer.DistinctKmers(b, t.k, t.options);
            double j = KmerAnalyzer.JaccardSimilarity(a, b, t.k, t.options);
            double jSwap = KmerAnalyzer.JaccardSimilarity(b, a, t.k, t.options);
            double jSelf = KmerAnalyzer.JaccardSimilarity(a, a, t.k, t.options);
            double cSelf = KmerAnalyzer.ContainmentIndex(a, a, t.k, t.options);
            double cAB = KmerAnalyzer.ContainmentIndex(a, b, t.k, t.options);
            double cBA = KmerAnalyzer.ContainmentIndex(b, a, t.k, t.options);
            double m = KmerAnalyzer.MashDistance(a, b, t.k, t.options);
            double mSwap = KmerAnalyzer.MashDistance(b, a, t.k, t.options);
            double mSelf = KmerAnalyzer.MashDistance(a, a, t.k, t.options);
            int shared = setA.Count(setB.Contains);
            int union = setA.Count + setB.Count - shared;
            double jExpected = union == 0 ? 0.0 : (double)shared / union;
            return (j == jSwap && j >= 0 && j <= 1 && Math.Abs(j - jExpected) < 1e-15
                    && jSelf == (setA.Count == 0 ? 0.0 : 1.0)
                    && cSelf == (setA.Count == 0 ? 0.0 : 1.0)
                    && Math.Abs(cAB * setA.Count - shared) < 1e-9 && Math.Abs(cBA * setB.Count - shared) < 1e-9
                    && m == mSwap && m >= 0 && m <= 1 && mSelf == 0.0)
                .Label($"k={t.k} {t.options} J={j} C={cAB}/{cBA} D={m}");
        });
    }

    #endregion

    #region MinHash and FracMinHash sketches (F23, F30, F32, F33)

    /// <summary>
    /// F23: a bottom-s sketch with s ≥ |A ∪ B| keeps every hash of both sets, so <c>mash dist</c>'s x/s is the exact
    /// Jaccard index of the canonical (or, with canonical = false, ACGT-only forward) k-mer sets and the distance equals
    /// <see cref="KmerAnalyzer.MashDistance"/> (k ≥ 17: 64-bit hashes, collisions negligible).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property MinHash_SketchCoveringUnion_EqualsExactJaccard()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 300)
                  from k in Gen.Choose(17, 24)
                  from canonical in Gen.Elements(true, false)
                  from slack in Gen.Choose(0, 5)
                  select (seed, len, k, canonical, slack);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var rng = new Random(t.seed);
            string a = Mixed(rng, t.len, "ACGTACGTacgtN");
            string b = Mutate(rng, a, 25, "ACGT");
            var options = t.canonical ? new KmerCountingOptions(Canonical: true) : new KmerCountingOptions(AcgtOnly: true);
            var setA = KmerAnalyzer.DistinctKmers(a, t.k, options);
            var setB = KmerAnalyzer.DistinctKmers(b, t.k, options);
            int union = setA.Union(setB).Count();
            int s = Math.Max(1, union + t.slack);
            var cmp = KmerAnalyzer.CompareMinHashSketches(
                KmerAnalyzer.CreateMinHashSketch(a, t.k, s, t.canonical),
                KmerAnalyzer.CreateMinHashSketch(b, t.k, s, t.canonical));
            double exact = KmerAnalyzer.JaccardSimilarity(a, b, t.k, options);
            double mash = KmerAnalyzer.MashDistance(a, b, t.k, options);
            bool distanceOk = union == 0 ? cmp.Distance == 0.0 : Math.Abs(cmp.Distance - mash) < 1e-12;
            return (cmp.Denominator == union && Math.Abs(cmp.Jaccard - exact) < 1e-12 && distanceOk
                    && cmp.PValue >= 0 && cmp.PValue <= 1)
                .Label($"k={t.k} canonical={t.canonical} union={union} x/s={cmp.SharedHashes}/{cmp.Denominator} exact={exact}");
        });
    }

    /// <summary>
    /// F30: with scaled = 1 the FracMinHash sketch keeps every hash (max_hash = 2^64 − 1) and the bias factor is 1, so
    /// sourmash's jaccard and contained_by equal the exact canonical Jaccard and containment indices.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 200)]
    public Property FracMinHash_Scaled1_EqualsExactJaccardAndContainment()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 300)
                  from k in Gen.Choose(1, 31)
                  select (seed, len, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var rng = new Random(t.seed);
            string a = Mixed(rng, t.len, "ACGTACGTacgtN");
            string b = Mutate(rng, a.Substring(0, rng.Next(a.Length + 1)), 20, "ACGT");
            var canonical = new KmerCountingOptions(Canonical: true);
            var cmp = KmerAnalyzer.CompareFracMinHashSketches(
                KmerAnalyzer.CreateFracMinHashSketch(a, t.k, 1),
                KmerAnalyzer.CreateFracMinHashSketch(b, t.k, 1));
            double j = KmerAnalyzer.JaccardSimilarity(a, b, t.k, canonical);
            double cAB = KmerAnalyzer.ContainmentIndex(a, b, t.k, canonical);
            double cBA = KmerAnalyzer.ContainmentIndex(b, a, t.k, canonical);
            return (Math.Abs(cmp.Jaccard - j) < 1e-12 && Math.Abs(cmp.ContainmentAInB - cAB) < 1e-12
                    && Math.Abs(cmp.ContainmentBInA - cBA) < 1e-12
                    && Math.Abs(cmp.MaxContainment - Math.Max(cAB, cBA)) < 1e-12)
                .Label($"k={t.k} J={cmp.Jaccard}/{j} C={cmp.ContainmentAInB}/{cAB}");
        });
    }

    /// <summary>
    /// F32/F33: downsampling a sketch from S to S′ ≥ S keeps exactly the hashes ≤ max_hash(S′) with their abundances —
    /// the sketch that sketching the sequence at S′ gives (sourmash <c>MinHash.downsample</c>).
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 150)]
    public Property FracMinHash_Downsample_EqualsDirectSketch()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 3000)
                  from k in Gen.Choose(1, 31)
                  from scaled in Gen.Choose(1, 40)
                  from factor in Gen.Choose(1, 25)
                  from track in Gen.Elements(true, false)
                  select (seed, len, k, scaled: (long)scaled, newScaled: (long)scaled * factor, track);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            string s = Mixed(new Random(t.seed), t.len, "ACGTACGTacgtN");
            var down = KmerAnalyzer.DownsampleFracMinHash(
                KmerAnalyzer.CreateFracMinHashSketch(s, t.k, t.scaled, trackAbundance: t.track), t.newScaled);
            var direct = KmerAnalyzer.CreateFracMinHashSketch(s, t.k, t.newScaled, trackAbundance: t.track);
            bool abundOk = t.track
                ? down.Abundances!.SequenceEqual(direct.Abundances!)
                : down.Abundances is null && direct.Abundances is null;
            return (down.Scaled == direct.Scaled && down.MaxHash == direct.MaxHash
                    && down.Hashes.SequenceEqual(direct.Hashes) && abundOk)
                .Label($"k={t.k} {t.scaled}->{t.newScaled} track={t.track}");
        });
    }

    #endregion

    #region KMP positions and the k-mer generator (F4, F5)

    /// <summary>
    /// F5: the KMP scan reports exactly the overlapping occurrences found by a naive window-by-window comparison, in
    /// ascending order, case-insensitively — stressed on a two-letter alphabet with self-overlapping (periodic) patterns
    /// where the prefix-function fallback matters, and on absent and longer-than-text patterns.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 400)]
    public Property KmpPositions_EqualNaiveScan()
    {
        var gen = from seed in Seeds
                  from len in Gen.Choose(0, 200)
                  from plen in Gen.Choose(1, 12)
                  select (seed, len, plen);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var rng = new Random(t.seed);
            string s = Mixed(rng, t.len, "AAAAAACaac");
            string pattern = rng.Next(3) == 0 && s.Length >= t.plen
                ? s.Substring(rng.Next(s.Length - t.plen + 1), t.plen).ToLowerInvariant()
                : Mixed(rng, t.plen, "AAACa");
            string su = s.ToUpperInvariant(), pu = pattern.ToUpperInvariant();
            var naive = Enumerable.Range(0, Math.Max(0, su.Length - pu.Length + 1))
                .Where(i => string.CompareOrdinal(su, i, pu, 0, pu.Length) == 0).ToList();
            return KmerAnalyzer.FindKmerPositions(s, pattern).SequenceEqual(naive)
                .Label($"s={s} pattern={pattern}");
        });
    }

    /// <summary>
    /// F4: <c>GenerateAllKmers</c> is the k-fold Cartesian product in alphabet order (itertools.product): |Σ|^k strings,
    /// equal to an independent recursive enumeration, for arbitrary distinct alphabets.
    /// </summary>
    [FsCheck.NUnit.Property(MaxTest = 100)]
    public Property GenerateAllKmers_EqualsRecursiveProduct()
    {
        var gen = from seed in Seeds
                  from size in Gen.Choose(1, 5)
                  from k in Gen.Choose(1, 5)
                  select (seed, size, k);
        return Prop.ForAll(gen.ToArbitrary(), t =>
        {
            var rng = new Random(t.seed);
            string alphabet = new("ACGTNXYZacgt".OrderBy(_ => rng.Next()).Take(t.size).ToArray());
            IEnumerable<string> Product(int depth) =>
                depth == 0 ? [string.Empty] : Product(depth - 1).SelectMany(p => alphabet.Select(c => p + c));
            var expected = Product(t.k).ToList();
            var actual = KmerAnalyzer.GenerateAllKmers(t.k, alphabet).ToList();
            return (actual.Count == (int)Math.Pow(t.size, t.k) && actual.SequenceEqual(expected))
                .Label($"alphabet={alphabet} k={t.k}");
        });
    }

    #endregion
}
