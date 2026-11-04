using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Review;

/// <summary>
/// Independent review of 1o.2 (docs/refactoring/suffix-tree-upgrade/adr/REVIEW-S1a-opt.md): the software prefetch of the induction passes and of the naming has no
/// architectural effect, so a wrong distance, side or address passes every output test (mutants P04-P07 of 1o.2 and R20, R24, R40, R43 of the review survived all of
/// them). These tests PIN the prefetch targets: a recording <see cref="ITypedText"/> logs every <c>Prefetch</c> and every access of a pass, and the position prefetched in
/// iteration <c>i</c> must be exactly the one that iteration <c>i + SaisPrefetch.Induce</c> (<c>SaisPrefetch.Name</c> for the naming) reads. The pass is run on the FINAL suffix
/// array of a text (a state in which the look-ahead entries are final, and which the induction pass reproduces unchanged), so the match is exact, not statistical.
/// The prefetches of the LCP kernels go through the concrete <see cref="Buf{T}"/> and cannot be intercepted this way (their mutants Q06-Q08, D23, R69 stay "equivalent by nature").
/// </summary>
[TestFixture]
public class ReviewS1aOptPrefetchTests
{
    private sealed class Log
    {
        internal readonly List<int> Prefetches = [];

        internal readonly List<int> Accesses = [];

        /// <summary>For the naming: the iteration (number of Prefetch calls so far minus 1) in which each access happened.</summary>
        internal readonly List<int> AccessIteration = [];
    }

    private readonly struct Recording(PackedBytes inner, Log log) : ITypedText
    {
        public void Get(int i, out int symbol, out int isS)
        {
            log.Accesses.Add(i);
            log.AccessIteration.Add(log.Prefetches.Count - 1);
            inner.Get(i, out symbol, out isS);
        }

        public int Symbol(int i)
        {
            log.Accesses.Add(i);
            log.AccessIteration.Add(log.Prefetches.Count - 1);
            return inner.Symbol(i);
        }

        public int IsS(int i)
        {
            log.Accesses.Add(i);
            log.AccessIteration.Add(log.Prefetches.Count - 1);
            return inner.IsS(i);
        }

        public bool Same(int x, int y)
        {
            log.Accesses.Add(x);
            log.AccessIteration.Add(log.Prefetches.Count - 1);
            log.Accesses.Add(y);
            log.AccessIteration.Add(log.Prefetches.Count - 1);
            return inner.Same(x, y);
        }

        public void Prefetch(int i) => log.Prefetches.Add(i);
    }

    /// <summary>Everything a pass needs, built from the definition: final suffix array (n + 1 layout), packed symbol+type bytes, bucket heads and exclusive tails.</summary>
    private sealed class Level
    {
        internal required int N { get; init; }

        internal required int[] Sa { get; init; }

        internal required byte[] Packed { get; init; }

        internal required int[] Heads { get; init; }

        internal required int[] Tails { get; init; }

        internal required int M { get; init; }
    }

    private static Level Build(string text)
    {
        PreparedText p = TextPreparation.Prepare(text);
        int n = p.Length;
        int[] sa = SaisReference.BuildSuffixArray(text);
        var rank = new int[n + 1];
        for (int r = 0; r <= n; r++) rank[sa[r]] = r;
        var isS = new int[n + 1];
        isS[n] = 1;
        for (int i = n - 1; i >= 0; i--) isS[i] = rank[i] < rank[i + 1] ? 1 : 0; // S-type iff the suffix is smaller than its successor
        var packed = new byte[n + 1];
        for (int i = 0; i < n; i++) packed[i] = (byte)((p.Bytes[i] << 1) | isS[i]);
        var counts = new int[p.AlphabetSize];
        for (int i = 0; i < n; i++) counts[p.Bytes[i]]++;
        var heads = new int[counts.Length];
        var tails = new int[counts.Length];
        int sum = 1;
        for (int c = 0; c < counts.Length; c++)
        {
            heads[c] = sum;
            sum += counts[c];
            tails[c] = sum;
        }

        int m = 0;
        for (int i = 1; i < n; i++) if (isS[i] == 1 && isS[i - 1] == 0) m++;
        return new Level { N = n, Sa = sa, Packed = packed, Heads = heads, Tails = tails, M = m };
    }

    private static IEnumerable<string> Texts()
    {
        yield return StringFamilies.RandomDna(3000, 5);
        yield return StringFamilies.Random(2500, 20, 6);
        yield return StringFamilies.Random(2000, 100, 7);
        yield return StringFamilies.Fibonacci(2100);
        yield return StringFamilies.Periodic(2200, 7, 8);
        yield return StringFamilies.ThueMorse(1900);
    }

    [Test]
    public void LPass_PrefetchesThePredecessorOfTheEntryScannedInduceIterationsLater()
    {
        foreach (string text in Texts())
        {
            Level lv = Build(text);
            using var arena = new NativeArena(ArenaOptions.Guarded());
            Buf<int> sa = arena.Allocate<int>(lv.N + 1, zero: false, "sa");
            lv.Sa.CopyTo(sa.AsSpan(lv.N + 1));
            Buf<int> bkt = arena.Allocate<int>(lv.Heads.Length, zero: false, "bkt");
            lv.Heads.CopyTo(bkt.AsSpan(lv.Heads.Length));
            Buf<byte> packed = arena.Allocate<byte>(lv.N + 1, zero: false, "packed");
            lv.Packed.CopyTo(packed.AsSpan(lv.N + 1));
            var log = new Log();
            SaisInduction.InduceL(new Recording(new PackedBytes(packed), log), sa, bkt, lv.N);
            arena.VerifyCanaries();

            Assert.That(sa.AsSpan(lv.N + 1).ToArray(), Is.EqualTo(lv.Sa), "inducing from the final suffix array reproduces it");
            Assert.That(log.Prefetches, Has.Count.EqualTo(lv.N + 1), "one prefetch per scanned entry");
            Assert.That(log.Accesses, Has.Count.EqualTo(lv.N + 1), "one text access per scanned entry");
            int d = SaisPrefetch.Induce;
            Assert.That(d, Is.InRange(8, 128), "the distance must stay on the measured plateau (PROGRESS-S1.md 1o.2: flat from 12 to 64 iterations)");
            for (int i = 0; i + d <= lv.N; i++)
                if (log.Prefetches[i] != log.Accesses[i + d])
                    Assert.Fail($"L pass, iteration {i}: prefetched position {log.Prefetches[i]}, iteration {i + d} reads position {log.Accesses[i + d]}");
            int tail = Math.Max(lv.Sa[lv.N] - 1, 0); // the look-ahead is clamped to SA[n]: its predecessor (0 for an empty or position-0 entry)
            for (int i = Math.Max(0, lv.N - d + 1); i <= lv.N; i++)
                Assert.That(log.Prefetches[i], Is.EqualTo(tail), $"L pass, iteration {i}: the look-ahead is clamped to SA[n]");
        }
    }

    [Test]
    public void SPass_PrefetchesThePredecessorOfTheEntryScannedInduceIterationsLater()
    {
        foreach (string text in Texts())
        {
            Level lv = Build(text);
            using var arena = new NativeArena(ArenaOptions.Guarded());
            Buf<int> sa = arena.Allocate<int>(lv.N + 1, zero: false, "sa");
            lv.Sa.CopyTo(sa.AsSpan(lv.N + 1));
            Buf<int> bkt = arena.Allocate<int>(lv.Tails.Length, zero: false, "bkt");
            lv.Tails.CopyTo(bkt.AsSpan(lv.Tails.Length));
            Buf<byte> packed = arena.Allocate<byte>(lv.N + 1, zero: false, "packed");
            lv.Packed.CopyTo(packed.AsSpan(lv.N + 1));
            var log = new Log();
            SaisInduction.InduceS(new Recording(new PackedBytes(packed), log), sa, bkt, lv.N);
            arena.VerifyCanaries();

            Assert.That(sa.AsSpan(lv.N + 1).ToArray(), Is.EqualTo(lv.Sa), "inducing from the final suffix array reproduces it");
            Assert.That(log.Prefetches, Has.Count.EqualTo(lv.N), "one prefetch per scanned entry (i = n .. 1)");
            int d = SaisPrefetch.Induce;
            Assert.That(d, Is.InRange(8, 128), "the distance must stay on the measured plateau (PROGRESS-S1.md 1o.2: flat from 12 to 64 iterations)");
            // iteration k scans index i = n - k; the entry that iteration k + d scans is SA[n - k - d] (clamped to index 1 by the look-ahead)
            for (int k = 0; (lv.N - k) - d >= 1; k++)
                if (log.Prefetches[k] != log.Accesses[k + d])
                    Assert.Fail($"S pass, iteration {k}: prefetched position {log.Prefetches[k]}, iteration {k + d} reads position {log.Accesses[k + d]}");
            int tail = Math.Max(lv.Sa[1] - 1, 0); // clamped to SA[1]
            for (int k = 0; k < lv.N; k++)
                if (lv.N - k - d < 1)
                    Assert.That(log.Prefetches[k], Is.EqualTo(tail), $"S pass, iteration {k}: the look-ahead is clamped to SA[1]");
        }
    }

    [Test]
    public void Naming_PrefetchesTheTextAtTheLmsPositionComparedNameIterationsLater()
    {
        foreach (string text in Texts())
        {
            Level lv = Build(text);
            int n = lv.N, m = lv.M;
            using var arena = new NativeArena(ArenaOptions.Guarded());
            Buf<int> sa = arena.Allocate<int>(n + 1, zero: false, "sa");
            // state before the naming: SA[1..m] = the LMS positions in suffix order (LMS substrings are sorted in this order, equal ones are adjacent)
            sa[0] = n;
            int k = 0;
            for (int r = 1; r <= n; r++)
            {
                int pos = lv.Sa[r];
                bool lms = pos > 0 && (lv.Packed[pos] & 1) == 1 && (lv.Packed[pos - 1] & 1) == 0;
                if (lms) sa[1 + k++] = pos;
            }

            Assert.That(k, Is.EqualTo(m));
            Buf<byte> packed = arena.Allocate<byte>(n + 1, zero: false, "packed");
            lv.Packed.CopyTo(packed.AsSpan(n + 1));
            var log = new Log();
            int names = SaisInduction.NameLmsSubstrings(new Recording(new PackedBytes(packed), log), sa, n, m);
            arena.VerifyCanaries();

            Assert.That(names, Is.InRange(1, m));
            Assert.That(log.Prefetches, Has.Count.EqualTo(m), "one prefetch per sorted LMS position");
            int d = SaisPrefetch.Name;
            Assert.That(d, Is.InRange(4, 64), "the naming distance must stay on the measured plateau (PROGRESS-S1.md 1o.2: 8-64 iterations)");
            int lastLms = lv.Sa.Skip(1).Last(pos => pos > 0 && (lv.Packed[pos] & 1) == 1 && (lv.Packed[pos - 1] & 1) == 0);
            for (int i = Math.Max(0, m - d); i < m; i++)
                Assert.That(log.Prefetches[i], Is.EqualTo(lastLms), $"naming, iteration {i}: the look-ahead is clamped to SA[m], the last sorted LMS position");
            var accessedIn = new List<HashSet<int>>();
            for (int i = 0; i < m; i++) accessedIn.Add([]);
            for (int a = 0; a < log.Accesses.Count; a++) accessedIn[log.AccessIteration[a]].Add(log.Accesses[a]);
            for (int i = 0; i + d < m; i++)
                if (!accessedIn[i + d].Contains(log.Prefetches[i]))
                    Assert.Fail($"naming, iteration {i}: prefetched position {log.Prefetches[i]}, iteration {i + d} does not read it");
        }
    }
}
