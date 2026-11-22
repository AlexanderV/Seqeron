using SuffixTree.Esa.Bench.Links;

namespace SuffixTree.Esa.Tests.Handoff;

/// <summary>
/// Work package 1b.9f: the pieces of the hand-off query check (<c>sesa handoff queries</c> of the Bench application, compiled into this project from <c>apps/SuffixTree.Esa.Bench/Links</c>):
/// the deterministic query set, the streamed oracle (a rolling hash over one pass of the input text, every candidate verified byte by byte) and the binary search over a suffix array that
/// answers from the container's SA and the input text without the engine. Every piece is compared with the naive definition on random texts.
/// </summary>
[TestFixture]
public sealed class HandoffOracleTests
{
    private static byte[] RandomText(int seed, int n, string alphabet)
    {
        var rng = new Random(seed);
        var t = new byte[n];
        for (int i = 0; i < n; i++) t[i] = (byte)alphabet[rng.Next(alphabet.Length)];
        return t;
    }

    private static List<long> NaivePositions(byte[] text, byte[] q)
    {
        var r = new List<long>();
        for (int i = 0; i + q.Length <= text.Length; i++)
        {
            if (text.AsSpan(i, q.Length).SequenceEqual(q)) r.Add(i);
        }
        return r;
    }

    private static Func<long, int, byte[]> Reader(byte[] text) => (pos, len) => text.AsSpan((int)pos, len).ToArray();

    private static long[] NaiveSa(byte[] text)
    {
        // the suffix array with the empty suffix first: positions 0..n sorted by the suffix, a shorter suffix before a longer one that it prefixes
        long[] sa = [.. Enumerable.Range(0, text.Length + 1).Select(i => (long)i)];
        Array.Sort(sa, (a, b) => text.AsSpan((int)a).SequenceCompareTo(text.AsSpan((int)b)));
        return sa;
    }

    // ---------------------------------------------------------------- the query set

    [Test]
    public void TheAlphabet_IsTheSortedDistinctBytes()
    {
        Assert.That(HandoffQuerySet.AlphabetOf("GATTACA"u8), Is.EqualTo("ACGT"u8.ToArray()));
        Assert.That(HandoffQuerySet.AlphabetOf([]), Is.Empty);
    }

    [Test]
    public void TheSet_IsDeterministic_AndTheSeedChangesIt()
    {
        byte[] text = RandomText(1, 5000, "ACGT");
        var a = HandoffQuerySet.Build(text.Length, Reader(text), "ACGT"u8, [12, 20], 30, 20, 7);
        var b = HandoffQuerySet.Build(text.Length, Reader(text), "ACGT"u8, [12, 20], 30, 20, 7);
        var c = HandoffQuerySet.Build(text.Length, Reader(text), "ACGT"u8, [12, 20], 30, 20, 8);
        Assert.That(a, Has.Count.EqualTo(100));
        for (int i = 0; i < a.Count; i++)
        {
            Assert.That(b[i].Bytes, Is.EqualTo(a[i].Bytes));
            Assert.That(b[i].Origin, Is.EqualTo(a[i].Origin));
        }
        Assert.That(a.Select(q => Convert.ToHexString(q.Bytes)).SequenceEqual(c.Select(q => Convert.ToHexString(q.Bytes))), Is.False);
    }

    [Test]
    public void ThePresentQueries_AreTheTextAtTheirOrigin_AndTheAbsentOnesUseTheAlphabet()
    {
        byte[] text = RandomText(2, 3000, "ACGTN");
        var set = HandoffQuerySet.Build(text.Length, Reader(text), "ACGT"u8, [8, 16, 32], 40, 25, 99);
        Assert.That(set, Has.Count.EqualTo(3 * (40 + 25)));
        Assert.That(set.Select(q => q.Index), Is.EqualTo(Enumerable.Range(0, set.Count)));
        foreach (HandoffQuery q in set)
        {
            Assert.That(q.Bytes.Length, Is.AnyOf(8, 16, 32));
            if (q.Present)
            {
                Assert.That(q.Origin, Is.InRange(0L, text.Length - q.Bytes.Length));
                Assert.That(text.AsSpan((int)q.Origin, q.Bytes.Length).SequenceEqual(q.Bytes), Is.True);
            }
            else
            {
                Assert.That(q.Origin, Is.EqualTo(-1));
                Assert.That(q.Bytes.All(b => "ACGT"u8.Contains(b)), Is.True);
            }
        }
        Assert.That(set.Count(q => q.Present), Is.EqualTo(3 * 40));
    }

    [Test]
    public void TheCutPositions_ReachTheLastWindow_AndNeverLeaveTheText()
    {
        byte[] text = RandomText(3, 40, "AC");
        var set = HandoffQuerySet.Build(text.Length, Reader(text), "AC"u8, [40], 50, 0, 5);
        // (an empty set would pass Is.All: the audit's mutant "a window as long as the text is no present query" survived until the count was asserted)
        Assert.That(set, Has.Count.EqualTo(50), "a window as long as the text is a present query");
        Assert.That(set.Select(q => q.Origin), Is.All.EqualTo(0L), "a window as long as the text can only start at 0");
        set = HandoffQuerySet.Build(text.Length, Reader(text), "AC"u8, [39], 200, 0, 5);
        Assert.That(set.Select(q => q.Origin).Distinct().Order(), Is.EqualTo(new long[] { 0, 1 }));
    }

    [Test]
    public void ALengthAboveTheText_GivesNoPresentQuery_ButTheAbsentOnesStay()
    {
        byte[] text = RandomText(4, 10, "AC");
        var set = HandoffQuerySet.Build(text.Length, Reader(text), "AC"u8, [20], 5, 3, 1);
        Assert.That(set.Count(q => q.Present), Is.Zero);
        Assert.That(set.Count(q => !q.Present), Is.EqualTo(3));
    }

    // ---------------------------------------------------------------- the streamed oracle

    [Test]
    public void TheScan_EqualsTheNaiveCountAndPositions_AtEveryBlockSize()
    {
        foreach ((int seed, string alphabet, int n) in new[] { (11, "ACGT", 4000), (12, "AB", 3000), (13, "A", 500), (14, "ACGTN", 2500) })
        {
            byte[] text = RandomText(seed, n, alphabet);
            var set = HandoffQuerySet.Build(n, Reader(text), [.. alphabet.Select(c => (byte)c)], [1, 2, 5, 12, 31], 40, 20, (ulong)seed);
            foreach (int block in new[] { 1, 7, 33, 100, 4096, 1 << 20 })
            {
                ScanResult r = RollingScan.Run(new MemoryStream(text), set, 50, block);
                foreach (HandoffQuery q in set)
                {
                    List<long> naive = NaivePositions(text, q.Bytes);
                    Assert.That(r.Counts[q.Index], Is.EqualTo(naive.Count), $"count of query {q.Index} (len {q.Bytes.Length}) at block {block}, text {alphabet}");
                    Assert.That(r.Positions[q.Index], Is.EqualTo(naive.Take(50).ToList()), $"positions of query {q.Index} at block {block}");
                }
            }
        }
    }

    [Test]
    public void TheScan_CountsOverlappingOccurrences_AndCapsOnlyThePositions()
    {
        byte[] text = new byte[1000];
        Array.Fill(text, (byte)'a');
        var q = new HandoffQuery(0, true, 0, "aaaa"u8.ToArray());
        ScanResult r = RollingScan.Run(new MemoryStream(text), [q], 10, 64);
        Assert.That(r.Counts[0], Is.EqualTo(997));
        Assert.That(r.Positions[0], Is.EqualTo(Enumerable.Range(0, 10).Select(i => (long)i).ToList()));
    }

    [Test]
    public void TheScan_FindsNothing_ForAQueryLongerThanTheText_AndForNoQueries()
    {
        byte[] text = "ACGT"u8.ToArray();
        ScanResult r = RollingScan.Run(new MemoryStream(text), [new HandoffQuery(0, false, -1, "ACGTA"u8.ToArray())], 5);
        Assert.That(r.Counts[0], Is.Zero);
        Assert.That(RollingScan.Run(new MemoryStream(text), [], 5).Counts, Is.Empty);
    }

    [Test]
    public void TheScan_NeverReportsAWindowThatOnlyHasTheSameHash()
    {
        // two different windows of one length: the hash is of the whole window, the verification compares the bytes; a query that is not in the text must count 0 whatever the table holds
        byte[] text = RandomText(21, 20000, "ACGT");
        var qs = new List<HandoffQuery>();
        for (int i = 0; i < 300; i++) qs.Add(new HandoffQuery(i, false, -1, RandomText(1000 + i, 9, "ACGT")));
        ScanResult r = RollingScan.Run(new MemoryStream(text), qs, 5, 512);
        foreach (HandoffQuery q in qs) Assert.That(r.Counts[q.Index], Is.EqualTo(NaivePositions(text, q.Bytes).Count));
    }

    [Test]
    public void TheScan_NeverCountsAWindowWithTheSameHash_ForTwoStringsThatReallyCollide()
    {
        // The test above cannot fail: no two 9-letter strings share a 64-bit hash by chance. Two that do, by construction: the Thue-Morse string and its complement, 2^11 letters 'a' / 'b'.
        // The hash is a polynomial of an odd base modulo 2^64, the difference of the two hashes is the product of (1 - base^(2^j)) for j < 11, and every factor is even: its 2-adic valuation is at
        // least j + 2 for j >= 1 and 2 for j = 0 (the base ends in binary ...10101), 2 + (3 + 4 + ... + 12) = 77 > 64 in all, so the difference is 0 modulo 2^64. The text holds the complement twice; the query is the string itself.
        const int k = 11;
        int len = 1 << k;
        var tm = new byte[len];
        var complement = new byte[len];
        for (int i = 0; i < len; i++)
        {
            bool odd = (System.Numerics.BitOperations.PopCount((uint)i) & 1) == 1;
            tm[i] = odd ? (byte)'b' : (byte)'a';
            complement[i] = odd ? (byte)'a' : (byte)'b';
        }
        byte[] text = [.. complement, (byte)'x', .. complement];
        var same = new HandoffQuery(0, false, -1, complement);
        var collides = new HandoffQuery(1, false, -1, tm);
        Assert.That(tm, Is.Not.EqualTo(complement));
        ScanResult r = RollingScan.Run(new MemoryStream(text), [same, collides], 5, 1000);
        Assert.That(r.Counts[0], Is.EqualTo(2), "the complement itself is found twice");
        Assert.That(r.Counts[1], Is.Zero, "the string with the same hash is not in the text: the byte comparison must reject both windows");
        Assert.That(r.Positions[1], Is.Empty);
    }

    [Test]
    public void TheScan_WorksOnAStreamThatReturnsShortReads()
    {
        byte[] text = RandomText(22, 3000, "ACGT");
        var set = HandoffQuerySet.Build(text.Length, Reader(text), "ACGT"u8, [6, 14], 30, 10, 3);
        ScanResult whole = RollingScan.Run(new MemoryStream(text), set, 20, 256);
        ScanResult trickle = RollingScan.Run(new TrickleStream(text, 5), set, 20, 256);
        Assert.That(trickle.Counts, Is.EqualTo(whole.Counts));
        for (int i = 0; i < set.Count; i++) Assert.That(trickle.Positions[i], Is.EqualTo(whole.Positions[i]));
    }

    private sealed class TrickleStream(byte[] data, int most) : MemoryStream(data)
    {
        public override int Read(byte[] buffer, int offset, int count) => base.Read(buffer, offset, Math.Min(count, most));

        public override int Read(Span<byte> buffer) => base.Read(buffer[..Math.Min(buffer.Length, most)]);
    }

    // ---------------------------------------------------------------- the search over a suffix array

    [Test]
    public void TheSearch_GivesTheRankRangeOfTheSuffixesThatStartWithTheQuery()
    {
        foreach ((int seed, string alphabet, int n) in new[] { (31, "ACGT", 600), (32, "AB", 400), (33, "A", 60), (34, "ACGTN", 300) })
        {
            byte[] text = RandomText(seed, n, alphabet);
            long[] sa = NaiveSa(text);
            var set = HandoffQuerySet.Build(n, Reader(text), [.. alphabet.Select(c => (byte)c)], [1, 2, 3, 6, 11], 60, 40, (ulong)seed);
            foreach (HandoffQuery q in set)
            {
                (long lo, long hi) = SaSearch.Range(q.Bytes, n, r => sa[r], (pos, buf) => ReadAt(text, pos, buf));
                List<long> expected = NaivePositions(text, q.Bytes);
                Assert.That(hi - lo, Is.EqualTo(expected.Count), $"count of query {q.Index} on {alphabet}");
                Assert.That(sa[(int)lo..(int)hi].Order(), Is.EqualTo(expected), $"positions of query {q.Index}");
            }
        }
    }

    [Test]
    public void TheSearch_HandlesASuffixShorterThanTheQuery_TheEmptyTextAndAnEmptyRange()
    {
        byte[] text = "banana"u8.ToArray();
        long[] sa = NaiveSa(text);
        (long lo, long hi) = SaSearch.Range("nana"u8, text.Length, r => sa[r], (pos, buf) => ReadAt(text, pos, buf));
        Assert.That(hi - lo, Is.EqualTo(1));
        (lo, hi) = SaSearch.Range("nanab"u8, text.Length, r => sa[r], (pos, buf) => ReadAt(text, pos, buf));
        Assert.That(hi - lo, Is.Zero, "a query longer than every suffix that it could match");
        (lo, hi) = SaSearch.Range("a"u8, text.Length, r => sa[r], (pos, buf) => ReadAt(text, pos, buf));
        Assert.That(hi - lo, Is.EqualTo(3));
        (lo, hi) = SaSearch.Range("z"u8, text.Length, r => sa[r], (pos, buf) => ReadAt(text, pos, buf));
        Assert.That(hi - lo, Is.Zero);
        long[] saEmpty = [0];
        (lo, hi) = SaSearch.Range("a"u8, 0, r => saEmpty[r], (pos, buf) => 0);
        Assert.That(hi - lo, Is.Zero);
    }

    [Test]
    public void TheSearch_OfASwappedSuffixArray_DoesNotAgreeWithTheOracle()
    {
        byte[] text = RandomText(41, 800, "ACGT");
        long[] sa = NaiveSa(text);
        (sa[100], sa[500]) = (sa[500], sa[100]);
        var set = HandoffQuerySet.Build(text.Length, Reader(text), "ACGT"u8, [3, 5, 8], 100, 100, 41);
        ScanResult oracle = RollingScan.Run(new MemoryStream(text), set, 1000);
        int different = 0;
        foreach (HandoffQuery q in set)
        {
            (long lo, long hi) = SaSearch.Range(q.Bytes, text.Length, r => sa[r], (pos, buf) => ReadAt(text, pos, buf));
            if (hi - lo != oracle.Counts[q.Index]) different++;
        }
        Assert.That(different, Is.GreaterThan(0), "a binary search over an unsorted array must go wrong somewhere in 600 queries");
    }

    private static int ReadAt(byte[] text, long pos, byte[] buffer)
    {
        int count = (int)Math.Min(buffer.Length, text.Length - pos);
        text.AsSpan((int)pos, count).CopyTo(buffer);
        return count;
    }
}
