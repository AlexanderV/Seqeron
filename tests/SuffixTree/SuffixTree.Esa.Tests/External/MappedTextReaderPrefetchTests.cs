using SuffixTree.Esa.External;

namespace SuffixTree.Esa.Tests.External;

/// <summary>The cache hint that PH4 sends a few iterations ahead (work package 1b.7c): it never throws, never changes a result and refuses nothing.</summary>
[TestFixture]
public class MappedTextReaderPrefetchTests
{
    [Test]
    public void Prefetch_AcceptsAnyPosition_AndChangesNothing()
    {
        var rng = new SplitMix64Local(5);
        foreach (int bits in new[] { 1, 2, 3, 5, 8 })
        {
            int n = 500;
            var symbols = new byte[n];
            for (int i = 0; i < n; i++) symbols[i] = (byte)(rng.Next() % (1u << bits));
            using var text = new PinnedText(symbols, bits);
            MappedTextReader r = text.Reader;
            byte[] before = (byte[])text.FileBytes.Clone();
            foreach (long i in new long[] { 0, 1, 7, 8, 63, 64, n - 1, n, n + 1, n + 100, long.MaxValue, long.MinValue, -1, -1000, uint.MaxValue, 1L << 40 })
            {
                r.Prefetch(i);
            }
            for (int k = 0; k < 2000; k++) r.Prefetch((long)(rng.Next() % 1000) - 200);
            Assert.That(text.FileBytes, Is.EqualTo(before), $"bits {bits}: the file is untouched");
            for (int i = 0; i < n; i++) Assert.That(r.Symbol(i), Is.EqualTo((uint)symbols[i]), $"bits {bits} symbol {i}");
        }
    }

    [Test]
    public void Prefetch_AfterDispose_IsHarmless()
    {
        var text = new PinnedText(new byte[] { 0, 1, 1, 0 }, 1);
        MappedTextReader r = text.Reader;
        text.Dispose();
        Assert.DoesNotThrow(() =>
        {
            r.Prefetch(0);
            r.Prefetch(2);
        });
    }

    private sealed class SplitMix64Local(ulong seed)
    {
        private ulong _s = seed;

        internal uint Next()
        {
            _s += 0x9E3779B97F4A7C15UL;
            ulong z = _s;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9UL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBUL;
            return (uint)((z ^ (z >> 31)) >> 7);
        }
    }
}
