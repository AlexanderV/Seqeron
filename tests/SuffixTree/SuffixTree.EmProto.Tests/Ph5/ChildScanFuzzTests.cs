using System.Globalization;
using SuffixTree.EmProto.Ph5;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph5;

/// <summary>
/// Differential fuzz of PH5 with tiny spill blocks (block sizes 2 to 50 entries, so that every spill, refill and eviction path is exercised on texts of 10^3 to 10^5 symbols), tiny windows (16 bytes to 4 KiB) and one to four
/// window slots, both codecs, the sequential and the concurrent NSV pass, against the in-RAM structures of the product. Category <c>Fuzz</c>, <c>[Explicit]</c>; environment <c>EMPROTO_FUZZ_SEED</c> (default 1) and
/// <c>EMPROTO_FUZZ_SECONDS</c> (the budget, default 120).
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Differential fuzz of PH5; run with --filter Category=Fuzz (see PROGRESS-S1.md, Log E5)")]
public class ChildScanFuzzTests
{
    [Test]
    public void TinyBlocksAndWindows_Texts1kTo100k_EqualTheInRamStructures()
    {
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("EMPROTO_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double seconds = double.Parse(Environment.GetEnvironmentVariable("EMPROTO_FUZZ_SECONDS") ?? "120", CultureInfo.InvariantCulture);
        string dir = Path.Combine(Path.GetTempPath(), "emproto-ph5-fuzz-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(dir);
        try
        {
            var clock = System.Diagnostics.Stopwatch.StartNew();
            var rng = new SplitMix64(seed * 0x9E3779B97F4A7C15UL);
            long cases = 0, spilling = 0, parallel = 0, packed = 0;
            long maxStack = 0;
            while (clock.Elapsed.TotalSeconds < seconds)
            {
                string text;
                int pick = rng.NextInt(30);
                if (pick == 0)
                {
                    string x = StringFamilies.RandomDna(2000 + rng.NextInt(3000), rng.Next());
                    text = x + x + x;
                }
                else if (pick == 1)
                {
                    text = new string('a', 1000 + rng.NextInt(20_000)) + (rng.NextInt(2) == 0 ? "b" : string.Empty);
                }
                else
                {
                    int n = new[] { 1000, 3000, 10_000, 30_000, 100_000 }[rng.NextInt(5)];
                    var all = StringFamilies.Adversarial(n, rng.Next()).Where(f => f.Text.Length > 0).ToList();
                    text = all[rng.NextInt(all.Count)].Text;
                }
                ChildReference r = ChildReference.Of(text);
                ChildCodec codec = rng.NextInt(2) == 0 ? ChildCodec.Tagged : ChildCodec.Packed;
                bool par = rng.NextInt(3) == 0;
                ChildBuildOptions o = ChildChecks.Tiny(Path.Combine(dir, "scratch"), 2 + rng.NextInt(49), 16 << rng.NextInt(9), 1 + rng.NextInt(4), codec, flush: rng.NextInt(4) != 0) with
                {
                    ReadBufferBytes = 8 + rng.NextInt(900),
                    SpillFlushBytes = rng.NextInt(2) == 0 ? 0 : 50 + rng.NextInt(2000),
                    ParallelBackward = par,
                };
                ChildBuildStats st = ChildChecks.RunAndCompare(r, o, dir, $"fuzz case {cases} (seed {seed}): n={r.N} block={o.BlockEntries} window={o.WindowBytes} slots={o.WindowSlots} {codec} parallel={par}", both: codec == ChildCodec.Tagged);
                cases++;
                if (st.Forward.Spill.BlocksWritten > 0 || (st.Backward?.Spill.BlocksWritten ?? 0) > 0) spilling++;
                if (par && codec == ChildCodec.Tagged) parallel++;
                if (codec == ChildCodec.Packed) packed++;
                maxStack = Math.Max(maxStack, st.Forward.StackPeak);
            }
            TestContext.Out.WriteLine($"fuzz: {cases} cases ({spilling} with a spill, {parallel} concurrent NSV, {packed} packed), deepest stack {maxStack}, {clock.Elapsed.TotalSeconds:F0} s");
            Assert.That(cases, Is.GreaterThan(0));
            Assert.That(spilling, Is.GreaterThan(cases / 4), "a large share of the cases leaves the two blocks of RAM");
        }
        finally
        {
            Directory.Delete(dir, true);
        }
    }
}
