using System.Buffers.Binary;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// The corpus of the corruption rows of work package S2.5 (explicit): containers that carry the sections that the navigation pack reads in place, damaged. Two families, written to <c>$S25_CORPUS_DIR</c>:
/// <c>A_NNNN_Kind.sesa</c> are the cases of the differential harness (the same draw and the same damages as <see cref="CorruptionCorpusS24Tests"/>, another master seed) over a container that has the RMQ section in three cases out of four and the
/// PSI section in three out of four; <c>P_NNNN_What.sesa</c> are damages aimed at the two sections of this package: a byte of PSI or of RMQ flipped, a share of their words replaced by extreme values, two PSI entries exchanged, PSI rotated by
/// one rank, the RMQ levels cleared, one field of the RMQ table entry changed (each sealed, so that the checksums do not stop them). The runner of <c>analysis/bench/corpus_s24</c> opens every file through the public API of a build in Checked and
/// in Trusted mode (one process each); the rows are in PROGRESS-S2.md (work package S2.5).
/// </summary>
[TestFixture]
[Explicit("Writes 2,100 container files; run by analysis/bench/corpus_s24/run_corpus_s24.py (see PROGRESS-S2.md, work package S2.5)")]
public class CorruptionCorpusS25Tests
{
    private const int FamilyA = 1_500;
    private const int FamilyP = 600;

    [Test]
    public void WriteTheCorpus()
    {
        string? dir = Environment.GetEnvironmentVariable("S25_CORPUS_DIR");
        if (string.IsNullOrEmpty(dir)) Assert.Ignore("S25_CORPUS_DIR names the directory of the corpus");
        Directory.CreateDirectory(dir);
        ulong master = ulong.TryParse(Environment.GetEnvironmentVariable("S25_CORPUS_SEED"), out ulong s) ? s : 20261125UL;
        int a = 0;
        for (long i = 0; a < FamilyA && i < 40_000; i++)
        {
            (FuzzCase c, DamageKind kind, ulong damageSeed) = CorruptionFuzz.Draw(master, i);
            string text = c.Text();
            if (text.Length == 0 || text.Any(ch => ch > 255)) continue;
            var rng = new SplitMix64(damageSeed);
            bool psi = rng.NextInt(4) != 0;
            bool rmq = rng.NextInt(4) != 0;
            LcpChildCase clean = LcpChildCase.Of(text);
            byte[] image = SesaTestFile.Image(clean, psi, codec: (damageSeed & 4) != 0 ? ChildCodec.Tagged : ChildCodec.Packed, rmq: rmq);
            (byte[] damaged, _, _) = CorruptionFuzz.Damage(image, kind, rng, clean);
            File.WriteAllBytes(Path.Combine(dir, $"A_{i:D4}_{kind}.sesa"), damaged);
            a++;
        }
        int p = 0;
        for (long i = 0; p < FamilyP && i < 40_000; i++)
        {
            (FuzzCase c, _, ulong damageSeed) = CorruptionFuzz.Draw(master ^ 0x5DEECE66DUL, i);
            string text = c.Text();
            if (text.Length < 3 || text.Any(ch => ch > 255)) continue;
            var rng = new SplitMix64(damageSeed);
            LcpChildCase clean = LcpChildCase.Of(text);
            byte[] image = SesaTestFile.Image(clean, psi: true, codec: (damageSeed & 4) != 0 ? ChildCodec.Tagged : ChildCodec.Packed, rmq: true);
            (byte[] damaged, string what) = DamageNavigationSections(image, p % 8, rng);
            File.WriteAllBytes(Path.Combine(dir, $"P_{p:D4}_{what}.sesa"), damaged);
            p++;
        }
        TestContext.Out.WriteLine($"CORPUS written={a + p} (A={a}, P={p}) master={master}");
        Assert.That(a + p, Is.EqualTo(FamilyA + FamilyP));
    }

    /// <summary>One damage of the second family, sealed (the checksums are recomputed after it).</summary>
    private static (byte[] Image, string What) DamageNavigationSections(byte[] clean, int variant, SplitMix64 rng)
    {
        byte[] image = (byte[])clean.Clone();
        SesaOracle o = SesaOracle.FromBytes(image);
        OracleSection psi = o.Get(SesaFormat.KindPsi);
        OracleSection rmq = o.Get(SesaFormat.KindRmq);
        long n = o.TextLen;
        string what;
        switch (variant)
        {
            case 0:
            {
                int at = (int)psi.Offset + rng.NextInt((int)psi.Length);
                image[at] ^= (byte)(1 + rng.NextInt(255));
                what = "PsiByteFlip";
                break;
            }
            case 1:
            {
                WordsOf(image, psi, n, rng);
                what = "PsiWords";
                break;
            }
            case 2:
            {
                int r = rng.NextInt((int)Math.Max(1, psi.Length / 4 - 1));
                int at = (int)psi.Offset + (4 * r);
                int x = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(at));
                int y = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan(at + 4));
                BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(at), y);
                BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan(at + 4), x);
                what = "PsiSwap";
                break;
            }
            case 3:
            {
                int words = (int)(psi.Length / 4);
                var rotated = new int[words];
                for (int w = 0; w < words; w++) rotated[(w + 1) % words] = BinaryPrimitives.ReadInt32LittleEndian(image.AsSpan((int)psi.Offset + (4 * w)));
                for (int w = 0; w < words; w++) BinaryPrimitives.WriteInt32LittleEndian(image.AsSpan((int)psi.Offset + (4 * w)), rotated[w]);
                what = "PsiRotated";
                break;
            }
            case 4:
            {
                int at = (int)rmq.Offset + rng.NextInt((int)rmq.Length);
                image[at] ^= (byte)(1 + rng.NextInt(255));
                what = "RmqByteFlip";
                break;
            }
            case 5:
            {
                WordsOf(image, rmq, n, rng);
                what = "RmqWords";
                break;
            }
            case 6:
            {
                image.AsSpan((int)rmq.Offset, (int)rmq.Length).Clear();
                what = "RmqCleared";
                break;
            }
            default:
            {
                int entry = 128 + (64 * rmq.Index);
                int field = rng.NextInt(5) switch { 0 => 24, 1 => 32, 2 => 40, 3 => 44, _ => 48 };
                image[entry + field + rng.NextInt(4)] ^= (byte)(1 + rng.NextInt(255));
                what = "RmqEntryField";
                break;
            }
        }
        SesaTestFile.Reseal(image);
        return (image, what);
    }

    private static void WordsOf(byte[] image, OracleSection s, long n, SplitMix64 rng)
    {
        int words = (int)(s.Length / 4);
        if (words == 0) return;
        uint[] extreme = [0, 1, (uint)n, (uint)n + 1, (uint)n + 2, 0x7FFFFFFF, 0x80000000, 0xFFFFFFFF, 0xFFFFFFFE, 0x40000000, 0x3FFFFFFF, 0xC0000001];
        int hits = 1 + rng.NextInt(Math.Max(1, words / 5));
        for (int h = 0; h < hits; h++)
        {
            int w = rng.NextInt(words);
            uint v = rng.NextInt(3) switch { 0 => extreme[rng.NextInt(extreme.Length)], 1 => (uint)rng.Next(), _ => (uint)rng.NextInt((int)Math.Max(2, n + 2)) };
            BinaryPrimitives.WriteUInt32LittleEndian(image.AsSpan((int)s.Offset + (4 * w)), v);
        }
    }
}
