using System.Buffers.Binary;
using System.Globalization;
using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// Seeded fuzz of the validation levels of the reader (work package 1b.7f). Category <c>Fuzz</c>, <c>[Explicit]</c>; environment <c>ESA_FUZZ_SEED</c>, <c>ESA_FUZZ_SCALE</c> as for the other fuzz classes. A case is a text from the generator
/// families, its container, and one damage drawn from the kinds of the corruption tests: a flip of a byte anywhere (the levels must name it: the header page at <c>Header</c>, the body at <c>Sections</c>), a body word replaced and sealed again
/// (<c>Structure</c> or a plain answer), a truncation. An outcome other than "opens", <see cref="InvalidDataException"/> or <see cref="NotSupportedException"/>, and a flip that no level sees, fails the case.
/// </summary>
[TestFixture]
[Category("Fuzz")]
[Explicit("Fuzz of the validation levels; run with --filter FullyQualifiedName~SesaReaderFuzzTests (see PROGRESS-S1.md, Log 1b.7f)")]
public class SesaReaderFuzzTests
{
    [Test]
    public void Flips_UpTo4000Symbols_3000Cases() => Run(1, 4_000, 3_000, 81);

    [Test]
    public void Flips_UpTo60000Symbols_150Cases() => Run(4_001, 60_000, 150, 82);

    private static void Run(int minN, int maxN, long cases, ulong salt)
    {
        ulong seed = ulong.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SEED") ?? "1", CultureInfo.InvariantCulture);
        double scale = double.Parse(Environment.GetEnvironmentVariable("ESA_FUZZ_SCALE") ?? "1", CultureInfo.InvariantCulture);
        long count = Math.Max(1, (long)(cases * scale));
        int ran = 0;
        using var dir = new TempDir();
        string path = dir.File("f.sesa");
        for (long i = 0; i < count; i++)
        {
            ulong h = seed ^ (salt * 0xA24BAED4963EE407UL) ^ ((ulong)i * 0x9FB21C651E98DF25UL);
            var rng = new SplitMix64(new SplitMix64(h).Next());
            int n = minN + rng.NextInt(maxN - minN + 1);
            FuzzCase c = FuzzTexts.Pick(rng, n);
            string text = c.Text();
            if (text.Any(ch => ch > 255) || text.Length == 0) continue;
            byte[] clean = SesaTestFile.Image(LcpChildCase.Of(text), psi: rng.NextInt(4) != 0);
            byte[] image = (byte[])clean.Clone();
            int kind = rng.NextInt(4);
            string what;
            EsaValidationLevel level;
            switch (kind)
            {
                case 0:
                    {
                        int at = rng.NextInt(SesaFormat.PageBytes);
                        image[at] ^= (byte)(1 + rng.NextInt(255));
                        what = $"flip of page byte {at}";
                        level = EsaValidationLevel.Header;
                        break;
                    }
                case 1:
                    {
                        int at = SesaFormat.PageBytes + rng.NextInt(image.Length - SesaFormat.PageBytes);
                        image[at] ^= (byte)(1 + rng.NextInt(255));
                        what = $"flip of body byte {at}";
                        level = EsaValidationLevel.Sections;
                        break;
                    }
                case 2:
                    {
                        SesaOracle o = SesaOracle.FromBytes(image);
                        OracleSection s = o.Get(new[] { SesaFormat.KindSa, SesaFormat.KindLcp, SesaFormat.KindChild }[rng.NextInt(3)]);
                        int at = (int)s.Offset + rng.NextInt((int)s.Length);
                        image[at] ^= (byte)(1 + rng.NextInt(255));
                        SesaTestFile.Reseal(image);
                        what = $"sealed flip of byte {at} of section {s.Kind}";
                        level = EsaValidationLevel.Structure;
                        break;
                    }
                default:
                    {
                        int length = rng.NextInt(image.Length);
                        image = image[..length];
                        what = $"truncation to {length}";
                        level = EsaValidationLevel.Header;
                        break;
                    }
            }
            File.WriteAllBytes(path, image);
            string label = $"case {i} (n {text.Length}, {c.Id}): {what}";
            if (kind == 2)
            {
                // a sealed flip of the suffix array, the LCP array or the child table is seen by Structure, except where the flipped byte is one that the format leaves free (a child cell that holds no field changes to a
                // value, which is a difference too; an LCP escape slot is not free): every one of the three sections is fully determined by the text, so every flip must be seen
                Assert.That(() => SesaTestFile.Open(path, level).Dispose(), Throws.TypeOf<InvalidDataException>(), label);
            }
            else if (kind == 3)
            {
                Assert.That(() => SesaTestFile.Open(path, EsaValidationLevel.Header).Dispose(), Throws.TypeOf<InvalidDataException>(), label);
            }
            else
            {
                Assert.That(() => SesaTestFile.Open(path, level).Dispose(), Throws.TypeOf<InvalidDataException>().Or.TypeOf<NotSupportedException>(), label);
            }
            ran++;
        }
        Assert.That(ran, Is.GreaterThan(0));
    }
}
