using SuffixTree.Esa.External;
using SuffixTree.Esa.Tests.External;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Fuzz;

/// <summary>
/// The corpus of the corruption rows of work package S2.4 (explicit): the first 1,500 corruption cases of the differential harness for a master seed (<c>S24_CORPUS_SEED</c>, default 20261124): the text, the kind of damage and the damage seed
/// are those that <see cref="CorruptionFuzz.Draw"/> gives, the damaged container is the one that <see cref="CorruptionFuzz.Run"/> would write; each is written to <c>$S24_CORPUS_DIR/NNNN_Kind.sesa</c>. The runner of
/// <c>analysis/bench/corpus_s24</c> opens every file through the public API of a build, in Checked and in Trusted mode, one process each (a crash in Trusted mode is the documented behaviour of that mode).
/// </summary>
[TestFixture]
[Explicit("Writes 1,500 container files; run by analysis/bench/corpus_s24/run_corpus_s24.py (see PROGRESS-S2.md, work package S2.4)")]
public class CorruptionCorpusS24Tests
{
    [Test]
    public void WriteTheCorpus()
    {
        string? dir = Environment.GetEnvironmentVariable("S24_CORPUS_DIR");
        if (string.IsNullOrEmpty(dir)) Assert.Ignore("S24_CORPUS_DIR names the directory of the corpus");
        Directory.CreateDirectory(dir);
        ulong master = ulong.TryParse(Environment.GetEnvironmentVariable("S24_CORPUS_SEED"), out ulong s) ? s : 20261124UL;
        int written = 0;
        for (long i = 0; i < 1_500; i++)
        {
            (FuzzCase c, DamageKind kind, ulong damageSeed) = CorruptionFuzz.Draw(master, i);
            string text = c.Text();
            if (text.Length == 0 || text.Any(ch => ch > 255)) continue; // the harness skips them too
            var rng = new SplitMix64(damageSeed);
            bool psi = rng.NextInt(4) != 0;
            LcpChildCase clean = LcpChildCase.Of(text);
            byte[] image = SesaTestFile.Image(clean, psi, codec: (damageSeed & 4) != 0 ? ChildCodec.Tagged : ChildCodec.Packed);
            (byte[] damaged, _, _) = CorruptionFuzz.Damage(image, kind, rng, clean);
            File.WriteAllBytes(Path.Combine(dir, $"{i:D4}_{kind}.sesa"), damaged);
            written++;
        }
        TestContext.Out.WriteLine($"CORPUS written={written} master={master}");
        Assert.That(written, Is.GreaterThan(1_000));
    }
}
