using System.Runtime.InteropServices;
using System.Security.Cryptography;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Ph6;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;
using Prod = global::SuffixTree.Esa.External;

namespace SuffixTree.EmProto.Tests.Ph6;

/// <summary>
/// Gate G2 of work package 1b.7e for PH6: the production BWT pass and Psi distribution write the same bytes as the prototype of E9 (<c>BwtPass</c>, <c>PsiStream</c>, code unchanged) with the same counted figures,
/// on generated texts and on the data files of E9 (<c>/tmp/stbench/e9</c>, skipped when absent), which also carry the in-RAM reference tables of the E9 oracle.
/// </summary>
[TestFixture]
public unsafe class ProductionPh6Tests
{
    private static string E9Dir => Path.Combine(Environment.GetEnvironmentVariable("SEQERON_BENCH_DIR") ?? "/tmp/stbench", "e9");

    /// <summary>Runs the production stages on a suffix array at offset 64 of a container-like file; returns the BWT bytes, the result and the PSI entries.</summary>
    internal static (byte[] Bwt, Prod.BwtResult Result, Prod.PsiStats Psi, byte[] PsiBytes) RunProduction(Ph6Dir dir, string saPath, string packedPath, long n, int bits, int sigma, Prod.BwtOptions? bwt = null, Prod.PsiOptions? psi = null)
    {
        byte[] saBytes = File.ReadAllBytes(saPath);
        long psiOffset = (64 + saBytes.Length + 63) / 64 * 64;
        string container = dir.File("container.bin");
        using Prod.MappedOutputFile file = Prod.MappedOutputFile.Create(container, psiOffset + (4 * (n + 1)), new Prod.MappedFileOptions { Reservation = Prod.ReservationMode.Sparse });
        using (Prod.MappedWindowWriter w = file.CreateWriter(64, saBytes.Length, 1 << 20)) w.Append(saBytes);
        file.Complete(false);
        using Prod.MappedTextReader text = Prod.MappedTextReader.Open(packedPath, n, bits);
        string bwtPath = dir.File("prod.bwt");
        Prod.BwtResult r = Prod.BwtPass.Build(container, 64, text, sigma, bwtPath, bwt);
        Prod.PsiStats stats = Prod.PsiStream.Build(bwtPath, n, r.PrimaryRow, r.Counts, file, psiOffset, psi);
        file.Complete(false);
        byte[] psiBytes = new byte[4 * (n + 1)];
        using (var h = File.OpenHandle(container)) RandomAccess.Read(h, psiBytes, psiOffset);
        return (File.ReadAllBytes(bwtPath), r, stats, psiBytes);
    }

    private static void CompareWithPrototype(Ph6Case c, string label, int window = 1 << 20, int group = 64)
    {
        using var dir = new Ph6Dir();
        string sa = c.WriteSa(dir);
        string pk = c.WritePacked(dir);
        int bits = Math.Max(1, PackedText.BitsFor(c.Sigma));
        string protoBwt = dir.File("proto.bwt");
        BwtMeta meta;
        using (var mapped = new MappedFile(pk, writable: false))
        using (var text = new PackedText(mapped.Pointer, c.N, bits, mapped, ownsOwner: false))
        {
            meta = BwtPass.Build(sa, text, c.Sigma, protoBwt, new BwtOptions { WindowBytes = window });
        }
        string protoPsi = dir.File("proto.psi");
        PsiStats protoStats = PsiStream.Build(protoBwt, meta, protoPsi, new PsiOptions { WindowBytes = window, SymbolsPerPass = group });
        (byte[] bwt, Prod.BwtResult r, Prod.PsiStats stats, byte[] psi) = RunProduction(dir, sa, pk, c.N, bits, c.Sigma, new Prod.BwtOptions { WindowBytes = window }, new Prod.PsiOptions { WindowBytes = window, SymbolsPerPass = group });
        Assert.That(bwt, Is.EqualTo(File.ReadAllBytes(protoBwt)), $"{label}: BWT");
        Assert.That(psi, Is.EqualTo(File.ReadAllBytes(protoPsi)), $"{label}: PSI");
        Assert.That((r.N, r.PrimaryRow), Is.EqualTo((meta.N, meta.PrimaryRow)), $"{label}: meta");
        Assert.That(r.Counts, Is.EqualTo(meta.Counts), $"{label}: counts");
        Assert.That((stats.Passes, stats.SymbolsPerPass, stats.Windows, stats.MappedBoundBytes, stats.BwtBytesRead, stats.PsiBytesWritten), Is.EqualTo((protoStats.Passes, protoStats.SymbolsPerPass, protoStats.Windows, protoStats.MappedBoundBytes, protoStats.BwtBytesRead, protoStats.PsiBytesWritten)), $"{label}: the figures of the pass");
        Assert.That(psi, Is.EqualTo(c.ExpectedPsiBytes(false)), $"{label}: PSI equals the in-RAM tables of P1");
    }

    [Test]
    public void RandomTexts_AllAlphabets_EqualThePrototype_AtTheDefaultAndSmallWindows()
    {
        int v = 0;
        foreach (int sigma in new[] { 1, 2, 4, 5, 20, 64, 65, 200, 256 })
        {
            foreach ((int window, int group) in new[] { (1 << 20, 64), (64, 7), (4, 1) })
            {
                v++;
                int n = sigma + 1 + ((v * 53) % 400);
                Ph6Case c = Ph6Case.Of(Ph6Case.RandomText(n, sigma, (ulong)(77 + v)));
                CompareWithPrototype(c, $"sigma {sigma} n {n} window {window} group {group}", window, group);
            }
        }
    }

    [Test]
    public void TheFamilies_AtTenToTheFive_EqualThePrototype()
    {
        foreach ((string name, string text) in new[] { ("a^n", StringFamilies.AllSame(100_000)), ("fibonacci", StringFamilies.Fibonacci(100_000)), ("thue-morse", StringFamilies.ThueMorse(100_000)), ("dna", StringFamilies.RandomDna(100_000, 3)) })
        {
            CompareWithPrototype(Ph6Case.Of(text), $"{name} 10^5");
        }
    }

    [Test]
    public void TheDataFilesOfE9_GiveTheReferenceBwtAndPsi()
    {
        int checkedFiles = 0;
        foreach (string name in new[] { "a_1m", "fib_1m", "tm_1m", "prot_1m", "rnd256_10m", "rep_10m", "dna_10m" })
        {
            string pk = Path.Combine(E9Dir, name + ".pk"), sa = Path.Combine(E9Dir, name + ".sa"), bwtRef = Path.Combine(E9Dir, name + ".bwt"), psiRef = Path.Combine(E9Dir, name + ".ref.psi");
            if (!File.Exists(pk) || !File.Exists(sa) || !File.Exists(bwtRef) || !File.Exists(psiRef)) continue;
            string[] meta = File.ReadAllLines(pk + ".meta");
            long n = long.Parse(meta.First(l => l.StartsWith("n=", StringComparison.Ordinal))[2..], System.Globalization.CultureInfo.InvariantCulture);
            int bits = int.Parse(meta.First(l => l.StartsWith("bits=", StringComparison.Ordinal))[5..], System.Globalization.CultureInfo.InvariantCulture);
            int sigma = int.Parse(meta.First(l => l.StartsWith("sigma=", StringComparison.Ordinal))[6..], System.Globalization.CultureInfo.InvariantCulture);
            using var dir = new Ph6Dir();
            (byte[] bwt, _, _, byte[] psi) = RunProduction(dir, sa, pk, n, bits, sigma);
            Assert.That(SHA256.HashData(bwt), Is.EqualTo(SHA256.HashData(File.ReadAllBytes(bwtRef))), $"{name}: BWT");
            Assert.That(SHA256.HashData(psi), Is.EqualTo(SHA256.HashData(File.ReadAllBytes(psiRef))), $"{name}: PSI");
            checkedFiles++;
        }
        if (checkedFiles == 0) Assert.Ignore($"no reference files in {E9Dir}");
    }
}
