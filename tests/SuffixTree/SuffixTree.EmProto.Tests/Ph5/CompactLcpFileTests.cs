using SuffixTree.EmProto.Ph5;
using SuffixTree.Esa;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph5;

/// <summary>The reader of the compact LCP of PH4 (base file + overflow pairs): resolved values forwards (the scan of layout A, PSV) and backwards (the scan of the NSV).</summary>
[TestFixture]
public class CompactLcpFileTests
{
    private string _dir = "";

    [SetUp]
    public void SetUp()
    {
        _dir = Path.Combine(Path.GetTempPath(), "clf-" + Guid.NewGuid().ToString("N"));
        Directory.CreateDirectory(_dir);
    }

    [TearDown]
    public void TearDown() => Directory.Delete(_dir, true);

    private static IEnumerable<(string Name, string Text, LcpKind Kind)> Cases()
    {
        yield return ("random-dna", StringFamilies.RandomDna(3000, 5), LcpKind.U8);
        yield return ("u8-with-escapes", ChildSamples.U8WithEscapes(), LcpKind.U8);
        yield return ("fibonacci", StringFamilies.Fibonacci(4000), LcpKind.U8);
        yield return ("runs", StringFamilies.Runs(5000, 256, 3), LcpKind.U8);
        yield return ("all-same-70k", StringFamilies.AllSame(70_000), LcpKind.U16);
        yield return ("all-same-150k", StringFamilies.AllSame(150_000), LcpKind.U32);
        yield return ("period-300", string.Concat(Enumerable.Repeat(StringFamilies.RandomDna(300, 9), 40)), LcpKind.U16);
    }

    [Test]
    public void Forward_And_Backward_YieldTheDenseValues_InEveryWidth()
    {
        var kinds = new HashSet<LcpKind>();
        foreach ((string name, string text, _) in Cases())
        {
            ChildReference r = ChildReference.Of(text);
            kinds.Add(r.Lcp.Kind);
            (string b, string o) = r.WriteLcpFiles(_dir, name);
            uint[] dense = ChildFiles.Dense(r.Lcp);
            foreach (int buffer in new[] { 8, 12, 100, 1 << 16 })
            {
                using var file = new CompactLcpFile(b, o, r.N, buffer);
                Assert.That(file.OverflowCount, Is.EqualTo(r.Lcp.OverflowCount), $"{name}: overflow pairs");
                var forward = new List<uint>();
                long expectedFirst = 0;
                file.ForEachForward((first, values) =>
                {
                    Assert.That(first, Is.EqualTo(expectedFirst), $"{name}: forward blocks are contiguous");
                    forward.AddRange(values.ToArray());
                    expectedFirst += values.Length;
                });
                Assert.That(forward, Is.EqualTo(dense), $"{name}, buffer {buffer}: forward");
                var backward = new List<uint>();
                long expectedEnd = r.N + 1;
                file.ForEachBackward((first, values) =>
                {
                    Assert.That(first + values.Length, Is.EqualTo(expectedEnd), $"{name}: backward blocks are contiguous");
                    expectedEnd = first;
                    backward.InsertRange(0, values.ToArray());
                });
                Assert.That(backward, Is.EqualTo(dense), $"{name}, buffer {buffer}: backward");
            }
        }
        Assert.That(kinds, Is.EquivalentTo(new[] { LcpKind.U8, LcpKind.U16, LcpKind.U32 }), "the cases cover the three widths");
    }

    [Test]
    public void WrongLengths_AndCorruptOverflow_AreRefused()
    {
        ChildReference r = ChildReference.Of(ChildSamples.U8WithEscapes());
        (string b, string o) = r.WriteLcpFiles(_dir, "bad");
        Assert.That(r.Lcp.Kind, Is.EqualTo(LcpKind.U8));
        Assert.That(r.Lcp.OverflowCount, Is.GreaterThan(2), "the case has escapes");
        Assert.Throws<ArgumentException>(() => new CompactLcpFile(b, o, r.N + 7));
        byte[] pairs = File.ReadAllBytes(o);
        pairs[0] ^= 1; // the first pair names another rank
        string o2 = Path.Combine(_dir, "bad2.ovf");
        File.WriteAllBytes(o2, pairs);
        using var file = new CompactLcpFile(b, o2, r.N);
        Assert.Throws<InvalidDataException>(() => file.ForEachForward((_, _) => { }));
        Assert.Throws<InvalidDataException>(() => file.ForEachBackward((_, _) => { }));
        string o3 = Path.Combine(_dir, "bad3.ovf");
        File.WriteAllBytes(o3, pairs.AsSpan(0, pairs.Length - 8).ToArray()); // one pair missing
        using var missing = new CompactLcpFile(b, o3, r.N);
        Assert.Throws<InvalidDataException>(() => missing.ForEachForward((_, _) => { }));
    }
}
