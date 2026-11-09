using SuffixTree.Esa;
using SuffixTree.EmProto.Build;
using SuffixTree.EmProto.Harness;
using SuffixTree.EmProto.Text;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Build;

/// <summary>Dense symbols of a string family plus the oracle suffix array (SaisReference, which is itself pinned to the libsais vectors and to a naive sort).</summary>
internal sealed record DenseCase(string Name, byte[] Symbols, int Sigma)
{
    internal int MinBits => PackedText.BitsFor(Math.Max(Sigma, 2));

    internal uint[] ReferenceSa()
    {
        int[] sa = SaisReference.Build(new ByteSymbols(Symbols), Math.Max(Sigma, 1));
        return [.. sa.Select(x => (uint)x)];
    }
}

internal static class BuilderTestSupport
{
    internal static DenseCase Dense(string name, string text)
    {
        PreparedText p = TextPreparation.Prepare(text);
        if (p.Width != SymbolWidth.Byte) throw new ArgumentException("alphabet above 256 symbols");
        return new DenseCase(name, p.Bytes.ToArray(), p.AlphabetSize);
    }

    internal static IEnumerable<DenseCase> Families(int n, ulong seed)
    {
        foreach ((string name, string text) in StringFamilies.Adversarial(n, seed))
        {
            if (text.Any(c => c > 255) || text.Length == 0) continue;
            yield return Dense($"{name}-n{n}", text);
        }
    }

    internal static uint[] Build(DenseCase c, int bits, BuildOptions options, DcsRankTable? dcs = null)
    {
        using PackedText t = PackedText.Pack(c.Symbols, bits);
        using var sink = new MemorySaSink(c.Symbols.Length + 1L);
        _ = EmBuilder.Run(t, options, sink, dcs);
        return sink.ToArray();
    }

    internal static string Describe(uint[] expected, uint[] actual)
    {
        if (expected.Length != actual.Length) return $"length {actual.Length} instead of {expected.Length}";
        for (int i = 0; i < expected.Length; i++)
        {
            if (expected[i] != actual[i]) return $"first difference at rank {i}: {actual[i]} instead of {expected[i]}";
        }
        return "equal";
    }
}
