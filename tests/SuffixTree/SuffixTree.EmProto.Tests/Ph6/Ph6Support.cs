using System.Runtime.InteropServices;
using SuffixTree.EmProto.Ph6;
using SuffixTree.EmProto.Text;
using SuffixTree.EmProto.Tests.WholePipeline;
using SuffixTree.Esa;
using SuffixTree.Esa.Bench.Links;
using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.EmProto.Tests.Ph6;

/// <summary>A scratch directory of one test.</summary>
internal sealed class Ph6Dir : IDisposable
{
    internal Ph6Dir() => Path = System.IO.Path.Combine(System.IO.Path.GetTempPath(), "emproto-e9-" + Guid.NewGuid().ToString("N"));

    internal string Path { get; }

    internal string File(string name) => System.IO.Path.Combine(EnsureCreated(), name);

    private string EnsureCreated()
    {
        Directory.CreateDirectory(Path);
        return Path;
    }

    public void Dispose()
    {
        if (Directory.Exists(Path)) Directory.Delete(Path, true);
    }
}

/// <summary>A text in the form the streaming code consumes, plus the in-RAM oracle of everything it must produce.</summary>
internal sealed class Ph6Case
{
    internal required string Text { get; init; }

    internal required byte[] Symbols { get; init; }

    internal required int Sigma { get; init; }

    internal required int[] Sa { get; init; }

    internal required int[] IsaRam { get; init; }

    /// <summary>P1's <c>psi</c>: entry 0 is -1.</summary>
    internal required int[] PsiRam { get; init; }

    internal int N => Text.Length;

    /// <summary>Row of the suffix 0 (its BWT symbol is the sentinel).</summary>
    internal int Primary => Array.IndexOf(Sa, 0);

    internal long[] Counts()
    {
        var c = new long[Sigma];
        foreach (byte s in Symbols) c[s]++;
        return c;
    }

    /// <summary>The expected PSI file: P1's table with entry 0 replaced by 0 (the convention) or by the row of the suffix 0 (the cyclic Psi).</summary>
    internal byte[] ExpectedPsiBytes(bool cyclic)
    {
        var psi = (int[])PsiRam.Clone();
        psi[0] = cyclic ? Primary : 0;
        return MemoryMarshal.AsBytes(psi.AsSpan()).ToArray();
    }

    internal byte[] ExpectedIsaBytes() => MemoryMarshal.AsBytes(IsaRam.AsSpan()).ToArray();

    /// <summary>BWT by the definition, one byte per row; the row of the suffix 0 holds the dummy 0.</summary>
    internal byte[] ExpectedBwt()
    {
        var bwt = new byte[N + 1];
        for (int j = 0; j <= N; j++) bwt[j] = Sa[j] == 0 ? (byte)0 : Symbols[Sa[j] - 1];
        return bwt;
    }

    /// <summary>The text with its in-RAM suffix array (<see cref="EsaIndex"/>, SA-IS) and the ISA / Psi of the prototype P1 (<c>LinkTables</c>).</summary>
    internal static Ph6Case Of(string text)
    {
        (byte[] symbols, char[] letters) = PipelineChecks.Dense(text);
        EsaIndex index = EsaIndex.Build(text, ChildLayoutKind.Tagged);
        int[] sa = index.SuffixArray;
        int[] isa = LinkTables.BuildIsa(sa);
        return new Ph6Case { Text = text, Symbols = symbols, Sigma = letters.Length, Sa = sa, IsaRam = isa, PsiRam = LinkTables.BuildPsi(sa, isa) };
    }

    /// <summary>Same, with the definitional suffix array (sorted suffixes): only for texts of a few hundred symbols.</summary>
    internal static Ph6Case OfNaive(string text)
    {
        (byte[] symbols, char[] letters) = PipelineChecks.Dense(text);
        int[] sa = NaiveSuffixArray.Build(text);
        int[] isa = LinkTables.BuildIsa(sa);
        return new Ph6Case { Text = text, Symbols = symbols, Sigma = letters.Length, Sa = sa, IsaRam = isa, PsiRam = LinkTables.BuildPsi(sa, isa) };
    }

    internal string WriteSa(Ph6Dir dir, string name = "x.sa")
    {
        string path = dir.File(name);
        File.WriteAllBytes(path, MemoryMarshal.AsBytes(Sa.AsSpan()).ToArray());
        return path;
    }

    /// <summary>The packed text file and its <c>.meta</c> as <c>pack</c> writes them.</summary>
    internal string WritePacked(Ph6Dir dir, string name = "x.pk")
    {
        string path = dir.File(name);
        int bits = Math.Max(1, PackedText.BitsFor(Sigma));
        var bytes = new byte[PackedText.TotalBytes(N, bits)];
        unsafe
        {
            fixed (byte* p = bytes) PackedText.PackInto(Symbols, bits, p);
        }
        File.WriteAllBytes(path, bytes);
        File.WriteAllText(path + ".meta", $"n={N}\nbits={bits}\nsigma={Sigma}\n");
        return path;
    }

    /// <summary>A text with exactly <paramref name="sigma"/> distinct letters (char codes 0..sigma-1) of length n &gt;= sigma: every letter once in a random order first, then random letters.</summary>
    internal static string RandomText(int n, int sigma, ulong seed)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(n, sigma);
        var rng = new SplitMix64(seed);
        var letters = Enumerable.Range(0, sigma).ToArray();
        for (int i = sigma - 1; i > 0; i--)
        {
            int j = rng.NextInt(i + 1);
            (letters[i], letters[j]) = (letters[j], letters[i]);
        }
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = (char)(i < sigma ? letters[i] : rng.NextInt(sigma));
        return new string(c);
    }
}
