using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Vectors;

/// <summary>
/// Cheap, default-run checks of the committed libsais vector manifest (<c>analysis/sa_vectors.tsv</c>, work package 1a.6): it is well formed,
/// lists exactly the expected corpora, and the hashing convention (libsais <c>SA[0..n-1]</c> = our <c>SA[1..n]</c>, libsais <c>LCP[0..n-1]</c> =
/// our <c>LCP[1..n]</c>, little-endian int32) reproduces hashes that the libsais driver printed for a tiny text. The heavy comparison against
/// the inputs is <see cref="LibsaisVectorTests"/> (category <c>Vectors</c>, explicit).
/// </summary>
[TestFixture]
public class SaVectorsManifestTests
{
    [Test]
    public void Manifest_ListsExactlyTheExpectedCorpora()
    {
        var manifest = VectorManifest.Load();
        var expected = new List<(string Name, VectorManifest.CorpusKind Kind, int N)>();
        foreach (var (name, n) in VectorManifest.ExpectedFiles) expected.Add((name, VectorManifest.CorpusKind.File, n));
        foreach (var (name, n) in VectorManifest.ExpectedGenerated()) expected.Add((name, VectorManifest.CorpusKind.Generated, n));

        Assert.That(manifest.Entries.Select(e => (e.Corpus, e.Kind, e.N)), Is.EquivalentTo(expected));
        Assert.That(expected.Count(e => e.Kind == VectorManifest.CorpusKind.Generated), Is.EqualTo(2 * 37 + 1), "37 families at two sizes + all-distinct");
    }

    [Test]
    public void Manifest_WidthMatchesTheGeneratedText_AndEqualHashesMeanOrderIsomorphicTexts()
    {
        var manifest = VectorManifest.Load();
        var preparedByHash = new Dictionary<string, (string Corpus, byte[] Symbols)>();
        foreach (var e in manifest.Entries.Where(e => e.Kind == VectorManifest.CorpusKind.Generated && e.N == 100_000))
        {
            string text = VectorManifest.GenerateText(e.Corpus, e.N);
            Assert.That(text.Length, Is.EqualTo(e.N), e.Corpus);
            Assert.That(text.All(c => c <= 255) ? "u8" : "u16", Is.EqualTo(e.Width), e.Corpus);

            // Two corpora may share a suffix array only when their texts are the same up to an order-preserving remap
            // (random-protein and random-sigma20 draw the same ranks; periodic-p2-s1, periodic-aab-p2, akb-k1 and ab-repeat are all "abab...").
            PreparedText prepared = TextPreparation.Prepare(text);
            byte[] symbols = prepared.Width == SymbolWidth.Byte ? prepared.Bytes.ToArray() : [];
            if (preparedByHash.TryGetValue(e.SaSha256, out var other))
                Assert.That(prepared.Width == SymbolWidth.Byte && symbols.AsSpan().SequenceEqual(other.Symbols), Is.True,
                    $"{e.Corpus} and {other.Corpus} have the same SA hash but different texts");
            else
                preparedByHash.Add(e.SaSha256, (e.Corpus, symbols));
        }
        foreach (var e in manifest.Entries)
            Assert.That(e.SaSha256, Is.Not.EqualTo(e.LcpSha256), e.Corpus);
        Assert.That(manifest.Entries.Count(e => e.Kind == VectorManifest.CorpusKind.File), Is.EqualTo(4));
    }

    [Test]
    public void Parse_RejectsMalformedRows()
    {
        string h = new('a', 64);
        const string header = "corpus\tkind\tn\twidth\tsa_sha256\tlcp_sha256\n";
        Assert.DoesNotThrow(() => VectorManifest.Parse($"# c\n{header}x\tfile\t5\tu8\t{h}\t{h}\n"));
        Assert.Throws<FormatException>(() => VectorManifest.Parse($"x\tfile\t5\tu8\t{h}\t{h}\n"), "header missing");
        Assert.Throws<FormatException>(() => VectorManifest.Parse($"{header}x\tfile\t5\tu8\t{h}\n"), "field count");
        Assert.Throws<FormatException>(() => VectorManifest.Parse($"{header}x\tother\t5\tu8\t{h}\t{h}\n"), "kind");
        Assert.Throws<FormatException>(() => VectorManifest.Parse($"{header}x\tfile\t-5\tu8\t{h}\t{h}\n"), "n");
        Assert.Throws<FormatException>(() => VectorManifest.Parse($"{header}x\tfile\t5\tu32\t{h}\t{h}\n"), "width");
        Assert.Throws<FormatException>(() => VectorManifest.Parse($"{header}x\tfile\t5\tu8\t{h.ToUpperInvariant()}\t{h}\n"), "hash case");
        Assert.Throws<FormatException>(() => VectorManifest.Parse($"{header}x\tfile\t5\tu8\t{h}\t{h}\nx\tfile\t5\tu8\t{h}\t{h}\n"), "duplicate");
    }

    /// <summary>
    /// "banana": libsais gives SA = [5, 3, 1, 0, 4, 2] and LCP = [0, 1, 3, 0, 0, 2]; the hashes below were printed by the libsais driver
    /// (<c>libsais + libsais_plcp + libsais_lcp</c>, see the comment block of the manifest) for the 6-byte file.
    /// </summary>
    [Test]
    public void HashingConvention_ReproducesTheLibsaisDriverOnBanana()
    {
        var (sa, lcp) = EsaPipeline.Build("banana");
        Assert.That(sa, Is.EqualTo(new[] { 6, 5, 3, 1, 0, 4, 2 }), "n + 1 layout: SA[0] = n, SA[1..n] = the libsais SA");
        Assert.That(lcp.ToArray(), Is.EqualTo(new[] { 0, 0, 1, 3, 0, 0, 2 }), "LCP[0] = LCP[1] = 0, LCP[2..n] = libsais LCP[1..n-1]");
        Assert.That(EsaPipeline.HashSa(sa), Is.EqualTo("b2aab8610e2695af5a3dc5f079aa6e91215a77e56aef3b6bb678fcde3ea0983d"));
        Assert.That(EsaPipeline.HashLcp(lcp, 6), Is.EqualTo("a34ee68dd19d130c6668beb56b20879ae92f78bc98823a8fa8073768122795fe"));
    }
}
