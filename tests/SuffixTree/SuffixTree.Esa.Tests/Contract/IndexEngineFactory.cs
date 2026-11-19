using SuffixTree.Esa.Tests.External;
using SuffixTree.Indexing;

namespace SuffixTree.Esa.Tests.Contract;

/// <summary>
/// A way to build one engine over a text, for the engine conformance suite (<see cref="IndexConformanceTests"/>): the suite is parameterised
/// over these factories, so a new engine joins by adding one entry to <see cref="IndexEngineFactories.All"/>.
/// </summary>
/// <param name="Name">Display name, also the NUnit fixture name.</param>
/// <param name="Create">Builds the index of a text.</param>
/// <param name="Accepts">Whether the engine can hold the text (null: every text); the container of 1b.7f holds code units below 256 only, and the suite skips the other texts for it.</param>
public sealed record IndexEngineFactory(string Name, Func<string, ISuffixIndex> Create, Func<string, bool>? Accepts = null)
{
    /// <summary>Whether <see cref="Create"/> can be called for <paramref name="text"/>.</summary>
    public bool AcceptsText(string text) => Accepts?.Invoke(text) ?? true;

    /// <inheritdoc />
    public override string ToString() => Name;
}

/// <summary>The registered engines. The Ukkonen tree behind a test-only adapter is the oracle engine (L2 of the oracle hierarchy).</summary>
public static class IndexEngineFactories
{
    /// <summary>Every engine the conformance suite runs against.</summary>
    public static IEnumerable<IndexEngineFactory> All { get; } =
    [
        new IndexEngineFactory("ukkonen-adapter", text => new UkkonenIndexAdapter(text)),
        // the enhanced suffix array of work package 1b.5, with each of the three child-table layouts (the default is Tagged)
        new IndexEngineFactory("esa", text => EsaSuffixIndex.Build(text)),
        new IndexEngineFactory("esa-packed", text => EsaSuffixIndex.Build(text, ChildLayoutKind.Packed)),
        new IndexEngineFactory("esa-reference3", text => EsaSuffixIndex.Build(text, ChildLayoutKind.Reference3)),
        // the facade of work package 1b.6 (ISuffixTree + ISuffixTreeNavigation + ISuffixArrayAccess) as a contract engine, default and one other layout
        new IndexEngineFactory("esa-facade", text => EsaSuffixTree.Build(text, null)),
        new IndexEngineFactory("esa-facade-packed", text => EsaSuffixTree.Build(text, ChildLayoutKind.Packed)),
        // the container reader of work package 1b.7f: the tree reads a SESA file in place (Checked mode, Header validation); the file holds code units below 256 only
        new IndexEngineFactory("esa-file", text => SesaTestFile.FileTree(text), SesaTestFile.IsByteText),
        new IndexEngineFactory("esa-file-trusted", text => SesaTestFile.FileTree(text, psi: false, access: EsaAccessMode.Trusted, level: EsaValidationLevel.Structure), SesaTestFile.IsByteText),
    ];
}
