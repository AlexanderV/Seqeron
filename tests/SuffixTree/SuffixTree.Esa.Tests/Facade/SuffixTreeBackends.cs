using SuffixTree.Esa;
using SuffixTree.Persistent;
using UkkonenTree = global::SuffixTree.SuffixTree;

namespace SuffixTree.Esa.Tests.Facade;

/// <summary>
/// One way to build an <see cref="ISuffixTree"/> over a text: the seam through which the existing ISuffixTree-level assertions are re-run for every backend
/// (work package 1b.6). The cross-backend selection lives here, in the tests, because a production factory that offers Ukkonen, v6 and the ESA together would
/// have to reference all three projects (the ESA factory of the product is <see cref="EsaSuffixTreeFactory"/>).
/// </summary>
/// <param name="Name">Display name, also the NUnit fixture name.</param>
/// <param name="FromString">Builds the tree of a string.</param>
/// <param name="FromSource">Builds the tree of a text source.</param>
/// <param name="IsOracle">True for the backends that are the oracles of the ESA (the in-memory Ukkonen tree and the persistent v6 tree).</param>
public sealed record SuffixTreeBackend(string Name, Func<string, ISuffixTree> FromString, Func<ITextSource, ISuffixTree> FromSource, bool IsOracle)
{
    /// <inheritdoc />
    public override string ToString() => Name;
}

/// <summary>The registered backends.</summary>
public static class SuffixTreeBackends
{
    /// <summary>The in-memory Ukkonen tree (oracle L2).</summary>
    public static readonly SuffixTreeBackend Ukkonen = new("ukkonen", UkkonenTree.Build, UkkonenTree.Build, true);

    /// <summary>The persistent v6 tree in heap mode (oracle L3).</summary>
    public static readonly SuffixTreeBackend PersistentHeap = new(
        "persistent-heap",
        text => PersistentSuffixTreeFactory.Create(new StringTextSource(text)),
        source => PersistentSuffixTreeFactory.Create(source),
        true);

    /// <summary>The ESA facade with the default child-table layout.</summary>
    public static readonly SuffixTreeBackend Esa = new("esa", EsaSuffixTreeFactory.Create, EsaSuffixTreeFactory.Create, false);

    /// <summary>The ESA facade with the LCP-disambiguated single-array child table.</summary>
    public static readonly SuffixTreeBackend EsaPacked = new("esa-packed", text => EsaSuffixTree.Build(text, ChildLayoutKind.Packed), source => EsaSuffixTree.Build(source.ToString()!, ChildLayoutKind.Packed), false);

    /// <summary>The ESA facade with the three-array reference child table.</summary>
    public static readonly SuffixTreeBackend EsaReference3 = new("esa-reference3", text => EsaSuffixTree.Build(text, ChildLayoutKind.Reference3), source => EsaSuffixTree.Build(source.ToString()!, ChildLayoutKind.Reference3), false);

    /// <summary>Every backend.</summary>
    public static IReadOnlyList<SuffixTreeBackend> All { get; } = [Ukkonen, PersistentHeap, Esa, EsaPacked, EsaReference3];

    /// <summary>The ESA backends (every child-table layout).</summary>
    public static IReadOnlyList<SuffixTreeBackend> EsaBackends { get; } = [Esa, EsaPacked, EsaReference3];

    /// <summary>The oracle backends.</summary>
    public static IReadOnlyList<SuffixTreeBackend> Oracles { get; } = [Ukkonen, PersistentHeap];

    /// <summary>The backends that run the derived copies of the existing tests: both oracles (they prove each ported assertion is backend-agnostic) and the default ESA.</summary>
    public static IEnumerable<TestFixtureData> PortedFixtures()
    {
        foreach (SuffixTreeBackend backend in new[] { Ukkonen, PersistentHeap, Esa })
            yield return new TestFixtureData(backend).SetArgDisplayNames(backend.Name);
    }

    /// <summary>All five backends, for fixtures that are cheap.</summary>
    public static IEnumerable<TestFixtureData> AllFixtures()
    {
        foreach (SuffixTreeBackend backend in All)
            yield return new TestFixtureData(backend).SetArgDisplayNames(backend.Name);
    }

    /// <summary>The ESA layouts only.</summary>
    public static IEnumerable<TestFixtureData> EsaFixtures()
    {
        foreach (SuffixTreeBackend backend in EsaBackends)
            yield return new TestFixtureData(backend).SetArgDisplayNames(backend.Name);
    }
}
