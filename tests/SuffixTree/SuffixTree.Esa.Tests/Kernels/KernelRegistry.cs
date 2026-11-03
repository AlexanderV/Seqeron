using SuffixTree.Esa.Tests.Fuzz;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>What a suffix-array engine returns: the suffix array, the recursion levels it went through and, when it can report them, the intermediate states.</summary>
internal sealed record SaisRun(int[] Sa, IReadOnlyList<SaisLevel> Levels, IReadOnlyList<SaisLevelState>? States);

/// <summary>What an LCP engine returns: the compact LCP, and the permuted LCP array when the engine exposes it.</summary>
internal sealed record LcpRun(CompactLcp Lcp, int[]? Plcp);

/// <summary>A registered suffix-array engine. <see cref="Run"/> gets the prepared text and the guarded arena the harness will verify afterwards.</summary>
internal sealed record SaisEngine(string Name, Func<PreparedText, NativeArena, SaisRun> Run);

/// <summary>A registered LCP engine; the suffix array is the (reference) input.</summary>
internal sealed record LcpEngine(string Name, Func<PreparedText, int[], NativeArena, LcpRun> Run);

/// <summary>
/// The kernel-vs-reference registry (safety infrastructure, work package 1o.0). Every optimized kernel is registered here once and is then run
/// through the same harness (<see cref="KernelHarness"/>) on every corpus: bit-identical results with <see cref="SaisReference"/> /
/// <see cref="LcpReference"/>, equal recursion trace, equal intermediate states against the definitional oracle, intact canaries after the call,
/// and the same result for several poison patterns in the uninitialised workspace. The reference engines are registered too, so the harness
/// itself is exercised against the permanent oracle.
/// </summary>
internal static class KernelRegistry
{
    internal static readonly SaisEngine ReferenceSais = new("SaisReference", RunReferenceSais);

    internal static readonly LcpEngine ReferenceLcp = new("LcpReference", RunReferenceLcp);

    /// <summary>All suffix-array engines, the reference first.</summary>
    internal static IReadOnlyList<SaisEngine> SaisEngines { get; } = [ReferenceSais];

    /// <summary>All LCP engines, the reference first.</summary>
    internal static IReadOnlyList<LcpEngine> LcpEngines { get; } = [ReferenceLcp];

    /// <summary>The suffix-array engines that are not the reference (the optimized kernels).</summary>
    internal static IEnumerable<SaisEngine> OptimizedSais => SaisEngines.Where(e => e != ReferenceSais);

    /// <summary>The LCP engines that are not the reference (the optimized kernels).</summary>
    internal static IEnumerable<LcpEngine> OptimizedLcp => LcpEngines.Where(e => e != ReferenceLcp);

    /// <summary>The reference SA with its trace, as the harness compares against.</summary>
    internal static SaisRun RunReferenceSais(PreparedText prepared, NativeArena arena)
    {
        var trace = new SaisTrace();
        int[] sa = prepared.Width == SymbolWidth.Byte
            ? SaisReference.Build(prepared.ByteReader(), prepared.AlphabetSize, trace)
            : SaisReference.Build(prepared.UShortReader(), prepared.AlphabetSize, trace);
        return new SaisRun(sa, trace.Levels, null);
    }

    /// <summary>The reference LCP (and PLCP).</summary>
    internal static LcpRun RunReferenceLcp(PreparedText prepared, int[] sa, NativeArena arena)
    {
        if (prepared.Width == SymbolWidth.Byte)
            return new LcpRun(LcpReference.Build(prepared.ByteReader(), sa), LcpReference.ComputePlcp(prepared.ByteReader(), sa));
        return new LcpRun(LcpReference.Build(prepared.UShortReader(), sa), LcpReference.ComputePlcp(prepared.UShortReader(), sa));
    }

    /// <summary>The pipelines of the optimized kernels, each combined with the reference of the other stage and (when both exist) with each other.</summary>
    internal static IEnumerable<(string Name, EsaBuilder Build)> OptimizedPipelines()
    {
        foreach (SaisEngine s in OptimizedSais) yield return ($"{s.Name}+{ReferenceLcp.Name}", Pipeline(s, ReferenceLcp));
        foreach (LcpEngine l in OptimizedLcp) yield return ($"{ReferenceSais.Name}+{l.Name}", Pipeline(ReferenceSais, l));
        foreach (SaisEngine s in OptimizedSais)
            foreach (LcpEngine l in OptimizedLcp)
                yield return ($"{s.Name}+{l.Name}", Pipeline(s, l));
    }

    /// <summary>An <see cref="EsaBuilder"/> (text to SA + LCP) from an engine pair, for the fuzz harness and the vector tests.</summary>
    internal static EsaBuilder Pipeline(SaisEngine sais, LcpEngine lcp) => text =>
    {
        PreparedText prepared = TextPreparation.Prepare(text);
        using var arena = new NativeArena(ArenaOptions.Guarded());
        SaisRun s = sais.Run(prepared, arena);
        arena.VerifyCanaries();
        LcpRun l = lcp.Run(prepared, s.Sa, arena);
        arena.VerifyCanaries();
        return (s.Sa, l.Lcp);
    };
}
