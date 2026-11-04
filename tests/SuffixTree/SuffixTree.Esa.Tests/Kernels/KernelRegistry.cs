using SuffixTree.Esa.Tests.Fuzz;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>What a suffix-array engine returns: the suffix array, the recursion levels it went through and, when it can report them, the intermediate states.</summary>
internal sealed record SaisRun(int[] Sa, IReadOnlyList<SaisLevel> Levels, IReadOnlyList<SaisLevelState>? States);

/// <summary>What an LCP engine returns: the compact LCP, and the permuted LCP array when the engine exposes it.</summary>
internal sealed record LcpRun(CompactLcp Lcp, int[]? Plcp);

/// <summary>A registered suffix-array engine. <see cref="Run"/> gets the prepared text and the guarded arena the harness will verify afterwards.</summary>
internal sealed record SaisEngine(string Name, Func<PreparedText, NativeArena, SaisRun> Run, bool InHeavyTests = true);

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

    /// <summary>1o.1b: the SA-IS pointer kernel as it ships (<see cref="NoProbe"/>), run in a guarded arena with a poisoned suffix-array block.</summary>
    internal static readonly SaisEngine PointerSais = new("SaisKernel.Pointer", RunPointerSais);

    /// <summary>The same kernel with a recording probe: the intermediate states (bucket arrays, LMS order, names) are compared with the oracle. Not part of the heavy categories (copies every level).</summary>
    internal static readonly SaisEngine PointerSaisWithStates = new("SaisKernel.Pointer+states", RunPointerSaisWithStates, InHeavyTests: false);

    /// <summary>1o.1a: the pointer kernel for Phi, PLCP and the compact LCP (<see cref="LcpKernel"/>), run on guarded copies of its inputs.</summary>
    internal static readonly LcpEngine PointerLcp = new("LcpKernel.Pointer", RunPointerLcp);

    /// <summary>1o.3b: the suffix-array-order scan (<see cref="LcpKernel.BuildDirect{TSym}"/>), with the pointer kernel as its fallback; it has no PLCP array.</summary>
    internal static readonly LcpEngine DirectLcp = new("LcpKernel.Direct", (p, sa, arena) => RunDirectLcp(p, sa, arena, forceTextOrder: false));

    /// <summary>1o.3c: the same kernel with the per-rank resolution switched off, i.e. every long rank is resolved in text order (the path a long repeat takes); a scan that gives up still falls back to the three-pass kernel.</summary>
    internal static readonly LcpEngine DirectTextOrderLcp = new("LcpKernel.DirectTextOrder", (p, sa, arena) => RunDirectLcp(p, sa, arena, forceTextOrder: true));

    /// <summary>All suffix-array engines, the reference first.</summary>
    internal static IReadOnlyList<SaisEngine> SaisEngines { get; } = [ReferenceSais, PointerSais, PointerSaisWithStates];

    /// <summary>All LCP engines, the reference first.</summary>
    internal static IReadOnlyList<LcpEngine> LcpEngines { get; } = [ReferenceLcp, PointerLcp, DirectLcp, DirectTextOrderLcp];

    /// <summary>The suffix-array engines that are not the reference (the optimized kernels).</summary>
    internal static IEnumerable<SaisEngine> OptimizedSais => SaisEngines.Where(e => e != ReferenceSais && e.InHeavyTests);

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

    private static LcpRun RunPointerLcp(PreparedText prepared, int[] sa, NativeArena arena)
    {
        Buf<int> saBuf = arena.Allocate<int>(sa.Length, zero: false, "sa");
        sa.CopyTo(saBuf.AsSpan(sa.Length));
        return prepared.Width == SymbolWidth.Byte
            ? RunPointerLcp(prepared.Bytes, sa, saBuf, arena)
            : RunPointerLcp(prepared.UShorts, sa, saBuf, arena);
    }

    private static LcpRun RunDirectLcp(PreparedText prepared, int[] sa, NativeArena arena, bool forceTextOrder)
    {
        Buf<int> saBuf = arena.Allocate<int>(sa.Length, zero: false, "sa");
        sa.CopyTo(saBuf.AsSpan(sa.Length));
        return prepared.Width == SymbolWidth.Byte
            ? RunDirectLcp(prepared.Bytes, sa, saBuf, arena, forceTextOrder)
            : RunDirectLcp(prepared.UShorts, sa, saBuf, arena, forceTextOrder);
    }

    private static LcpRun RunDirectLcp<TSym>(ReadOnlySpan<TSym> symbols, int[] original, Buf<int> sa, NativeArena arena, bool forceTextOrder)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int n = symbols.Length;
        Buf<TSym> text = arena.Allocate<TSym>(n, zero: false, "text");
        symbols.CopyTo(text.AsSpan(n));
        CompactLcp lcp = LcpKernel.BuildDirect(text, n, sa, arena, report: forceTextOrder ? new DirectReport { ForceTextOrder = true } : null);
        if (!sa.AsSpan(n + 1).SequenceEqual(original) || !text.AsSpan(n).SequenceEqual(symbols))
            throw new InvalidOperationException("the LCP kernel modified its input");
        return new LcpRun(lcp, null);
    }

    private static LcpRun RunPointerLcp<TSym>(ReadOnlySpan<TSym> symbols, int[] original, Buf<int> sa, NativeArena arena)
        where TSym : unmanaged, IEquatable<TSym>
    {
        int n = symbols.Length;
        Buf<TSym> text = arena.Allocate<TSym>(n, zero: false, "text");
        symbols.CopyTo(text.AsSpan(n));
        CompactLcp lcp = LcpKernel.Build(text, n, sa, arena);
        Buf<int> plcp = LcpKernel.ComputePlcp(text, n, sa, arena, out _, out _);
        if (!sa.AsSpan(n + 1).SequenceEqual(original) || !text.AsSpan(n).SequenceEqual(symbols))
            throw new InvalidOperationException("the LCP kernel modified its input");
        return new LcpRun(lcp, plcp.AsSpan(n).ToArray());
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

    private static SaisRun RunPointerSais(PreparedText prepared, NativeArena arena) => RunPointerSais(prepared, arena, withStates: false);

    private static SaisRun RunPointerSaisWithStates(PreparedText prepared, NativeArena arena) => RunPointerSais(prepared, arena, withStates: true);

    private static SaisRun RunPointerSais(PreparedText prepared, NativeArena arena, bool withStates)
    {
        int n = prepared.Length;
        Buf<int> sa = arena.Allocate<int>(n + 1, zero: false, "sa"); // poisoned in a guarded arena: every entry must be written by the kernel
        var trace = new SaisTrace();
        var recorder = new SaisStateRecorder();
        if (prepared.Width == SymbolWidth.Byte)
        {
            Buf<byte> text = arena.Allocate<byte>(n, zero: false, "text");
            prepared.Bytes.CopyTo(text.AsSpan(n));
            if (withStates) SaisKernel.Build(text, n, prepared.AlphabetSize, sa, arena, new RecordingSaisProbe(recorder), trace);
            else SaisKernel.Build(text, n, prepared.AlphabetSize, sa, arena, default(NoProbe), trace);
            if (!text.AsSpan(n).SequenceEqual(prepared.Bytes)) throw new InvalidOperationException("the SA-IS kernel modified its text");
        }
        else
        {
            Buf<ushort> text = arena.Allocate<ushort>(n, zero: false, "text");
            prepared.UShorts.CopyTo(text.AsSpan(n));
            if (withStates) SaisKernel.Build(text, n, prepared.AlphabetSize, sa, arena, new RecordingSaisProbe(recorder), trace);
            else SaisKernel.Build(text, n, prepared.AlphabetSize, sa, arena, default(NoProbe), trace);
            if (!text.AsSpan(n).SequenceEqual(prepared.UShorts)) throw new InvalidOperationException("the SA-IS kernel modified its text");
        }
        return new SaisRun(sa.AsSpan(n + 1).ToArray(), trace.Levels, withStates && n > 0 ? recorder.ToStates() : null);
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
