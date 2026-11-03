using SuffixTree.Esa.Tests.TestSupport;

namespace SuffixTree.Esa.Tests.Kernels;

/// <summary>
/// Runs a registered engine against <see cref="SaisReference"/> / <see cref="LcpReference"/> and the independent oracles. Every check returns
/// <c>null</c> when all is equal and otherwise a one-line description of the FIRST difference, so that the self-tests can assert which class of
/// defect is reported. The workspace of the engine is a guarded <see cref="NativeArena"/> whose canaries are verified after every call, and
/// the call is repeated with different poison patterns in the uninitialised memory (a kernel that reads a cell before it writes it gives
/// different results).
/// </summary>
internal static class KernelHarness
{
    /// <summary>Poison bytes of the uninitialised workspace; the first is the one used for large texts.</summary>
    internal static readonly byte[] PoisonPatterns = [0xCD, 0x00, 0xFF, 0x5A];

    /// <summary>The independent oracle limit: naive SA up to this length, the O(n) verifier above.</summary>
    private const int NaiveLimit = 2000;

    internal static void AssertSais(SaisEngine engine, string text, bool allPoisons = true)
        => Assert.That(CheckSais(engine, text, allPoisons), Is.Null, $"{engine.Name} on {Describe(text)}");

    internal static void AssertLcp(LcpEngine engine, string text, bool allPoisons = true)
        => Assert.That(CheckLcp(engine, text, allPoisons), Is.Null, $"{engine.Name} on {Describe(text)}");

    /// <summary>SA equality with the reference and the independent oracle, equal recursion trace, equal intermediate states, intact canaries, poison independence.</summary>
    internal static string? CheckSais(SaisEngine engine, string text, bool allPoisons = true)
    {
        PreparedText prepared = TextPreparation.Prepare(text);
        using var refArena = new NativeArena();
        SaisRun reference = KernelRegistry.RunReferenceSais(prepared, refArena);
        List<SaisLevelState>? oracle = null;
        foreach (byte poison in allPoisons ? PoisonPatterns : PoisonPatterns.Take(1))
        {
            using var arena = new NativeArena(ArenaOptions.Guarded(poison));
            SaisRun run;
            try
            {
                run = engine.Run(prepared, arena);
                arena.VerifyCanaries();
            }
            catch (EsaMemoryViolationException ex)
            {
                return $"memory violation ({ex.Kind}): {ex.Message}";
            }
            catch (Exception ex)
            {
                return $"exception {ex.GetType().Name}: {ex.Message}";
            }

            string? d = FirstDifference("SA", run.Sa, reference.Sa);
            if (d is not null) return $"poison 0x{poison:X2}: {d}";
            if (text.Length <= NaiveLimit)
            {
                d = FirstDifference("SA vs naive", run.Sa, NaiveSuffixArray.Build(text));
                if (d is not null) return d;
            }
            else
            {
                string? verdict = SuffixArrayVerifier.Check(text, run.Sa);
                if (verdict is not null) return "SA fails the independent verifier: " + verdict;
            }

            if (run.Levels.Count != reference.Levels.Count)
                return $"trace: {run.Levels.Count} recursion levels, reference {reference.Levels.Count}";
            for (int i = 0; i < run.Levels.Count; i++)
                if (run.Levels[i] != reference.Levels[i])
                    return $"trace level {i}: {run.Levels[i]}, reference {reference.Levels[i]}";

            if (run.States is not null)
            {
                oracle ??= SaisStateOracle.Compute(Symbols(prepared), prepared.AlphabetSize);
                d = CompareStates(run.States, oracle);
                if (d is not null) return d;
            }
        }
        return null;
    }

    /// <summary>LCP (kind, escape count, every entry) and PLCP equality with the reference and the independent naive oracle, intact canaries, poison independence.</summary>
    internal static string? CheckLcp(LcpEngine engine, string text, bool allPoisons = true)
    {
        PreparedText prepared = TextPreparation.Prepare(text);
        int[] sa = SaisReference.BuildSuffixArray(text);
        using var refArena = new NativeArena();
        LcpRun reference = KernelRegistry.RunReferenceLcp(prepared, sa, refArena);
        int[] expected = reference.Lcp.ToArray();
        foreach (byte poison in allPoisons ? PoisonPatterns : PoisonPatterns.Take(1))
        {
            using var arena = new NativeArena(ArenaOptions.Guarded(poison));
            LcpRun run;
            try
            {
                run = engine.Run(prepared, sa, arena);
                arena.VerifyCanaries();
            }
            catch (EsaMemoryViolationException ex)
            {
                return $"memory violation ({ex.Kind}): {ex.Message}";
            }
            catch (Exception ex)
            {
                return $"exception {ex.GetType().Name}: {ex.Message}";
            }

            if (run.Lcp.Kind != reference.Lcp.Kind) return $"LCP kind {run.Lcp.Kind}, reference {reference.Lcp.Kind}";
            if (run.Lcp.OverflowCount != reference.Lcp.OverflowCount) return $"LCP escapes {run.Lcp.OverflowCount}, reference {reference.Lcp.OverflowCount}";
            if (run.Lcp.Length != reference.Lcp.Length) return $"LCP length {run.Lcp.Length}, reference {reference.Lcp.Length}";
            string? d = FirstDifference("LCP", run.Lcp.ToArray(), expected);
            if (d is not null) return $"poison 0x{poison:X2}: {d}";
            if (run.Plcp is not null && reference.Plcp is not null)
            {
                d = FirstDifference("PLCP", run.Plcp, reference.Plcp);
                if (d is not null) return d;
            }
            if (text.Length <= NaiveLimit)
            {
                d = FirstDifference("LCP vs naive", run.Lcp.ToArray(), NaiveLcp.Build(text, sa));
                if (d is not null) return d;
            }
            else
            {
                d = FirstDifference("LCP vs Kasai", run.Lcp.ToArray(), KasaiLcp.Build(text, sa));
                if (d is not null) return d;
            }
        }
        return null;
    }

    /// <summary>The first difference between two per-level state lists (bucket arrays, LMS positions, LMS order, names, reduced string), or null.</summary>
    internal static string? CompareStates(IReadOnlyList<SaisLevelState> actual, IReadOnlyList<SaisLevelState> expected)
    {
        if (actual.Count != expected.Count) return $"states: {actual.Count} levels, oracle {expected.Count}";
        for (int i = 0; i < actual.Count; i++)
        {
            SaisLevelState a = actual[i], e = expected[i];
            if (a.Length != e.Length || a.Sigma != e.Sigma) return $"level {i}: (n, sigma) = ({a.Length}, {a.Sigma}), oracle ({e.Length}, {e.Sigma})";
            string? d = FirstDifference($"level {i} symbol counts", a.Counts, e.Counts)
                ?? FirstDifference($"level {i} bucket heads", a.Heads, e.Heads)
                ?? FirstDifference($"level {i} bucket tails", a.Tails, e.Tails)
                ?? FirstDifference($"level {i} LMS positions", a.LmsPositions, e.LmsPositions);
            if (d is not null) return d;
            if (a.NameCount != e.NameCount) return $"level {i}: {a.NameCount} names, oracle {e.NameCount}";
            d = FirstDifference($"level {i} reduced string (names in text order)", a.Reduced, e.Reduced);
            if (d is not null) return d;

            // LMS order: a permutation of the LMS positions whose names are non-decreasing and equal to the oracle's sorted names
            int[] nameOf = SaisStateOracle.NameOfPosition(e);
            if (a.LmsSorted.Length != e.LmsSorted.Length) return $"level {i}: LMS order has {a.LmsSorted.Length} entries, oracle {e.LmsSorted.Length}";
            if (!a.LmsSorted.Order().SequenceEqual(e.LmsPositions)) return $"level {i}: LMS order is not a permutation of the LMS positions";
            for (int k = 0; k < a.LmsSorted.Length; k++)
            {
                int got = nameOf[a.LmsSorted[k]], want = nameOf[e.LmsSorted[k]];
                if (got != want) return $"level {i}: LMS order differs at {k}: position {a.LmsSorted[k]} has name {got}, oracle order has name {want}";
            }
        }
        return null;
    }

    /// <summary>The symbols of a prepared text as a dense int array (for the oracle).</summary>
    internal static int[] Symbols(PreparedText prepared)
    {
        var symbols = new int[prepared.Length];
        for (int i = 0; i < symbols.Length; i++) symbols[i] = prepared.SymbolAt(i);
        return symbols;
    }

    private static string? FirstDifference(string what, int[] actual, int[] expected)
    {
        int k = Math.Min(actual.Length, expected.Length);
        for (int i = 0; i < k; i++)
            if (actual[i] != expected[i]) return $"{what}[{i}] = {actual[i]}, expected {expected[i]}";
        return actual.Length == expected.Length ? null : $"{what} has {actual.Length} entries, expected {expected.Length}";
    }

    private static string Describe(string text)
        => text.Length <= 40 ? $"'{text.Replace("\0", "\\0", StringComparison.Ordinal)}' (n = {text.Length})" : $"n = {text.Length}, starts '{text[..24]}'";
}
