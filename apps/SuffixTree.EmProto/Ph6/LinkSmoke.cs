using SuffixTree.Esa;
using SuffixTree.Esa.Bench.Links;

namespace SuffixTree.EmProto.Ph6;

/// <summary>One kernel x endpoint-source cell of the smoke test.</summary>
/// <param name="Mismatches">Links on which the run with the tables under test differs from the reference run (the flag, the step count, the target, or an exception).</param>
/// <param name="Threw">Links on which the kernel threw with the tables under test (counted in <paramref name="Mismatches"/> as well).</param>
/// <param name="TargetOnly">Links that differ in the target alone: same success flag, same step count (counted in <paramref name="Mismatches"/> as well).</param>
internal sealed record SmokeCell(LinkKernel Kernel, EndpointSource Source, long Nodes, long Mismatches, long GaveUp, long DefinitionMismatches, long Threw = 0, long TargetOnly = 0);

internal sealed record SmokeResult(long NodeCount, long SampleCount, SmokeCell[] Cells)
{
    internal long Mismatches => Cells.Sum(c => c.Mismatches);

    internal long DefinitionMismatches => Cells.Sum(c => c.DefinitionMismatches);

    internal bool Ok => Mismatches == 0 && DefinitionMismatches == 0;
}

/// <summary>
/// The suffix-link smoke test of E9: every kernel of prototype P1 with every endpoint source on a sample of internal nodes, once with the in-RAM ISA / Psi of <c>LinkTables</c> and once with the tables of the streaming code
/// (ISA from DIST, Psi from the BWT, Psi from the ISA); the results (success flag, target, step count) must be identical, and, where a definition is given, equal to it.
/// </summary>
internal static class LinkSmoke
{
    /// <param name="index">The in-RAM index (tagged child table) of the text.</param>
    /// <param name="isaRam">P1's ISA (n + 1 entries).</param>
    /// <param name="psiRam">P1's Psi (n + 1 entries, entry 0 = -1).</param>
    /// <param name="isaOther">The ISA under test.</param>
    /// <param name="psiOther">The Psi under test, entry 0 = -1 (<see cref="PsiFiles.LoadForP1"/>).</param>
    /// <param name="sample">Number of nodes: with <paramref name="p1Sampling"/> the <c>--sample</c> of the P1 benchmark (sample A of 2 x sample + 4000 drawn nodes), else the number of nodes drawn (all when fewer exist).</param>
    /// <param name="seed">Seed of the draw (P1 default 1).</param>
    /// <param name="budget">Step budget of a link (P1: 3,000,000); a link that exceeds it in the reference run is counted as given up and must give up in the other run too.</param>
    /// <param name="fanout">Fan-out of the range helper (P1 default 16).</param>
    /// <param name="definition">Optional definitional target of a node.</param>
    internal static SmokeResult Run(EsaIndex index, int[] isaRam, int[] psiRam, int[] isaOther, int[] psiOther, int sample, ulong seed, long budget, int fanout, bool p1Sampling, Func<LinkNode, LinkTarget>? definition)
    {
        ArgumentNullException.ThrowIfNull(index);
        var helper = new LcpRangeHelper(index.Lcp, fanout);
        var reference = new LinkContext(index, isaRam, psiRam, helper);
        var other = new LinkContext(index, isaOther, psiOther, helper);
        LinkNodeSet nodes = LinkNodeSet.Enumerate(index.Lcp);
        int[] sa = index.SuffixArray;
        int[] chosen = [];
        if (nodes.Count > 0 && p1Sampling)
        {
            int[] drawn = nodes.Sample((2 * sample) + 4000, seed);
            int half = Math.Min(sample, Math.Max(1, (drawn.Length - 4000) / 2));
            chosen = drawn[..half];
        }
        else if (nodes.Count > 0)
        {
            chosen = nodes.Sample(sample, seed);
        }

        LinkNode[] picked = [.. chosen.Select(i => nodes.Node(i, sa))];
        LinkTarget[]? expected = definition is null ? null : [.. picked.Select(definition)];
        var cells = new List<SmokeCell>();
        foreach (LinkKernel kernel in Enum.GetValues<LinkKernel>())
            foreach (EndpointSource source in Enum.GetValues<EndpointSource>())
            {
                long mismatches = 0, gaveUp = 0, defMismatches = 0, threwCount = 0, targetOnly = 0;
                for (int k = 0; k < picked.Length; k++)
                {
                    LinkNode v = picked[k];
                    bool okA = LinkKernels.TryLink(kernel, source, reference, v, budget, out LinkTarget a, out long stepsA);
                    bool threw = false;
                    LinkTarget b = default;
                    long stepsB = 0;
                    bool okB = false;
                    try
                    {
                        okB = LinkKernels.TryLink(kernel, source, other, v, budget, out b, out stepsB);
                    }
                    catch (Exception ex) when (ex is IndexOutOfRangeException or ArgumentException or InvalidOperationException)
                    {
                        threw = true; // wrong tables can send a kernel outside its arrays: that is a mismatch, not a crash of the test
                        threwCount++;
                    }
                    bool sameRun = okA == okB && stepsA == stepsB;
                    if (threw || !sameRun || (okA && a != b)) mismatches++;
                    if (!threw && sameRun && okA && a != b) targetOnly++;
                    if (!okA) gaveUp++;
                    else if (expected is not null && a != expected[k]) defMismatches++;
                    if (!threw && okB && expected is not null && b != expected[k] && !(okA && a != expected[k])) defMismatches++;
                }
                cells.Add(new SmokeCell(kernel, source, picked.Length, mismatches, gaveUp, defMismatches, threwCount, targetOnly));
            }
        return new SmokeResult(nodes.Count, picked.Length, [.. cells]);
    }
}
