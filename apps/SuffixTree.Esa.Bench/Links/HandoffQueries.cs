// STUB (tests first): the runner of the hand-off query check. Replaced by the implementation in the next commit.
namespace SuffixTree.Esa.Bench.Links;

internal sealed record HandoffQueryOptions
{
    internal required string Container { get; init; }

    internal required string Input { get; init; }

    internal ulong Seed { get; init; } = 20261122;

    internal IReadOnlyList<int> Lengths { get; init; } = [12, 20, 32];

    internal int Present { get; init; } = 300;

    internal int Absent { get; init; } = 200;

    internal int Cap { get; init; } = 1000;

    internal string Engine { get; init; } = "auto";
}

internal static class HandoffQueryRun
{
    internal static int Execute(HandoffQueryOptions options, Action<string, string, string> row) => throw new NotImplementedException();
}
