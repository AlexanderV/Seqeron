namespace SuffixTree.Esa.Tests.Safety;

/// <summary>
/// The STUB of the tests-first commit of work package S2.2: the rule does not exist yet (every configuration passes), so that
/// <see cref="BuildExpectationTests"/> is red against it. Replaced by the rule in the next commit.
/// </summary>
internal static class BuildExpectation
{
    internal const string VariableName = "ESA_EXPECT_BOUNDS";

    internal static string? Verdict(string? value, bool boundsChecked) => null;

    internal static string? VerdictOfTheEnvironment() => Verdict(Environment.GetEnvironmentVariable(VariableName), EsaBuild.BoundsChecked);
}
