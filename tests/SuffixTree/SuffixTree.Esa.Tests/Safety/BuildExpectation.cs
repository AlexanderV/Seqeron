namespace SuffixTree.Esa.Tests.Safety;

/// <summary>
/// The rule of the guard of the checked build (work package S2.2, review finding R04). A CI job that is meant to run the checked build
/// (<c>dotnet test ... -p:EsaBounds=true</c>) sets <c>ESA_EXPECT_BOUNDS=true</c>; this rule fails the run when the <c>SuffixTree.Esa</c> assembly it tests was NOT compiled with
/// <c>ESA_BOUNDS</c> (<see cref="EsaBuild.BoundsChecked"/> is false): 19 tests call <c>Assert.Ignore</c> in that build, and a skipped test is green.
/// <list type="bullet">
/// <item>Not set: no expectation. Both builds are honest (the default CI runs the default build without the variable; a developer's own <c>-p:EsaBounds=true</c> run is checked
/// without it, and is the stronger of the two: nothing is skipped by the build mode).</item>
/// <item><c>true</c> / <c>false</c> (any case, blanks around are ignored): the build must be the one named. <c>false</c> in a checked build fails as well: a variable that
/// states the expectation states it both ways.</item>
/// <item>Anything else, the empty string included (a workflow whose expression expanded to nothing): a failure in either build, never read as "not set" or as "false".</item>
/// </list>
/// Test code, not production code: the variable and the rule belong to the test run.
/// </summary>
internal static class BuildExpectation
{
    internal const string VariableName = "ESA_EXPECT_BOUNDS";

    /// <summary>Null when the configuration is honest; otherwise the message of the failure.</summary>
    /// <param name="value">The value of the variable, null when it is not set.</param>
    /// <param name="boundsChecked">Whether the assembly under test was compiled with <c>ESA_BOUNDS</c>.</param>
    internal static string? Verdict(string? value, bool boundsChecked)
    {
        if (value is null) return null;
        bool expectChecked;
        string trimmed = value.Trim();
        if (string.Equals(trimmed, "true", StringComparison.OrdinalIgnoreCase)) expectChecked = true;
        else if (string.Equals(trimmed, "false", StringComparison.OrdinalIgnoreCase)) expectChecked = false;
        else return $"{VariableName} is '{value}': it must be unset, 'true' or 'false'";

        if (expectChecked == boundsChecked) return null;
        return expectChecked
            ? $"{VariableName}=true but SuffixTree.Esa was compiled WITHOUT ESA_BOUNDS: this is a default build, in which the 19 tests that need the checked accessor are skipped (and a skipped test is green). " +
              "Build with -p:EsaBounds=true (dotnet build ... -p:EsaBounds=true, then dotnet test --no-build), or unset the variable."
            : $"{VariableName}=false but SuffixTree.Esa was compiled WITH ESA_BOUNDS: this is the checked build. Unset the variable, or build without -p:EsaBounds=true.";
    }

    /// <summary>The verdict of the process: the variable as it is set now, the build of the assembly under test.</summary>
    internal static string? VerdictOfTheEnvironment() => Verdict(Environment.GetEnvironmentVariable(VariableName), EsaBuild.BoundsChecked);
}
