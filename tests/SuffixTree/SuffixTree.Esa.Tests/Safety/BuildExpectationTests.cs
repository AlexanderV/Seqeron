namespace SuffixTree.Esa.Tests.Safety;

/// <summary>
/// The guard of the checked build (work package S2.2, review finding R04): the CI job that is meant to run the checked build sets
/// <c>ESA_EXPECT_BOUNDS=true</c>, and a run whose <c>SuffixTree.Esa</c> assembly was NOT compiled with <c>ESA_BOUNDS</c> (a default build: the 19 tests that call
/// <c>Assert.Ignore</c> unless <see cref="EsaBuild.BoundsChecked"/> then skip, and a skipped test is green) must fail. The rule is in
/// <see cref="BuildExpectation"/>; this class pins every cell of it, the reading of the variable by its name, and (through
/// <see cref="BufAndArenaTests.BuildMode_IsReported_AndMatchesTheExpectedMode"/>) the live check of the run it is in.
///
/// The three honest configurations all pass: a default build without the variable (the default CI, dotnet.yml), a checked build with <c>true</c>
/// (the deep workflow), and a checked build WITHOUT the variable (a developer's <c>dotnet test -p:EsaBounds=true</c>: every test is stronger there,
/// nothing is skipped by the build mode, and the 19 tests report their own mode, so no expectation is needed). The dishonest one is the default build with <c>true</c>.
/// </summary>
[TestFixture]
public class BuildExpectationTests
{
    private const string Variable = "ESA_EXPECT_BOUNDS";

    /// <summary>(value of the variable, whether the assembly is checked, whether the verdict is a failure).</summary>
    public static IEnumerable<TestCaseData> Cells()
    {
        // not set: both builds are honest
        yield return new TestCaseData(null, false, false).SetName("Verdict_Unset_DefaultBuild_Holds");
        yield return new TestCaseData(null, true, false).SetName("Verdict_Unset_CheckedBuild_Holds");
        // true: only the checked build holds
        yield return new TestCaseData("true", true, false).SetName("Verdict_True_CheckedBuild_Holds");
        yield return new TestCaseData("true", false, true).SetName("Verdict_True_DefaultBuild_Fails");
        yield return new TestCaseData("True", true, false).SetName("Verdict_TrueCapitalised_CheckedBuild_Holds");
        yield return new TestCaseData("True", false, true).SetName("Verdict_TrueCapitalised_DefaultBuild_Fails");
        yield return new TestCaseData("TRUE", false, true).SetName("Verdict_TrueUpperCase_DefaultBuild_Fails");
        yield return new TestCaseData(" true ", false, true).SetName("Verdict_TrueWithBlanks_DefaultBuild_Fails");
        yield return new TestCaseData(" true ", true, false).SetName("Verdict_TrueWithBlanks_CheckedBuild_Holds");
        // false: only the default build holds
        yield return new TestCaseData("false", false, false).SetName("Verdict_False_DefaultBuild_Holds");
        yield return new TestCaseData("False", false, false).SetName("Verdict_FalseCapitalised_DefaultBuild_Holds");
        yield return new TestCaseData("FALSE", false, false).SetName("Verdict_FalseUpperCase_DefaultBuild_Holds");
        yield return new TestCaseData(" false ", false, false).SetName("Verdict_FalseWithBlanks_DefaultBuild_Holds");
        yield return new TestCaseData("false", true, true).SetName("Verdict_False_CheckedBuild_Fails");
        yield return new TestCaseData("FALSE", true, true).SetName("Verdict_FalseUpperCase_CheckedBuild_Fails");
        // anything else is a mistake of the workflow, and fails whatever the build is (a typo must not read as "not set" or as "false")
        foreach (string value in new[] { "", " ", "1", "0", "yes", "no", "tru", "truee", "t", "on", "checked", "true;" })
        {
            yield return new TestCaseData(value, false, true).SetName($"Verdict_Garbage_{Escape(value)}_DefaultBuild_Fails");
            yield return new TestCaseData(value, true, true).SetName($"Verdict_Garbage_{Escape(value)}_CheckedBuild_Fails");
        }
    }

    private static string Escape(string s) => s.Length == 0 ? "empty" : s.Replace(" ", "blank").Replace(";", "semicolon");

    [TestCaseSource(nameof(Cells))]
    public void Verdict_FollowsTheTableOfTheRule(string? value, bool boundsChecked, bool fails)
    {
        string? verdict = BuildExpectation.Verdict(value, boundsChecked);
        if (!fails)
        {
            Assert.That(verdict, Is.Null, "an honest configuration must pass");
            return;
        }
        Assert.That(verdict, Is.Not.Null.And.Not.Empty, "a dishonest configuration must fail");
        Assert.That(verdict, Does.Contain(Variable), "the message names the variable");
    }

    [Test]
    public void TheFailureOfADefaultBuild_SaysWhatToDo_AndTheFailureOfACheckedBuild_SaysWhatToUnset()
    {
        string unchecked_ = BuildExpectation.Verdict("true", boundsChecked: false)!;
        Assert.That(unchecked_, Does.Contain("-p:EsaBounds=true").And.Contain("ESA_BOUNDS"), "the fix of the default build is named");
        Assert.That(unchecked_, Does.Contain("skipped"), "the reason (a skipped test is green) is named");
        string checked_ = BuildExpectation.Verdict("false", boundsChecked: true)!;
        Assert.That(checked_, Does.Contain("false").And.Contain("ESA_BOUNDS"));
        string garbage = BuildExpectation.Verdict("1", boundsChecked: true)!;
        Assert.That(garbage, Does.Contain("'1'").And.Contain("true").And.Contain("false").And.Contain("unset"), "the accepted values are listed and the value is quoted");
    }

    [Test]
    public void TheVariableIsReadByItsName_AndTheVerdictUsesTheBuildOfTheEsaAssembly()
    {
        Assert.That(BuildExpectation.VariableName, Is.EqualTo(Variable));
        string? old = Environment.GetEnvironmentVariable(Variable);
        try
        {
            Environment.SetEnvironmentVariable(Variable, null);
            Assert.That(BuildExpectation.VerdictOfTheEnvironment(), Is.Null, "unset: no expectation");

            Environment.SetEnvironmentVariable(Variable, "true");
            Assert.That(BuildExpectation.VerdictOfTheEnvironment() is null, Is.EqualTo(EsaBuild.BoundsChecked), "true holds exactly in the checked build");
            Environment.SetEnvironmentVariable(Variable, "false");
            Assert.That(BuildExpectation.VerdictOfTheEnvironment() is null, Is.EqualTo(!EsaBuild.BoundsChecked), "false holds exactly in the default build");
            Environment.SetEnvironmentVariable(Variable, "maybe");
            Assert.That(BuildExpectation.VerdictOfTheEnvironment(), Is.Not.Null, "garbage fails in either build");
        }
        finally
        {
            Environment.SetEnvironmentVariable(Variable, old);
        }
    }

    [Test]
    public void TheLiveGuardOfTheRun_FailsExactlyWhenTheVariableContradictsTheBuild()
    {
        // BufAndArenaTests.BuildMode_IsReported_AndMatchesTheExpectedMode is the check that the CI job relies on: called here with the variable set in this process, so that the
        // glue between the variable, the build and the failure is pinned in a run that has no variable at all (the default CI)
        string? old = Environment.GetEnvironmentVariable(Variable);
        var live = new BufAndArenaTests();
        try
        {
            Environment.SetEnvironmentVariable(Variable, null);
            Assert.DoesNotThrow(live.BuildMode_IsReported_AndMatchesTheExpectedMode, "not set");

            Environment.SetEnvironmentVariable(Variable, "true");
            if (EsaBuild.BoundsChecked)
                Assert.DoesNotThrow(live.BuildMode_IsReported_AndMatchesTheExpectedMode, "true in the checked build");
            else
                Assert.That(Assert.Throws<AssertionException>(live.BuildMode_IsReported_AndMatchesTheExpectedMode)!.Message, Does.Contain("-p:EsaBounds=true"), "true in a default build");

            Environment.SetEnvironmentVariable(Variable, "false");
            if (EsaBuild.BoundsChecked)
                Assert.Throws<AssertionException>(live.BuildMode_IsReported_AndMatchesTheExpectedMode, "false in the checked build");
            else
                Assert.DoesNotThrow(live.BuildMode_IsReported_AndMatchesTheExpectedMode, "false in a default build");
        }
        finally
        {
            Environment.SetEnvironmentVariable(Variable, old);
        }
    }

    [Test]
    public void ABuildThatIsExpectedToBeCheckedButIsNot_IsTheOnlyConfigurationThatCannotPass()
    {
        // the whole space of the rule: 2 builds x {unset, true, false}; exactly the two crossed cells fail
        var failing = new List<string>();
        foreach (bool built in new[] { false, true })
            foreach (string? value in new[] { null, "true", "false" })
                if (BuildExpectation.Verdict(value, built) is not null) failing.Add($"{(built ? "checked" : "default")} build with {value ?? "no variable"}");
        Assert.That(failing, Is.EquivalentTo(new[] { "default build with true", "checked build with false" }));
    }
}
