namespace Seqeron.Genomics.Tests.Architecture;

/// <summary>
/// Audit A2-44 (B09 F87): a doc-lint for the files of the alignment / genomic-analysis batch. One member must carry ONE <c>&lt;summary&gt;</c>: a contiguous
/// <c>///</c> block with two summaries (the stale text of a removed member left above the real one: <c>SelectCenterByPairwiseScore</c>) is invisible to the compiler
/// (CS1587 / CS1570 only fire for misplaced / malformed XML) and misleads IntelliSense. The scan is textual and needs no build option.
/// </summary>
[TestFixture]
public class XmlDocLintTests
{
    private static string RepoRoot()
    {
        var dir = new DirectoryInfo(TestContext.CurrentContext.TestDirectory);
        while (dir is not null && !File.Exists(Path.Combine(dir.FullName, "Seqeron.sln")))
            dir = dir.Parent;
        return dir?.FullName ?? throw new InvalidOperationException("Seqeron.sln not found above the test directory.");
    }

    private static IEnumerable<string> BatchFiles()
    {
        string algorithms = Path.Combine(RepoRoot(), "src", "Seqeron", "Algorithms");
        foreach (var f in Directory.EnumerateFiles(Path.Combine(algorithms, "Seqeron.Genomics.Alignment"), "SequenceAligner*.cs"))
            yield return f;
        yield return Path.Combine(algorithms, "Seqeron.Genomics.Alignment", "AnchorBasedAligner.cs");
        yield return Path.Combine(algorithms, "Seqeron.Genomics.Analysis", "GenomicAnalyzer.cs");
        yield return Path.Combine(algorithms, "Seqeron.Genomics.Infrastructure", "AlignmentTypes.cs");
        yield return Path.Combine(algorithms, "Seqeron.Genomics.Core", "OrfScanner.cs");
    }

    private static List<string> Scan(string name, string[] lines)
    {
        var violations = new List<string>();
        int summaries = 0, blockStart = 0;
        for (int i = 0; i <= lines.Length; i++)
        {
            bool isDoc = i < lines.Length && lines[i].TrimStart().StartsWith("///", StringComparison.Ordinal);
            if (isDoc)
            {
                if (i == 0 || !lines[i - 1].TrimStart().StartsWith("///", StringComparison.Ordinal))
                    blockStart = i + 1;
                if (lines[i].Contains("<summary>", StringComparison.Ordinal) || lines[i].Contains("<summary ", StringComparison.Ordinal))
                    summaries++;
                continue;
            }

            if (summaries > 1)
                violations.Add($"{name}:{blockStart}: {summaries} <summary> elements in one doc block");
            summaries = 0;
        }

        return violations;
    }

    [Test]
    public void EveryDocBlock_HasAtMostOneSummary()
    {
        var violations = BatchFiles().SelectMany(f => Scan(Path.GetFileName(f), File.ReadAllLines(f))).ToList();
        Assert.That(violations, Is.Empty);
    }

    [Test]
    public void TheLint_FindsTheOriginalDefect_AndAcceptsOneSummary()
    {
        // the pattern of SelectCenterByPairwiseScore before the fix
        string[] defect =
        [
            "    /// <summary>", "    /// Gusfield (1993) center ...", "    /// </summary>",
            "    /// <summary>Gusfield center for an arbitrary optimal pair score ...</summary>", "    private static int SelectCenter() => 0;",
        ];
        Assert.That(Scan("x.cs", defect), Has.Count.EqualTo(1));
        string[] fine =
        [
            "    /// <summary>One.</summary>", "    int A();", "", "    /// <summary>Two.</summary>", "    int B();",
        ];
        Assert.That(Scan("x.cs", fine), Is.Empty);
    }
}
