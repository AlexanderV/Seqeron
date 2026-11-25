using System.Reflection;
using System.Text.RegularExpressions;
using System.Xml.Linq;

namespace Seqeron.Genomics.Tests.Architecture;

/// <summary>
/// Licence separation guard (B09 WP79, docs/Validation/review-2026-09/LICENSES.md). Seqeron's main assemblies are MIT-licensed; the ports of copyleft programs (GPL / LGPL)
/// live in separate assemblies named by their subject, each with its own LICENSE file and an SPDX header in every source file. The rules:
/// (1) no MIT project, MIT test project or MIT assembly references a copyleft assembly - by project reference, by InternalsVisibleTo, by assembly reference, by namespace or by a
/// string that names it (reflection); (2) every .cs file of a copyleft project (and of its test project) carries the SPDX header of the project's licence;
/// (3) every copyleft project has a LICENSE file that names the licence.
/// The copyleft projects are the ones whose csproj declares a <c>PackageLicenseExpression</c> that is not MIT.
/// </summary>
[TestFixture]
[Category("Architecture")]
public class LicenceSeparationTests
{
    private sealed record Project(string Name, string Dir, string CsprojPath, string? LicenseExpression, string[] ProjectReferences, string[] InternalsVisibleTo);

    private static readonly string Root = FindRepoRoot();

    private static readonly Project[] Projects = LoadProjects();

    private static bool IsCopyleft(Project p) => p.LicenseExpression is { } e && !e.StartsWith("MIT", StringComparison.Ordinal);

    private static IEnumerable<Project> CopyleftProjects => Projects.Where(IsCopyleft);

    /// <summary>A test project belongs to a copyleft assembly when its name is the assembly's name + ".Tests" (or it is the test project of a copyleft MCP server).</summary>
    private static bool IsCopyleftTestProject(Project p) => p.Name.EndsWith(".Tests", StringComparison.Ordinal)
        && CopyleftProjects.Any(c => c.Name + ".Tests" == p.Name);

    private static string FindRepoRoot()
    {
        string? dir = AppContext.BaseDirectory;
        while (dir != null && !File.Exists(Path.Combine(dir, "Seqeron.sln")))
            dir = Path.GetDirectoryName(dir);
        return dir ?? throw new DirectoryNotFoundException("Seqeron.sln not found above " + AppContext.BaseDirectory);
    }

    private static Project[] LoadProjects()
    {
        var list = new List<Project>();
        foreach (var path in Directory.EnumerateFiles(Path.Combine(Root, "src"), "*.csproj", SearchOption.AllDirectories)
                     .Concat(Directory.EnumerateFiles(Path.Combine(Root, "tests"), "*.csproj", SearchOption.AllDirectories)))
        {
            if (path.Contains($"{Path.DirectorySeparatorChar}obj{Path.DirectorySeparatorChar}", StringComparison.Ordinal))
                continue;
            var doc = XDocument.Load(path);
            string? expression = doc.Descendants("PackageLicenseExpression").FirstOrDefault()?.Value.Trim();
            var refs = doc.Descendants("ProjectReference").Select(e => Path.GetFileNameWithoutExtension((string?)e.Attribute("Include") ?? "")).ToArray();
            var ivt = doc.Descendants("InternalsVisibleTo").Select(e => (string?)e.Attribute("Include") ?? "").ToArray();
            list.Add(new Project(Path.GetFileNameWithoutExtension(path), Path.GetDirectoryName(path)!, path, expression, refs, ivt));
        }

        return [.. list];
    }

    private static IEnumerable<string> SourceFiles(Project p) => Directory.EnumerateFiles(p.Dir, "*.cs", SearchOption.AllDirectories)
        .Where(f => !f.Contains($"{Path.DirectorySeparatorChar}obj{Path.DirectorySeparatorChar}", StringComparison.Ordinal)
                    && !f.Contains($"{Path.DirectorySeparatorChar}bin{Path.DirectorySeparatorChar}", StringComparison.Ordinal));

    [Test]
    public void CopyleftAssemblies_Exist_AndTheirNamesDoNotAdvertiseTheLicence()
    {
        var names = CopyleftProjects.Select(p => p.Name).ToList();
        Assert.That(names, Is.Not.Empty, "the ports of the copyleft programs live in their own assemblies");
        Assert.That(names.Where(n => Regex.IsMatch(n, "gpl|copyleft", RegexOptions.IgnoreCase)), Is.Empty,
            "assemblies are named by what they hold; the licence is stated in LICENSE and in the file headers");
    }

    /// <summary>The copyleft assemblies of the ports (one per subject area, B09 WP79 / WP80) as LICENSES.md section 2 lists them: the set the owner decision "variant 1" asked for. Read as data - this guard never names an assembly literally.</summary>
    private static IEnumerable<string> ListedCopyleftAssemblies()
    {
        string text = File.ReadAllText(Path.Combine(Root, "docs", "Validation", "review-2026-09", "LICENSES.md"));
        int start = text.IndexOf("## 2. Copyleft assemblies", StringComparison.Ordinal);
        int end = text.IndexOf("## 3. ", start, StringComparison.Ordinal);
        foreach (Match row in Regex.Matches(text[start..end], @"(?m)^\| ((?:`Seqeron\.[\w.]+`(?:, )?)+)[^|]*\|"))
        {
            foreach (Match name in Regex.Matches(row.Groups[1].Value, @"`(Seqeron\.[\w.]+)`"))
                yield return name.Groups[1].Value;
        }
    }

    [Test]
    public void EveryExpectedCopyleftAssembly_ExistsWithItsOwnTestProject()
    {
        var copyleft = CopyleftProjects.Select(p => p.Name).ToHashSet();
        var all = Projects.Select(p => p.Name).ToHashSet();
        var violations = new List<string>();
        var listed = ListedCopyleftAssemblies().ToList();
        Assert.That(listed.Count, Is.GreaterThanOrEqualTo(20), "LICENSES.md section 2 lists the copyleft assemblies and servers");
        foreach (var name in listed)
        {
            if (!copyleft.Contains(name))
                violations.Add(name + ": no copyleft project of this name (PackageLicenseExpression missing or MIT)");
            if (!all.Contains(name + ".Tests"))
                violations.Add(name + ": no test project " + name + ".Tests");
        }

        Assert.That(violations, Is.Empty);
    }

    [Test]
    public void NoMitMcpServer_ExposesAToolOfACopyleftMcpServer()
    {
        // tool names stay identical across the separation, but each tool lives on exactly one server: the copyleft servers' tools (their docs/mcp/tools/<server> documents) are not declared by an MIT server.
        var copyleftTools = new Dictionary<string, string>(StringComparer.Ordinal);
        foreach (var p in CopyleftProjects.Where(p => p.Name.StartsWith("Seqeron.Mcp.", StringComparison.Ordinal)))
        {
            string docs = Path.Combine(Root, "docs", "mcp", "tools", p.Name["Seqeron.Mcp.".Length..].ToLowerInvariant());
            if (!Directory.Exists(docs))
                continue;
            foreach (var f in Directory.EnumerateFiles(docs, "*.md"))
                copyleftTools[Path.GetFileNameWithoutExtension(f)] = p.Name;
        }

        Assert.That(copyleftTools, Is.Not.Empty, "the copyleft MCP servers have tool documents");
        var violations = new List<string>();
        foreach (var p in Projects.Where(p => !IsCopyleft(p) && p.Name.StartsWith("Seqeron.Mcp.", StringComparison.Ordinal) && !p.Name.EndsWith(".Tests", StringComparison.Ordinal)))
        {
            foreach (var f in SourceFiles(p))
            {
                string text = File.ReadAllText(f);
                foreach (Match m in Regex.Matches(text, "McpServerTool\\(Name = \"(\\w+)\""))
                {
                    if (copyleftTools.TryGetValue(m.Groups[1].Value, out var owner))
                        violations.Add($"{p.Name}: tool '{m.Groups[1].Value}' is a tool of {owner}");
                }
            }
        }

        Assert.That(violations, Is.Empty);
    }

    [Test]
    public void NoMitProject_ReferencesACopyleftProject()
    {
        var copyleft = CopyleftProjects.Select(p => p.Name).ToHashSet();
        var violations = Projects.Where(p => !IsCopyleft(p) && !IsCopyleftTestProject(p) && !IsCopyleftMcpTestProject(p))
            .SelectMany(p => p.ProjectReferences.Where(copyleft.Contains).Select(r => $"{p.Name} -> {r}"))
            .ToList();
        Assert.That(violations, Is.Empty, "MIT projects must not depend on copyleft projects (the copyleft projects may depend on the MIT ones)");
    }

    /// <summary>The test project of a copyleft MCP server (Seqeron.Mcp.X.Tests tests Seqeron.Mcp.X).</summary>
    private static bool IsCopyleftMcpTestProject(Project p) => IsCopyleftTestProject(p);

    [Test]
    public void NoMitProject_NamesACopyleftAssembly_InInternalsVisibleTo()
    {
        var copyleft = CopyleftProjects.Select(p => p.Name).ToHashSet();
        var violations = Projects.Where(p => !IsCopyleft(p))
            .SelectMany(p => p.InternalsVisibleTo.Where(i => copyleft.Any(c => i == c || i.StartsWith(c + ".", StringComparison.Ordinal))).Select(i => $"{p.Name} -> {i}"))
            .ToList();
        Assert.That(violations, Is.Empty, "an InternalsVisibleTo of an MIT assembly must not name a copyleft assembly");
    }

    [Test]
    public void NoMitSourceFile_NamesACopyleftAssemblyOrNamespace()
    {
        var copyleft = CopyleftProjects.Select(p => p.Name).Where(n => !n.EndsWith(".Tests", StringComparison.Ordinal)).ToArray();
        var rx = new Regex(@"\b(" + string.Join("|", copyleft.Select(Regex.Escape)) + @")\b");
        var violations = new List<string>();
        foreach (var p in Projects.Where(p => !IsCopyleft(p) && !IsCopyleftTestProject(p)))
            foreach (var f in SourceFiles(p))
            {
                // this guard names the copyleft assemblies only as data (it discovers them from the csproj files), never literally
                string text = Regex.Replace(File.ReadAllText(f), @"/\*.*?\*/|//[^\n]*", "", RegexOptions.Singleline); // documentation may mention a program; code may not name an assembly
                if (rx.IsMatch(text))
                    violations.Add(Path.GetRelativePath(Root, f));
            }

        Assert.That(violations, Is.Empty, "an MIT source file (namespace, using, reflection string) must not name a copyleft assembly");
    }

    [Test]
    public void MitAssemblies_DoNotReferenceACopyleftAssembly()
    {
        var copyleft = CopyleftProjects.Select(p => p.Name).ToHashSet();
        var mit = new[]
        {
            typeof(DnaSequence).Assembly, typeof(ScoringMatrix).Assembly, typeof(DisorderPredictor).Assembly, typeof(SequenceAligner).Assembly, typeof(FastaParser).Assembly,
            typeof(VariantCaller).Assembly, typeof(CrisprDesigner).Assembly, typeof(PhylogeneticAnalyzer).Assembly, typeof(PopulationGeneticsAnalyzer).Assembly,
            typeof(GenomeAssemblyAnalyzer).Assembly, typeof(PanGenomeAnalyzer).Assembly, typeof(ImmuneAnalyzer).Assembly, typeof(ReportGenerator).Assembly,
        };
        var violations = mit.SelectMany(a => a.GetReferencedAssemblies().Where(r => copyleft.Contains(r.Name ?? "")).Select(r => $"{a.GetName().Name} -> {r.Name}")).ToList();
        Assert.That(violations, Is.Empty);
        // the MIT assemblies carry no InternalsVisibleTo attribute that names a copyleft assembly
        var ivt = mit.SelectMany(a => a.GetCustomAttributes<System.Runtime.CompilerServices.InternalsVisibleToAttribute>().Select(x => x.AssemblyName))
            .Where(n => copyleft.Any(c => n == c || n.StartsWith(c + ".", StringComparison.Ordinal))).ToList();
        Assert.That(ivt, Is.Empty);
    }

    [Test]
    public void EveryCsFileOfACopyleftProject_CarriesTheSpdxHeaderOfItsLicence()
    {
        var violations = new List<string>();
        foreach (var p in Projects.Where(p => IsCopyleft(p) || IsCopyleftTestProject(p)))
        {
            var owner = IsCopyleft(p) ? p : CopyleftProjects.First(c => c.Name + ".Tests" == p.Name);
            string[] allowed = LicensesOf(owner);
            foreach (var f in SourceFiles(p))
            {
                string firstLine = File.ReadLines(f).FirstOrDefault()?.TrimStart('﻿') ?? "";
                var m = Regex.Match(firstLine, @"^// SPDX-License-Identifier: (\S+)$");
                if (!m.Success || !allowed.Contains(m.Groups[1].Value))
                    violations.Add(Path.GetRelativePath(Root, f));
            }
        }

        Assert.That(violations, Is.Empty, "every .cs file of a copyleft assembly starts with `// SPDX-License-Identifier: <the licence of the assembly>`");
    }

    /// <summary>The identifiers a file of the project may carry: the project's own; a GPL project may also hold files of an LGPL-derived origin (and the reverse is not allowed).</summary>
    private static string[] LicensesOf(Project p) => p.LicenseExpression switch
    {
        "GPL-3.0-or-later" => ["GPL-3.0-or-later", "GPL-2.0-or-later", "LGPL-2.1-or-later", "LGPL-3.0-or-later"],
        "AGPL-3.0-only" => ["AGPL-3.0-only"],
        "GPL-2.0-or-later" => ["GPL-2.0-or-later", "LGPL-2.1-or-later"],
        var e => [e ?? ""],
    };

    [Test]
    public void EveryCopyleftProject_HasALicenseFileThatNamesItsLicence()
    {
        var violations = new List<string>();
        foreach (var p in CopyleftProjects)
        {
            string path = Path.Combine(p.Dir, "LICENSE");
            if (!File.Exists(path))
            {
                violations.Add(p.Name + ": LICENSE missing");
                continue;
            }

            string text = File.ReadAllText(path);
            if (!text.Contains(p.LicenseExpression!, StringComparison.Ordinal) || !Regex.IsMatch(text, "GNU (Lesser |Library )?General Public License"))
                violations.Add(p.Name + ": LICENSE does not state " + p.LicenseExpression);
        }

        Assert.That(violations, Is.Empty);
    }
}
