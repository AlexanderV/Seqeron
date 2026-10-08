using System.Reflection;
using System.Text.RegularExpressions;
using System.Text.Json;
using ModelContextProtocol.Server;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Pins the machine-readable docs (<c>docs/mcp/tools/moltools/{tool}.mcp.json</c>) of every MolTools tool that wraps a
/// PRIMER-* / PROBE-* class (PrimerDesigner, ProbeDesigner, ntthal) to the wrapper itself (B07 audit round 3, A3-22):
/// <c>inputSchema</c> lists exactly the parameters of the schema the MCP SDK advertises for the method — same names in
/// the same order, same JSON types, descriptions and defaults, same required set — and <c>outputSchema</c> lists exactly
/// the public properties of the result type. A parameter added to a wrapper without updating its .mcp.json fails here.
/// </summary>
[TestFixture]
public sealed class B07ToolDocsContractTests
{
    private static readonly string[] B07Tools =
    {
        "design_primers", "evaluate_primer", "primer_melting_temperature", "primer_melting_temperature_salt",
        "longest_homopolymer", "longest_dinucleotide_repeat", "hairpin_potential", "primer_dimer", "three_prime_stability",
        "generate_primer_candidates", "design_probes", "design_tiling_probes", "design_antisense_probes",
        "design_molecular_beacon", "validate_probe", "design_probes_primer3", "analyze_oligo",
        "oligo_extinction_coefficient", "oligo_concentration_from_absorbance",
    };

    private static string DocsDirectory()
    {
        for (var dir = new DirectoryInfo(TestContext.CurrentContext.TestDirectory); dir is not null; dir = dir.Parent)
        {
            string candidate = Path.Combine(dir.FullName, "docs", "mcp", "tools", "moltools");
            if (Directory.Exists(candidate))
                return candidate;
        }
        throw new DirectoryNotFoundException("docs/mcp/tools/moltools not found above the test directory.");
    }

    private static MethodInfo ToolMethod(string name) =>
        typeof(MolToolsTools).GetMethods(BindingFlags.Public | BindingFlags.Static)
            .Single(m => m.GetCustomAttribute<McpServerToolAttribute>()?.Name == name);

    // JSON values equal up to formatting (numbers compared numerically).
    private static bool JsonEqual(JsonElement a, JsonElement b)
    {
        if (a.ValueKind == JsonValueKind.Number && b.ValueKind == JsonValueKind.Number)
            return a.GetDouble() == b.GetDouble();
        if (a.ValueKind != b.ValueKind)
            return false;
        return a.ValueKind switch
        {
            JsonValueKind.Array => a.GetArrayLength() == b.GetArrayLength()
                && a.EnumerateArray().Zip(b.EnumerateArray()).All(p => JsonEqual(p.First, p.Second)),
            JsonValueKind.String => a.GetString() == b.GetString(),
            JsonValueKind.Object => a.EnumerateObject().Count() == b.EnumerateObject().Count()
                && a.EnumerateObject().All(p => b.TryGetProperty(p.Name, out var q) && JsonEqual(p.Value, q)),
            _ => true,
        };
    }

    private static string[] Required(JsonElement schema) =>
        schema.TryGetProperty("required", out var r) ? r.EnumerateArray().Select(x => x.GetString()!).ToArray() : Array.Empty<string>();

    [TestCaseSource(nameof(B07Tools))]
    public void McpJson_InputSchema_MatchesWrapperParameters(string tool)
    {
        var method = ToolMethod(tool);
        var advertised = McpServerTool.Create(method, target: null, new McpServerToolCreateOptions { Name = tool })
            .ProtocolTool.InputSchema;
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.Combine(DocsDirectory(), tool + ".mcp.json")));
        var documented = doc.RootElement.GetProperty("inputSchema");

        var want = advertised.GetProperty("properties").EnumerateObject().ToList();
        var have = documented.GetProperty("properties").EnumerateObject().ToList();
        Assert.That(have.Select(p => p.Name), Is.EqualTo(want.Select(p => p.Name)), $"{tool}: parameter names / order");
        Assert.That(Required(documented), Is.EquivalentTo(Required(advertised)), $"{tool}: required parameters");
        Assert.Multiple(() =>
        {
            foreach (var (w, h) in want.Zip(have))
            {
                Assert.That(JsonEqual(h.Value.GetProperty("type"), w.Value.GetProperty("type")), Is.True, $"{tool}.{w.Name}: type");
                Assert.That(h.Value.GetProperty("description").GetString(), Is.EqualTo(w.Value.GetProperty("description").GetString()),
                    $"{tool}.{w.Name}: description");
                bool wd = w.Value.TryGetProperty("default", out var wDefault);
                bool hd = h.Value.TryGetProperty("default", out var hDefault);
                Assert.That(hd, Is.EqualTo(wd), $"{tool}.{w.Name}: default present");
                if (wd && hd)
                    Assert.That(JsonEqual(hDefault, wDefault), Is.True, $"{tool}.{w.Name}: default {hDefault} vs {wDefault}");
            }
        });
    }

    // docRef "Seqeron.Genomics.MolTools/{File}.cs#L{n}" and the .md "Source: [File.cs#L{n}](…)" must both point at a
    // declaration line of the methodId's method (B07 audit round 5, A5-1: anchors had rotted after code edits).
    [TestCaseSource(nameof(B07Tools))]
    public void DocRef_LineAnchor_PointsAtMethodDeclaration(string tool)
    {
        string docsDir = DocsDirectory();
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.Combine(docsDir, tool + ".mcp.json")));
        string methodId = doc.RootElement.GetProperty("methodId").GetString()!;
        string docRef = doc.RootElement.GetProperty("docRef").GetString()!;
        string method = methodId.Split('.')[^1];

        var m = Regex.Match(docRef, @"^Seqeron\.Genomics\.MolTools/(?<file>\w+\.cs)#L(?<line>\d+)\b");
        Assert.That(m.Success, Is.True, $"{tool}: docRef '{docRef}' has no File.cs#L<n> anchor");
        string file = m.Groups["file"].Value;
        int line = int.Parse(m.Groups["line"].Value);

        string repoRoot = Path.GetFullPath(Path.Combine(docsDir, "..", "..", "..", ".."));
        string[] source = File.ReadAllLines(Path.Combine(repoRoot, "src", "Seqeron", "Algorithms", "Seqeron.Genomics.MolTools", file));
        Assert.That(line, Is.InRange(1, source.Length), $"{tool}: {file}#L{line} out of range");
        string target = source[line - 1];
        string previous = line >= 2 ? source[line - 2] : string.Empty;
        Assert.That(Regex.IsMatch(target, $@"\b{Regex.Escape(method)}\(") &&
                    (target.Contains("public static") || previous.Contains("public static")),
            Is.True, $"{tool}: {file}#L{line} is '{target.Trim()}', not the declaration of {method}");

        string md = File.ReadAllText(Path.Combine(docsDir, tool + ".md"));
        var src = Regex.Match(md, @"Source: \[(?<file>\w+\.cs)#L(?<a>\d+)\]\([^)]*/(?<file2>\w+\.cs)#L(?<b>\d+)\)");
        Assert.That(src.Success, Is.True, $"{tool}.md: no anchored 'Source:' link");
        Assert.That((src.Groups["file"].Value, src.Groups["a"].Value, src.Groups["file2"].Value, src.Groups["b"].Value),
            Is.EqualTo((file, line.ToString(), file, line.ToString())), $"{tool}.md: Source anchor differs from docRef");
    }

    [TestCaseSource(nameof(B07Tools))]
    public void McpJson_OutputSchema_ListsResultProperties(string tool)
    {
        var resultType = ToolMethod(tool).ReturnType;
        var want = resultType.GetProperties(BindingFlags.Public | BindingFlags.Instance)
            .Select(p => JsonNamingPolicy.CamelCase.ConvertName(p.Name));
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.Combine(DocsDirectory(), tool + ".mcp.json")));
        var have = doc.RootElement.GetProperty("outputSchema").GetProperty("properties").EnumerateObject().Select(p => p.Name);
        Assert.That(have, Is.EquivalentTo(want), $"{tool}: outputSchema properties vs {resultType.Name}");
    }
}
