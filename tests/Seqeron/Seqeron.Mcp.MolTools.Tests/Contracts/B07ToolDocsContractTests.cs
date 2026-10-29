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

    // Source files whose tool-doc line anchors are pinned (B08 round 4, H29). The tool list is derived from the docs
    // directory - every <tool>.mcp.json whose docRef points into one of these files - so a new CRISPR / restriction /
    // primer / probe tool is guarded the day its doc is added; the old hardcoded list silently exempted all of them.
    // A batch adds the stems of the files it has audited. A failure here means line numbers drifted:
    //   python3 scripts/docs/sync_moltools_docref.py --fix && python3 scripts/skills/gen-catalog.py
    private static readonly string[] AuditedSourceStems = { "PrimerDesigner", "ProbeDesigner", "CrisprDesigner", "RestrictionAnalyzer" };

    // docRef "Seqeron.Genomics.MolTools/{File}.cs#L{n}"; the file name may carry partial-class dots
    // (RestrictionAnalyzer.PartialDigest.cs), which the former `\w+\.cs` rejected.
    private static readonly Regex DocRefPattern =
        new(@"^Seqeron\.Genomics\.MolTools/(?<file>[\w.]+\.cs)#L(?<line>\d+)\b", RegexOptions.Compiled);

    // Source files whose tools document their examples in machine-comparable form (a JSON subset of the real result).
    // The primer / probe / restriction docs describe some outputs in prose ("10 enzymes: AluI, RsaI, ..."); a batch adds a
    // stem here once its examples are executable.
    private static readonly string[] ExecutableExampleStems = { "CrisprDesigner" };

    public static IEnumerable<string> ExecutableExampleTools() => ToolsWithDocRefIn(ExecutableExampleStems);

    public static IEnumerable<string> AnchoredTools() => ToolsWithDocRefIn(AuditedSourceStems);

    private static IEnumerable<string> ToolsWithDocRefIn(string[] stems)
    {
        string docsDir = DocsDirectory();
        foreach (string path in Directory.GetFiles(docsDir, "*.mcp.json").OrderBy(p => p, StringComparer.Ordinal))
        {
            using var doc = JsonDocument.Parse(File.ReadAllText(path));
            var m = DocRefPattern.Match(doc.RootElement.GetProperty("docRef").GetString()!);
            if (m.Success && stems.Any(stem => m.Groups["file"].Value.StartsWith(stem, StringComparison.Ordinal)))
                yield return doc.RootElement.GetProperty("toolName").GetString()!;
        }
    }

    [Test]
    public void AnchoredTools_AreDerivedFromTheDocsDirectory_AndCoverEveryGuardedFamily()
    {
        var tools = AnchoredTools().ToList();

        Assert.That(tools, Is.SupersetOf(B07Tools), "the primer / probe tools stay guarded");
        Assert.That(tools, Does.Contain("find_pam_sites").And.Contain("crispr_system_info").And.Contain("restriction_digest"),
            "CRISPR and restriction tools are guarded now");
        Assert.That(tools, Does.Contain("partial_digest"), "a partial-class file name (RestrictionAnalyzer.PartialDigest.cs) parses");
    }

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

    [TestCaseSource(nameof(AnchoredTools))]
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

    // The docRef anchor of every guarded tool must point at a declaration line of the methodId's method, and the .md must
    // link the same anchor (B07 audit round 5, A5-1: anchors had rotted after code edits; B08 round 4, H29: 18 of 55 were
    // stale because this test only ever saw 19 primer / probe tools).
    [TestCaseSource(nameof(AnchoredTools))]
    public void DocRef_LineAnchor_PointsAtMethodDeclaration(string tool)
    {
        string docsDir = DocsDirectory();
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.Combine(docsDir, tool + ".mcp.json")));
        string methodId = doc.RootElement.GetProperty("methodId").GetString()!;
        string docRef = doc.RootElement.GetProperty("docRef").GetString()!;
        string method = methodId.Split('.')[^1];

        var m = DocRefPattern.Match(docRef);
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
            Is.True, $"{tool}: {file}#L{line} is '{target.Trim()}', not the declaration of {method} " +
                     "(python3 scripts/docs/sync_moltools_docref.py --fix)");

        // The md carries the same anchor: the first "Source: [File.cs#Ln](...File.cs#Ln)" link must equal the docRef; a
        // doc that states its sources in another shape must still link the docRef anchor somewhere.
        string md = File.ReadAllText(Path.Combine(docsDir, tool + ".md"));
        var src = Regex.Match(md, @"Source: \[(?<file>[\w.]+\.cs)#L(?<a>\d+)\]\([^)]*/(?<file2>[\w.]+\.cs)#L(?<b>\d+)\)");
        if (src.Success)
        {
            Assert.That((src.Groups["file"].Value, src.Groups["a"].Value, src.Groups["file2"].Value, src.Groups["b"].Value),
                Is.EqualTo((file, line.ToString(), file, line.ToString())), $"{tool}.md: Source anchor differs from docRef");
        }
        else
        {
            Assert.That(md, Does.Contain($"{file}#L{line})"), $"{tool}.md: no link to the docRef anchor {file}#L{line}");
        }

        // Secondary anchors in the same doc (e.g. get_enzyme's "Overhang:" / "Isoschizomers:" links) must also sit on a
        // declaration, whatever method it is.
        foreach (Match a in Regex.Matches(md, @"(?<f>[\w.]+\.cs)#L(?<n>\d+)\)"))
        {
            string[] lines = File.ReadAllLines(Path.Combine(repoRoot, "src", "Seqeron", "Algorithms", "Seqeron.Genomics.MolTools", a.Groups["f"].Value));
            int n = int.Parse(a.Groups["n"].Value);
            Assert.That(n, Is.InRange(1, lines.Length), $"{tool}.md: {a.Value} out of range");
            Assert.That(lines[n - 1].Contains("public static") || (n >= 2 && lines[n - 2].Contains("public static")), Is.True,
                $"{tool}.md: {a.Groups["f"].Value}#L{n} is '{lines[n - 1].Trim()}', not a public static declaration");
        }
    }

    // Every documented example is an executable claim: its input, run through the real wrapper, must reproduce the
    // documented output (B08 round 5 A5-2: the design_guide_rnas example said "score 100" while the tool returns 90, and
    // nothing ran the CRISPR examples). Outputs may be abridged: an object lists a subset of the properties, an array a
    // prefix of the elements; numbers compare to 1e-9, an "error" output expects an exception whose type name or message
    // matches.
    [TestCaseSource(nameof(ExecutableExampleTools))]
    public void McpJson_Examples_ReproduceTheTool(string tool)
    {
        var method = ToolMethod(tool);
        var parameters = method.GetParameters();
        var options = ModelContextProtocol.McpJsonUtilities.DefaultOptions;
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.Combine(DocsDirectory(), tool + ".mcp.json")));

        Assert.Multiple(() =>
        {
            foreach (var example in doc.RootElement.GetProperty("examples").EnumerateArray())
            {
                string name = $"{tool}/{example.GetProperty("name").GetString()}";
                var input = example.GetProperty("input");
                var output = example.GetProperty("output");
                object?[] args = parameters.Select(p =>
                    input.TryGetProperty(p.Name!, out var value) ? value.Deserialize(p.ParameterType, options)
                    : p.HasDefaultValue ? p.DefaultValue
                    : null).ToArray();

                if (output.TryGetProperty("error", out var error))
                {
                    var thrown = Assert.Catch(() =>
                    {
                        try { method.Invoke(null, args); }
                        catch (TargetInvocationException e) { System.Runtime.ExceptionServices.ExceptionDispatchInfo.Capture(e.InnerException!).Throw(); }
                    }, name);
                    string expected = error.GetString()!;
                    Assert.That(thrown!.GetType().Name == expected || thrown.Message.Contains(expected, StringComparison.OrdinalIgnoreCase),
                        Is.True, $"{name}: documented error '{expected}' vs {thrown.GetType().Name}: {thrown.Message}");
                    continue;
                }

                object? result;
                try { result = method.Invoke(null, args); }
                catch (TargetInvocationException e) { Assert.Fail($"{name}: threw {e.InnerException}"); continue; }
                var actual = JsonSerializer.SerializeToElement(result, result!.GetType(), options);
                string? difference = FirstDifference(output, actual, "$");
                Assert.That(difference, Is.Null, $"{name}: {difference}");
            }
        });
    }

    private static string? FirstDifference(JsonElement documented, JsonElement actual, string path)
    {
        switch (documented.ValueKind)
        {
            case JsonValueKind.Object:
                if (actual.ValueKind != JsonValueKind.Object)
                    return $"{path}: documented object, tool returned {actual.ValueKind}";
                foreach (var property in documented.EnumerateObject())
                {
                    if (!actual.TryGetProperty(property.Name, out var counterpart))
                        return $"{path}.{property.Name}: not in the tool's result";
                    if (FirstDifference(property.Value, counterpart, $"{path}.{property.Name}") is { } inner)
                        return inner;
                }
                return null;
            case JsonValueKind.Array:
                if (actual.ValueKind != JsonValueKind.Array)
                    return $"{path}: documented array, tool returned {actual.ValueKind}";
                if (documented.GetArrayLength() > actual.GetArrayLength())
                    return $"{path}: documented {documented.GetArrayLength()} element(s), tool returned {actual.GetArrayLength()}";
                var elements = actual.EnumerateArray().ToArray();
                int i = 0;
                foreach (var element in documented.EnumerateArray())
                {
                    if (FirstDifference(element, elements[i], $"{path}[{i}]") is { } inner)
                        return inner;
                    i++;
                }
                return null;
            case JsonValueKind.Number:
                return actual.ValueKind == JsonValueKind.Number && Math.Abs(documented.GetDouble() - actual.GetDouble()) <= 1e-9 * Math.Max(1, Math.Abs(actual.GetDouble()))
                    ? null : $"{path}: documented {documented.GetRawText()}, tool returned {actual.GetRawText()}";
            default:
                return documented.ValueKind == actual.ValueKind && documented.GetRawText() == actual.GetRawText()
                    ? null : $"{path}: documented {documented.GetRawText()}, tool returned {actual.GetRawText()}";
        }
    }

    [TestCaseSource(nameof(AnchoredTools))]
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
