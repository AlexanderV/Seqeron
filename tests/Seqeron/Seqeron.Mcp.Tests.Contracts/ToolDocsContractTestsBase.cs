using System.ComponentModel;
using DescriptionAttribute = System.ComponentModel.DescriptionAttribute;
using System.Reflection;
using System.Text.Json;
using NUnit.Framework;

namespace Seqeron.Mcp.Tests.Contracts;

/// <summary>
/// Shared docs contract (audit A5-5, B09 F170): the documents of a server's tools (<c>docs/mcp/tools/&lt;server&gt;/&lt;tool&gt;.mcp.json</c> and <c>.md</c>) cannot drift from the code.
/// For every <c>[McpServerTool]</c> read by REFLECTION the test compares the tool description, the parameter names (and their order), type, default and <c>[Description]</c> with the
/// <c>description</c> / <c>inputSchema</c> of the JSON and the parameter names / Required cell of the <c>## Input Schema</c> table of the markdown document. The fix is always
/// <c>python3 scripts/docs/sync_mcp_tool_docs.py</c> (it regenerates the same fields from the same reflection). A derived fixture names the tool class, the server directory and
/// (for a server shared with other batches) the <c>methodId</c> prefixes of the tools it owns.
/// </summary>
public abstract class ToolDocsContractTestsBase
{
    protected abstract Type ToolsType { get; }

    /// <summary>The directory of the server under docs/mcp/tools, e.g. "alignment".</summary>
    protected abstract string DocsServer { get; }

    /// <summary>Only the tools whose document's <c>methodId</c> starts with one of these are checked (empty = every tool of the server).</summary>
    protected virtual string[] OwnedMethodIdPrefixes => [];

    [Test]
    public void ToolDocuments_MatchTheToolMetadataOfTheCode()
    {
        string dir = Path.Combine(FindRepoRoot(), "docs", "mcp", "tools", DocsServer);
        var problems = new List<string>();
        int checkedTools = 0;
        foreach (var method in ToolsType.GetMethods(BindingFlags.Public | BindingFlags.Static | BindingFlags.Instance | BindingFlags.DeclaredOnly))
        {
            var attribute = method.GetCustomAttributes(inherit: false).FirstOrDefault(a => a.GetType().Name == "McpServerToolAttribute");
            if (attribute is null) continue;
            string? name = attribute.GetType().GetProperty("Name")?.GetValue(attribute) as string;
            if (string.IsNullOrWhiteSpace(name)) continue;
            string jsonPath = Path.Combine(dir, name + ".mcp.json");
            if (!File.Exists(jsonPath))
            {
                problems.Add($"{name}: no document {name}.mcp.json");
                continue;
            }

            using var doc = JsonDocument.Parse(File.ReadAllText(jsonPath));
            var root = doc.RootElement;
            string methodId = root.TryGetProperty("methodId", out var mid) ? mid.GetString() ?? "" : "";
            if (OwnedMethodIdPrefixes.Length > 0 && !OwnedMethodIdPrefixes.Any(p => methodId.StartsWith(p, StringComparison.Ordinal)))
                continue;
            checkedTools++;
            CompareJson(name, method, root, problems);
            CompareMarkdown(name, method, Path.Combine(dir, name + ".md"), problems);
        }

        Assert.That(checkedTools, Is.GreaterThan(0), "no tool of the server was checked");
        Assert.That(problems, Is.Empty, $"{problems.Count} difference(s) between the tool documents and the code ({checkedTools} tools checked); run `python3 scripts/docs/sync_mcp_tool_docs.py`:\n" + string.Join("\n", problems.Take(25)));
    }

    private static IEnumerable<ParameterInfo> InputParameters(MethodInfo m) =>
        m.GetParameters().Where(p => p.ParameterType != typeof(CancellationToken) && p.ParameterType.Namespace?.StartsWith("ModelContextProtocol", StringComparison.Ordinal) != true);

    private static string JsonType(Type t)
    {
        t = Nullable.GetUnderlyingType(t) ?? t;
        if (t == typeof(string) || t == typeof(char) || t.IsEnum) return "string";
        if (t == typeof(bool)) return "boolean";
        if (t == typeof(int) || t == typeof(long) || t == typeof(short) || t == typeof(byte) || t == typeof(sbyte) || t == typeof(uint) || t == typeof(ulong) || t == typeof(ushort)) return "integer";
        if (t == typeof(double) || t == typeof(float) || t == typeof(decimal)) return "number";
        if (t.IsArray) return "array";
        if (t.IsGenericType && t.GetGenericTypeDefinition() is var g && (g == typeof(List<>) || g == typeof(IEnumerable<>) || g == typeof(IReadOnlyList<>) || g == typeof(IList<>) || g == typeof(IReadOnlyCollection<>) || g == typeof(ICollection<>)))
            return "array";
        return "object";
    }

    private static void CompareJson(string name, MethodInfo method, JsonElement root, List<string> problems)
    {
        string? description = method.GetCustomAttribute<DescriptionAttribute>()?.Description;
        if (description is not null && (!root.TryGetProperty("description", out var d) || d.GetString() != description))
            problems.Add($"{name}: description differs");
        var parameters = InputParameters(method).ToList();
        if (!root.TryGetProperty("inputSchema", out var schema) || !schema.TryGetProperty("properties", out var props))
        {
            if (parameters.Count > 0) problems.Add($"{name}: no inputSchema.properties");
            return;
        }

        var docNames = props.EnumerateObject().Select(p => p.Name).ToList();
        var codeNames = parameters.Select(p => p.Name!).ToList();
        foreach (var missing in codeNames.Except(docNames)) problems.Add($"{name}: parameter {missing} is missing in the document");
        foreach (var extra in docNames.Except(codeNames)) problems.Add($"{name}: the document has the parameter {extra}, the code has not");
        if (!docNames.Intersect(codeNames).SequenceEqual(codeNames.Intersect(docNames))) problems.Add($"{name}: parameter order differs");
        var required = new List<string>();
        foreach (var p in parameters)
        {
            if (!p.HasDefaultValue) required.Add(p.Name!);
            if (!props.TryGetProperty(p.Name!, out var node)) continue;
            if (node.TryGetProperty("type", out var type) && type.GetString() != JsonType(p.ParameterType))
                problems.Add($"{name}: type of {p.Name} is {type.GetString()}, the code says {JsonType(p.ParameterType)}");
            string? pd = p.GetCustomAttribute<DescriptionAttribute>()?.Description;
            if (pd is not null && (!node.TryGetProperty("description", out var dd) || dd.GetString() != pd))
                problems.Add($"{name}: description of {p.Name} differs");
            bool codeHasDefault = p.HasDefaultValue && p.DefaultValue is not null && p.DefaultValue != DBNull.Value;
            bool docHasDefault = node.TryGetProperty("default", out var docDefault);
            if (codeHasDefault != docHasDefault)
                problems.Add($"{name}: default of {p.Name}: the code {(codeHasDefault ? "has " + p.DefaultValue : "has none")}, the document {(docHasDefault ? "has " + docDefault : "has none")}");
            else if (codeHasDefault && !DefaultEquals(p.DefaultValue!, docDefault))
                problems.Add($"{name}: default of {p.Name} is {docDefault}, the code says {p.DefaultValue}");
        }

        var docRequired = schema.TryGetProperty("required", out var req) ? req.EnumerateArray().Select(e => e.GetString()!).ToList() : [];
        if (!docRequired.SequenceEqual(required)) problems.Add($"{name}: required is [{string.Join(", ", docRequired)}], the code says [{string.Join(", ", required)}]");
    }

    private static bool DefaultEquals(object code, JsonElement doc) => code switch
    {
        bool b => doc.ValueKind is JsonValueKind.True or JsonValueKind.False && doc.GetBoolean() == b,
        string s => doc.ValueKind == JsonValueKind.String && doc.GetString() == s,
        char c => doc.ValueKind == JsonValueKind.String && doc.GetString() == c.ToString(),
        Enum e => doc.ValueKind == JsonValueKind.String && doc.GetString() == e.ToString(),
        IConvertible n when doc.ValueKind == JsonValueKind.Number => Math.Abs(n.ToDouble(System.Globalization.CultureInfo.InvariantCulture) - doc.GetDouble()) < 1e-12,
        _ => false,
    };

    private static void CompareMarkdown(string name, MethodInfo method, string path, List<string> problems)
    {
        if (!File.Exists(path))
        {
            problems.Add($"{name}: no document {name}.md");
            return;
        }

        string text = File.ReadAllText(path);
        int at = text.IndexOf("## Input Schema", StringComparison.Ordinal);
        var parameters = InputParameters(method).ToList();
        if (at < 0)
        {
            if (parameters.Count > 0) problems.Add($"{name}: the .md has no Input Schema section");
            return;
        }

        int end = text.IndexOf("\n## ", at + 4, StringComparison.Ordinal);
        var rows = text[at..(end < 0 ? text.Length : end)].Split('\n').Where(l => l.StartsWith("| `", StringComparison.Ordinal)).ToList();
        var docNames = rows.Select(r => r.Split('|')[1].Trim().Trim('`')).ToList();
        var codeNames = parameters.Select(p => p.Name!).ToList();
        if (!docNames.SequenceEqual(codeNames))
            problems.Add($"{name}: the .md Input Schema table lists [{string.Join(", ", docNames)}], the code has [{string.Join(", ", codeNames)}]");
        else
            for (int i = 0; i < rows.Count; i++)
            {
                string cell = rows[i].Split('|')[3].Trim();
                if (cell != (parameters[i].HasDefaultValue ? "No" : "Yes"))
                    problems.Add($"{name}: Required cell of {codeNames[i]} in the .md is {cell}");
            }
    }

    private static string FindRepoRoot()
    {
        for (var dir = new DirectoryInfo(AppContext.BaseDirectory); dir is not null; dir = dir.Parent)
            if (Directory.Exists(Path.Combine(dir.FullName, "docs", "mcp", "tools")))
                return dir.FullName;
        throw new DirectoryNotFoundException("docs/mcp/tools not found above " + AppContext.BaseDirectory);
    }
}
