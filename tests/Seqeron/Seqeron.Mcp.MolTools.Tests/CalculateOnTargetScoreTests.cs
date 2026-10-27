using System.Reflection;
using System.Text.Json;
using ModelContextProtocol.Server;
using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Pins the two published on-target efficacy scorers exposed by the MolTools server. The expected
/// values are the reference ones locked by Seqeron.Genomics.Tests (Doench 2014: the two worked
/// examples distributed with CRISPOR's doenchScore.py; Rule Set 2: a row of the verified Azimuth
/// oracle on which the reference and Microsoft's own fixture agree).
/// </summary>
[TestFixture]
public class CalculateOnTargetScoreTests
{
    private const string Doench2014Example1 = "TATAGCTGCGATCTGAGGTAGGGAGGGACC"; // reference 0.713089368437
    private const string Doench2014Example2 = "TCCGCACCTGTCACGGTCGGGGCTTGGCGC"; // reference 0.0189838463593
    private const string Oracle30Mer = "AAAAAAAAAAAAAAAGCTAACAGCAGGAGT"; // nopos_oracle ref_score 0.529037

    [Test]
    public void CalculateOnTargetDoench2014_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.calculate_on_target_doench2014(Doench2014Example1));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_doench2014(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_doench2014(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_doench2014("ACGTACGT"));
        // 30-mer without an NGG PAM at offsets 25-26.
        Assert.Throws<ArgumentException>(
            () => MolToolsTools.calculate_on_target_doench2014("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"));
    }

    [Test]
    public void CalculateOnTargetDoench2014_Binding_ReproducesTheReferenceExamples()
    {
        Assert.That(MolToolsTools.calculate_on_target_doench2014(Doench2014Example1).Score,
            Is.EqualTo(71.30893469547244).Within(1e-9)); // doenchScore.py calcDoenchScore, executed (its docstring literal 71.3089368437 differs)
        Assert.That(MolToolsTools.calculate_on_target_doench2014(Doench2014Example2).Score,
            Is.EqualTo(1.8983843193336125).Within(1e-9)); // executed (docstring literal 1.89838463593 differs)
    }

    [Test]
    public void CalculateOnTargetRuleSet2_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer));
        Assert.DoesNotThrow(() => MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer, 10, 25.0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_rule_set2(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_rule_set2(null!));
        Assert.Throws<ArgumentException>(() => MolToolsTools.calculate_on_target_rule_set2("ACGTACGT"));
        // The gene-context arguments must be supplied together.
        Assert.Throws<ArgumentException>(
            () => MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer, 10, null));
        Assert.Throws<ArgumentException>(
            () => MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer, null, 25.0));
    }

    [Test]
    public void CalculateOnTargetRuleSet2_Binding_ReproducesTheVerifiedOracleRow()
    {
        Assert.That(MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer).Score,
            Is.EqualTo(0.529037).Within(1e-5));

        // The gene-context (full) model is a different model and must give a different value.
        Assert.That(MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer, 10, 25.0).Score,
            Is.Not.EqualTo(MolToolsTools.calculate_on_target_rule_set2(Oracle30Mer).Score));
    }

    // B08 WP5/D14: the doc pair (docs/mcp/tools/moltools/{tool}.md + .mcp.json) is pinned to the tool — the documented
    // input schema is the one the MCP SDK advertises, the output schema is the result record, and every documented
    // example is a real tool output (numeric examples bit-exact, error examples by message).
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

    private static OnTargetScoreResultInvoker Invoker(string tool) => tool switch
    {
        "calculate_on_target_doench2014" => input => MolToolsTools.calculate_on_target_doench2014(
            input.GetProperty("context_30mer").GetString()!).Score,
        "calculate_on_target_rule_set2" => input => MolToolsTools.calculate_on_target_rule_set2(
            input.GetProperty("context_30mer").GetString()!,
            input.TryGetProperty("amino_acid_cut_position", out var aa) ? aa.GetInt32() : null,
            input.TryGetProperty("percent_peptide", out var pp) ? pp.GetDouble() : null).Score,
        _ => throw new ArgumentOutOfRangeException(nameof(tool)),
    };

    private delegate double OnTargetScoreResultInvoker(JsonElement input);

    [TestCase("calculate_on_target_doench2014")]
    [TestCase("calculate_on_target_rule_set2")]
    public void DocPair_SchemaAndExamples_MatchTheTool(string tool)
    {
        string dir = DocsDirectory();
        Assert.That(File.Exists(Path.Combine(dir, tool + ".md")), Is.True, $"{tool}.md missing");
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.Combine(dir, tool + ".mcp.json")));
        var root = doc.RootElement;
        Assert.That(root.GetProperty("toolName").GetString(), Is.EqualTo(tool));

        var method = typeof(MolToolsTools).GetMethods(BindingFlags.Public | BindingFlags.Static)
            .Single(m => m.GetCustomAttribute<McpServerToolAttribute>()?.Name == tool);
        var advertised = McpServerTool.Create(method, target: null, new McpServerToolCreateOptions { Name = tool })
            .ProtocolTool.InputSchema;
        var documented = root.GetProperty("inputSchema");
        Assert.That(documented.GetProperty("properties").EnumerateObject().Select(p => p.Name),
            Is.EqualTo(advertised.GetProperty("properties").EnumerateObject().Select(p => p.Name)), "parameter names / order");
        Assert.That(documented.GetProperty("required").EnumerateArray().Select(x => x.GetString()),
            Is.EquivalentTo(advertised.GetProperty("required").EnumerateArray().Select(x => x.GetString())), "required");
        Assert.That(root.GetProperty("description").GetString(),
            Is.EqualTo(method.GetCustomAttribute<System.ComponentModel.DescriptionAttribute>()!.Description), "description");
        Assert.That(root.GetProperty("outputSchema").GetProperty("properties").EnumerateObject().Select(p => p.Name),
            Is.EqualTo(new[] { "score" }));

        var invoke = Invoker(tool);
        Assert.Multiple(() =>
        {
            foreach (var example in root.GetProperty("examples").EnumerateArray())
            {
                string name = example.GetProperty("name").GetString()!;
                var input = example.GetProperty("input");
                var output = example.GetProperty("output");
                if (output.TryGetProperty("score", out var score))
                    Assert.That(invoke(input), Is.EqualTo(score.GetDouble()), name);
                else
                    Assert.That(() => invoke(input), Throws.InstanceOf<ArgumentException>()
                        .With.Message.StartsWith(output.GetProperty("error").GetString()!), name);
            }
        });
    }
}
