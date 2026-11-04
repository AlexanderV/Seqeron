using System.Reflection;
using System.Text.Json;
using ModelContextProtocol.Server;
using NUnit.Framework;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Review 2026-09 B08 D8 / F53: <c>partial_digest</c>. Reference values: Rosalind PDPL ("Creating a
/// Restriction Map") sample X = {0,2,4,7,10} → ΔX = {2,2,3,3,4,5,6,7,8,10}; Biopython 1.88
/// <c>RestrictionBatch([AluI, Sau3AI, TaqI]).search(Seq("AGCTGATCGA"))</c> → AluI [3], Sau3AI [5],
/// TaqI [8] (1-based), i.e. 0-based cuts {2,4,7} on a 10 bp molecule.
/// </summary>
[TestFixture]
public class PartialDigestTests
{
    // The MCP SDK's own serializer settings (camelCase, string enums, nulls omitted) — what a client receives.
    private static readonly JsonSerializerOptions Web = ModelContextProtocol.McpJsonUtilities.DefaultOptions;

    [Test]
    public void PartialDigest_Schema_ValidatesCorrectly()
    {
        Assert.Multiple(() =>
        {
            Assert.DoesNotThrow(() => MolToolsTools.partial_digest("AGCTGATCGA", new[] { "AluI" }));
            Assert.Throws<ArgumentException>(() => MolToolsTools.partial_digest("", new[] { "AluI" }));
            Assert.Throws<ArgumentException>(() => MolToolsTools.partial_digest("AGCTGATCGA", Array.Empty<string>()));
            Assert.Throws<ArgumentException>(() => MolToolsTools.partial_digest("AGCTGATCGA", null!));
            Assert.Throws<ArgumentOutOfRangeException>(() =>
                MolToolsTools.partial_digest("AGCTGATCGA", new[] { "AluI" }, max_missed_cuts: -1));
        });
    }

    [Test]
    public void PartialDigest_Binding_LengthsAreTheRosalindDifferenceMultiset()
    {
        var fragments = MolToolsTools.partial_digest("AGCTGATCGA", new[] { "AluI", "Sau3AI", "TaqI" }).Fragments;

        Assert.Multiple(() =>
        {
            Assert.That(fragments.Select(f => f.Fragment.Length).OrderBy(x => x),
                Is.EqualTo(new[] { 2, 2, 3, 3, 4, 5, 6, 7, 8, 10 }));
            Assert.That(fragments.Select(f => f.MissedCuts), Is.EqualTo(new[] { 0, 0, 0, 0, 1, 1, 1, 2, 2, 3 }));
            // The missedCuts = 0 prefix is restriction_digest's output, member for member.
            Assert.That(fragments.Where(f => f.MissedCuts == 0).Select(f => f.Fragment),
                Is.EqualTo(MolToolsTools.restriction_digest("AGCTGATCGA", new[] { "AluI", "Sau3AI", "TaqI" }).Fragments));
        });
    }

    [Test]
    public void PartialDigest_TopologyAndLimit_AreForwarded()
    {
        var circular = MolToolsTools.partial_digest(
            "AGCTGATCGA", new[] { "AluI", "Sau3AI", "TaqI" }, MoleculeTopology.Circular).Fragments;
        var limited = MolToolsTools.partial_digest(
            "AGCTGATCGA", new[] { "AluI", "Sau3AI", "TaqI" }, MoleculeTopology.Circular, max_missed_cuts: 1).Fragments;
        var library = RestrictionAnalyzer.PartialDigest(
            new DnaSequence("AGCTGATCGA"), MoleculeTopology.Circular, "AluI", "Sau3AI", "TaqI").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(circular, Is.EqualTo(library));
            Assert.That(circular, Has.Count.EqualTo(10), "k² + 1 for k = 3");
            Assert.That(circular[^1].Fragment.Sequence, Is.EqualTo("AGCTGATCGA"), "intact circle last");
            Assert.That(limited, Is.EqualTo(library.Take(6)));
        });
    }

    [Test]
    public void PartialDigest_DocPair_SchemaDescriptionAndExamples_MatchTheTool()
    {
        string dir = DocsDirectory();
        Assert.That(File.Exists(Path.Combine(dir, "partial_digest.md")), Is.True, "partial_digest.md missing");
        using var doc = JsonDocument.Parse(File.ReadAllText(Path.Combine(dir, "partial_digest.mcp.json")));
        var root = doc.RootElement;
        var method = typeof(MolToolsTools).GetMethod(nameof(MolToolsTools.partial_digest))!;
        var advertised = McpServerTool.Create(method, target: null, new McpServerToolCreateOptions { Name = "partial_digest" })
            .ProtocolTool.InputSchema;

        Assert.Multiple(() =>
        {
            Assert.That(root.GetProperty("toolName").GetString(), Is.EqualTo("partial_digest"));
            Assert.That(root.GetProperty("description").GetString(),
                Is.EqualTo(method.GetCustomAttribute<System.ComponentModel.DescriptionAttribute>()!.Description), "description");
            var have = root.GetProperty("inputSchema").GetProperty("properties").EnumerateObject().ToList();
            var want = advertised.GetProperty("properties").EnumerateObject().ToList();
            Assert.That(have.Select(p => p.Name), Is.EqualTo(want.Select(p => p.Name)), "parameter names / order");
            foreach (var (w, h) in want.Zip(have))
            {
                Assert.That(h.Value.GetProperty("description").GetString(),
                    Is.EqualTo(w.Value.GetProperty("description").GetString()), $"{w.Name}: description");
            }
            Assert.That(root.GetProperty("inputSchema").GetProperty("required").EnumerateArray().Select(x => x.GetString()),
                Is.EquivalentTo(advertised.GetProperty("required").EnumerateArray().Select(x => x.GetString())), "required");
            Assert.That(root.GetProperty("outputSchema").GetProperty("properties").EnumerateObject().Select(p => p.Name),
                Is.EqualTo(new[] { "fragments" }));

            foreach (var example in root.GetProperty("examples").EnumerateArray())
            {
                string name = example.GetProperty("name").GetString()!;
                var input = example.GetProperty("input");
                var output = example.GetProperty("output");
                if (output.TryGetProperty("error", out var error))
                {
                    Assert.That(() => Invoke(input), Throws.InstanceOf<ArgumentException>()
                        .With.Message.StartsWith(error.GetString()!), name);
                    continue;
                }
                Assert.That(JsonSerializer.Serialize(Invoke(input), Web), Is.EqualTo(JsonSerializer.Serialize(output)), name);
            }
        });
    }

    private static PartialDigestResult Invoke(JsonElement input) =>
        MolToolsTools.partial_digest(
            input.GetProperty("sequence").GetString()!,
            input.GetProperty("enzyme_names").EnumerateArray().Select(e => e.GetString()!).ToArray(),
            input.TryGetProperty("topology", out var t) ? Enum.Parse<MoleculeTopology>(t.GetString()!) : MoleculeTopology.Linear,
            input.TryGetProperty("max_missed_cuts", out var m) ? m.GetInt32() : null);

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
}
