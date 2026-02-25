using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using NUnit.Framework;
using SuffixTree.Mcp.Core.Tools;

namespace SuffixTree.Mcp.Core.Tests;

[TestFixture]
public class ToolRegistrationContractTests
{
    [Test]
    public void SuffixTreeCoreTools_RegisteredToolNames_AreStable()
    {
        var registered = GetRegisteredToolNames(typeof(SuffixTreeCoreTools));
        var expected = new[]
        {
            "suffix_tree_contains",
            "suffix_tree_count",
            "suffix_tree_find_all",
            "suffix_tree_lrs",
            "suffix_tree_lcs",
            "suffix_tree_stats"
        };

        Assert.That(registered.Values, Is.EquivalentTo(expected));
        Assert.That(registered.Count, Is.EqualTo(expected.Length));
    }

    [Test]
    public void SuffixTreeCompatibilityWrappers_AreNotDirectlyRegistered()
    {
        var suffixWrapperMethods = new[]
        {
            nameof(SuffixTreeTools.SuffixTreeContains),
            nameof(SuffixTreeTools.SuffixTreeCount),
            nameof(SuffixTreeTools.SuffixTreeFindAll),
            nameof(SuffixTreeTools.SuffixTreeLrs),
            nameof(SuffixTreeTools.SuffixTreeLcs),
            nameof(SuffixTreeTools.SuffixTreeStats)
        };

        foreach (string methodName in suffixWrapperMethods)
        {
            var method = typeof(SuffixTreeTools).GetMethod(methodName, BindingFlags.Public | BindingFlags.Static);
            Assert.That(method, Is.Not.Null, $"Method not found: {methodName}");
            Assert.That(GetRegisteredToolName(method!), Is.Null, $"Compatibility wrapper must not be MCP-registered: {methodName}");
        }
    }

    [Test]
    public void GenomicsTools_RegisteredToolNames_AreStable()
    {
        var registered = GetRegisteredToolNames(typeof(SuffixTreeTools));
        var expected = new[]
        {
            "find_longest_repeat",
            "find_longest_common_region",
            "calculate_similarity",
            "hamming_distance",
            "edit_distance",
            "count_approximate_occurrences"
        };

        Assert.That(registered.Values, Is.EquivalentTo(expected));
        Assert.That(registered.Count, Is.EqualTo(expected.Length));
    }

    private static Dictionary<string, string> GetRegisteredToolNames(Type toolsType)
    {
        var result = new Dictionary<string, string>(StringComparer.Ordinal);
        var methods = toolsType.GetMethods(BindingFlags.Public | BindingFlags.Static);
        foreach (var method in methods)
        {
            string? toolName = GetRegisteredToolName(method);
            if (toolName == null)
                continue;

            result[method.Name] = toolName;
        }

        return result;
    }

    private static string? GetRegisteredToolName(MethodInfo method)
    {
        object? attribute = method
            .GetCustomAttributes(inherit: false)
            .FirstOrDefault(a => string.Equals(a.GetType().Name, "McpServerToolAttribute", StringComparison.Ordinal));

        if (attribute == null)
            return null;

        var nameProperty = attribute.GetType().GetProperty("Name", BindingFlags.Public | BindingFlags.Instance);
        var explicitName = nameProperty?.GetValue(attribute) as string;
        return string.IsNullOrWhiteSpace(explicitName) ? method.Name : explicitName;
    }
}
