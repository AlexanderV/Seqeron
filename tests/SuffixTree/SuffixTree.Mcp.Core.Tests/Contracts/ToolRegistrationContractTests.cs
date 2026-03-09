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
    public void SuffixTreeGenomicsTools_RegisteredToolNames_AreStable()
    {
        var registered = GetRegisteredToolNames(typeof(SuffixTreeGenomicsTools));
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

    [Test]
    public void AllRegisteredTools_HaveTitleAndReadOnlyAnnotations()
    {
        var allToolTypes = new[] { typeof(SuffixTreeCoreTools), typeof(SuffixTreeGenomicsTools) };

        foreach (var toolType in allToolTypes)
        {
            var methods = toolType.GetMethods(BindingFlags.Public | BindingFlags.Static);
            foreach (var method in methods)
            {
                var attr = method.GetCustomAttributes(inherit: false)
                    .FirstOrDefault(a => string.Equals(a.GetType().Name, "McpServerToolAttribute", StringComparison.Ordinal));

                if (attr == null) continue;

                var titleProp = attr.GetType().GetProperty("Title", BindingFlags.Public | BindingFlags.Instance);
                var title = titleProp?.GetValue(attr) as string;
                Assert.That(title, Is.Not.Null.And.Not.Empty,
                    $"Tool {method.Name} on {toolType.Name} must have a non-empty Title");

                var readOnlyProp = attr.GetType().GetProperty("ReadOnly", BindingFlags.Public | BindingFlags.Instance);
                var readOnly = readOnlyProp?.GetValue(attr);
                Assert.That(readOnly, Is.Not.Null,
                    $"Tool {method.Name} on {toolType.Name} must have ReadOnly set");
            }
        }
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
