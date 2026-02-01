using NUnit.Framework;
using Seqeron.Genomics;
using System;
using System.Collections.Generic;
using System.Linq;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Canonical tests for PHYLO-NEWICK-001: Newick I/O.
/// Tests ToNewick() and ParseNewick() methods in PhylogeneticAnalyzer.
/// 
/// Evidence: Wikipedia (Newick format), PHYLIP documentation (Felsenstein)
/// TestSpec: TestSpecs/PHYLO-NEWICK-001.md
/// Algorithm: docs/algorithms/Phylogenetics/Newick_Format.md
/// </summary>
[TestFixture]
[Category("PHYLO-NEWICK-001")]
[Description("Newick format parsing and export tests")]
public class PhylogeneticAnalyzer_NewickIO_Tests
{
    #region ToNewick Tests

    [Test]
    [Description("MUST: ToNewick output ends with semicolon (Wikipedia Grammar: Tree → Subtree ';')")]
    public void ToNewick_SimpleTree_EndsWithSemicolon()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);

        // Assert
        Assert.That(newick, Does.EndWith(";"),
            "Newick format must end with semicolon (Wikipedia Grammar)");
    }

    [Test]
    [Description("MUST: ToNewick with branch lengths includes colons (Wikipedia: 'A:0.1,B:0.2')")]
    public void ToNewick_WithBranchLengths_IncludesColons()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root, includeBranchLengths: true);

        // Assert
        Assert.That(newick, Does.Contain(":"),
            "Branch lengths use colon notation (Wikipedia Examples)");
    }

    [Test]
    [Description("MUST: ToNewick without branch lengths has no colons")]
    public void ToNewick_WithoutBranchLengths_NoColons()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGT",
            ["B"] = "TCGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root, includeBranchLengths: false);

        // Assert
        Assert.That(newick, Does.Not.Contain(":"),
            "Without branch lengths, no colons should appear");
    }

    [Test]
    [Description("MUST: ToNewick output contains all leaf names")]
    public void ToNewick_ContainsAllLeafNames()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["Human"] = "ACGTACGT",
            ["Chimp"] = "ACGTACGA",
            ["Mouse"] = "TCGTACGT",
            ["Rat"] = "TCGTACGA"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(newick, Does.Contain("Human"), "Must contain 'Human'");
            Assert.That(newick, Does.Contain("Chimp"), "Must contain 'Chimp'");
            Assert.That(newick, Does.Contain("Mouse"), "Must contain 'Mouse'");
            Assert.That(newick, Does.Contain("Rat"), "Must contain 'Rat'");
        });
    }

    [Test]
    [Description("MUST: ToNewick on null node returns empty string")]
    public void ToNewick_NullNode_ReturnsEmptyString()
    {
        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(null!);

        // Assert
        Assert.That(newick, Is.Empty,
            "ToNewick(null) should return empty string");
    }

    [Test]
    [Description("SHOULD: ToNewick handles large trees correctly")]
    public void ToNewick_LargeTree_ProducesValidFormat()
    {
        // Arrange - 8 taxa
        var sequences = new Dictionary<string, string>
        {
            ["Taxon1"] = "ACGTACGTACGT",
            ["Taxon2"] = "ACGTACGTACGA",
            ["Taxon3"] = "ACGTACGTACGC",
            ["Taxon4"] = "TCGTACGTACGT",
            ["Taxon5"] = "TCGTACGTACGA",
            ["Taxon6"] = "GCGTACGTACGT",
            ["Taxon7"] = "GCGTACGTACGA",
            ["Taxon8"] = "CCGTACGTACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(newick, Does.EndWith(";"), "Must end with semicolon");
            Assert.That(newick.Count(c => c == '('), Is.EqualTo(newick.Count(c => c == ')')),
                "Parentheses must be balanced");
            for (int i = 1; i <= 8; i++)
            {
                Assert.That(newick, Does.Contain($"Taxon{i}"), $"Must contain Taxon{i}");
            }
        });
    }

    #endregion

    #region ParseNewick Tests

    [Test]
    [Description("MUST: ParseNewick parses simple binary tree (Wikipedia: '(A,B);')")]
    public void ParseNewick_SimpleBinaryTree_ParsesCorrectly()
    {
        // Arrange - Wikipedia example
        string newick = "(A,B);";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(node, Is.Not.Null, "Root should not be null");
            Assert.That(node.IsLeaf, Is.False, "Root should be internal node");
            Assert.That(node.Left, Is.Not.Null, "Left child should exist");
            Assert.That(node.Right, Is.Not.Null, "Right child should exist");
            Assert.That(node.Left!.IsLeaf, Is.True, "Left should be leaf");
            Assert.That(node.Right!.IsLeaf, Is.True, "Right should be leaf");
        });
    }

    [Test]
    [Description("MUST: ParseNewick extracts branch lengths (Wikipedia: '(A:0.1,B:0.2);')")]
    public void ParseNewick_WithBranchLengths_ExtractsValues()
    {
        // Arrange - Wikipedia popular format example
        string newick = "(A:0.1,B:0.2);";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(node.Left!.BranchLength, Is.EqualTo(0.1).Within(0.0001),
                "Left branch length should be 0.1");
            Assert.That(node.Right!.BranchLength, Is.EqualTo(0.2).Within(0.0001),
                "Right branch length should be 0.2");
            Assert.That(node.Left.Name, Is.EqualTo("A"), "Left name should be A");
            Assert.That(node.Right.Name, Is.EqualTo("B"), "Right name should be B");
        });
    }

    [Test]
    [Description("MUST: ParseNewick handles nested trees (Wikipedia: '((A,B),(C,D));')")]
    public void ParseNewick_NestedBinaryTree_ParsesRecursively()
    {
        // Arrange - Wikipedia example
        string newick = "((A,B),(C,D));";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(node.IsLeaf, Is.False, "Root should be internal");
            Assert.That(node.Left!.IsLeaf, Is.False, "Left subtree should be internal");
            Assert.That(node.Right!.IsLeaf, Is.False, "Right subtree should be internal");
        });

        var leaves = PhylogeneticAnalyzer.GetLeaves(node).Select(l => l.Name).ToList();
        Assert.That(leaves, Has.Count.EqualTo(4), "Should have 4 leaves");
    }

    [Test]
    [Description("MUST: ParseNewick preserves leaf count")]
    public void ParseNewick_LeafCountMatchesInput()
    {
        // Arrange - test various sizes
        var testCases = new[]
        {
            ("(A,B);", 2),
            ("((A,B),C);", 3),
            ("((A,B),(C,D));", 4),
            ("(((A,B),C),(D,E));", 5)
        };

        foreach (var (newick, expectedCount) in testCases)
        {
            // Act
            var node = PhylogeneticAnalyzer.ParseNewick(newick);
            var leaves = PhylogeneticAnalyzer.GetLeaves(node).ToList();

            // Assert
            Assert.That(leaves, Has.Count.EqualTo(expectedCount),
                $"Tree '{newick}' should have {expectedCount} leaves");
        }
    }

    [Test]
    [Description("SHOULD: ParseNewick extracts internal node names")]
    public void ParseNewick_InternalNodeNames_ExtractsName()
    {
        // Arrange - simple binary tree with internal node name
        string newick = "(A,B)Root;";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.That(node.Name, Is.EqualTo("Root"),
            "Internal node name should be 'Root'");
    }

    [Test]
    [Description("SHOULD: ParseNewick handles binary tree with all features")]
    public void ParseNewick_FullFormat_ParsesAllFields()
    {
        // Arrange - Binary tree with names and branch lengths
        // Note: Implementation supports binary trees only
        string newick = "((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)Root;";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(node.Name, Is.EqualTo("Root"), "Root name should be 'Root'");
            Assert.That(node.Left!.Name, Is.EqualTo("AB"), "Left subtree should be 'AB'");
            Assert.That(node.Right!.Name, Is.EqualTo("CD"), "Right subtree should be 'CD'");
        });

        var leaves = PhylogeneticAnalyzer.GetLeaves(node).Select(l => l.Name).OrderBy(n => n).ToList();
        Assert.That(leaves, Is.EquivalentTo(new[] { "A", "B", "C", "D" }),
            "Should contain all four leaf names");
    }

    [Test]
    [Description("COULD: ParseNewick handles single taxon (PHYLIP: 'A;')")]
    public void ParseNewick_SingleTaxon_ParsesSingleNode()
    {
        // Arrange - PHYLIP example
        string newick = "A;";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(node, Is.Not.Null, "Node should not be null");
            Assert.That(node.Name, Is.EqualTo("A"), "Name should be A");
            Assert.That(node.IsLeaf, Is.True, "Single taxon is a leaf");
        });
    }

    #endregion

    #region Round-Trip Tests

    [Test]
    [Description("MUST: Round-trip preserves all leaf names")]
    public void RoundTrip_PreservesLeafNames()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["Human"] = "ACGTACGT",
            ["Chimp"] = "ACGTACGA",
            ["Gorilla"] = "ACGTACGC",
            ["Mouse"] = "TCGTACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);
        var parsed = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        var originalLeaves = PhylogeneticAnalyzer.GetLeaves(tree.Root)
            .Select(l => l.Name).OrderBy(n => n).ToList();
        var parsedLeaves = PhylogeneticAnalyzer.GetLeaves(parsed)
            .Select(l => l.Name).OrderBy(n => n).ToList();

        Assert.That(parsedLeaves, Is.EquivalentTo(originalLeaves),
            "Round-trip must preserve all leaf names");
    }

    [Test]
    [Description("MUST: Round-trip preserves tree topology (leaf count)")]
    public void RoundTrip_PreservesTopology()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "AAAAAAAA",
            ["B"] = "AAAAAAAC",
            ["C"] = "CCCCCCCC",
            ["D"] = "CCCCCCCA",
            ["E"] = "GGGGGGGG"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root);
        var parsed = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        int originalLeafCount = PhylogeneticAnalyzer.GetLeaves(tree.Root).Count();
        int parsedLeafCount = PhylogeneticAnalyzer.GetLeaves(parsed).Count();

        Assert.That(parsedLeafCount, Is.EqualTo(originalLeafCount),
            "Round-trip must preserve leaf count (topology indicator)");
    }

    [Test]
    [Description("SHOULD: Round-trip preserves branch lengths with tolerance")]
    public void RoundTrip_PreservesBranchLengths()
    {
        // Arrange
        var sequences = new Dictionary<string, string>
        {
            ["A"] = "ACGTACGT",
            ["B"] = "TCGTACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);
        double originalLength = PhylogeneticAnalyzer.CalculateTreeLength(tree.Root);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root, includeBranchLengths: true);
        var parsed = PhylogeneticAnalyzer.ParseNewick(newick);
        double parsedLength = PhylogeneticAnalyzer.CalculateTreeLength(parsed);

        // Assert - allow for formatting precision loss
        Assert.That(parsedLength, Is.EqualTo(originalLength).Within(0.001),
            "Round-trip should preserve total branch length within tolerance");
    }

    #endregion

    #region Edge Cases

    [Test]
    [Description("MUST: ParseNewick throws on empty string")]
    public void ParseNewick_EmptyString_ThrowsArgumentException()
    {
        // Act & Assert
        var ex = Assert.Throws<ArgumentException>(() =>
            PhylogeneticAnalyzer.ParseNewick(""));

        Assert.That(ex!.Message, Does.Contain("empty").IgnoreCase,
            "Exception message should indicate empty input");
    }

    [Test]
    [Description("MUST: ParseNewick throws on null string")]
    public void ParseNewick_NullString_ThrowsArgumentException()
    {
        // Act & Assert
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticAnalyzer.ParseNewick(null!));
    }

    [Test]
    [Description("SHOULD: ParseNewick handles whitespace-only as empty")]
    public void ParseNewick_WhitespaceOnly_ThrowsArgumentException()
    {
        // Act & Assert
        Assert.Throws<ArgumentException>(() =>
            PhylogeneticAnalyzer.ParseNewick("   "));
    }

    [Test]
    [Description("SHOULD: ParseNewick handles trailing whitespace gracefully")]
    public void ParseNewick_TrailingWhitespace_ParsesCorrectly()
    {
        // Arrange
        string newick = "(A,B);  ";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.That(node, Is.Not.Null);
        Assert.That(PhylogeneticAnalyzer.GetLeaves(node).Count(), Is.EqualTo(2));
    }

    [Test]
    [Description("SHOULD: ParseNewick handles missing semicolon gracefully")]
    public void ParseNewick_MissingSemicolon_ParsesCorrectly()
    {
        // Arrange - implementation strips trailing semicolon, should handle absence
        string newick = "(A,B)";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.That(node, Is.Not.Null);
        Assert.That(PhylogeneticAnalyzer.GetLeaves(node).Count(), Is.EqualTo(2));
    }

    #endregion
}
