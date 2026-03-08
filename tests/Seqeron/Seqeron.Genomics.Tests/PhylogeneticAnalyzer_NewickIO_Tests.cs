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
        // Arrange - known Newick input with exact branch lengths
        string input = "(A:0.1,B:0.2);";
        var tree = PhylogeneticAnalyzer.ParseNewick(input);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree, includeBranchLengths: true);

        // Assert - verify Name:number format per Wikipedia grammar (Length → ":" number)
        Assert.Multiple(() =>
        {
            Assert.That(newick, Does.Contain(":"),
                "Branch lengths use colon notation (Wikipedia Examples)");
            Assert.That(newick, Does.Match(@"A:\d+\.\d+"),
                "Leaf A must have Name:number format");
            Assert.That(newick, Does.Match(@"B:\d+\.\d+"),
                "Leaf B must have Name:number format");
        });
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

        // Assert - verify structure AND leaf names per Wikipedia grammar
        Assert.Multiple(() =>
        {
            Assert.That(node, Is.Not.Null, "Root should not be null");
            Assert.That(node.IsLeaf, Is.False, "Root should be internal node");
            Assert.That(node.Left, Is.Not.Null, "Left child should exist");
            Assert.That(node.Right, Is.Not.Null, "Right child should exist");
            Assert.That(node.Left!.IsLeaf, Is.True, "Left should be leaf");
            Assert.That(node.Right!.IsLeaf, Is.True, "Right should be leaf");
            Assert.That(node.Left!.Name, Is.EqualTo("A"), "Left leaf name should be 'A'");
            Assert.That(node.Right!.Name, Is.EqualTo("B"), "Right leaf name should be 'B'");
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

        // Assert - verify full tree structure per Wikipedia grammar
        Assert.Multiple(() =>
        {
            Assert.That(node.IsLeaf, Is.False, "Root should be internal");
            Assert.That(node.Left!.IsLeaf, Is.False, "Left subtree should be internal");
            Assert.That(node.Right!.IsLeaf, Is.False, "Right subtree should be internal");

            // Left subtree: (A,B)
            Assert.That(node.Left!.Left!.Name, Is.EqualTo("A"), "Left-Left should be 'A'");
            Assert.That(node.Left!.Right!.Name, Is.EqualTo("B"), "Left-Right should be 'B'");

            // Right subtree: (C,D)
            Assert.That(node.Right!.Left!.Name, Is.EqualTo("C"), "Right-Left should be 'C'");
            Assert.That(node.Right!.Right!.Name, Is.EqualTo("D"), "Right-Right should be 'D'");
        });

        var leaves = PhylogeneticAnalyzer.GetLeaves(node).Select(l => l.Name).OrderBy(n => n).ToList();
        Assert.That(leaves, Is.EquivalentTo(new[] { "A", "B", "C", "D" }), "Should have exactly leaves A, B, C, D");
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
        // Arrange - Binary tree with names and branch lengths (Wikipedia full format)
        // Note: Implementation supports binary trees only
        string newick = "((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)Root;";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert - verify ALL fields: names AND branch lengths
        Assert.Multiple(() =>
        {
            // Node names
            Assert.That(node.Name, Is.EqualTo("Root"), "Root name should be 'Root'");
            Assert.That(node.Left!.Name, Is.EqualTo("AB"), "Left subtree should be 'AB'");
            Assert.That(node.Right!.Name, Is.EqualTo("CD"), "Right subtree should be 'CD'");

            // Branch lengths for internal nodes
            Assert.That(node.Left!.BranchLength, Is.EqualTo(0.3).Within(0.0001),
                "AB branch length should be 0.3");
            Assert.That(node.Right!.BranchLength, Is.EqualTo(0.6).Within(0.0001),
                "CD branch length should be 0.6");

            // Leaf names
            Assert.That(node.Left!.Left!.Name, Is.EqualTo("A"), "Left-Left should be 'A'");
            Assert.That(node.Left!.Right!.Name, Is.EqualTo("B"), "Left-Right should be 'B'");
            Assert.That(node.Right!.Left!.Name, Is.EqualTo("C"), "Right-Left should be 'C'");
            Assert.That(node.Right!.Right!.Name, Is.EqualTo("D"), "Right-Right should be 'D'");

            // Leaf branch lengths
            Assert.That(node.Left!.Left!.BranchLength, Is.EqualTo(0.1).Within(0.0001),
                "A branch length should be 0.1");
            Assert.That(node.Left!.Right!.BranchLength, Is.EqualTo(0.2).Within(0.0001),
                "B branch length should be 0.2");
            Assert.That(node.Right!.Left!.BranchLength, Is.EqualTo(0.4).Within(0.0001),
                "C branch length should be 0.4");
            Assert.That(node.Right!.Right!.BranchLength, Is.EqualTo(0.5).Within(0.0001),
                "D branch length should be 0.5");
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

    [Test]
    [Description("MUST: ParseNewick handles root branch length (Olsen: tree ==> descendant_list [root_label] [:branch_length] ;)")]
    public void ParseNewick_RootBranchLength_ExtractsValue()
    {
        // Arrange - Olsen grammar: tree ==> descendant_list [root_label] [:branch_length] ;
        // TreeAlign convention: root node with branch length 0.0
        // Source: Olsen (1990), also PHYLIP notes (TreeAlign writes root branch length)
        string newick = "(A:0.1,B:0.2):0.0;";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(node.BranchLength, Is.EqualTo(0.0).Within(0.0001),
                "Root branch length should be 0.0");
            Assert.That(node.Left!.BranchLength, Is.EqualTo(0.1).Within(0.0001),
                "Left branch length should be 0.1");
            Assert.That(node.Right!.BranchLength, Is.EqualTo(0.2).Within(0.0001),
                "Right branch length should be 0.2");
        });
    }

    [Test]
    [Description("MUST: ParseNewick handles root with name and branch length (Olsen grammar)")]
    public void ParseNewick_RootNameAndBranchLength_ParsesCorrectly()
    {
        // Arrange - Wikipedia example: '(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;' (all have distance to parent)
        // Adapted for binary: root with name and branch length
        string newick = "(A:0.1,B:0.2)Root:0.5;";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert
        Assert.Multiple(() =>
        {
            Assert.That(node.Name, Is.EqualTo("Root"), "Root name should be 'Root'");
            Assert.That(node.BranchLength, Is.EqualTo(0.5).Within(0.0001),
                "Root branch length should be 0.5");
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
    [Description("MUST: Round-trip preserves tree topology (sibling relationships)")]
    public void RoundTrip_PreservesTopology()
    {
        // Arrange - known topology: (A,B) siblings, (C,D) siblings
        // Use known Newick to avoid coupling with UPGMA implementation
        string input = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);";
        var tree = PhylogeneticAnalyzer.ParseNewick(input);

        // Act - round-trip: ToNewick → ParseNewick
        string newick = PhylogeneticAnalyzer.ToNewick(tree, includeBranchLengths: true);
        var parsed = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert - verify actual topology, not just leaf count
        var leftLeaves = PhylogeneticAnalyzer.GetLeaves(parsed.Left!)
            .Select(l => l.Name).OrderBy(n => n).ToList();
        var rightLeaves = PhylogeneticAnalyzer.GetLeaves(parsed.Right!)
            .Select(l => l.Name).OrderBy(n => n).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(parsed.IsLeaf, Is.False, "Root should be internal");
            Assert.That(leftLeaves, Is.EquivalentTo(new[] { "A", "B" }),
                "Left subtree must contain siblings A and B");
            Assert.That(rightLeaves, Is.EquivalentTo(new[] { "C", "D" }),
                "Right subtree must contain siblings C and D");
        });
    }

    [Test]
    [Description("SHOULD: Round-trip preserves branch lengths with tolerance")]
    public void RoundTrip_PreservesBranchLengths()
    {
        // Arrange - known Newick with exact branch lengths (not BuildTree)
        // Tests round-trip fidelity independent of UPGMA implementation
        string input = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);";
        var tree = PhylogeneticAnalyzer.ParseNewick(input);

        // Act - round-trip: ToNewick → ParseNewick
        string newick = PhylogeneticAnalyzer.ToNewick(tree, includeBranchLengths: true);
        var parsed = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert - verify INDIVIDUAL branch lengths, not just total
        // Tolerance = ±0.0001 (F4 format precision: ±0.00005 per round-trip)
        Assert.Multiple(() =>
        {
            Assert.That(parsed.Left!.Left!.BranchLength, Is.EqualTo(0.1).Within(0.0001),
                "A branch length should survive round-trip");
            Assert.That(parsed.Left!.Right!.BranchLength, Is.EqualTo(0.2).Within(0.0001),
                "B branch length should survive round-trip");
            Assert.That(parsed.Left!.BranchLength, Is.EqualTo(0.3).Within(0.0001),
                "Left subtree branch length should survive round-trip");
            Assert.That(parsed.Right!.Left!.BranchLength, Is.EqualTo(0.4).Within(0.0001),
                "C branch length should survive round-trip");
            Assert.That(parsed.Right!.Right!.BranchLength, Is.EqualTo(0.5).Within(0.0001),
                "D branch length should survive round-trip");
            Assert.That(parsed.Right!.BranchLength, Is.EqualTo(0.6).Within(0.0001),
                "Right subtree branch length should survive round-trip");
        });
    }

    [Test]
    [Description("MUST: Round-trip preserves internal node names (Wikipedia grammar: Internal → '(' BranchSet ')' Name)")]
    public void RoundTrip_FullFormat_PreservesInternalNames()
    {
        // Arrange - Wikipedia full format adapted to binary tree:
        // Original Wikipedia: (A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F; (trifurcating)
        // Binary equivalent with all features: names + branch lengths on all nodes
        string input = "((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)Root;";

        // Act
        var parsed = PhylogeneticAnalyzer.ParseNewick(input);
        string output = PhylogeneticAnalyzer.ToNewick(parsed, includeBranchLengths: true);
        var reparsed = PhylogeneticAnalyzer.ParseNewick(output);

        // Assert - Internal node names must survive round-trip
        Assert.Multiple(() =>
        {
            Assert.That(reparsed.Name, Is.EqualTo("Root"),
                "Root name must be preserved in round-trip");
            Assert.That(reparsed.Left!.Name, Is.EqualTo("AB"),
                "Left subtree name must be preserved in round-trip");
            Assert.That(reparsed.Right!.Name, Is.EqualTo("CD"),
                "Right subtree name must be preserved in round-trip");
        });

        // Also verify leaf names
        var leaves = PhylogeneticAnalyzer.GetLeaves(reparsed).Select(l => l.Name).OrderBy(n => n).ToList();
        Assert.That(leaves, Is.EquivalentTo(new[] { "A", "B", "C", "D" }),
            "Leaf names must be preserved in round-trip");
    }

    [Test]
    [Description("MUST: ToNewick emits internal node names for valid labels (Wikipedia grammar)")]
    public void ToNewick_WithValidInternalNames_EmitsNames()
    {
        // Arrange - Build a tree with valid internal node names (not UPGMA auto-generated)
        string input = "((A:0.1,B:0.2)AB:0.3,(C:0.4,D:0.5)CD:0.6)Root;";
        var tree = PhylogeneticAnalyzer.ParseNewick(input);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree, includeBranchLengths: true);

        // Assert - Internal node names should appear in output
        Assert.Multiple(() =>
        {
            Assert.That(newick, Does.Contain(")AB:"),
                "Output must contain internal name 'AB' after closing paren");
            Assert.That(newick, Does.Contain(")CD:"),
                "Output must contain internal name 'CD' after closing paren");
            Assert.That(newick, Does.EndWith("Root;"),
                "Output must end with root name 'Root' before semicolon");
        });
    }

    [Test]
    [Description("MUST: ToNewick omits invalid internal node names (Olsen: unquoted labels prohibit metacharacters)")]
    public void ToNewick_WithMetacharacterNames_OmitsInvalidNames()
    {
        // Arrange - UPGMA/NJ generates names like '(Human,Chimp)' containing Newick metacharacters.
        // Per Olsen spec, unquoted labels may not contain parentheses, commas, etc.
        var sequences = new Dictionary<string, string>
        {
            ["Human"] = "ACGTACGT",
            ["Chimp"] = "TCGTACGT"
        };
        var tree = PhylogeneticAnalyzer.BuildTree(sequences);

        // Act
        string newick = PhylogeneticAnalyzer.ToNewick(tree.Root, includeBranchLengths: true);

        // Assert - Auto-generated names with metacharacters should NOT appear as labels
        // The output should have ')' immediately followed by ':' or ';', not by a name containing '('
        Assert.That(newick, Does.Not.Match(@"\)[^:;,)]+\("),
            "Invalid metacharacter names must not appear in output");
        // But leaf names should still be present
        Assert.Multiple(() =>
        {
            Assert.That(newick, Does.Contain("Human"), "Leaf name 'Human' must be present");
            Assert.That(newick, Does.Contain("Chimp"), "Leaf name 'Chimp' must be present");
        });
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
    [Description("SHOULD: ParseNewick handles missing semicolon gracefully (Wikipedia grammar requires ';' but parser is lenient)")]
    public void ParseNewick_MissingSemicolon_ParsesCorrectly()
    {
        // Arrange - Wikipedia grammar: Tree → Subtree ";". Semicolon is required.
        // However, parser is lenient: it strips trailing ';' if present and parses regardless.
        // Source: Wikipedia Newick format § Grammar rules.
        string newick = "(A,B)";

        // Act
        var node = PhylogeneticAnalyzer.ParseNewick(newick);

        // Assert - verify full parse result, not just non-null
        Assert.Multiple(() =>
        {
            Assert.That(node, Is.Not.Null, "Root should not be null");
            Assert.That(node.IsLeaf, Is.False, "Root should be internal");
            Assert.That(node.Left!.Name, Is.EqualTo("A"), "Left leaf should be 'A'");
            Assert.That(node.Right!.Name, Is.EqualTo("B"), "Right leaf should be 'B'");
        });
    }

    #endregion
}
