using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class RestrictionDigestTests
{
    [Test]
    public void RestrictionDigest_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.restriction_digest("AAAGAATTCAAA", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest("", new[] { "EcoRI" }));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest("AAAGAATTCAAA", Array.Empty<string>()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.restriction_digest("AAAGAATTCAAA", null!));
    }

    /// <summary>
    /// Review 2026-09 B08 F48: isoschizomers cutting at the same positions (Biopython 1.88
    /// <c>MspI.search(Seq("CCGGAAACCGG")) == HpaII.search(...) == [2, 9]</c>) are all listed in
    /// <c>leftEnzymes</c>/<c>rightEnzymes</c> (ordinal order); the single-name fields are the first of
    /// the list, independent of the order of <c>enzyme_names</c>; both fields reach the JSON output.
    /// </summary>
    [Test]
    public void RestrictionDigest_CoincidingCuts_ListAllEnzymes()
    {
        var a = MolToolsTools.restriction_digest("CCGGAAACCGG", new[] { "MspI", "HpaII" }).Fragments;
        var b = MolToolsTools.restriction_digest("CCGGAAACCGG", new[] { "HpaII", "MspI" }).Fragments;
        string json = System.Text.Json.JsonSerializer.Serialize(
            a[1], new System.Text.Json.JsonSerializerOptions(System.Text.Json.JsonSerializerDefaults.Web));

        Assert.Multiple(() =>
        {
            Assert.That(a.Select(f => f.Sequence), Is.EqualTo(new[] { "C", "CGGAAAC", "CGG" }));
            Assert.That(a[0].LeftEnzymes, Is.Empty);
            Assert.That(a[0].RightEnzymes, Is.EqualTo(new[] { "HpaII", "MspI" }));
            Assert.That(a[1].LeftEnzymes, Is.EqualTo(new[] { "HpaII", "MspI" }));
            Assert.That(a[1].LeftEnzyme, Is.EqualTo("HpaII"));
            Assert.That(a[1].LeftOverhang, Is.EqualTo("CG"), "Biopython MspI/HpaII ovhgseq");
            Assert.That(b, Is.EqualTo(a));
            Assert.That(json, Does.Contain("\"leftEnzymes\":[\"HpaII\",\"MspI\"]"));
            Assert.That(json, Does.Contain("\"leftEnzyme\":\"HpaII\""));
        });
    }

    [Test]
    public void RestrictionDigest_Binding_InvokesSuccessfully()
    {
        // "AAAGAATTCAAA" (len 12), EcoRI cuts forward at position 4.
        // Fragments: [0,4)="AAAG", [4,12)="AATTCAAA".
        var fragments = MolToolsTools.restriction_digest("AAAGAATTCAAA", new[] { "EcoRI" }).Fragments;

        Assert.Multiple(() =>
        {
            Assert.That(fragments, Has.Count.EqualTo(2));

            Assert.That(fragments[0].Sequence, Is.EqualTo("AAAG"));
            Assert.That(fragments[0].StartPosition, Is.EqualTo(0));
            Assert.That(fragments[0].Length, Is.EqualTo(4));
            Assert.That(fragments[0].LeftEnzyme, Is.Null);
            Assert.That(fragments[0].RightEnzyme, Is.EqualTo("EcoRI"));

            Assert.That(fragments[1].Sequence, Is.EqualTo("AATTCAAA"));
            Assert.That(fragments[1].StartPosition, Is.EqualTo(4));
            Assert.That(fragments[1].Length, Is.EqualTo(8));
            Assert.That(fragments[1].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(fragments[1].RightEnzyme, Is.Null);
        });

        // No cut site -> single fragment equal to the whole sequence.
        var uncut = MolToolsTools.restriction_digest("AAAAAAAA", new[] { "EcoRI" }).Fragments;
        Assert.Multiple(() =>
        {
            Assert.That(uncut, Has.Count.EqualTo(1));
            Assert.That(uncut[0].Sequence, Is.EqualTo("AAAAAAAA"));
        });
    }

    /// <summary>
    /// The tool exposes <c>topology</c> and defaults it to <c>Linear</c>, so a circular digest is
    /// reachable over MCP. Reference (Biopython 1.88):
    /// <c>EcoRI.catalyse(Seq("TTCAAAAAAAGAA"), linear=False) == (Seq("AATTCAAAAAAAG"),)</c> while
    /// <c>EcoRI.search(Seq("TTCAAAAAAAGAA")) == []</c> — the <c>GAATTC</c> site spans the origin.
    /// </summary>
    [Test]
    public void RestrictionDigest_Topology_CircularFindsOriginSpanningSite()
    {
        var linear = MolToolsTools.restriction_digest("TTCAAAAAAAGAA", new[] { "EcoRI" }).Fragments;
        var circular = MolToolsTools
            .restriction_digest("TTCAAAAAAAGAA", new[] { "EcoRI" }, MoleculeTopology.Circular)
            .Fragments;

        Assert.Multiple(() =>
        {
            // Default is Linear (backward compatible): the site does not exist on a linear molecule.
            Assert.That(linear, Has.Count.EqualTo(1));
            Assert.That(linear[0].Sequence, Is.EqualTo("TTCAAAAAAAGAA"));
            Assert.That(linear[0].LeftEnzyme, Is.Null);

            Assert.That(circular, Has.Count.EqualTo(1));
            Assert.That(circular[0].Sequence, Is.EqualTo("AATTCAAAAAAAG"));
            Assert.That(circular[0].StartPosition, Is.EqualTo(11));
            Assert.That(circular[0].Length, Is.EqualTo(13));
            Assert.That(circular[0].LeftEnzyme, Is.EqualTo("EcoRI"));
            Assert.That(circular[0].RightEnzyme, Is.EqualTo("EcoRI"));
        });
    }

    /// <summary>
    /// Fragment ends created by a cut carry their blunt/5'/3' chemistry and overhang bases.
    /// Reference (Biopython 1.88): <c>EcoRI.ovhgseq == "AATT"</c> with <c>is_5overhang()</c>;
    /// <c>PstI.ovhgseq == "TGCA"</c> with <c>is_3overhang()</c>. An original terminus of a linear
    /// molecule is not a cut end, so its overhang fields are null.
    /// </summary>
    [Test]
    public void RestrictionDigest_ReportsFragmentEndOverhangs()
    {
        var ecoRI = MolToolsTools.restriction_digest("AAAGAATTCAAA", new[] { "EcoRI" }).Fragments;
        var pstI = MolToolsTools.restriction_digest("CTGCAGAAACTGCAG", new[] { "PstI" }).Fragments;

        Assert.Multiple(() =>
        {
            Assert.That(ecoRI[0].LeftOverhangType, Is.Null);
            Assert.That(ecoRI[0].LeftOverhang, Is.Null);
            Assert.That(ecoRI[0].RightOverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That(ecoRI[0].RightOverhang, Is.EqualTo("AATT"));
            Assert.That(ecoRI[1].LeftOverhangType, Is.EqualTo(OverhangType.FivePrime));
            Assert.That(ecoRI[1].LeftOverhang, Is.EqualTo("AATT"));
            Assert.That(ecoRI[1].RightOverhangType, Is.Null);

            Assert.That(pstI, Has.Count.EqualTo(3));
            Assert.That(pstI[0].RightOverhangType, Is.EqualTo(OverhangType.ThreePrime));
            Assert.That(pstI[0].RightOverhang, Is.EqualTo("TGCA"));
            Assert.That(pstI[1].LeftOverhang, Is.EqualTo("TGCA"));
        });
    }
}
