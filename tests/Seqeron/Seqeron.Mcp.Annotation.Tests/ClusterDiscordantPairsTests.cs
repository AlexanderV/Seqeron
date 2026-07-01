using NUnit.Framework;
using Seqeron.Mcp.Annotation.Tools;

namespace Seqeron.Mcp.Annotation.Tests;

[TestFixture]
public class ClusterDiscordantPairsTests
{
    private static ReadPairSignatureDto Pair(string id, int p1, int p2) =>
        new(id, "chr1", p1, '+', "chr2", p2, '-', 400, true);

    // Three interchromosomal (chr1->chr2) pairs at nearby anchors -> one translocation cluster.
    private static List<ReadPairSignatureDto> Cluster() =>
        new() { Pair("r1", 100, 500), Pair("r2", 110, 505), Pair("r3", 120, 510) };

    [Test]
    public void ClusterDiscordantPairs_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnnotationTools.ClusterDiscordantPairs(Cluster()));
        Assert.Throws<ArgumentNullException>(() => AnnotationTools.ClusterDiscordantPairs(null!));
    }

    [Test]
    public void ClusterDiscordantPairs_Binding_InvokesSuccessfully()
    {
        // StructuralVariantAnalyzer.ClusterDiscordantPairs: nearby same-chromosome-pair anchors within
        // clusterDistance and >= minSupport (3) form one SV. Interchromosomal -> Translocation.
        var result = AnnotationTools.ClusterDiscordantPairs(Cluster());

        Assert.That(result.Variants, Has.Count.EqualTo(1));
        Assert.Multiple(() =>
        {
            Assert.That(result.Variants[0].Type, Is.EqualTo("Translocation"));
            Assert.That(result.Variants[0].SupportingReads, Is.EqualTo(3));
        });
    }

    [Test]
    public void ClusterDiscordantPairs_BelowMinSupport_ReturnsEmpty()
    {
        // Only two supporting pairs (< default minSupport 3) -> no SV emitted.
        var pairs = new List<ReadPairSignatureDto> { Pair("r1", 100, 500), Pair("r2", 110, 505) };
        var result = AnnotationTools.ClusterDiscordantPairs(pairs);
        Assert.That(result.Variants, Is.Empty);
    }
}
