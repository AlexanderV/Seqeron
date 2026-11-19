using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>multiple_align_clustalw: type='auto' (ClustalW's guess from the second sequence) and clustering='upgma' (audit A2-14, B09 F92).</summary>
[TestFixture]
public class ClustalWAutoTypeUpgmaTests
{
    private static readonly string[] Dna = ["ACGTACGTAAGCT", "ACGTACGTCAGCT", "ACGTTCGTCAGGT", "TCGTTCGTCAGGT"];

    [Test]
    public void Auto_ResolvesFromTheSecondSequence_LikeClustalW()
    {
        var auto = AlignmentTools.MultipleAlignClustalW(Dna, type: "auto");
        var dna = AlignmentTools.MultipleAlignClustalW(Dna, type: "dna");
        Assert.That(auto.AlignedSequences, Is.EqualTo(dna.AlignedSequences));

        // protein first sequence, DNA second: ClustalW decides DNA from the second sequence (it printed 'bp' for all three), so the first row is then read as (invalid) DNA
        string[] mixed = ["MKVLAAGIV", "ACGTACGTAC", "MKVLAAGIV"];
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(mixed, type: "auto"), "decided DNA by the second sequence; 'M', 'K', 'V', 'L' ... of the first are not DNA here");
        string[] protein = ["ACGTACGTAC", "MKVLAAGIVL", "MKVLAGIVLL"];
        Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignClustalW(protein, type: "auto"), "second sequence is protein -> PROTEIN");
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(Dna, type: "rna"));
    }

    [Test]
    public void Clustering_Upgma_IsAcceptedAndKeepsTheSequences_UnknownIsRejected()
    {
        var r = AlignmentTools.MultipleAlignClustalW(Dna, clustering: "upgma");
        Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(Dna));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(Dna, clustering: "ward"));
    }
}
