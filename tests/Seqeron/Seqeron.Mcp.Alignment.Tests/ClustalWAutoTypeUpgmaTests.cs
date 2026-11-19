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

    [Test]
    public void Iteration_MatchesTheLibrary_AndIsValidated()
    {
        // a set on which the 'remove first' iteration of clustalw 2.1 changes the rows (clustalw_iteration_oracle.tsv row 350 family): the tool must equal the library call
        string[] seqs = ["AACCCCCAAACAACCCCAA", "AAACCCCCAACAACCCCAA", "AAACCCCAAACAACCCCAA", "CATACTTGACT", "AATCACCGAAACAACCCCAA"];
        var plain = AlignmentTools.MultipleAlignClustalW(seqs, distance: "clustalw-quick");
        foreach (var mode in new[] { "tree", "alignment" })
        {
            var tool = AlignmentTools.MultipleAlignClustalW(seqs, distance: "clustalw-quick", iteration: mode, numIterations: 3);
            var lib = Seqeron.Genomics.Alignment.SequenceAligner.MultipleAlignClustalW(seqs, new Seqeron.Genomics.Infrastructure.ClustalWOptions
            {
                PairwiseMode = Seqeron.Genomics.Infrastructure.ClustalWPairwiseMode.Quick,
                Iteration = mode == "tree" ? Seqeron.Genomics.Infrastructure.ClustalWIteration.Tree : Seqeron.Genomics.Infrastructure.ClustalWIteration.Alignment,
            });
            Assert.That(tool.AlignedSequences, Is.EqualTo(lib.AlignedSequences), mode);
            Assert.That(tool.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(seqs));
        }

        Assert.That(plain.AlignedSequences[3], Is.EqualTo("-----CATACTTGACT----"), "this set is the one on which the program's rows (iteration alignment: set 350) keep the plain alignment");
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(seqs, iteration: "ALIGNMENTS"));
        Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignClustalW(seqs, iteration: "tree", numIterations: 0));
        Assert.That(AlignmentTools.MultipleAlignClustalW(seqs, iteration: "none").AlignedSequences, Is.EqualTo(AlignmentTools.MultipleAlignClustalW(seqs).AlignedSequences));
    }
}
