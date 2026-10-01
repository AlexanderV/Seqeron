using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// generate_emboss_consensus — delegation to <c>MotifFinder.GenerateEmbossConsensus</c>.
/// Expected values: EMBOSS 6.6.0 cons outputs locked in MotifFinder_AlignmentConsensus_Tests (Emboss_000/003/010/013).
/// </summary>
[TestFixture]
public class GenerateEmbossConsensusTests
{
    private static readonly string[] Dna = { "ACGTAC-T", "ACGTTCAT", "AGGTAC-T", "tCGAAG-T" };
    private static readonly string[] Protein = { "MKVLAAGIVG", "MKVLSAGIVA", "MRVLAAG-VG", "MKILTAGLVG" };

    [Test]
    public void GenerateEmbossConsensus_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.GenerateEmbossConsensus(Dna));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateEmbossConsensus(null!));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateEmbossConsensus(new[] { "ACGT" }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateEmbossConsensus(Dna, "rna"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.GenerateEmbossConsensus(new[] { "ACG", "AC" }));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.GenerateEmbossConsensus(Dna, identity: -1));
    }

    [Test]
    public void GenerateEmbossConsensus_Defaults_EqualsCons()
    {
        // cons -sequence a.fa -snucleotide
        Assert.That(AnalysisTools.GenerateEmbossConsensus(Dna).Consensus, Is.EqualTo("ACGTACnT"));
        // cons -sequence a.fa -sprotein
        Assert.That(AnalysisTools.GenerateEmbossConsensus(Protein, "protein").Consensus, Is.EqualTo("MKVLAAGIVG"));
    }

    [Test]
    public void GenerateEmbossConsensus_PluralitySetcase_EqualsCons_AndEqualsLibrary()
    {
        // cons -plurality 1.0 -setcase 3.0
        Assert.That(AnalysisTools.GenerateEmbossConsensus(Dna, plurality: 1.0, setcase: 3.0).Consensus, Is.EqualTo("acGtacaT"));
        Assert.That(AnalysisTools.GenerateEmbossConsensus(Protein, "Protein", 1.0, 0, 3.0).Consensus, Is.EqualTo("MKVLaAGiVg"));
        Assert.That(AnalysisTools.GenerateEmbossConsensus(Dna, identity: 4).Consensus,
            Is.EqualTo(global::Seqeron.Genomics.Analysis.MotifFinder.GenerateEmbossConsensus(Dna, identity: 4)));
    }
}
