using NUnit.Framework;
using Seqeron.Genomics.MolTools;
using Seqeron.Mcp.MolTools.Models;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

[TestFixture]
public class ReduceSecondaryStructureTests
{
    private static CodonUsageTableInput EColi() => new(Preset: "EColiK12");

    [Test]
    public void ReduceSecondaryStructure_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.reduce_secondary_structure("ATGATG", EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.reduce_secondary_structure("", EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.reduce_secondary_structure(null!, EColi()));
        Assert.Throws<ArgumentException>(() => MolToolsTools.reduce_secondary_structure("ATGATG", EColi(), 0));
        Assert.Throws<ArgumentException>(() => MolToolsTools.reduce_secondary_structure(
            "ATGATG", EColi(), 40, (CodonOptimizer.SecondaryStructureScoring)42));
    }

    [Test]
    public void ReduceSecondaryStructure_Scoring_DefaultsToPairCount_AndExposesTurner2004Mfe()
    {
        // review 2026-09, CODON-OPT-001 F41 (B02 A1-6): the Turner 2004 window-MFE mode is
        // opt-in; the default stays the canonical-pair-count heuristic.
        const string gene = "AUGAAACUGGCAGAUCUGGCGCAGCUGGCAUUUCGUGAAGAUCUGCGUCAGGCGCUGGAUCGU"
            + "UUCAUGGCGCAGCUGGAUAAAGCGCUGGCGGAAUAA";

        var byDefault = MolToolsTools.reduce_secondary_structure(gene, EColi());
        var pairCount = MolToolsTools.reduce_secondary_structure(
            gene, EColi(), 40, CodonOptimizer.SecondaryStructureScoring.PairCount);
        var turner = MolToolsTools.reduce_secondary_structure(
            gene, EColi(), 40, CodonOptimizer.SecondaryStructureScoring.Turner2004Mfe);

        Assert.Multiple(() =>
        {
            Assert.That(byDefault.OptimizedSequence, Is.EqualTo(pairCount.OptimizedSequence));
            Assert.That(byDefault.OptimizedSequence, Is.EqualTo(gene), "no window exceeds the 0.5 pair fraction");
            Assert.That(turner.OptimizedSequence, Is.EqualTo(
                "AUGAAGCUGGCAGAUCUUGCUCAACUUGCAUUUCGUGAAGACCUUCGGCAAGCUUUAGAUCGU"
                + "UUCAUGGCACAAUUAGAUAAAGCAUUAGCAGAAUAA"),
                "Turner 2004 MFE scoring rewrites 19 codons (locked in CodonOptimizer_SecondaryStructureScoring_Tests)");
            Assert.That(turner.OptimizedSequence.Length, Is.EqualTo(gene.Length));
        });
    }

    [Test]
    public void ReduceSecondaryStructure_Binding_InvokesSuccessfully()
    {
        // A sequence shorter than the window keeps its codons, but the output is always the
        // normalised RNA form (review 2026-09, CODON-OPT-001 F24).
        var shortSeq = MolToolsTools.reduce_secondary_structure("ATGATG", EColi(), window_size: 40);
        Assert.That(shortSeq.OptimizedSequence, Is.EqualTo("AUGAUG"));

        // A sequence at/above the window is processed; synonymous swaps preserve the codon count.
        // 15 codons (45 nt) > 40-nt window.
        string coding = string.Concat(System.Linq.Enumerable.Repeat("GCG", 15)); // 15x Ala
        var result = MolToolsTools.reduce_secondary_structure(coding, EColi(), window_size: 40);
        Assert.That(result.OptimizedSequence.Length, Is.EqualTo(coding.Length));
    }
}
