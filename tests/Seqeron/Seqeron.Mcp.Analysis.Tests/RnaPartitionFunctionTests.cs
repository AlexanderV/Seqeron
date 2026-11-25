using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>rna_partition_function</c> MCP tool and the <c>temperature</c> parameter of the MFE tools (RNA-PARTITION-001).
/// Expected values: ViennaRNA 2.7.2 <c>fc.pf()</c>/<c>bpp()</c>/<c>pr_structure()</c>/<c>mean_bp_distance()</c>/<c>mfe()</c> (dangles = 2,
/// NNDB-precision Turner 2004 parameter file), NOT the wrapper output.
/// </summary>
[TestFixture]
public class RnaPartitionFunctionTests
{
    [Test]
    public void RnaPartitionFunction_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.RnaPartitionFunction("GGGAAAACCC"));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaPartitionFunction(""));
        Assert.Throws<ArgumentException>(() => AnalysisTools.RnaPartitionFunction(null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.RnaPartitionFunction("GGGAAAACCC", temperature: 0));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.RnaPartitionFunction("GGGAAAACCC", minProbability: 1.5));
    }

    [Test]
    public void RnaPartitionFunction_Binding_InvokesSuccessfully()
    {
        var r = AnalysisTools.RnaPartitionFunction("gggaaaaccc");
        Assert.Multiple(() =>
        {
            Assert.That(r.MfeStructure, Is.EqualTo("(((....)))"));
            Assert.That(r.MfeFreeEnergy, Is.EqualTo(-2.02).Within(1e-9));
            Assert.That(r.EnsembleFreeEnergy, Is.EqualTo(-2.1529).Within(0.005));           // fc.pf()
            Assert.That(r.MfeProbability, Is.EqualTo(0.8060).Within(0.005));                // pr_structure
            Assert.That(r.EnsembleDiversity, Is.EqualTo(1.3195).Within(0.01));              // mean_bp_distance
            Assert.That(r.CentroidStructure, Is.EqualTo("(((....)))"));
            Assert.That(r.UnpairedProbabilities, Has.Length.EqualTo(10));
            Assert.That(r.UnpairedProbabilities[0], Is.EqualTo(0.0593).Within(0.005));      // 1 - sum bpp
            Assert.That(r.PartitionFunction, Is.EqualTo(Math.Exp(-r.EnsembleFreeEnergy / (1.98717e-3 * 310.15))).Within(1e-9));
            var p18 = r.BasePairProbabilities.Single(p => p.Position1 == 1 && p.Position2 == 8);
            Assert.That(p18.Probability, Is.EqualTo(0.8376).Within(0.005));                 // bpp
            Assert.That(r.Temperature, Is.EqualTo(310.15));
        });
    }

    [Test]
    public void RnaPartitionFunction_MinProbability_FiltersListedPairs()
    {
        var all = AnalysisTools.RnaPartitionFunction("GGGAAAACCC", minProbability: 0);
        var filtered = AnalysisTools.RnaPartitionFunction("GGGAAAACCC", minProbability: 0.5);
        Assert.Multiple(() =>
        {
            Assert.That(all.BasePairProbabilities, Has.Length.EqualTo(9));
            Assert.That(filtered.BasePairProbabilities, Has.Length.EqualTo(3));
            Assert.That(filtered.BasePairProbabilities.All(p => p.Probability >= 0.5));
        });
    }

    [Test]
    public void TemperatureParameter_RescalesMfeTools()
    {
        // ViennaRNA fc.mfe() / eval_structure at 20 °C (NNDB-precision file): -3.39 (int-truncated parameters; ours unrounded)
        Assert.Multiple(() =>
        {
            Assert.That(AnalysisTools.MinimumFreeEnergy("GGGAAAACCC", 3, 293.15).Mfe, Is.EqualTo(-3.39).Within(0.02));
            Assert.That(AnalysisTools.RnaMfeStructure("GGGAAAACCC", 3, 293.15).MinimumFreeEnergy, Is.EqualTo(-3.39).Within(0.02));
            Assert.That(AnalysisTools.RnaStructureEnergy("GGGAAAACCC", "(((....)))", 293.15).Energy, Is.EqualTo(-3.39).Within(0.02));
            Assert.That(AnalysisTools.MinimumFreeEnergy("GGGAAAACCC", 3, 310.15).Mfe, Is.EqualTo(-2.02).Within(1e-9));
            Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MinimumFreeEnergy("GGGAAAACCC", 3, -1));
        });
    }
}
