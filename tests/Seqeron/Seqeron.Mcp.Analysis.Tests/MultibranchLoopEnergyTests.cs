using NUnit.Framework;
using Seqeron.Mcp.Analysis.Tools;

namespace Seqeron.Mcp.Analysis.Tests;

/// <summary>
/// Tests for the <c>multibranch_loop_energy</c> MCP tool.
/// Expected values: linear Turner 2004 model 9.3 + 0·unpaired − 0.9·helices (ViennaRNA 2.7.2 default.c,
/// eval_structure_verbose "Multi loop" 660 for a GC-closed 3-way loop, dangles = 0) and the NNDB efn2 model
/// (rna_2004_multibranch_loops.html §9.1), NOT the wrapper output.
/// </summary>
[TestFixture]
public class MultibranchLoopEnergyTests
{
    [Test]
    public void MultibranchLoopEnergy_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => AnalysisTools.MultibranchLoopEnergy(3, 6));
        Assert.Throws<ArgumentOutOfRangeException>(() => AnalysisTools.MultibranchLoopEnergy(2, 6));
        // Strain is an efn2 term: rejected without the per-gap counts, and must match them when given.
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultibranchLoopEnergy(3, 1, hasStrain: true));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultibranchLoopEnergy(3, 1, unpairedBetweenHelices: new[] { 1, 0, 0 }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultibranchLoopEnergy(3, 2, unpairedBetweenHelices: new[] { 1, 0, 0 }));
        Assert.Throws<ArgumentException>(() => AnalysisTools.MultibranchLoopEnergy(4, 1, unpairedBetweenHelices: new[] { 1, 0, 0 }));
    }

    [Test]
    public void MultibranchLoopEnergy_Binding_InvokesSuccessfully()
    {
        Assert.Multiple(() =>
        {
            // Linear (engine) model: 9.3 − 0.9·3 = 6.6, independent of the unpaired count (b = 0).
            Assert.That(AnalysisTools.MultibranchLoopEnergy(3, 6).Energy, Is.EqualTo(6.6).Within(1e-9));
            Assert.That(AnalysisTools.MultibranchLoopEnergy(3, 3, stackingEnergy: -4.1).Energy, Is.EqualTo(2.5).Within(1e-9));
            // NNDB efn2: gaps (1,0,0) → mean asymmetry 2/3, strained 3-way: 9.25 + 0.61 − 1.89 − 2.0 + 3.14 = 9.11.
            Assert.That(AnalysisTools.MultibranchLoopEnergy(3, 1, hasStrain: true, stackingEnergy: -2.0,
                unpairedBetweenHelices: new[] { 1, 0, 0 }).Energy, Is.EqualTo(9.11).Within(1e-9));
            // NNDB §9.3 example: 4 helices, mean asymmetry capped at 2.0 → 9.25 + 1.82 − 2.52 = 8.55 (NNDB prints 8.6).
            Assert.That(AnalysisTools.MultibranchLoopEnergy(4, 8, unpairedBetweenHelices: new[] { 2, 0, 1, 5 }).Energy,
                Is.EqualTo(8.55).Within(1e-9));
        });
    }
}
