using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Tests for the <c>analyze_oligo</c> MCP tool (<see cref="MolToolsTools.analyze_oligo"/>),
/// a thin wrapper over <c>ProbeDesigner.AnalyzeOligo</c>.
///
/// Expected values are derived from the algorithm's own definitions, NOT from whatever the
/// wrapper currently returns:
///  - Tm (short, &lt; 14 bp): Wallace rule Tm = 2·(A+T) + 4·(G+C)
///    (ThermoConstants.WallaceAtContribution=2, WallaceGcContribution=4, WallaceMaxLength=14).
///  - Tm (&gt;= 14 bp): salt-adjusted Tm = 81.5 + 16.6·log10(0.05) + 41·gcFraction − 600/length.
///  - GC: fraction in [0,1] (SequenceExtensions.CalculateGcFractionFast).
///  - MW (Da): Σ base weights (A=331.2, C=307.2, G=347.2, T=322.2, U=308.2) − (length−1)·18.0.
///  - ε260 (M⁻¹·cm⁻¹): Σ base contributions (A=15400, C=7400, G=11500, T=8700, U=9900).
/// Source: src/Seqeron/Algorithms/Seqeron.Genomics.MolTools/ProbeDesigner.cs#L1299 (AnalyzeOligo)
/// and src/Seqeron/Algorithms/Seqeron.Genomics.Infrastructure/ThermoConstants.cs.
/// </summary>
[TestFixture]
public class AnalyzeOligoTests
{
    [Test]
    public void AnalyzeOligo_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => MolToolsTools.analyze_oligo("ATGC"));
        Assert.Throws<ArgumentException>(() => MolToolsTools.analyze_oligo(""));
        Assert.Throws<ArgumentException>(() => MolToolsTools.analyze_oligo(null!));
    }

    [Test]
    public void AnalyzeOligo_Binding_InvokesSuccessfully()
    {
        // "ATGC": length 4 (< 14) -> Wallace rule.
        //   Tm  = 2*(A+T) + 4*(G+C) = 2*2 + 4*2 = 12
        //   GC  = 2/4 = 0.5
        //   MW  = 331.2 + 322.2 + 347.2 + 307.2 - 3*18 = 1253.8
        //   eps = 15400 + 8700 + 11500 + 7400 = 43000
        var r = MolToolsTools.analyze_oligo("ATGC");
        Assert.Multiple(() =>
        {
            Assert.That(r.Tm, Is.EqualTo(12.0).Within(1e-9));
            Assert.That(r.GcContent, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(r.MolecularWeight, Is.EqualTo(1253.8).Within(1e-6));
            Assert.That(r.ExtinctionCoefficient, Is.EqualTo(43000.0).Within(1e-9));
        });
    }

    [Test]
    public void AnalyzeOligo_20mer_UsesSaltAdjustedTm()
    {
        // "ACGTACGTACGTACGTACGT": length 20 (>= 14) -> salt-adjusted Tm.
        //   GC  = 0.5
        //   Tm  = 81.5 + 16.6*log10(0.05) + 41*0.5 - 600/20 = 50.40290207...
        //   MW  = 5*(331.2+307.2+347.2+322.2) - 19*18 = 6197.0
        //   eps = 5*(15400+7400+11500+8700) = 215000
        var r = MolToolsTools.analyze_oligo("ACGTACGTACGTACGTACGT");
        Assert.Multiple(() =>
        {
            Assert.That(r.GcContent, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(r.Tm, Is.EqualTo(50.402902071977906).Within(1e-6));
            Assert.That(r.MolecularWeight, Is.EqualTo(6197.0).Within(1e-6));
            Assert.That(r.ExtinctionCoefficient, Is.EqualTo(215000.0).Within(1e-9));
        });
    }

    [Test]
    public void AnalyzeOligo_IsCaseInsensitive()
    {
        var upper = MolToolsTools.analyze_oligo("ATGC");
        var lower = MolToolsTools.analyze_oligo("atgc");
        Assert.Multiple(() =>
        {
            Assert.That(lower.Tm, Is.EqualTo(upper.Tm).Within(1e-9));
            Assert.That(lower.GcContent, Is.EqualTo(upper.GcContent).Within(1e-9));
            Assert.That(lower.MolecularWeight, Is.EqualTo(upper.MolecularWeight).Within(1e-9));
            Assert.That(lower.ExtinctionCoefficient, Is.EqualTo(upper.ExtinctionCoefficient).Within(1e-9));
        });
    }
}
