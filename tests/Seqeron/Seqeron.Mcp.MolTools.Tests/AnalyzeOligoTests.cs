using NUnit.Framework;
using Seqeron.Mcp.MolTools.Tools;

namespace Seqeron.Mcp.MolTools.Tests;

/// <summary>
/// Tests for the <c>analyze_oligo</c> MCP tool (<see cref="MolToolsTools.analyze_oligo"/>),
/// a thin wrapper over <c>ProbeDesigner.AnalyzeOligo</c>.
///
/// Expected values come from the reference implementations, NOT from whatever the wrapper returns:
///  - Tm: primer3-py 2.3.1 calc_tm(seq, mv_conc=50, dv_conc=0, dntp_conc=0, dna_conc=50)
///    (Primer3 seqtm at the Primer3 hybridization-probe defaults; null when not computable).
///  - GC: fraction in [0,1] (SequenceExtensions.CalculateGcFractionFast).
///  - MW (Da): Biopython 1.88 Bio.SeqUtils.molecular_weight(seq, "DNA" | "RNA") (single-stranded).
///  - ε260 (M⁻¹·cm⁻¹): Σ base contributions (A=15400, C=7400, G=11500, T=8700, U=9900).
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
        // "ATGC": primer3.calc_tm(mv=50, dv=0, dntp=0, dna=50) = -54.53620887888866;
        //   GC = 0.5; Biopython molecular_weight("ATGC", "DNA") = 1253.8027;
        //   eps = 15400 + 8700 + 11500 + 7400 = 43000
        var r = MolToolsTools.analyze_oligo("ATGC");
        Assert.Multiple(() =>
        {
            Assert.That(r.Tm, Is.EqualTo(-54.53620887888866).Within(1e-9));
            Assert.That(r.GcContent, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(r.MolecularWeight, Is.EqualTo(1253.8027).Within(1e-6));
            Assert.That(r.ExtinctionCoefficient, Is.EqualTo(43000.0).Within(1e-9));
        });
    }

    [Test]
    public void AnalyzeOligo_20mer_UsesPrimer3NearestNeighborTm()
    {
        // "ACGTACGTACGTACGTACGT": GC 0.5; primer3.calc_tm(mv=50, dv=0, dntp=0, dna=50) = 53.99351583691026;
        //   Biopython molecular_weight(..., "DNA") = 6196.9523; eps = 5*(15400+7400+11500+8700) = 215000
        var r = MolToolsTools.analyze_oligo("ACGTACGTACGTACGTACGT");
        Assert.Multiple(() =>
        {
            Assert.That(r.GcContent, Is.EqualTo(0.5).Within(1e-9));
            Assert.That(r.Tm, Is.EqualTo(53.99351583691026).Within(1e-9));
            Assert.That(r.MolecularWeight, Is.EqualTo(6196.9523).Within(1e-6));
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

    [Test]
    public void AnalyzeOligo_Rna_TmNullAndRnaMolecularWeight()
    {
        // RNA: Primer3 seqtm has no U (DNA nearest-neighbour model) -> Tm null;
        // Biopython molecular_weight("ACGU", "RNA") = 1303.7737 (UMP 324.1813).
        var r = MolToolsTools.analyze_oligo("ACGU");
        Assert.Multiple(() =>
        {
            Assert.That(r.Tm, Is.Null);
            Assert.That(r.MolecularWeight, Is.EqualTo(1303.7737).Within(1e-6));
        });
    }
}
