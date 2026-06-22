using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.EpigeneticsAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// EPIGEN-* mutation killers: exact-value and inclusive-boundary tests pinning rules whose
/// canonical tests only used range assertions, leaving relational / arithmetic mutants alive.
///
/// Evidence:
///  - ChromHMM present/absent binarization at an inclusive threshold (Ernst &amp; Kellis 2012).
///  - Roadmap Epigenomics single-mark → state mapping.
///  - Allele-specific-methylation imprinting score and DMR call (Barlow &amp; Bartolomei 2014).
///  - CpG / CHG / CHH context window definition (Krueger &amp; Andrews 2011, Bismark).
/// </summary>
[TestFixture]
public class EpigeneticsAnalyzer_MutationKillers_Tests
{
    private const double Tol = 1e-9;

    #region PredictChromatinState — inclusive presence threshold (signal == threshold)

    [Test]
    public void PredictChromatinState_K4me1ExactlyAtThreshold_IsWeakEnhancer()
    {
        // Signal exactly equal to the default 0.5 presence call must count as PRESENT (>=).
        var state = PredictChromatinState(0, 0.5, 0, 0, 0, 0);
        Assert.That(state, Is.EqualTo(ChromatinState.WeakEnhancer)); // a '>' mutant ⇒ LowSignal
    }

    [Test]
    public void PredictChromatinState_K36me3ExactlyAtThreshold_IsTranscribed()
    {
        Assert.That(PredictChromatinState(0, 0, 0, 0.5, 0, 0), Is.EqualTo(ChromatinState.Transcribed));
    }

    [Test]
    public void PredictChromatinState_K27me3ExactlyAtThreshold_IsRepressed()
    {
        Assert.That(PredictChromatinState(0, 0, 0, 0, 0.5, 0), Is.EqualTo(ChromatinState.Repressed));
    }

    [Test]
    public void PredictChromatinState_K9me3ExactlyAtThreshold_IsHeterochromatin()
    {
        Assert.That(PredictChromatinState(0, 0, 0, 0, 0, 0.5), Is.EqualTo(ChromatinState.Heterochromatin));
    }

    [Test]
    public void AnnotateHistoneModifications_MarkExactlyAtThreshold_IsPresent()
    {
        // signal == threshold is present (guard is strict 'signal < threshold' ⇒ LowSignal only below).
        var mods = AnnotateHistoneModifications(new[] { (0, 10, "H3K4me3", 0.5) }).ToList();
        Assert.That(mods, Has.Count.EqualTo(1));
        Assert.That(mods[0].PredictedState, Is.EqualTo(ChromatinState.ActivePromoter));
    }

    #endregion

    #region PredictImprintedGenes — score formula & inclusive call boundaries

    [Test]
    public void PredictImprintedGenes_DifferenceExactlyAtCutoff_IsReported()
    {
        // |maternal − paternal| == minDifference (0.4): the inclusive '>=' reports it.
        var genes = PredictImprintedGenes(new[] { ("g", 0, 10, 0.5, 0.1) }).ToList();
        Assert.That(genes, Has.Count.EqualTo(1)); // a '>' mutant would drop it
    }

    [Test]
    public void PredictImprintedGenes_ExactImprintingScoreAndOrigin()
    {
        // maternal 0.9, paternal 0.1: diff = 0.8; score = diff/(m+p+0.01) = 0.8/1.01.
        var gene = PredictImprintedGenes(new[] { ("g", 0, 10, 0.9, 0.1) }).Single();
        Assert.That(gene.ImprintingScore, Is.EqualTo(0.8 / 1.01).Within(Tol)); // pins the score formula
        Assert.That(gene.ParentalOrigin, Is.EqualTo("Maternal"));
        Assert.That(gene.HasDMR, Is.True); // 0.8 > 0.5
    }

    [Test]
    public void PredictImprintedGenes_DmrThresholdIsStrict()
    {
        // diff exactly 0.5: HasDMR uses strict '> 0.5' ⇒ false. A '>=' mutant flips it to true.
        var gene = PredictImprintedGenes(new[] { ("g", 0, 10, 0.6, 0.1) }).Single();
        Assert.That(gene.HasDMR, Is.False);
        Assert.That(gene.ParentalOrigin, Is.EqualTo("Maternal"));
    }

    #endregion

    #region GetMethylationContext — index guards

    [Test]
    public void GetMethylationContext_IndexEqualsLength_ReturnsNull()
    {
        // index == length must be rejected by the index>=length guard (a '>' mutant would
        // read past the end and throw).
        Assert.That(GetMethylationContext("ACGT", 4), Is.Null);
    }

    [Test]
    public void GetMethylationContext_TerminalCytosineNoDownstream_ReturnsNull()
    {
        // Last base is C with no downstream base: the 'index+1 >= length' guard returns null
        // (a '>' or 'index-1' mutant would read sequence[index+1] out of range).
        Assert.That(GetMethylationContext("AAAC", 3), Is.Null);
    }

    #endregion
}
