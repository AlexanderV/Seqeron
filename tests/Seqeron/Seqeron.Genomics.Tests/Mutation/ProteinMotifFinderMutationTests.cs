using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// PROTMOTIF-* mutation killers: exact-value and boundary tests pinning published formulas
/// and parse rules the canonical tests only checked loosely.
///
/// Evidence:
///  - PROSITE PA-line → regex grammar incl. N/C-terminus '&lt;','&gt;' inside brackets
///    (PROSITE User Manual §IV.E; De Castro et al. 2006).
///  - Coiled-coil heptad a/d occupancy (Lupas 1991; Mason &amp; Arndt 2004).
///  - Disorder propensity sliding window normalised to [0,1].
///  - von Heijne (1986) signal-peptide weight-matrix acceptance threshold.
///  - Kyte &amp; Doolittle (1982) transmembrane hydropathy window.
///  - Information content IC = Σ log₂(20/allowed) (Schneider &amp; Stephens 1990).
/// </summary>
[TestFixture]
public class ProteinMotifFinderMutationTests
{
    private const double Tol = 1e-9;

    #region ConvertPrositeToRegex — terminus inside character class

    [Test]
    public void ConvertPrositeToRegex_NTerminusInsideBracket_EmitsAnchoredAlternation()
    {
        // PROSITE '[<G]' = "G at the N-terminus OR any G": (?:^|G).
        Assert.That(ConvertPrositeToRegex("[<G]"), Is.EqualTo("(?:^|G)"));
    }

    [Test]
    public void ConvertPrositeToRegex_CTerminusOnlyInsideBracket_EmitsEndAnchor()
    {
        // '[>]' with no letters = pure C-terminus anchor '$'.
        Assert.That(ConvertPrositeToRegex("[>]"), Is.EqualTo("$"));
    }

    [Test]
    public void ConvertPrositeToRegex_NTerminusOnlyInsideBracket_EmitsStartAnchor()
    {
        // '[<]' with no letters = pure N-terminus anchor '^'.
        Assert.That(ConvertPrositeToRegex("[<]"), Is.EqualTo("^"));
    }

    [Test]
    public void ConvertPrositeToRegex_StandardPattern_ExactRegex()
    {
        Assert.That(ConvertPrositeToRegex("N-{P}-[ST]-{P}"), Is.EqualTo("N[^P][ST][^P]"));
        Assert.That(ConvertPrositeToRegex("[ST]-x(2)-[DE]"), Is.EqualTo("[ST].{2}[DE]"));
        Assert.That(ConvertPrositeToRegex("[AG]-x(4)-G-K-[ST]"), Is.EqualTo("[AG].{4}GK[ST]"));
    }

    #endregion

    #region PredictCoiledCoils — heptad a/d occupancy & minimum length boundary

    [Test]
    public void PredictCoiledCoils_PerfectHeptad_ExactRegionAtMinLength()
    {
        // Heptad (abcdefg) with I at a(0) and d(3): "IAAIAAA" ×3 = 21 residues.
        // windowSize 21 (== MinCoiledCoilRegion 3 heptads): one window, occupancy 1.0,
        // region exactly 21 residues ⇒ kept by the strict '< MinCoiledCoilRegion' length guard.
        string seq = string.Concat(Enumerable.Repeat("IAAIAAA", 3)); // 21 nt
        var regions = PredictCoiledCoils(seq, windowSize: 21, threshold: 0.5).ToList();

        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(0));
        Assert.That(regions[0].End, Is.EqualTo(20));        // 0 + 21 - 1
        Assert.That(regions[0].Score, Is.EqualTo(1.0).Within(Tol));
    }

    #endregion

    #region PredictDisorderedRegions — region end coordinate (i-1 + window/2)

    [Test]
    public void PredictDisorderedRegions_ProlineThenOrdered_ExactRegionEnd()
    {
        // 40×P (disorder propensity 0.41 ⇒ normalized 0.91) then 30×I (ordered).
        // window 21, threshold 0.5: disordered windows are i = 0..28; the region closes at
        // i = 29 with End = (29-1) + 21/2 = 38, peak score = 0.91.
        string seq = new string('P', 40) + new string('I', 30);
        var regions = PredictDisorderedRegions(seq, windowSize: 21, threshold: 0.5).ToList();

        Assert.That(regions, Is.Not.Empty);
        Assert.That(regions[0].Start, Is.EqualTo(0));
        Assert.That(regions[0].End, Is.EqualTo(38));     // kills i+1 and window*2 mutants
        Assert.That(regions[0].Score, Is.EqualTo(0.91).Within(1e-9));
    }

    #endregion

    #region PredictSignalPeptide — acceptance threshold is inclusive

    [Test]
    public void PredictSignalPeptide_MinWeightEqualToScore_IsAccepted()
    {
        // Classic eukaryotic signal peptide. Whatever its score, setting minWeight == score
        // must accept it (inclusive '>='); a '>' mutant would reject the exact-boundary case.
        const string seq = "MKWVTFISLLLLFSSAYSRGVFRR";
        var baseline = PredictSignalPeptide(seq);
        Assert.That(baseline, Is.Not.Null);

        var atThreshold = PredictSignalPeptide(seq, prokaryote: false, minWeight: baseline!.Value.Score);
        Assert.That(atThreshold!.Value.IsLikelySignalPeptide, Is.True);
    }

    #endregion

    #region PredictTransmembraneHelices — threshold inclusive & exact peak

    [Test]
    public void PredictTransmembraneHelices_ThresholdEqualToWindowMean_IsDetected()
    {
        // 25×I: every 19-residue window has mean hydropathy 4.5 (KD(I)=4.5). With threshold
        // exactly 4.5 the inclusive '>=' detects it; a '>' mutant finds nothing.
        var regions = PredictTransmembraneHelices(new string('I', 25), windowSize: 19, threshold: 4.5).ToList();
        Assert.That(regions, Is.Not.Empty);
    }

    [Test]
    public void PredictTransmembraneHelices_AllHydrophobic_ExactRegionAndPeak()
    {
        var regions = PredictTransmembraneHelices(new string('I', 25)).ToList();
        Assert.That(regions, Has.Count.EqualTo(1));
        Assert.That(regions[0].Start, Is.EqualTo(0));
        Assert.That(regions[0].End, Is.EqualTo(24));
        Assert.That(regions[0].Score, Is.EqualTo(4.5).Within(Tol));
    }

    #endregion

    #region Motif scoring — information content (bits) & quantifier parsing

    [Test]
    public void FindMotifByPattern_SingleResidues_ExactInformationContent()
    {
        // "RGD": three fully-specified positions ⇒ IC = 3·log₂(20/1) = 3·log₂(20).
        // E-value = (N − L + 1)·2^(−IC) = 5 · 20⁻³ = 5/8000.
        var match = FindMotifByPattern("AARGDAA", "RGD").Single();
        Assert.That(match.Score, Is.EqualTo(3 * Math.Log2(20)).Within(Tol));
        Assert.That(match.EValue, Is.EqualTo(5.0 / 8000.0).Within(1e-12));
    }

    [Test]
    public void FindMotifByPattern_QuantifiedClass_ExactInformationContent()
    {
        // "[RK]{2}.[ST]": two [RK] (allowed 2 each) + one '.' (20, contributes 0) + one [ST] (2).
        // IC = 2·log₂(20/2) + 0 + log₂(20/2) = 3·log₂(10). A mis-parsed {2} quantifier
        // (repeat ≠ 2) drops one [RK] term ⇒ 2·log₂(10) ≠ 3·log₂(10).
        var match = FindMotifByPattern("RKAST", "[RK]{2}.[ST]").First();
        Assert.That(match.Score, Is.EqualTo(3 * Math.Log2(10)).Within(Tol));
    }

    #endregion
}
