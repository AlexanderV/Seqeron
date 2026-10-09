// SEQ-GC-ANALYSIS-001 — Comprehensive GC Analysis
// Evidence: docs/Evidence/SEQ-GC-ANALYSIS-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-GC-ANALYSIS-001.md
// Source: Lobry JR (1996) Mol Biol Evol 13(5):660-665; Madigan & Martinko (2003) Brock Biology of Microorganisms;
//         Biopython Bio.SeqUtils v1.84 GC_skew; Cuemath Population Variance.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GcSkewCalculator_AnalyzeGcContent_Tests
{
    #region AnalyzeGcContent(DnaSequence) — Canonical

    // M1/M2/M3 — overall scalar metrics on "GGGCCAT" (G=3,C=2,A=1,T=1,n=7).
    // GC% = (3+2)/7*100 = 71.4285...; GC skew = (3-2)/5 = 0.2; AT skew = (1-1)/2 = 0.0.
    // An off-by-one or wrong formula (e.g. fraction not %, or (G-C)/n) would not produce these exact values.
    [Test]
    public void AnalyzeGcContent_KnownSequence_ComputesOverallContentSkewAndAtSkew()
    {
        var seq = new DnaSequence("GGGCCAT");

        var result = GcSkewCalculator.AnalyzeGcContent(seq);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverallGcContent, Is.EqualTo(5.0 / 7.0 * 100.0).Within(1e-10),
                "GC% = (G+C)/(A+T+G+C)*100 = (3+2)/7*100 = 71.42857142857143 (Brock/Wikipedia GC-content)");
            Assert.That(result.OverallGcSkew, Is.EqualTo(0.2).Within(1e-10),
                "GC skew = (G-C)/(G+C) = (3-2)/(3+2) = 0.2 (Lobry 1996 / Biopython GC_skew)");
            Assert.That(result.OverallAtSkew, Is.EqualTo(0.0).Within(1e-10),
                "AT skew = (A-T)/(A+T) = (1-1)/2 = 0.0 (balanced A/T)");
            Assert.That(result.SequenceLength, Is.EqualTo(7),
                "SequenceLength equals the input length (INV-05)");
        });
    }

    // M4/M5 — windowed population variance on "GGCC" w=2 step=2 -> windows GG(+1,100%), CC(-1,100%).
    // GcSkewVariance = ((1-0)^2+(-1-0)^2)/2 = 1.0; GcContentVariance = ((100-100)^2+(100-100)^2)/2 = 0.0.
    // Sample variance (/N-1) would give 2.0, not 1.0 — this pins the population estimator.
    [Test]
    public void AnalyzeGcContent_TwoWindows_ComputesPopulationVariance()
    {
        var seq = new DnaSequence("GGCC");

        var result = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 2, stepSize: 2);

        Assert.Multiple(() =>
        {
            Assert.That(result.WindowedGcSkew.Count, Is.EqualTo(2), "Two non-overlapping windows GG and CC");
            Assert.That(result.GcSkewVariance, Is.EqualTo(1.0).Within(1e-10),
                "Population variance of {+1,-1}: mean 0, (1+1)/2 = 1.0 (Cuemath population variance, /N not /N-1)");
            Assert.That(result.GcContentVariance, Is.EqualTo(0.0).Within(1e-10),
                "Both windows are 100% GC: variance of {100,100} = 0.0");
        });
    }

    // M6 — windowed GC% population variance on "AAAGGGCCCTTT" w=3 step=3 -> windows 0,100,100,0.
    // mean = 50; var = ((0-50)^2*2 + (100-50)^2*2)/4 = (2500*4)/4 = 2500.
    [Test]
    public void AnalyzeGcContent_FourWindows_ComputesGcContentVariance()
    {
        var seq = new DnaSequence("AAAGGGCCCTTT");

        var result = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 3, stepSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(result.WindowedGcContent.Select(w => w.GcContent).ToArray(),
                Is.EqualTo(new[] { 0.0, 100.0, 100.0, 0.0 }).Within(1e-10),
                "Windows AAA,GGG,CCC,TTT have GC% 0,100,100,0");
            Assert.That(result.GcContentVariance, Is.EqualTo(2500.0).Within(1e-10),
                "Population variance of {0,100,100,0}: mean 50, (2500*4)/4 = 2500 (Cuemath /N)");
        });
    }

    // M7 — windowing geometry on "ACGTACGTAC" w=4 step=2 -> 4 windows (starts 0,2,4,6).
    // First window: start=0, end=3 (inclusive), position = start + w/2 = 2 (Biopython multi-window; INV-05).
    [Test]
    public void AnalyzeGcContent_SlidingWindows_HasExactCountAndBoundaries()
    {
        var seq = new DnaSequence("ACGTACGTAC");

        var result = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 4, stepSize: 2);

        Assert.Multiple(() =>
        {
            Assert.That(result.WindowedGcSkew.Count, Is.EqualTo(4),
                "floor((10-4)/2)+1 = 4 full windows (INV-05)");
            Assert.That(result.WindowedGcContent.Count, Is.EqualTo(4),
                "GC-content profile has the same window count as the skew profile");
            Assert.That(result.WindowedGcContent[0].WindowStart, Is.EqualTo(0), "First window starts at index 0");
            Assert.That(result.WindowedGcContent[0].WindowEnd, Is.EqualTo(3), "First window end is inclusive: 0+4-1 = 3");
            Assert.That(result.WindowedGcContent[0].Position, Is.EqualTo(2), "Window position is midpoint start + w/2 = 2");
        });
    }

    // M8 — pure-G upper bound: GC skew = (4-0)/4 = +1; GC% = 100 (skew range bound, Wikipedia GC skew).
    [Test]
    public void AnalyzeGcContent_PureGuanine_SkewIsPlusOneContentIsHundred()
    {
        var seq = new DnaSequence("GGGG");

        var result = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 4, stepSize: 4);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverallGcSkew, Is.EqualTo(1.0).Within(1e-10),
                "Pure-G: (G-C)/(G+C) = (4-0)/4 = +1 (upper bound)");
            Assert.That(result.OverallGcContent, Is.EqualTo(100.0).Within(1e-10),
                "All bases are G: GC% = 4/4*100 = 100");
        });
    }

    // M8b — pure-C lower bound: GC skew = (0-4)/4 = -1; GC% = 100 (skew range lower bound, Wikipedia GC skew:
    // "-1 corresponds to G = 0"). Complements M8 so both documented bounds (+1 and -1) are pinned directly.
    [Test]
    public void AnalyzeGcContent_PureCytosine_SkewIsMinusOneContentIsHundred()
    {
        var seq = new DnaSequence("CCCC");

        var result = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 4, stepSize: 4);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverallGcSkew, Is.EqualTo(-1.0).Within(1e-10),
                "Pure-C: (G-C)/(G+C) = (0-4)/4 = -1 (lower bound; Wikipedia GC skew '-1 corresponds to G = 0')");
            Assert.That(result.OverallGcContent, Is.EqualTo(100.0).Within(1e-10),
                "All bases are C: GC% = 4/4*100 = 100");
        });
    }

    // M9 — no G/C: skew zero-division -> 0; GC% numerator 0 -> 0; AT skew of ATATAT = 0 (balanced).
    [Test]
    public void AnalyzeGcContent_NoGcBases_SkewAndContentAreZero()
    {
        var seq = new DnaSequence("ATATAT");

        var result = GcSkewCalculator.AnalyzeGcContent(seq);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverallGcSkew, Is.EqualTo(0.0).Within(1e-10),
                "G+C = 0 => GC skew defined as 0 (Biopython zero-division handling)");
            Assert.That(result.OverallGcContent, Is.EqualTo(0.0).Within(1e-10),
                "No G/C => GC% = 0");
            Assert.That(result.OverallAtSkew, Is.EqualTo(0.0).Within(1e-10),
                "A=3,T=3 => AT skew (3-3)/6 = 0");
        });
    }

    // M3b — asymmetric AT skew pins the (A-T)/(A+T) formula and its sign: "AAAT" has A=3,T=1 =>
    // AT skew = (3-1)/(3+1) = +0.5. A balanced-only test (M3) cannot distinguish (A-T)/(A+T) from
    // (T-A)/(A+T) or (A-T)/n; this asymmetric case does.
    [Test]
    public void AnalyzeGcContent_AdenineRich_AtSkewIsPositiveHalf()
    {
        var seq = new DnaSequence("AAAT");

        var result = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 4, stepSize: 4);

        Assert.That(result.OverallAtSkew, Is.EqualTo(0.5).Within(1e-10),
            "AT skew = (A-T)/(A+T) = (3-1)/(3+1) = +0.5 (Charneski et al. 2011 / Lobry 1996)");
    }

    #endregion

    #region Edge Cases and Failure Modes

    // S1 — sequence shorter than window: no full window => empty profiles, variances 0; scalars still computed.
    [Test]
    public void AnalyzeGcContent_SequenceShorterThanWindow_EmptyWindowsButScalarsComputed()
    {
        var seq = new DnaSequence("ACGT");

        var result = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 10, stepSize: 1);

        Assert.Multiple(() =>
        {
            Assert.That(result.WindowedGcSkew, Is.Empty, "No full window fits => empty skew profile");
            Assert.That(result.WindowedGcContent, Is.Empty, "No full window fits => empty content profile");
            Assert.That(result.GcSkewVariance, Is.EqualTo(0.0).Within(1e-10), "No windows => variance 0");
            Assert.That(result.GcContentVariance, Is.EqualTo(0.0).Within(1e-10), "No windows => variance 0");
            Assert.That(result.OverallGcContent, Is.EqualTo(50.0).Within(1e-10),
                "Whole-sequence GC% = (G+C)/n*100 = 2/4*100 = 50 (ACGT)");
            Assert.That(result.SequenceLength, Is.EqualTo(4), "Length is still the whole sequence");
        });
    }

    // S2 — null DnaSequence throws ArgumentNullException (parity with sibling methods).
    [Test]
    public void AnalyzeGcContent_NullDnaSequence_Throws()
    {
        Assert.That(() => GcSkewCalculator.AnalyzeGcContent((DnaSequence)null!),
            NUnit.Framework.Throws.ArgumentNullException, "Null DnaSequence must throw ArgumentNullException");
    }

    // S3 — null/empty string overload returns a zero result with empty profiles and length 0.
    [Test]
    public void AnalyzeGcContent_NullOrEmptyString_ReturnsZeroResult()
    {
        var fromNull = GcSkewCalculator.AnalyzeGcContent((string)null!);
        var fromEmpty = GcSkewCalculator.AnalyzeGcContent("");

        Assert.Multiple(() =>
        {
            Assert.That(fromNull.SequenceLength, Is.EqualTo(0), "Null string => length 0");
            Assert.That(fromNull.WindowedGcSkew, Is.Empty, "Null string => no windows");
            Assert.That(fromNull.OverallGcContent, Is.EqualTo(0.0).Within(1e-10), "Null string => GC% 0");
            Assert.That(fromEmpty.SequenceLength, Is.EqualTo(0), "Empty string => length 0");
            Assert.That(fromEmpty.WindowedGcContent, Is.Empty, "Empty string => no windows");
        });
    }

    #endregion

    #region AnalyzeGcContent(string) — Delegate (overload equivalence)

    // C1 — string and DnaSequence overloads produce identical results for the same sequence.
    [Test]
    public void AnalyzeGcContent_StringOverload_MatchesDnaSequenceOverload()
    {
        const string s = "AAAGGGCCCTTT";
        var fromString = GcSkewCalculator.AnalyzeGcContent(s, windowSize: 3, stepSize: 3);
        var fromDna = GcSkewCalculator.AnalyzeGcContent(new DnaSequence(s), windowSize: 3, stepSize: 3);

        Assert.Multiple(() =>
        {
            Assert.That(fromString.OverallGcContent, Is.EqualTo(fromDna.OverallGcContent).Within(1e-10),
                "Both overloads compute the same overall GC content");
            Assert.That(fromString.OverallGcSkew, Is.EqualTo(fromDna.OverallGcSkew).Within(1e-10),
                "Both overloads compute the same overall GC skew");
            Assert.That(fromString.OverallAtSkew, Is.EqualTo(fromDna.OverallAtSkew).Within(1e-10),
                "Both overloads compute the same overall AT skew");
            Assert.That(fromString.GcContentVariance, Is.EqualTo(fromDna.GcContentVariance).Within(1e-10),
                "Both overloads compute the same GC-content variance");
            Assert.That(fromString.GcSkewVariance, Is.EqualTo(fromDna.GcSkewVariance).Within(1e-10),
                "Both overloads compute the same GC-skew variance");
            Assert.That(fromString.WindowedGcSkew.Count, Is.EqualTo(fromDna.WindowedGcSkew.Count),
                "Both overloads produce the same number of windows");
        });
    }

    // C1b — string overload lowercases input: "gggccat" must match uppercase "GGGCCAT".
    [Test]
    public void AnalyzeGcContent_StringOverload_IsCaseInsensitive()
    {
        var lower = GcSkewCalculator.AnalyzeGcContent("gggccat");

        Assert.Multiple(() =>
        {
            Assert.That(lower.OverallGcContent, Is.EqualTo(5.0 / 7.0 * 100.0).Within(1e-10),
                "Lowercase input is uppercased: GC% = 71.42857142857143");
            Assert.That(lower.OverallGcSkew, Is.EqualTo(0.2).Within(1e-10),
                "Lowercase input is uppercased: GC skew = 0.2");
        });
    }

    #endregion

    #region Review 2026-09 — RNA U in the GC denominator, window/step guards, variance canonical

    // F1 — GC content counts U (RNA) in the denominator, like T. Source: GC-content is defined for
    // "a DNA or RNA molecule … adenine and uracil in RNA" (Wikipedia "GC-content"); reference
    // Biopython 1.88 gc_fraction (denominator counts "ATWUatwu"), docstring example
    // gc_fraction("GGAUCUUCGGAUCU") = 0.50. Before the fix U was dropped: 7/9·100 = 77.78.
    [Test]
    public void AnalyzeGcContent_RnaString_CountsUracilInGcDenominator_MatchesBiopython()
    {
        var result = GcSkewCalculator.AnalyzeGcContent("GGAUCUUCGGAUCU", windowSize: 7, stepSize: 7);
        var frac = GcSkewCalculator.AnalyzeGcContent("ACGU", windowSize: 4, stepSize: 4, fraction: true);

        Assert.Multiple(() =>
        {
            Assert.That(result.OverallGcContent, Is.EqualTo(50.0).Within(1e-10),
                "Biopython gc_fraction('GGAUCUUCGGAUCU')*100 = 50.0");
            Assert.That(result.WindowedGcContent.Select(w => w.GcContent).ToArray(),
                Is.EqualTo(new[] { 42.857142857142854, 57.14285714285714 }).Within(1e-10),
                "Biopython gc_fraction of windows GGAUCUU / CGGAUCU ×100 = 3/7, 4/7");
            Assert.That(result.GcContentVariance, Is.EqualTo(51.0204081632653).Within(1e-9),
                "numpy.var([300/7, 400/7]) (ddof=0) = 51.0204081632653");
            Assert.That(frac.OverallGcContent, Is.EqualTo(0.5).Within(1e-12),
                "Biopython gc_fraction('ACGU') = 0.5");
        });
    }

    // F2 — windowSize/stepSize < 1 are rejected eagerly on both overloads (a zero step never
    // terminated; window 0 emitted empty windows). Biopython GC_skew(seq, 0) raises ValueError.
    [TestCase(0, 1)]
    [TestCase(-1, 1)]
    [TestCase(4, 0)]
    [TestCase(4, -2)]
    public void AnalyzeGcContent_WindowOrStepBelowOne_Throws(int window, int step)
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => GcSkewCalculator.AnalyzeGcContent(new DnaSequence("ACGTACGT"), window, step),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => GcSkewCalculator.AnalyzeGcContent("ACGTACGT", window, step),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => GcSkewCalculator.AnalyzeGcContent("", window, step),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>(),
                "validated eagerly, before the null/empty short-circuit");
        });
    }

    // Canonical population variance: Cuemath worked example {12,13,12,14,19} -> 6.8 = numpy.var(ddof=0).
    [Test]
    public void StatisticsHelper_PopulationVariance_MatchesNumpyVar()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Seqeron.Genomics.Infrastructure.StatisticsHelper.PopulationVariance(
                new[] { 12.0, 13.0, 12.0, 14.0, 19.0 }), Is.EqualTo(6.8).Within(1e-12));
            Assert.That(Seqeron.Genomics.Infrastructure.StatisticsHelper.PopulationVariance(
                Array.Empty<double>()), Is.EqualTo(0.0));
        });
    }

    #endregion

    #region CalculateWindowedGcContent — public sliding-GC driver (review 2026-09 B03 R18, F27)

    // M13 — per-window GC fraction locked to Biopython 1.88:
    //   s = "ATGCGCGATTACGGCCATatgcgt"; [gc_fraction(s[i:i+5]) for i in range(0, len(s)-4, 3)]
    //   = [0.6, 0.8, 0.2, 0.6, 0.8, 0.2, 0.6]  (lower case folded; trailing partial window dropped).
    // Position = start + 5/2, WindowEnd = start + 4 (inclusive).
    [Test]
    public void CalculateWindowedGcContent_MatchesBiopythonPerWindowGcFraction()
    {
        const string s = "ATGCGCGATTACGGCCATatgcgt";
        double[] expected = { 0.6, 0.8, 0.2, 0.6, 0.8, 0.2, 0.6 };

        var fractions = GcSkewCalculator.CalculateWindowedGcContent(s, 5, 3, fraction: true).ToList();
        var percents = GcSkewCalculator.CalculateWindowedGcContent(s, 5, 3).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(fractions.Select(p => p.GcContent), Is.EqualTo(expected).Within(1e-12));
            Assert.That(percents.Select(p => p.GcContent), Is.EqualTo(expected.Select(f => f * 100.0)).Within(1e-10));
            Assert.That(fractions.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 3, 6, 9, 12, 15, 18 }));
            Assert.That(fractions.Select(p => p.Position), Is.EqualTo(new[] { 2, 5, 8, 11, 14, 17, 20 }));
            Assert.That(fractions.Select(p => p.WindowEnd), Is.EqualTo(new[] { 4, 7, 10, 13, 16, 19, 22 }));
        });
    }

    // M13 (cont.) — it is the driver behind AnalyzeGcContent: identical windows and values; the
    // DnaSequence overload agrees with the string overload.
    [Test]
    public void CalculateWindowedGcContent_IsTheAnalyzeGcContentDriver()
    {
        const string s = "ATGCGCGATTACGGCCATATGCGT";
        var analysis = GcSkewCalculator.AnalyzeGcContent(s, 5, 3, fraction: true).WindowedGcContent;

        Assert.Multiple(() =>
        {
            Assert.That(GcSkewCalculator.CalculateWindowedGcContent(s, 5, 3, fraction: true), Is.EqualTo(analysis));
            Assert.That(GcSkewCalculator.CalculateWindowedGcContent(new DnaSequence(s), 5, 3, fraction: true),
                Is.EqualTo(analysis));
            // numpy.var(ddof=0) of the Biopython per-window fractions = 0.053877551020408164.
            Assert.That(GcSkewCalculator.AnalyzeGcContent(s, 5, 3, fraction: true).GcContentVariance,
                Is.EqualTo(0.053877551020408164).Within(1e-12));
        });
    }

    // M14 — ambiguity overload locked to Biopython 1.88 gc_fraction(window, ambiguous=…):
    //   "GGSNNBWAAC", w=4, step=2 → windows GGSN, SNNB, NBWA, WAAC
    //   remove   [1.0, 1.0, 0.0, 0.25]
    //   ignore   [0.75, 0.25, 0.0, 0.25]
    //   weighted [0.875, 0.6666666666666666, 0.29166666666666663, 0.25]
    [TestCase(SequenceExtensions.GcAmbiguityMode.Remove, new[] { 1.0, 1.0, 0.0, 0.25 })]
    [TestCase(SequenceExtensions.GcAmbiguityMode.Ignore, new[] { 0.75, 0.25, 0.0, 0.25 })]
    [TestCase(SequenceExtensions.GcAmbiguityMode.Weighted, new[] { 0.875, 0.6666666666666666, 0.29166666666666663, 0.25 })]
    public void CalculateWindowedGcContent_AmbiguityMode_MatchesBiopython(
        SequenceExtensions.GcAmbiguityMode mode, double[] expected)
    {
        var points = GcSkewCalculator.CalculateWindowedGcContent("GGSNNBWAAC", 4, 2, true, mode).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(points.Select(p => p.GcContent), Is.EqualTo(expected).Within(1e-12));
            Assert.That(points.Select(p => p.WindowStart), Is.EqualTo(new[] { 0, 2, 4, 6 }));
            Assert.That(GcSkewCalculator.CalculateWindowedGcContent("GGSNNBWAAC", 4, 2, false, mode)
                .Select(p => p.GcContent), Is.EqualTo(expected.Select(f => f * 100.0)).Within(1e-10));
        });
    }

    // M15 — guards as the sibling windowed methods: window/step < 1 throw eagerly (also for ""),
    // null DnaSequence throws, null/empty string or window > length yields no windows; a huge step
    // terminates after the first window (no int overflow).
    [Test]
    public void CalculateWindowedGcContent_Guards()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => GcSkewCalculator.CalculateWindowedGcContent("ACGT", 0, 1),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => GcSkewCalculator.CalculateWindowedGcContent("", 4, 0),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => GcSkewCalculator.CalculateWindowedGcContent("ACGT", 2, 0, true,
                    SequenceExtensions.GcAmbiguityMode.Remove),
                NUnit.Framework.Throws.InstanceOf<ArgumentOutOfRangeException>());
            Assert.That(() => GcSkewCalculator.CalculateWindowedGcContent((DnaSequence)null!),
                NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(GcSkewCalculator.CalculateWindowedGcContent((string)null!), Is.Empty);
            Assert.That(GcSkewCalculator.CalculateWindowedGcContent(""), Is.Empty);
            Assert.That(GcSkewCalculator.CalculateWindowedGcContent("ACGT", 5, 1), Is.Empty);
            Assert.That(GcSkewCalculator.CalculateWindowedGcContent("GCAT", 2, int.MaxValue, true)
                .Select(p => p.GcContent), Is.EqualTo(new[] { 1.0 }));
        });
    }

    #endregion
}
