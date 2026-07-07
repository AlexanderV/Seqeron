// CONVENTIONS — opt-in Biopython/VCF compatibility modes
// Verifies the opt-in compatibility surfaces that retire docs/Validation/LIMITATIONS.md §3:
//   - SequenceExtensions.CalculateGcFraction(GcAmbiguityMode)  (Biopython gc_fraction ambiguous=...)
//   - GcSkewCalculator.AnalyzeGcContent(..., fraction: true)   (Biopython [0,1] vs default %)
//   - SequenceStatistics.CalculateGcContentProfile(..., fraction: true)
//   - Variant.VcfPosition                                      (VCF 1-based POS)
// Source (gc_fraction + _gc_values, ambiguous="remove"/"ignore"/"weighted"):
//   Biopython Bio/SeqUtils/__init__.py (master), retrieved 2026-06-23 from
//   https://raw.githubusercontent.com/biopython/biopython/master/Bio/SeqUtils/__init__.py
// Source (VCF 1-based POS): VCF v4.3 spec §1.4.1, https://samtools.github.io/hts-specs/VCFv4.3.pdf
namespace Seqeron.Genomics.Tests.Unit.Core;

[TestFixture]
public class ConventionCompatibility_OptIn_Tests
{
    private const double Tol = 1e-10;

    #region CalculateGcFraction(GcAmbiguityMode) — Biopython gc_fraction

    // Remove mode (default Biopython): gc = count(C,G,S); length = gc + count(A,T,W,U).
    // "GCAT" -> gc=2, length=4 -> 0.5. No ambiguity codes present, so identical to default.
    [Test]
    public void CalculateGcFraction_Remove_UnambiguousSequence_MatchesGcOverLength()
    {
        double result = "GCAT".CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Remove);

        Assert.That(result, Is.EqualTo(0.5).Within(Tol),
            "Remove: gc=count(C,G,S)=2, length=2+count(A,T,W,U)=4 -> 0.5 (Biopython gc_fraction).");
    }

    // Remove: S counts toward GC numerator AND length; W counts toward length only.
    // "GCSW": gc=count(C,G,S)=3, length=3+count(W)=4 -> 0.75.
    [Test]
    public void CalculateGcFraction_Remove_SCountsAsGc_WCountsAsLengthOnly()
    {
        double result = "GCSW".CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Remove);

        Assert.That(result, Is.EqualTo(0.75).Within(Tol),
            "Remove: S in numerator (G|C), W in denominator only -> gc=3, length=4 -> 0.75.");
    }

    // Remove: other ambiguity codes (N) are excluded from numerator AND denominator.
    // "GCATN": gc=2, length=4 (N dropped) -> 0.5.
    [Test]
    public void CalculateGcFraction_Remove_OtherAmbiguityExcludedFromBoth()
    {
        double result = "GCATN".CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Remove);

        Assert.That(result, Is.EqualTo(0.5).Within(Tol),
            "Remove: N excluded from both numerator and length -> gc=2, length=4 -> 0.5.");
    }

    // Ignore mode: numerator = count(C,G,S); denominator = FULL length (every char).
    // "GCATN": gc=2, length=5 -> 0.4.
    [Test]
    public void CalculateGcFraction_Ignore_DenominatorIsFullLength()
    {
        double result = "GCATN".CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Ignore);

        Assert.That(result, Is.EqualTo(0.4).Within(Tol),
            "Ignore: gc=count(C,G,S)=2 over full length 5 -> 0.4 (Biopython ambiguous='ignore').");
    }

    // Weighted mode: ambiguity codes add their mean GC value; denominator = full length.
    // "GCATN": gc = 2 (G,C) + N*0.5 = 2.5; length = 5 -> 0.5.
    [Test]
    public void CalculateGcFraction_Weighted_NContributesHalf()
    {
        double result = "GCATN".CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Weighted);

        Assert.That(result, Is.EqualTo(0.5).Within(Tol),
            "Weighted: N=0.5 -> gc=2.5 over length 5 -> 0.5 (Biopython _gc_values[N]=0.5).");
    }

    // Weighted: V = 2/3, H = 1/3 per Biopython _gc_values.
    // "VH": gc = 2/3 + 1/3 = 1.0; length = 2 -> 0.5.
    [Test]
    public void CalculateGcFraction_Weighted_VAndHUseExactGcValues()
    {
        double result = "VH".CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Weighted);

        Assert.That(result, Is.EqualTo(0.5).Within(Tol),
            "Weighted: V=2/3, H=1/3 -> gc=1.0 over length 2 -> 0.5 (Biopython _gc_values).");
    }

    // Empty input returns 0 in every mode ("Note that this will return zero for an empty sequence").
    [Test]
    public void CalculateGcFraction_EmptyOrNull_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That("".CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Remove),
                Is.EqualTo(0.0).Within(Tol), "Empty -> 0 (Remove).");
            Assert.That(((string?)null)!.CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Weighted),
                Is.EqualTo(0.0).Within(Tol), "Null -> 0 (Weighted).");
        });
    }

    // The opt-in overload must NOT change the existing parameterless default (A/T/G/C/U only).
    [Test]
    public void CalculateGcFraction_DefaultOverload_Unchanged()
    {
        Assert.That("GCAT".CalculateGcFraction(), Is.EqualTo(0.5).Within(Tol),
            "Default parameterless overload remains (G+C)/(A+T+G+C) fraction.");
    }

    #endregion

    #region AnalyzeGcContent(fraction:) — GcSkewCalculator

    // Default keeps percentage; opt-in fraction divides by 100. GCATGCAT -> 50% / 0.5.
    [Test]
    public void AnalyzeGcContent_FractionTrue_ReportsZeroToOne()
    {
        string seq = "GCATGCAT";
        var pct = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 4, stepSize: 4, fraction: false);
        var frac = GcSkewCalculator.AnalyzeGcContent(seq, windowSize: 4, stepSize: 4, fraction: true);

        Assert.Multiple(() =>
        {
            Assert.That(pct.OverallGcContent, Is.EqualTo(50.0).Within(Tol),
                "Default still emits percentage (50%).");
            Assert.That(frac.OverallGcContent, Is.EqualTo(0.5).Within(Tol),
                "fraction:true emits Biopython [0,1] (0.5).");
            Assert.That(frac.WindowedGcContent.All(p => p.GcContent <= 1.0), Is.True,
                "Windowed GC content is in [0,1] under fraction:true.");
            Assert.That(frac.WindowedGcContent[0].GcContent, Is.EqualTo(0.5).Within(Tol),
                "First window 'GCAT' -> 0.5 fraction.");
        });
    }

    #endregion

    #region CalculateGcContentProfile(fraction:) — SequenceStatistics

    [Test]
    public void CalculateGcContentProfile_FractionTrue_ReportsZeroToOne()
    {
        string seq = "GCATGC";
        List<double> pct = SequenceStatistics.CalculateGcContentProfile(seq, windowSize: 2, stepSize: 2, fraction: false).ToList();
        List<double> frac = SequenceStatistics.CalculateGcContentProfile(seq, windowSize: 2, stepSize: 2, fraction: true).ToList();

        Assert.Multiple(() =>
        {
            // Windows: "GC"=100%, "AT"=0%, "GC"=100%.
            Assert.That(pct, Is.EqualTo(new[] { 100.0, 0.0, 100.0 }).AsCollection,
                "Default percentage profile unchanged.");
            Assert.That(frac, Is.EqualTo(new[] { 1.0, 0.0, 1.0 }).AsCollection,
                "fraction:true emits the same windows in [0,1].");
        });
    }

    #endregion

    #region Variant.VcfPosition — VCF 1-based

    [Test]
    public void VcfPosition_IsOneBasedOffsetOfInternalPosition()
    {
        var v = new Variant(Position: 0, ReferenceAllele: "A", AlternateAllele: "G",
            Type: VariantType.SNP, QueryPosition: 0);

        Assert.Multiple(() =>
        {
            Assert.That(v.Position, Is.EqualTo(0), "Internal Position stays 0-based.");
            Assert.That(v.VcfPosition, Is.EqualTo(1),
                "VcfPosition = Position + 1 (VCF v4.3 §1.4.1: first base is POS 1).");
        });
    }

    [Test]
    public void VcfPosition_MatchesToVcfLinesPosColumn()
    {
        var variants = new[]
        {
            new Variant(Position: 41, ReferenceAllele: "C", AlternateAllele: "T",
                Type: VariantType.SNP, QueryPosition: 41)
        };

        // The POS column emitted by ToVcfLines must equal Variant.VcfPosition (both 1-based).
        string dataLine = VariantCaller.ToVcfLines(variants, chromosome: "chr1").Last();
        string posColumn = dataLine.Split('\t')[1];

        Assert.That(posColumn, Is.EqualTo(variants[0].VcfPosition.ToString()),
            "ToVcfLines POS column equals Variant.VcfPosition (42).");
    }

    #endregion
}
