// SEQ-SUMMARY-001 — Sequence Summary (aggregation of length, GC, entropy, complexity, Tm, composition)
// Evidence: docs/Evidence/SEQ-SUMMARY-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-SUMMARY-001.md
// Source: Biopython Bio.SeqUtils gc_fraction / MeltingTemp (Cock et al. 2009);
//         Shannon C.E. (1948) A Mathematical Theory of Communication;
//         Trifonov E.N. (1990) linguistic sequence complexity.

using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

[TestFixture]
public class SequenceStatistics_SummarizeNucleotideSequence_Tests
{
    private const double Tolerance = 1e-10;

    #region SummarizeNucleotideSequence

    // M1 — Worked example "ATGCATGC": composition A=T=G=C=2 -> GcContent = 4/8 = 0.5;
    // four equally frequent symbols -> Shannon H = log2(4) = 2.0 bits; length 8 < 14 ->
    // Wallace Tm = 2*(A+T) + 4*(G+C) = 2*4 + 4*4 = 24.0. Values derived from the cited formulas.
    [Test]
    public void SummarizeNucleotideSequence_BalancedTetramer_ReturnsExactValues()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence("ATGCATGC");

        Assert.Multiple(() =>
        {
            Assert.That(summary.Length, Is.EqualTo(8),
                "Length is the raw character count (INV-01)");
            Assert.That(summary.GcContent, Is.EqualTo(0.5).Within(Tolerance),
                "GC fraction = (G+C)/total = 4/8 = 0.5 per Biopython gc_fraction");
            Assert.That(summary.Entropy, Is.EqualTo(2.0).Within(Tolerance),
                "uniform over 4 symbols -> H = log2(4) = 2.0 bits per Shannon 1948");
            Assert.That(summary.MeltingTemperature, Is.EqualTo(24.0).Within(Tolerance),
                "len<14 -> Wallace Tm = 2*(A+T)+4*(G+C) = 8+16 = 24.0");
            // Complexity is the mean of per-word-size vocabulary-usage ratios U_k = observed/possible,
            // k=1..6 (Trifonov 1990 vocabulary usage; this method uses the mean, not the product).
            // Hand-computed independently: k=1: 4/4=1; k=2: 4/7; k=3: 4/6; k=4: 4/5; k=5: 4/4=1; k=6: 3/3=1
            // mean = (1 + 4/7 + 2/3 + 0.8 + 1 + 1)/6 = 0.83968253968253968...
            Assert.That(summary.Complexity, Is.EqualTo(0.8396825396825397).Within(Tolerance),
                "Complexity = mean of vocabulary-usage ratios over k=1..6 = 0.8396825396825397");
        });
    }

    // M2 — Every summary field must equal the canonical per-metric method's value on the
    // same input (aggregation-consistency, INV-02..INV-06). Uses a 16-mer (GC branch for Tm).
    [Test]
    public void SummarizeNucleotideSequence_AllFields_EqualCanonicalMethods()
    {
        const string seq = "ATGCATGCATGCATGC";
        var summary = SequenceStatistics.SummarizeNucleotideSequence(seq);
        var comp = SequenceStatistics.CalculateNucleotideComposition(seq);

        Assert.Multiple(() =>
        {
            Assert.That(summary.Length, Is.EqualTo(comp.Length),
                "INV-01: Length equals composition length");
            Assert.That(summary.GcContent, Is.EqualTo(comp.GcContent).Within(Tolerance),
                "INV-02: GcContent equals CalculateNucleotideComposition.GcContent");
            Assert.That(summary.Entropy,
                Is.EqualTo(SequenceStatistics.CalculateShannonEntropy(seq)).Within(Tolerance),
                "INV-03: Entropy equals CalculateShannonEntropy");
            Assert.That(summary.Complexity,
                Is.EqualTo(SequenceStatistics.CalculateLinguisticComplexity(seq)).Within(Tolerance),
                "INV-04: Complexity equals CalculateLinguisticComplexity");
            Assert.That(summary.MeltingTemperature,
                Is.EqualTo(SequenceStatistics.CalculateMeltingTemperature(seq, useWallaceRule: seq.Length < 14)).Within(Tolerance),
                "INV-05: MeltingTemperature equals CalculateMeltingTemperature with the len<14 flag");
        });
    }

    // M3 — Tm uses the GC/Marmur-Doty branch for length >= 14. For "ATGCATGCATGCATGC"
    // (len 16, GC 8): Tm = 64.9 + 41*(8-16.4)/16 = 43.375 (repo Marmur-Doty variant, SEQ-TM-001).
    [Test]
    public void SummarizeNucleotideSequence_LengthAtLeast14_UsesGcFormula()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence("ATGCATGCATGCATGC");

        Assert.That(summary.MeltingTemperature, Is.EqualTo(43.375).Within(Tolerance),
            "len>=14 -> GC formula Tm = 64.9 + 41*(GC-16.4)/N = 64.9 + 41*(8-16.4)/16 = 43.375");
    }

    // M4 — Tm uses the Wallace branch for length < 14. "ATGC" (A+T=2, G+C=2):
    // Tm = 2*2 + 4*2 = 12.0; must equal the canonical method with useWallaceRule:true.
    [Test]
    public void SummarizeNucleotideSequence_ShortSequence_UsesWallaceRule()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence("ATGC");

        Assert.Multiple(() =>
        {
            Assert.That(summary.MeltingTemperature, Is.EqualTo(12.0).Within(Tolerance),
                "len<14 -> Wallace Tm = 2*(A+T)+4*(G+C) = 4+8 = 12.0");
            Assert.That(summary.MeltingTemperature,
                Is.EqualTo(SequenceStatistics.CalculateMeltingTemperature("ATGC", useWallaceRule: true)).Within(Tolerance),
                "Wallace branch matches the canonical method with useWallaceRule:true");
        });
    }

    // M5 — Composition dictionary counts (A,T,G,C,U,N) must equal CalculateNucleotideComposition.
    // RNA + ambiguous input "AUGCNNA": A=2, U=1, G=1, C=1, N=2, T=0.
    [Test]
    public void SummarizeNucleotideSequence_RnaWithN_CompositionMatchesCounts()
    {
        const string seq = "AUGCNNA";
        var summary = SequenceStatistics.SummarizeNucleotideSequence(seq);
        var comp = SequenceStatistics.CalculateNucleotideComposition(seq);

        Assert.Multiple(() =>
        {
            Assert.That(summary.Composition['A'], Is.EqualTo(2), "two A's; equals comp.CountA");
            Assert.That(summary.Composition['A'], Is.EqualTo(comp.CountA), "INV-06: A count matches composition");
            Assert.That(summary.Composition['U'], Is.EqualTo(1), "one U; RNA base counted");
            Assert.That(summary.Composition['U'], Is.EqualTo(comp.CountU), "INV-06: U count matches composition");
            Assert.That(summary.Composition['G'], Is.EqualTo(1), "one G");
            Assert.That(summary.Composition['C'], Is.EqualTo(1), "one C");
            Assert.That(summary.Composition['N'], Is.EqualTo(2), "two N (ambiguous) counted");
            Assert.That(summary.Composition['N'], Is.EqualTo(comp.CountN), "INV-06: N count matches composition");
            Assert.That(summary.Composition['T'], Is.EqualTo(0), "no T in an RNA sequence");
        });
    }

    // M6 — Empty input returns the degenerate summary (all zero), per the empty-sequence
    // handling of every per-metric method (gc_fraction returns 0 on empty).
    [Test]
    public void SummarizeNucleotideSequence_EmptyString_ReturnsZeroSummary()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence("");

        Assert.Multiple(() =>
        {
            Assert.That(summary.Length, Is.EqualTo(0), "empty length is 0");
            Assert.That(summary.GcContent, Is.EqualTo(0.0).Within(Tolerance), "empty GcContent is 0");
            Assert.That(summary.Entropy, Is.EqualTo(0.0).Within(Tolerance), "empty entropy is 0");
            Assert.That(summary.Complexity, Is.EqualTo(0.0).Within(Tolerance), "empty complexity is 0");
            Assert.That(summary.MeltingTemperature, Is.EqualTo(0.0).Within(Tolerance), "empty Tm is 0");
            Assert.That(summary.Composition['A'], Is.EqualTo(0), "empty composition counts are 0");
        });
    }

    // S1 — Null input is guarded identically to empty input (no throw).
    [Test]
    public void SummarizeNucleotideSequence_Null_ReturnsZeroSummary()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence(null!);

        Assert.Multiple(() =>
        {
            Assert.That(summary.Length, Is.EqualTo(0), "null length is 0");
            Assert.That(summary.GcContent, Is.EqualTo(0.0).Within(Tolerance), "null GcContent is 0");
            Assert.That(summary.Entropy, Is.EqualTo(0.0).Within(Tolerance), "null entropy is 0");
            Assert.That(summary.MeltingTemperature, Is.EqualTo(0.0).Within(Tolerance), "null Tm is 0");
        });
    }

    // S2 — Case-insensitivity: lowercase input yields the identical summary as uppercase,
    // because each per-metric method uppercases internally.
    [Test]
    public void SummarizeNucleotideSequence_LowercaseInput_MatchesUppercase()
    {
        var lower = SequenceStatistics.SummarizeNucleotideSequence("atgcatgc");
        var upper = SequenceStatistics.SummarizeNucleotideSequence("ATGCATGC");

        Assert.Multiple(() =>
        {
            Assert.That(lower.GcContent, Is.EqualTo(upper.GcContent).Within(Tolerance), "GcContent case-insensitive");
            Assert.That(lower.Entropy, Is.EqualTo(upper.Entropy).Within(Tolerance), "Entropy case-insensitive");
            Assert.That(lower.Complexity, Is.EqualTo(upper.Complexity).Within(Tolerance), "Complexity case-insensitive");
            Assert.That(lower.MeltingTemperature, Is.EqualTo(upper.MeltingTemperature).Within(Tolerance), "Tm case-insensitive");
            Assert.That(lower.Composition['G'], Is.EqualTo(upper.Composition['G']), "composition case-insensitive");
        });
    }

    // C1 — Bounds invariant (INV-07): 0 <= GcContent <= 1 and 0 <= Complexity < 1 for a DNA fragment.
    [Test]
    public void SummarizeNucleotideSequence_DnaFragment_RespectsBounds()
    {
        var summary = SequenceStatistics.SummarizeNucleotideSequence("ATGGCCATTGCATAGCTAGCT");

        Assert.Multiple(() =>
        {
            Assert.That(summary.GcContent, Is.InRange(0.0, 1.0), "GC fraction is bounded in [0,1]");
            Assert.That(summary.Complexity, Is.GreaterThan(0.0), "linguistic complexity is positive for a real fragment");
            Assert.That(summary.Complexity, Is.LessThan(1.0), "linguistic complexity stays below 1 for a non-maximal fragment");
        });
    }

    #endregion
}
