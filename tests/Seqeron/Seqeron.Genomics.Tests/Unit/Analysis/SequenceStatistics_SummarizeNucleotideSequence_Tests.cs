// SEQ-SUMMARY-001 — Sequence Summary (aggregation of length, GC, entropy, complexity, Tm, composition)
// Evidence: docs/Evidence/SEQ-SUMMARY-001-Evidence.md
// TestSpec: tests/TestSpecs/SEQ-SUMMARY-001.md
// Source: Biopython Bio.SeqUtils gc_fraction / MeltingTemp (Cock et al. 2009);
//         Shannon C.E. (1948) A Mathematical Theory of Communication;
//         Orlov Y.L. & Potapov V.N. (2004) NAR 32:W628; Troyanskaya O.G. et al. (2002) linguistic complexity.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

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
            // Complexity = Σ V_k / Σ min(4^k, N−k+1), k=1..6 (Orlov & Potapov 2004 summation form,
            // canonical SequenceComplexity; B03 F21 replaced the former unsourced mean of U_k = 529/630).
            // Hand-computed: V = 4,4,4,4,4,3 (23); V_max = 4,7,6,5,4,3 (29) -> 23/29.
            Assert.That(summary.Complexity, Is.EqualTo(23.0 / 29.0).Within(Tolerance),
                "Complexity = Σ V_k / Σ V_max,k over k=1..6 = 23/29");
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
                Is.EqualTo(SequenceStatistics.CalculateMeltingTemperature(seq, useWallaceRule: true)).Within(Tolerance),
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

    // 2026-09 B03 review (SEQ-SUMMARY-001): every field locked to an executed Python reference —
    // GcContent = Biopython 1.88 gc_fraction(seq, "remove"); Entropy = scipy.stats.entropy(counts, base=2);
    // Tm = Biopython Tm_Wallace (< 14 A/C/G/T/U) or Tm_GC(userset=(64.9, 0.41, 672.4, 0), saltcorr=0)
    // (= 64.9 + 41·(G+C − 16.4)/N; both back-transcribe U, F20); Complexity = exact-fraction
    // Σ V_k / Σ min(4^k, N − k + 1), k = 1..6 (Orlov & Potapov 2004; Python reference, B03 F21).
    [TestCase("ATGCATGC", 0.5, 2.0, 24.0, 23.0 / 29.0)]
    [TestCase("ACGTACGGTACCAGTTAGCA", 0.5, 1.9854752972273346, 51.78000000000001, 38.0 / 43.0)]
    [TestCase("AUGCAUGC", 0.5, 2.0, 24.0, 23.0 / 29.0)]
    [TestCase("GGGAAAUUUCCCAAAUGC", 0.4444444444444444, 1.974937501201927, 45.76666666666668, 23.0 / 26.0)]
    [TestCase("ATTTGGATT", 0.2222222222222222, 1.4355205042826666, 22.0, 29.0 / 34.0)]
    public void SummarizeNucleotideSequence_MatchesPythonReferences(
        string seq, double gc, double entropy, double tm, double complexity)
    {
        var s = SequenceStatistics.SummarizeNucleotideSequence(seq);
        Assert.Multiple(() =>
        {
            Assert.That(s.Length, Is.EqualTo(seq.Length));
            Assert.That(s.GcContent, Is.EqualTo(gc).Within(1e-12), "Biopython gc_fraction");
            Assert.That(s.Entropy, Is.EqualTo(entropy).Within(1e-12), "scipy entropy base 2");
            Assert.That(s.MeltingTemperature, Is.EqualTo(tm).Within(1e-9), "Biopython Tm_Wallace / Tm_GC");
            Assert.That(s.Complexity, Is.EqualTo(complexity).Within(1e-12), "Σ V_k / Σ V_max,k, k=1..6");
        });
    }

    // RNA spelling: U is read as T by every Tm/GC component (Biopython back-transcription, F20), so an
    // RNA and its DNA spelling have identical GC, entropy, complexity and Tm; only the T/U counts move.
    [TestCase("AUGCAUGC", "ATGCATGC")]
    [TestCase("GGGAAAUUUCCCAAAUGC", "GGGAAATTTCCCAAATGC")]
    [TestCase("acguacgu", "ACGTACGT")]
    public void SummarizeNucleotideSequence_RnaSpelling_EqualsDnaSpelling(string rna, string dna)
    {
        var r = SequenceStatistics.SummarizeNucleotideSequence(rna);
        var d = SequenceStatistics.SummarizeNucleotideSequence(dna);
        Assert.Multiple(() =>
        {
            Assert.That(r.GcContent, Is.EqualTo(d.GcContent));
            Assert.That(r.Entropy, Is.EqualTo(d.Entropy).Within(1e-12));
            Assert.That(r.Complexity, Is.EqualTo(d.Complexity).Within(1e-12));
            Assert.That(r.MeltingTemperature, Is.EqualTo(d.MeltingTemperature), "Tm reads U as T");
            Assert.That(r.Composition['U'], Is.EqualTo(d.Composition['T']));
            Assert.That(r.Composition['T'], Is.EqualTo(0));
        });
    }

    // LINGUISTIC (B03 F21): SequenceStatistics.CalculateLinguisticComplexity delegates to the canonical
    // SequenceComplexity sum form (Orlov & Potapov 2004; Troyanskaya 2002 / Rosalind LING at m ≥ N).
    // Values are exact fractions from an executed Python reference of Σ V_k / Σ min(4^k, N−k+1).
    [TestCase("ATTTGGATT", 9, 7.0 / 8.0)]          // Rosalind LING sample, full m: 0.875
    [TestCase("ATTTGGATT", int.MaxValue, 7.0 / 8.0)]
    [TestCase("ATTTGGATT", 6, 29.0 / 34.0)]         // default m = 6 (was 293/336 = mean of U_k)
    [TestCase("AAAAAAAAAA", 6, 2.0 / 13.0)]
    [TestCase("ATGCATGCATGC", 6, 24.0 / 49.0)]
    [TestCase("ATATATAT", 6, 12.0 / 29.0)]
    [TestCase("GGGAAAUUUCCC", 6, 45.0 / 49.0)]
    [TestCase("ACGTN", 6, 15.0 / 14.0)]             // N is a fifth symbol: > 1 (canonical, R23)
    [TestCase("ATGC", 0, 0.0)]
    [TestCase("", 6, 0.0)]
    public void CalculateLinguisticComplexity_MatchesSumFormReference(string seq, int maxK, double expected)
    {
        Assert.That(SequenceStatistics.CalculateLinguisticComplexity(seq, maxK),
            Is.EqualTo(expected).Within(1e-12));
    }

    // Differential lock: identical to the canonical implementation for the same m (maxK ≡ maxWordLength).
    [Test]
    public void CalculateLinguisticComplexity_EqualsCanonicalSequenceComplexity()
    {
        var rng = new Random(20260928);
        const string alphabet = "ACGTNacgtu-";
        for (int t = 0; t < 1000; t++)
        {
            int n = rng.Next(0, 40);
            var chars = new char[n];
            for (int i = 0; i < n; i++) chars[i] = alphabet[rng.Next(t % 2 == 0 ? 4 : alphabet.Length)];
            string seq = new(chars);
            int m = rng.Next(-1, 20);
            Assert.That(SequenceStatistics.CalculateLinguisticComplexity(seq, m),
                Is.EqualTo(SequenceComplexity.CalculateLinguisticComplexity(seq, m)), $"seq={seq} m={m}");
        }
        Assert.That(SequenceStatistics.CalculateLinguisticComplexity(null!), Is.EqualTo(0.0));
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
