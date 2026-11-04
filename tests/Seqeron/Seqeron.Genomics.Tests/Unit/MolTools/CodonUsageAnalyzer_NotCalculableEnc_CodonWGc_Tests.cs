// CODON-ENC-001 (A1-11) — explicit not-calculable ENC API (TryCalculateEnc, IsEncCalculable)
// CODON-STATS-001 (A1-14) — CodonW GC convention for GC1/GC2/GC3/OverallGc (stop codons excluded)
// Evidence: docs/Evidence/CODON-ENC-001-Evidence.md, docs/Evidence/CODON-STATS-001-Evidence.md
// TestSpec: tests/TestSpecs/CODON-ENC-001.md (M30–M34), tests/TestSpecs/CODON-STATS-001.md (R11–R16)
// Reference: CodonW 1.4.4 (Ubuntu codonw_1.4.4.orig.tar.gz, compiled locally) run as
//   `codonw genes.fa -base -enc -gc -gc3s -L_aa -code <c> -silent -nomenu -nowarn -machine`
//   (codon_us.c enc_out prints "*****" when Nc is not calculated; gc_out skips stop codons,
//   `if (pcu->ca[id] == 11) continue;`, denominator Len_aa); EMBOSS 6.6.0.0 `cusp` binary.
//   CodonW prints 3 dp (GC indices, fractions) / 2 dp (Nc); the exact rationals below are the
//   per-position counts of each gene and agree with CodonW's printed values.

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonUsageAnalyzer_NotCalculableEnc_CodonWGc_Tests
{
    private const double Tol = 1e-9;

    private static string Repeat(string s, int n) => string.Concat(Enumerable.Repeat(s, n));

    // Fully populated biased gene (CODON-ENC-001 M3, CodonW -enc 41.29): Nc = 41.288461538461526.
    private static readonly string BiasedGene =
        Repeat("TTT", 4) + "TTC" +
        Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
        Repeat("ATT", 3) + Repeat("ATC", 2) + "ATA" +
        Repeat("GTG", 4) + "GTC" +
        Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
        Repeat("CGC", 4) + Repeat("CGT", 2) +
        Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

    // Random gene (seed 20261014, "r53"): 52 triplets + 1 trailing base; one triplet with N
    // (skipped), 51 valid codons: 48 sense codons under table 1 (three TAG stops), 46 under
    // NCBI 2 (three TAG + two AGA stops).
    // CodonW 1.4.4: -code 0 Nc 42.85, Len_aa 48, GC 0.521, GC1 0.562, GC2 0.521, GC3 0.479, GC3s 0.479;
    //               -code 1 (NCBI 2) Nc 40.38, Len_aa 46, GC 0.529, GC1 0.587, GC2 0.500, GC3 0.500, GC3s 0.500.
    private const string R53 =
        "TTGGGTACACCCAACGGCGGTGGATAGCGCTTGCGTTATGAAAACCCCGTGCCCATCTTTGCTTAGGCAAAGTCCGCAAGACTAATACGCC" +
        "AAGTNGGATTTAGTGACCACCCGCTGTAGGAGGTGTGTACAAGTAAAACCGTTAGAGTACCGTACA";

    // Random gene (same seed, "r155"): CodonW -code 0 / -code 1 both print Nc "*****".
    private const string R155 =
        "CACTTCCCATTTGTCTAGCATGTAACRGGACATCTGCCGAGCAATNTAGTTAAAATGTATACCGTTCGAYAGGTAAAATCGGGAATCAA" +
        "GGTGGTCGGATTTCGCGCCAAGTGCATAGTRCTTCCAGAGGTCGAGATCCACGGGAAGGAATGAATGTAGGTC";

    #region TryCalculateEnc — explicit not-calculable result (A1-11)

    // CodonW prints "*****" (Nc not calculated) for each of these (CODON-ENC-001 F15 cases + r155).
    [TestCase("TTTTTTTTTTTC", 1, TestName = "TryCalculateEnc_PheOnly_CodonWStars_False")]
    [TestCase("ATGAAAGAGCTGTTCGCCAAA", 1, TestName = "TryCalculateEnc_ShortGene_CodonWStars_False")]
    [TestCase("ATGAAACTGTTCGCCAAATAA", 1, TestName = "TryCalculateEnc_ShortGeneWithStop_CodonWStars_False")]
    [TestCase("ATGTGGATGTGG", 1, TestName = "TryCalculateEnc_MetTrpOnly_CodonWStars_False")]
    [TestCase(R155, 1, TestName = "TryCalculateEnc_RandomGeneR155_Table1_CodonWStars_False")]
    [TestCase(R155, 2, TestName = "TryCalculateEnc_RandomGeneR155_Table2_CodonWStars_False")]
    public void TryCalculateEnc_CodonWNotCalculated_ReturnsFalseAndZero(string seq, int table)
    {
        var code = GeneticCode.GetByTableNumber(table);
        bool ok = CodonUsageAnalyzer.TryCalculateEnc(seq, code, out double enc);

        Assert.Multiple(() =>
        {
            Assert.That(ok, Is.False);
            Assert.That(enc, Is.EqualTo(0.0));
            Assert.That(CodonUsageAnalyzer.CalculateEnc(seq, code), Is.EqualTo(0.0), "old API keeps the 0 sentinel");
            Assert.That(CodonUsageAnalyzer.GetStatistics(seq, code).IsEncCalculable, Is.False);
        });
    }

    // Calculable genes: value equal to CodonW (locked values: M3 41.2885 / 45.00 for tables 1/2,
    // r53 42.85 / 40.38).
    [TestCase(1, 41.288461538461526)]
    [TestCase(2, 45.0)]
    public void TryCalculateEnc_BiasedGene_TrueAndCodonWValue(int table, double expected)
    {
        var code = GeneticCode.GetByTableNumber(table);
        bool ok = CodonUsageAnalyzer.TryCalculateEnc(BiasedGene, code, out double enc);

        Assert.Multiple(() =>
        {
            Assert.That(ok, Is.True);
            Assert.That(enc, Is.EqualTo(expected).Within(Tol));
            Assert.That(CodonUsageAnalyzer.CalculateEnc(BiasedGene, code), Is.EqualTo(enc));
            var stats = CodonUsageAnalyzer.GetStatistics(BiasedGene, code);
            Assert.That(stats.IsEncCalculable, Is.True);
            Assert.That(stats.Enc, Is.EqualTo(enc));
        });
    }

    [TestCase(1, 42.85340746733685, 42.85)]
    [TestCase(2, 40.37735849056604, 40.38)]
    public void TryCalculateEnc_RandomGeneR53_MatchesCodonW(int table, double exact, double codonW2dp)
    {
        bool ok = CodonUsageAnalyzer.TryCalculateEnc(R53, GeneticCode.GetByTableNumber(table), out double enc);

        Assert.Multiple(() =>
        {
            Assert.That(ok, Is.True);
            Assert.That(enc, Is.EqualTo(exact).Within(Tol));
            Assert.That(Math.Round(enc, 2), Is.EqualTo(codonW2dp));
        });
    }

    [Test]
    public void TryCalculateEnc_DefaultCodeOverloads_EqualTable1()
    {
        Assert.That(CodonUsageAnalyzer.TryCalculateEnc(BiasedGene, out double s), Is.True);
        Assert.That(CodonUsageAnalyzer.TryCalculateEnc(new DnaSequence(BiasedGene), out double d), Is.True);
        Assert.That(CodonUsageAnalyzer.TryCalculateEnc(new DnaSequence(BiasedGene), GeneticCode.VertebrateMitochondrial, out double d2), Is.True);
        Assert.Multiple(() =>
        {
            Assert.That(s, Is.EqualTo(41.288461538461526).Within(Tol));
            Assert.That(d, Is.EqualTo(s));
            Assert.That(d2, Is.EqualTo(45.0).Within(Tol));
            Assert.That(CodonUsageAnalyzer.TryCalculateEnc("TTTTTTTTTTTC", out double z), Is.False);
            Assert.That(z, Is.EqualTo(0.0));
        });
    }

    [TestCase(null)]
    [TestCase("")]
    [TestCase("A")]
    [TestCase("NNNRYK")]
    public void TryCalculateEnc_NoCodons_False(string? seq)
    {
        Assert.That(CodonUsageAnalyzer.TryCalculateEnc(seq!, GeneticCode.Standard, out double enc), Is.False);
        Assert.That(enc, Is.EqualTo(0.0));
    }

    [Test]
    public void TryCalculateEnc_NullArguments_Throw()
    {
        Assert.Throws<ArgumentNullException>(() => CodonUsageAnalyzer.TryCalculateEnc("ATG", null!, out _));
        Assert.Throws<ArgumentNullException>(() => CodonUsageAnalyzer.TryCalculateEnc((DnaSequence)null!, out _));
        Assert.Throws<ArgumentNullException>(() =>
            CodonUsageAnalyzer.TryCalculateEnc(new DnaSequence("ATG"), null!, out _));
    }

    [Test]
    public void GetStatistics_EmptyInput_EncNotCalculable()
    {
        Assert.That(CodonUsageAnalyzer.GetStatistics("").IsEncCalculable, Is.False);
        Assert.That(default(CodonUsageStatistics).IsEncCalculable, Is.False);
    }

    #endregion

    #region GetStatistics GC convention — CodonW gc_out (A1-14)

    // `ATGAAACTGTTCGCCAAATAA` (terminal TAA). CodonW -code 0: Len_aa 6, GC 0.389, GC1 0.333,
    // GC2 0.167, GC3 0.667, GC3s 0.600. EMBOSS cusp includes the stop: 7 codons, GC 7/21.
    [Test]
    public void GetStatistics_CodonWConvention_TerminalStop_MatchesCodonW()
    {
        const string gene = "ATGAAACTGTTCGCCAAATAA";
        var w = CodonUsageAnalyzer.GetStatistics(gene, GeneticCode.Standard, CodonGcConvention.CodonW);
        var cusp = CodonUsageAnalyzer.GetStatistics(gene, GeneticCode.Standard);

        Assert.Multiple(() =>
        {
            Assert.That(w.Gc1, Is.EqualTo(100.0 * 2 / 6).Within(Tol));
            Assert.That(w.Gc2, Is.EqualTo(100.0 * 1 / 6).Within(Tol));
            Assert.That(w.Gc3, Is.EqualTo(100.0 * 4 / 6).Within(Tol));
            Assert.That(w.OverallGc, Is.EqualTo(100.0 * 7 / 18).Within(Tol));
            Assert.That(Math.Round(w.OverallGc / 100, 3), Is.EqualTo(0.389));
            Assert.That(w.Gc3s, Is.EqualTo(60.0).Within(Tol));
            Assert.That(w.TotalCodons, Is.EqualTo(7), "TotalCodons does not depend on the convention");
            Assert.That(cusp.OverallGc, Is.EqualTo(100.0 * 7 / 21).Within(Tol));
            Assert.That(cusp.Gc3s, Is.EqualTo(w.Gc3s));
        });
    }

    // Internal stops + ambiguous triplet + trailing partial base, tables 1 and 2 (stops differ:
    // under NCBI 2 AGA/AGG are also stops). Exact counts over the sense codons of r53.
    [TestCase(1, 48, 27, 25, 23, 0.521, 0.562, 0.521, 0.479, 0.479)]
    [TestCase(2, 46, 27, 23, 23, 0.529, 0.587, 0.500, 0.500, 0.500)]
    public void GetStatistics_CodonWConvention_RandomGeneR53_MatchesCodonW(
        int table, int lenAa, int g1, int g2, int g3,
        double cwGc, double cwGc1, double cwGc2, double cwGc3, double cwGc3s)
    {
        var s = CodonUsageAnalyzer.GetStatistics(R53, GeneticCode.GetByTableNumber(table), CodonGcConvention.CodonW);

        Assert.Multiple(() =>
        {
            Assert.That(s.Gc1, Is.EqualTo(100.0 * g1 / lenAa).Within(Tol));
            Assert.That(s.Gc2, Is.EqualTo(100.0 * g2 / lenAa).Within(Tol));
            Assert.That(s.Gc3, Is.EqualTo(100.0 * g3 / lenAa).Within(Tol));
            Assert.That(s.OverallGc, Is.EqualTo(100.0 * (g1 + g2 + g3) / (3 * lenAa)).Within(Tol));
            Assert.That(Math.Round(s.OverallGc / 100, 3), Is.EqualTo(cwGc).Within(1e-12));
            Assert.That(Math.Round(s.Gc1 / 100, 3), Is.EqualTo(cwGc1).Within(1e-12));
            Assert.That(Math.Round(s.Gc2 / 100, 3), Is.EqualTo(cwGc2).Within(1e-12));
            Assert.That(Math.Round(s.Gc3 / 100, 3), Is.EqualTo(cwGc3).Within(1e-12));
            Assert.That(Math.Round(s.Gc3s / 100, 3), Is.EqualTo(cwGc3s).Within(1e-12));
            Assert.That(s.TotalCodons, Is.EqualTo(51));
        });
    }

    // `ATGGCTTAAGCCNNNGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA`: 15 triplets, NNN skipped, 4 stops in table 1.
    // CodonW -code 0: Len_aa 10, GC 0.600, GC1 0.600, GC2 0.500, GC3 0.700, GC3s 0.625;
    // -code 1 (NCBI 2, TGA = Trp): Len_aa 11, GC 0.576, GC1 0.545, GC2 0.545, GC3 0.636, GC3s 0.636.
    [TestCase(1, 10, 6, 5, 7, 62.5)]
    [TestCase(2, 11, 6, 6, 7, 100.0 * 7 / 11)]
    public void GetStatistics_CodonWConvention_InternalStopsAndAmbiguous_MatchesCodonW(
        int table, int lenAa, int g1, int g2, int g3, double gc3s)
    {
        const string gene = "ATGGCTTAAGCCNNNGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA";
        var s = CodonUsageAnalyzer.GetStatistics(gene, GeneticCode.GetByTableNumber(table), CodonGcConvention.CodonW);

        Assert.Multiple(() =>
        {
            Assert.That(s.Gc1, Is.EqualTo(100.0 * g1 / lenAa).Within(Tol));
            Assert.That(s.Gc2, Is.EqualTo(100.0 * g2 / lenAa).Within(Tol));
            Assert.That(s.Gc3, Is.EqualTo(100.0 * g3 / lenAa).Within(Tol));
            Assert.That(s.OverallGc, Is.EqualTo(100.0 * (g1 + g2 + g3) / (3 * lenAa)).Within(Tol));
            Assert.That(s.Gc3s, Is.EqualTo(gc3s).Within(Tol));
            Assert.That(s.TotalCodons, Is.EqualTo(14));
        });
    }

    // EMBOSS 6.6.0.0 `cusp` binary on `ATGGCTTAAGCCGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA` (14 codons,
    // 4 stops): Coding GC 47.62%, 1st 42.86%, 2nd 42.86%, 3rd 57.14% — the default convention.
    [Test]
    public void GetStatistics_DefaultConvention_MatchesEmbossCusp()
    {
        const string gene = "ATGGCTTAAGCCGCGTGAGGCCTGCTTTAGAAAAAGTGGTAA";
        var d = CodonUsageAnalyzer.GetStatistics(gene);
        var e = CodonUsageAnalyzer.GetStatistics(gene, GeneticCode.Standard, CodonGcConvention.EmbossCusp);

        Assert.Multiple(() =>
        {
            Assert.That(d.TotalCodons, Is.EqualTo(14));
            Assert.That(d.Gc1, Is.EqualTo(100.0 * 6 / 14).Within(Tol));
            Assert.That(d.Gc2, Is.EqualTo(100.0 * 6 / 14).Within(Tol));
            Assert.That(d.Gc3, Is.EqualTo(100.0 * 8 / 14).Within(Tol));
            Assert.That(Math.Round(d.OverallGc, 2), Is.EqualTo(47.62));
            Assert.That(Math.Round(d.Gc1, 2), Is.EqualTo(42.86));
            Assert.That(Math.Round(d.Gc3, 2), Is.EqualTo(57.14));
            Assert.That(e.Gc1, Is.EqualTo(d.Gc1));
            Assert.That(e.Gc2, Is.EqualTo(d.Gc2));
            Assert.That(e.Gc3, Is.EqualTo(d.Gc3));
            Assert.That(e.TotalCodons, Is.EqualTo(d.TotalCodons));
        });
    }

    // Stop-free gene: the two conventions coincide (CodonW == cusp, existing Evidence).
    [Test]
    public void GetStatistics_StopFreeGene_ConventionsAgree()
    {
        var a = CodonUsageAnalyzer.GetStatistics(BiasedGene, GeneticCode.Standard, CodonGcConvention.CodonW);
        var b = CodonUsageAnalyzer.GetStatistics(BiasedGene, GeneticCode.Standard, CodonGcConvention.EmbossCusp);
        Assert.Multiple(() =>
        {
            Assert.That(a.Gc1, Is.EqualTo(b.Gc1));
            Assert.That(a.Gc2, Is.EqualTo(b.Gc2));
            Assert.That(a.Gc3, Is.EqualTo(b.Gc3));
        });
    }

    // Only stop codons: CodonW writes no GC value ("appear to be too short"); 0 here.
    [Test]
    public void GetStatistics_CodonWConvention_StopsOnly_Zero()
    {
        var s = CodonUsageAnalyzer.GetStatistics("TAATAGTGA", GeneticCode.Standard, CodonGcConvention.CodonW);
        Assert.Multiple(() =>
        {
            Assert.That(s.Gc1, Is.EqualTo(0.0));
            Assert.That(s.Gc2, Is.EqualTo(0.0));
            Assert.That(s.Gc3, Is.EqualTo(0.0));
            Assert.That(s.TotalCodons, Is.EqualTo(3));
        });
    }

    [Test]
    public void GetStatistics_CodonWConvention_DnaSequenceOverload_EqualsString()
    {
        const string gene = "ATGAAACTGTTCGCCAAATAA";
        var a = CodonUsageAnalyzer.GetStatistics(new DnaSequence(gene), GeneticCode.Standard, CodonGcConvention.CodonW);
        var b = CodonUsageAnalyzer.GetStatistics(gene, GeneticCode.Standard, CodonGcConvention.CodonW);
        Assert.That(a.Gc1, Is.EqualTo(b.Gc1));
        Assert.That(a.OverallGc, Is.EqualTo(b.OverallGc));
    }

    [Test]
    public void GetStatistics_GcConvention_InvalidArguments_Throw()
    {
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            CodonUsageAnalyzer.GetStatistics("ATG", GeneticCode.Standard, (CodonGcConvention)7));
        Assert.Throws<ArgumentNullException>(() =>
            CodonUsageAnalyzer.GetStatistics("ATG", null!, CodonGcConvention.CodonW));
        Assert.Throws<ArgumentNullException>(() =>
            CodonUsageAnalyzer.GetStatistics((DnaSequence)null!, GeneticCode.Standard, CodonGcConvention.CodonW));
    }

    #endregion
}
