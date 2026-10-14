// CODON-ENC-001 (audit A2-7, A2-8) — opt-in Nc estimators and the CodonW cap option:
//   EncMethod.Novembre2002 (background-corrected Nc′), EncMethod.SunYangXia2013, EncOptions.CodonWCap.
// Evidence: docs/Evidence/CODON-ENC-001-Evidence.md; TestSpec: tests/TestSpecs/CODON-ENC-001.md (M40–M47)
// References:
//   codonbias 0.5.0 (PyPI codon-bias, scores.py EffectiveNumberOfCodons, k_mer=1), run as
//     Novembre 2002:   EffectiveNumberOfCodons(bg_correction=True, robust=False, pseudocount=0, mean='unweighted', genetic_code=t)
//     Sun–Yang–Xia:    EffectiveNumberOfCodons(genetic_code=t)  (robust=True, pseudocount=1, mean='weighted')
//     + variants (bg_correction=True; pseudocount=0.5 weighted; pseudocount=0 unweighted; get_score(background=…)).
//   CodonW 1.4.4 (Ubuntu codonw_1.4.4.orig.tar.gz, compiled locally): `codonw genes.fa -enc -code c -silent -nomenu
//     -nowarn -machine`, c = 0..7 = NCBI 1, 2, 3, 4, 5, 6, 10, 9; enc_out prints any Nc > 61 as "61.00".
//   All codonbias values below agree with the C# result to ≤ 2.3e-15 relative (57351 gene × table × option cases).

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonUsageAnalyzer_EncVariants_Tests
{
    private const double Tol = 1e-9;

    private static string Repeat(string s, int n) => string.Concat(Enumerable.Repeat(s, n));

    private static readonly string[] AllCodons =
        (from a in "ACGT" from b in "ACGT" from c in "ACGT" select $"{a}{b}{c}").ToArray();

    // CODON-ENC-001 M3 (CodonW 41.29).
    private static readonly string M3 =
        Repeat("TTT", 4) + "TTC" +
        Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
        Repeat("ATT", 3) + Repeat("ATC", 2) + "ATA" +
        Repeat("GTG", 4) + "GTC" +
        Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
        Repeat("CGC", 4) + Repeat("CGT", 2) +
        Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

    // CODON-ENC-001 M10: codon_i × ((3i mod 7) + 1), codons in ACGT lexicographic order.
    private static readonly string Det =
        string.Concat(AllCodons.Select((c, i) => Repeat(c, (3 * i) % 7 + 1)));

    // CODON-ENC-001 M12: codon_i × ((7i mod 5) + 1) — Wright Nc overshoots every sense-codon cap.
    private static readonly string Over =
        string.Concat(AllCodons.Select((c, i) => Repeat(c, (7 * i) % 5 + 1)));

    // Random gene r53 (seed 20261014; one N triplet, in-frame TAG stops, trailing base).
    private const string R53 =
        "TTGGGTACACCCAACGGCGGTGGATAGCGCTTGCGTTATGAAAACCCCGTGCCCATCTTTGCTTAGGCAAAGTCCGCAAGACTAATACGCC" +
        "AAGTNGGATTTAGTGACCACCCGCTGTAGGAGGTGTGTACAAGTAAAACCGTTAGAGTACCGTACA";

    // Random genes r203 / r15 (seed 20261014, WP16 generator, 120 codons each).
    private const string R203 =
        "GACGACTTCTAAAAGTAGCTGTATCGTTAAGAAAAACCTCGTCCGGCGACTGGGAGTTCGGCACTTGAAGTGTATCAGGTTATCGAATTGATGCACCCGCCACACGAG" +
        "TTTATTGGGCTTTGATTCCACGAACTACACTCGAGCAGTATCCGCAACACTTAATGTGCCGCCGCGTCTCTAAACACTCACGCAATCCCACCGTCAGTCAATTCACTG" +
        "GGCACCCTGCATATAAAAGGCCGCGAGAAGTCGATGGAGCAGCGCGAGGCTATGCTGGTGCAATTGTACAGGTCGTTAACAGACCTCTCATCCATAGGCCATGCACGCC" +
        "ACCACAACCAGTCCATGAAAACCAATGAGGGAAGT";

    private const string R15 =
        "TAGCTTGCGGGGATTGCCGTGATGTATCGCATCGTGGTGGAACTCACTCGGGTCGGTTAGCCAAGAAGTCGACAGGATTGAAGTTTTCGGGGACCCTCGATACTCCTG" +
        "AATAACTATTGTGGGTGCTACTTTCTCGCATATAGCCTCAGGTTGCTTGTGACCCCGGGGAAGCCCTTAGCGGTGGCCAGACTGTGTGTAGCGAACTTCCCGCACACCT" +
        "TCCTCTTTGCGAGGATTGCTAGGATAGCCGCACTGGCTTGGCTTTAGGTGAAAGCCCTTCATCATACCCGAGCGACTGCTTGCCTGTCTGAAATCCGGCATTTATTTGG" +
        "ATTGGCCCATTACGCGTAAGCTCGTCAGTACTTG";

    private const string PheOnly = "TTTTTTTTTTTC";

    private const string Background = "GGGCCCGCGATATGCGCGGCCGCAAT";

    private static string Gene(string name) => name switch
    {
        "M3" => M3,
        "Det" => Det,
        "Over" => Over,
        "R53" => R53,
        "R203" => R203,
        "R15" => R15,
        "Phe" => PheOnly,
        _ => throw new ArgumentOutOfRangeException(nameof(name)),
    };

    private static double Enc(string gene, int table, EncOptions options)
    {
        bool ok = CodonUsageAnalyzer.TryCalculateEnc(gene, GeneticCode.GetByTableNumber(table), options, out double enc);
        Assert.That(ok, Is.True, "a gene with sense codons is always calculable under Novembre / Sun–Yang–Xia");
        return enc;
    }

    #region Novembre 2002 (A2-7)

    // codonbias EffectiveNumberOfCodons(bg_correction=True, robust=False, pseudocount=0, mean='unweighted').
    [TestCase("M3", 1, 45.49951134064658)]
    [TestCase("R53", 1, 59.23386199093016)]
    [TestCase("R203", 1, 60.502183497673535)]
    [TestCase("Det", 1, 61.0)] // re-adjusted to the 61 sense codons
    [TestCase("Phe", 1, 61.0)] // empty classes imputed with F̄ = 1/z
    [TestCase("M3", 2, 46.946215313358586)]
    [TestCase("R53", 2, 59.32323189503357)]
    [TestCase("Det", 2, 59.95683413399217)]
    [TestCase("R15", 2, 60.0)] // table 2: 60 sense codons
    public void Novembre2002_MatchesCodonbias(string gene, int table, double expected)
    {
        Assert.That(Enc(Gene(gene), table, new EncOptions { Method = EncMethod.Novembre2002 }),
            Is.EqualTo(expected).Within(Tol));
    }

    // codonbias get_score(seq, background=…): the background sequence's base composition replaces the gene's.
    [TestCase("M3", 1, 38.480446976398305)]
    [TestCase("Det", 1, 49.56413413981868)]
    [TestCase("Over", 1, 53.63482734039479)]
    [TestCase("R53", 1, 34.22859075786503)]
    [TestCase("R203", 1, 60.69278318528453)]
    [TestCase("R15", 1, 51.54164369125411)]
    [TestCase("M3", 2, 37.925611282065496)]
    [TestCase("R15", 2, 50.44213241219982)]
    public void Novembre2002_ExplicitBackground_MatchesCodonbias(string gene, int table, double expected)
    {
        var options = new EncOptions { Method = EncMethod.Novembre2002, Background = Background };
        Assert.That(Enc(Gene(gene), table, options), Is.EqualTo(expected).Within(Tol));
    }

    // codonbias (bg_correction=True, robust=False, pseudocount=0.5, mean='weighted').
    [TestCase("M3", 1, 54.84871855749488)]
    [TestCase("Phe", 1, 43.80149130574615)]
    [TestCase("M3", 2, 50.34127053389516)]
    [TestCase("Det", 2, 59.95795016993473)]
    public void Novembre2002_PseudocountWeighted_MatchesCodonbias(string gene, int table, double expected)
    {
        var options = new EncOptions { Method = EncMethod.Novembre2002, Pseudocount = 0.5, WeightedAverage = true };
        Assert.That(Enc(Gene(gene), table, options), Is.EqualTo(expected).Within(Tol));
    }

    [Test]
    public void Novembre2002_WithoutBackgroundCorrection_IsWrightHomozygosity()
    {
        // With uniform expected frequencies χ² = kΣp² − 1, so F = (N·χ² + N − k)/(k(N − 1)) = (NΣp² − 1)/(N − 1),
        // Wright's Eq. 1: on a gene where every class is estimable both estimators agree (M3: CodonW 41.29).
        var options = new EncOptions { Method = EncMethod.Novembre2002, BackgroundCorrection = false };
        Assert.That(Enc(M3, 1, options), Is.EqualTo(41.288461538461526).Within(Tol));
        Assert.That(Enc(M3, 1, options), Is.EqualTo(CodonUsageAnalyzer.CalculateEnc(M3)).Within(1e-12));
    }

    #endregion

    #region Sun, Yang & Xia 2013 (A2-7)

    // codonbias EffectiveNumberOfCodons() defaults (robust=True, pseudocount=1, mean='weighted', bg_correction=False).
    [TestCase("M3", 1, 50.957809330628805)]
    [TestCase("Det", 1, 54.99679076828325)]
    [TestCase("Over", 1, 55.376384357027874)]
    [TestCase("R53", 1, 52.22708180943929)]
    [TestCase("R203", 1, 53.55409422696348)]
    [TestCase("R15", 1, 53.479948012949244)]
    [TestCase("M3", 2, 50.70834399908201)]
    [TestCase("Det", 2, 53.5004850251249)]
    [TestCase("R53", 2, 51.40115206789308)]
    [TestCase("R203", 2, 52.696068170034266)]
    public void SunYangXia2013_MatchesCodonbias(string gene, int table, double expected)
    {
        Assert.That(Enc(Gene(gene), table, new EncOptions { Method = EncMethod.SunYangXia2013 }),
            Is.EqualTo(expected).Within(Tol));
    }

    [Test]
    public void SunYangXia2013_PheOnly_HandDerivation()
    {
        // Pseudocount 1: Phe (TTT 3+1, TTC 1+1) N = 6, F = (16 + 4)/36 = 5/9; every other 2-fold amino acid
        // N = 2, F = 1/2. Weighted F̄₂ = (6·5/9 + 8·2·1/2)/(6 + 16) = 17/33; F̄₃ = 1/3, F̄₄ = 1/4, F̄₆ = 1/6.
        // Nc = 2 + 9·33/17 + 3 + 20 + 18 = 1028/17 (codonbias 60.470588235294116).
        Assert.That(Enc(PheOnly, 1, new EncOptions { Method = EncMethod.SunYangXia2013 }),
            Is.EqualTo(1028.0 / 17.0).Within(Tol));
    }

    // codonbias (bg_correction=True) on top of the Sun–Yang–Xia defaults.
    [TestCase("M3", 1, 47.56338039024118)]
    [TestCase("Det", 1, 54.96527947242783)]
    [TestCase("Phe", 1, 26.608219785444312)]
    [TestCase("R15", 2, 52.798967973279275)]
    public void SunYangXia2013_BackgroundCorrected_MatchesCodonbias(string gene, int table, double expected)
    {
        var options = new EncOptions { Method = EncMethod.SunYangXia2013, BackgroundCorrection = true };
        Assert.That(Enc(Gene(gene), table, options), Is.EqualTo(expected).Within(Tol));
    }

    // codonbias (pseudocount=0, mean='unweighted').
    [TestCase("M3", 1, 33.91223204458498)]
    [TestCase("R53", 1, 34.396744240129955)]
    [TestCase("R203", 2, 44.45168033800177)]
    public void SunYangXia2013_NoPseudocountUnweighted_MatchesCodonbias(string gene, int table, double expected)
    {
        var options = new EncOptions { Method = EncMethod.SunYangXia2013, Pseudocount = 0, WeightedAverage = false };
        Assert.That(Enc(Gene(gene), table, options), Is.EqualTo(expected).Within(Tol));
    }

    #endregion

    #region Contracts (A2-7)

    [Test]
    public void DefaultOptions_BitIdenticalToWright()
    {
        foreach (string gene in new[] { M3, Det, Over, R53, R203, R15, PheOnly, "", "A" })
            foreach (int table in new[] { 1, 2, 3, 9 })
            {
                var code = GeneticCode.GetByTableNumber(table);
                bool expectedOk = CodonUsageAnalyzer.TryCalculateEnc(gene, code, out double expected);
                bool ok = CodonUsageAnalyzer.TryCalculateEnc(gene, code, EncOptions.Default, out double enc);
                Assert.Multiple(() =>
                {
                    Assert.That(ok, Is.EqualTo(expectedOk));
                    Assert.That(BitConverter.DoubleToInt64Bits(enc), Is.EqualTo(BitConverter.DoubleToInt64Bits(expected)));
                    Assert.That(CodonUsageAnalyzer.CalculateEnc(gene, code, new EncOptions()), Is.EqualTo(expected));
                });
            }
    }

    [TestCase(EncMethod.Novembre2002)]
    [TestCase(EncMethod.SunYangXia2013)]
    public void NoSenseCodon_NotCalculable(EncMethod method)
    {
        var options = new EncOptions { Method = method };
        Assert.Multiple(() =>
        {
            foreach (string gene in new[] { "", "A", "NNNRYK", "TAATAGTGA" })
            {
                Assert.That(CodonUsageAnalyzer.TryCalculateEnc(gene, GeneticCode.Standard, options, out double enc), Is.False, gene);
                Assert.That(enc, Is.EqualTo(0.0));
                Assert.That(CodonUsageAnalyzer.CalculateEnc(gene, GeneticCode.Standard, options), Is.EqualTo(0.0));
            }
            // AGG is a stop codon in table 2 only.
            Assert.That(CodonUsageAnalyzer.TryCalculateEnc("AGG", GeneticCode.GetByTableNumber(2), options, out _), Is.False);
            Assert.That(CodonUsageAnalyzer.TryCalculateEnc("AGG", GeneticCode.Standard, options, out _), Is.True);
        });
    }

    [TestCase(EncMethod.Novembre2002)]
    [TestCase(EncMethod.SunYangXia2013)]
    public void RnaLowerCaseAndDnaSequenceOverloads_Agree(EncMethod method)
    {
        var options = new EncOptions { Method = method };
        double expected = CodonUsageAnalyzer.CalculateEnc(M3, GeneticCode.Standard, options);
        string rnaLower = M3.Replace('T', 'U').ToLowerInvariant();
        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.CalculateEnc(rnaLower, GeneticCode.Standard, options), Is.EqualTo(expected).Within(1e-12));
            Assert.That(CodonUsageAnalyzer.CalculateEnc(new DnaSequence(M3), GeneticCode.Standard, options), Is.EqualTo(expected));
            Assert.That(CodonUsageAnalyzer.TryCalculateEnc(new DnaSequence(M3), GeneticCode.Standard, options, out double enc), Is.True);
            Assert.That(enc, Is.EqualTo(expected));
        });
    }

    [Test]
    public void InvalidOptions_Throw()
    {
        var code = GeneticCode.Standard;
        var invalid = new[]
        {
            new EncOptions { Pseudocount = 1 },
            new EncOptions { BackgroundCorrection = true },
            new EncOptions { WeightedAverage = true },
            new EncOptions { Background = Background },
            new EncOptions { Method = EncMethod.Novembre2002, CodonWCap = true },
            new EncOptions { Method = EncMethod.SunYangXia2013, CodonWCap = true },
            new EncOptions { Method = EncMethod.SunYangXia2013, Pseudocount = -1 },
            new EncOptions { Method = EncMethod.SunYangXia2013, Pseudocount = double.NaN },
            new EncOptions { Method = EncMethod.SunYangXia2013, Pseudocount = double.PositiveInfinity },
            new EncOptions { Method = EncMethod.SunYangXia2013, Background = Background },
            new EncOptions { Method = EncMethod.Novembre2002, BackgroundCorrection = false, Background = Background },
            new EncOptions { Method = (EncMethod)99 },
        };
        Assert.Multiple(() =>
        {
            foreach (var options in invalid)
            {
                Assert.That(() => CodonUsageAnalyzer.TryCalculateEnc(M3, code, options, out _),
                    NUnit.Framework.Throws.ArgumentException.With.Property("ParamName").EqualTo("options"), options.ToString());
                Assert.That(() => CodonUsageAnalyzer.CalculateEnc(M3, code, options), NUnit.Framework.Throws.ArgumentException);
            }
        });
    }

    [Test]
    public void NullArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => CodonUsageAnalyzer.TryCalculateEnc(M3, null!, EncOptions.Default, out _), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => CodonUsageAnalyzer.TryCalculateEnc(M3, GeneticCode.Standard, null!, out _), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => CodonUsageAnalyzer.CalculateEnc(M3, GeneticCode.Standard, (EncOptions)null!), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => CodonUsageAnalyzer.CalculateEnc((DnaSequence)null!, GeneticCode.Standard, EncOptions.Default), NUnit.Framework.Throws.ArgumentNullException);
            Assert.That(() => CodonUsageAnalyzer.TryCalculateEnc((DnaSequence)null!, GeneticCode.Standard, EncOptions.Default, out _), NUnit.Framework.Throws.ArgumentNullException);
        });
    }

    #endregion

    #region CodonW cap (A2-8)

    // CodonW 1.4.4 -enc prints 61.00 for gene M12 under every code (c = 0..7); the sense-codon cap gives
    // 61/60/62/62/62/63/62/62 (NCBI 1, 2, 3, 4, 5, 6, 9, 10).
    [TestCase(1, 61.0, 61.0)]
    [TestCase(2, 60.0, 61.0)]
    [TestCase(3, 62.0, 61.0)]
    [TestCase(4, 62.0, 61.0)]
    [TestCase(5, 62.0, 61.0)]
    [TestCase(6, 63.0, 61.0)]
    [TestCase(9, 62.0, 61.0)]
    [TestCase(10, 62.0, 61.0)]
    public void CodonWCap_OvershootGene_MatchesCodonW(int table, double senseCodonCap, double codonW)
    {
        var code = GeneticCode.GetByTableNumber(table);
        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.CalculateEnc(Over, code), Is.EqualTo(senseCodonCap), "default unchanged");
            Assert.That(CodonUsageAnalyzer.CalculateEnc(Over, code, new EncOptions { CodonWCap = true }), Is.EqualTo(codonW));
        });
    }

    [Test]
    public void CodonWCap_Table2_ValueBetween60And61_MatchesCodonW()
    {
        // r203 under NCBI 2: Wright Nc 60.71 > 60 sense codons; CodonW -code 1 prints 60.71.
        var code = GeneticCode.GetByTableNumber(2);
        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.CalculateEnc(R203, code), Is.EqualTo(60.0));
            Assert.That(CodonUsageAnalyzer.CalculateEnc(R203, code, new EncOptions { CodonWCap = true }),
                Is.EqualTo(60.71285921278657).Within(Tol));
        });
    }

    [Test]
    public void CodonWCap_Table3_ValueBetween61And62_MatchesCodonW()
    {
        // r15 under NCBI 3: Wright Nc 61.89 ≤ 62 sense codons; CodonW -code 2 prints 61.00.
        var code = GeneticCode.GetByTableNumber(3);
        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.CalculateEnc(R15, code), Is.EqualTo(61.89427123695478).Within(Tol));
            Assert.That(CodonUsageAnalyzer.CalculateEnc(R15, code, new EncOptions { CodonWCap = true }), Is.EqualTo(61.0));
        });
    }

    [Test]
    public void CodonWCap_BelowCap_AndNotCalculable_Unchanged()
    {
        var capped = new EncOptions { CodonWCap = true };
        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.CalculateEnc(M3, GeneticCode.Standard, capped), Is.EqualTo(41.288461538461526).Within(Tol));
            Assert.That(CodonUsageAnalyzer.TryCalculateEnc(PheOnly, GeneticCode.Standard, capped, out double enc), Is.False, "CodonW *****");
            Assert.That(enc, Is.EqualTo(0.0));
        });
    }

    #endregion
}
