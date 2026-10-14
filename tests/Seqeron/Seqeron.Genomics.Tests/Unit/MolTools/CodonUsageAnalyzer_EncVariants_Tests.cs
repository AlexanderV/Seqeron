// CODON-ENC-001 (audit A2-7, A2-8, A3-2) — opt-in Nc estimators and the CodonW cap option:
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

    #region ENCprime (A3-2)

    // ENCprime (Novembre 2002 reference program; github.com/jnovembre/ENCprime @ 709ace50, src/ENCprime.c with
    // the JN 2/28/06 3-fold fix) compiled locally and run as `ENCprime <codcnt> <acgtfreq> <code> out 0 -q`.
    // It prints Nc / Ncp with %.4f, so values are locked to ±5e-5 of the printed number.
    private const double EncPrimeTol = 5e-5;

    // example/ExSeqs.fasta and example/ExBackgroundSeqs.fasta of the ENCprime repository (shipped results:
    // example/ExResults, Nc/Ncp computed with the background file as the SeqCount -n input).
    private const string ExHmaif =
        "ATGAACTGTGTTTGCCGCCTGGTCCTGGTCGTGCTGAGCCTGTGGCCAGATACAGCTGTCGCCCCTGGGCCACCACCTGGCCCCCCTCGAGTTTCCCCAGACCCTCGGGC" +
        "CGAGCTGGACAGCACCGTGCTCCTGACCCGCTCTCTCCTGGCGGACACGCGGCAGCTGGCTGCACAGCTGAGGGACAAATTCCCAGCTGACGGGGACCACAACCTGGATT" +
        "CCCTGCCCACCCTGGCCATGAGTGCGGGGGCACTGGGAGCTCTACAGCTCCCAGGTGTGCTGACAAGGCTGCGAGCGGACCTACTGTCCTACCTGCGGCACGTGCAGTGG" +
        "CTGCGCCGGGCAGGTGGCTCTTCCCTGAAGACCCTGGAGCCCGAGCTGGGCACCCTGCAGGCCCGACTGGACCGGCTGCTGCGCCGGCTGCAGCTCCTGATGTCCCGCCT" +
        "GGCCCTGCCCCAGCCACCCCCGGACCCGCCGGCGCCCCCGCTGGCGCCCCCCTCCTCAGCCTGGGGGGGCATCAGGGCCGCCCACGCCATCCTGGGGGGGCTGCACCTGA" +
        "CACTTGACTGGGCCGTGAGGGGACTGCTGCTGCTGAAGACTCGGCTGTGA";
    private const string ExHmaifBackground =
        "AGTAGAGACGGGGTTTCACCATTCAAGACCGGCCAGGCTGGTCTTGAACTCCTGACCTCCTGAGCCACCCACCTTGGCCTCTCAAAGTGTTGGGATTACAGGCATGAGCC" +
        "ACTGCGCCCAGCCTACCTGTTTATTTTTATGTGTTCTTTAGTACAATGTAAATTTCCATGAGGATAAGGACCTTACATGTCTTCCTCATGATCACTACATTGATATGTGT" +
        "GTTGTGTGTGTGTGCAAAAGGGGAGCTGGAAAAACAGTTGAGTCCAAAATACAGGTAAAATATTGCACTAAAGAAAAACGAGGCTGGGTGTAGTGGCTCACAACTGTAAT" +
        "CCCAGCATTTTGGGAGGCTGAGGCAGGCAGATCACGAGGTCAGGAGTTCAAGACCAGCCTGGCCAACATGGTGAAACCGTGTCTCTACTAAAAATACAAAAATTAGCTGG" +
        "GCACGGTGGCACACGCCTGTAATTCCAGCTACTCGGGAGGCTGAGGCAGGAGAATCACCTGAACCCAGGAGGTGGAGGTTGCAATAAGCGGAGATCAAGCCACTGCACTC" +
        "TGGCCCGGGCGACAGAGCAAGACTCTGTCTCGGAAAAAAAAAAAAAAAAAAAAAAAGAAAAGAAAAAAGAAAGCCCCATCACTGAAAACAGGTTTGAAATATAAAAGAAA" +
        "TAAGGTGATTGGTAAATGGGGTAAAGTTGAAGCTTCTTCCATAGCTTATGCTACGATCTTTATTCCCAACCAGGTTGTAGGCCTTTGCTCCACAAGTCTGGGTGCCTTGG" +
        "AGGAAGGTTTCACAAACAACACAAAAATGTCTCCTTTCTCTTACCACATGGGCTGTGGCCGGCCCTGATGGTGTCCCCTGCCTTCCTTACCTGTGAACTGTGTTTCTTCC" +
        "CTCAAAATCATCCAGGGCTCTTTTCCTCGCTCCAGTAAATCAACCACAACTGGTTTAGAAATGGAATGTCCTCCTTATAAGAAAAGAAATGCCACATGGTATAAAAATAA" +
        "AAACAAAAGT";
    private const string ExHs14kdapt =
        "ATGTCGTCCTTGATCAGAAGGGTGATCAGCACCGCGAAAGCCCCAGGGGCCATTGGACCCTACAGTCAAGCTGTATTAGTCGACAGGACCATTTACATTTCAGGACAGAT" +
        "AGGCATGGACCCTTCAAGTGGACAGCTTGTGTCAGGAGGGGTAGCAGAAGAAGCTAAACAAGCTCTTAAAAACATGGGTGAAATTCTGAAAGCTGCAGGCTGTGACTTCA" +
        "CTAACGTGGTGAAAACAACTGTTCTTCTGGCTGACATAAATGACTTCAATACTGTCAATGAAATCTACAAACAGTATTTCAAGAGTAATTTTCCTGCTAGAGCTGCTTAC" +
        "CAAGTTGCTGCTTTACCCAAAGGCAGCCGAATTGAAATTGAAGCAGTAGCTATCCAAGGACCACTGACAACGGCATCACTATAA";
    private const string ExHs14kdaptBackground =
        "TCCCAACCTCTGCCCTCTGACACACCCCAGTGTGTGTTGTTCCCCTCTATGTGTCCATGAGTCCTCATCATTTACCTTCCACTTATAAGTAAGAACATGTGGTATTTGGT" +
        "TTTCTGTTCCTGCATTAGTTTGCTAAGGATAGTGGCCTCCAGCTCCATCCACATCCCTGAAAAGAACATAATCTCATTCTTTTTTGTGGCTGCATAATATTCCATGGTAT" +
        "ATATGTACCACATTTTCTTTATCCAGTCTATCATTGGTGGACATTTGGGTTGATTTCATGTGTTTGCTATTGTGAATAGTTGAAATGGCATATCTTTCATTGTTCATTAA" +
        "CTAAATTTATGATTTTTATTTTTAAAACACAAGTATTTTCACTTGTGTAAAAAACATTTTTGTCCTTATTTTTTTGGAAAAAACTTTCCAAATGAAATAAATCAAATGAA" +
        "TGCAGATTTCTTAAATAATCTGTGTTTTTTTGTAAAACAGAAATCTATTTTAATGTGTTGAAGCAGAACTTTATCAAACTAGTGCTTCCAGTAAGATAATTTCACAACTG" +
        "ACAGGTTAAATGAAGAGGAATTCTGACTGATGGCAAGTCAGAAAACTTTGGTCATTTATCTTTCCTTTCACACAGATGGCTCTAATTCCAGGATAAAATGACTGCTGGCA" +
        "GTAGGAAGAGGAAGACAAAAAGAGACAAGAAATTGAGATACTTCAAATGTCAATGGGACTGGAAAGGATTTGATAAAGGTGGGTGGCCGCTGTTCCCCTACCTTTGGGCT" +
        "ACGATATTTCCTTTGGAAAGGAGGATTTGGCTGAAAGGTGGTGAAGTAAAATTCCTCAATACCTCAGTTAATGAGCACAGTGCCTTTCCTGGAGGGTCTCTGTGGTATTC" +
        "TTCCTATTCTTTCTTATACTTCTGTATTACATGTTGCCTCACCATCAGTTTCCTCATCAGAACATTTCCCCATCTCAACTCATACTTCTGTACTAGCTCTCATCCTATTC" +
        "TCACTTTCTA";

    private static double EncPrime(string gene, int table, EncMethod method, string? background = null)
    {
        var options = new EncOptions { Method = method, Background = background };
        bool ok = CodonUsageAnalyzer.TryCalculateEnc(gene, GeneticCode.GetByTableNumber(table), options, out double enc);
        Assert.That(ok, Is.True, "ENCprime reports a value for every non-empty gene");
        return enc;
    }

    // example/ExResults (shipped with ENCprime): "HMAIF.PE1: 40.0709 36.2391 …", "HS14KDAPT.PE1: 56.3938 57.2395 …".
    [TestCase("HMAIF", 40.0709, 36.2391)]
    [TestCase("HS14KDAPT", 56.3938, 57.2395)]
    public void EncPrime_ShippedExample_MatchesExResults(string gene, double nc, double ncPrime)
    {
        (string seq, string bg) = gene == "HMAIF" ? (ExHmaif, ExHmaifBackground) : (ExHs14kdapt, ExHs14kdaptBackground);
        Assert.Multiple(() =>
        {
            Assert.That(EncPrime(seq, 1, EncMethod.EncPrimeNc), Is.EqualTo(nc).Within(EncPrimeTol));
            Assert.That(EncPrime(seq, 1, EncMethod.EncPrimeNcPrime, bg), Is.EqualTo(ncPrime).Within(EncPrimeTol));
        });
    }

    // ENCprime with the gene itself as the nucleotide-composition input (SeqCount -n on the same FASTA).
    [TestCase("HMAIF", 1, 40.0709, 51.1290)]
    [TestCase("HMAIF", 2, 39.4889, 50.9255)]
    [TestCase("HS14KDAPT", 1, 56.3938, 58.0541)]
    [TestCase("HS14KDAPT", 2, 55.1371, 57.4080)]
    [TestCase("M3", 1, 50.7500, 48.9509)]
    [TestCase("M3", 2, 47.8636, 45.5472)]
    [TestCase("Det", 1, 58.5287, 58.3944)]
    [TestCase("Det", 2, 56.6344, 56.5795)]
    [TestCase("Over", 1, 60.5455, 60.5030)]
    [TestCase("Over", 2, 60.0000, 60.0000)]
    [TestCase("R53", 1, 61.0000, 61.0000)] // every amino acid n ≤ 5: all classes imputed with 1/k
    [TestCase("R203", 1, 56.7452, 57.5932)]
    [TestCase("R203", 2, 56.4075, 56.7819)]
    [TestCase("R203", 11, 56.7452, 57.5932)]
    [TestCase("R15", 1, 55.4554, 56.4604)]
    [TestCase("R15", 2, 55.3108, 56.5900)]
    [TestCase("Phe", 1, 61.0000, 61.0000)]
    public void EncPrime_GeneAsBackground_MatchesProgram(string gene, int table, double nc, double ncPrime)
    {
        string seq = gene switch { "HMAIF" => ExHmaif, "HS14KDAPT" => ExHs14kdapt, _ => Gene(gene) };
        Assert.Multiple(() =>
        {
            Assert.That(EncPrime(seq, table, EncMethod.EncPrimeNc), Is.EqualTo(nc).Within(EncPrimeTol));
            Assert.That(EncPrime(seq, table, EncMethod.EncPrimeNcPrime), Is.EqualTo(ncPrime).Within(EncPrimeTol));
        });
    }

    // Separate background (nucleotide-composition file from the BG sequence); Nc is unaffected.
    [TestCase("M3", 1, 46.5174)]
    [TestCase("M3", 2, 46.2959)]
    [TestCase("Det", 1, 42.4160)]
    [TestCase("Over", 1, 48.2205)]
    [TestCase("Over", 2, 47.7853)]
    [TestCase("R203", 1, 53.5760)]
    [TestCase("R15", 1, 50.3005)]
    [TestCase("R15", 2, 48.6160)]
    public void EncPrime_ExplicitBackground_MatchesProgram(string gene, int table, double ncPrime)
    {
        Assert.That(EncPrime(Gene(gene), table, EncMethod.EncPrimeNcPrime, Background), Is.EqualTo(ncPrime).Within(EncPrimeTol));
    }

    // Rule probes (ENCprime values; first column Nc, second Ncp; table 1 / table 2).
    [TestCase("TTTTTTTTTTTTTTTTTC", 1, 56.5000, 59.6304)]  // Phe n = 6 > 5: F = 2/3 enters; F̄₃ = 1/3 (no 2+4 imputation)
    [TestCase("TTTTTTTTTTTTTTTTTC", 2, 54.0000, 58.1739)]
    [TestCase("TTTTTTTTTTTTTTC", 1, 61.0000, 61.0000)]     // Phe n = 5: not estimated (threshold n > 5)
    [TestCase("TTTTTTTTTTTTTTTTTCGCTGCTGCTGCTGCTGCTGCCGCC", 1, 43.8654, 61.0000)] // no Ile: F̄₃ = (F̄₂ + F̄₄)/2
    [TestCase("TTTTTTTTTTTTTTTTTCGCTGCTGCTGCTGCTGCTGCCGCC", 2, 40.5000, 60.0000)]
    [TestCase("TTTTTTTTTTTCTTCTTCTATTATTATTACTACTAC", 1, 61.0000, 55.2034)] // F = 0.4 < 1/2: class mean floored at 1/2
    [TestCase("ATGATGATGATGATGATGATGATGTGGTGGTGGTGGTGGTGGTGGTGG", 2, 48.0000, 55.3846)] // table 2: Met/Trp 2-fold
    [TestCase("TAATAGTGA", 1, 61.0000, 61.0000)] // no sense codon: uninformative Σ K_z·z
    public void EncPrime_RuleProbes_MatchProgram(string gene, int table, double nc, double ncPrime)
    {
        Assert.Multiple(() =>
        {
            Assert.That(EncPrime(gene, table, EncMethod.EncPrimeNc), Is.EqualTo(nc).Within(EncPrimeTol));
            Assert.That(EncPrime(gene, table, EncMethod.EncPrimeNcPrime), Is.EqualTo(ncPrime).Within(EncPrimeTol));
        });
    }

    [Test]
    public void EncPrime_RuleProbe_HandDerivations()
    {
        // Phe TTT×5 TTC×1 (n = 6): F = (6·26/36 − 1)/5 = 2/3; Nc = 2 + 9·3/2 + 3 + 5·4 + 3·6 = 56.5.
        Assert.That(EncPrime("TTTTTTTTTTTTTTTTTC", 1, EncMethod.EncPrimeNc), Is.EqualTo(56.5).Within(1e-12));
        // Floor: Phe 3/3 and Tyr 3/3 (n = 6): F = (6·½ − 1)/5 = 0.4 < ½ → F̄₂ = ½ → Nc = 61.
        Assert.That(EncPrime("TTTTTTTTTTTCTTCTTCTATTATTATTACTACTAC", 1, EncMethod.EncPrimeNc), Is.EqualTo(61.0));
    }

    [Test]
    public void EncPrime_ExpectedBelowEps_UsesProgramSentinel()
    {
        // Background without T: e(GCT) = 0 < 1e-6, so test_null resets χ² to −1 and keeps accumulating
        // GCC, GCA, GCG (e = 1/3 each): χ² = −1 + 6·(4/3 + 1/3 + 1/3) = 11; F′ = (11 + 6 − 4)/(4·5) = 0.65;
        // Nc′ = 2 + 18 + 3 + 5/0.65 + 18 = 48.6923 (ENCprime prints 48.6923; Nc 46.0000).
        const string gene = "GCCGCCGCCGCCGCCGCC";
        Assert.Multiple(() =>
        {
            Assert.That(EncPrime(gene, 1, EncMethod.EncPrimeNcPrime, "CCCAAAGGG"), Is.EqualTo(41.0 + 5.0 / 0.65).Within(1e-12));
            Assert.That(EncPrime(gene, 1, EncMethod.EncPrimeNcPrime, "CCCAAAGGG"), Is.EqualTo(48.6923).Within(EncPrimeTol));
            Assert.That(EncPrime(gene, 1, EncMethod.EncPrimeNc), Is.EqualTo(46.0).Within(1e-12));
            Assert.That(EncPrime(gene + "TTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGG", 1, EncMethod.EncPrimeNcPrime, "CCCAAAGGG"),
                Is.EqualTo(38.8784).Within(EncPrimeTol));
        });
    }

    [Test]
    public void EncPrime_Table13_DeviatesFromProgramParsingBug()
    {
        // ENCprime get_code sets n_aa to the index of the last codon's amino acid (GGG → G); in NCBI table 13
        // G is first met at AGA, so V, A, D, E fall outside n_aa and are dropped from Nc (M3: 37.6786 / 35.1829).
        // The library keeps all 21 amino acids: values equal ENCprime with n_aa = number of distinct letters.
        Assert.Multiple(() =>
        {
            Assert.That(EncPrime(M3, 13, EncMethod.EncPrimeNc), Is.EqualTo(45.9643).Within(EncPrimeTol));
            Assert.That(EncPrime(M3, 13, EncMethod.EncPrimeNcPrime), Is.EqualTo(43.8972).Within(EncPrimeTol));
        });
    }

    [Test]
    public void EncPrime_RnaLowerCaseAndDnaSequence_Agree()
    {
        foreach (var method in new[] { EncMethod.EncPrimeNc, EncMethod.EncPrimeNcPrime })
        {
            double expected = EncPrime(M3, 1, method);
            Assert.Multiple(() =>
            {
                Assert.That(EncPrime(M3.Replace('T', 'U').ToLowerInvariant(), 1, method), Is.EqualTo(expected).Within(1e-12));
                Assert.That(CodonUsageAnalyzer.CalculateEnc(new DnaSequence(M3), GeneticCode.Standard, new EncOptions { Method = method }),
                    Is.EqualTo(expected));
            });
        }
    }

    [Test]
    public void EncPrime_EmptyInput_NotCalculable_AndInvalidOptionsThrow()
    {
        Assert.Multiple(() =>
        {
            foreach (var method in new[] { EncMethod.EncPrimeNc, EncMethod.EncPrimeNcPrime })
            {
                Assert.That(CodonUsageAnalyzer.TryCalculateEnc("", GeneticCode.Standard, new EncOptions { Method = method }, out double enc), Is.False);
                Assert.That(enc, Is.EqualTo(0.0));
                foreach (var bad in new[]
                {
                    new EncOptions { Method = method, CodonWCap = true },
                    new EncOptions { Method = method, Pseudocount = 0 },
                    new EncOptions { Method = method, BackgroundCorrection = true },
                    new EncOptions { Method = method, WeightedAverage = false },
                })
                    Assert.That(() => CodonUsageAnalyzer.TryCalculateEnc(M3, GeneticCode.Standard, bad, out _),
                        NUnit.Framework.Throws.ArgumentException.With.Property("ParamName").EqualTo("options"), bad.ToString());
            }
            Assert.That(() => CodonUsageAnalyzer.TryCalculateEnc(M3, GeneticCode.Standard,
                new EncOptions { Method = EncMethod.EncPrimeNc, Background = Background }, out _),
                NUnit.Framework.Throws.ArgumentException.With.Property("ParamName").EqualTo("options"));
        });
    }

    #endregion
}
