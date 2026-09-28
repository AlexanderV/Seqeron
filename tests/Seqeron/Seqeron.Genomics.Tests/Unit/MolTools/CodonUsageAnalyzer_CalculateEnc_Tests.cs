// CODON-ENC-001 — Effective Number of Codons (ENC / Nc)
// Evidence: docs/Evidence/CODON-ENC-001-Evidence.md
// TestSpec: tests/TestSpecs/CODON-ENC-001.md
// Source: Wright F (1990). Gene 87(1):23-29; reproduced verbatim in
//         Fuglsang A (2004). Biochem Biophys Res Commun 317:957-964 (Eqs. 1-5a).

using System.Text;

namespace Seqeron.Genomics.Tests.Unit.MolTools;

[TestFixture]
public class CodonUsageAnalyzer_CalculateEnc_Tests
{
    // Standard genetic code (NCBI table 1) sense codons, grouped by amino-acid degeneracy.
    private static readonly string[] AllSenseCodons =
    {
        "TTT","TTC","TTA","TTG","CTT","CTC","CTA","CTG",
        "ATT","ATC","ATA","ATG","GTT","GTC","GTA","GTG",
        "TCT","TCC","TCA","TCG","CCT","CCC","CCA","CCG",
        "ACT","ACC","ACA","ACG","GCT","GCC","GCA","GCG",
        "TAT","TAC","CAT","CAC","CAA","CAG",
        "AAT","AAC","AAA","AAG","GAT","GAC","GAA","GAG",
        "TGT","TGC","TGG","CGT","CGC","CGA","CGG",
        "AGT","AGC","AGA","AGG","GGT","GGC","GGA","GGG"
        // TAA, TAG, TGA (stop) deliberately excluded — not counted in Nc.
    };

    // One representative codon per amino acid (degeneracy-collapsing set).
    private static readonly string[] OneCodonPerAminoAcid =
    {
        "TTT", // Phe
        "CTG", // Leu (6-fold)
        "ATT", // Ile (3-fold)
        "ATG", // Met (single)
        "GTG", // Val (4-fold)
        "TCT", // Ser (6-fold)
        "CCG", // Pro (4-fold)
        "ACC", // Thr (4-fold)
        "GCT", // Ala (4-fold)
        "TAT", // Tyr (2-fold)
        "CAT", // His (2-fold)
        "CAG", // Gln (2-fold)
        "AAT", // Asn (2-fold)
        "AAA", // Lys (2-fold)
        "GAT", // Asp (2-fold)
        "GAA", // Glu (2-fold)
        "TGT", // Cys (2-fold)
        "TGG", // Trp (single)
        "CGT", // Arg (6-fold)
        "GGT"  // Gly (4-fold)
    };

    // Gene M3 (see CalculateEnc_FullyPopulatedBiasedGene_MatchesIndependentReference).
    private static readonly string M3Gene =
        Repeat("TTT", 4) + "TTC" +
        Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
        Repeat("ATT", 3) + Repeat("ATC", 2) + "ATA" +
        Repeat("GTG", 4) + "GTC" +
        Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
        Repeat("CGC", 4) + Repeat("CGT", 2) +
        Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

    // gene = concat over the 64 codons in TCAG order (index i) of codon × ((multiplier·i mod modulus) + 1).
    private static string DeterministicGene(int modulus, int multiplier)
    {
        const string bases = "TCAG";
        var sb = new StringBuilder();
        int i = 0;
        foreach (char b1 in bases)
            foreach (char b2 in bases)
                foreach (char b3 in bases)
                {
                    sb.Append(Repeat(new string(new[] { b1, b2, b3 }), (multiplier * i % modulus) + 1));
                    i++;
                }
        return sb.ToString();
    }

    private static string Repeat(string codon, int times)
    {
        var sb = new StringBuilder(codon.Length * times);
        for (int i = 0; i < times; i++) sb.Append(codon);
        return sb.ToString();
    }

    #region CalculateEnc(string) — canonical

    // M1 — Maximally biased gene: exactly one codon per amino acid (each used twice so
    // F̂ is defined). Wright Eq. (1) gives F̂ = 1 for every class ⇒ N̂c(aa) = 1, and
    // Eq. (3) sums to 9 + 1 + 5 + 3 + 2 = 20 (Fuglsang 2004 extreme-bias limit).
    [Test]
    public void CalculateEnc_OneCodonPerAminoAcid_ReturnsTwenty()
    {
        var seq = string.Concat(OneCodonPerAminoAcid.Select(c => Repeat(c, 2)));

        double enc = CodonUsageAnalyzer.CalculateEnc(seq);

        Assert.That(enc, Is.EqualTo(20.0).Within(1e-9),
            "Extreme bias (one codon per amino acid) must yield Nc = 20 per Wright/Fuglsang Eq. (3).");
    }

    // M2 — Near-uniform usage: every sense codon present in equal counts (2 each).
    // Each class F̂ is well below its asymptotic value at this small count, so raw Eq. (3)
    // overshoots 61 and must be re-adjusted down to exactly 61 (Fuglsang 2004 cap rule).
    [Test]
    public void CalculateEnc_NearUniformUsage_CapsAtSixtyOne()
    {
        var seq = string.Concat(AllSenseCodons.Select(c => Repeat(c, 2)));

        double enc = CodonUsageAnalyzer.CalculateEnc(seq);

        Assert.That(enc, Is.EqualTo(61.0).Within(1e-9),
            "Near-uniform codon usage overshoots Eq. (3) and is re-adjusted to exactly 61.");
    }

    // M3 — Fully-populated biased gene (every degeneracy class estimable, no Eq. 5a / no
    // absent-class fallback). Exercises Wright Eq. (1) homozygosity + Eq. (4) within-class
    // averaging + Eq. (3) aggregation on a realistic input. Codon counts:
    //   Phe(2):  TTTx4 TTCx1        Leu(6):  CTGx3 CTCx2 TTAx1
    //   Ile(3):  ATTx3 ATCx2 ATAx1  Val(4):  GTGx4 GTCx1
    //   Ser(6):  AGCx3 TCTx2 TCAx1  Arg(6):  CGCx4 CGTx2
    //   Gly(4):  GGCx3 GGTx2 GGAx1
    // Per-class average homozygosities (Eq. 1, then class mean):
    //   F̂₂=0.6, F̂₃=0.2666666667, F̂₄=0.4333333333, F̂₆=0.3333333333.
    // Eq. (3): Nc = 2 + 9/0.6 + 1/0.266… + 5/0.433… + 3/0.333… = 41.288461538461526.
    // Expected value is computed by an INDEPENDENT reference implementation of the Wright/
    // codonW algorithm (Python, /tmp/enc_ref.py during validation), NOT read back from the code.
    [Test]
    public void CalculateEnc_FullyPopulatedBiasedGene_MatchesIndependentReference()
    {
        string seq =
            Repeat("TTT", 4) + "TTC" +
            Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
            Repeat("ATT", 3) + Repeat("ATC", 2) + "ATA" +
            Repeat("GTG", 4) + "GTC" +
            Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
            Repeat("CGC", 4) + Repeat("CGT", 2) +
            Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

        double enc = CodonUsageAnalyzer.CalculateEnc(seq);

        Assert.That(enc, Is.EqualTo(41.288461538461526).Within(1e-9),
            "Fully-populated biased gene: F̂₂=0.6, F̂₃=0.2666…, F̂₄=0.4333…, F̂₆=0.3333… (Eq.1) ⇒ "
            + "Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆ = 41.288461538461526 (independent reference).");
    }

    // M4 — Invariant INV-01: 20 ≤ Nc ≤ 61 whenever Nc is calculable (property test on
    // deterministic genes in which every synonymous class is estimable).
    // (Validation 2026-09: the former cases were genes with empty synonymous classes, e.g.
    // Lys-only "AAA…", for which CodonW enc_out does not calculate Nc — they now return 0,
    // see CalculateEnc_EmptySynonymousClass_NotCalculated.)
    [TestCase(7, 3)]
    [TestCase(5, 7)]
    [TestCase(11, 3)]
    [TestCase(13, 5)]
    [TestCase(9, 2)]
    public void CalculateEnc_AnyValidSequence_StaysWithinRange(int modulus, int multiplier)
    {
        string seq = DeterministicGene(modulus, multiplier);

        double enc = CodonUsageAnalyzer.CalculateEnc(seq);

        Assert.Multiple(() =>
        {
            Assert.That(enc, Is.GreaterThanOrEqualTo(20.0),
                "INV-01: Nc cannot fall below the extreme-bias limit 20.");
            Assert.That(enc, Is.LessThanOrEqualTo(61.0),
                "INV-01: Nc cannot exceed 61 (Eq. 3 re-adjustment).");
        });
    }

    // M5 — Isoleucine (the only 3-fold amino acid) absent, but the 2-, 4- and 6-fold classes
    // are all populated, so Wright Eq. (5a) F̂₃ = (F̂₂ + F̂₄)/2 genuinely fires (codonW: "an
    // exception is made … isoleucine 3-fold absent ⇒ F̂₃ = average of F̂₂ and F̂₄", Peden thesis).
    // Same gene as M3 with Ile removed; counts:
    //   Phe(2): TTTx4 TTCx1   Leu(6): CTGx3 CTCx2 TTAx1   Val(4): GTGx4 GTCx1
    //   Ser(6): AGCx3 TCTx2 TCAx1   Arg(6): CGCx4 CGTx2   Gly(4): GGCx3 GGTx2 GGAx1
    // F̂₂=0.6, F̂₄=0.4333333333, F̂₆=0.3333333333; Eq. (5a) ⇒ F̂₃=(0.6+0.4333…)/2=0.5166666667.
    // Eq. (3): Nc = 2 + 9/F̂₂ + 1/F̂₃ + 5/F̂₄ + 3/F̂₆ = 39.47394540942927 (independent reference).
    [Test]
    public void CalculateEnc_IsoleucineAbsent_UsesEq5aFallback()
    {
        string seq =
            Repeat("TTT", 4) + "TTC" +
            Repeat("CTG", 3) + Repeat("CTC", 2) + "TTA" +
            Repeat("GTG", 4) + "GTC" +
            Repeat("AGC", 3) + Repeat("TCT", 2) + "TCA" +
            Repeat("CGC", 4) + Repeat("CGT", 2) +
            Repeat("GGC", 3) + Repeat("GGT", 2) + "GGA";

        double enc = CodonUsageAnalyzer.CalculateEnc(seq);

        Assert.That(enc, Is.EqualTo(39.47394540942927).Within(1e-9),
            "Ile absent, F̂₂/F̂₄/F̂₆ estimable ⇒ Eq. (5a) F̂₃=(F̂₂+F̂₄)/2=0.5166…; "
            + "Nc = 2 + 9/0.6 + 1/0.5166… + 5/0.4333… + 3/0.3333… = 39.47394540942927 (independent reference).");
    }

    // M5b — Whole synonymous class absent ⇒ Nc NOT calculated (returns 0).
    // Source: CodonW README_indices.txt (Peden 1999): "When there are no amino acids in a
    // synonymous family, Nc is not calculated as the gene is either too short or has extremely
    // skewed amino acid usage (Wright 1990). An exception to this is made for genetic codes where
    // isoleucine is the only 3-fold synonymous amino acid". Confirmed with the compiled CodonW
    // 1.4.4 binary (`-enc`), which prints "*****" for every case below.
    // Validation 2026-09 (finding F15): this test previously pinned a library-specific 29.0 for
    // the Phe-only gene (absent classes contributing their full codon count) — corrected to the
    // sourced behaviour.
    [TestCase("TTTTTTTTTTTC", TestName = "CalculateEnc_EmptySynonymousClass_PheOnly_NotCalculated")]
    [TestCase("TTTTTTTTTTTCCTGCTGCTGCTCCTC", TestName = "CalculateEnc_EmptySynonymousClass_IleAnd4FoldAbsent_NotCalculated")]
    [TestCase("AAAAAAAAAAAAAAAAAA", TestName = "CalculateEnc_EmptySynonymousClass_LysOnly_NotCalculated")]
    [TestCase("ATGAAAGAGCTGTTCGCCAAA", TestName = "CalculateEnc_EmptySynonymousClass_ShortGene_NotCalculated")]
    [TestCase("ATGTGGATGTGG", TestName = "CalculateEnc_EmptySynonymousClass_MetTrpOnly_NotCalculated")]
    [TestCase("A", TestName = "CalculateEnc_EmptySynonymousClass_NoCompleteCodon_NotCalculated")]
    public void CalculateEnc_EmptySynonymousClass_NotCalculated(string seq)
    {
        Assert.That(CodonUsageAnalyzer.CalculateEnc(seq), Is.EqualTo(0.0),
            "CodonW enc_out: a synonymous class with no estimable amino acid (other than the "
            + "single 3-fold Ile class) ⇒ Nc not calculated (\"*****\"); the library returns 0.");
    }

    // F̂ = 0 (every observed codon of the amino acid used exactly once) is not an estimate:
    // CodonW enc_out adds an amino acid to its class average only "if (bb > 0.0000001)".
    // Gene M3 plus His as CAT+CAC (n=2, Σp²=0.5 ⇒ F̂ = 0): His is left out, so Nc equals M3's
    // 41.288461538461526 (CodonW 1.4.4 binary: 41.29 for both genes). Before the 2026-09 fix His
    // lowered F̄₂ to (0.6+0)/2 = 0.3 and gave 9/0.3 = 30 for the 2-fold class (Nc = 56.29).
    [Test]
    public void CalculateEnc_ZeroHomozygosityAminoAcid_ExcludedFromClassAverage()
    {
        double enc = CodonUsageAnalyzer.CalculateEnc(M3Gene + "CATCAC");

        Assert.That(enc, Is.EqualTo(41.288461538461526).Within(1e-9));
    }

    // Every sense codon exactly once: every F̂ = 0, so no class is estimable and CodonW prints
    // "*****". (The pre-2026-09 code returned 20 — "extreme bias" — for this maximally even gene.)
    [Test]
    public void CalculateEnc_AllSenseCodonsOnce_NotCalculated()
    {
        Assert.That(CodonUsageAnalyzer.CalculateEnc(string.Concat(AllSenseCodons)), Is.EqualTo(0.0));
    }

    // Standard-code gene with every class estimable: gene(i) = codon_i × ((3·i mod 7) + 1) over the
    // 64 codons in TCAG order. CodonW 1.4.4 `-enc -code 0` prints 57.56; the Python port of
    // enc_out and codonbias 0.5.0 EffectiveNumberOfCodons(robust=False, pseudocount=0,
    // mean='unweighted') give 57.5614857446809.
    [Test]
    public void CalculateEnc_DeterministicGene_MatchesCodonW()
    {
        Assert.That(CodonUsageAnalyzer.CalculateEnc(DeterministicGene(7, 3)),
            Is.EqualTo(57.5614857446809).Within(1e-9));
    }

    #region CalculateEnc(string, GeneticCode) — genetic-code-aware classes

    // Synonymous classes come from the genetic code (CodonW `-enc` honours `-code`).
    // References: CodonW 1.4.4 binary (2 dp) / enc_out port and codonbias 0.5.0 (Wright mode):
    //   table 2 (CodonW code 1, vertebrate mito: 12 two-fold, 6 four-fold, 2 six-fold, no 1/3-fold)
    //   table 3 (CodonW code 2, yeast mito: 8-fold Thr = ACN + CTN), table 9 (CodonW code 7,
    //   echinoderm mito: Ile and Asn both 3-fold).
    [TestCase(1, 57.5614857446809, TestName = "CalculateEnc_Table1_DeterministicGene_CodonW_57_56")]
    [TestCase(2, 55.38187053526956, TestName = "CalculateEnc_Table2_DeterministicGene_CodonW_55_38")]
    [TestCase(3, 57.4784046236798, TestName = "CalculateEnc_Table3_DeterministicGene_CodonW_57_48")]
    [TestCase(9, 58.152584900842236, TestName = "CalculateEnc_Table9_DeterministicGene_CodonW_58_15")]
    public void CalculateEnc_AlternativeGeneticCode_MatchesCodonW(int table, double expected)
    {
        double enc = CodonUsageAnalyzer.CalculateEnc(DeterministicGene(7, 3), GeneticCode.GetByTableNumber(table));

        Assert.That(enc, Is.EqualTo(expected).Within(1e-9));
    }

    // Gene M3 under other codes (CodonW 1.4.4: table 2 → 45.00, table 3 → 43.19, table 9 → 43.83).
    [TestCase(2, 45.0)]
    [TestCase(3, 43.18531468531467)]
    [TestCase(9, 43.83333333333333)]
    public void CalculateEnc_M3Gene_AlternativeGeneticCode_MatchesCodonW(int table, double expected)
    {
        Assert.That(CodonUsageAnalyzer.CalculateEnc(M3Gene, GeneticCode.GetByTableNumber(table)),
            Is.EqualTo(expected).Within(1e-9));
    }

    // Overshoot is re-adjusted down to the number of sense codons of the code (Wright 1990:
    // uniform usage gives Nc = Σ_z K_z·z = number of sense codons; 61 in table 1, CodonW's cap).
    // codonbias 0.5.0 caps at the same value: gene(i) = codon_i × ((7·i mod 5) + 1) → 61 (table 1),
    // 60 (table 2), 62 (tables 3, 9). CodonW 1.4.4 hard-codes 61 for every code (documented divergence
    // for codes whose sense-codon count is not 61).
    [TestCase(1, 61.0)]
    [TestCase(2, 60.0)]
    [TestCase(3, 62.0)]
    [TestCase(9, 62.0)]
    public void CalculateEnc_Overshoot_CappedAtSenseCodonCount(int table, double expected)
    {
        Assert.That(CodonUsageAnalyzer.CalculateEnc(DeterministicGene(5, 7), GeneticCode.GetByTableNumber(table)),
            Is.EqualTo(expected));
    }

    // RNA spelling and lower case give the same Nc (CodonW ident_codon reads U as T).
    [Test]
    public void CalculateEnc_RnaAndLowerCase_EqualDna()
    {
        string dna = DeterministicGene(7, 3);
        double expected = CodonUsageAnalyzer.CalculateEnc(dna);

        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.CalculateEnc(dna.Replace('T', 'U')), Is.EqualTo(expected));
            Assert.That(CodonUsageAnalyzer.CalculateEnc(dna.ToLowerInvariant().Replace('t', 'u')), Is.EqualTo(expected));
        });
    }

    [Test]
    public void CalculateEnc_GeneticCodeOverloads_NullArguments_Throw()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentNullException>(() => CodonUsageAnalyzer.CalculateEnc("ATG", null!));
            Assert.Throws<ArgumentNullException>(() => CodonUsageAnalyzer.CalculateEnc(new DnaSequence("ATG"), null!));
            Assert.Throws<ArgumentNullException>(() => CodonUsageAnalyzer.CalculateEnc((DnaSequence)null!, GeneticCode.Standard));
        });
    }

    [Test]
    public void CalculateEnc_DnaSequenceWithCode_EqualsStringOverload()
    {
        string seq = DeterministicGene(7, 3);
        var code = GeneticCode.GetByTableNumber(2);

        Assert.That(CodonUsageAnalyzer.CalculateEnc(new DnaSequence(seq), code),
            Is.EqualTo(CodonUsageAnalyzer.CalculateEnc(seq, code)));
    }

    #endregion

    // M7 — Empty / null string returns 0 (degenerate input contract).
    [Test]
    public void CalculateEnc_EmptyString_ReturnsZero()
    {
        Assert.Multiple(() =>
        {
            Assert.That(CodonUsageAnalyzer.CalculateEnc(""), Is.EqualTo(0.0),
                "Empty sequence has no codons; contract returns 0.");
            Assert.That(CodonUsageAnalyzer.CalculateEnc((string)null!), Is.EqualTo(0.0),
                "Null string is treated as empty; contract returns 0.");
        });
    }

    // S1 — Case insensitivity: lowercase input is normalized to upper case.
    [Test]
    public void CalculateEnc_LowercaseInput_EqualsUppercase()
    {
        string seq = Repeat("TTT", 3) + "TTC" + Repeat("GCT", 2) + Repeat("GCC", 2);

        double upper = CodonUsageAnalyzer.CalculateEnc(seq);
        double lower = CodonUsageAnalyzer.CalculateEnc(seq.ToLowerInvariant());

        Assert.That(lower, Is.EqualTo(upper).Within(1e-12),
            "Lowercase input must be normalized and produce the identical Nc.");
    }

    // S2 — Codons containing non-ACGT characters are skipped (consistent with CountCodons).
    [Test]
    public void CalculateEnc_InvalidCodonsSkipped_EqualsCleanSequence()
    {
        string clean = Repeat("TTT", 3) + "TTC";
        string withN = clean + "NNN" + "TAN"; // two non-ACGT codons appended

        double encClean = CodonUsageAnalyzer.CalculateEnc(clean);
        double encWithN = CodonUsageAnalyzer.CalculateEnc(withN);

        Assert.That(encWithN, Is.EqualTo(encClean).Within(1e-12),
            "Non-ACGT codons are skipped, so they must not change Nc.");
    }

    #endregion

    #region CalculateEnc(DnaSequence) — delegate

    // M6 — Null DnaSequence overload throws ArgumentNullException (contract).
    [Test]
    public void CalculateEnc_NullDnaSequence_Throws()
    {
        Assert.Throws<ArgumentNullException>(() => CodonUsageAnalyzer.CalculateEnc((DnaSequence)null!),
            "Null DnaSequence must raise ArgumentNullException per the contract.");
    }

    // S3 — DnaSequence overload delegates to the string overload (identical result).
    [Test]
    public void CalculateEnc_DnaSequenceOverload_DelegatesToStringOverload()
    {
        string seq = Repeat("TTT", 3) + "TTC" + Repeat("GCT", 2) + Repeat("GCC", 2);

        double viaString = CodonUsageAnalyzer.CalculateEnc(seq);
        double viaDna = CodonUsageAnalyzer.CalculateEnc(new DnaSequence(seq));

        Assert.That(viaDna, Is.EqualTo(viaString).Within(1e-12),
            "DnaSequence overload must delegate to the string computation and return the same Nc.");
    }

    #endregion

    #region Intermediate-bias exact derivation (COULD)

    // C1 — Every multi-codon amino acid uses its first two codons at a fixed 2:1 ratio
    // (counts 2 and 1, n=3); singlets contribute count 1. By Eq. (1) each represented
    // amino acid then has Σp² = (2/3)² + (1/3)² = 5/9 and F̂ = (3·5/9 − 1)/(3 − 1) = 1/3,
    // so the average homozygosity of every degeneracy class is exactly 1/3. By Eq. (3):
    //   Nc = 2 + 9/(1/3) + 1/(1/3) + 5/(1/3) + 3/(1/3) = 2 + 27 + 3 + 15 + 9 = 56.
    // This exercises an exact intermediate value (no clamping, all classes estimable).
    [Test]
    public void CalculateEnc_TwoToOneBiasAllAminoAcids_EqualsFiftySix()
    {
        var sb = new StringBuilder();
        foreach (var codons in MultiCodonFamilies)
        {
            sb.Append(Repeat(codons[0], 2)); // major codon, count 2
            sb.Append(codons[1]);            // minor codon, count 1
        }
        sb.Append("ATG"); // Met (singlet)
        sb.Append("TGG"); // Trp (singlet)

        double enc = CodonUsageAnalyzer.CalculateEnc(sb.ToString());

        Assert.That(enc, Is.EqualTo(56.0).Within(1e-9),
            "Uniform 2:1 bias gives F̂ = 1/3 in every class; Eq. (3) yields Nc = 2 + 27 + 3 + 15 + 9 = 56.");
    }

    #endregion

    // First two codons of each multi-codon amino acid family (standard genetic code).
    private static readonly string[][] MultiCodonFamilies =
    {
        new[]{"TTT","TTC"}, // Phe
        new[]{"TTA","TTG"}, // Leu (6-fold)
        new[]{"ATT","ATC"}, // Ile (3-fold)
        new[]{"GTT","GTC"}, // Val
        new[]{"TCT","TCC"}, // Ser (6-fold)
        new[]{"CCT","CCC"}, // Pro
        new[]{"ACT","ACC"}, // Thr
        new[]{"GCT","GCC"}, // Ala
        new[]{"TAT","TAC"}, // Tyr
        new[]{"CAT","CAC"}, // His
        new[]{"CAA","CAG"}, // Gln
        new[]{"AAT","AAC"}, // Asn
        new[]{"AAA","AAG"}, // Lys
        new[]{"GAT","GAC"}, // Asp
        new[]{"GAA","GAG"}, // Glu
        new[]{"TGT","TGC"}, // Cys
        new[]{"CGT","CGC"}  // Arg (6-fold)
    };
}
