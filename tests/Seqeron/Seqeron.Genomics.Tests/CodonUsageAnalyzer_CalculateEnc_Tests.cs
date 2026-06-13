// CODON-ENC-001 — Effective Number of Codons (ENC / Nc)
// Evidence: docs/Evidence/CODON-ENC-001-Evidence.md
// TestSpec: tests/TestSpecs/CODON-ENC-001.md
// Source: Wright F (1990). Gene 87(1):23-29; reproduced verbatim in
//         Fuglsang A (2004). Biochem Biophys Res Commun 317:957-964 (Eqs. 1-5a).

using System;
using System.Linq;
using System.Text;
using NUnit.Framework;
using Seqeron.Genomics;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.Tests;

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

    // M3 — Single two-fold amino acid Phe (TTT x3, TTC x1). Hand derivation by Eq. (1):
    // n=4, p=(0.75,0.25), Σp²=0.625, F̂=(4*0.625-1)/3=0.5, so 9/F̂₂=18. The 3-,4-,6-fold
    // classes have no estimable amino acid and contribute their full counts 1+5+3. Met+Trp=2.
    // Nc = 2 + 18 + 1 + 5 + 3 = 29.0.
    [Test]
    public void CalculateEnc_SinglePhenylalanineTwoFold_MatchesHandDerivation()
    {
        string seq = Repeat("TTT", 3) + "TTC"; // F̂₂ = 0.5

        double enc = CodonUsageAnalyzer.CalculateEnc(seq);

        Assert.That(enc, Is.EqualTo(29.0).Within(1e-9),
            "Phe TTTx3/TTCx1 gives F̂=0.5 (Eq.1); Nc = 2 + 9/0.5 + 1 + 5 + 3 = 29 (Eq.3 with absent classes at full count).");
    }

    // M4 — Invariant INV-01: 20 ≤ Nc ≤ 61 for any non-empty coding sequence (property test).
    [TestCase("ATGAAAGAGCTGTTCGCCAAA")]
    [TestCase("ATGGCTGCAGCTGCAGGTGGCGGAGGG")]
    [TestCase("TTTTTCTTATTGCTTCTCCTACTG")]
    [TestCase("ATGTGGATGTGGATGTGG")]
    [TestCase("AAAAAAAAAAAAAAAAAA")]
    public void CalculateEnc_AnyValidSequence_StaysWithinRange(string seq)
    {
        double enc = CodonUsageAnalyzer.CalculateEnc(seq);

        Assert.Multiple(() =>
        {
            Assert.That(enc, Is.GreaterThanOrEqualTo(20.0),
                "INV-01: Nc cannot fall below the extreme-bias limit 20.");
            Assert.That(enc, Is.LessThanOrEqualTo(61.0),
                "INV-01: Nc cannot exceed 61 (Eq. 3 re-adjustment).");
        });
    }

    // M5 — Isoleucine absent, but a 2-fold (Phe) and a 4-fold (Ala) amino acid present.
    // Wright Eq. (5a) sets F̂₃ = (F̂₂ + F̂₄)/2. Hand derivation:
    //   Phe TTTx3,TTCx1: F̂₂ = 0.5 ⇒ 9/F̂₂ = 18.
    //   Ala GCTx2,GCCx2: n=4, Σp²=0.5, F̂₄ = (4*0.5-1)/3 = 1/3 ⇒ 5/F̂₄ = 15.
    //   F̂₃ = (0.5 + 1/3)/2 = 0.416666...; 1/F̂₃ = 2.4.
    //   6-fold absent ⇒ full count 3. Met+Trp = 2.
    //   Nc = 2 + 18 + 2.4 + 15 + 3 = 40.4.
    [Test]
    public void CalculateEnc_IsoleucineAbsent_UsesEq5aFallback()
    {
        string seq = Repeat("TTT", 3) + "TTC" + Repeat("GCT", 2) + Repeat("GCC", 2);

        double enc = CodonUsageAnalyzer.CalculateEnc(seq);

        Assert.That(enc, Is.EqualTo(40.4).Within(1e-9),
            "With no isoleucine, F̂₃ = (F̂₂+F̂₄)/2 (Eq. 5a) gives Nc = 2 + 18 + 2.4 + 15 + 3 = 40.4.");
    }

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
