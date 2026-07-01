// ANNOT-CODONUSAGE-001 — Relative Synonymous Codon Usage (RSCU)
// Evidence: docs/Evidence/ANNOT-CODONUSAGE-001-Evidence.md
// TestSpec: tests/TestSpecs/ANNOT-CODONUSAGE-001.md
// Source: Sharp PM, Li W-H (1986). Nucleic Acids Res 14(19):7737-7749. https://doi.org/10.1093/nar/14.19.7737
//         Formula (verbatim): RSCU = n_i * x_ij / Σ_j x_ij — https://www.lirmm.fr/~rivals/rscu/
//         Reference impl: SouradiptoC/CodonU internal_comp.py `rscu` (sense codons only, pooled counts).

namespace Seqeron.Genomics.Tests.Unit.Annotation;

using System;
using System.Collections.Generic;
using System.Linq;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Core;

[TestFixture]
public class GenomeAnnotator_GetCodonUsage_Tests
{
    // Leucine family codons (NCBI table 1, n_i = 6).
    private static readonly string[] LeuFamily = { "TTA", "TTG", "CTT", "CTC", "CTA", "CTG" };

    #region GetCodonUsage (RSCU, Standard code)

    // M1 — Leu family: CDS CTTCTTCTGTTA → CTT=2, CTG=1, TTA=1, Σ=4, n_i=6.
    //      RSCU = 6*x/4: CTT=3.0, CTG=1.5, TTA=1.5, TTG/CTC/CTA=0. (LIRMM formula, NCBI Leu family.)
    [Test]
    public void GetCodonUsage_LeucineWorkedExample_ReturnsExactRscuValues()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "CTTCTTCTGTTA" });

        Assert.Multiple(() =>
        {
            Assert.That(rscu["CTT"], Is.EqualTo(3.0).Within(1e-10), "CTT: 6*2/4 = 3.0 (preferred Leu codon).");
            Assert.That(rscu["CTG"], Is.EqualTo(1.5).Within(1e-10), "CTG: 6*1/4 = 1.5.");
            Assert.That(rscu["TTA"], Is.EqualTo(1.5).Within(1e-10), "TTA: 6*1/4 = 1.5.");
            Assert.That(rscu["TTG"], Is.EqualTo(0.0).Within(1e-10), "TTG unobserved: 6*0/4 = 0.0.");
            Assert.That(rscu["CTC"], Is.EqualTo(0.0).Within(1e-10), "CTC unobserved: 0.0.");
            Assert.That(rscu["CTA"], Is.EqualTo(0.0).Within(1e-10), "CTA unobserved: 0.0.");
        });
    }

    // M2 — Uniform usage = no bias: Phe TTT=1, TTC=1, n_i=2, Σ=2 → RSCU 1.0 each. (PMC2528880: RSCU=1.0 ⇒ no bias.)
    [Test]
    public void GetCodonUsage_UniformPhenylalanineUsage_ReturnsOnePointZero()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "TTTTTC" });

        Assert.Multiple(() =>
        {
            Assert.That(rscu["TTT"], Is.EqualTo(1.0).Within(1e-10), "Uniform usage ⇒ RSCU = 1.0 (no bias).");
            Assert.That(rscu["TTC"], Is.EqualTo(1.0).Within(1e-10), "Uniform usage ⇒ RSCU = 1.0 (no bias).");
        });
    }

    // M3 — Single-codon amino acid: Met (ATG) n_i=1 ⇒ RSCU always 1.0 regardless of count. (NCBI table 1, INV-02.)
    [Test]
    public void GetCodonUsage_SingleCodonAminoAcid_ReturnsOnePointZero()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "ATGATG" });

        Assert.That(rscu["ATG"], Is.EqualTo(1.0).Within(1e-10),
            "Methionine has one codon (n_i=1), so RSCU = 1*2/2 = 1.0.");
    }

    // M3b — Single-codon amino acid Trp (TGG), the other n_i=1 family (INV-02). NCBI table 1: Trp = TGG only.
    [Test]
    public void GetCodonUsage_TryptophanSingleCodon_ReturnsOnePointZero()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "TGGTGGTGG" });

        Assert.That(rscu["TGG"], Is.EqualTo(1.0).Within(1e-10),
            "Tryptophan has one codon (n_i=1), so RSCU = 1*3/3 = 1.0 regardless of count.");
    }

    // INV-04 / table-1 cardinality: the output enumerates exactly the 61 sense codons of the
    // Standard code (64 codons − 3 stops), and no stop codon ever appears. (NCBI translation table 1.)
    [Test]
    public void GetCodonUsage_OutputContainsExactlyThe61SenseCodons_NoStops()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "ATGTGG" });

        Assert.Multiple(() =>
        {
            Assert.That(rscu.Count, Is.EqualTo(61),
                "Standard code has 64 codons − 3 stops = 61 sense codons (NCBI table 1).");
            Assert.That(rscu.ContainsKey("TAA"), Is.False, "TAA stop excluded.");
            Assert.That(rscu.ContainsKey("TAG"), Is.False, "TAG stop excluded.");
            Assert.That(rscu.ContainsKey("TGA"), Is.False, "TGA stop excluded.");
        });
    }

    // M4 — Pooling across sequences: ["CTTCTT","CTGTTA"] pools to the M1 Leu counts (CodonU pools the reference set).
    [Test]
    public void GetCodonUsage_PoolsCountsAcrossSequences_MatchesAggregate()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "CTTCTT", "CTGTTA" });

        Assert.Multiple(() =>
        {
            Assert.That(rscu["CTT"], Is.EqualTo(3.0).Within(1e-10), "Pooled CTT=2 across both sequences ⇒ 3.0.");
            Assert.That(rscu["CTG"], Is.EqualTo(1.5).Within(1e-10), "Pooled CTG=1 ⇒ 1.5.");
            Assert.That(rscu["TTA"], Is.EqualTo(1.5).Within(1e-10), "Pooled TTA=1 ⇒ 1.5.");
        });
    }

    // M5 — Stop codons excluded: ATGTAA (Met + TAA stop). Output has ATG, no stop codon keys. (CodonU forward_table.)
    [Test]
    public void GetCodonUsage_StopCodons_AreExcludedFromOutput()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "ATGTAA" });

        Assert.Multiple(() =>
        {
            Assert.That(rscu["ATG"], Is.EqualTo(1.0).Within(1e-10), "Met still reported with RSCU 1.0.");
            Assert.That(rscu.ContainsKey("TAA"), Is.False, "TAA is a stop codon, excluded from RSCU.");
            Assert.That(rscu.ContainsKey("TAG"), Is.False, "TAG stop codon excluded.");
            Assert.That(rscu.ContainsKey("TGA"), Is.False, "TGA stop codon excluded.");
        });
    }

    // M6 — INV-01: Σ RSCU over an observed synonymous family = n_i (= 6 for Leu).
    [Test]
    public void GetCodonUsage_FamilySumInvariant_EqualsFamilySize()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "CTTCTTCTGTTA" });

        double leuSum = LeuFamily.Sum(c => rscu[c]);

        Assert.That(leuSum, Is.EqualTo(6.0).Within(1e-10),
            "Σ RSCU over the 6-codon Leu family equals n_i = 6 (INV-01).");
    }

    // S1 — Case-insensitive input: lower-case yields identical M1 result.
    [Test]
    public void GetCodonUsage_LowerCaseInput_TreatedCaseInsensitively()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "cttcttctgtta" });

        Assert.Multiple(() =>
        {
            Assert.That(rscu["CTT"], Is.EqualTo(3.0).Within(1e-10), "Lower-case input uppercased ⇒ CTT 3.0.");
            Assert.That(rscu["CTG"], Is.EqualTo(1.5).Within(1e-10), "Lower-case input uppercased ⇒ CTG 1.5.");
        });
    }

    // S2 — Partial trailing codon ignored: ATGAT (5 nt) reads ATG only; trailing AT dropped.
    [Test]
    public void GetCodonUsage_PartialTrailingCodon_IsIgnored()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { "ATGAT" });

        Assert.That(rscu["ATG"], Is.EqualTo(1.0).Within(1e-10),
            "Only the complete in-frame codon ATG is counted; trailing 'AT' is ignored.");
    }

    // C1 — Null input throws ArgumentNullException.
    [Test]
    public void GetCodonUsage_NullSequences_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => GenomeAnnotator.GetCodonUsage((IEnumerable<string>)null!),
            "Null coding-sequence collection must throw ArgumentNullException.");
    }

    // C2 — Empty / all-empty input yields all-zero RSCU (no codon observed; families total 0 ⇒ 0.0).
    [Test]
    public void GetCodonUsage_EmptyInput_ReturnsAllZeroRscu()
    {
        var fromEmptyString = GenomeAnnotator.GetCodonUsage(new[] { "" });
        var fromEmptyList = GenomeAnnotator.GetCodonUsage(Array.Empty<string>());

        Assert.Multiple(() =>
        {
            Assert.That(fromEmptyString.Count, Is.EqualTo(61), "Output enumerates all 61 sense codons of the code.");
            Assert.That(fromEmptyString.Values.All(v => v == 0.0), Is.True,
                "With nothing observed, every synonymous family totals 0 ⇒ RSCU 0.0 for all codons.");
            Assert.That(fromEmptyList.Count, Is.EqualTo(61), "Empty enumerable also yields the 61-codon zero map.");
            Assert.That(fromEmptyList.Values.All(v => v == 0.0), Is.True,
                "Empty enumerable: nothing observed ⇒ RSCU 0.0 for all codons.");
        });
    }

    #endregion

    #region GetCodonUsage (GeneticCode overload — delegation smoke)

    // S3 — Explicit GeneticCode.Standard overload equals the default overload (delegation).
    [Test]
    public void GetCodonUsage_StandardCodeOverload_MatchesDefaultOverload()
    {
        var defaultResult = GenomeAnnotator.GetCodonUsage(new[] { "CTTCTTCTGTTA" });
        var explicitResult = GenomeAnnotator.GetCodonUsage(new[] { "CTTCTTCTGTTA" }, GeneticCode.Standard);

        Assert.Multiple(() =>
        {
            Assert.That(explicitResult["CTT"], Is.EqualTo(defaultResult["CTT"]).Within(1e-10),
                "Explicit Standard-code overload delegates to the same computation.");
            Assert.That(explicitResult["CTG"], Is.EqualTo(defaultResult["CTG"]).Within(1e-10),
                "Explicit Standard-code overload delegates to the same computation.");
        });
    }

    #endregion
}
