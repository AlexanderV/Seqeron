// ONCO-HLA-001 — HLA nomenclature parsing + allele-specific HLA LOH (LOHHLA)
// Evidence: docs/Evidence/ONCO-HLA-001-Evidence.md
// TestSpec: tests/TestSpecs/ONCO-HLA-001.md
// Source: Marsh SGE et al. (2010). Tissue Antigens 75(4):291–455 (WHO HLA Nomenclature; hla.alleles.org "Naming Alleles").
//         McGranahan N et al. (2017). Cell 171(6):1259–1271 (LOHHLA): allele CN < 0.5 => loss; allelic imbalance paired-t p < 0.01.

using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Oncology;

[TestFixture]
public class OncologyAnalyzer_HlaAnalysis_Tests
{
    #region ParseHlaAllele

    // M1 — hla.alleles.org: two-field name HLA-A*02:01 (minimum valid; F1=type, F2=protein).
    [Test]
    public void ParseHlaAllele_TwoFieldName_ParsesGeneAndFields()
    {
        HlaAllele allele = OncologyAnalyzer.ParseHlaAllele("HLA-A*02:01");

        Assert.Multiple(() =>
        {
            Assert.That(allele.Gene, Is.EqualTo("A"), "Gene token between 'HLA-' and '*' is A.");
            Assert.That(allele.Fields, Is.EqualTo(new[] { "02", "01" }), "Two colon-separated digit fields, verbatim with leading zeros.");
            Assert.That(allele.AlleleGroup, Is.EqualTo("02"), "Field 1 is the allele group / type.");
            Assert.That(allele.Protein, Is.EqualTo("01"), "Field 2 is the specific HLA protein.");
            Assert.That(allele.Suffix, Is.EqualTo(HlaExpressionSuffix.None), "No trailing expression suffix.");
            Assert.That(allele.Name, Is.EqualTo("HLA-A*02:01"), "Normalized name round-trips the input (INV-02).");
        });
    }

    // M2 — three-field name HLA-B*07:02:01 (F3 = synonymous coding substitutions).
    [Test]
    public void ParseHlaAllele_ThreeFieldName_ParsesAllThreeFields()
    {
        HlaAllele allele = OncologyAnalyzer.ParseHlaAllele("HLA-B*07:02:01");

        Assert.Multiple(() =>
        {
            Assert.That(allele.Gene, Is.EqualTo("B"), "Gene is B.");
            Assert.That(allele.Fields, Is.EqualTo(new[] { "07", "02", "01" }), "Three fields per WHO nomenclature (type:protein:synonymous).");
            Assert.That(allele.Suffix, Is.EqualTo(HlaExpressionSuffix.None), "No suffix.");
        });
    }

    // M3 — four-field name HLA-C*07:02:01:03 (F4 = non-coding differences; max field count).
    [Test]
    public void ParseHlaAllele_FourFieldName_ParsesAllFourFields()
    {
        HlaAllele allele = OncologyAnalyzer.ParseHlaAllele("HLA-C*07:02:01:03");

        Assert.Multiple(() =>
        {
            Assert.That(allele.Gene, Is.EqualTo("C"), "Gene is C.");
            Assert.That(allele.Fields, Is.EqualTo(new[] { "07", "02", "01", "03" }), "Four fields (type:protein:synonymous:non-coding) is the maximum.");
            Assert.That(allele.Name, Is.EqualTo("HLA-C*07:02:01:03"), "Round-trips four fields (INV-02).");
        });
    }

    // M4 — expression suffix HLA-A*24:02:01:02L (L = Low cell-surface expression per hla.alleles.org).
    [Test]
    public void ParseHlaAllele_WithExpressionSuffix_ParsesSuffix()
    {
        HlaAllele allele = OncologyAnalyzer.ParseHlaAllele("HLA-A*24:02:01:02L");

        Assert.Multiple(() =>
        {
            Assert.That(allele.Fields, Is.EqualTo(new[] { "24", "02", "01", "02" }), "Suffix letter is stripped before field splitting.");
            Assert.That(allele.Suffix, Is.EqualTo(HlaExpressionSuffix.Low), "Trailing 'L' maps to Low cell-surface expression.");
            Assert.That(allele.Name, Is.EqualTo("HLA-A*24:02:01:02L"), "Normalized name re-appends the suffix (INV-02).");
        });
    }

    // M4b — every WHO expression suffix {N,L,S,C,A,Q} maps to its documented status and round-trips
    // in the normalized Name. Letters/meanings verbatim from hla.alleles.org "Assigning Suffixes":
    // N=Null(not expressed), L=Low cell-surface, S=Secreted (soluble, not on surface),
    // C=Cytoplasm (not on surface), A=Aberrant (doubt whether expressed), Q=Questionable.
    [TestCase("HLA-A*01:01:01:02N", HlaExpressionSuffix.Null)]
    [TestCase("HLA-A*24:02:01:02L", HlaExpressionSuffix.Low)]
    [TestCase("HLA-B*44:02:01:02S", HlaExpressionSuffix.Secreted)]
    [TestCase("HLA-A*01:01:01:01C", HlaExpressionSuffix.Cytoplasm)]
    [TestCase("HLA-A*01:01:01:01A", HlaExpressionSuffix.Aberrant)]
    [TestCase("HLA-A*32:11Q", HlaExpressionSuffix.Questionable)]
    public void ParseHlaAllele_EachExpressionSuffix_MapsAndRoundTrips(string name, HlaExpressionSuffix expected)
    {
        HlaAllele allele = OncologyAnalyzer.ParseHlaAllele(name);

        Assert.Multiple(() =>
        {
            Assert.That(allele.Suffix, Is.EqualTo(expected), "Suffix letter maps to the WHO expression status.");
            Assert.That(allele.Name, Is.EqualTo(name), "Normalized Name re-appends the suffix verbatim (INV-02).");
        });
    }

    // M4c — the suffix letter is matched case-insensitively (lowercase 'l' => Low).
    [Test]
    public void ParseHlaAllele_LowercaseSuffix_MapsCaseInsensitively()
    {
        HlaAllele allele = OncologyAnalyzer.ParseHlaAllele("HLA-A*24:02:01:02l");

        Assert.That(allele.Suffix, Is.EqualTo(HlaExpressionSuffix.Low), "Suffix is matched case-insensitively.");
    }

    // M5 — missing 'HLA-' prefix is invalid.
    [Test]
    public void ParseHlaAllele_MissingPrefix_ThrowsFormatException()
    {
        Assert.That(() => OncologyAnalyzer.ParseHlaAllele("A*02:01"), NUnit.Framework.Throws.TypeOf<FormatException>(),
            "WHO nomenclature requires the 'HLA-' prefix; 'A*02:01' has none.");
    }

    // M6 — single field HLA-A*02 is incomplete (two-field minimum: "at least a four digit name").
    [Test]
    public void ParseHlaAllele_SingleField_ThrowsFormatException()
    {
        Assert.That(() => OncologyAnalyzer.ParseHlaAllele("HLA-A*02"), NUnit.Framework.Throws.TypeOf<FormatException>(),
            "All alleles receive at least a four-digit (two-field) name; one field is invalid.");
    }

    // M7 — five fields exceeds the four-field maximum.
    [Test]
    public void ParseHlaAllele_FiveFields_ThrowsFormatException()
    {
        Assert.That(() => OncologyAnalyzer.ParseHlaAllele("HLA-A*02:01:01:01:01"), NUnit.Framework.Throws.TypeOf<FormatException>(),
            "WHO nomenclature allows up to four fields; five is invalid.");
    }

    // M8 — trailing letter outside {N,L,S,C,A,Q} is not a valid expression suffix.
    [Test]
    public void ParseHlaAllele_InvalidSuffix_ThrowsFormatException()
    {
        Assert.That(() => OncologyAnalyzer.ParseHlaAllele("HLA-A*02:01X"), NUnit.Framework.Throws.TypeOf<FormatException>(),
            "'X' is not in the allowed expression-status suffix set {N,L,S,C,A,Q}.");
    }

    // M8b — a non-numeric field is invalid (fields are digit groups).
    [Test]
    public void ParseHlaAllele_NonNumericField_ThrowsFormatException()
    {
        Assert.That(() => OncologyAnalyzer.ParseHlaAllele("HLA-A*0A:01"), NUnit.Framework.Throws.TypeOf<FormatException>(),
            "HLA fields must be digit groups; '0A' contains a letter.");
    }

    // M8c — the gene/field separator '*' is mandatory and must bound a gene and a field block:
    // missing '*', empty gene ('HLA-*...'), and empty field block ('HLA-A*') are all invalid.
    [TestCase("HLA-A02:01")]    // no '*' at all
    [TestCase("HLA-*02:01")]    // '*' present but gene token empty
    [TestCase("HLA-A*")]        // '*' present but field block empty
    public void ParseHlaAllele_MissingOrEmptyGeneFieldSeparator_ThrowsFormatException(string name)
    {
        Assert.That(() => OncologyAnalyzer.ParseHlaAllele(name), NUnit.Framework.Throws.TypeOf<FormatException>(),
            "WHO nomenclature requires a gene and a field block separated by '*'.");
    }

    // C2 — lowercase gene/prefix is normalized to upper-case.
    [Test]
    public void ParseHlaAllele_LowercaseGene_NormalizesToUpper()
    {
        HlaAllele allele = OncologyAnalyzer.ParseHlaAllele("hla-a*02:01");

        Assert.That(allele.Gene, Is.EqualTo("A"), "Gene name is upper-cased; the 'HLA-' prefix is matched case-insensitively.");
    }

    // S3 — null / empty / whitespace are rejected by argument validation.
    [Test]
    public void ParseHlaAllele_NullEmptyOrWhitespace_Throws()
    {
        Assert.Multiple(() =>
        {
            Assert.That(() => OncologyAnalyzer.ParseHlaAllele(null!), NUnit.Framework.Throws.TypeOf<ArgumentNullException>(),
                "Null name is an ArgumentNullException.");
            Assert.That(() => OncologyAnalyzer.ParseHlaAllele(""), NUnit.Framework.Throws.TypeOf<ArgumentException>(),
                "Empty name is an ArgumentException.");
            Assert.That(() => OncologyAnalyzer.ParseHlaAllele("   "), NUnit.Framework.Throws.TypeOf<ArgumentException>(),
                "Whitespace-only name is an ArgumentException.");
        });
    }

    #endregion

    #region TryParseHlaAllele (delegate — smoke)

    // S1 — TryParse returns true and the parsed allele for a valid name (delegates to ParseHlaAllele).
    [Test]
    public void TryParseHlaAllele_ValidName_ReturnsTrue()
    {
        bool ok = OncologyAnalyzer.TryParseHlaAllele("HLA-A*02:01", out HlaAllele allele);

        Assert.Multiple(() =>
        {
            Assert.That(ok, Is.True, "Valid name parses successfully.");
            Assert.That(allele.Gene, Is.EqualTo("A"), "Out value is populated from the delegate.");
        });
    }

    // S2 — TryParse returns false and default for an invalid name (no throw).
    [Test]
    public void TryParseHlaAllele_InvalidName_ReturnsFalse()
    {
        bool ok = OncologyAnalyzer.TryParseHlaAllele("bad", out HlaAllele allele);

        Assert.Multiple(() =>
        {
            Assert.That(ok, Is.False, "Invalid name returns false instead of throwing.");
            Assert.That(allele, Is.EqualTo(default(HlaAllele)), "Out value is default on failure.");
        });
    }

    // S2b — TryParse returns false for null (no throw).
    [Test]
    public void TryParseHlaAllele_Null_ReturnsFalse()
    {
        bool ok = OncologyAnalyzer.TryParseHlaAllele(null, out _);

        Assert.That(ok, Is.False, "Null name returns false rather than throwing.");
    }

    #endregion

    #region DetectHlaLoh

    // M9 — LOHHLA: allele 2 CN 0.30 < 0.5 with significant imbalance (p=0.001 < 0.01) => LOH, lose allele 2.
    [Test]
    public void DetectHlaLoh_Allele2LowWithSignificantImbalance_CallsLohOnAllele2()
    {
        HlaLohResult result = OncologyAnalyzer.DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 1.8, "HLA-A*11:01", 0.30, 0.001));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsLoh, Is.True, "CN(allele2)=0.30 < 0.5 and p=0.001 < 0.01 => HLA LOH (McGranahan 2017).");
            Assert.That(result.LostAllele, Is.EqualTo(HlaLostAllele.Allele2), "The sub-0.5-copy allele (allele 2) is the lost allele (INV-04).");
            Assert.That(result.AllelicImbalanceSignificant, Is.True, "p=0.001 < 0.01 => significant imbalance.");
        });
    }

    // M10 — symmetric: allele 1 CN 0.10 < 0.5, p=0.0005 < 0.01 => LOH, lose allele 1.
    [Test]
    public void DetectHlaLoh_Allele1LowWithSignificantImbalance_CallsLohOnAllele1()
    {
        HlaLohResult result = OncologyAnalyzer.DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 0.10, "HLA-A*11:01", 1.50, 0.0005));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsLoh, Is.True, "CN(allele1)=0.10 < 0.5 and p=0.0005 < 0.01 => HLA LOH.");
            Assert.That(result.LostAllele, Is.EqualTo(HlaLostAllele.Allele1), "Allele 1 is the lost allele (INV-04).");
        });
    }

    // M11 — both alleles CN >= 0.5 => heterozygous retained, no LOH.
    [Test]
    public void DetectHlaLoh_BothAllelesRetained_NoLoh()
    {
        HlaLohResult result = OncologyAnalyzer.DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 1.10, "HLA-A*11:01", 0.90, 0.30));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsLoh, Is.False, "Both CN >= 0.5 => no allele lost => no LOH.");
            Assert.That(result.LostAllele, Is.EqualTo(HlaLostAllele.None), "No lost allele.");
        });
    }

    // M12 — over-calling guard: CN 0.40 < 0.5 but p=0.05 >= 0.01 => NOT LOH.
    [Test]
    public void DetectHlaLoh_LowCopyButNonSignificantImbalance_NoLoh()
    {
        HlaLohResult result = OncologyAnalyzer.DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 1.60, "HLA-A*11:01", 0.40, 0.05));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsLoh, Is.False, "p=0.05 >= 0.01 => imbalance not significant => no LOH (explicit LOHHLA over-calling guard).");
            Assert.That(result.AllelicImbalanceSignificant, Is.False, "p=0.05 is not below the 0.01 threshold.");
            Assert.That(result.LostAllele, Is.EqualTo(HlaLostAllele.None), "No LOH call despite low copy number.");
        });
    }

    // M13 — boundary: CN exactly 0.5 is retained (strict '< 0.5').
    [Test]
    public void DetectHlaLoh_CopyNumberExactlyHalf_NoLoh()
    {
        HlaLohResult result = OncologyAnalyzer.DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 1.50, "HLA-A*11:01", 0.50, 0.001));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsLoh, Is.False, "CN=0.5 is NOT < 0.5; the threshold is strict (INV-05).");
            Assert.That(result.LostAllele, Is.EqualTo(HlaLostAllele.None), "No allele is below 0.5.");
        });
    }

    // M14 — boundary: p exactly 0.01 is not significant (strict 'p < 0.01').
    [Test]
    public void DetectHlaLoh_PValueExactlyThreshold_NoLoh()
    {
        HlaLohResult result = OncologyAnalyzer.DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 1.70, "HLA-A*11:01", 0.40, 0.01));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsLoh, Is.False, "p=0.01 is NOT < 0.01; the threshold is strict (INV-05).");
            Assert.That(result.AllelicImbalanceSignificant, Is.False, "p=0.01 is not below 0.01.");
        });
    }

    // C1 — both alleles < 0.5 with significant imbalance => homozygous loss, not allele-specific LOH (ASM-01).
    [Test]
    public void DetectHlaLoh_BothAllelesLost_ReportsHomozygousLossNotLoh()
    {
        HlaLohResult result = OncologyAnalyzer.DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 0.20, "HLA-A*11:01", 0.30, 0.001));

        Assert.Multiple(() =>
        {
            Assert.That(result.IsLoh, Is.False, "Both CN < 0.5 is homozygous deletion, not allele-specific LOH (Evidence assumption ASM-01).");
            Assert.That(result.LostAllele, Is.EqualTo(HlaLostAllele.Both), "Both homologs below threshold => Both.");
        });
    }

    // S4 — negative copy number is rejected.
    [Test]
    public void DetectHlaLoh_NegativeCopyNumber_ThrowsArgumentException()
    {
        Assert.That(() => OncologyAnalyzer.DetectHlaLoh(
                new HlaAlleleCopyNumber("HLA-A*02:01", -1.0, "HLA-A*11:01", 1.0, 0.001)),
            NUnit.Framework.Throws.InstanceOf<ArgumentException>(), "Copy numbers must be non-negative.");
    }

    // S5 — p value outside [0,1] is rejected.
    [Test]
    public void DetectHlaLoh_PValueOutOfRange_ThrowsArgumentException()
    {
        Assert.That(() => OncologyAnalyzer.DetectHlaLoh(
                new HlaAlleleCopyNumber("HLA-A*02:01", 1.0, "HLA-A*11:01", 0.30, 1.5)),
            NUnit.Framework.Throws.InstanceOf<ArgumentException>(), "Allelic-imbalance p value must be in [0, 1].");
    }

    #endregion
}
