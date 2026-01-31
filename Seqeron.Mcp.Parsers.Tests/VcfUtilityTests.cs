using NUnit.Framework;
using Seqeron.Mcp.Parsers.Tools;

namespace Seqeron.Mcp.Parsers.Tests;

[TestFixture]
public class VcfUtilityTests
{
    // ========================
    // vcf_classify Tests
    // ========================

    [Test]
    public void VcfClassify_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfClassify("A", "G"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfClassify("", "G"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfClassify("A", ""));
    }

    [Test]
    public void VcfClassify_Binding_IdentifiesSNP()
    {
        var result = ParsersTools.VcfClassify("A", "G");
        Assert.That(result.VariantType, Is.EqualTo("SNP"));
        Assert.That(result.RefLength, Is.EqualTo(1));
        Assert.That(result.AltLength, Is.EqualTo(1));
        Assert.That(result.LengthDifference, Is.EqualTo(0));
    }

    [Test]
    public void VcfClassify_Binding_IdentifiesInsertion()
    {
        var result = ParsersTools.VcfClassify("A", "ATG");
        Assert.That(result.VariantType, Is.EqualTo("Insertion"));
        Assert.That(result.LengthDifference, Is.EqualTo(2));
    }

    [Test]
    public void VcfClassify_Binding_IdentifiesDeletion()
    {
        var result = ParsersTools.VcfClassify("ATG", "A");
        Assert.That(result.VariantType, Is.EqualTo("Deletion"));
        Assert.That(result.LengthDifference, Is.EqualTo(2));
    }

    [Test]
    public void VcfClassify_Binding_IdentifiesMNP()
    {
        var result = ParsersTools.VcfClassify("AT", "GC");
        Assert.That(result.VariantType, Is.EqualTo("MNP"));
    }

    // ========================
    // vcf_is_snp Tests
    // ========================

    [Test]
    public void VcfIsSNP_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfIsSNP("A", "G"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfIsSNP("", "G"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfIsSNP("A", ""));
    }

    [Test]
    public void VcfIsSNP_Binding_DetectsSNP()
    {
        var result = ParsersTools.VcfIsSNP("A", "G");
        Assert.That(result.IsSNP, Is.True);
        Assert.That(result.RefAllele, Is.EqualTo("A"));
        Assert.That(result.AltAllele, Is.EqualTo("G"));
    }

    [Test]
    public void VcfIsSNP_Binding_NotSNPForIndel()
    {
        var result = ParsersTools.VcfIsSNP("A", "ATG");
        Assert.That(result.IsSNP, Is.False);
    }

    // ========================
    // vcf_is_indel Tests
    // ========================

    [Test]
    public void VcfIsIndel_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfIsIndel("A", "ATG"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfIsIndel("", "ATG"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfIsIndel("A", ""));
    }

    [Test]
    public void VcfIsIndel_Binding_DetectsInsertion()
    {
        var result = ParsersTools.VcfIsIndel("A", "ATG");
        Assert.That(result.IsIndel, Is.True);
        Assert.That(result.IsInsertion, Is.True);
        Assert.That(result.IsDeletion, Is.False);
    }

    [Test]
    public void VcfIsIndel_Binding_DetectsDeletion()
    {
        var result = ParsersTools.VcfIsIndel("ATG", "A");
        Assert.That(result.IsIndel, Is.True);
        Assert.That(result.IsInsertion, Is.False);
        Assert.That(result.IsDeletion, Is.True);
    }

    [Test]
    public void VcfIsIndel_Binding_NotIndelForSNP()
    {
        var result = ParsersTools.VcfIsIndel("A", "G");
        Assert.That(result.IsIndel, Is.False);
    }

    // ========================
    // vcf_variant_length Tests
    // ========================

    [Test]
    public void VcfVariantLength_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfVariantLength("A", "ATG"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfVariantLength("", "ATG"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfVariantLength("A", ""));
    }

    [Test]
    public void VcfVariantLength_Binding_CalculatesCorrectly()
    {
        var result = ParsersTools.VcfVariantLength("A", "ATGC");
        Assert.That(result.Length, Is.EqualTo(3));
        Assert.That(result.RefLength, Is.EqualTo(1));
        Assert.That(result.AltLength, Is.EqualTo(4));
    }

    [Test]
    public void VcfVariantLength_Binding_ZeroForSNP()
    {
        var result = ParsersTools.VcfVariantLength("A", "G");
        Assert.That(result.Length, Is.EqualTo(0));
    }

    // ========================
    // vcf_is_hom_ref Tests
    // ========================

    [Test]
    public void VcfIsHomRef_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfIsHomRef("0/0"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfIsHomRef(""));
    }

    [Test]
    public void VcfIsHomRef_Binding_DetectsHomRef()
    {
        var result1 = ParsersTools.VcfIsHomRef("0/0");
        Assert.That(result1.Result, Is.True);
        Assert.That(result1.CheckType, Is.EqualTo("HomozygousReference"));

        var result2 = ParsersTools.VcfIsHomRef("0|0");
        Assert.That(result2.Result, Is.True);
    }

    [Test]
    public void VcfIsHomRef_Binding_NotHomRefForHet()
    {
        var result = ParsersTools.VcfIsHomRef("0/1");
        Assert.That(result.Result, Is.False);
    }

    // ========================
    // vcf_is_hom_alt Tests
    // ========================

    [Test]
    public void VcfIsHomAlt_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfIsHomAlt("1/1"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfIsHomAlt(""));
    }

    [Test]
    public void VcfIsHomAlt_Binding_DetectsHomAlt()
    {
        var result1 = ParsersTools.VcfIsHomAlt("1/1");
        Assert.That(result1.Result, Is.True);
        Assert.That(result1.CheckType, Is.EqualTo("HomozygousAlternate"));

        var result2 = ParsersTools.VcfIsHomAlt("2|2");
        Assert.That(result2.Result, Is.True);
    }

    [Test]
    public void VcfIsHomAlt_Binding_NotHomAltForHet()
    {
        var result = ParsersTools.VcfIsHomAlt("0/1");
        Assert.That(result.Result, Is.False);
    }

    // ========================
    // vcf_is_het Tests
    // ========================

    [Test]
    public void VcfIsHet_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => ParsersTools.VcfIsHet("0/1"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfIsHet(""));
    }

    [Test]
    public void VcfIsHet_Binding_DetectsHet()
    {
        var result1 = ParsersTools.VcfIsHet("0/1");
        Assert.That(result1.Result, Is.True);
        Assert.That(result1.CheckType, Is.EqualTo("Heterozygous"));

        var result2 = ParsersTools.VcfIsHet("1|2");
        Assert.That(result2.Result, Is.True);
    }

    [Test]
    public void VcfIsHet_Binding_NotHetForHom()
    {
        var result1 = ParsersTools.VcfIsHet("0/0");
        Assert.That(result1.Result, Is.False);

        var result2 = ParsersTools.VcfIsHet("1/1");
        Assert.That(result2.Result, Is.False);
    }

    // ========================
    // vcf_has_flag Tests
    // ========================

    [Test]
    public void VcfHasFlag_Schema_ValidatesCorrectly()
    {
        var vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t30\tPASS\tDB";
        Assert.DoesNotThrow(() => ParsersTools.VcfHasFlag(vcf, "DB"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfHasFlag("", "DB"));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfHasFlag(vcf, ""));
    }

    [Test]
    public void VcfHasFlag_Binding_FindsFlagPresent()
    {
        var vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t30\tPASS\tDB\nchr1\t200\t.\tC\tT\t40\tPASS\t.";
        var result = ParsersTools.VcfHasFlag(vcf, "DB");

        Assert.That(result.Flag, Is.EqualTo("DB"));
        Assert.That(result.RecordsWithFlag, Is.EqualTo(1));
        Assert.That(result.TotalRecords, Is.EqualTo(2));
        Assert.That(result.Percentage, Is.EqualTo(50.0));
    }

    [Test]
    public void VcfHasFlag_Binding_NoFlagFound()
    {
        var vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t30\tPASS\t.";
        var result = ParsersTools.VcfHasFlag(vcf, "DB");

        Assert.That(result.RecordsWithFlag, Is.EqualTo(0));
    }

    // ========================
    // vcf_write Tests
    // ========================

    [Test]
    public void VcfWrite_Schema_ValidatesCorrectly()
    {
        var vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t30\tPASS\t.";
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfWrite("", vcf));
        Assert.Throws<ArgumentException>(() => ParsersTools.VcfWrite("test.vcf", ""));
    }

    [Test]
    public void VcfWrite_Binding_WritesToFile()
    {
        var tempFile = Path.GetTempFileName();
        try
        {
            var vcf = "##fileformat=VCFv4.3\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tA\tG\t30\tPASS\t.\nchr1\t200\t.\tC\tT\t40\tPASS\t.";
            var result = ParsersTools.VcfWrite(tempFile, vcf);

            Assert.That(result.FilePath, Is.EqualTo(tempFile));
            Assert.That(result.RecordsWritten, Is.EqualTo(2));
            Assert.That(File.Exists(tempFile), Is.True);

            var content = File.ReadAllText(tempFile);
            Assert.That(content, Does.Contain("#CHROM"));
            Assert.That(content, Does.Contain("chr1"));
        }
        finally
        {
            if (File.Exists(tempFile))
                File.Delete(tempFile);
        }
    }
}
