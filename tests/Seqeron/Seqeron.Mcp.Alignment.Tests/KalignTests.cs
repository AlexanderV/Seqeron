using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;

namespace Seqeron.Mcp.Alignment.Tests;

/// <summary>multiple_align_kalign / kalign_command: Kalign 3.3.5 (audit A2-20, B09 F118). The tools delegate to SequenceAligner.KalignRun / KalignProgram.Run; the oracle comparison lives in the library tests.</summary>
[TestFixture]
public class KalignTests
{
    private static readonly string[] Dna = ["ACGTACGTTAGCTAGCTAGCATCGATCGATGCATGCATCG", "ACGTACGTAGCTAGCTAGCATCGATCGATGCATGCCATCG", "ACGTTCGTTAGCTAGCTTGCATCGATCGATGCATGCATCG"];

    [Test]
    public void MultipleAlignKalign_EqualsTheLibrary_AndReportsKeptSequencesAndAlphabet()
    {
        var tool = AlignmentTools.MultipleAlignKalign(["ACGTACGTAC", "----", "ACGTCGTAC", "ACGTACGAC"], type: "dna");
        var lib = Seqeron.Genomics.Alignment.SequenceAligner.KalignRun(["ACGTACGTAC", "----", "ACGTCGTAC", "ACGTACGAC"], new Seqeron.Genomics.Alignment.KalignOptions { Type = Seqeron.Genomics.Alignment.KalignType.Dna });
        Assert.Multiple(() =>
        {
            Assert.That(tool.AlignedSequences, Is.EqualTo(lib.Alignment.AlignedSequences));
            Assert.That(tool.KeptIndices, Is.EqualTo(new[] { 0, 2, 3 }));
            Assert.That(tool.Alphabet, Is.EqualTo("dna"));
            Assert.That(tool.GuideTree, Has.Length.EqualTo(2));
        });
    }

    [Test]
    public void MultipleAlignKalign_OneThreadIsTheSerialVariant_ThreeSetsOfThePrograms()
    {
        var one = AlignmentTools.MultipleAlignKalign(Dna, threads: 1);
        Assert.That(one.AlignedSequences[0], Is.EqualTo("ACGTACGTTAGCTAGCTAGCATCGATCGATGCATG-CATCG"), "kalign -n 1");
        Assert.That(AlignmentTools.MultipleAlignKalign(Dna, type: "dna").AlignedSequences[0], Is.EqualTo("ACGTACGTTAGCTAGCTAGCATCGATCGATGCATGC-ATCG"), "kalign -n 4 --type dna");
    }

    [Test]
    public void MultipleAlignKalign_Validation()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignKalign(Dna, type: "nucleotide"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignKalign(Dna, type: "protein"), "DNA with --type protein is refused by the program");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignKalign(Dna, threads: 0));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignKalign(["ACGT"]));
            Assert.Throws<ArgumentNullException>(() => AlignmentTools.MultipleAlignKalign(null!));
            Assert.DoesNotThrow(() => AlignmentTools.MultipleAlignKalign(Dna, type: "dna-internal"));
        });
    }

    [Test]
    public void KalignCommand_RunsTheCommandLine_FastaOutputWith60Columns()
    {
        var files = new Dictionary<string, string> { ["in.fa"] = ">s1 first\n" + Dna[0] + "\n>s2\n" + Dna[1] + "\n>s3\n" + Dna[2] + "\n" };
        var res = AlignmentTools.KalignCommand(["-q", "-i", "in.fa", "-o", "out.fa", "--type", "dna"], files);
        Assert.Multiple(() =>
        {
            Assert.That(res.ExitCode, Is.Zero);
            Assert.That(res.Files["out.fa"], Does.StartWith(">s1 first\nACGTACGTTAGCTAGCTAGCATCGATCGATGCATGC-ATCG\n"));
            Assert.That(AlignmentTools.KalignCommand(["-q", "--version"]).Stdout, Is.EqualTo("kalign 3.3.5\n"));
            Assert.That(AlignmentTools.KalignCommand(["-q", "--gp", "3", "-i", "in.fa"], files).ExitCode, Is.EqualTo(1), "--gp is ambiguous");
            Assert.That(AlignmentTools.KalignCommand(["-q", "-i", "missing.fa"], files).Messages, Has.Some.Contain("does not exist"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.KalignCommand(null!));
        });
    }
}

/// <summary>multiple_align_mafft (audit A2-19, B09 F119): the alignment of two sequences of MAFFT 7.505; delegates to SequenceAligner.MultipleAlignMafft (oracle comparison in the library tests).</summary>
[TestFixture]
public class MafftTests
{
    [Test]
    public void MultipleAlignMafft_HandRunOfTheProgram_ProteinPair()
    {
        // mafft --anysymbol: ">a MKVLAAGIVG-LLLAQ" ">b mkvlagivallaq" gives (installed 7.505) MKVLAAgIVGLLLAQ / mkvl-agiva-llaq; here the original letters are restored
        var r = AlignmentTools.MultipleAlignMafft(["MKVLAAgIVGLLLAQ", "mkvlagivallaq"]);
        Assert.That(r.AlignedSequences, Is.EqualTo(new[] { "MKVLAAgIVGLLLAQ", "mkvl-agiva-llaq" }));
        var plain = AlignmentTools.MultipleAlignMafft(["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ"], anySymbol: false);
        Assert.That(plain.AlignedSequences, Is.EqualTo(new[] { "MKVLAAGIVGLLLAQ", "MKVL-AGIVA-LLAQ" }), "the plain command writes amino acids in upper case");
    }

    [Test]
    public void MultipleAlignMafft_Validation()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<NotSupportedException>(() => AlignmentTools.MultipleAlignMafft(["ACGT", "ACGA", "ACGG"]), "three sequences: the progressive stage is registered (A2-19b)");
            Assert.That(AlignmentTools.MultipleAlignMafft(["ACGT"]).AlignedSequences, Is.EqualTo(new[] { "ACGT" }), "one sequence is written back (A3-4; --anysymbol keeps the letters)");
            Assert.That(AlignmentTools.MultipleAlignMafft([]).AlignedSequences, Is.Empty, "mafft exits 0 with an empty output (A3-4)");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(["ACGT", "ACGA"], sequenceType: "rna"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(["ACGT", "ACGA"], blosum: 63));
            Assert.Throws<ArgumentNullException>(() => AlignmentTools.MultipleAlignMafft(null!));
        });
    }
}
