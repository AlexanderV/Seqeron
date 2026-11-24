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

/// <summary>multiple_align_mafft (audit A2-19 / A2-19b, B09 F119 / F283): MAFFT 7.505, progressive strategies; delegates to SequenceAligner.RunMafft (oracle comparison in the library tests).</summary>
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
    public void MultipleAlignMafft_ThreeSequences_CarriesTheProgramsStandardError()
    {
        var r = AlignmentTools.MultipleAlignMafft(["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ"], anySymbol: false);
        Assert.Multiple(() =>
        {
            Assert.That(r.StandardError, Does.Contain("Strategy:\n FFT-NS-2 (Fast but rough)"));
            Assert.That(AlignmentTools.MultipleAlignMafft(["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ"], quiet: true).StandardError, Is.Empty);
            Assert.That(AlignmentTools.MultipleAlignMafft(["MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ"], format: "fasta", names: ["a", "b", "c"]).Formatted, Does.StartWith(">a"));
        });
    }

    [Test]
    public void MultipleAlignMafft_MaxIterate_RefinesAndRefusesThreads()
    {
        var seqs = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ", "MKVLAAGIVGLLAQ" };
        var r = AlignmentTools.MultipleAlignMafft(seqs, anySymbol: false, maxIterate: 2);
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences, Has.Length.EqualTo(4));
            Assert.That(r.StandardError, Does.Contain("FFT-NS-i (Standard)").And.Contain("Iterative refinement method (max. 2 iterations)"));
            Assert.Throws<NotSupportedException>(() => AlignmentTools.MultipleAlignMafft(seqs, maxIterate: 2, threads: 4), "the program is non-deterministic with threads in the iterative step");
            Assert.That(AlignmentTools.MultipleAlignMafft(seqs, anySymbol: false, maxIterate: 2, threads: 4, threadsIterate: 0).AlignedSequences, Is.EqualTo(r.AlignedSequences));
        });
    }

    [Test]
    public void MultipleAlignMafft_PairLibrary_G_L_E_Ins1()
    {
        // hand run of the installed mafft 7.505: mafft --globalpair --maxiterate 0 (5 proteins of 60 residues, B09 F285): the rows of the library alignment
        var seqs = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ", "MKVLAAGIVGLLAQ" };
        foreach (var (tag, text) in new[] { ("global", "G-INS-1 (Not tested.)"), ("local", "L-INS-1 (Not tested.)"), ("genaf", "E-INS-1 (Not tested.)") })
        {
            var r = AlignmentTools.MultipleAlignMafft(seqs, anySymbol: false, pairLibrary: tag);
            Assert.Multiple(() =>
            {
                Assert.That(r.AlignedSequences, Has.Length.EqualTo(4), tag);
                Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(seqs), tag);
                Assert.That(r.StandardError, Does.Contain("Strategy:\n " + text), tag);
            });
        }
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.MultipleAlignMafft(seqs, anySymbol: false, pairLibrary: "global", maxIterate: 2).StandardError, Does.Contain("G-INS-i (Suitable for sequences of similar lengths, very slow)"), "G-INS-i");
            Assert.Throws<NotSupportedException>(() => AlignmentTools.MultipleAlignMafft(seqs, pairLibrary: "global", maxIterate: 2, threads: 4), "the program is non-deterministic with threads in the iterative step");
            Assert.That(AlignmentTools.MultipleAlignMafft(seqs, pairLibrary: "global", weightI: 0.0).AlignedSequences, Has.Length.EqualTo(seqs.Length), "--weighti 0 is ported (F317)");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(seqs, pairLibrary: "x"));
        });
    }

    [Test]
    public void MultipleAlignMafft_Auto_PicksTheStrategyOfTheScript()
    {
        var seqs = new[] { "MKVLAAGIVGLLLAQ", "MKVLAGIVALLAQ", "MKILAAGIVGLLAQ", "MKVLAAGIVGLLAQ" };
        var r = AlignmentTools.MultipleAlignMafft(seqs, anySymbol: false, auto: true);
        Assert.That(r.StandardError, Does.Contain("L-INS-i (Probably most accurate, very slow)"), "fewer than 100 sequences: L-INS-i");
    }

    [Test]
    public void MultipleAlignMafft_Validation()
    {
        Assert.Multiple(() =>
        {
            Assert.That(AlignmentTools.MultipleAlignMafft(["ACGTACGT", "ACGAACGT", "ACGGACGT"]).AlignedSequences, Has.Length.EqualTo(3), "three sequences: the progressive stage (A2-19b)");
            Assert.That(AlignmentTools.MultipleAlignMafft(["ACGT", "ACGA", "ACGG"], retree: 0, treeOut: true).TreeFile, Does.EndWith(";\n"), "--retree 0: the guide tree only (B09 F313)");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(["ACGT", "ACGA", "ACGG"], partTree: "bogus"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(["ACGT", "ACGA", "ACGG"], scoringModel: "pam"));
            Assert.That(AlignmentTools.MultipleAlignMafft(["ACGT"]).AlignedSequences, Is.EqualTo(new[] { "ACGT" }), "one sequence is written back (A3-4; --anysymbol keeps the letters)");
            Assert.That(AlignmentTools.MultipleAlignMafft([]).AlignedSequences, Is.Empty, "mafft exits 0 with an empty output (A3-4)");
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(["ACGT", "ACGA"], sequenceType: "rna"));
            Assert.Throws<ArgumentException>(() => AlignmentTools.MultipleAlignMafft(["ACGT", "ACGA"], blosum: 63));
            Assert.Throws<ArgumentNullException>(() => AlignmentTools.MultipleAlignMafft(null!));
        });
    }
}
