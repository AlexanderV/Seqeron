using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A3-9d (B09 F205): <c>promer</c> of MUMmer 3.23 - the perl script, <c>prepro</c>, <c>mummer</c> on the translations, <c>mgaps</c>, <c>postpro</c>, <c>show-coords -r</c> - on random protein-coding DNA
/// (TestData/Alignment/promer_oracle.jsonl.gz, scripts/alignment/generate_promer_oracle.py; every option family, the temporary files compared as well) and the delta programs on PROMER delta files
/// (TestData/Alignment/promer_show_oracle.jsonl.gz, scripts/alignment/generate_promer_show_oracle.py: <c>show-aligns</c> and <c>show-snps</c> in amino acids, <c>show-coords -k</c>, <c>show-tiling</c>, <c>show-diff</c>, <c>delta-filter</c>).
/// </summary>
[TestFixture]
public class MummerPromer_Oracle_Tests
{
    [Test]
    public void Promer_AllStagesEqualTheProgram_AllRows()
    {
        var (failures, c) = MummerPromerOracleCheck.CheckAll();
        Assert.That(failures, Is.Empty, $"{failures.Count} checks differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(c["rows"], Is.GreaterThanOrEqualTo(300));
        Assert.That(c["errors"], Is.GreaterThan(3));
        Assert.That(c["withAlignments"], Is.GreaterThan(250));
        Assert.That(c["reverse"], Is.GreaterThan(100));
        Assert.That(c["multi"], Is.GreaterThan(100));
        Assert.That(c["nodelta"], Is.GreaterThan(10));
        Assert.That(c["coords"], Is.GreaterThan(30));
        Assert.That(c["matrices"], Is.GreaterThan(50));
        Assert.That(c["warnings"], Is.GreaterThan(30));
        Assert.That(c["modes"], Is.EqualTo(3));
    }

    [Test]
    public void DeltaPrograms_OnPromerDeltaFiles_EqualTheProgram_AllRuns()
    {
        var (failures, cases, perProgram, f) = MummerShowOracleCheck.CheckAll("promer_show_oracle.jsonl.gz");
        Assert.That(failures, Is.Empty, $"{failures.Count} runs differ:\n" + string.Join("\n", failures.Take(6)));
        Assert.That(cases, Is.GreaterThanOrEqualTo(100));
        foreach (string program in new[] { "delta-filter", "show-coords", "show-aligns", "show-snps", "show-diff", "show-tiling" })
            Assert.That(perProgram.GetValueOrDefault(program), Is.GreaterThanOrEqualTo(200), program);
        Assert.That(f["tilingWithPath"], Is.GreaterThan(100));
        Assert.That(f["snpSelect"], Is.GreaterThan(5));
        Assert.That(f["filterLis"], Is.GreaterThan(100));
    }

    private static byte[] OneBased(string s) => [0, .. Encoding.ASCII.GetBytes(s), 0];

    [Test]
    public void Translate_SixFrames_AsTheProgram()
    {
        // prepro -m 0 -q on ">g2 atgaaatagttyrcnacgNNNtgg": the amino acids before the masking (the program writes the stop as J and the unknown X as O: MKJFOTOWO, JNSOOROO, EIVOOOOO)
        var warnings = new StringBuilder();
        var dna = OneBased("atgaaatagttyrcnacgnnntgg");
        string Frame(int f)
        {
            var (t, n) = MummerTranslate.Translate(dna, dna.Length - 2, f, warnings);
            return Encoding.ASCII.GetString(t, 1, n);
        }

        Assert.That(Frame(1), Is.EqualTo("MK*FXTXW"));
        Assert.That(Frame(2), Is.EqualTo("*NSXXRX"));
        Assert.That(Frame(3), Is.EqualTo("EIVXXXX"));
        Assert.That(MummerTranslate.Translate(dna, dna.Length - 2, 7, warnings).Length, Is.EqualTo(-1));
        Assert.That(warnings.ToString(), Is.Empty);
        // a letter that is no nucleotide: a warning in the forward frames, the word before compdna (0 = A) in the reverse frames
        var bad = OneBased("cccecccc");
        var w = new StringBuilder();
        MummerTranslate.Translate(bad, 8, 1, w);
        Assert.That(w.ToString(), Is.EqualTo("WARNING: Forcing unrecognized DNA char to N\n"));
        w.Clear();
        var (rev, rn) = MummerTranslate.Translate(bad, 8, 4, w);
        var (revT, revTn) = MummerTranslate.Translate(OneBased("ccctcccc"), 8, 4, w);
        Assert.That(w.ToString(), Is.Empty);
        Assert.That(Encoding.ASCII.GetString(rev, 1, rn), Is.EqualTo(Encoding.ASCII.GetString(revT, 1, revTn)));   // the character acts as t (its complement is read as A)
    }

    [Test]
    public void Prepro_TranslatesAndMasksLikeTheProgram()
    {
        // prepro -m 5 -q / -r g.fa, ">g1 desc GCTTAAGCTGCTGCTGCTGCTTAAGCT" ("A*AAAAA*A": the stretches between stops and the sequence ends of at most 5 residues are masked)
        const string fasta = ">g1 desc\nGCTTAAGCTGCTGCTGCTGCTTAAGCT\n>g2\natgaaatagttyrcnacgNNNtgg\n";
        string q = PromerPrepro.Run(fasta, false, 5, "g.fa", out string warnings);
        Assert.That(q, Does.StartWith(">g1.1\nOJOOOOOJOO\n>g1.2\nLKLLLLLKO\n>g1.3\nLSCCCCLSO\n>g1.4\nSLSSSSSLSO\n>g1.5\nOJOOOOOJO\n>g1.6\nLKQQQQLKO\n>g2.1\nOOJOOOOOO\n"));
        Assert.That(q, Does.EndWith(">g2.6\nOOOOTISO\n"));
        Assert.That(warnings, Is.Empty);
        Assert.That(PromerPrepro.Run(fasta, true, 5, "g.fa", out _),
            Is.EqualTo(">allcontigs g.fa\nXJXXXXXJXXLKLLLLLKXLSCCCCLSXSLSSSSSLSXXJXXXXXJXLKQQQQLKXXXJX\nXXXXXJNSXXRXXEIVXXXXXPXRXXLFHXXXVXNYFXXXXXTISX\n"));
        Assert.That(PromerPrepro.Run(fasta, true, 0, "g.fa", out _), Does.StartWith(">allcontigs g.fa\nAJAAAAAJAXLKLLLLLKXLSCCCCLSXSLSSSSSLSXAJAAAAAJXLKQQQQLKXMKJF\n"));
        PromerPrepro.Run(fasta, false, -3, "g.fa", out string negative);
        Assert.That(negative, Is.EqualTo("WARNING: Invalid maximum mask length -3, ignored\n"));
    }

    [Test]
    public void Postpro_ReadsTheMgapsTextOfThePrograms_AndPrintsTheClusterFile()
    {
        // the pipeline of the example in MummerProgramToolsTests with --nodelta: the cluster file of postpro -d (frames as signed numbers, three nucleotides per amino acid)
        const string reference = ">ref gene\nACGTTGCATATAAAAACACTGCACGGATTTGCGGAGAGAATTCGATTCCCGTCGAATGGATCGCACACCAGAAGGGTCCAGGAGATTTCTATTATCGATTCTTTGCAGCACATCACGGAGGATTGATGCATGTACACGGCCATCCTCTAGCGAATCAAGGATTCAAGGACATGTTATATTGGGCATATATGGATCTGAACTTATCCATGGGATCCTCTGGATC\n";
        const string query = ">qry\nTTTACAAGAATACGGCCAGAATCTGCGGAGAGTGGTCAATCCCGGTCGAATGGATAGCCCACCAGAAAGGTCCGGGAGACTTTTACTATCGATTCTTCGCGGGACACCACGGCGGTCTCATGCACGTACACGGACACCCACTTGCAAACCAGAAAAAGGGATTCAAAGATATGCTTTACTGGGCTTACATGGACCTTAACCTGTCTATGGGCAGTTCGACGTA\n";
        var r = PromerProgram.Run(reference, query, "ref.fa", "qry.fa", new PromerOptions { MinMatch = 5, MinCluster = 10, Delta = false });
        Assert.That(r.Delta, Is.Null);
        Assert.That(r.Cluster, Does.StartWith("ref.fa qry.fa\nPROMER\n>ref qry 223 223\n"));
        // the same mgaps text through postpro directly gives the same cluster file
        var again = PromerPostpro.Run("ref.fa", "qry.fa", reference, query, r.Mgaps, new PromerPostproOptions { Delta = false });
        Assert.That(again.Cluster, Is.EqualTo(r.Cluster));
        Assert.Throws<FormatException>(() => PromerPostpro.Run("r", "q", ">r\nacgt\n", ">q\nacgt\n", "> q.1\n       1        1      3    none      -      -\n   garbage\n"));
        Assert.Throws<FormatException>(() => PromerPostpro.Run("r", "q", "", "", ""));
    }
}
