using NUnit.Framework;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-19, B09 F119): MAFFT 7.505 stage 1 - the parts that the pair oracle (<see cref="SequenceAligner_MafftPair_Tests"/>) does not name: the scoring matrices and integer
/// penalties of <c>constants()</c> against what the installed <c>disttbfast -d</c> prints, the hand-run outputs of the program, the contracts.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_Mafft_Tests
{
    private const string Amino = "ARNDCQEGHILKMFPSTWYVBZX.-J";

    private static int Index(string alphabet, char c) => alphabet.LastIndexOf(c);

    [Test]
    public void AminoMatrix_DefaultBlosum62_EqualsTheProgramsPrintedMatrix()
    {
        // disttbfast -b 62 -f -1.53 -h 0 -d (rows A, R, C, W; columns in the order A R N D C Q E G H I L K M F P S T W Y V)
        var expected = new Dictionary<char, int[]>
        {
            ['A'] = [484, -40, -52, -74, 58, 19, 14, 114, -61, -31, -45, 26, 7, -119, 18, 208, 94, -150, -75, 80],
            ['R'] = [-40, 635, 55, -59, -234, 195, 87, -128, 74, -195, -113, 305, -36, -175, -109, 23, -12, -165, -68, -147],
            ['C'] = [58, -234, -163, -241, 940, -186, -256, -147, -195, -22, -27, -200, -41, -135, -176, 12, 13, -128, -138, 19],
            ['W'] = [-150, -165, -264, -315, -128, -93, -180, -146, -131, -155, -62, -192, -41, 188, -260, -172, -140, 1128, 310, -180],
        };
        var (m, pen, pex) = SequenceAligner.MafftEngine.ScoringMatrix(false, new MafftOptions());
        Assert.Multiple(() =>
        {
            Assert.That(pen, Is.EqualTo(-917));
            Assert.That(pex, Is.Zero);
            foreach (var (letter, row) in expected)
            {
                int i = Index(Amino, letter);
                for (int j = 0; j < 20; j++)
                    Assert.That(m[i, j], Is.EqualTo(row[j]), $"{letter} x {Amino[j]}");
            }
            for (int i = 0; i < 26; i++)
                for (int j = 0; j < 26; j++)
                {
                    Assert.That(m[i, j], Is.EqualTo(m[j, i]), $"symmetry {Amino[i]} {Amino[j]}");
                    if (i >= 20 || j >= 20)
                        Assert.That(m[i, j], Is.Zero, "B Z X . - J score 0");
                }
        });
    }

    [Test]
    public void AminoMatrix_Blosum80_WithOpAndEp_EqualsTheProgramsPrintedMatrix()
    {
        // disttbfast -b 80 -f -2.0 -h -0.1 -d: "parameters: -1199, 0, -59"; the script's --ep 0.1 is -h -0.1
        var (m, pen, pex) = SequenceAligner.MafftEngine.ScoringMatrix(false, new MafftOptions { Blosum = 80, GapOpen = 2.0, Offset = 0.1 });
        int[] rowA = [571, 9, 9, 9, 121, 65, 65, 177, 9, 9, 9, 121, 65, -47, 121, 290, 177, -104, -47];
        int[] rowN = [9, 121, 683, 290, -104, 177, 121, 121, 233, -160, -160, 177, -47, -160, -47, 233, 177, -216, -47];
        Assert.Multiple(() =>
        {
            Assert.That(pen, Is.EqualTo(-1199));
            Assert.That(pex, Is.Zero);
            for (int j = 0; j < 19; j++)
                Assert.That(m[0, j], Is.EqualTo(rowA[j]), "A x " + Amino[j]);
            for (int j = 0; j < 19; j++)
                Assert.That(m[2, j], Is.EqualTo(rowN[j]), "N x " + Amino[j]);
        });
    }

    [Test]
    public void NucleotideMatrix_Kimura200Pam_EqualsTheProgramsPrintedMatrix()
    {
        const string nuc = "agctuAGCTUnNbdhkmnrsvwyx-O";
        var (m, pen, pex) = SequenceAligner.MafftEngine.ScoringMatrix(true, new MafftOptions());
        var (m2, pen2, _) = SequenceAligner.MafftEngine.ScoringMatrix(true, new MafftOptions { GapOpen = 2.0, Offset = 0.1 });
        Assert.Multiple(() =>
        {
            // disttbfast -D -d: a-a 600, a-g 128 (transition), a-c / a-t -364
            Assert.That(new[] { m[0, 0], m[0, 1], m[0, 2], m[0, 3], m[2, 3], m[3, 3] }, Is.EqualTo(new[] { 600, 128, -364, -364, 128, 600 }));
            Assert.That(m[Index(nuc, 'u'), Index(nuc, 'u')], Is.EqualTo(600), "u scores like t");
            Assert.That(m[Index(nuc, 'u'), Index(nuc, 't')], Is.EqualTo(600).Or.EqualTo(0), "the table copies row and column of t");
            Assert.That(m[Index(nuc, 'r'), Index(nuc, 'a')], Is.EqualTo((int)Math.Round((600 + 128) / 2.0, MidpointRounding.AwayFromZero)), "an ambiguity code scores the mean of its bases");
            Assert.That(pen, Is.EqualTo(-2753));
            Assert.That(pex, Is.Zero);
            // -f -2.0 -h -0.1: a-a 659, a-g 187, a-c -305
            Assert.That(new[] { m2[0, 0], m2[0, 1], m2[0, 2] }, Is.EqualTo(new[] { 659, 187, -305 }));
            Assert.That(pen2, Is.EqualTo(-3599));
        });
    }

    [Test]
    public void HandRuns_OfTheInstalledProgram()
    {
        Assert.Multiple(() =>
        {
            // mafft --quiet --thread 1 d.fa: ">a ACGUACGTTAGC" ">b ACGTACGATAGC" (nucleotides come back in lower case)
            Assert.That(SequenceAligner.MultipleAlignMafft(["ACGUACGTTAGC", "ACGTACGATAGC"], new MafftOptions { AnySymbol = false }).AlignedSequences, Is.EqualTo(new[] { "acguacgttagc", "acgtacgatagc" }));
            // IUPAC codes: ACGTACGTTAGCNNRY / ACGTACGATAGCCC (the shorter sequence ends with a free end gap)
            Assert.That(SequenceAligner.MultipleAlignMafft(["ACGTACGTTAGCNNRY", "ACGTACGATAGCCC"], new MafftOptions { AnySymbol = false }).AlignedSequences, Is.EqualTo(new[] { "acgtacgttagcnnry", "acgtacgatagccc--" }));
            // --anysymbol keeps '.', '*' and the original case: "MKVLAA.GIVG-LLL*AQ" / "mkvlagivallaq"
            Assert.That(SequenceAligner.MultipleAlignMafft(["MKVLAA.GIVG-LLL*AQ", "mkvlagivallaq"]).AlignedSequences, Is.EqualTo(new[] { "MKVLAA.GIVGLLL*AQ", "mkvla--givall--aq" }));
            // --nuc on a protein pair
            Assert.That(SequenceAligner.MultipleAlignMafft(["ACGTACGTTAGCNNRY", "ACGTACGATAGCCC"], new MafftOptions { SequenceType = MafftSequenceType.Amino }).AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(new[] { "ACGTACGTTAGCNNRY", "ACGTACGATAGCCC" }));
        });
    }

    [Test]
    public void EndGaps_AreFree_TheShorterSequenceSitsWhereItScoresBest()
    {
        var rows = SequenceAligner.MultipleAlignMafft(["TTTTTTTTTTACGTACGTAGTTTTTTTTTT", "ACGTACGTAG"]).AlignedSequences;
        Assert.That(rows[1], Is.EqualTo("----------ACGTACGTAG----------"));
    }

    [Test]
    public void Contracts()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.MultipleAlignMafft([]).AlignedSequences, Is.Empty);
            Assert.That(SequenceAligner.MultipleAlignMafft(["ACGT"]).AlignedSequences, Is.EqualTo(new[] { "ACGT" }), "A3-4: one sequence is written back (mafft 7.505 'Only 1 sequence found' is a warning)");
            Assert.That(SequenceAligner.MultipleAlignMafft(["ACGTACGT", "ACGAACGT", "ACGGACGT"]).AlignedSequences, Has.Length.EqualTo(3), "three sequences: the progressive stage (A2-19b, SequenceAligner_MafftProgressive_Tests)");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft(["ACGT", "---"]), "no residues");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft(["ACGT", "AC=GT"]), "'=' stops the program with --anysymbol");
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft(["ACGT", "ACGé"]));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft(["ACGT", "ACGA"], new MafftOptions { Blosum = 63 }));
            Assert.Throws<ArgumentException>(() => SequenceAligner.MultipleAlignMafft(["ACGT", "ACGA"], new MafftOptions { GapOpen = double.NaN }));
            Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignMafft(null!));
            Assert.Throws<OperationCanceledException>(() => SequenceAligner.MultipleAlignMafft(["ACGT", "ACGA"], new MafftOptions { CancellationToken = new CancellationToken(true) }));
        });
    }

    // B09 A3-4: a single sequence is written back by the program (installed mafft 7.505 --quiet [--anysymbol] on one record): gaps dropped, plain mode nucleotides lower / amino acids upper case.
    [TestCase("MKVLAAGI", "MKVLAAGI", "MKVLAAGI")]
    [TestCase("mkvLAAgi", "MKVLAAGI", "mkvLAAgi")]
    [TestCase("AC-GT.AC", "acgt.ac", "ACGT.AC")]
    [TestCase("ACGTNNRYacgt", "acgtnnryacgt", "ACGTNNRYacgt")]
    [TestCase("AC1GT", "acgt", "AC1GT")]
    [TestCase("MKV*LA", "MKVLA", "MKV*LA")]
    [TestCase("----", "", "")]
    public void SingleSequence_IsWrittenBackLikeTheProgram(string input, string plain, string anySymbol)
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAligner.MultipleAlignMafft([input], new MafftOptions { AnySymbol = false }).AlignedSequences, Is.EqualTo(new[] { plain }), "mafft --quiet");
            Assert.That(SequenceAligner.MultipleAlignMafft([input]).AlignedSequences, Is.EqualTo(new[] { anySymbol }), "mafft --quiet --anysymbol");
        });
    }
}
