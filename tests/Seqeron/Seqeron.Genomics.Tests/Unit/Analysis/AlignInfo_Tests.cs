using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// B09 WP67 (F308, R14-13): <see cref="SummaryInfo"/> / <see cref="Pssm"/> of Biopython 1.79 (<c>Bio.Align.AlignInfo</c>). Expected values: the tests of Biopython 1.79 itself (<c>Tests/test_AlignInfo.py</c>: the nucleotide, protein and
/// pseudo-count alignments; the example of the PSSM theory notes of the University of Marseille), the docstring examples and values derived by hand (each derivation is in the test). The recorded run of 6 199 calls on the real class is
/// the package handoff/R14-13-aligninfo-biopython179 (STATUS in its README).
/// </summary>
[TestFixture]
public class AlignInfo_Tests
{
    private static MultipleSeqAlignment Alignment(params string[] rows) => new(rows.Select((s, i) => new SeqRecord(s, $"ID{i + 1:000}")));

    private static Dictionary<char, double> Uniform(string letters) => letters.ToDictionary(c => c, _ => 1.0 / letters.Length);

    [Test]
    public void Nucleotides_BiopythonTestSuite()
    {
        // GFF/multi.fna: test1 ACGTCGCG, test2 GGGGCCCC, test3 AAACACAC
        var summary = new SummaryInfo(Alignment("ACGTCGCG", "GGGGCCCC", "AAACACAC"));
        string consensus = summary.DumbConsensus(ambiguous: "N");
        Assert.Multiple(() =>
        {
            Assert.That(consensus, Is.EqualTo("NNNNNNNN"));
            Assert.That(summary.GapConsensus(ambiguous: "N"), Is.EqualTo("NNNNNNNN"));
            Assert.That(summary.PosSpecificScoreMatrix(consensus, ["-"]).ToString(), Is.EqualTo(
                "    A   C   G   T\nN  2.0 0.0 1.0 0.0\nN  1.0 1.0 1.0 0.0\nN  1.0 0.0 2.0 0.0\nN  0.0 1.0 1.0 1.0\nN  1.0 2.0 0.0 0.0\nN  0.0 2.0 1.0 0.0\nN  1.0 2.0 0.0 0.0\nN  0.0 2.0 1.0 0.0\n"));
            Assert.That(summary.InformationContent(eFreqTable: new Dictionary<char, double> { ['A'] = 0.25, ['G'] = 0.25, ['T'] = 0.25, ['C'] = 0.25 }, charsToIgnore: ["-"]), Is.EqualTo(7.32029999423075).Within(1e-6));
        });
    }

    [Test]
    public void Proteins_BiopythonTestSuite()
    {
        var a = Alignment("MHQAIFIYQIGYP*LKSGYIQSIRSPEYDNW-", "MH--IFIYQIGYAYLKSGYIQSIRSPEY-NW*", "MHQAIFIYQIGYPYLKSGYIQSIRSPEYDNW*");
        Assert.That(a.GetAlignmentLength(), Is.EqualTo(32));
        var s = new SummaryInfo(a);
        string consensus = s.DumbConsensus(ambiguous: "X");
        Assert.Multiple(() =>
        {
            Assert.That(consensus, Is.EqualTo("MHQAIFIYQIGYXXLKSGYIQSIRSPEYDNW*"));
            Assert.That(s.GapConsensus(ambiguous: "X"), Is.EqualTo("MHXXIFIYQIGYXXLKSGYIQSIRSPEYXNWX"));
        });

        string pssm = s.PosSpecificScoreMatrix(s.GapConsensus(ambiguous: "X"), ["-", "*"]).ToString();
        string[] lines = pssm.Split('\n');
        Assert.Multiple(() =>
        {
            Assert.That(lines[0], Is.EqualTo("    A   D   E   F   G   H   I   K   L   M   N   P   Q   R   S   W   Y"));
            Assert.That(lines[1], Is.EqualTo("M  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"));
            Assert.That(lines[3], Is.EqualTo("X  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0"));
            Assert.That(lines[13], Is.EqualTo("X  1.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.0 0.0 0.0 0.0 0.0 0.0"));
            Assert.That(lines[29], Is.EqualTo("X  0.0 2.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"));
            Assert.That(lines[31], Is.EqualTo("W  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 3.0 0.0"));
            Assert.That(lines[32], Is.EqualTo("X  0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0"));
            Assert.That(lines, Has.Length.EqualTo(34));
            Assert.That(s.InformationContent(eFreqTable: Uniform("ACDEFGHIKLMNPQRSTVWY"), charsToIgnore: ["-", "*"]), Is.EqualTo(133.061475107).Within(1e-6));
        });
    }

    [Test]
    public void PseudoCount_PssmTheoryExample()
    {
        // http://biologie.univ-mrs.fr/upload/p202/01.4.PSSM_theory.pdf: eight 12-mers
        var summary = new SummaryInfo(Alignment("AACCACGTTTAA", "CACCACGTGGGT", "CACCACGTTCGC", "GCGCACGTGGGG", "TCGCACGTTGTG", "TGGCACGTGTTT", "TGACACGTGGGA", "TTACACGTGCGC"));
        double ic = summary.InformationContent(eFreqTable: new Dictionary<char, double> { ['A'] = 0.325, ['G'] = 0.175, ['T'] = 0.325, ['C'] = 0.175 }, logBase: Math.E, pseudoCount: 1);
        double[] expected = [0.110, 0.090, 0.360, 1.290, 0.800, 1.290, 1.290, 0.80, 0.610, 0.390, 0.470, 0.040];
        Assert.Multiple(() =>
        {
            Assert.That(ic, Is.EqualTo(7.546).Within(1e-3));
            Assert.That(summary.IcVector, Has.Count.EqualTo(12));
            for (int i = 0; i < 12; i++)
                Assert.That(summary.IcVector[i], Is.EqualTo(expected[i]).Within(0.005), $"column {i}");
            Assert.That(summary.IcVector.Sum(), Is.EqualTo(ic).Within(1e-12));
        });
    }

    [Test]
    public void ReplacementDictionary_WeightedExampleOfTheDocstring()
    {
        // GTATC 0.5 / AT--C 0.8 / CTGTC 1.0 (letters ACGT; the letter of the earlier row first), derived column by column:
        // c0 G,A,C: (G,A) .4, (G,C) .5, (A,C) .8; c1 T,T,T: (T,T) .4 + .5 + .8; c2 A,-,G: (A,G) .5; c3 T,-,T: (T,T) .5; c4 C,C,C: (C,C) .4 + .5 + .8
        var rows = new[] { ("GTATC", 0.5), ("AT--C", 0.8), ("CTGTC", 1.0) };
        var records = rows.Select((r, i) => new SeqRecord(r.Item1, $"r{i}", annotations: new JsonObject { ["weight"] = r.Item2 })).ToList();
        var rep = new SummaryInfo(new MultipleSeqAlignment(records)).ReplacementDictionary("ACGT");
        Assert.Multiple(() =>
        {
            Assert.That(rep, Has.Count.EqualTo(16));
            Assert.That(rep.Keys.First(), Is.EqualTo(('A', 'A')));
            Assert.That(rep[('G', 'A')], Is.EqualTo(0.4).Within(1e-12));
            Assert.That(rep[('G', 'C')], Is.EqualTo(0.5).Within(1e-12));
            Assert.That(rep[('A', 'C')], Is.EqualTo(0.8).Within(1e-12));
            Assert.That(rep[('T', 'T')], Is.EqualTo(2.2).Within(1e-12));
            Assert.That(rep[('A', 'G')], Is.EqualTo(0.5).Within(1e-12));
            Assert.That(rep[('C', 'C')], Is.EqualTo(1.7).Within(1e-12));
            Assert.That(rep[('C', 'G')], Is.EqualTo(0.0));
            Assert.Throws<ArgumentException>(() => new SummaryInfo(new MultipleSeqAlignment(records)).ReplacementDictionary("ACGT", "-"), "skip_chars is deprecated");
            Assert.Throws<InvalidOperationException>(() => new SummaryInfo(new MultipleSeqAlignment(records)).ReplacementDictionary(), "letters is required");
        });
    }

    [Test]
    public void InformationContent_HandDerivedColumn()
    {
        // one column A, A, C with expected 0.5 / 0.5 in bits: 2/3 log2(4/3) + 1/3 log2(2/3) = 0.08170416594551044 (derived by hand: 0.6666667 x 0.4150375 - 0.3333333 x 0.5849625)
        var s = new SummaryInfo(Alignment("A", "A", "C"));
        double ic = s.InformationContent(eFreqTable: new Dictionary<char, double> { ['A'] = 0.5, ['C'] = 0.5 });
        Assert.Multiple(() =>
        {
            Assert.That(ic, Is.EqualTo(2.0 / 3 * Math.Log2(4.0 / 3) + 1.0 / 3 * Math.Log2(2.0 / 3)).Within(1e-14));
            Assert.That(ic, Is.EqualTo(0.08170416594551044).Within(1e-12));
            Assert.That(s.IcVector, Is.EqualTo(new[] { ic }));
        });

        // a gap is not counted for the frequencies of the letters but it is in the total when it is not ignored: A, A, - gives f(A) = 2/3
        var gapped = new SummaryInfo(Alignment("A", "A", "-"));
        Assert.That(gapped.InformationContent(eFreqTable: new Dictionary<char, double> { ['A'] = 0.5, ['C'] = 0.5 }), Is.EqualTo(2.0 / 3 * Math.Log2(4.0 / 3) + 0).Within(1e-12));
        Assert.That(gapped.InformationContent(eFreqTable: new Dictionary<char, double> { ['A'] = 0.5, ['C'] = 0.5 }, charsToIgnore: ["-"]), Is.EqualTo(Math.Log2(2.0)).Within(1e-12));
    }

    [Test]
    public void InformationContent_ErrorsAsBiopython()
    {
        var s = new SummaryInfo(Alignment("AC", "AG"));
        var table = new Dictionary<char, double> { ['A'] = 0.25, ['C'] = 0.25, ['G'] = 0.25 };
        Assert.Multiple(() =>
        {
            // no expected frequencies: random_expected is None in Biopython and the division fails
            Assert.That(Assert.Throws<InvalidOperationException>(() => s.InformationContent())!.Message, Is.EqualTo("unsupported operand type(s) for /: 'float' and 'NoneType'"));
            Assert.That(Assert.Throws<ArgumentException>(() => s.InformationContent(eFreqTable: table, end: 3))!.Message, Is.EqualTo("Start (0) and end (3) are not in the range 0 to 2"));
            Assert.That(Assert.Throws<ArgumentException>(() => s.InformationContent(eFreqTable: table, pseudoCount: -1))!.Message, Is.EqualTo("Positive value required for pseudo_count, -1 provided"));
            Assert.That(Assert.Throws<ArgumentException>(() => s.InformationContent(eFreqTable: new Dictionary<char, double> { ['A'] = 0.5 }))!.Message, Is.EqualTo("C not found in expected frequency table"));
            Assert.That(Assert.Throws<DivideByZeroException>(() => s.InformationContent(eFreqTable: table, logBase: 1))!.Message, Is.EqualTo("float division by zero"));
            Assert.Throws<InvalidOperationException>(() => s.PosSpecificScoreMatrix("ACG"), "the axis must have the alignment length (AssertionError)");
        });
    }

    [Test]
    public void PosSpecificScoreMatrix_MutatesCharsToIgnore_AndPrintsWeights()
    {
        // Biopython appends '-' to the list it is given
        var ignore = new List<string>();
        var records = new[] { new SeqRecord("A-", "a", annotations: new JsonObject { ["weight"] = 0.25 }), new SeqRecord("AC", "b") };
        var pssm = new SummaryInfo(new MultipleSeqAlignment(records)).PosSpecificScoreMatrix(null, ignore);
        Assert.Multiple(() =>
        {
            Assert.That(ignore, Is.EqualTo(new[] { "-" }));
            Assert.That(pssm.Count, Is.EqualTo(2));
            Assert.That(pssm[0]['A'], Is.EqualTo(1.25));
            Assert.That(pssm[1]['C'], Is.EqualTo(1.0));
            Assert.That(pssm.GetResidue(0), Is.EqualTo("A"));
            Assert.That(pssm.ToString(), Is.EqualTo("    A   C\nA  1.2 0.0\nC  0.0 1.0\n"), "1.25 prints as 1.2 (round half to even of the exact binary value)");
        });
    }

    [Test]
    public void PrintInfoContent_And_GetColumn()
    {
        var s = new SummaryInfo(Alignment("AC", "AG"));
        Assert.That(s.GetColumn(1), Is.EqualTo("CG"));
        s.InformationContent(eFreqTable: new Dictionary<char, double> { ['A'] = 0.25, ['C'] = 0.25, ['G'] = 0.25 });
        // column 0: A conserved: 1 x log2(1 / 0.25) = 2; column 1: C and G: 2 x 0.5 x log2(0.5 / 0.25) = 1
        Assert.That(s.PrintInfoContent(), Is.EqualTo("0 A 2.000\n1 C 1.000\n"));
        Assert.Throws<InvalidOperationException>(() => new SummaryInfo(Alignment("AC", "AG")).PrintInfoContent(), "ic_vector is empty and information_content() has no table");
    }
}
