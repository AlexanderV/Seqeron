using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit R14-8, B09 F301): all 30 named matrices of Biopython 1.88 (<c>Bio.Align.substitution_matrices.load</c>) are built in: BLASTN, BLASTP, FENG, HOXD70, LEVIN, MCLACHLAN, MEGABLAST, RAO, STR, TRANS
/// (integer) and BENNER6 / 22 / 74, DAYHOFF, GENETIC, GONNET1992, JOHNSON, JONES, MDM78, RISLER (fractional: stored exactly as integers in units of 1 / <see cref="SubstitutionMatrix.Denominator"/>, scored exactly by the aligners)
/// next to the BLOSUM / PAM / NUC.4.4 ones, plus the codon matrix SCHNEIDER and <see cref="SubstitutionMatrix.Parse"/>. Oracle TestData/Alignment/biopython_matrices_oracle.tsv.gz
/// (scripts/alignment/generate_matrix_oracle.py): every cell of the 29 residue matrices and 841 optimal global / local scores of seeded pairs with integer and fractional affine gap scores (Biopython <c>PairwiseAligner</c>).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class SubstitutionMatrix_Biopython_Tests
{
    private static readonly Lazy<List<string[]>> Oracle = new(() =>
    {
        using var stream = TestDataResource.Open("Seqeron.Genomics.Tests.TestData.Alignment.biopython_matrices_oracle.tsv")
            ?? throw new InvalidOperationException("biopython_matrices_oracle.tsv not embedded");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var rows = new List<string[]>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            rows.Add(line.Split('\t'));
        return rows;
    });

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    [Test]
    public void EveryCellOfEveryNamedMatrix_EqualsBiopython()
    {
        var cells = Oracle.Value.Where(r => r[0] == "M").ToList();
        Assert.That(cells, Has.Count.EqualTo(29));
        foreach (var r in cells)
        {
            Assert.That(SubstitutionMatrix.TryGet(r[1], out var m), Is.True, r[1]);
            string alphabet = r[2];
            var values = r[3].Split(' ');
            for (int i = 0; i < alphabet.Length; i++)
            {
                for (int j = 0; j < alphabet.Length; j++)
                {
                    Assert.That(m.TryGetValue(alphabet[i], alphabet[j], out decimal v), Is.True, $"{r[1]} {alphabet[i]}{alphabet[j]}");
                    Assert.That((double)v, Is.EqualTo(D(values[(i * alphabet.Length) + j])).Within(1e-9), $"{r[1]} {alphabet[i]}{alphabet[j]}");
                }
            }
        }
    }

    [Test]
    public void OptimalScores_WithTheNamedMatrices_EqualBiopythonPairwiseAligner()
    {
        var rows = Oracle.Value.Where(r => r[0] == "A").ToList();
        Assert.That(rows, Has.Count.EqualTo(812));
        foreach (var r in rows)
        {
            var scoring = new MatrixScoring(SubstitutionMatrix.Get(r[1]), D(r[3]), D(r[4]));
            var result = r[2] == "global" ? SequenceAligner.GlobalAlignAffine(r[5], r[6], scoring) : SequenceAligner.LocalAlignAffine(r[5], r[6], scoring);
            Assert.That(result.ExactScore, Is.EqualTo(D(r[7])).Within(1e-6), $"{r[1]} {r[2]} gap {r[3]}/{r[4]} {r[5]} / {r[6]}");
        }

        Assert.That(rows.Select(r => r[1]).Distinct().Count(), Is.EqualTo(29));
    }

    [Test]
    public void FractionalMatrices_AreStoredExactly_AndTheIntegerAccessorsRefuseThem()
    {
        var m = SubstitutionMatrix.Get("GONNET1992");
        Assert.That(m.IsFractional, Is.True);
        Assert.That(m.Denominator, Is.EqualTo(10));
        Assert.That(SubstitutionMatrix.Get("JOHNSON").Denominator, Is.EqualTo(100));
        Assert.That(SubstitutionMatrix.Get("BLOSUM62").Denominator, Is.EqualTo(1));
        Assert.That(m.GetValue('W', 'W'), Is.EqualTo(14.2m));
        Assert.That(m.IsPositive('W', 'W'), Is.True);
        Assert.That(m.MaxAbsValue, Is.EqualTo(14.2m));
        Assert.Throws<InvalidOperationException>(() => m.GetScore('A', 'A'));
        Assert.Throws<InvalidOperationException>(() => _ = m.MaxAbsScore);
        Assert.That(SubstitutionMatrix.Get("BLASTP").GetScore('U', 'C'), Is.EqualTo(9));
        // a whole-number table built with decimals is the integer matrix
        var whole = new SubstitutionMatrix("W", "AC", new decimal[,] { { 2m, -1m }, { -1m, 2m } });
        Assert.That((whole.IsFractional, whole.GetScore('A', 'C')), Is.EqualTo((false, -1)));
        Assert.Throws<ArgumentException>(() => new SubstitutionMatrix("T", "AC", new decimal[,] { { 0.0000001m, 0 }, { 0, 1 } }));
    }

    [Test]
    public void BlastpMatrix_ScoresUOJ_AsBlastpDoes_WhileBlosum62KeepsRejectingThem()
    {
        // blastp 2.12.0 -ungapped (measured): U is scored as C (U/C 9, U/A 0), O as X (-1 with every residue, -4 with *), J = the NCBI J column (J/I 3, J/L 3, J/M 2, J/A -1)
        var blastp = SubstitutionMatrix.Get("BLASTP");
        Assert.That((blastp.GetScore('U', 'U'), blastp.GetScore('U', 'A'), blastp.GetScore('U', 'C'), blastp.GetScore('O', 'W'), blastp.GetScore('O', '*'), blastp.GetScore('J', 'I'), blastp.GetScore('J', 'M'), blastp.GetScore('J', 'A')), Is.EqualTo((9, 0, 9, -1, -4, 3, 2, -1)));
        Assert.That(SubstitutionMatrix.Get("BLOSUM62").Contains('U'), Is.False, "Biopython / EMBOSS BLOSUM62 has no U, O, J: the named BLOSUM matrices are unchanged");
        Assert.That(SequenceAligner.GlobalAlignAffine("MKTAUYIAK", "MKTAAYIAK", new MatrixScoring(blastp, -11, -1)).Score, Is.EqualTo(39), "blastp scores MKTAUYIAK / MKTAAYIAK 39: U/A contributes 0");
    }

    [Test]
    public void SchneiderCodonMatrix_EqualsBiopython()
    {
        var m = CodonSubstitutionMatrix.Schneider;
        Assert.That(m.Codons, Has.Count.EqualTo(64));
        Assert.That(m.GetValue("AAA", "AAA"), Is.EqualTo(11.6m));
        Assert.That(m.GetValue("AAA", "AAC"), Is.EqualTo(-2.7m));
        Assert.That(m.GetValue("aaa", "aag"), Is.EqualTo(9.7m));
        Assert.That(m.GetValue("AAU", "AAT"), Is.EqualTo(m.GetValue("AAT", "AAT")), "U is read as T");
        Assert.That(m.TryGetValue("AAA", "AAN", out _), Is.False);
    }

    private const string NcbiBlosum62Fragment = "# comment\n   A  R  N\nA  4 -1 -2\nR -1  5  0\nN -2  0  6\n";

    [Test]
    public void Parse_ReadsTheNcbiEmbossAndBiopythonFormat_AndRoundTrips()
    {
        var m = SubstitutionMatrix.Parse(NcbiBlosum62Fragment, "FRAG");
        Assert.That((m.Name, m.Alphabet, m.IsFractional), Is.EqualTo(("FRAG", "ARN", false)));
        Assert.That((m.GetScore('A', 'A'), m.GetScore('R', 'N'), m.GetScore('N', 'A')), Is.EqualTo((4, 0, -2)));
        // rows in another order are placed by their symbol; lower case folds; decimals make a fractional matrix
        var f = SubstitutionMatrix.Parse("  a  c\nc 0.5 1.25\na 2 -0.5\n");
        Assert.That((f.IsFractional, f.Denominator, f.GetValue('A', 'C'), f.GetValue('C', 'A'), f.GetValue('C', 'C')), Is.EqualTo((true, 100, -0.5m, 0.5m, 1.25m)));
        // every built-in matrix survives text -> matrix -> text
        foreach (string name in SubstitutionMatrix.BiopythonNames.Concat(["BLOSUM62", "PAM250", "EDNAFULL", "BLOSUM45"]))
        {
            var original = SubstitutionMatrix.Get(name);
            var back = SubstitutionMatrix.Parse(original.ToText(), name);
            Assert.That(back.Alphabet, Is.EqualTo(original.Alphabet), name);
            foreach (char a in original.Alphabet)
            {
                foreach (char b in original.Alphabet)
                    Assert.That(back.GetValue(a, b), Is.EqualTo(original.GetValue(a, b)), $"{name} {a}{b}");
            }
        }
    }

    [Test]
    public void Parse_Refusals()
    {
        Assert.Throws<ArgumentException>(() => SubstitutionMatrix.Parse(""));
        Assert.Throws<ArgumentException>(() => SubstitutionMatrix.Parse("A C\nA 1 2\n"), "a row is missing");
        Assert.Throws<ArgumentException>(() => SubstitutionMatrix.Parse("A C\nA 1 2\nC 3\n"), "a short row");
        Assert.Throws<ArgumentException>(() => SubstitutionMatrix.Parse("A C\nA 1 x\nC 3 4\n"), "not a number");
        Assert.Throws<ArgumentException>(() => SubstitutionMatrix.Parse("A A\nA 1 2\nA 3 4\n"), "duplicate symbol");
        Assert.Throws<ArgumentException>(() => SubstitutionMatrix.Parse("AAA CCC\nAAA 1 2\nCCC 3 4\n"), "a codon header");
        Assert.Throws<ArgumentNullException>(() => SubstitutionMatrix.Parse(null!));
    }

    [Test]
    public void FractionalMatrix_WithTheAligners_AllModesAndTheIntegerOnlyConsumers()
    {
        var scoring = new MatrixScoring(SubstitutionMatrix.Get("JOHNSON"), -2.5, -0.5);
        var global = SequenceAligner.GlobalAlignAffine("ACDEFGHIK", "ACDFGHIKL", scoring);
        var local = SequenceAligner.LocalAlignAffine("ACDEFGHIK", "ACDFGHIKL", scoring);
        Assert.That(global.ExactScore, Is.Not.EqualTo(Math.Round(global.ExactScore)), "the score is fractional (two decimals)");
        Assert.That(local.ExactScore, Is.GreaterThanOrEqualTo(global.ExactScore));
        var stats = SequenceAligner.CalculateStatistics(global, scoring.Matrix);
        Assert.That(stats.Matches, Is.GreaterThan(0));
        // EMBOSS needle reads integer matrices only (its data files are integer): a clear refusal, not a wrong score
        Assert.Throws<InvalidOperationException>(() => SequenceAligner.NeedleAlign("ACDEFGHIK", "ACDFGHIKL", SubstitutionMatrix.Get("JOHNSON")));
    }
}
