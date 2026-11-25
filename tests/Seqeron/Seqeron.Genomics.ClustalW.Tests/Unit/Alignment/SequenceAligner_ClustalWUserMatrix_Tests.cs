// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10d, file part, B09 F64): user-defined matrix files of ClustalW 2.1 - <c>-MATRIX=file</c> (one matrix or a CLUSTAL_SERIES file),
/// <c>-PWMATRIX=file</c> - read by <see cref="ClustalWUserMatrix.Parse"/> / <see cref="ClustalWMatrixFile.Parse"/> (ports of <c>SubMatrix::readUserMatrix</c> /
/// <c>readMatrixSeries</c>) and used by <see cref="ClustalWAligner.MultipleAlignClustalW(IEnumerable{string}, ClustalWProteinOptions, MatrixScoring?)"/> and
/// <see cref="ClustalWAligner.ClustalWPairwiseAlign(string, string, ClustalWProteinPairwiseOptions)"/>.
///
/// Oracle (TestData/Alignment/clustalw_user_matrix_oracle.tsv, scripts/alignment/generate_clustalw_user_matrix_oracle.py): the installed <c>clustalw</c> 2.1
/// <c>-TYPE=PROTEIN -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT -MATRIX=file / -PWMATRIX=file</c> on 240 sets whose matrix files (embedded in the TSV) are written in
/// NCBI / EMBOSS layout (comments, labelled rows, '*' row / column), unlabelled rows, integer / "%.1f" / "%g" floating-point rows (x10, scale per row), tab separated,
/// header subsets, shuffled residue order, junk in the upper triangle, trailers, readUserMatrix quirks (J and lower-case header letters, 2-3 residue matrices) and
/// CLUSTAL_SERIES files with 2-4 entries (gaps / overlaps in the limits); every output row and every printed pairwise "Score:" must be identical.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWUserMatrix_Tests
{
    private sealed record Row(string Id, string[] Seqs, string[] Rows, string Options, Dictionary<string, string> Files, (int I, int J, int Score)[] PairScores);

    private static IEnumerable<string[]> ReadTsv(string resource)
    {
        string name = "Seqeron.Genomics.Tests.TestData.Alignment." + resource;
        using var stream = TestDataResource.Open(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        string? line;
        while ((line = reader.ReadLine()) is not null)
            if (line.Length > 0)
                yield return line.Split('\t');
    }

    private static readonly Lazy<List<Row>> Oracle = new(() =>
        ReadTsv("clustalw_user_matrix_oracle.tsv").Select(f => new Row(
            f[0], f[1].Split(';'), f[2].Split(';'), f[3],
            JsonSerializer.Deserialize<Dictionary<string, string>>(f[4])!,
            f[5].Length == 0 ? [] : f[5].Split(';').Select(p => p.Split(':')).Select(p => (int.Parse(p[0]), int.Parse(p[1]), int.Parse(p[2]))).ToArray())).ToList());

    /// <summary>The oracle's option string: the matrix options name files of the set, everything else is the shared clustalw option syntax.</summary>
    private static ClustalWProteinOptions ParseOptions(Row r)
    {
        var others = new List<string>();
        ClustalWUserMatrix? pw = null;
        ClustalWMatrixFile? file = null;
        foreach (string kv in r.Options.Length == 0 ? [] : r.Options.Split(','))
        {
            var p = kv.Split('=');
            if (p[0] == "matrix")
                file = ClustalWMatrixFile.Parse(r.Files[p[1]], name => r.Files[name], p[1]);
            else if (p[0] == "pwmatrix")
                pw = ClustalWUserMatrix.Parse(r.Files[p[1]], p[1]);
            else
                others.Add(kv);
        }

        var o = SequenceAligner_ClustalWProteinMultiple_Tests.Parse(string.Join(',', others));
        return o with
        {
            UserMatrix = file?.Matrix,
            UserMatrixSeries = file?.Series,
            Pairwise = o.Pairwise with { UserMatrix = pw },
        };
    }

    [Test]
    public void MultipleAlignClustalW_UserMatrixFiles_Rows_EqualClustalW21_OnAllOracleSets()
    {
        var rows = Oracle.Value;
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(240));
        var bad = new List<string>();
        foreach (var r in rows)
        {
            string[] got = ClustalWAligner.MultipleAlignClustalW(r.Seqs, ParseOptions(r)).AlignedSequences.ToArray();
            if (!got.SequenceEqual(r.Rows))
                bad.Add($"{r.Id} ({r.Options})");
        }

        Assert.That(bad, Is.Empty);
    }

    /// <summary>The pairwise stage with <c>-PWMATRIX=file</c>: the printed "Sequences (i:j) Aligned. Score:" of the binary (the full stage; not printed with -QUICKTREE).</summary>
    [Test]
    public void ClustalWPairwiseAlign_UserMatrixFile_Scores_EqualClustalW21_OnAllOraclePairs()
    {
        int pairs = 0, bad = 0, withUserPw = 0;
        string? first = null;
        foreach (var r in Oracle.Value.Where(r => !r.Options.Split(',').Contains("quick")))
        {
            var o = ParseOptions(r);
            if (o.Pairwise.UserMatrix is not null)
                withUserPw++;
            foreach (var (i, j, score) in r.PairScores)
            {
                pairs++;
                var got = ClustalWAligner.ClustalWPairwiseAlign(r.Seqs[i - 1], r.Seqs[j - 1], o.Pairwise);
                if (got.ClustalWScore != score && ++bad == 1)
                    first = $"set {r.Id} pair {i}:{j} ({r.Options}): expected {score}, got {got.ClustalWScore}";
            }
        }

        Assert.That(pairs, Is.GreaterThan(2000));
        Assert.That(withUserPw, Is.GreaterThan(40));
        Assert.That(bad, Is.Zero, first);
    }

    [Test]
    public void UserMatrixOracle_CoversEveryFormatAndPath()
    {
        var rows = Oracle.Value;
        var files = rows.SelectMany(r => r.Files.Where(f => f.Key != "series.txt").Select(f => f.Value)).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(rows.Count(r => r.Options.Contains("matrix=series.txt")), Is.GreaterThan(60), "series files");
            Assert.That(rows.Count(r => r.Options.Split(',').Any(o => o.StartsWith("matrix=m"))), Is.GreaterThan(30), "single -MATRIX files");
            Assert.That(rows.Count(r => r.Options.Contains("pwmatrix=")), Is.GreaterThan(60), "-PWMATRIX files");
            Assert.That(files.Count(f => f.Contains("0.5") || f.Contains(".0 ")), Is.GreaterThan(40), "floating-point matrices");
            Assert.That(files.Count(f => f.Contains('\t')), Is.GreaterThan(20), "tab separated");
            Assert.That(files.Count(f => f.StartsWith('#')), Is.GreaterThan(40), "comment header");
            Assert.That(rows.Count(r => r.Options.Contains("quick")), Is.GreaterThan(10), "fast pairwise stage");
            Assert.That(rows.Count(r => r.Options.Contains("negative")), Is.GreaterThan(10), "-NEGATIVE");
        });
    }

    /// <summary>The oracle discriminates: the same sets aligned with the built-in matrix give other rows on many sets, so the user matrix really drives the result.</summary>
    [Test]
    public void UserMatrixFiles_ChangeRowsAgainstBuiltInMatrices_OnManyOracleSets()
    {
        int differs = 0, total = 0;
        foreach (var r in Oracle.Value.Where(r => r.Options.Contains("matrix=series.txt") || r.Options.Contains("matrix=m")))
        {
            total++;
            var o = ParseOptions(r);
            var builtin = o with { UserMatrix = null, UserMatrixSeries = null };
            if (!ClustalWAligner.MultipleAlignClustalW(r.Seqs, builtin).AlignedSequences.SequenceEqual(r.Rows))
                differs++;
        }

        Assert.That(total, Is.GreaterThan(100));
        Assert.That(differs, Is.GreaterThan(total / 5), "user -MATRIX files change the rows");
    }

    // ---- reader unit tests ---------------------------------------------------------------------------------------------------------------------

    private const string Header = "   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *";

    private static string BlastText(SubstitutionMatrix m, bool labels = true, string comment = "# BLOSUM62\n")
    {
        var sb = new System.Text.StringBuilder(comment);
        sb.Append(Header).Append('\n');
        string letters = "ARNDCQEGHILKMFPSTWYVBZX*";
        foreach (char a in letters)
        {
            if (labels)
                sb.Append(a).Append(' ');
            sb.Append(string.Join(' ', letters.Select(b => a == '*' || b == '*' ? "-4" : m.GetScore(a, b).ToString())));
            sb.Append('\n');
        }

        return sb.ToString();
    }

    [Test]
    public void Parse_NcbiLayout_GivesTheCanonicalMatrix()
    {
        var u = ClustalWUserMatrix.Parse(BlastText(SubstitutionMatrix.Blosum62), "B62");
        var m = u.ToSubstitutionMatrix();
        Assert.Multiple(() =>
        {
            Assert.That(u.Residues, Is.EqualTo("ARNDCQEGHILKMFPSTWYVBZX"));
            Assert.That(u.RecognisedResidues, Is.EqualTo(23));
            Assert.That(m.Alphabet, Is.EqualTo("ARNDCQEGHILKMFPSTWYVBZX"));
            foreach (char a in m.Alphabet)
                foreach (char b in m.Alphabet)
                    Assert.That(m.GetScore(a, b), Is.EqualTo(SubstitutionMatrix.Blosum62.GetScore(a, b)), $"{a}{b}");
        });
    }

    [Test]
    public void Parse_LabelledAndUnlabelledRows_AreTheSameMatrix()
    {
        var withLabels = ClustalWUserMatrix.Parse(BlastText(SubstitutionMatrix.Blosum62, labels: true));
        var without = ClustalWUserMatrix.Parse(BlastText(SubstitutionMatrix.Blosum62, labels: false, comment: ""));
        Assert.That(without.Cells(), Is.EqualTo(withLabels.Cells()));
        Assert.That(without.Cross, Is.EqualTo(withLabels.Cross));
    }

    [Test]
    public void Parse_FloatingPointRows_AreScaledByTenAndTruncated_PerRow()
    {
        // row scale = 10 when the FIRST value token of the row contains a '.', else 1 (so "2 1.5" reads 2 and 1)
        const string text = "A C G *\n" +
                            "A 2.5 0 0 0\n" +
                            "C -0.5 3.0 0 0\n" +
                            "G 1 0.7 4 0\n" +
                            "* 0 0 0 0\n";
        var u = ClustalWUserMatrix.Parse(text);
        Assert.That(u.Values.Take(6), Is.EqualTo(new short[] { 25, -5, 30, 1, 0, 4 }));
        var m = u.ToSubstitutionMatrix();
        Assert.That(new[] { m.GetScore('A', 'A'), m.GetScore('C', 'A'), m.GetScore('A', 'C'), m.GetScore('G', 'C') }, Is.EqualTo(new[] { 25, -5, -5, 0 }));
    }

    [Test]
    public void Parse_OnlyTheLowerTriangleIsRead()
    {
        var a = ClustalWUserMatrix.Parse("A C G *\nA 5 9 9 9\nC 1 6 9 9\nG 2 3 7 9\n* 0 0 0 0\n").ToSubstitutionMatrix();
        var b = ClustalWUserMatrix.Parse("A C G *\nA 5 -7 -7 0\nC 1 6 -7 0\nG 2 3 7 0\n* 0 0 0 0\n").ToSubstitutionMatrix();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Assert.That(a.ScoreAt(i, j), Is.EqualTo(b.ScoreAt(i, j)));
        Assert.That(a.GetScore('A', 'C'), Is.EqualTo(1));
        Assert.That(a.GetScore('C', 'A'), Is.EqualTo(1));
    }

    [Test]
    public void Parse_CommentsBlankLinesAndTrailer_AreSkipped()
    {
        const string text = "# c1\n\n   \n# c2\nA C G *\n\nA 5 0 0 0\nC 1 6 0 0\n\nG 2 3 7 0\n* 0 0 0 0\n! trailer\nnot read\n";
        var u = ClustalWUserMatrix.Parse(text);
        Assert.That(u.Values.Take(6), Is.EqualTo(new short[] { 5, 1, 6, 2, 3, 7 }));
    }

    [Test]
    public void Parse_MalformedFiles_AreRejectedLikeClustalW()
    {
        // every message below is what the installed clustalw 2.1 reports for the same file ("ERROR: wrong format in matrix", "single matrix expected")
        Assert.Multiple(() =>
        {
            Assert.Throws<FormatException>(() => ClustalWUserMatrix.Parse("A C G\nA 5 0 0\nC 1 6 0\nG 2 3 7\n"), "no row for the non-letter header column (rows must be letters + 1)");
            Assert.Throws<FormatException>(() => ClustalWUserMatrix.Parse("A C G *\nA 5 0 0 0\nC 1 6 0 0\nG 2 3 7 0\n"), "no '*' row");
            Assert.Throws<FormatException>(() => ClustalWUserMatrix.Parse("# only a comment\n"), "no header");
            Assert.Throws<FormatException>(() => ClustalWUserMatrix.Parse("CLUSTAL_SERIES\nMATRIX 0 100 x\n"), "a series where one matrix is expected");
            Assert.Throws<FormatException>(() => ClustalWUserMatrix.Parse("A C G *\nA 5\nC 1 6 0 0\nG 2 3 7 0\n* 0 0 0 0\n"), "a row with fewer tokens than residues");
        });
    }

    [Test]
    public void Parse_UnknownHeaderLetters_AreAcceptedWithAShiftedTable_AsInClustalW()
    {
        // clustalw only warns "residue J in matrix ... not recognised"; Cross keeps -1 and the canonical form does not exist
        var u = ClustalWUserMatrix.Parse("A C J G *\nA 5 0 0 0 0\nC 1 6 0 0 0\nJ 2 3 8 0 0\nG 2 3 4 7 0\n* 0 0 0 0 0\n");
        Assert.That(u.Cross, Is.EqualTo(new[] { 0, 2, -1, 6 }));
        Assert.That(u.RecognisedResidues, Is.EqualTo(3));
        Assert.Throws<ArgumentException>(() => u.ToSubstitutionMatrix());
    }

    [Test]
    public void MatrixFile_Series_ParsesEntriesAndFollowsClustalWsRules()
    {
        var files = new Dictionary<string, string>
        {
            ["a.mat"] = BlastText(SubstitutionMatrix.Blosum62),
            ["b.mat"] = BlastText(SubstitutionMatrix.Pam250),
        };
        var f = ClustalWMatrixFile.Parse("# series\nCLUSTAL_SERIES\n# comment\nMATRIX 0 40 a.mat\nsome ignored line\nMATRIX 40 100 b.mat\n", n => files[n], "s.txt");
        Assert.Multiple(() =>
        {
            Assert.That(f.Matrix, Is.Null);
            Assert.That(f.Series!.Select(e => (e.LowerLimit, e.UpperLimit)), Is.EqualTo(new[] { (0, 40), (40, 100) }));
            Assert.That(f.Series![1].Matrix.ToSubstitutionMatrix().GetScore('A', 'A'), Is.EqualTo(SubstitutionMatrix.Pam250.GetScore('A', 'A')));
            Assert.Throws<FormatException>(() => ClustalWMatrixFile.Parse("CLUSTAL_SERIES\nMATRIX 50 40 a.mat\n", n => files[n]), "ulimit <= llimit");
            Assert.Throws<FormatException>(() => ClustalWMatrixFile.Parse("CLUSTAL_SERIES\nMATRIX 0 101 a.mat\n", n => files[n]), "limit above 100");
            Assert.Throws<FormatException>(() => ClustalWMatrixFile.Parse("CLUSTAL_SERIES\nMATRIX 0 x a.mat\n", n => files[n]), "not three fields");
            Assert.Throws<FormatException>(() => ClustalWMatrixFile.Parse("CLUSTAL_SERIES\n# nothing\n", n => files[n]), "no entry");
            // a file that is not a series is one matrix
            Assert.That(ClustalWMatrixFile.Parse(files["a.mat"]).Matrix, Is.Not.Null);
        });
    }

    [Test]
    public void MatrixFile_Series_ReadsAtMostTenMatrices()
    {
        var text = "CLUSTAL_SERIES\n" + string.Concat(Enumerable.Range(0, 12).Select(i => $"MATRIX {i} {i + 1} a.mat\n"));
        var f = ClustalWMatrixFile.Parse(text, _ => BlastText(SubstitutionMatrix.Blosum62));
        Assert.That(f.Series!.Count, Is.EqualTo(10));
    }

    /// <summary>A canonical <see cref="SubstitutionMatrix"/> passed through <see cref="ClustalWUserMatrix.FromSubstitutionMatrix"/> is the matrix ClustalW reads from the equivalent file.</summary>
    [Test]
    public void FromSubstitutionMatrix_EqualsTheParsedFile_AndDrivesTheAlignmentIdentically()
    {
        var fromObject = ClustalWUserMatrix.FromSubstitutionMatrix(SubstitutionMatrix.Blosum62);
        var fromFile = ClustalWUserMatrix.Parse(BlastText(SubstitutionMatrix.Blosum62));
        Assert.That(fromObject.Cells(), Is.EqualTo(fromFile.Cells()));
        Assert.That(fromObject.Cross, Is.EqualTo(fromFile.Cross));

        var set = Oracle.Value.First(r => r.Options.Length == 0 || !r.Options.Contains("quick"));
        var viaObject = ClustalWAligner.MultipleAlignClustalW(set.Seqs, new ClustalWProteinOptions { UserMatrix = fromObject, Pairwise = new ClustalWProteinPairwiseOptions { UserMatrix = fromObject } });
        var viaFile = ClustalWAligner.MultipleAlignClustalW(set.Seqs, new ClustalWProteinOptions { UserMatrix = fromFile, Pairwise = new ClustalWProteinPairwiseOptions { UserMatrix = fromFile } });
        Assert.That(viaObject.AlignedSequences, Is.EqualTo(viaFile.AlignedSequences));
    }

    /// <summary>
    /// Series selection and scale: the entry whose [llimit, ulimit] contains the percent identity (the first entry if none does) with scale
    /// 0.5 + (pcid - llimit) / ((ulimit - llimit) * 2): a one-entry series covering every identity equals the single matrix only for the scale 0.75 (pcid 50 in 0..100),
    /// so the two options are different code paths on the oracle - asserted through the rows of the oracle (above) and, here, by the entry choice.
    /// </summary>
    [Test]
    public void SeriesAndSingleMatrix_AreDistinctPaths_OnTheOracle()
    {
        int differs = 0;
        foreach (var r in Oracle.Value.Where(r => r.Options.Contains("matrix=series.txt")).Take(60))
        {
            var o = ParseOptions(r);
            var single = o with { UserMatrixSeries = null, UserMatrix = o.UserMatrixSeries![0].Matrix };
            if (!ClustalWAligner.MultipleAlignClustalW(r.Seqs, single).AlignedSequences.SequenceEqual(r.Rows))
                differs++;
        }

        Assert.That(differs, Is.GreaterThan(5));
    }
}
