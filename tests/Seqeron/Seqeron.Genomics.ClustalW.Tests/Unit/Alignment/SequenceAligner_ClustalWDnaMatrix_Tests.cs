// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using System.Globalization;
using System.Text.Json;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A1-10d part 3, B09 F68): the nucleotide matrices of ClustalW 2.1 - <c>-DNAMATRIX=IUB | CLUSTALW | file</c>, <c>-PWDNAMATRIX=IUB | CLUSTALW | file</c> and
/// <c>-TRANSWEIGHT</c> - in <see cref="ClustalWAligner.MultipleAlignClustalW(IEnumerable{string}, ClustalWOptions?, ScoringMatrix?)"/> and
/// <see cref="ClustalWAligner.ClustalWPairwiseAlign(string, string, ClustalWDnaPairwiseOptions)"/>.
///
/// Oracle (TestData/Alignment/clustalw_dna_matrix_oracle.tsv, scripts/alignment/generate_clustalw_dna_matrix_oracle.py): the installed <c>clustalw</c> 2.1, driven two ways - (cli) the
/// command line <c>-TYPE=DNA -ALIGN -OUTPUT=FASTA -OUTORDER=INPUT -DNAMATRIX=x -PWDNAMATRIX=y</c>: IUB / CLUSTALW are used, a FILE is read and validated but then registered under the PROTEIN
/// slot (a defect of <c>CommandLineParser</c>), so the DNA stages keep IUB (<see cref="ClustalWUserDnaMatrixUse.CommandLine"/>); and (menu) the interactive menu fed through stdin, which registers
/// the file for DNA, so it is used (<see cref="ClustalWUserDnaMatrixUse.Applied"/>). The files (embedded in the TSV) come in many layouts; every output row, every printed pairwise
/// "Score:" (full stage) and every printed "Group g: Sequences: n Score: s" / "Delayed" line of the first pass must be identical.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ClustalWDnaMatrix_Tests
{
    private sealed record Row(string Id, string Mode, string[] Seqs, string Options, Dictionary<string, string> Files, string[]? Rows,
        (int I, int J, int Score)[] PairScores, (int Group, int Entries, int? Score)[] Groups);

    private static readonly Lazy<List<Row>> Oracle = new(() =>
    {
        const string name = "Seqeron.Genomics.Tests.TestData.Alignment.clustalw_dna_matrix_oracle.tsv";
        using var stream = System.Reflection.Assembly.GetExecutingAssembly().GetManifestResourceStream(name)
            ?? throw new InvalidOperationException($"Embedded oracle '{name}' not found.");
        using var reader = new StreamReader(stream);
        reader.ReadLine();
        var list = new List<Row>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Row(f[0], f[1], f[2].Split(';'), f[3], JsonSerializer.Deserialize<Dictionary<string, string>>(f[4])!,
                f[5] == "EXIT" ? null : f[5].Split(';'),
                f[6].Length == 0 ? [] : f[6].Split(';').Select(p => p.Split(':')).Select(p => (int.Parse(p[0]), int.Parse(p[1]), int.Parse(p[2]))).ToArray(),
                f[7].Length == 0 ? [] : f[7].Split(';').Select(p => p.Split(':')).Select(p => (int.Parse(p[0]), p.Length == 2 ? 0 : int.Parse(p[1]), p.Length == 2 ? (int?)null : int.Parse(p[2]))).ToArray()));
        }

        return list;
    });

    private static double D(string s) => double.Parse(s, CultureInfo.InvariantCulture);

    private static ClustalWDnaMatrix Builtin(string s) => s.ToUpperInvariant() switch
    {
        "IUB" => ClustalWDnaMatrix.Iub,
        "CLUSTALW" => ClustalWDnaMatrix.Clustalw,
        _ => throw new ArgumentException(s),
    };

    /// <summary>The oracle's option string as <see cref="ClustalWOptions"/>; reading a malformed matrix file throws <see cref="FormatException"/> like the program's "wrong format" exit.</summary>
    private static ClustalWOptions ParseOptions(Row r)
    {
        var o = new ClustalWOptions
        {
            UserDnaMatrixUse = r.Mode == "cli" ? ClustalWUserDnaMatrixUse.CommandLine : ClustalWUserDnaMatrixUse.Applied,
        };
        if (r.Options.Length == 0)
            return o;
        foreach (string kv in r.Options.Split(','))
        {
            var p = kv.Split('=');
            switch (p[0])
            {
                case "gapopen": o = o with { GapOpen = D(p[1]) }; break;
                case "gapext": o = o with { GapExtend = D(p[1]) }; break;
                case "transweight": o = o with { TransitionWeight = float.Parse(p[1], CultureInfo.InvariantCulture) }; break;
                case "maxdiv": o = o with { DivergenceCutoff = int.Parse(p[1]) }; break;
                case "negative": o = o with { NegativeMatrix = true }; break;
                case "quick": o = o with { PairwiseMode = ClustalWPairwiseMode.Quick }; break;
                case "pwgapopen": o = o with { PairwiseGapOpen = D(p[1]) }; break;
                case "pwgapext": o = o with { PairwiseGapExtend = D(p[1]) }; break;
                case "dnamatrix":
                    o = r.Files.TryGetValue(p[1], out var text) ? o with { UserDnaMatrix = ClustalWUserMatrix.Parse(text, p[1]) } : o with { DnaMatrix = Builtin(p[1]) };
                    break;
                case "pwdnamatrix":
                    o = r.Files.TryGetValue(p[1], out var pwText) ? o with { PairwiseUserDnaMatrix = ClustalWUserMatrix.Parse(pwText, p[1]) } : o with { PairwiseDnaMatrix = Builtin(p[1]) };
                    break;
                default: throw new InvalidOperationException($"unknown oracle option {p[0]}");
            }
        }

        return o;
    }

    private static ClustalWDnaPairwiseOptions PairwiseOf(ClustalWOptions o) => new()
    {
        GapOpen = o.PairwiseGapOpen,
        GapExtend = o.PairwiseGapExtend,
        TransitionWeight = o.TransitionWeight,
        Matrix = o.PairwiseDnaMatrix,
        UserMatrix = o.PairwiseUserDnaMatrix,
        UserMatrixUse = o.UserDnaMatrixUse,
    };

    /// <summary>
    /// The profile stage indexes its matrix with Xref[0], [2], [4], [11] and [12] (the letters at those positions of the file's header: the transition weights). A header with fewer than
    /// 13 residues leaves some of them -1 and ClustalW writes outside its <c>matrix[NUMRES][NUMRES]</c> array: undefined behaviour (8 of the 37 such menu runs of the generator crashed with
    /// SIGSEGV and are not in the oracle; the port drops the writes). Applies only when the user matrix is used for the profile stage.
    /// </summary>
    private static readonly string[] UndefinedBehaviourIds = ["0", "6", "8", "18", "26", "40", "74", "154", "172", "181", "189", "194", "202", "219", "224", "228", "229", "230", "232", "285", "289", "303", "304", "314", "326", "343", "349", "356", "358"];

    private static readonly string[] UndefinedBehaviourDifferingIds = ["303"];

    private static bool IsUndefinedBehaviour(Row r)
    {
        var o = ParseOptions(r);
        if (o.UserDnaMatrix is null || o.UserDnaMatrixUse != ClustalWUserDnaMatrixUse.Applied)
            return false;
        return new[] { 0, 2, 4, 11, 12 }.Any(p => p >= o.UserDnaMatrix.Cross.Count || o.UserDnaMatrix.Cross[p] < 0);
    }

    [Test]
    public void MultipleAlignClustalW_DnaMatrices_Rows_EqualClustalW21_OnAllOracleSets()
    {
        var rows = Oracle.Value.Where(r => r.Rows is not null).ToList();
        Assert.That(rows.Count, Is.GreaterThanOrEqualTo(150));
        var bad = new List<string>();
        var undefinedIds = new List<string>();
        var undefinedDiffer = new List<string>();
        foreach (var r in rows)
        {
            string[] got = ClustalWAligner.MultipleAlignClustalW(r.Seqs, ParseOptions(r)).AlignedSequences.ToArray();
            bool equal = got.SequenceEqual(r.Rows!);
            if (IsUndefinedBehaviour(r))
            {
                undefinedIds.Add(r.Id);
                if (!equal)
                    undefinedDiffer.Add(r.Id);
            }
            else if (!equal)
            {
                bad.Add($"{r.Id} {r.Mode} ({r.Options})");
            }
        }

        Assert.Multiple(() =>
        {
            Assert.That(bad, Is.Empty);
            // undefined behaviour of the C++ (writes outside the matrix array): the port - which drops those writes - is checked on exactly these recorded runs, and differs from the binary on exactly these
            Assert.That(undefinedIds, Is.EqualTo(UndefinedBehaviourIds), "the recorded runs whose user matrix header leaves an Xref index at -1; observed: " + string.Join(", ", undefinedIds.Select(x => "\"" + x + "\"")));
            Assert.That(undefinedDiffer, Is.EqualTo(UndefinedBehaviourDifferingIds), "of those, the runs on which the dropped writes change the binary's rows; observed: " + string.Join(", ", undefinedDiffer.Select(x => "\"" + x + "\"")));
        });
    }

    /// <summary>The first-pass "Group g: Sequences: n Score: s" lines carry the profile-alignment optimum of every step: they pin the profile matrix (and gap scale) numerically.</summary>
    [Test]
    public void MultipleAlignClustalW_DnaMatrices_GroupScores_EqualClustalW21_OnAllOracleSets()
    {
        int checkedGroups = 0;
        var bad = new List<string>();
        foreach (var r in Oracle.Value.Where(r => r.Rows is not null && !IsUndefinedBehaviour(r)))
        {
            var (_, groups) = ClustalWAligner.MultipleAlignClustalWWithGroupScores(r.Seqs, ParseOptions(r));
            // ClustalW prints "Delayed" without the entries when no profile alignment was done or its score is not positive
            var got = groups.Select(g => g.Score is null ? (g.Group, 0, (int?)null) : (g.Group, g.Entries, g.Score)).ToArray();
            checkedGroups += got.Length;
            if (!got.SequenceEqual(r.Groups))
                bad.Add($"{r.Id} {r.Mode} ({r.Options}): expected {string.Join(';', r.Groups)} got {string.Join(';', got)}");
        }

        Assert.That(checkedGroups, Is.GreaterThan(800));
        Assert.That(bad, Is.Empty);
    }

    /// <summary>The pairwise stage: the printed "Sequences (i:j) Aligned. Score:" of the binary (full stage only; the fast stage has no matrix).</summary>
    [Test]
    public void ClustalWPairwiseAlign_DnaMatrices_Scores_EqualClustalW21_OnAllOraclePairs()
    {
        int pairs = 0, withMatrix = 0, bad = 0;
        string? first = null;
        foreach (var r in Oracle.Value.Where(r => r.Rows is not null && !r.Options.Split(',').Contains("quick")))
        {
            var o = PairwiseOf(ParseOptions(r));
            if (o.UserMatrix is not null && o.UserMatrixUse == ClustalWUserDnaMatrixUse.Applied || o.Matrix == ClustalWDnaMatrix.Clustalw)
                withMatrix++;
            foreach (var (i, j, score) in r.PairScores)
            {
                pairs++;
                var got = ClustalWAligner.ClustalWPairwiseAlign(r.Seqs[i - 1], r.Seqs[j - 1], o);
                if (got.ClustalWScore != score && ++bad == 1)
                    first = $"set {r.Id} {r.Mode} pair {i}:{j} ({r.Options}): expected {score}, got {got.ClustalWScore}";
            }
        }

        Assert.That(pairs, Is.GreaterThan(1500));
        Assert.That(withMatrix, Is.GreaterThan(60));
        Assert.That(bad, Is.Zero, first);
    }

    /// <summary>A malformed file ends clustalw ("wrong format in matrix", exit status 1, no output): the reader rejects it with a <see cref="FormatException"/>.</summary>
    [Test]
    public void MalformedMatrixFiles_AreRejectedLikeClustalW()
    {
        var exits = Oracle.Value.Where(r => r.Rows is null).ToList();
        Assert.That(exits.Count, Is.GreaterThanOrEqualTo(10));
        foreach (var r in exits)
            Assert.Throws<FormatException>(() => ParseOptions(r), r.Id);
    }

    private static bool IsShortHeader(string line)
    {
        if (line.Length == 0 || line.StartsWith('#'))
            return false;
        var words = line.Trim().Split(new[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
        return words.Length is >= 4 and <= 12 && words.All(w => w.Length == 1);
    }

    [Test]
    public void Oracle_CoversEveryModeMatrixAndLayout()
    {
        var rows = Oracle.Value;
        var files = rows.SelectMany(r => r.Files.Values).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(rows.Count, Is.GreaterThanOrEqualTo(330));
            Assert.That(rows.Count(r => r.Mode == "cli"), Is.GreaterThan(100));
            Assert.That(rows.Count(r => r.Mode == "menu"), Is.GreaterThan(150));
            Assert.That(rows.Count(r => r.Mode == "cli" && r.Options.Contains("dnamatrix=") && r.Options.Contains("CLUSTALW", StringComparison.OrdinalIgnoreCase)), Is.GreaterThan(30), "CLUSTALW matrix");
            Assert.That(rows.Count(r => r.Mode == "cli" && r.Rows is not null && r.Files.Count > 0), Is.GreaterThan(40), "files on the command line");
            Assert.That(rows.Count(r => r.Mode == "menu" && r.Options.Contains("dnamatrix=")), Is.GreaterThan(100), "user profile matrix (menu)");
            Assert.That(rows.Count(r => r.Mode == "menu" && r.Options.Contains("pwdnamatrix=")), Is.GreaterThan(100), "user pairwise matrix (menu)");
            foreach (string key in new[] { "gapopen", "gapext", "transweight", "maxdiv", "negative", "quick", "pwgapopen", "pwgapext" })
                Assert.That(rows.Count(r => r.Options.Split(',').Any(o => o.Split('=')[0] == key)), Is.GreaterThan(20), key);
            Assert.That(files.Count(t => t.Contains('\t')), Is.GreaterThan(10), "tab separated");
            Assert.That(files.Count(t => t.Contains('.') && t.Split('\n').Skip(3).Any(l => l.Contains(".0 ") || l.Contains(".5 "))), Is.GreaterThan(10), "floating-point rows");
            Assert.That(files.Count(t => t.StartsWith('#')), Is.GreaterThan(30), "comment lines");
            Assert.That(files.Count(t => t.Split('\n').Take(6).Any(IsShortHeader)), Is.GreaterThan(10), "short headers");
        });
    }

    /// <summary>
    /// The matrix choice decides the alignment: over the oracle sets the rows (or the group scores) differ from the IUB run of the same set on most
    /// user-matrix menu runs, and the command-line files change nothing.
    /// </summary>
    [Test]
    public void Oracle_DiscriminatesTheMatrixOptions()
    {
        int menuDiffers = 0, menuSets = 0, cliFileSame = 0, cliFileSets = 0, clustalwDiffers = 0, clustalwSets = 0, scoreDiffers = 0;
        foreach (var r in Oracle.Value.Where(r => r.Rows is not null))
        {
            var o = ParseOptions(r);
            var withoutFiles = o with { UserDnaMatrix = null, PairwiseUserDnaMatrix = null };
            var (rowsPlain, groupsPlain) = ClustalWAligner.MultipleAlignClustalWWithGroupScores(r.Seqs, withoutFiles);
            bool same = rowsPlain.SequenceEqual(r.Rows!);
            if (r.Mode == "menu")
            {
                menuSets++;
                if (!same)
                    menuDiffers++;
                if (!groupsPlain.Select(g => g.Score is null ? (g.Group, 0, (int?)null) : (g.Group, g.Entries, g.Score)).SequenceEqual(r.Groups))
                    scoreDiffers++;
            }
            else if (o.UserDnaMatrix is not null || o.PairwiseUserDnaMatrix is not null)
            {
                cliFileSets++;
                if (same)
                    cliFileSame++;   // the command line reads the file but keeps the built-in matrix: the rows are those of the run without the file
            }
            else if (o.DnaMatrix == ClustalWDnaMatrix.Clustalw || o.PairwiseDnaMatrix == ClustalWDnaMatrix.Clustalw)
            {
                clustalwSets++;
                var iub = o with { DnaMatrix = ClustalWDnaMatrix.Iub, PairwiseDnaMatrix = ClustalWDnaMatrix.Iub };
                if (!ClustalWAligner.MultipleAlignClustalW(r.Seqs, iub).AlignedSequences.SequenceEqual(r.Rows!))
                    clustalwDiffers++;
            }
        }

        Assert.Multiple(() =>
        {
            Assert.That(menuDiffers, Is.GreaterThan(menuSets / 3), "user matrices change the rows");
            Assert.That(scoreDiffers, Is.GreaterThan(menuSets * 2 / 3), "user matrices change the group scores");
            Assert.That(cliFileSame, Is.EqualTo(cliFileSets), "command-line files change nothing (the PROTEIN-slot defect)");
            Assert.That(clustalwDiffers, Is.GreaterThan(clustalwSets / 4), "CLUSTALW matrix changes the rows");
        });
    }

    [Test]
    public void CommandLineUse_IgnoresTheFile_AppliedUses_It()
    {
        var r = Oracle.Value.First(x => x.Mode == "menu" && x.Rows is not null && x.Options.Contains("dnamatrix=") && !x.Options.Contains("quick"));
        var applied = ParseOptions(r);
        var cli = applied with { UserDnaMatrixUse = ClustalWUserDnaMatrixUse.CommandLine };
        var plain = applied with { UserDnaMatrix = null, PairwiseUserDnaMatrix = null };
        Assert.Multiple(() =>
        {
            Assert.That(ClustalWAligner.MultipleAlignClustalW(r.Seqs, cli).AlignedSequences, Is.EqualTo(ClustalWAligner.MultipleAlignClustalW(r.Seqs, plain).AlignedSequences));
            Assert.That(ClustalWAligner.MultipleAlignClustalW(r.Seqs, applied).AlignedSequences, Is.EqualTo(r.Rows));
        });
    }

    /// <summary>A user matrix holding exactly the IUB table of <c>matrices.h</c> (in its packing order) behaves like the built-in IUB matrix (same cells, same Xref).</summary>
    [Test]
    public void UserMatrixEqualToIub_BehavesLikeTheBuiltInMatrix()
    {
        const string order = "ABCDGHKMNRSTUVWXY";
        string Identical(int match, int mismatch)
        {
            var sb = new System.Text.StringBuilder("# IUB\n");
            sb.Append(string.Join(' ', order.Select(c => c.ToString()))).Append(" *\n");
            for (int i = 0; i <= order.Length; i++)
            {
                sb.Append(i < order.Length ? order[i].ToString() : "*");
                for (int j = 0; j <= order.Length; j++)
                {
                    int v = i == order.Length || j == order.Length ? -4 : IubScore(order[i], order[j], match, mismatch);
                    sb.Append(' ').Append(v);
                }

                sb.Append('\n');
            }

            return sb.ToString();
        }

        var user = ClustalWUserMatrix.Parse(Identical(10, -9), "iub.mat");
        var rnd = new Random(7);
        for (int t = 0; t < 12; t++)
        {
            string[] seqs = Enumerable.Range(0, 5).Select(_ => string.Concat(Enumerable.Range(0, rnd.Next(30, 90)).Select(_ => "ACGT"[rnd.Next(4)]))).ToArray();
            var builtIn = new ClustalWOptions { TransitionWeight = 0.25f };
            var viaFile = builtIn with { UserDnaMatrix = user, PairwiseUserDnaMatrix = user };
            Assert.That(ClustalWAligner.MultipleAlignClustalW(seqs, viaFile).AlignedSequences, Is.EqualTo(ClustalWAligner.MultipleAlignClustalW(seqs, builtIn).AlignedSequences), $"set {t}");
        }

        // pairwise too
        var a = "ACGTACGGTACGTTAGC";
        var b = "ACGTTCGGACGTAGC";
        Assert.That(ClustalWAligner.ClustalWPairwiseAlign(a, b, new ClustalWDnaPairwiseOptions { UserMatrix = user }).ClustalWScore,
            Is.EqualTo(ClustalWAligner.ClustalWPairwiseAlign(a, b).ClustalWScore));
    }

    private static int IubScore(char a, char b, int match, int mismatch)
    {
        // swgapdnamt of matrices.h: +10 when the IUB sets overlap, -9 otherwise (verified against the table by the test above through equal behaviour)
        const string order = "ABCDGHKMNRSTUVWXY";
        int i = order.IndexOf(a), j = order.IndexOf(b);
        if (i < j)
            (i, j) = (j, i);
        return Triangle(i * (i + 1) / 2 + j) ? match : mismatch;
    }

    private static bool Triangle(int index) => IubTriangle[index] > 0;

    private static readonly short[] IubTriangle =
    {
        10,
        -9, 10,
        -9, 10, 10,
        10, 10, -9, 10,
        -9, 10, -9, 10, 10,
        10, 10, 10, 10, -9, 10,
        -9, 10, -9, 10, 10, 10, 10,
        10, 10, 10, 10, -9, 10, -9, 10,
        10, 10, 10, 10, 10, 10, 10, 10, 10,
        10, 10, -9, 10, 10, 10, 10, 10, 10, 10,
        -9, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
        -9, 10, -9, 10, -9, 10, 10, -9, 10, -9, -9, 10,
        -9, 10, -9, 10, -9, 10, 10, -9, 10, -9, -9, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, -9, -9, 10,
        10, 10, -9, 10, -9, 10, 10, 10, 10, 10, -9, 10, 10, 10, 10,
        10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
        -9, 10, 10, 10, -9, 10, 10, 10, 10, -9, 10, 10, 10, 10, 10, 10, 10,
    };

    [Test]
    public void PairwiseOptions_Defaults_AndContracts()
    {
        var o = new ClustalWDnaPairwiseOptions();
        Assert.Multiple(() =>
        {
            Assert.That(o.GapOpen, Is.EqualTo(15.0));
            Assert.That(o.GapExtend, Is.EqualTo(6.66));
            Assert.That(o.TransitionWeight, Is.EqualTo(0.5f));
            Assert.That(o.Matrix, Is.EqualTo(ClustalWDnaMatrix.Iub));
            Assert.That(o.UserMatrixUse, Is.EqualTo(ClustalWUserDnaMatrixUse.Applied));
            Assert.That(new ClustalWOptions().DnaMatrix, Is.EqualTo(ClustalWDnaMatrix.Iub));
            Assert.That(new ClustalWOptions().UserDnaMatrixUse, Is.EqualTo(ClustalWUserDnaMatrixUse.Applied));
            // the options overload with defaults = the classic overload
            Assert.That(ClustalWAligner.ClustalWPairwiseAlign("ACGTACGTTAGC", "ACGTCGTTAGGC", o).ClustalWScore, Is.EqualTo(ClustalWAligner.ClustalWPairwiseAlign("ACGTACGTTAGC", "ACGTCGTTAGGC").ClustalWScore));
            Assert.Throws<ArgumentNullException>(() => ClustalWAligner.ClustalWPairwiseAlign("ACGT", "ACGT", (ClustalWDnaPairwiseOptions)null!));
            Assert.Throws<ArgumentException>(() => ClustalWAligner.ClustalWPairwiseAlign("ACGT", "ACGT", o with { GapOpen = -1 }));
            Assert.Throws<ArgumentException>(() => ClustalWAligner.ClustalWPairwiseAlign("ACGT", "ACGT", o with { TransitionWeight = float.NaN }));
            Assert.Throws<ArgumentOutOfRangeException>(() => ClustalWAligner.ClustalWPairwiseAlign("ACGT", "ACGT", o with { Matrix = (ClustalWDnaMatrix)9 }));
            Assert.That(ClustalWAligner.ClustalWPairwiseAlign(new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT"), o with { Matrix = ClustalWDnaMatrix.Clustalw }).ClustalWScore, Is.EqualTo(100));
        });
    }
}
