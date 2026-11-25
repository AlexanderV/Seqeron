// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;

namespace Seqeron.Genomics.TCoffee;

using Seqeron.Genomics.ClustalW;
using static Seqeron.Genomics.ClustalW.ClustalWAligner;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

/// <summary>
/// The result of <see cref="TCoffeeAligner.TCoffeeEvaluate"/>: the consistency scores of the T-Coffee 13.41 evaluation (<c>t_coffee -infile aln -mode evaluate</c>, evaluate mode <c>triplet</c>): a score 0 - 9 per residue,
/// a score per column (the <c>cons</c> row), a score per sequence (0 - 100) and a score of the alignment (0 - 100).
/// </summary>
public sealed class TCoffeeEvaluation
{
    /// <summary>The sequence names in the order of the alignment.</summary>
    public required IReadOnlyList<string> Names { get; init; }

    /// <summary>The aligned rows, as given.</summary>
    public required IReadOnlyList<string> Rows { get; init; }

    /// <summary>The score 0 - 9 of every residue; -1 for a gap (the program's <c>NO_COLOR_RESIDUE</c>), per row and column.</summary>
    public required int[][] ResidueScores { get; init; }

    /// <summary>The score 0 - 9 of every column (the <c>cons</c> row); -1 when the column has no score.</summary>
    public required int[] ColumnScores { get; init; }

    /// <summary>The score 0 - 100 of every sequence.</summary>
    public required int[] SequenceScores { get; init; }

    /// <summary>The score 0 - 100 of the alignment (<c>SCORE=</c>).</summary>
    public required int AlignmentScore { get; init; }

    private const string Header = "T-COFFEE, Version_13.41.0.28bdc39 (2019-11-30 10:21:53 - Revision 5d5a1c1 - Build 465)\nCedric Notredame \nCPU TIME:{0} sec.\n";

    private static char Digit(int score) => score switch
    {
        < 0 => '-',
        > 9 => '#',
        _ => (char)('0' + score),
    };

    /// <summary>The <c>score_ascii</c> file (<c>output_reliability_format</c> with the ascii callbacks of <c>io_func.c</c>): the header, <c>SCORE=</c>, the per-sequence scores, the blocks of the residue scores with the <c>cons</c> row.</summary>
    /// <param name="lineLength">The number of characters per line of the program's ascii format (50 by default, cut to a multiple of three).</param>
    /// <param name="cpuSeconds">The CPU time printed in the header (the program prints its own).</param>
    public string ToScoreAscii(int lineLength = 50, int cpuSeconds = 0)
    {
        int nseq = Names.Count, len = ColumnScores.Length;
        int maxLen = Names.Max(n => n.Length);
        int width = lineLength - (lineLength % 3);
        if (width <= 0)
            throw new ArgumentOutOfRangeException(nameof(lineLength));
        var sb = new StringBuilder();
        sb.Append(string.Format(CultureInfo.InvariantCulture, Header, cpuSeconds));
        sb.Append("SCORE=").Append(AlignmentScore.ToString(CultureInfo.InvariantCulture)).Append("\n*\n BAD AVG GOOD\n*\n");
        string Label(string name) => (name.Length > maxLen ? name[..maxLen] : name).PadRight(maxLen + 2) + " ";
        for (int a = 0; a < nseq; a++)
            sb.Append(Label(Names[a])).Append(": ").Append(SequenceScores[a].ToString(CultureInfo.InvariantCulture).PadLeft(3)).Append('\n');
        sb.Append(Label("cons")).Append(": ").Append((AlignmentScore / 10).ToString(CultureInfo.InvariantCulture).PadLeft(3)).Append('\n');
        sb.Append('\n');
        for (int start = 0; start < len; start += width)
        {
            int end = Math.Min(start + width, len);
            for (int b = 0; b <= nseq; b++)
            {
                sb.Append(Label(b == nseq ? "cons" : Names[b]));
                for (int c = start; c < end; c++)
                    sb.Append(b == nseq ? Digit(ColumnScores[c]) : Digit(ResidueScores[b][c]));
                sb.Append('\n');
            }

            sb.Append("\n\n");
        }

        sb.Append("\n\n\n");
        return sb.ToString();
    }
}

public static partial class TCoffeeAligner
{
    /// <summary>
    /// The evaluation of an alignment that <c>t_coffee -infile aln -mode evaluate -output score_ascii</c> (T-Coffee 13.41, evaluate mode <c>triplet</c>, <c>-n_core 1</c>) computes: the library of the sequences of the alignment
    /// (<c>proba_pair</c> unless <paramref name="options"/> say otherwise) is relaxed; a residue scores <c>trunc(10 x S / M)</c> (at most 9) where <c>M</c> is 1000 plus the sum of <c>min(w(r, r2), w(r2, r3))</c> over the two-step paths of the library
    /// starting at the residue and <c>S</c> is 1000 plus the part of that sum whose end residue lies in the column of the residue; the column and the sequence scores are the same ratios over the column and the sequence, the alignment
    /// score is the ratio over all residues (x 100). The rows equal the file the program writes (<see cref="TCoffeeEvaluation.ToScoreAscii"/>).
    /// </summary>
    /// <param name="names">The sequence names (distinct); the program sorts the sequences by name before it builds the library.</param>
    /// <param name="rows">The aligned rows: letters and <c>-</c>, equal lengths, no empty sequence, no column of gaps only.</param>
    /// <param name="options">The library options of the program (default: <c>proba_pair</c>, one relaxation).</param>
    public static TCoffeeEvaluation TCoffeeEvaluate(IReadOnlyList<string> names, IReadOnlyList<string> rows, TCoffeeProgramOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(names);
        ArgumentNullException.ThrowIfNull(rows);
        options ??= new TCoffeeProgramOptions();
        int k = rows.Count;
        if (names.Count != k)
            throw new ArgumentException("Names must have one entry per row.", nameof(names));
        if (k == 0)
            throw new ArgumentException("The alignment has no sequence.", nameof(rows));
        int len = rows[0].Length;
        foreach (var r in rows)
        {
            ArgumentNullException.ThrowIfNull(r);
            if (r.Length != len)
                throw new ArgumentException("The rows of an alignment must have the same length.", nameof(rows));
            if (r.Any(c => c != '-' && !char.IsAsciiLetter(c)))
                throw new ArgumentException("Rows must consist of letters and '-'.", nameof(rows));
        }

        if (names.Distinct(StringComparer.Ordinal).Count() != k)
            throw new ArgumentException("Sequence names must be distinct.", nameof(names));
        var ungapped = rows.Select(r => r.Replace("-", string.Empty, StringComparison.Ordinal)).ToArray();
        RequireProgramMembers(ungapped, "t_coffee", nameof(rows));
        if (k == 1)
            throw new TCoffeeAbortedException(TCoffeeAbortedException.NoMethodMessage);
        for (int c = 0; c < len; c++)
        {
            if (rows.All(r => r[c] == '-'))
                throw new ArgumentException("The alignment has a column of gaps only.", nameof(rows));
        }

        TcRefuseModes(options, k);

        // -inorder aligned (the program's default): the sequences are sorted by name (byte order) before anything is computed; the alignment keeps its order
        var perm = Enumerable.Range(0, k).OrderBy(i => names[i], StringComparer.Ordinal).ToArray();
        var inverse = new int[k];
        for (int i = 0; i < k; i++)
            inverse[perm[i]] = i;
        var seqs = perm.Select(i => ungapped[i]).ToList();
        var kind = TcKindOf(seqs, options.Protein);
        var cl = TcBuildChecked(seqs, options.Methods, options.LalignTopAlignments, kind, options.Matrix, options.MethodParameters);
        if (cl.EntryCount == 0)
            throw new TCoffeeAbortedException(TCoffeeAbortedException.NoMethodMessage + " (no library entry: a method without suitable sequences produced no pair)");
        var index = TcIndex.FromConstraintList(cl, seqs.Select(s => s.Length).ToArray());
        if (options.FilterLib > 0)
            index.Filter(options.FilterLib);
        int relax = options.RelaxLibIterations ?? (options.RelaxLibrary ? 1 : 0);
        for (int a = 0; a < relax; a++)
            index.Relax();

        // aln2pos_simple: the 1-based residue number of every cell, 0 for a gap
        var pos = new int[k][];
        for (int b = 0; b < k; b++)
        {
            pos[b] = new int[len];
            int n = 0;
            for (int c = 0; c < len; c++)
                pos[b][c] = rows[b][c] == '-' ? 0 : ++n;
        }

        const double NormF = 1000;
        var lookup = new bool[k][];
        for (int s = 0; s < k; s++)
            lookup[s] = new bool[len + 2];
        var maxSeq = new double[k];
        var scoreSeq = new double[k];
        double scoreAln = 0, maxAln = 0;
        var residueScores = Enumerable.Range(0, k).Select(_ => new int[len]).ToArray();
        var columnScores = new int[len];
        for (int a = 0; a < len; a++)
        {
            double scoreCol = 0, maxCol = 0;
            for (int b = 0; b < k; b++)
            {
                if (pos[b][a] > 0)
                    lookup[inverse[b]][pos[b][a]] = true;
            }

            for (int b = 0; b < k; b++)
            {
                double scoreRes = NormF, maxRes = NormF;
                int s1 = inverse[b], r1 = pos[b][a];
                if (r1 <= 0)
                {
                    residueScores[b][a] = -1;
                    continue;
                }

                foreach (var e2 in index.Edges[s1][r1])
                {
                    foreach (var e3 in index.Edges[e2.Seq][e2.Res])
                    {
                        double score = Math.Min(e2.Weight, e3.Weight);
                        maxRes += score;
                        maxCol += score;
                        maxSeq[b] += score;
                        maxAln += score;
                        if (lookup[e3.Seq][e3.Res])
                        {
                            scoreRes += score;
                            scoreCol += score;
                            scoreSeq[b] += score;
                            scoreAln += score;
                        }
                    }
                }

                residueScores[b][a] = Math.Min((int)(scoreRes * 10 / maxRes), 9);
            }

            for (int b = 0; b < k; b++)
            {
                if (pos[b][a] > 0)
                    lookup[inverse[b]][pos[b][a]] = false;
            }

            columnScores[a] = maxCol == 0 ? -1 : Math.Min((int)(scoreCol * 10 / maxCol), 9);
        }

        return new TCoffeeEvaluation
        {
            Names = [.. names],
            Rows = [.. rows],
            ResidueScores = residueScores,
            ColumnScores = columnScores,
            SequenceScores = [.. Enumerable.Range(0, k).Select(a => maxSeq[a] == 0 ? 0 : (int)(scoreSeq[a] * 100 / maxSeq[a]))],
            AlignmentScore = maxAln == 0 ? 0 : (int)(scoreAln * 100 / maxAln),
        };
    }
}
