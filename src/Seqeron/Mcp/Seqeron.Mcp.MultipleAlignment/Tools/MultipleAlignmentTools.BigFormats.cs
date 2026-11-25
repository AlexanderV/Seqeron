// SPDX-License-Identifier: GPL-3.0-or-later
// This file is part of Seqeron.Mcp.MultipleAlignment, licensed under the terms stated in LICENSE (this directory).
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Alignment.Tools;
using static Seqeron.Mcp.Alignment.Tools.AlignmentTools;

namespace Seqeron.Mcp.MultipleAlignment.Tools;

public partial class MultipleAlignmentTools
{
    [McpServerTool(Name = "tcoffee_evaluate", Title = "MSA - T-Coffee 13.41 evaluation of an alignment (score_ascii)", ReadOnly = true)]
    [Description("The consistency evaluation of an alignment that t_coffee -infile aln -mode evaluate -output score_ascii (T-Coffee 13.41, Notredame, Higgins & Heringa 2000; evaluate mode triplet, -n_core 1) computes: the library of the sequences of the alignment (proba_pair by default) is relaxed, and a residue scores trunc(10 x S / M), at most 9, where M is 1000 plus the sum of min(w(r, r2), w(r2, r3)) over the two-step library paths starting at the residue and S is 1000 plus the part of that sum whose end residue lies in the column of the residue; the column scores (the cons row), the sequence scores (0 - 100) and the alignment score (SCORE=, 0 - 100) are the same ratios over the column, the sequence and the whole alignment. Returns the scores and the score_ascii file (header, SCORE=, the per-sequence scores, blocks of 48 columns with the residue scores 0 - 9, - for a gap, and the cons row). The sequences are sorted by name before the library is built, as the program does (-inorder aligned). Library options as in multiple_align_tcoffee ('options' is a t_coffee flag string, e.g. '-method slow_pair -relax_lib 0'). VERIFICATION: the score_ascii text equals the installed t_coffee on 300 / 300 committed alignments (DNA and protein, 2 - 8 rows, t_coffee and random alignments, lower case; scripts/alignment/generate_tcoffee_score_oracle.py).")]
    public static TCoffeeEvaluateResult TCoffeeEvaluate(
        [Description("The sequence names (distinct), one per row.")] string[] names,
        [Description("The aligned rows: letters and '-', equal lengths, no column of gaps only.")] string[] rows,
        [Description("Library method: proba_pair (default), slow_pair, lalign_id_pair or any combination.")] string[]? methods = null,
        [Description("Further t_coffee flags (library options such as -relax_lib N, -filter_lib T, -matrix NAME, -lalign_n_top N).")] string? options = null,
        [Description("The line length of the ascii format (50 by default, cut to a multiple of three).")] int lineLength = 50)
    {
        ArgumentNullException.ThrowIfNull(names);
        ArgumentNullException.ThrowIfNull(rows);
        var tcOptions = new global::Seqeron.Genomics.TCoffee.TCoffeeProgramOptions();
        if (methods is { Length: > 0 })
            tcOptions = tcOptions.WithCommandLine("-method " + string.Join(",", methods.Select(m => (m ?? "").Trim().ToLowerInvariant())));
        if (!string.IsNullOrWhiteSpace(options))
            tcOptions = tcOptions.WithCommandLine(options);
        var e = global::Seqeron.Genomics.TCoffee.TCoffeeAligner.TCoffeeEvaluate(names, rows, tcOptions);
        return new TCoffeeEvaluateResult(e.AlignmentScore, e.SequenceScores, e.ColumnScores, e.ResidueScores, e.ToScoreAscii(lineLength));
    }
}
