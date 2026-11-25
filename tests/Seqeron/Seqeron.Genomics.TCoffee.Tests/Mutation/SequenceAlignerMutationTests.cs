// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using Seqeron.Genomics.TCoffee;
namespace Seqeron.Genomics.Tests.Mutation;

/// <summary>
/// Targeted mutation-killing tests for SequenceAligner.cs (checklist 04 rows 35/36/37/38/226:
/// ALIGN-GLOBAL/LOCAL/SEMI/MULTI-001 + ALIGN-STATS-001).
///
/// The canonical suite under-covered the cancellable Needleman-Wunsch overload (a full
/// re-implementation), the EMBOSS srspair formatter, the statistics percentages, and the
/// two multiple-sequence aligners (star + progressive). These pin the published algorithm
/// outputs with independent ground truth so operator/relational/logical mutants diverge.
/// </summary>
[TestFixture]
public class SequenceAlignerMutationTests
{

    // ── Star MSA (MultipleAlign): center select, reconcile, consensus, sum-of-pairs ───

    private static int SumOfPairs(string[] rows, ScoringMatrix s)
    {
        int len = rows[0].Length, total = 0;
        for (int p = 0; p < len; p++)
            for (int i = 0; i < rows.Length; i++)
                for (int j = i + 1; j < rows.Length; j++)
                {
                    char ci = rows[i][p], cj = rows[j][p];
                    if (ci != '-' || cj != '-')
                    {
                        if (ci == '-' || cj == '-') total += s.GapExtend;
                        else if (ci == cj) total += s.Match;
                        else total += s.Mismatch;
                    }
                }
        return total;
    }

    private static string DeGap(string s) => s.Replace("-", "");

    // ── Progressive MSA (MultipleAlignProgressive): guide tree + profile alignment ────

    [Test]
    public void MultipleAlignProgressive_IdenticalSequences_ExactConsensusAndSumOfPairs()
    {
        var seqs = new[] { new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT"), new DnaSequence("ACGTACGT") };
        var r = Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs, SequenceAligner.SimpleDna);

        r.AlignedSequences.Should().OnlyContain(s => s == "ACGTACGT");
        r.Consensus.Should().Be("ACGTACGT");
        r.TotalScore.Should().Be(24);
    }

    [Test]
    public void MultipleAlignProgressive_DivergentSequences_RowsValid_AndTotalScoreIsTrueSumOfPairs()
    {
        var inputs = new[] { "ACGTACGT", "ACGTTACGT", "ACGACGT", "ACGTACGA" };
        var seqs = inputs.Select(s => new DnaSequence(s)).ToArray();
        var r = Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(seqs, SequenceAligner.SimpleDna);

        r.AlignedSequences.Should().HaveCount(4);
        r.AlignedSequences.Select(x => x.Length).Distinct().Should().ContainSingle();
        // rows are reprojected to input order
        for (int i = 0; i < inputs.Length; i++)
            DeGap(r.AlignedSequences[i]).Should().Be(inputs[i]);
        r.TotalScore.Should().Be(SumOfPairs(r.AlignedSequences, SequenceAligner.SimpleDna));
    }

    [Test]
    public void MultipleAlignProgressive_EmptyInput_ReturnsEmpty()
    {
        Seqeron.Genomics.ClustalW.ClustalWAligner.MultipleAlignProgressive(System.Array.Empty<DnaSequence>())
            .Should().Be(MultipleAlignmentResult.Empty);
    }}
