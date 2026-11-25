// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>
/// MCP <c>global_align_co_optimal</c> (B09 F59, audit A1-7b) and <c>multiple_align_exact</c> (B09 F60, audit A1-11). Expected values: Biopython 1.88
/// <c>PairwiseAligner(mode='global').align</c> set and exhaustive enumeration (core tests, scripts/alignment/generate_global_cooptimal_oracle.py);
/// exact optimum of two independent Python programs and exhaustive enumeration (scripts/alignment/generate_msa_exact_oracle.py).
/// </summary>
[TestFixture]
public class GlobalCoOptimalAndExactMsaTests
{

    [Test]
    public void MultipleAlignExact_BeatsOrEqualsHeuristics_AndReportsTheProblemSize()
    {
        string[] seqs = ["ACGTTGCA", "ACGTGCA", "ACTTGCAA", "ACGTTCA"];
        foreach (var model in new[] { "linear", "affine" })
        {
            var exact = AlignmentTools.MultipleAlignExact(seqs, 2, -3, -5, -2, model);
            var star = AlignmentTools.MultipleAlign(seqs, 2, -3, -5, -2, model);
            var prog = MultipleAlignmentTools.MultipleAlignProgressive(seqs, 2, -3, -5, -2, gapModel: model);
            Assert.That(exact.TotalScore, Is.GreaterThanOrEqualTo(star.TotalScore), model);
            Assert.That(exact.TotalScore, Is.GreaterThanOrEqualTo(prog.TotalScore), model);
            Assert.That(exact.AlignedSequences.Select(r => r.Replace("-", "")), Is.EqualTo(seqs));
            Assert.That(exact.LowerBound, Is.LessThanOrEqualTo(exact.TotalScore));
            Assert.That(exact.LatticeCells, Is.EqualTo(9L * 8 * 9 * 8));
            Assert.That(exact.VisitedCells, Is.LessThanOrEqualTo(exact.LatticeCells));
            // the full lattice gives the same optimum
            Assert.That(AlignmentTools.MultipleAlignExact(seqs, 2, -3, -5, -2, model, carrilloLipman: false).TotalScore, Is.EqualTo(exact.TotalScore), model);
        }
    }}
