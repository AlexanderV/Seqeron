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
using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>MCP <c>write_big_alignment</c> / <c>parse_big_alignment</c> (B09 A2-27d, F234 - F236) and <c>tcoffee_evaluate</c> (F237).</summary>
[TestFixture]
public class BigFormatTests
{

    [Test]
    public void TCoffeeEvaluate_AsTheProgram()
    {
        var r = MultipleAlignmentTools.TCoffeeEvaluate(["a", "b", "c"], ["ACGTACGTTAGCTAGCTAGGATCGATCGA-T", "ACGTACGATAGCTAGCTA-GATCGATCGATC", "ACGTTCGTTAGCTCGCTAGGATCGATCGA-T"]);
        const string expected = "T-COFFEE, Version_13.41.0.28bdc39 (2019-11-30 10:21:53 - Revision 5d5a1c1 - Build 465)\nCedric Notredame \nCPU TIME:0 sec.\nSCORE=97\n*\n BAD AVG GOOD\n*\n"
            + "a   :  97\nb   :  97\nc   :  97\nc   :   9\n\n"
            + "a   99999999999999999967999999999-9\nb   999999999999999999-799999999979\nc   99999999999999999967999999999-9\nc   9999999999999999995699999999939\n\n\n\n\n\n";
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignmentScore, Is.EqualTo(97));
            Assert.That(r.SequenceScores, Is.EqualTo(new[] { 97, 97, 97 }));
            Assert.That(r.ScoreAscii, Is.EqualTo(expected));
            Assert.That(() => MultipleAlignmentTools.TCoffeeEvaluate(["a"], ["ACGT"]), Throws.Exception);
        });
    }}
