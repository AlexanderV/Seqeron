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

/// <summary>supermatcher_align / needleall_align / ksw2_extend and the `format` parameter of the pairwise and multiple-alignment tools (audit A2-27b / A2-28c / A2-29b, B09 F124 - F126). The tools delegate to
/// EmbossPairwisePrograms.Supermatcher / NeedleAll, Ksw2Aligner and AlignmentFormats / EmbossAlignmentWriter; the oracle comparisons live in the library tests.</summary>
[TestFixture]
public class EmbossProgramToolTests
{

    [Test]
    public void FormatParameter_MultipleAlignmentTools_WriteTheCanonicalFormats()
    {
        var r = AlignmentTools.MultipleAlign(["ACGT", "ACGA", "ACGT"], format: "phylip");
        Assert.That(r.Formatted, Does.StartWith("3 4\nseq1      "));
        var s = AlignmentTools.MultipleAlign(["ACGT", "ACGA", "ACGT"], format: "stockholm", names: ["a", "b", "c"]);
        Assert.That(s.Formatted, Does.StartWith("# STOCKHOLM 1.0\n#=GF SQ   3\n"));
        var m = AlignmentTools.MultipleAlign(["ACGT", "ACGA", "ACGT"], format: "msf");
        Assert.That(m.Formatted, Does.Contain("MSF: 4"));
        Assert.That(AlignmentTools.MultipleAlign(["ACGT", "ACGA"]).Formatted, Is.Null);
    }}
