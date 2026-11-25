// SPDX-License-Identifier: GPL-3.0-or-later
// Wrappers of the copyleft multiple-alignment assemblies (Kalign 3.3.5, GPL-3.0-or-later and others).
// Part of Seqeron.Mcp.MultipleAlignment, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.PartialOrderAlignment;
using Seqeron.Genomics.Kalign;
using Seqeron.Genomics.TCoffee;
using Seqeron.Genomics.ClustalW;
using NUnit.Framework;
using Seqeron.Mcp.MultipleAlignment.Tools;

namespace Seqeron.Mcp.MultipleAlignment.Tests;

/// <summary>The Kalign share of <c>Seqeron.Mcp.Alignment.Tests.NullContractTests</c> (moved with the tool, B09 WP79).</summary>
[TestFixture]
public class KalignNullContractTests
{
    [Test]
    public void KalignCommand_NullArgumentElement_IsAnArgumentNullException()
    {
        string[] args = ["-i", null!];
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => MultipleAlignmentTools.KalignCommand(args))!.ParamName, Is.EqualTo("arguments"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.KalignCommand(null!));
        });
    }
}
