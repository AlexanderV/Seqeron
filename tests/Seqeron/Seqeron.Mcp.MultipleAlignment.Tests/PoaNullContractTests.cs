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

/// <summary>
/// Audit A5-3 (B09 F169): the MCP tools over the secondary entry points refuse a null element / a nameless sequence / a matrix beyond the char-indexed table with a typed
/// <see cref="ArgumentException"/> (the shared message style of the library), not a <see cref="NullReferenceException"/> / <see cref="IndexOutOfRangeException"/>.
/// </summary>
[TestFixture]
public class PoaNullContractTests
{

    [Test]
    public void MultipleAlignPoa_NullName_AndWideMatrix_AreArgumentExceptions()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => MultipleAlignmentTools.MultipleAlignPoa(["ACGT", "ACGA"], names: ["a", null!]))!.Message, Does.Contain("Sequence 1 has no name"));
            Assert.That(Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignPoa(["ACGT", "ACGA"], matrix: new string('A', 129) + "\n"))!.Message, Does.Contain("at most 128"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.MultipleAlignPoa(["ACGT", "ACGA"], matrix: "Aé\n"));
        });
    }
}
