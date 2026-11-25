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
/// Audit A5-3 (B09 F169): the MCP tools over the secondary entry points refuse a null element / a nameless sequence / a matrix beyond the char-indexed table with a typed
/// <see cref="ArgumentException"/> (the shared message style of the library), not a <see cref="NullReferenceException"/> / <see cref="IndexOutOfRangeException"/>.
/// </summary>
[TestFixture]
public class NullContractTests
{
    [Test]
    public void CommandTools_NullArgumentElement_IsAnArgumentNullException()
    {
        string[] args = ["-i", null!];
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => MultipleAlignmentTools.ProbConsCommand(args))!.ParamName, Is.EqualTo("arguments"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.ProbConsCommand(null!));
        });
    }

    [Test]
    public void TCoffeeStageTools_NullElementAndForeignCharacters_AreArgumentExceptions()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => MultipleAlignmentTools.TCoffeeKtupSimilarity(["ACGTACGT", null!]))!.ParamName, Is.EqualTo("sequences"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeLalignTopAlignments("ACGT*ACGT", "ACGTACGT"));
            Assert.Throws<ArgumentException>(() => MultipleAlignmentTools.TCoffeeLalignTopAlignments("ACGTACGT", "ACGéACGT"));
        });
    }

    [Test]
    public void MultipleAlignTools_NullSequenceList_AndNullStructureOutput_AreArgumentExceptions()
    {
        // found by the reflection fuzzer of audit A5-3: ParseDna(null) and a null element of structureAlignerOutputs were NullReferenceExceptions
        Assert.Multiple(() =>
        {
            Assert.That(Assert.Throws<ArgumentNullException>(() => MultipleAlignmentTools.MultipleAlignClustalW(null!))!.ParamName, Is.EqualTo("sequences"));
            Assert.Throws<ArgumentNullException>(() => MultipleAlignmentTools.MultipleAlignTCoffee(["ACGTACGT", "ACGTACGA"], structureAlignerOutputs: [null!]));
        });
    }
}
