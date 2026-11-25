// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from POA 2.0 (Partial Order Alignment, https://sourceforge.net/projects/poamsa/), Copyright (C) 2004-2006 Christopher Lee; GPL-2.0-or-later.
// Part of Seqeron.Genomics.PartialOrderAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.PartialOrderAlignment;
using System.Text.Json.Nodes;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 / ALIGN-STATS-001 (audit A5-3, B09 F169): the SECONDARY public entry points of the alignment classes refuse a null list / null element / out-of-range value with a typed
/// <see cref="ArgumentException"/> family exception naming the parameter (the library's contract; the main entry points were fixed in F81 / F82). Before: <see cref="NullReferenceException"/> /
/// <see cref="IndexOutOfRangeException"/> / <see cref="KeyNotFoundException"/> found by a reflection fuzzer over every public static method.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_NullContract_Tests
{

    // ------------------------------------------------------------------ (8) POA: a nameless sequence and a matrix beyond the char-indexed table

    private static void ThrowsNull(TestDelegate call, string paramName) =>
        Assert.That(Assert.Throws<ArgumentNullException>(call)!.ParamName, Is.EqualTo(paramName));

    [Test]
    public void EmptyList_Poa_RefusesLikeTheProgram()
    {
        var ex = Assert.Throws<ArgumentException>(() => PoaAligner.MultipleAlignPoa(Array.Empty<string>()));
        Assert.That(ex!.Message, Does.StartWith("POA cannot align an empty list of sequences: the program refuses an empty input ("));
        Assert.That(ex.Message, Does.Contain("MultipleAlignMafft and the Seqeron aligners").And.Contain("return an empty result."));
    }

    [Test]
    public void PoaRunner_NullArgumentOrNullFileText_IsAnArgumentNullException()
    {
        string?[] args = ["-i", null];
        var files = new Dictionary<string, string> { ["in.fa"] = null! };
        Assert.Multiple(() =>
        {
            ThrowsNull(() => PoaProgram.Run(args!), "arguments");
            ThrowsNull(() => PoaProgram.Run(["-read_fasta", "in.fa"], files), "files");
            // the program's own argument handling is untouched: an empty list is the usage screen
            Assert.That(PoaProgram.Run([]).ExitCode, Is.EqualTo(255));
        });
    }

    [Test]
    public void Poa_NamelessSequence_AndMatrixBeyondTheCharIndexedTable_AreArgumentExceptions()
    {
        ThrowsNull(() => PoaAligner.PoaRun([new PoaSequence(null!, null, "ACGT")]), "sequences");
        ThrowsNull(() => PoaAligner.PoaRun([new PoaSequence("a", null, null!)]), "sequences");
        ThrowsNull(() => PoaAligner.PoaRun([null!]), "sequences");
        var tooMany = Assert.Throws<ArgumentException>(() => PoaScoreMatrix.Parse(new string('A', 129) + "\n"))!;
        Assert.That(tooMany.Message, Does.Contain("129 symbols").And.Contain("at most 128"));
        var wide = Assert.Throws<ArgumentException>(() => PoaScoreMatrix.Parse("Aé\n"))!;
        Assert.That(wide.Message, Does.Contain("U+00E9"));
        Assert.That(PoaScoreMatrix.Parse(new string('A', 128) + "\n").Symbols, Has.Length.EqualTo(128));
    }
}
