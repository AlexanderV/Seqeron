// SPDX-License-Identifier: GPL-3.0-or-later
// Derived from Kalign 3.3.5 (https://github.com/TimoLassmann/kalign), Copyright (C) 2006, 2019, 2020 Timo Lassmann; GPL-3.0-or-later.
// Part of Seqeron.Genomics.Kalign, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Kalign;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The Kalign share of the program-port contracts of <c>SequenceAligner_InputContract_Tests</c> and <c>SequenceAligner_NullContract_Tests</c> (moved with the port to this assembly, B09 WP79).</summary>
[TestFixture]
public class SequenceAligner_KalignContract_Tests
{
    private static void ThrowsNull(TestDelegate call, string paramName) =>
        Assert.That(Assert.Throws<ArgumentNullException>(call)!.ParamName, Is.EqualTo(paramName));

    [Test]
    public void EmptyList_Kalign_RefusesLikeTheProgram()
    {
        var ex = Assert.Throws<ArgumentException>(() => KalignAligner.MultipleAlignKalign(Array.Empty<string>()));
        Assert.That(ex!.Message, Does.StartWith("Kalign cannot align an empty list of sequences: the program refuses an empty input ("));
        Assert.That(ex.Message, Does.Contain("MultipleAlignMafft and the Seqeron aligners").And.Contain("return an empty result."));
    }

    [Test]
    public void KalignRunner_NullArgumentOrNullFileText_IsAnArgumentNullException()
    {
        string?[] args = ["-i", null];
        var files = new Dictionary<string, string> { ["in.fa"] = null! };
        Assert.Multiple(() =>
        {
            ThrowsNull(() => KalignProgram.Parse(args!), "arguments");
            ThrowsNull(() => KalignProgram.Run(args!), "arguments");
            ThrowsNull(() => KalignProgram.Run(["-i", "in.fa"], files), "files");
            Assert.That(KalignProgram.Run(["--no-such-option"]).ExitCode, Is.EqualTo(1));
        });
    }
}
