// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Emboss.SequenceFormats;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// A5-8c3h (B09 WP59): the nibble / nib readers of EMBOSS 6.6.0 (<c>seqReadNibble</c>) fail on every file (measured with the installed <c>seqret -sformat1 nibble</c> and <c>nib</c> on FASTA text and on a nibble file written by
/// <c>seqret -osformat2 nibble</c>, valgrind clean): "Unable to read sequence" and "Died: seqret terminated: Bad value for '-sequence' with -auto defined".
/// </summary>
[TestFixture]
public class EmbossNibbleReader_Tests
{
    private static readonly string[] Expected = ["Error: Unable to read sequence 'in.dat'", "Died: seqret terminated: Bad value for '-sequence' with -auto defined"];

    [TestCase("nibble", false)]
    [TestCase("nib", false)]
    [TestCase("nibble", true)]
    [TestCase("nib", true)]
    public void ExplicitNibbleFormat_FailsLikeTheLibrary(string format, bool nibbleFile)
    {
        // the file the library's writer makes for ten bases: the magic number 3aede96b, the length, 4-bit codes
        string text = nibbleFile ? new string([(char)0x3a, (char)0xed, (char)0xe9, (char)0x6b, (char)10, (char)0, (char)0, (char)0, (char)0x33, (char)0x11, (char)0x22, (char)0]) : ">s1\nACGT\n";
        var o = EmbossSeqretPrograms.Seqret(new EmbossSeqretRequest("in.dat", new Dictionary<string, string> { ["in.dat"] = text }, new EmbossSequenceQualifiers { Format = format }));
        Assert.That(o.Died, Is.True);
        Assert.That(o.Messages, Is.EqualTo(Expected));
    }
}
