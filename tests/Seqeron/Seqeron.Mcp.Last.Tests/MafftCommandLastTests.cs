// SPDX-License-Identifier: GPL-3.0-or-later
// Wrapper of the LAST 1256 port (GPL-3.0-or-later).
// Part of Seqeron.Mcp.Last, licensed under the GNU General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Last.Tools;

namespace Seqeron.Mcp.Last.Tests;

/// <summary>The tool <c>mafft_command_last</c>: the script mafft with the port of LAST behind it.</summary>
[TestFixture]
public sealed class MafftCommandLastTests
{
    private const string Two = ">a\nACGTTGCATGCATTAGCATCGATCGGATCGATTGCATCGATCGGATTAGCATCGAGGATCGATTAGGCTAGCTAGCTAGGATCGATC\n>b\nACGTTGCATGCATTAGCATCGATCGGATCGATTGCATCGATCGGATTAGCATCGAGGATCGATTAGGCTAGCTAGCTAGGATCGATC\n";

    [Test]
    public void LastPair_TwoSequences_IsAnswered()
    {
        var r = LastTools.MafftCommandLast(["--lastpair", "in.fa"], new Dictionary<string, string> { ["in.fa"] = Two });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.Zero);
            Assert.That(r.Stdout, Does.StartWith(">a\n"));
            Assert.That(r.Stderr, Does.Contain("Pairwise alignments were computed by LAST"));
        });
    }

    [Test]
    public void LastMultiPair_WithoutFragments_StopsWithTheMessageOfTheScript()
    {
        var r = LastTools.MafftCommandLast(["--lastmultipair", "in.fa"], new Dictionary<string, string> { ["in.fa"] = Two });
        Assert.Multiple(() =>
        {
            Assert.That(r.ExitCode, Is.EqualTo(1));
            Assert.That(r.Stderr, Is.EqualTo("Specify --addfragments too\n"));
        });
    }

    [Test]
    public void NoArguments_AreRefused()
    {
        Assert.Throws<ArgumentException>(() => LastTools.MafftCommandLast([]));
    }
}
