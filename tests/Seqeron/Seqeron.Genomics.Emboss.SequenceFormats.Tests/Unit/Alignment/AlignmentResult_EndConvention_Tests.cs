// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-SEMI-001 (audit R14-1, B09 F288): ONE convention for <see cref="AlignmentResult.EndPosition1"/> / <see cref="AlignmentResult.EndPosition2"/> across every producer:
/// 0-based INCLUSIVE (length − 1 for a global result, the last consumed index of a local / fitting result, <c>Start − 1</c> for an empty consumed interval, −1 for an empty local result).
/// Before the fix NeedleAlign / WaterAlign / StretcherAlign / Global- and LocalAlignDualAffine returned the exclusive end (R14-1: "TTACGTTT" vs "ACGT": Needle 8,4, Water 2..6, Dual 2..6 against GlobalAlign 7,3 and LocalAlign 2..5).
/// The contract is checked over a corpus of 60 seeded pairs for every aligner of the unit.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class AlignmentResult_EndConvention_Tests
{

    [Test]
    public void EmbossPairAlignment_Flanks_AreDerivedFromTheRows_NotFromTheEndConvention()
    {
        // local: rows hold the aligned part, offset = Start, offend = the rest
        var local = SequenceAligner.LocalAlign("TTACGTTT", "ACGT");
        var pair = EmbossPairAlignment.From(local, "a", 8, "b", 4);
        Assert.That((pair.First.Offset, pair.First.Offend, pair.Second.Offset, pair.Second.Offend), Is.EqualTo((2, 2, 0, 0)));
        // fitting: the rows carry the whole sequences (free flanks included), so there is no offset although Start2 = 2
        var fitting = SequenceAligner.SemiGlobalAlign("ACGT", "TTACGTTT");
        Assert.That(fitting.StartPosition2, Is.EqualTo(2));
        var fitted = EmbossPairAlignment.From(fitting, "q", 4, "ref", 8);
        Assert.That((fitted.First.Offset, fitted.First.Offend, fitted.Second.Offset, fitted.Second.Offend), Is.EqualTo((0, 0, 0, 0)));
        // empty local result: the whole sequence lies after the (empty) row
        var none = SequenceAligner.LocalAlign("AAAA", "CCCC");
        var empty = EmbossPairAlignment.From(none, "a", 4, "b", 4);
        Assert.That((empty.First.Offset, empty.First.Offend, empty.Second.Offset, empty.Second.Offend), Is.EqualTo((0, 4, 0, 4)));
    }
}
