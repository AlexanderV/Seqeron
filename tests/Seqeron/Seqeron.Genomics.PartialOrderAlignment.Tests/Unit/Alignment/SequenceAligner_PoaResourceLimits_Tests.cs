// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from POA 2.0 (Partial Order Alignment, https://sourceforge.net/projects/poamsa/), Copyright (C) 2004-2006 Christopher Lee; GPL-2.0-or-later.
// Part of Seqeron.Genomics.PartialOrderAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.PartialOrderAlignment;
using NUnit.Framework;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>The POA share of <c>SequenceAligner_ResourceLimits_Tests</c> (moved with the port to this assembly, B09 WP79).</summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_PoaResourceLimits_Tests
{
    private const long MaxAllocatedBytes = 64L << 20;

    private static T Refused<T>(Action call) where T : Exception
    {
        long before = GC.GetAllocatedBytesForCurrentThread();
        var ex = Assert.Throws<T>(() => call());
        long allocated = GC.GetAllocatedBytesForCurrentThread() - before;
        Assert.That(allocated, Is.LessThan(MaxAllocatedBytes), "the limit is checked before the matrices are allocated");
        return ex!;
    }

    [Test]
    public void Poa_Above2Pow27Cells_IsRefusedWithoutAllocating()
    {
        string a = new('A', 12_000), b = new('C', 12_000);
        var ex = Refused<ArgumentException>(() => PoaAligner.MultipleAlignPoa([a, b]));
        Assert.That(ex.Message, Does.StartWith("POA keeps").And.Contain("PoaAligner.PoaMaxMatrixCells").And.Contain("2^31"));
    }

    [Test]
    public void Poa_LimitMessage_IsCultureInvariant()
    {
        var odd = new System.Globalization.CultureInfo("en-US");
        odd.NumberFormat.NumberGroupSeparator = ".";
        odd.NumberFormat.NumberDecimalSeparator = ",";
        var saved = System.Globalization.CultureInfo.CurrentCulture;
        string poa;
        try
        {
            System.Globalization.CultureInfo.CurrentCulture = odd;
            string pa = new('A', 12_000), pb = new('C', 12_000);
            poa = Refused<ArgumentException>(() => PoaAligner.MultipleAlignPoa([pa, pb])).Message;
        }
        finally { System.Globalization.CultureInfo.CurrentCulture = saved; }
        Assert.That(poa, Does.Contain("12,000 graph nodes").And.Not.Contain("12.000"));
    }
}
