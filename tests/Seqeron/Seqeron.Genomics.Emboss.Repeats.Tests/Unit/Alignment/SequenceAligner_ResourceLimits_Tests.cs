// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 / ALIGN-MULTI-001 (audit A3-2, A3-3; B09 F131, F132): the resource limits of the program ports (one <see cref="ArgumentException"/> style, 64-bit arithmetic, raised
/// BEFORE anything large is allocated: every limit test asserts the exception AND that the call allocated less than 64 MB) and the NaN / Infinity / overflowing penalty contract of the EMBOSS ports.
/// The programs' own limit is "none": installed needle / water die with "Sequences too big. Try 'stretcher'" / "Try 'matcher' or 'supermatcher'" when the allocation fails (measured under a 4 GB
/// address-space limit with 30 000 x 30 000 residues), installed stretcher wraps its 32-bit score and warns "stretcher CheckScore failed" (1e9 penalties: score 1 949 672 837).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ResourceLimits_Tests
{

    [Test]
    public void WordMatches_AboveTheMatchLimit_AreRefused_BySeqMatchAll()
    {
        // wordsize 1 over two random 5 200-base sequences: about 0.1875 x 27 million = 5.1 million maximal matches (> 2^22 = 4 194 304); the engine behind seqmatchall and megamerger is shared (PairWordMatcher), every later list is bounded
        string a = Dna(5_200, 31), b = Dna(5_200, 32);
        long before = GC.GetAllocatedBytesForCurrentThread();
        var ex = Assert.Throws<ArgumentException>(() => EmbossRepeatPrograms.SeqMatchAll([a, b], 1))!;
        Assert.That(GC.GetAllocatedBytesForCurrentThread() - before, Is.LessThan(1L << 30), "only the engine's own list is built");
        Assert.That(ex.Message, Does.StartWith("EMBOSS seqmatchall / megamerger keep every maximal word match").And.Contain("SequenceAligner.WordMatchMaxMatches"));
        Assert.That(SequenceAligner.WordMatchMaxMatches, Is.EqualTo(4_194_304L));
        Assert.That(EmbossRepeatPrograms.SeqMatchAll([Dna(300, 33), Dna(300, 34)], 4), Has.Count.EqualTo(282), "an ordinary pair is untouched: the 282 match lines of installed seqmatchall 6.6.0 -wordsize 4 on these sequences");
    }

    private static string Dna(int n, int seed)
    {
        var r = new Random(seed);
        var c = new char[n];
        for (int i = 0; i < n; i++) c[i] = "ACGT"[r.Next(4)];
        return new string(c);
    }
}
