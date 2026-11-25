// SPDX-License-Identifier: GPL-3.0-only
// Derived from DECIPHER 3.9.4 (Bioconductor), ConsensusSequence, Copyright (C) Erik Wright; GPL-3.0-only.
// Part of Seqeron.Genomics.AlignmentConsensus, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Tests.Fuzzing;

using Seqeron.Genomics.AlignmentConsensus;

/// <summary>
/// Heavy-tier fuzz tests for the review 2026-09 B05 audit additions (docs/Validation/review-2026-09/B05.md F27–F35):
/// weighted edit / Damerau distances and alignments, PWM p-values (DNA, Markov, alphabet) with random options,
/// generic-alphabet PWMs, DECIPHER / EMBOSS Auto consensus, RSAT -seqtype oligo analysis, dyad analysis and σ70
/// promoter pairing / Promoter Calculator. Random inputs (fixed seeds) must either produce a well-formed result or
/// throw one of the documented exceptions (ArgumentException and subtypes; OverflowException for a weighted distance
/// beyond int.MaxValue) — never anything else. σ70 pairing is also checked against a brute-force enumeration and its
/// reverse-complement relation.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
[Category("Matching")]
public class B05AuditFuzzTests
{
    private static string RandomString(Random rng, int length, string pool) =>
        new(Enumerable.Range(0, length).Select(_ => pool[rng.Next(pool.Length)]).ToArray());

    private static void Documented(Action call, bool allowOverflow = false)
    {
        try
        {
            call();
        }
        catch (ArgumentException)
        {
        }
        catch (OverflowException) when (allowOverflow)
        {
        }
    }

    [Test]
    public void DecipherAndEmbossAuto_RandomInput_OnlyDocumentedExceptions()
    {
        var rng = new Random(350_104);
        const string pool = "ACGTUNRYKMSWBDHVacgtu-.+~*?XEFILPQZJO 1é";
        for (int trial = 0; trial < 1000; trial++)
        {
            var rows = Enumerable.Range(0, rng.Next(0, 8)).Select(_ => RandomString(rng, rng.Next(0, 20), pool)).ToArray();
            int maxLen = rows.Length == 0 ? 0 : rows.Max(r => r.Length);
            Documented(() =>
            {
                string c = DecipherConsensus.GenerateDecipherConsensus(rows, (DecipherSequenceType)rng.Next(0, 4), rng.NextDouble() * 1.2 - 0.1,
                    rng.Next(2) == 0, "+-NX#"[rng.Next(5)], rng.Next(3) == 0 ? null : rng.NextDouble() * 1.2, rng.Next(2) == 0, rng.Next(2) == 0);
                Assert.That(c.Length, Is.EqualTo(maxLen));
            });
            Documented(() =>
            {
                string c = MotifFinder.GenerateEmbossConsensus(rows, true, ConsensusResidueType.Auto, rng.Next(3) == 0 ? null : rng.Next(-1, 5));
                // the rows are read as EMBOSS reads an alignment row (A4-3): digits, blanks and other characters are dropped before the padding
                Assert.That(c.Length, Is.EqualTo(rows.Length == 0 ? 0 : rows.Max(r => Seqeron.Genomics.Infrastructure.EmbossSequenceReader.CleanGapped(r).Length)));
            });
        }
    }}
