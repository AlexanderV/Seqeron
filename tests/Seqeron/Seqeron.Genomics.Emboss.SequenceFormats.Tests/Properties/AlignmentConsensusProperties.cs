// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.SequenceFormats, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property tests for the B05 follow-up consensus algorithms (MOTIF-CONS-001 / MOTIF-GENERATE-001):
///   C1 (default threshold)  GenerateConsensus(rows, 0.25) is bit-identical to GenerateConsensus(rows).
///   C2 (threshold nesting)  for t1 ≤ t2 the base set of GenerateConsensus(rows, t2) is a subset of the
///                           base set at t1 (column by column); every output symbol is an IUPAC code.
///   C3 (EMBOSS shape)       GenerateEmbossConsensus returns one symbol per column; every symbol is the
///                           (case-folded) residue of some row in that column or the no-consensus symbol;
///                           upper case ⇔ positive matches above setcase, so setcase = −1 gives all upper case.
///   C4 (dumb shape)         GenerateDumbConsensus symbols are either the ambiguous symbol or the strict
///                           majority residue of the column with fraction ≥ threshold.
/// Seeded random inputs, deterministic.
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Motif")]
public class AlignmentConsensusProperties
{

    [Test]
    public void C3_Emboss_SymbolsComeFromColumnOrNoConsensus()
    {
        var rng = new Random(5503);
        for (int trial = 0; trial < 300; trial++)
        {
            bool protein = trial % 2 == 1;
            string alphabet = protein ? "ACDEFGHIKLMNPQRSTVWYX-" : "ACGTRYN-";
            var rows = RandomRows(rng, rng.Next(2, 10), rng.Next(0, 25), alphabet);
            var type = protein ? ConsensusResidueType.Protein : ConsensusResidueType.Nucleotide;
            char nocon = protein ? 'X' : 'N';
            string cons = EmbossConsensus.GenerateEmbossConsensus(rows, type, plurality: (float)rng.Next(0, 6), identity: rng.Next(0, 3));
            Assert.That(cons.Length, Is.EqualTo(rows[0].Length));
            for (int k = 0; k < cons.Length; k++)
            {
                char u = char.ToUpperInvariant(cons[k]);
                Assert.That(u == nocon || rows.Any(r => char.ToUpperInvariant(r[k]) == u), Is.True, $"col {k} '{cons[k]}'");
            }

            string upper = EmbossConsensus.GenerateEmbossConsensus(rows, type, plurality: 0f, setcase: -1f);
            Assert.That(upper, Is.EqualTo(upper.ToUpperInvariant()));
        }
    }

    private static string[] RandomRows(Random rng, int n, int length, string alphabet) =>
        Enumerable.Range(0, n).Select(_ => new string(Enumerable.Range(0, length)
            .Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray())).ToArray();
}
