// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP server Seqeron.Mcp.OligoThermodynamics (Primer3-style primer / probe design, ntthal); licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the server project).
using NUnit.Framework;
using Seqeron.Mcp.OligoThermodynamics.Tools;

namespace Seqeron.Mcp.OligoThermodynamics.Tests;

[TestFixture]
public class PrimerDimerTests
{
    [Test]
    public void PrimerDimer_Schema_ValidatesCorrectly()
    {
        Assert.DoesNotThrow(() => OligoThermodynamicsTools.primer_dimer("AAAAAAAA", "AAAAAAAA"));
        Assert.Throws<ArgumentException>(() => OligoThermodynamicsTools.primer_dimer("", "AAAA"));
        Assert.Throws<ArgumentException>(() => OligoThermodynamicsTools.primer_dimer("AAAA", null!));
    }

    [Test]
    public void PrimerDimer_Binding_InvokesSuccessfully()
    {
        // Primer3 alignment-mode PRIMER_PAIR_COMPL_END (primer3-py 2.3.1 check_primers,
        // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0): primer2 ending CGATCGAT = revcomp(ATCGATCG) → 8.
        var strong = OligoThermodynamicsTools.primer_dimer("AACCGGTTAACCATCGATCG", "AACCGGTTAACGATCGAT");
        // TTCAGTCAGTCAGTGGCC + ACTGACTGACTGAGGCC: 3' GGCC/GGCC overlap → 4.
        var ggcc = OligoThermodynamicsTools.primer_dimer("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC", min_complementarity: 5);
        // Identical poly-A primers cannot pair → 0.
        var polyA = OligoThermodynamicsTools.primer_dimer("AAAAAAAA", "AAAAAAAA");
        Assert.Multiple(() =>
        {
            Assert.That(strong.ComplementaryBases, Is.EqualTo(8));
            Assert.That(strong.ComplEndScore, Is.EqualTo(8.0));
            // PRIMER_PAIR_0_COMPL_ANY (check_primers, alignment mode): 8.0 and 13.0; poly-A 0.
            Assert.That(strong.ComplAnyScore, Is.EqualTo(8.0));
            Assert.That(ggcc.ComplAnyScore, Is.EqualTo(13.0));
            Assert.That(polyA.ComplAnyScore, Is.EqualTo(0.0));
            Assert.That(strong.HasDimer, Is.True);
            Assert.That(ggcc.ComplementaryBases, Is.EqualTo(4));
            Assert.That(ggcc.HasDimer, Is.False, "4 < min 5");
            Assert.That(OligoThermodynamicsTools.primer_dimer("TTCAGTCAGTCAGTGGCC", "ACTGACTGACTGAGGCC").HasDimer, Is.True);
            Assert.That(polyA.ComplementaryBases, Is.EqualTo(0));
            Assert.That(polyA.HasDimer, Is.False);
        });
    }
}
