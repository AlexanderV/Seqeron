// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of Primer3 ntthal (thal.c) and of the Primer3-style primer / probe design; part of Seqeron.Genomics.OligoThermodynamics.Tests, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the production project).
using FsCheck;
using FsCheck.Fluent;

namespace Seqeron.Genomics.Tests.Algebraic;

/// <summary>
/// Algebraic-law tests for the MolTools area (CRISPR off-target scoring,
/// restriction digestion).
///
/// Algebraic testing pins the formal equations the scoring/transform functions
/// must obey: the identity value on a perfect/empty input and the monotonic /
/// conservation structure as the input degrades or is partitioned.
/// — docs/checklists/06_ALGEBRAIC_TESTING.md §Description, rows 20, 27.
/// </summary>
[TestFixture]
[Category("Algebraic")]
[Category("MolTools")]
public class MolToolsAlgebraicTests
{

    [Test]
    public void Dimer_Idempotent_Deterministic()
    {
        double a = ThermodynamicPrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", DimerNa, DimerCt);
        double b = ThermodynamicPrimerDesigner.CalculateSelfDimerMeltingTemperature("GCGCGCGC", DimerNa, DimerCt);
        b.Should().Be(a);
    }


    // ═══════════════════════════════════════════════════════════════════════
    // Unit: PRIMER-DIMER-001 — ntthal self/hetero-dimer (MolTools), row 242.
    //
    // Model: the most stable duplex between two strands is the contiguous Watson–Crick
    //        alignment whose ΔH/ΔS is the nearest-neighbour SUM over its consecutive
    //        base-pair steps (ntthal thermodynamic alignment).
    //   — PrimerDesigner.FindMostStableDimer / CalculateSelfDimerMeltingTemperature.
    //     TestSpec tests/TestSpecs/PRIMER-DIMER-001.md.
    //
    // Laws (row 242): ID — the contiguous-WC optimum equals the NN sum.
    //                 IDEMP — dimer search/Tm is a pure, deterministic function.
    // ═══════════════════════════════════════════════════════════════════════

    private const double DimerNa = 0.05;    // 50 mM monovalent
    private const double DimerCt = 50e-9;
}
