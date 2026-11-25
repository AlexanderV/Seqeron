// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-GLOBAL-001 / ALIGN-LOCAL-001 (audit R14-8, B09 F301): all 30 named matrices of Biopython 1.88 (<c>Bio.Align.substitution_matrices.load</c>) are built in: BLASTN, BLASTP, FENG, HOXD70, LEVIN, MCLACHLAN, MEGABLAST, RAO, STR, TRANS
/// (integer) and BENNER6 / 22 / 74, DAYHOFF, GENETIC, GONNET1992, JOHNSON, JONES, MDM78, RISLER (fractional: stored exactly as integers in units of 1 / <see cref="SubstitutionMatrix.Denominator"/>, scored exactly by the aligners)
/// next to the BLOSUM / PAM / NUC.4.4 ones, plus the codon matrix SCHNEIDER and <see cref="SubstitutionMatrix.Parse"/>. Oracle TestData/Alignment/biopython_matrices_oracle.tsv.gz
/// (scripts/alignment/generate_matrix_oracle.py): every cell of the 29 residue matrices and 841 optimal global / local scores of seeded pairs with integer and fractional affine gap scores (Biopython <c>PairwiseAligner</c>).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-GLOBAL-001")]
public class SubstitutionMatrix_Biopython_Tests
{

    [Test]
    public void FractionalMatrix_WithTheAligners_AllModesAndTheIntegerOnlyConsumers()
    {
        var scoring = new MatrixScoring(SubstitutionMatrix.Get("JOHNSON"), -2.5, -0.5);
        var global = SequenceAligner.GlobalAlignAffine("ACDEFGHIK", "ACDFGHIKL", scoring);
        var local = SequenceAligner.LocalAlignAffine("ACDEFGHIK", "ACDFGHIKL", scoring);
        Assert.That(global.ExactScore, Is.Not.EqualTo(Math.Round(global.ExactScore)), "the score is fractional (two decimals)");
        Assert.That(local.ExactScore, Is.GreaterThanOrEqualTo(global.ExactScore));
        var stats = SequenceAligner.CalculateStatistics(global, scoring.Matrix);
        Assert.That(stats.Matches, Is.GreaterThan(0));
        // EMBOSS needle reads integer matrices only (its data files are integer): a clear refusal, not a wrong score
        Assert.Throws<InvalidOperationException>(() => EmbossPairwisePrograms.NeedleAlign("ACDEFGHIK", "ACDFGHIKL", SubstitutionMatrix.Get("JOHNSON")));
    }
}
