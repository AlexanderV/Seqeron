// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.ClustalW;
using Seqeron.Genomics.TCoffee;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (audit A2-16, B09 F93): the paper path of <see cref="Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(IEnumerable{DnaSequence}, ScoringMatrix?, ConsistencyMsaOptions?)"/> now has the two options the T-Coffee
/// program path had since F58 / F66: the top-N non-intersecting local alignments of Lalign in the primary library (<see cref="ConsistencyMsaOptions.PrimaryLocalAlignments"/>, built from
/// <see cref="SequenceAligner.LocalAlignTopN(string, string, ScoringMatrix?, int, int)"/>, itself oracle-verified against EMBOSS matcher: F47) and the redundancy weights
/// (<see cref="ConsistencyMsaOptions.SequenceWeights"/> = <c>-seq_weight t_coffee</c>, the weight function of the program path, oracle-verified in F66, applied after the extension as
/// <c>re_weight_constraint_list</c> does). The combination in the paper-formula engine cannot be compared with a t_coffee run (the program's library is SIM + integer weights + relaxation, a different pipeline).
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_ConsistencyOptions_Tests
{
    // two sequences that share TWO separated motifs (so the second local alignment is not the first one), plus a distinct third
    private static readonly string[] Seqs = ["GGCATTACGTCCCCCCCCGATCGGATTACAGG", "TTCATTACGTAAAAAAAAAGATCGGATTACATT", "CCGATCGGATTCACGTAAAGTCAT"];

    private static double PercentId(string a, string b)
    {
        int aligned = 0, identical = 0;
        for (int c = 0; c < a.Length; c++)
        {
            if (a[c] == '-' || b[c] == '-')
                continue;
            aligned++;
            if (a[c] == b[c])
                identical++;
        }

        return aligned == 0 ? 0 : (int)(100.0 * identical) / aligned;
    }

    [Test]
    public void Defaults_AreTheSingleBestLocalAlignment_AndNoWeights()
    {
        var o = new ConsistencyMsaOptions();
        Assert.Multiple(() =>
        {
            Assert.That(o.PrimaryLocalAlignments, Is.EqualTo(1));
            Assert.That(o.SequenceWeights, Is.EqualTo(MsaSequenceWeights.None));
            Assert.That(Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(Seqs.Select(s => new DnaSequence(s)), null, o).AlignedSequences,
                Is.EqualTo(Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(Seqs.Select(s => new DnaSequence(s))).AlignedSequences));
        });
    }

    [Test]
    public void TopN_PrimaryLibrary_IsTheGlobalLibraryPlusEveryTopNLocalAlignmentWeightedByItsOwnIdentity()
    {
        var o = new ConsistencyMsaOptions { PrimaryLocalAlignments = 3 };
        var (lib, _) = Seqeron.Genomics.TCoffee.TCoffeeAligner.GetPrimaryLibrary(Seqs, null, o);
        var one = Seqeron.Genomics.TCoffee.TCoffeeAligner.GetPrimaryLibrary(Seqs, null, new ConsistencyMsaOptions()).Primary;

        // expected: for every pair the NW alignment (weight = its identity) + each of the top-3 Waterman-Eggert alignments (each its own identity), signal-added
        var expected = new Dictionary<(int, int, int, int), double>();
        void Add(int i, int pi, int j, int pj, double w) => expected[(i, pi, j, pj)] = expected.GetValueOrDefault((i, pi, j, pj)) + w;
        void AddAlignment(int i, int j, AlignmentResult a, int start1, int start2)
        {
            double w = PercentId(a.AlignedSequence1, a.AlignedSequence2);
            if (w <= 0)
                return;
            int p1 = start1, p2 = start2;
            for (int c = 0; c < a.AlignedSequence1.Length; c++)
            {
                char x = a.AlignedSequence1[c], y = a.AlignedSequence2[c];
                if (x != '-' && y != '-')
                {
                    Add(i, p1, j, p2, w);
                    p1++;
                    p2++;
                }
                else if (x != '-')
                    p1++;
                else if (y != '-')
                    p2++;
            }
        }

        int topAlignmentsUsed = 0;
        for (int i = 0; i < Seqs.Length; i++)
            for (int j = i + 1; j < Seqs.Length; j++)
            {
                AddAlignment(i, j, SequenceAligner.GlobalAlign(Seqs[i], Seqs[j]), 0, 0);
                foreach (var la in SequenceAligner.LocalAlignTopN(Seqs[i], Seqs[j], null, 3))
                {
                    AddAlignment(i, j, la, la.StartPosition1, la.StartPosition2);
                    topAlignmentsUsed++;
                }
            }

        Assert.Multiple(() =>
        {
            Assert.That(topAlignmentsUsed, Is.GreaterThan(3), "more than one local alignment per pair exists for these sequences");
            Assert.That(lib, Is.EquivalentTo(expected));
            Assert.That(lib.Count, Is.GreaterThan(one.Count), "the extra alignments add residue pairs");
            Assert.That(lib.Values.Sum(), Is.GreaterThan(one.Values.Sum()));
        });
    }

    [Test]
    public void TopN_One_IsExactlyTheOldSingleBestLibrary_AndInvalidValuesAreRejected()
    {
        var a = Seqeron.Genomics.TCoffee.TCoffeeAligner.GetPrimaryLibrary(Seqs, null, new ConsistencyMsaOptions { PrimaryLocalAlignments = 1 }).Primary;
        var b = Seqeron.Genomics.TCoffee.TCoffeeAligner.GetPrimaryLibrary(Seqs, null, new ConsistencyMsaOptions()).Primary;
        Assert.That(a, Is.EquivalentTo(b));
        Assert.Throws<ArgumentOutOfRangeException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(Seqs, MatrixScoring.Emboss(SubstitutionMatrix.EdnaFull, 10, 0.5), new ConsistencyMsaOptions { PrimaryLocalAlignments = 0 }));
        Assert.Throws<ArgumentOutOfRangeException>(() => Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(Seqs, MatrixScoring.Emboss(SubstitutionMatrix.EdnaFull, 10, 0.5), new ConsistencyMsaOptions { SequenceWeights = (MsaSequenceWeights)7 }));
    }

    [Test]
    public void SequenceWeights_TCoffee_ScaleTheExtendedWeightByTheSmallerSequenceWeight()
    {
        // a redundant pair of identical sequences and one distinct sequence: the duplicates share the weight mass, the distinct sequence weighs more; mean weight 1
        string[] seqs = [Seqs[0], Seqs[0], Seqs[2]];
        var (_, w) = Seqeron.Genomics.TCoffee.TCoffeeAligner.GetPrimaryLibrary(seqs, null, new ConsistencyMsaOptions { SequenceWeights = MsaSequenceWeights.TCoffee });
        Assert.That(w, Is.Not.Null);
        Assert.Multiple(() =>
        {
            Assert.That(w!.Average(), Is.EqualTo(1.0).Within(1e-5), "compute_t_coffee_weight rescales to a mean of 1");
            Assert.That(w[2], Is.GreaterThan(w[0]), "the distinct sequence outweighs the duplicated one");
            Assert.That(w[0], Is.EqualTo(w[1]).Within(1e-6));
        });

        var plain = new ConsistencyMsaOptions();
        var weighted = new ConsistencyMsaOptions { SequenceWeights = MsaSequenceWeights.TCoffee };
        foreach (var (si, pi, sj, pj) in new[] { (0, 3, 1, 3), (0, 5, 2, 5), (1, 8, 2, 7) })
        {
            double e0 = Seqeron.Genomics.TCoffee.TCoffeeAligner.GetExtendedWeight(seqs, si, pi, sj, pj, null, plain);
            double e1 = Seqeron.Genomics.TCoffee.TCoffeeAligner.GetExtendedWeight(seqs, si, pi, sj, pj, null, weighted);
            Assert.That(e1, Is.EqualTo(e0 * Math.Min(w![si], w[sj])).Within(1e-9), $"({si}.{pi}, {sj}.{pj})");
        }

        // and the whole pipeline accepts the option for both engines
        foreach (var engine in new[] { MsaConsistencyEngine.PaperExtension, MsaConsistencyEngine.TCoffeeLinked })
        {
            var r = Seqeron.Genomics.TCoffee.TCoffeeAligner.MultipleAlignConsistency(seqs.Select(s => new DnaSequence(s)), null, new ConsistencyMsaOptions { SequenceWeights = MsaSequenceWeights.TCoffee, PrimaryLocalAlignments = 3, ConsistencyEngine = engine });
            Assert.That(r.AlignedSequences.Select(x => x.Replace("-", "")), Is.EqualTo(seqs), engine.ToString());
        }
    }
}
