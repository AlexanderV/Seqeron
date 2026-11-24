using System.Text;
using NUnit.Framework;
using Seqeron.Genomics.Alignment;
using Seqeron.Genomics.Core;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// ALIGN-MULTI-001 (B09 WP78, R14-16, F362): the speed-ups of the guide-tree multiple aligners (<c>MultipleAlignProgressive</c> / <c>MultipleAlignIterative</c> / <c>MultipleAlignConsistency</c>).
/// <list type="bullet">
/// <item>The exact ones change no result: the identity count of the distance of a pair from the forward recurrence alone (no matrix, no traceback strings), the profile-profile column scores from the compositions of the columns, and the
/// consistency column scores from the support of the extended library are compared with the alignments written by the code before them (commit 12059b65a, <c>msa_baseline_wp78.tsv.gz</c>: 117 families, DNA with a scalar
/// matrix, proteins with BLOSUM62, linear and affine gaps, UPGMA and neighbour-joining trees, progressive / iterative / consistency).</item>
/// <item>The k-tuple estimate (<see cref="MsaOptions.FastDistanceCells"/>) is used only above the work limit; its distance is on the identity scale of the alignment distance (d = 1 − F^(1/k)).</item>
/// </list>
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MsaGuideDistance_Tests
{
    private sealed record Case(string Kind, string Options, string[] Sequences, string[] Rows, string Score);

    private static readonly Lazy<List<Case>> Baseline = new(() =>
    {
        using var stream = TestDataResource.Open("Seqeron.Genomics.Tests.TestData.Alignment.msa_baseline_wp78.tsv") ?? throw new InvalidOperationException("baseline resource missing");
        using var reader = new StreamReader(stream, Encoding.UTF8);
        var list = new List<Case>();
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0)
                continue;
            var f = line.Split('\t');
            list.Add(new Case(f[0], f[1], f[2].Split('|'), f[3].Split('|'), f[4]));
        }

        return list;
    });

    private static MultipleAlignmentResult Run(Case c)
    {
        var options = new MsaOptions
        {
            GapModel = c.Options == "affine" ? MsaGapModel.Affine : MsaGapModel.Linear,
            GuideTree = c.Options == "nj" ? MsaGuideTree.NeighborJoining : MsaGuideTree.Upgma,
        };
        var protein = MatrixScoring.Emboss(SubstitutionMatrix.Blosum62);
        return c.Kind switch
        {
            "progressive" => SequenceAligner.MultipleAlignProgressive(c.Sequences.Select(s => new DnaSequence(s)).ToList(), null, options),
            "iterative" => SequenceAligner.MultipleAlignIterative(c.Sequences.Select(s => new DnaSequence(s)).ToList(), null, options),
            "consistency" => SequenceAligner.MultipleAlignConsistency(c.Sequences.Select(s => new DnaSequence(s)).ToList()),
            "progressive-protein" => SequenceAligner.MultipleAlignProgressive(c.Sequences, protein, options),
            "iterative-protein" => SequenceAligner.MultipleAlignIterative(c.Sequences, protein, options),
            _ => throw new InvalidOperationException(c.Kind),
        };
    }

    [Test]
    public void Baseline_HasEveryKind()
    {
        var kinds = Baseline.Value.GroupBy(c => c.Kind + "/" + c.Options).Select(g => g.Key).ToList();
        Assert.That(Baseline.Value, Has.Count.EqualTo(117));
        Assert.That(kinds, Is.SupersetOf(new[] { "progressive/linear", "progressive/affine", "progressive/nj", "iterative/linear", "iterative/affine", "consistency/linear", "progressive-protein/linear", "progressive-protein/affine", "progressive-protein/nj", "iterative-protein/linear" }));
    }

    [Test]
    public void ExactSpeedUps_GiveTheAlignmentsOfTheCodeBeforeThem()
    {
        var failures = new List<string>();
        for (int k = 0; k < Baseline.Value.Count; k++)
        {
            var c = Baseline.Value[k];
            var result = Run(c);
            if (!result.AlignedSequences.SequenceEqual(c.Rows) || result.TotalScore.ToString("R", System.Globalization.CultureInfo.InvariantCulture) != c.Score)
                failures.Add($"case {k} {c.Kind} {c.Options}");
        }

        Assert.That(failures, Is.Empty, string.Join("; ", failures.Take(10)));
    }

    private static string Random(Random rng, string alphabet, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = alphabet[rng.Next(alphabet.Length)];
        return new string(chars);
    }

    [Test]
    public void IdentityCountWithoutTraceback_EqualsTheIdentitiesOfTheAlignment()
    {
        var rng = new Random(7808);
        var scalar = new PairScoring(new ScoringMatrix(2, -1, 0, -2));
        for (int round = 0; round < 400; round++)
        {
            string alphabet = round % 3 == 0 ? "ACGT" : round % 3 == 1 ? "AC" : "ACGTN-";
            string a = Random(rng, alphabet, rng.Next(0, 80));
            string b = rng.Next(2) == 0 ? Random(rng, alphabet, rng.Next(0, 80)) : a.Length > 3 ? a.Remove(rng.Next(0, a.Length - 1), 2) + Random(rng, alphabet, rng.Next(0, 4)) : a;
            var alignment = SequenceAligner.LinearAlignCore(a, b, scalar, EndGaps.None);
            int expected = 0;
            for (int i = 0; i < alignment.AlignedSequence1.Length; i++)
            {
                char c1 = alignment.AlignedSequence1[i];
                if (c1 != '-' && c1 == alignment.AlignedSequence2[i])
                    expected++;
            }

            Assert.That(SequenceAligner.LinearIdentityCount(a, b, scalar), Is.EqualTo(expected), $"a={a} b={b}");
        }
    }

    [Test]
    public void IdentityCountWithoutTraceback_EqualsTheIdentitiesOfTheAlignment_UnderAProteinMatrix()
    {
        var rng = new Random(7809);
        for (int round = 0; round < 200; round++)
        {
            string a = Random(rng, "ACDEFGHIKLMNPQRSTVWY", rng.Next(1, 70));
            string b = rng.Next(2) == 0 ? Random(rng, "ACDEFGHIKLMNPQRSTVWY", rng.Next(1, 70)) : new string(a.Select(c => rng.Next(5) == 0 ? 'W' : c).ToArray());
            var scoring = PairScoring.Create(MatrixScoring.Emboss(SubstitutionMatrix.Blosum62), a, b);
            var alignment = SequenceAligner.LinearAlignCore(a, b, scoring, EndGaps.None);
            int expected = 0;
            for (int i = 0; i < alignment.AlignedSequence1.Length; i++)
            {
                char c1 = alignment.AlignedSequence1[i];
                if (c1 != '-' && c1 == alignment.AlignedSequence2[i])
                    expected++;
            }

            Assert.That(SequenceAligner.LinearIdentityCount(a, b, scoring), Is.EqualTo(expected), $"a={a} b={b}");
        }
    }

    [Test]
    public void KTupleDistance_IsZeroForIdenticalAndGrowsWithTheDivergence()
    {
        var rng = new Random(7810);
        string basis = Random(rng, "ACGT", 1000);
        double previous = -1;
        foreach (double divergence in new[] { 0.0, 0.05, 0.15, 0.30 })
        {
            var mutated = new string(basis.Select(c => rng.NextDouble() < divergence ? "ACGT"[rng.Next(4)] : c).ToArray());
            var codes1 = SequenceAligner.KTupleCodes(basis, SequenceAligner.KTupleDnaLength)!;
            var codes2 = SequenceAligner.KTupleCodes(mutated, SequenceAligner.KTupleDnaLength)!;
            double d = SequenceAligner.KTupleDistance(codes1, codes2, basis.Length, mutated.Length, SequenceAligner.KTupleDnaLength);
            if (divergence == 0.0)
                Assert.That(d, Is.EqualTo(0.0));
            Assert.That(d, Is.GreaterThan(previous));
            // the substituted fraction is divergence * 3/4 of the residues, the estimate is on that identity scale (within sampling noise)
            Assert.That(d, Is.EqualTo(divergence * 0.75).Within(0.06), $"divergence {divergence}");
            previous = d;
        }

        Assert.That(SequenceAligner.KTupleCodes("ACGT", SequenceAligner.KTupleDnaLength), Is.Null, "shorter than 2k: the alignment distance is kept");
    }

    [Test]
    public void BelowTheLimit_TheDistancesAreTheAlignments_AboveItTheKTupleEstimateGuidesTheTree()
    {
        var rng = new Random(7811);
        string basis = Random(rng, "ACGT", 300);
        var family = Enumerable.Range(0, 8).Select(_ => new DnaSequence(new string(basis.Select(c => rng.NextDouble() < 0.08 ? "ACGT"[rng.Next(4)] : c).ToArray()))).ToList();

        var never = SequenceAligner.MultipleAlignProgressive(family, null, new MsaOptions { FastDistanceCells = 0 });
        var defaultOptions = SequenceAligner.MultipleAlignProgressive(family);
        Assert.That(defaultOptions.AlignedSequences, Is.EqualTo(never.AlignedSequences), "28 pairs of 90 000 cells are far below the default limit of 400 million cells: the alignment distances");
        Assert.That(new MsaOptions().FastDistanceCells, Is.EqualTo(400_000_000L));

        var ktuple = SequenceAligner.MultipleAlignProgressive(family, null, new MsaOptions { FastDistanceCells = 1 });
        for (int i = 0; i < family.Count; i++)
            Assert.That(ktuple.AlignedSequences[i].Replace("-", ""), Is.EqualTo(family[i].Sequence), "every row is its sequence");
        Assert.That(ktuple.TotalScore, Is.GreaterThanOrEqualTo((int)(never.TotalScore * 0.97)), "on related sequences the k-tuple tree gives an alignment of (nearly) the same sum-of-pairs score");
        Assert.That(ktuple.TotalScore, Is.LessThanOrEqualTo(never.TotalScore + 1000));
    }
}
