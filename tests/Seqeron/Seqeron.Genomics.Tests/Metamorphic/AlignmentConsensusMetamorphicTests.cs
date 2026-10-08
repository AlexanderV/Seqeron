namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic relations for the B05 follow-up consensus algorithms:
///   • EMB-CASE     lower-casing the input does not change GenerateEmbossConsensus (EMBOSS upper-cases input).
///   • EMB-GAPS     rewriting gaps '-' as '.' or '~' does not change it (EMBOSS gap normalisation).
///   • EMB-SCALE    multiplying every weight, plurality and setcase by 2 does not change it (exact in float).
///   • EMB-DUP      over gap-free A/C/G/T rows, duplicating every row equals weight 2 per row: the extra
///                  self-pair adds the same EDNAFULL diagonal (+5) to every score, so the argmax, the
///                  positive-match weights and the positional tie-break are unchanged.
///   • DUMB-PERM    row permutation does not change GenerateDumbConsensus (counts only; ties → ambiguous).
///   • DUMB-MONO    raising the threshold only turns residues into the ambiguous symbol.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class AlignmentConsensusMetamorphicTests
{
    private static string[] RandomRows(Random rng, int n, int length, string alphabet) =>
        Enumerable.Range(0, n).Select(_ => new string(Enumerable.Range(0, length)
            .Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray())).ToArray();

    [Test]
    public void Emboss_CaseGapAndWeightScaling_Invariant()
    {
        var rng = new Random(6601);
        for (int trial = 0; trial < 200; trial++)
        {
            bool protein = trial % 2 == 1;
            var type = protein ? ConsensusResidueType.Protein : ConsensusResidueType.Nucleotide;
            var rows = RandomRows(rng, rng.Next(2, 9), rng.Next(1, 30), protein ? "ACDEFGHIKLMNPQRSTVWY-" : "ACGTN-");
            float[] weights = rows.Select(_ => new[] { 0.25f, 0.5f, 1f, 1.5f, 2f }[rng.Next(5)]).ToArray();
            float total = weights.Sum();
            float plurality = MathF.Round(total * (float)rng.NextDouble() * 4) / 4, setcase = MathF.Round(total * (float)rng.NextDouble() * 4) / 4;
            int identity = rng.Next(0, 3);

            string baseline = MotifFinder.GenerateEmbossConsensus(rows, type, plurality, identity, setcase, weights);
            Assert.That(MotifFinder.GenerateEmbossConsensus(rows.Select(r => r.ToLowerInvariant()), type, plurality, identity, setcase, weights),
                Is.EqualTo(baseline), "case");
            Assert.That(MotifFinder.GenerateEmbossConsensus(rows.Select((r, i) => r.Replace('-', i % 2 == 0 ? '.' : '~')), type, plurality, identity, setcase, weights),
                Is.EqualTo(baseline), "gaps");
            Assert.That(MotifFinder.GenerateEmbossConsensus(rows, type, 2 * plurality, identity, 2 * setcase, weights.Select(w => 2 * w).ToArray()),
                Is.EqualTo(baseline), "scale");
        }
    }

    [Test]
    public void Emboss_DuplicatedRows_EqualDoubleWeight()
    {
        var rng = new Random(6602);
        for (int trial = 0; trial < 200; trial++)
        {
            var rows = RandomRows(rng, rng.Next(2, 7), rng.Next(1, 25), "ACGT");
            float plurality = rng.Next(0, 2 * rows.Length + 1), setcase = rng.Next(0, 2 * rows.Length + 1);
            string weighted = MotifFinder.GenerateEmbossConsensus(rows, plurality: plurality, setcase: setcase,
                weights: rows.Select(_ => 2f).ToArray());
            // Interleave duplicates so every row keeps its relative order (the highest-score tie-break is positional).
            string duplicated = MotifFinder.GenerateEmbossConsensus(rows.SelectMany(r => new[] { r, r }), plurality: plurality, setcase: setcase);
            Assert.That(duplicated, Is.EqualTo(weighted), string.Join(",", rows));
        }
    }

    [Test]
    public void Dumb_RowPermutationInvariant_AndThresholdMonotone()
    {
        var rng = new Random(6603);
        for (int trial = 0; trial < 300; trial++)
        {
            var rows = RandomRows(rng, rng.Next(1, 10), rng.Next(0, 25), trial % 2 == 0 ? "ACGT-" : "ACDEFGHIKL-.");
            double t1 = rng.NextDouble(), t2 = rng.NextDouble();
            if (t1 > t2) (t1, t2) = (t2, t1);
            string low = MotifFinder.GenerateDumbConsensus(rows, t1, '?');
            var shuffled = rows.OrderBy(_ => rng.Next()).ToArray();
            Assert.That(MotifFinder.GenerateDumbConsensus(shuffled, t1, '?'), Is.EqualTo(low), "permutation");
            string high = MotifFinder.GenerateDumbConsensus(rows, t2, '?');
            for (int k = 0; k < low.Length; k++)
                Assert.That(high[k] == low[k] || high[k] == '?', Is.True, $"monotone col {k}");
        }
    }
}
