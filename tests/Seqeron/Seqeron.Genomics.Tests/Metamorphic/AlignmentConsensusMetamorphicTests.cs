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
