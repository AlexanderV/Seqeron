namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Heavy-tier metamorphic relations for the review 2026-09 B05 audit additions on MOTIF-GENERATE-001 /
/// MOTIF-CONS-001 (docs/Validation/review-2026-09/B05.md F28: DECIPHER <c>ConsensusSequence</c> port, EMBOSS
/// <c>cons</c> Auto residue type + ragged-row padding).
///   • DEC-PERM   permuting the rows does not change GenerateDecipherConsensus (DECIPHER tallies per-column
///                fractions; rows over A/C/G/T and gaps only, so every increment is the same 1/n and the sums are
///                order-independent in floating point).
///   • DEC-CASE   lower-casing the input does not change it (Biostrings upper-cases).
///   • DEC-RNA    T → U with DecipherSequenceType.Rna gives the DNA consensus with T → U.
///   • DEC-ORDER  is unaffected by appending a row consisting only of gaps when terminal gaps are excluded
///                (the row contributes nothing; every column keeps its counted characters) — checked only on
///                columns that the other rows reach.
///   • EMB-AUTO   ConsensusResidueType.Auto == Nucleotide for rows over A/C/G/T/N and gaps, == Protein when the
///                first row contains a non-nucleotide amino-acid letter and no row contains X/? (EMBOSS
///                ajSeqTypeGapnucS on the first sequence).
///   • EMB-PAD    padRaggedRows = true on rectangular input == the original overload; on ragged input ==
///                the original overload on rows explicitly right-padded with '-' (ajSeqsetFill).
/// Fixed seeds, bounded sizes.
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ConsensusAuditMetamorphicTests
{
    private static string RandomRow(Random rng, int length, string alphabet) =>
        new(Enumerable.Range(0, length).Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray());

    private static string[] RandomRows(Random rng, int n, int length, string alphabet) =>
        Enumerable.Range(0, n).Select(_ => RandomRow(rng, length, alphabet)).ToArray();

    private static string Decipher(IEnumerable<string> rows, (double Threshold, bool Ambiguity, bool NonLetters, bool Terminal) o,
        DecipherSequenceType type = DecipherSequenceType.Dna) =>
        MotifFinder.GenerateDecipherConsensus(rows, type, o.Threshold, o.Ambiguity, '+', null, o.NonLetters, o.Terminal);

    private static (double, bool, bool, bool) RandomOptions(Random rng) =>
        (new[] { 0.0, 0.05, 0.1, 0.25, 0.4, 0.6 }[rng.Next(6)], rng.Next(2) == 0, rng.Next(2) == 0, rng.Next(2) == 0);

    [Test]
    public void Decipher_RowPermutation_Case_AndRna_Invariant()
    {
        var rng = new Random(280_001);
        for (int trial = 0; trial < 400; trial++)
        {
            var rows = RandomRows(rng, rng.Next(1, 10), rng.Next(0, 30), trial % 2 == 0 ? "ACGT-" : "AAACG-");
            var o = RandomOptions(rng);
            string baseline = Decipher(rows, o);
            string msg = $"{string.Join(",", rows)} {o}";

            var perm = rows.OrderBy(_ => rng.Next()).ToArray();
            Assert.That(Decipher(perm, o), Is.EqualTo(baseline), "perm " + msg);
            Assert.That(Decipher(rows.Select(r => r.ToLowerInvariant()), o), Is.EqualTo(baseline), "case " + msg);
            Assert.That(Decipher(rows.Select(r => r.Replace('T', 'U')), o, DecipherSequenceType.Rna),
                Is.EqualTo(baseline.Replace('T', 'U')), "rna " + msg);
        }
    }

    [Test]
    public void Decipher_RaggedRows_PermutationInvariant_LengthOfLongest()
    {
        var rng = new Random(280_002);
        for (int trial = 0; trial < 300; trial++)
        {
            var rows = Enumerable.Range(0, rng.Next(1, 8)).Select(_ => RandomRow(rng, rng.Next(0, 25), "ACGT-")).ToArray();
            var o = RandomOptions(rng);
            string baseline = Decipher(rows, o);
            Assert.That(baseline.Length, Is.EqualTo(rows.Max(r => r.Length)));
            Assert.That(Decipher(rows.Reverse(), o), Is.EqualTo(baseline), $"{string.Join(",", rows)} {o}");
        }
    }

    [Test]
    public void Decipher_AllGapRow_WithoutTerminalGaps_DoesNotChangeConsensus()
    {
        var rng = new Random(280_003);
        for (int trial = 0; trial < 300; trial++)
        {
            int length = rng.Next(1, 25);
            var rows = RandomRows(rng, rng.Next(1, 8), length, "ACGT-");
            var o = RandomOptions(rng);
            o.Item4 = false; // terminal gaps excluded: an all-gap row has nothing but terminal gaps
            string baseline = Decipher(rows, o);
            string withGapRow = Decipher(rows.Append(new string('-', length)), o);
            Assert.That(withGapRow, Is.EqualTo(baseline), $"{string.Join(",", rows)} {o}");
        }
    }

    [Test]
    public void Emboss_Auto_EqualsExplicitType_WhenUnambiguous()
    {
        var rng = new Random(280_004);
        const string protein = "ACDEFGHIKLMNPQRSTVWY-";
        const string proteinOnly = "EFILPQ";
        for (int trial = 0; trial < 300; trial++)
        {
            int length = rng.Next(1, 30);
            float? plurality = rng.Next(3) == 0 ? null : rng.Next(0, 6);
            int identity = rng.Next(0, 3);

            var dna = RandomRows(rng, rng.Next(2, 8), length, "ACGTN-");
            Assert.That(MotifFinder.GenerateEmbossConsensus(dna, ConsensusResidueType.Auto, plurality, identity),
                Is.EqualTo(MotifFinder.GenerateEmbossConsensus(dna, ConsensusResidueType.Nucleotide, plurality, identity)),
                "dna " + string.Join(",", dna));

            var prot = RandomRows(rng, rng.Next(2, 8), length, protein);
            var first = prot[0].ToCharArray();
            first[rng.Next(length)] = proteinOnly[rng.Next(proteinOnly.Length)];
            prot[0] = new string(first);
            Assert.That(MotifFinder.GenerateEmbossConsensus(prot, ConsensusResidueType.Auto, plurality, identity),
                Is.EqualTo(MotifFinder.GenerateEmbossConsensus(prot, ConsensusResidueType.Protein, plurality, identity)),
                "protein " + string.Join(",", prot));
        }
    }

    [Test]
    public void Emboss_PadRaggedRows_EqualsExplicitPadding()
    {
        var rng = new Random(280_005);
        for (int trial = 0; trial < 300; trial++)
        {
            bool protein = trial % 3 == 0;
            var type = protein ? ConsensusResidueType.Protein : ConsensusResidueType.Nucleotide;
            string alphabet = protein ? "ACDEFGHIKLMNPQRSTVWY-" : "ACGTN-";
            float? plurality = rng.Next(3) == 0 ? null : rng.Next(0, 6);

            var rect = RandomRows(rng, rng.Next(2, 8), rng.Next(1, 25), alphabet);
            Assert.That(MotifFinder.GenerateEmbossConsensus(rect, true, type, plurality),
                Is.EqualTo(MotifFinder.GenerateEmbossConsensus(rect, type, plurality)), "rect");

            var ragged = Enumerable.Range(0, rng.Next(2, 8)).Select(_ => RandomRow(rng, rng.Next(1, 25), alphabet)).ToArray();
            int max = ragged.Max(r => r.Length);
            var padded = ragged.Select(r => r.PadRight(max, '-')).ToArray();
            Assert.That(MotifFinder.GenerateEmbossConsensus(ragged, true, type, plurality),
                Is.EqualTo(MotifFinder.GenerateEmbossConsensus(padded, type, plurality)), "ragged " + string.Join(",", ragged));
            if (ragged.Any(r => r.Length != max))
                Assert.That(() => MotifFinder.GenerateEmbossConsensus(ragged, type, plurality), NUnit.Framework.Throws.InstanceOf<ArgumentException>());
        }
    }
}
