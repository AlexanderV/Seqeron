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
    private static readonly Dictionary<char, string> IupacSets = new()
    {
        ['A'] = "A", ['C'] = "C", ['G'] = "G", ['T'] = "T", ['R'] = "AG", ['Y'] = "CT", ['S'] = "CG", ['W'] = "AT",
        ['K'] = "GT", ['M'] = "AC", ['B'] = "CGT", ['D'] = "AGT", ['H'] = "ACT", ['V'] = "ACG", ['N'] = "ACGT",
    };

    private static string[] RandomRows(Random rng, int n, int length, string alphabet) =>
        Enumerable.Range(0, n).Select(_ => new string(Enumerable.Range(0, length)
            .Select(_ => alphabet[rng.Next(alphabet.Length)]).ToArray())).ToArray();

    [Test]
    public void C1_DefaultThreshold_BitIdentical()
    {
        var rng = new Random(5501);
        for (int trial = 0; trial < 500; trial++)
        {
            var rows = RandomRows(rng, rng.Next(1, 15), rng.Next(0, 30), trial % 2 == 0 ? "ACGT" : "ACGTacgtN-");
            Assert.That(MotifFinder.GenerateConsensus(rows, 0.25), Is.EqualTo(MotifFinder.GenerateConsensus(rows)));
        }
    }

    [Test]
    public void C2_ThresholdNesting_SubsetOfBaseSets()
    {
        var rng = new Random(5502);
        for (int trial = 0; trial < 300; trial++)
        {
            var rows = RandomRows(rng, rng.Next(1, 12), rng.Next(1, 20), "ACGTN-");
            double t1 = rng.NextDouble(), t2 = rng.NextDouble();
            if (t1 > t2) (t1, t2) = (t2, t1);
            string low = MotifFinder.GenerateConsensus(rows, t1), high = MotifFinder.GenerateConsensus(rows, t2);
            for (int i = 0; i < low.Length; i++)
            {
                Assert.That(IupacSets.ContainsKey(low[i]) && IupacSets.ContainsKey(high[i]), Is.True);
                bool lowHasBase = rows.Any(r => "ACGT".Contains(r[i]));
                if (lowHasBase)
                    Assert.That(IupacSets[high[i]].All(b => IupacSets[low[i]].Contains(b)), Is.True,
                        $"col {i}: {high[i]}@{t2} ⊄ {low[i]}@{t1}");
            }
        }
    }

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
            string cons = MotifFinder.GenerateEmbossConsensus(rows, type, plurality: (float)rng.Next(0, 6), identity: rng.Next(0, 3));
            Assert.That(cons.Length, Is.EqualTo(rows[0].Length));
            for (int k = 0; k < cons.Length; k++)
            {
                char u = char.ToUpperInvariant(cons[k]);
                Assert.That(u == nocon || rows.Any(r => char.ToUpperInvariant(r[k]) == u), Is.True, $"col {k} '{cons[k]}'");
            }

            string upper = MotifFinder.GenerateEmbossConsensus(rows, type, plurality: 0f, setcase: -1f);
            Assert.That(upper, Is.EqualTo(upper.ToUpperInvariant()));
        }
    }

    [Test]
    public void C4_Dumb_MajorityOrAmbiguous()
    {
        var rng = new Random(5504);
        for (int trial = 0; trial < 300; trial++)
        {
            var rows = RandomRows(rng, rng.Next(1, 10), rng.Next(0, 25), "ACGT-.");
            double threshold = rng.NextDouble();
            string cons = MotifFinder.GenerateDumbConsensus(rows, threshold, '?');
            for (int k = 0; k < cons.Length; k++)
            {
                if (cons[k] == '?') continue;
                var column = rows.Select(r => r[k]).Where(c => c != '-' && c != '.').ToList();
                int count = column.Count(c => c == cons[k]);
                Assert.That(column.GroupBy(c => c).Count(g => g.Count() >= count), Is.EqualTo(1), "unique maximum");
                Assert.That((double)count / column.Count, Is.GreaterThanOrEqualTo(threshold));
            }
        }
    }
}
