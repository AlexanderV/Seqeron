// GENOMIC-MOTIFS-001 (B09 A2-35) — string overloads of FindMotif / FindKnownMotifs / FindKnownMotifHits with N / IUPAC sequences.
// Contract: a sequence symbol other than A/C/G/T matches nothing (not even a motif N) and costs a mismatch when a mismatch budget is given — the B05 contract of
// MotifFinder.FindDegenerateMotif(string); on an A/C/G/T string every method equals the DnaSequence overload. Oracle: in-test brute force + the DnaSequence overloads.

namespace Seqeron.Genomics.Tests.Unit.Analysis;

[TestFixture]
public class GenomicAnalyzer_StringOverloads_Tests
{
    private static readonly Dictionary<char, string> Iupac = new()
    {
        ['A'] = "A", ['C'] = "C", ['G'] = "G", ['T'] = "T", ['R'] = "AG", ['Y'] = "CT", ['S'] = "CG", ['W'] = "AT", ['K'] = "GT", ['M'] = "AC",
        ['B'] = "CGT", ['D'] = "AGT", ['H'] = "ACT", ['V'] = "ACG", ['N'] = "ACGT",
    };

    private static int Mismatches(string seq, int start, string motif, bool iupac)
    {
        int mm = 0;
        for (int j = 0; j < motif.Length; j++)
        {
            char s = seq[start + j];
            bool ok = iupac ? Iupac[motif[j]].Contains(s) && "ACGT".Contains(s) : s == motif[j] && "ACGT".Contains(s);
            if (!ok) mm++;
        }
        return mm;
    }

    private static string[] Brute(string seq, string motif, int k, bool iupac)
        => Enumerable.Range(0, Math.Max(0, seq.Length - motif.Length + 1))
            .Where(i => Mismatches(seq, i, motif, iupac) <= k)
            .Select(i => $"{i}:{Mismatches(seq, i, motif, iupac)}").ToArray();

    private static string Reverse(string motif) => DnaSequence.GetReverseComplementString(motif);

    [Test]
    public void FindKnownMotifHits_String_MatchesBruteForce_OnNSequences()
    {
        var rnd = new Random(20261124);
        string[] motifs = ["ACG", "GAATTC", "TATAWAW", "GGNCC", "RYR", "AAAA", "CAT", "NN"];
        int rows = 0, withN = 0;
        for (int iter = 0; iter < 400; iter++)
        {
            string alpha = new[] { "ACGT", "ACGTN", "acgtnr", "AC" }[rnd.Next(4)];
            string seq = new string(Enumerable.Range(0, rnd.Next(0, 40)).Select(_ => alpha[rnd.Next(alpha.Length)]).ToArray());
            string upper = seq.ToUpperInvariant();
            if (upper.Any(c => !"ACGT".Contains(c))) withN++;
            foreach (string motif in motifs.OrderBy(_ => rnd.Next()).Take(3))
            {
                foreach (bool iupac in new[] { false, true })
                {
                    if (!iupac && motif.Any(c => !"ACGT".Contains(c))) continue;
                    int k = rnd.Next(0, 3);
                    var hits = GenomicAnalyzer.FindKnownMotifHits(seq, [motif], bothStrands: true, allowIupac: iupac, maxMismatches: k);
                    var forward = hits.TryGetValue(motif, out var h) ? h.Where(x => !x.IsReverseStrand).Select(x => $"{x.Position}:{x.Mismatches}").ToArray() : [];
                    var reverse = hits.TryGetValue(motif, out var h2) ? h2.Where(x => x.IsReverseStrand).Select(x => $"{x.Position}:{x.Mismatches}").ToArray() : [];
                    Assert.That(forward, Is.EqualTo(Brute(upper, motif, k, iupac)), $"{seq} {motif} iupac={iupac} k={k} forward");
                    Assert.That(reverse, Is.EqualTo(Brute(upper, Reverse(motif), k, iupac)), $"{seq} {motif} iupac={iupac} k={k} reverse");
                    rows++;
                }
            }
        }

        Assert.That(rows, Is.GreaterThan(1500));
        Assert.That(withN, Is.GreaterThan(100));
    }

    [Test]
    public void FindMotifAndKnownMotifs_String_AreTheDnaSequenceResults_OnAcgt_AndSkipNWindows()
    {
        var rnd = new Random(77);
        for (int iter = 0; iter < 300; iter++)
        {
            string seq = new string(Enumerable.Range(0, rnd.Next(1, 60)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
            var dna = new DnaSequence(seq);
            string[] motifs = Enumerable.Range(0, 4).Select(_ => new string(Enumerable.Range(0, rnd.Next(1, 5)).Select(_ => "ACGT"[rnd.Next(4)]).ToArray())).ToArray();
            Assert.That(GenomicAnalyzer.FindMotif(seq.ToLowerInvariant(), motifs[0]), Is.EqualTo(GenomicAnalyzer.FindMotif(dna, motifs[0])), seq);
            var a = GenomicAnalyzer.FindKnownMotifs(seq, motifs);
            var b = GenomicAnalyzer.FindKnownMotifs(dna, motifs);
            Assert.That(a.Keys, Is.EquivalentTo(b.Keys));
            foreach (var key in a.Keys)
                Assert.That(a[key], Is.EqualTo(b[key]), $"{seq} {key}");
        }

        // N never matches: ANA is not found in ACANAC even though the symbols are equal
        Assert.That(GenomicAnalyzer.FindMotif("ACANAC", "ANA"), Is.Empty);
        Assert.That(GenomicAnalyzer.FindMotif("ACANACAGA", "AGA"), Is.EqualTo(new[] { 6 }));
        Assert.That(GenomicAnalyzer.FindKnownMotifs("ACNACGACG", new[] { "ACG", "CNA" }).Keys, Is.EqualTo(new[] { "ACG" }));
        Assert.That(GenomicAnalyzer.FindKnownMotifs("ACNACGACG", new[] { "ACG" })["ACG"], Is.EqualTo(new[] { 3, 6 }));
    }

    [Test]
    public void StringOverloads_Validation()
    {
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindMotif((string)null!, "A"));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindKnownMotifs((string)null!, new[] { "A" }));
        Assert.Throws<ArgumentNullException>(() => GenomicAnalyzer.FindKnownMotifs("ACGT", null!));
        Assert.Throws<ArgumentOutOfRangeException>(() => GenomicAnalyzer.FindKnownMotifHits("ACGT", new[] { "A" }, false, false, -1));
        Assert.That(GenomicAnalyzer.FindMotif("ACGT", string.Empty), Is.Empty);
        Assert.That(GenomicAnalyzer.FindMotif(string.Empty, "A"), Is.Empty);
        Assert.That(GenomicAnalyzer.FindKnownMotifHits("", new[] { "A" }), Is.Empty);
        Assert.Throws<ArgumentException>(() => GenomicAnalyzer.FindKnownMotifHits("ACGT", new[] { "AQA" }, allowIupac: true));
    }
}
