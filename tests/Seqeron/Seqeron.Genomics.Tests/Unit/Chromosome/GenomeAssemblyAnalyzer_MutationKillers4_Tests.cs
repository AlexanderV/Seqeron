using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// ASSEMBLY-STATS-001 mutation killers (batch 4): exact-value pins for the internal arithmetic of the
/// comparison / completeness / repeat / suspicious-region methods (shared-kmer fractions, fragmented vs
/// complete classification, contig length boundaries, k-mer size/peak edges, imperfect tandem purity,
/// low-complexity score formula).
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_MutationKillers4_Tests
{
    private const double Tol = 1e-9;

    private static string RandomDna(int length, int seed)
    {
        var rng = new Random(seed);
        const string bases = "ACGT";
        return string.Concat(Enumerable.Range(0, length).Select(_ => bases[rng.Next(4)]));
    }

    [Test]
    public void CompareAssemblies_PartialOverlap_ExactFractions()
    {
        // canonical 4-mers: a1 {ACGT, CGTA (=TACG), GTAC}; a2 {ACGT, CGTA, GTAA (=TTAC), TAAA (=TTTA), AAAA};
        // shared = 2, union = 6 ⇒ J = 1/3 ⇒ Mash D = −ln(0.5)/4 (Ondov 2016 eq. 4).
        var a1 = new[] { ("s", "ACGTACGT") };
        var a2 = new[] { ("s", "ACGTAAAA") };
        var cmp = CompareAssemblies(a1, a2, kmerSize: 4);
        Assert.That(cmp.Containment1, Is.EqualTo(2.0 / 3).Within(Tol));
        Assert.That(cmp.Containment2, Is.EqualTo(2.0 / 5).Within(Tol));
        Assert.That(cmp.JaccardIndex, Is.EqualTo(1.0 / 3).Within(Tol));
        Assert.That(cmp.MashDistance, Is.EqualTo(Math.Log(2) / 4).Within(Tol));
    }

    [Test]
    public void AssessCompleteness_PartialCoverage_IsFragmented()
    {
        // 100-bp unique marker; assembly carries only its first 80 bp ⇒ ~71% k-mer coverage:
        // above the relaxed hit threshold but below the 0.9 completeness cutoff ⇒ Fragmented.
        string marker = RandomDna(100, 11);
        var r = AssessCompleteness(new[] { ("scaf", marker.Substring(0, 80)) }, new[] { ("g", marker) });
        Assert.That(r.Fragmented, Is.EqualTo(1));
        Assert.That(r.Complete, Is.EqualTo(0));
        Assert.That(r.Missing, Is.EqualTo(0));
    }

    [Test]
    public void ExtractContigs_LengthBoundaryIsInclusive()
    {
        // A 3-bp contig is kept at minContigLength 3 (inclusive '>=') but dropped at 4.
        // minGapLength 1: the single N is a gap (QUAST split_by_ns with Ns_break_threshold = 1).
        Assert.That(ExtractContigs(new[] { ("s", "AAAN") }, minContigLength: 3, minGapLength: 1).Single(),
            Is.EqualTo(("s_contig1", "AAA")));
        Assert.That(ExtractContigs(new[] { ("s", "AAAN") }, minContigLength: 4, minGapLength: 1).ToList(), Is.Empty);
        // Trailing contig after a gap also honours the inclusive bound.
        Assert.That(ExtractContigs(new[] { ("s", "NAAA") }, minContigLength: 3, minGapLength: 1).Single(),
            Is.EqualTo(("s_contig1", "AAA")));
    }

    [Test]
    public void EstimateCompletenessFromKmers_AllSingletons_AndExpectedCoverage()
    {
        // A one-line histogram has no valley (Merqury filt = 0): nothing is an error k-mer, the peak is
        // count 1 ⇒ error 0, size 2/1 = 2; no assembly ⇒ completeness NaN (oracle; B10 F44–F46).
        var allSingleton = EstimateCompletenessFromKmers(new[] { ("A", 1), ("B", 1) });
        Assert.That(allSingleton.Completeness, Is.NaN);
        Assert.That(allSingleton.ErrorRate, Is.EqualTo(0.0).Within(Tol));
        Assert.That(allSingleton.EstimatedGenomeSize, Is.EqualTo(2L));

        // expectedCoverage overrides the computed peak: size = round(11 / 2) = round(5.5) = 6 (R round, half-even).
        var forced = EstimateCompletenessFromKmers(
            new[] { ("AAA", 4), ("CCC", 4), ("GGG", 2), ("TTT", 1) }, expectedCoverage: 2);
        Assert.That(forced.EstimatedGenomeSize, Is.EqualTo(6L));
    }

    [Test]
    public void FindTandemRepeats_SingleSubstitution_ExactTrfPurity()
    {
        // (ACGTT)×20 with one T→A substitution in the 11th copy and 100 bp: TRF (match 2, mismatch −7)
        // aligns the whole array, 20 copies, one mismatch touching two adjacent-copy comparisons of 95:
        // purity 93/95; score 2·99 − 7 = 191.
        string unit = "ACGTT";
        string seq = string.Concat(Enumerable.Repeat(unit, 10)) + "ACGTA" + string.Concat(Enumerable.Repeat(unit, 9));
        var r = FindTandemRepeats(new[] { ("s", seq) }, minUnitLength: 5, minCopies: 3).First();
        Assert.That(r.Unit, Is.EqualTo("ACGTT"));
        Assert.That(r.Copies, Is.EqualTo(20.0).Within(Tol));
        Assert.That(r.Purity, Is.EqualTo(93.0 / 95.0).Within(Tol));
    }

    [Test]
    public void FindSuspiciousRegions_PolyA_ExactLowComplexityScore()
    {
        // 600-bp poly(A): a 500-bp window has canonical linguistic complexity (Troyanskaya 2002, m = 10, a = 4)
        // ΣV_i / Σmin(4^i, 501−i) = 10 / (4+16+64+256+496+495+494+493+492+491) = 10/3301,
        // giving the peak low-complexity score 1 − (10/3301)/0.3.
        var r = FindSuspiciousRegions(new[] { ("s", new string('A', 600)) }).Single();
        Assert.That(r.Score, Is.EqualTo(1.0 - (10.0 / 3301.0) / 0.3).Within(1e-9));
    }
}
