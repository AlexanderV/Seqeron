using Seqeron.Genomics.Chromosome;
using static Seqeron.Genomics.Chromosome.GenomeAssemblyAnalyzer;

namespace Seqeron.Genomics.Tests.Unit.Chromosome;

/// <summary>
/// B10 R1-WP1 (A15-7 / A11-6): the private <c>ToUpperInvariant</c> folds of GenomeAssemblyAnalyzer were replaced by the
/// canonical ASCII-only <c>SequenceCase.ToUpperAscii</c> (U+017F 'ſ' must not become 'S'; B01 F21 / B02), the
/// per-position <c>Substring</c> unit lookup of <see cref="FindRepetitiveRegions"/> by <c>KmerEnumerator</c> + span
/// dictionary lookup, and the private <c>IsAcgt</c> by the canonical <c>IsValidDna</c> (same on upper-case input).
/// Each test carries a non-ASCII letter ('ſ' U+017F, 'İ' U+0130).
/// </summary>
[TestFixture]
public class GenomeAssemblyAnalyzer_AsciiFold_Tests
{
    private static string Rep(string u, int n) => string.Concat(Enumerable.Repeat(u, n));

    // WindowMasker: a unit with a non-ACGT symbol is invalid (its window is skipped). A non-ASCII letter is a non-ACGT
    // symbol exactly like 'N', so replacing one base by 'ſ' or 'İ' gives the same regions as replacing it by 'N'.
    [TestCase('ſ')]
    [TestCase('İ')]
    public void FindRepetitiveRegions_NonAsciiLetter_ActsAsNonAcgtSymbol(char letter)
    {
        string clean = Rep("ACGGTCATTGCA", 40);
        string withN = clean.Remove(240, 1).Insert(240, "N");
        string withLetter = clean.Remove(240, 1).Insert(240, letter.ToString());

        var expected = FindRepetitiveRegions(new[] { ("s", withN) }, 12, 3, 40).ToList();
        var actual = FindRepetitiveRegions(new[] { ("s", withLetter) }, 12, 3, 40).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(expected, Is.Not.Empty);
            Assert.That(actual, Is.EqualTo(expected));
        });
    }

    // Unit count equals the pooled forward + reverse-complement count (WindowMasker both strands); a lower-case sequence
    // is folded first (control for the fold itself).
    [Test]
    public void FindRepetitiveRegions_LowerCase_EqualsUpperCase()
    {
        string seq = Rep("ACGGTCATTGCA", 30);

        Assert.That(
            FindRepetitiveRegions(new[] { ("s", seq.ToLowerInvariant()) }, 12, 3, 40).ToList(),
            Is.EqualTo(FindRepetitiveRegions(new[] { ("s", seq) }, 12, 3, 40).ToList()));
    }

    // K-mer spectrum ids are exact: 'ſ' is not merged with 'S' (the old invariant fold made "ſ" and "S" a duplicate k-mer).
    [Test]
    public void EstimateCompletenessFromKmers_LongSAndS_AreDistinctKmers()
    {
        Assert.DoesNotThrow(() => EstimateCompletenessFromKmers(
            new[] { ("ſ", 3), ("S", 3) }, null));
        Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(
            new[] { ("S", 3), ("s", 3) }, null), "ASCII lower-case 's' is still the duplicate of 'S'");
    }

    // Canonical evaluation: a non-ACGT spectrum k-mer is rejected and reported verbatim (not as 'ASA').
    [Test]
    public void EstimateCompletenessFromKmers_Canonical_NonAsciiKmer_RejectedAsNotAcgt()
    {
        var ex = Assert.Throws<ArgumentException>(() => EstimateCompletenessFromKmers(
            new[] { ("AſA", 4) }, new[] { ("a", "ACGTACGTACGT") }));

        Assert.That(ex!.Message, Does.Contain("'AſA'").And.Contain("not ACGT"));
    }

    // Marker search: a non-ACGT symbol (here 'ſ' in the marker, 'S' in the assembly) does not seed and does not break the
    // full-length local alignment; the marker is upper-cased through the ASCII fold at the call site. (SequenceAligner,
    // a B09 file, still folds with ToUpperInvariant internally, so ſ vs S scores as a match there — cross-batch note in
    // the WP report; this test therefore only pins the hit geometry, not the identity.)
    [Test]
    public void FindMarkerHits_NonAsciiMarkerBase_FullLengthHitFound()
    {
        const string core = "TGCCTGGTACATCCGCGAAATGCAGTCAAACCAGGCACTAGTCAATATTCAACATTTGAATGCCAAACTATACC";
        string left = "CCCAAAGGCCCGCCGATAAAAGCACGGCACGACACTCCGTCGAGTTGTTGTGCAATCG";
        string right = "ATGCGCCAGTATATTGCTCTTACTAACCGAAGAGTGGGGCACAATTAATCACTTAAAC";
        string assembly = left + core + "S" + core + right;
        string marker = core + "\u017F" + core;

        var hits = FindMarkerHits(marker, new[] { ("a", assembly) });

        Assert.That(hits.Any(h => h.Strand == '+' && h.Start == left.Length && h.End == left.Length + marker.Length - 1
                                  && h.MarkerStart == 0 && h.MarkerEnd == marker.Length - 1 && h.Coverage == 1.0));
    }

    // A lower-case marker over a lower-case assembly is the exact-match control (ASCII fold).
    [Test]
    public void FindMarkerHits_LowerCaseMarkerAndAssembly_FullIdentity()
    {
        const string core = "TGCCTGGTACATCCGCGAAATGCAGTCAAACCAGGCACTAGTCAATATTCAACATTTGAATGCCAAACTATACC";
        var hit = FindMarkerHits(core.ToLowerInvariant(), new[] { ("a", "ccc" + core.ToLowerInvariant() + "ggg") }).Single();

        Assert.That(hit.Identity, Is.EqualTo(1.0));
    }

    // MUM anchors (nucmer): a query symbol outside A/C/G/T becomes the never-matching '$' sentinel; a non-ASCII letter is
    // such a symbol, so 'ſ' (as 'N') splits the match identically.
    [Test]
    public void FindSyntenicBlocks_NonAsciiQueryLetter_ActsAsNonAcgt()
    {
        var rnd = new Random(8123);
        string refSeq = new string(Enumerable.Range(0, 400).Select(_ => "ACGT"[rnd.Next(4)]).ToArray());
        string withN = refSeq.Remove(200, 1).Insert(200, "N");
        string withLetter = refSeq.Remove(200, 1).Insert(200, "ſ");

        var expected = FindSyntenicBlocks(new[] { ("r", refSeq) }, new[] { ("q", withN) }, minBlockSize: 20, kmerSize: 12).ToList();
        var actual = FindSyntenicBlocks(new[] { ("r", refSeq) }, new[] { ("q", withLetter) }, minBlockSize: 20, kmerSize: 12).ToList();

        Assert.Multiple(() =>
        {
            Assert.That(expected, Is.Not.Empty);
            Assert.That(actual, Is.EqualTo(expected));
        });
    }

    // Local quality: the non-ASCII window symbol is excluded from G+C and A+C+G+T (Biopython gc_fraction "remove"), so
    // the window "ACGTſ" has GC 0.5 (fraction 0.5) and no N.
    [Test]
    public void CalculateLocalQuality_NonAsciiLetter_ExcludedFromGc()
    {
        var q = CalculateLocalQuality(new[] { ("s", "ACGTſ") }, 5, 5).Single();

        Assert.Multiple(() =>
        {
            Assert.That(q.GcContent, Is.EqualTo(0.5).Within(1e-12));
            Assert.That(q.NCount, Is.EqualTo(0));
        });
    }
}
