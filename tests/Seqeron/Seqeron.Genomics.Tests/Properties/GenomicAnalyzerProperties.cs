using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the high-level GenomicAnalyzer façade: common regions, known motifs,
/// ORFs, repeats, similarity, and tandem repeats.
///
/// Test Units: GENOMIC-COMMON-001, GENOMIC-MOTIFS-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Analysis")]
public class GenomicAnalyzerProperties
{
    private static string RandDna(Random rng, int len)
    {
        const string bases = "ACGT";
        var c = new char[len];
        for (int i = 0; i < len; i++) c[i] = bases[rng.Next(4)];
        return new string(c);
    }

    #region GENOMIC-COMMON-001: P: common region occurs in both inputs; R: positions valid; D: deterministic

    // FindCommonRegions / FindLongestCommonRegion locate substrings shared by two sequences via the
    // generalized suffix tree (Gusfield 1997). A reported region must occur at its stated positions
    // in BOTH sequences.

    /// <summary>Two sequences sharing a random 6-mer, plus that shared substring.</summary>
    private static Arbitrary<(string s1, string s2, string shared)> SharedRegionArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            string shared = RandDna(rng, 6);
            string s1 = RandDna(rng, 1 + rng.Next(5)) + shared + RandDna(rng, 1 + rng.Next(5));
            string s2 = RandDna(rng, 1 + rng.Next(5)) + shared + RandDna(rng, 1 + rng.Next(5));
            return (s1, s2, shared);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P + R): every reported common region occurs at its stated position in BOTH sequences,
    /// with valid coordinates and length ≥ minLength.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CommonRegions_OccurInBothInputs()
    {
        return Prop.ForAll(SharedRegionArbitrary(), input =>
        {
            var (s1, s2, _) = input;
            const int minLength = 3;
            var regions = GenomicAnalyzer.FindCommonRegions(new DnaSequence(s1), new DnaSequence(s2), minLength).ToList();
            bool ok = regions.All(r =>
                r.Length >= minLength &&
                r.PositionInFirst >= 0 && r.PositionInFirst + r.Length <= s1.Length &&
                r.PositionInSecond >= 0 && r.PositionInSecond + r.Length <= s2.Length &&
                s1.Substring(r.PositionInFirst, r.Length) == r.Sequence &&
                s2.Substring(r.PositionInSecond, r.Length) == r.Sequence);
            return ok.Label("a common region did not occur at its stated position in both inputs");
        });
    }

    /// <summary>
    /// INV-2 (P, positive control): the longest common region is at least as long as the embedded
    /// shared substring and occurs in both sequences.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CommonRegions_LongestCoversSharedSubstring()
    {
        return Prop.ForAll(SharedRegionArbitrary(), input =>
        {
            var (s1, s2, shared) = input;
            var lcr = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence(s1), new DnaSequence(s2));
            bool ok = !lcr.IsEmpty
                      && lcr.Length >= shared.Length
                      && s1.Contains(lcr.Sequence) && s2.Contains(lcr.Sequence);
            return ok.Label($"longest common region '{lcr.Sequence}' shorter than shared '{shared}' or not common");
        });
    }

    /// <summary>
    /// INV-3 (D): Common-region detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property CommonRegions_IsDeterministic()
    {
        return Prop.ForAll(SharedRegionArbitrary(), input =>
        {
            var (s1, s2, _) = input;
            var a = GenomicAnalyzer.FindCommonRegions(new DnaSequence(s1), new DnaSequence(s2), 3)
                .Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();
            var b = GenomicAnalyzer.FindCommonRegions(new DnaSequence(s1), new DnaSequence(s2), 3)
                .Select(r => (r.Sequence, r.PositionInFirst, r.PositionInSecond)).ToList();
            return a.SequenceEqual(b).Label("FindCommonRegions must be deterministic");
        });
    }

    /// <summary>
    /// INV-4 (boundary): sequences with no shared substring yield no common region.
    /// </summary>
    [Test]
    [Category("Property")]
    public void CommonRegions_NoOverlap_IsEmpty()
    {
        var lcr = GenomicAnalyzer.FindLongestCommonRegion(new DnaSequence("AAAAAA"), new DnaSequence("CCCCCC"));
        Assert.Multiple(() =>
        {
            Assert.That(lcr.IsEmpty, Is.True, "disjoint alphabets share no region");
            Assert.That(GenomicAnalyzer.FindCommonRegions(new DnaSequence("AAAAAA"), new DnaSequence("CCCCCC"), 3),
                Is.Empty);
        });
    }

    #endregion

    #region GENOMIC-MOTIFS-001: R: positions valid; P: motif matches queried set; D: deterministic

    // FindKnownMotifs returns, for each queried motif that occurs, all overlapping start positions
    // (sorted ascending) via exact suffix-tree matching (Gusfield 1997).

    /// <summary>A sequence with a motif list: two of its own 3-mers plus two random 3-mers.</summary>
    private static Arbitrary<(string seq, string[] motifs)> SeqMotifsArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            string seq = RandDna(rng, 20);
            var motifs = new[]
            {
                seq.Substring(rng.Next(seq.Length - 3 + 1), 3),
                seq.Substring(rng.Next(seq.Length - 3 + 1), 3),
                RandDna(rng, 3),
                RandDna(rng, 3),
            };
            return (seq, motifs);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R + P): every returned motif is one of the queried motifs, and its positions are exactly
    /// all overlapping occurrences of it in the sequence (valid, ascending, genuine substrings).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KnownMotifs_PositionsAreExactOccurrences()
    {
        return Prop.ForAll(SeqMotifsArbitrary(), input =>
        {
            var (seq, motifs) = input;
            var queried = motifs.Select(m => m.ToUpperInvariant()).ToHashSet();
            var result = GenomicAnalyzer.FindKnownMotifs(new DnaSequence(seq), motifs);

            bool ok = result.All(kv =>
            {
                if (!queried.Contains(kv.Key)) return false;
                var expected = Enumerable.Range(0, seq.Length - kv.Key.Length + 1)
                    .Where(i => seq.Substring(i, kv.Key.Length) == kv.Key).ToList();
                return kv.Value.SequenceEqual(expected) && kv.Value.Count > 0;
            });
            return ok.Label("a known-motif entry was not an exact occurrence set of a queried motif");
        });
    }

    /// <summary>
    /// INV-2 (P, positive control): a motif embedded in the sequence is reported at its position; an
    /// absent motif is omitted.
    /// </summary>
    [Test]
    [Category("Property")]
    public void KnownMotifs_PresentReported_AbsentOmitted()
    {
        var result = GenomicAnalyzer.FindKnownMotifs(new DnaSequence("GGGATCGGG"), new[] { "ATC", "TTT", "" });
        Assert.Multiple(() =>
        {
            Assert.That(result.ContainsKey("ATC"), Is.True);
            Assert.That(result["ATC"], Is.EqualTo(new[] { 3 }));
            Assert.That(result.ContainsKey("TTT"), Is.False, "absent motif omitted");
            Assert.That(result.ContainsKey(""), Is.False, "empty motif skipped");
        });
    }

    /// <summary>
    /// INV-3 (D): Known-motif search is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property KnownMotifs_IsDeterministic()
    {
        return Prop.ForAll(SeqMotifsArbitrary(), input =>
        {
            var (seq, motifs) = input;
            var a = GenomicAnalyzer.FindKnownMotifs(new DnaSequence(seq), motifs);
            var b = GenomicAnalyzer.FindKnownMotifs(new DnaSequence(seq), motifs);
            return (a.Count == b.Count && a.All(kv => b.TryGetValue(kv.Key, out var v) && v.SequenceEqual(kv.Value)))
                .Label("FindKnownMotifs must be deterministic");
        });
    }

    #endregion
}
