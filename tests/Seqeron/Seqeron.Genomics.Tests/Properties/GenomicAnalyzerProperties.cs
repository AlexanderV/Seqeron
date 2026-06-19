using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for the high-level GenomicAnalyzer façade: common regions, known motifs,
/// ORFs, repeats, similarity, and tandem repeats.
///
/// Test Units: GENOMIC-COMMON-001, GENOMIC-MOTIFS-001, GENOMIC-ORF-001, GENOMIC-REPEAT-001, GENOMIC-SIMILARITY-001, GENOMIC-TANDEM-001
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

    #region GENOMIC-ORF-001: P: starts ATG, ends stop, no internal stop; R: length %3 = 0, positions valid; D: deterministic

    // FindOpenReadingFrames reports each ATG→(first in-frame stop) span across all six frames
    // (Rosalind ORF semantics). The reported sequence spans the start through the stop inclusive.

    private static readonly string[] StopCodons = { "TAA", "TAG", "TGA" };

    private static bool IsValidOrf(string orf, int minLength)
    {
        if (orf.Length < minLength || orf.Length % 3 != 0) return false;
        if (!orf.StartsWith("ATG", StringComparison.Ordinal)) return false;
        if (!StopCodons.Contains(orf[^3..])) return false;
        // No in-frame stop before the final (terminating) codon.
        for (int i = 0; i <= orf.Length - 6; i += 3)
            if (StopCodons.Contains(orf.Substring(i, 3))) return false;
        return true;
    }

    /// <summary>
    /// INV-1 (P + R): every ORF starts with ATG, ends with an in-frame stop with no earlier in-frame
    /// stop, has length divisible by 3 and ≥ minLength, and valid coordinates on its strand.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Orfs_AreWellFormed()
    {
        return Prop.ForAll(SeqArbitrary(30), seq =>
        {
            const int minLength = 6;
            int len = seq.Length;
            var orfs = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(seq), minLength).ToList();
            bool ok = orfs.All(o =>
                IsValidOrf(o.Sequence, minLength) &&
                o.Position >= 0 && o.Position + o.Length <= len &&
                o.Frame is >= 1 and <= 3);
            return ok.Label("a reported ORF was malformed");
        });
    }

    /// <summary>
    /// INV-2 (P, positive control): an embedded forward ORF is found exactly.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Orfs_EmbeddedOrf_IsFound()
    {
        var orfs = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence("ATGAAATAA"), minLength: 9).ToList();
        Assert.That(orfs.Any(o => o.Sequence == "ATGAAATAA" && !o.IsReverseComplement && o.Position == 0), Is.True,
            "the ATG-AAA-TAA ORF must be reported on the forward strand");
    }

    /// <summary>
    /// INV-3 (D): ORF detection is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Orfs_IsDeterministic()
    {
        return Prop.ForAll(SeqArbitrary(30), seq =>
        {
            var a = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(seq), 6)
                .Select(o => (o.Sequence, o.Position, o.Frame, o.IsReverseComplement)).ToList();
            var b = GenomicAnalyzer.FindOpenReadingFrames(new DnaSequence(seq), 6)
                .Select(o => (o.Sequence, o.Position, o.Frame, o.IsReverseComplement)).ToList();
            return a.SequenceEqual(b).Label("FindOpenReadingFrames must be deterministic");
        });
    }

    /// <summary>Generates a DNA sequence of at least <paramref name="minLen"/> bases.</summary>
    private static Arbitrary<string> SeqArbitrary(int minLen) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= minLen)
            .Select(a => new string(a)).ToArbitrary();

    #endregion

    #region GENOMIC-REPEAT-001: R: positions valid; M: lower minLen → ≥ repeats; D: deterministic

    // FindRepeats returns every distinct substring of length ≥ minLength occurring at least twice
    // (suffix-tree / LCP repeats; CMU 15-451).

    /// <summary>
    /// INV-1 (R): every repeat has length ≥ minLength, occurs ≥ 2 times at distinct ascending valid
    /// positions whose substrings equal the repeat.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Repeats_AreValid()
    {
        return Prop.ForAll(SeqArbitrary(15), seq =>
        {
            const int minLength = 2;
            var repeats = GenomicAnalyzer.FindRepeats(new DnaSequence(seq), minLength).ToList();
            bool ok = repeats.All(r =>
                r.Length >= minLength && r.Count >= 2 &&
                r.Positions.SequenceEqual(r.Positions.OrderBy(p => p)) &&
                r.Positions.Distinct().Count() == r.Positions.Count &&
                r.Positions.All(p => p >= 0 && p + r.Length <= seq.Length && seq.Substring(p, r.Length) == r.Sequence));
            return ok.Label("a repeat was invalid");
        });
    }

    /// <summary>
    /// INV-2 (M): a lower minimum length reports at least as many repeats — minLen 3 ⊆ minLen 2.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Repeats_LowerMinLength_IsSuperset()
    {
        return Prop.ForAll(SeqArbitrary(15), seq =>
        {
            var loose = GenomicAnalyzer.FindRepeats(new DnaSequence(seq), 2).Select(r => r.Sequence).ToHashSet();
            var strict = GenomicAnalyzer.FindRepeats(new DnaSequence(seq), 3).Select(r => r.Sequence).ToHashSet();
            return strict.IsSubsetOf(loose).Label("minLen 3 repeats not ⊆ minLen 2 repeats");
        });
    }

    /// <summary>
    /// INV-3 (D + positive control): repeats are deterministic; a tandem "ATGATGATG" reports the ATG
    /// repeat at positions {0,3,6}.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Repeats_DeterministicAndGolden()
    {
        var dna = new DnaSequence("ATGATGATG");
        var a = GenomicAnalyzer.FindRepeats(dna, 3).Select(r => r.Sequence).OrderBy(s => s).ToList();
        var b = GenomicAnalyzer.FindRepeats(dna, 3).Select(r => r.Sequence).OrderBy(s => s).ToList();
        var atg = GenomicAnalyzer.FindRepeats(dna, 3).FirstOrDefault(r => r.Sequence == "ATG");
        Assert.Multiple(() =>
        {
            Assert.That(b, Is.EqualTo(a), "deterministic");
            Assert.That(atg.Sequence, Is.EqualTo("ATG"));
            Assert.That(atg.Positions, Is.EqualTo(new[] { 0, 3, 6 }));
        });
    }

    #endregion

    #region GENOMIC-SIMILARITY-001: R: similarity ∈ [0,100]; S: sym; I: sim(x,x)=100; D: deterministic

    // CalculateSimilarity is the k-mer Jaccard index reported as a percentage in [0,100]
    // (Jaccard 1901; Ondov et al. Mash 2016). NOTE: this implementation returns a PERCENTAGE,
    // so the self-identity value is 100 (not 1) — the checklist's [0,1]/sim(x,x)=1 is in fractions.

    /// <summary>INV-1 (R): similarity is a percentage in [0,100].</summary>
    [FsCheck.NUnit.Property]
    public Property Similarity_InZeroToHundred()
    {
        return Prop.ForAll(SeqArbitrary(8), SeqArbitrary(8), (s1, s2) =>
        {
            double sim = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(s1), new DnaSequence(s2));
            return (sim is >= 0.0 and <= 100.0).Label($"similarity {sim} outside [0,100]");
        });
    }

    /// <summary>INV-2 (S): similarity is symmetric.</summary>
    [FsCheck.NUnit.Property]
    public Property Similarity_IsSymmetric()
    {
        return Prop.ForAll(SeqArbitrary(8), SeqArbitrary(8), (s1, s2) =>
        {
            double ab = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(s1), new DnaSequence(s2));
            double ba = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(s2), new DnaSequence(s1));
            return (Math.Abs(ab - ba) < 1e-9).Label($"sim(a,b)={ab} ≠ sim(b,a)={ba}");
        });
    }

    /// <summary>INV-3 (I): a sequence (with at least one k-mer) is 100% similar to itself.</summary>
    [FsCheck.NUnit.Property]
    public Property Similarity_SelfIsHundred()
    {
        return Prop.ForAll(SeqArbitrary(8), seq =>
        {
            double sim = GenomicAnalyzer.CalculateSimilarity(new DnaSequence(seq), new DnaSequence(seq));
            return (Math.Abs(sim - 100.0) < 1e-9).Label($"sim(x,x)={sim}, expected 100");
        });
    }

    /// <summary>INV-4 (D + boundary): deterministic; invalid k-mer size rejected.</summary>
    [Test]
    [Category("Property")]
    public void Similarity_DeterministicAndBoundary()
    {
        var a = GenomicAnalyzer.CalculateSimilarity(new DnaSequence("ACGTACGT"), new DnaSequence("ACGTTTGG"));
        var b = GenomicAnalyzer.CalculateSimilarity(new DnaSequence("ACGTACGT"), new DnaSequence("ACGTTTGG"));
        Assert.Multiple(() =>
        {
            Assert.That(b, Is.EqualTo(a), "deterministic");
            Assert.Throws<ArgumentOutOfRangeException>(
                () => GenomicAnalyzer.CalculateSimilarity(new DnaSequence("ACGT"), new DnaSequence("ACGT"), 0));
        });
    }

    #endregion

    #region GENOMIC-TANDEM-001: R: repetitions ≥ minReps; P: unit repeated contiguously; D: deterministic

    // FindTandemRepeats locates consecutive repeating units (e.g. ATGATGATG = ATG × 3).

    /// <summary>
    /// INV-1 (R + P): every tandem repeat has Repetitions ≥ minReps, unit length ≥ minUnit, valid
    /// position, and the unit appears contiguously Repetitions times at that position.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Tandem_UnitsRepeatContiguously()
    {
        return Prop.ForAll(SeqArbitrary(12), seq =>
        {
            const int minUnit = 2;
            const int minReps = 2;
            var tandems = GenomicAnalyzer.FindTandemRepeats(new DnaSequence(seq), minUnit, minReps).ToList();
            bool ok = tandems.All(t =>
                t.Repetitions >= minReps && t.Unit.Length >= minUnit &&
                t.Position >= 0 && t.Position + t.TotalLength <= seq.Length &&
                seq.Substring(t.Position, t.TotalLength) == string.Concat(Enumerable.Repeat(t.Unit, t.Repetitions)));
            return ok.Label("a tandem repeat unit did not repeat contiguously");
        });
    }

    /// <summary>
    /// INV-2 (D + positive control): tandem detection is deterministic; "ATGATGATG" yields the ATG×3
    /// tandem at position 0.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Tandem_DeterministicAndGolden()
    {
        var dna = new DnaSequence("ATGATGATG");
        var a = GenomicAnalyzer.FindTandemRepeats(dna, 3, 2).Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
        var b = GenomicAnalyzer.FindTandemRepeats(dna, 3, 2).Select(t => (t.Unit, t.Position, t.Repetitions)).ToList();
        Assert.Multiple(() =>
        {
            Assert.That(b, Is.EqualTo(a), "deterministic");
            Assert.That(a.Any(t => t.Unit == "ATG" && t.Position == 0 && t.Repetitions == 3), Is.True,
                "ATG×3 tandem expected");
        });
    }

    #endregion
}
