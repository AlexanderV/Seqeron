using FsCheck;
using FsCheck.Fluent;
using FsCheck.NUnit;

namespace Seqeron.Genomics.Tests.Properties;

/// <summary>
/// Property-based tests for sequence/genome assembly algorithms (SequenceAssembler).
/// Verifies invariants from the literature each algorithm implements.
///
/// Test Units: ASSEMBLY-CONSENSUS-001, ASSEMBLY-CORRECT-001, ASSEMBLY-COVER-001, ASSEMBLY-DBG-001, ASSEMBLY-MERGE-001, ASSEMBLY-OLC-001
/// </summary>
[TestFixture]
[Category("Property")]
[Category("Assembly")]
public class AssemblyProperties
{
    /// <summary>Generates 2..5 equal-length DNA reads (gap-free) over {A,C,G,T}.</summary>
    private static Arbitrary<string[]> AlignedReadsArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int rows = 2 + rng.Next(4);   // 2..5 reads
            int len = 3 + rng.Next(6);    // length 3..8
            const string bases = "ACGT";
            var reads = new string[rows];
            for (int r = 0; r < rows; r++)
            {
                var chars = new char[len];
                for (int i = 0; i < len; i++) chars[i] = bases[rng.Next(4)];
                reads[r] = new string(chars);
            }
            return reads;
        }).ToArbitrary();

    #region ASSEMBLY-CONSENSUS-001: R: length = longest read; P: each committed position = majority base; D: deterministic

    // ComputeConsensus follows Biopython dumb_consensus / EMBOSS cons: per column, tally non-gap
    // residues; emit the unique maximum residue when its frequency ≥ threshold, otherwise the
    // ambiguous symbol 'N'. Consensus length equals the longest read.

    /// <summary>
    /// INV-1 (R): The consensus length equals the longest input read (the full alignment width).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Consensus_Length_EqualsLongestRead()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            int longest = reads.Max(r => r.Length);
            int len = SequenceAssembler.ComputeConsensus(reads).Length;
            return (len == longest).Label($"consensus length {len} ≠ longest read {longest}");
        });
    }

    /// <summary>
    /// INV-2 (P): Every emitted character is either the ambiguous symbol or a residue actually
    /// observed in that column — the consensus never invents a base.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Consensus_EachChar_IsAmbiguousOrObserved()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            string consensus = SequenceAssembler.ComputeConsensus(reads);
            for (int pos = 0; pos < consensus.Length; pos++)
            {
                if (consensus[pos] == 'N') continue;
                var column = reads.Where(r => pos < r.Length).Select(r => char.ToUpperInvariant(r[pos])).ToHashSet();
                if (!column.Contains(consensus[pos]))
                    return false.Label($"position {pos}: '{consensus[pos]}' not present in column");
            }
            return true.Label("all columns consistent");
        });
    }

    /// <summary>
    /// INV-3 (P, majority rule): when one base occupies a strict majority of a column (count &gt; half
    /// the non-gap residues), the default-threshold (0.5) consensus must emit exactly that base.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Consensus_StrictMajorityColumn_EmitsMajorityBase()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            string consensus = SequenceAssembler.ComputeConsensus(reads);
            for (int pos = 0; pos < consensus.Length; pos++)
            {
                var column = reads.Where(r => pos < r.Length).Select(r => char.ToUpperInvariant(r[pos])).ToList();
                var grouped = column.GroupBy(c => c).Select(g => (Base: g.Key, Count: g.Count())).ToList();
                var top = grouped.OrderByDescending(g => g.Count).First();
                bool strictMajority = top.Count * 2 > column.Count;
                if (strictMajority && consensus[pos] != top.Base)
                    return false.Label($"position {pos}: majority '{top.Base}' but consensus '{consensus[pos]}'");
            }
            return true.Label("all columns consistent");
        });
    }

    /// <summary>
    /// INV-4 (D): Consensus computation is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Consensus_IsDeterministic()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
            (SequenceAssembler.ComputeConsensus(reads) == SequenceAssembler.ComputeConsensus(reads))
                .Label("ComputeConsensus must be deterministic"));
    }

    /// <summary>
    /// INV-5 (positive/boundary controls): unanimous columns reproduce the read; tied columns and
    /// the empty input behave per the dumb_consensus rule.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Consensus_GoldenAndBoundaryCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.ComputeConsensus(new[] { "ACGT", "ACGT", "ACGT" }), Is.EqualTo("ACGT"),
                "unanimous columns yield the agreed base");
            // Two reads disagree at every column → ties → all ambiguous.
            Assert.That(SequenceAssembler.ComputeConsensus(new[] { "AA", "CC" }), Is.EqualTo("NN"),
                "tied columns emit the ambiguous symbol");
            Assert.That(SequenceAssembler.ComputeConsensus(Array.Empty<string>()), Is.EqualTo(""),
                "empty input yields an empty consensus");
            // Gap symbols are skipped: column 1 has only 'A' (the dash is ignored) → committed.
            Assert.That(SequenceAssembler.ComputeConsensus(new[] { "A-", "AA" }), Is.EqualTo("AA"),
                "gaps are excluded from the column tally");
        });
    }

    #endregion

    #region ASSEMBLY-CORRECT-001: P: corrected reads keep length/count; M: more coverage → fewer errors; D: deterministic

    // ErrorCorrectReads is Musket/Quake two-sided k-mer-spectrum correction (Liu et al. 2013; Kelley
    // et al. 2010): a k-mer is trusted when its multiplicity ≥ cut-off; an untrusted position is
    // changed to the unique base that makes all covering k-mers trusted. Corrections are single-base
    // substitutions, so length and read count are preserved.

    private static int Hamming(string a, string b)
    {
        int d = 0;
        for (int i = 0; i < a.Length; i++) if (a[i] != b[i]) d++;
        return d;
    }

    /// <summary>
    /// INV-1 (P): Correction preserves the read count and every read's length (substitution-only model).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Correct_PreservesCountAndLengths()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);
            bool ok = corrected.Count == reads.Length
                      && corrected.Select((c, i) => c.Length == reads[i].Length).All(x => x);
            return ok.Label("read count or a read length changed during correction");
        });
    }

    /// <summary>
    /// INV-2 (P, no spurious corrections): error-free data (all reads identical, fully solid) is
    /// returned unchanged — correction never edits a base covered by a trusted k-mer.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Correct_ErrorFreeData_IsUnchanged()
    {
        return Prop.ForAll(DnaReadArbitrary(20), seq =>
        {
            var reads = Enumerable.Repeat(seq, 4).ToArray();
            var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 5, minKmerFrequency: 2);
            return corrected.All(c => c == seq).Label("a correction was applied to error-free data");
        });
    }

    /// <summary>
    /// INV-3 (D): Error correction is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Correct_IsDeterministic()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            var a = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);
            var b = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);
            return a.SequenceEqual(b).Label("ErrorCorrectReads must be deterministic");
        });
    }

    /// <summary>
    /// INV-4 (M, coverage → correction): a single substitution error is recovered when many correct
    /// copies make the true k-mers solid, but is left uncorrected at coverage too low for the true
    /// k-mers to be trusted — more coverage yields fewer residual errors.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Correct_HigherCoverage_RemovesMoreErrors()
    {
        // Non-repetitive (distinct 5-mers) so a true k-mer's multiplicity tracks coverage, not
        // internal periodicity — the coverage effect is then isolated.
        const string truth = "TTGACCAGTCGAATGCCTAG";
        // One erroneous copy: substitute the base at position 10 (A→C).
        char[] e = truth.ToCharArray();
        e[10] = e[10] == 'C' ? 'G' : 'C';
        string erroneous = new string(e);

        // Low coverage: 1 correct + 1 erroneous (true k-mers seen once < cut-off 2).
        var low = new[] { truth, erroneous };
        var lowCorrected = SequenceAssembler.ErrorCorrectReads(low, kmerSize: 5, minKmerFrequency: 2);
        int lowErrors = Hamming(lowCorrected[1], truth);

        // High coverage: 4 correct + 1 erroneous (true k-mers solid, error k-mers rare).
        var high = new[] { truth, truth, truth, truth, erroneous };
        var highCorrected = SequenceAssembler.ErrorCorrectReads(high, kmerSize: 5, minKmerFrequency: 2);
        int highErrors = Hamming(highCorrected[4], truth);

        Assert.Multiple(() =>
        {
            Assert.That(highErrors, Is.LessThan(lowErrors), "more coverage must reduce residual errors");
            Assert.That(highErrors, Is.EqualTo(0), "sufficient coverage must recover the true sequence");
            Assert.That(lowErrors, Is.EqualTo(1), "insufficient coverage must leave the error uncorrected");
        });
    }

    /// <summary>
    /// INV-5 (boundary): k &lt; 1 throws; reads shorter than k are returned unchanged (upper-cased).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Correct_Boundaries()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<ArgumentOutOfRangeException>(
                () => SequenceAssembler.ErrorCorrectReads(new[] { "ACGT" }, kmerSize: 0));
            var shortReads = SequenceAssembler.ErrorCorrectReads(new[] { "acg" }, kmerSize: 5, minKmerFrequency: 2);
            Assert.That(shortReads, Is.EqualTo(new[] { "ACG" }), "reads shorter than k are upper-cased and unchanged");
        });
    }

    /// <summary>Generates a single DNA read of the given length.</summary>
    private static Arbitrary<string> DnaReadArbitrary(int len) =>
        Gen.Elements('A', 'C', 'G', 'T').ArrayOf().Where(a => a.Length >= len)
            .Select(a => new string(a, 0, len)).ToArbitrary();

    #endregion

    #region ASSEMBLY-COVER-001: R: depth ≥ 0; P: Σ depth = total placed bases (mean = bases/ref length); D: deterministic

    // CalculateCoverage places each read at its best ungapped match (≥ minOverlap matches) and
    // increments per-base depth over the spanned, end-clipped interval (SAMtools depth model). The
    // mean depth is the total placed bases divided by the reference length.

    private const int CoverMinOverlap = 5;

    /// <summary>Generates a 30 bp reference and 1..5 reads that are exact substrings of it (so they place).</summary>
    private static Arbitrary<(string reference, string[] reads)> CoverageInputArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            const string bases = "ACGT";
            int refLen = 30;
            var refChars = new char[refLen];
            for (int i = 0; i < refLen; i++) refChars[i] = bases[rng.Next(4)];
            string reference = new string(refChars);

            int n = 1 + rng.Next(5);
            var reads = new string[n];
            for (int r = 0; r < n; r++)
            {
                int len = 5 + rng.Next(6);               // 5..10
                int start = rng.Next(refLen - len + 1);   // fits fully within the reference
                reads[r] = reference.Substring(start, len);
            }
            return (reference, reads);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (R): The depth array has one entry per reference base and every entry is non-negative.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Coverage_DepthsAreNonNegative()
    {
        return Prop.ForAll(DnaReadArbitrary(30), AlignedReadsArbitrary(), (reference, reads) =>
        {
            int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: CoverMinOverlap);
            return (depth.Length == reference.Length && depth.All(d => d >= 0))
                .Label("depth array length wrong or a negative depth produced");
        });
    }

    /// <summary>
    /// INV-2 (P): For reads that place fully (here, exact substrings), the total depth equals the sum
    /// of read lengths, so the mean depth equals total placed bases ÷ reference length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Coverage_TotalDepth_EqualsPlacedBases()
    {
        return Prop.ForAll(CoverageInputArbitrary(), input =>
        {
            var (reference, reads) = input;
            int[] depth = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: CoverMinOverlap);
            long totalDepth = depth.Sum();
            long placedBases = reads.Sum(r => (long)r.Length);
            double meanDepth = (double)totalDepth / reference.Length;
            return (totalDepth == placedBases
                    && Math.Abs(meanDepth - (double)placedBases / reference.Length) < 1e-9)
                .Label($"Σdepth={totalDepth} ≠ placed bases={placedBases}");
        });
    }

    /// <summary>
    /// INV-3 (D): Coverage computation is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Coverage_IsDeterministic()
    {
        return Prop.ForAll(CoverageInputArbitrary(), input =>
        {
            var (reference, reads) = input;
            var a = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: CoverMinOverlap);
            var b = SequenceAssembler.CalculateCoverage(reference, reads, minOverlap: CoverMinOverlap);
            return a.SequenceEqual(b).Label("CalculateCoverage must be deterministic");
        });
    }

    /// <summary>
    /// INV-4 (boundary): no reads → all-zero depth; a read longer than the reference cannot place.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Coverage_Boundaries()
    {
        const string reference = "ACGTACGTACGTACGT";
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.CalculateCoverage(reference, Array.Empty<string>(), minOverlap: CoverMinOverlap),
                Is.All.Zero, "no reads must give zero depth everywhere");
            var tooLong = SequenceAssembler.CalculateCoverage(reference, new[] { reference + "ACGT" }, minOverlap: CoverMinOverlap);
            Assert.That(tooLong.Sum(), Is.EqualTo(0), "a read longer than the reference cannot be placed");
        });
    }

    #endregion

    #region ASSEMBLY-DBG-001: P: reconstructs an unambiguous genome; M: larger k → ≤ branching; R: result fields consistent; D: deterministic

    // AssembleDeBruijn decomposes reads into k-mers (edges of a (k-1)-mer graph) and spells one
    // Eulerian walk per component (Compeau et al. 2011; Langmead DBG notes). Reconstruction is exact
    // when the walk is unique — i.e. when no (k-1)-mer repeats. Larger k makes (k-1)-mers more
    // specific, removing the repeat-induced branch nodes.

    private static SequenceAssembler.AssemblyResult Dbg(string[] reads, int k) =>
        SequenceAssembler.AssembleDeBruijn(reads,
            new SequenceAssembler.AssemblyParameters(KmerSize: k, MinContigLength: 1));

    private static int BranchNodeCount(string seq, int k) =>
        SequenceAssembler.BuildDeBruijnGraph(new[] { seq }, k).Values.Count(v => v.Distinct().Count() >= 2);

    /// <summary>
    /// INV-1 (P): When a sequence has no repeated (k-1)-mer its de Bruijn graph is a simple path, so
    /// the unique Eulerian walk reconstructs exactly that sequence as a single contig.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Dbg_UnambiguousSequence_IsReconstructedExactly()
    {
        return Prop.ForAll(DnaReadArbitrary(25), seq =>
        {
            const int k = 10; // (k-1)=9-mer nodes
            var km1mers = Enumerable.Range(0, seq.Length - (k - 1) + 1).Select(i => seq.Substring(i, k - 1));
            bool unambiguous = km1mers.Distinct().Count() == seq.Length - (k - 1) + 1;
            if (!unambiguous)
                return true.Label("skip: repeated (k-1)-mer");

            var contigs = Dbg(new[] { seq }, k).Contigs;
            return (contigs.Count == 1 && contigs[0] == seq)
                .Label($"reconstruction failed: got [{string.Join(",", contigs)}], expected {seq}");
        });
    }

    /// <summary>
    /// INV-2 (M): Increasing k removes repeat-induced branching. A sequence with a repeated 4-mer
    /// branches at k=5 (4-mer nodes) but the branch disappears once the (k-1)-mer is longer than the
    /// repeat (k=9, 8-mer nodes).
    /// </summary>
    [Test]
    [Category("Property")]
    public void Dbg_LargerK_ReducesBranching()
    {
        const string seq = "AAAATTTTAAAACCCC"; // 4-mer "AAAA" repeats at offsets 0 and 8
        int branchSmallK = BranchNodeCount(seq, 5);
        int branchLargeK = BranchNodeCount(seq, 9);

        Assert.Multiple(() =>
        {
            Assert.That(branchSmallK, Is.GreaterThan(0), "the repeated 4-mer must branch at k=5");
            Assert.That(branchLargeK, Is.LessThan(branchSmallK), "larger k must not increase branching");
            Assert.That(branchLargeK, Is.EqualTo(0), "an 8-mer node no longer captures the 4-mer repeat");
        });
    }

    /// <summary>
    /// INV-3 (R): Reported assembly statistics are internally consistent — total length equals the
    /// sum of contig lengths and the longest-contig length equals the maximum contig length.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Dbg_ResultFields_AreConsistent()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            var result = Dbg(reads, 4);
            int sum = result.Contigs.Sum(c => c.Length);
            int max = result.Contigs.Count == 0 ? 0 : result.Contigs.Max(c => c.Length);
            return (result.TotalLength == sum && result.LongestContig == max
                    && result.Contigs.All(c => c.Length >= 1))
                .Label($"inconsistent stats: TotalLength={result.TotalLength} vs {sum}, LongestContig={result.LongestContig} vs {max}");
        });
    }

    /// <summary>
    /// INV-4 (D): De Bruijn assembly is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Dbg_IsDeterministic()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
            Dbg(reads, 4).Contigs.SequenceEqual(Dbg(reads, 4).Contigs)
                .Label("AssembleDeBruijn must be deterministic"));
    }

    /// <summary>
    /// INV-5 (boundary): empty input yields an empty assembly; k &lt; 2 is rejected by the graph builder.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Dbg_Boundaries()
    {
        Assert.Multiple(() =>
        {
            Assert.That(Dbg(Array.Empty<string>(), 4).Contigs, Is.Empty, "no reads → no contigs");
            Assert.Throws<ArgumentOutOfRangeException>(() => SequenceAssembler.BuildDeBruijnGraph(new[] { "ACGT" }, 1));
        });
    }

    #endregion

    #region ASSEMBLY-MERGE-001: P: merged length ≥ longer contig; P: result is a superstring (overlap collapsed once); D: deterministic

    // MergeContigs collapses a known suffix/prefix overlap of length l so the shared region appears
    // once: merged = contig1 + contig2[l:], length |c1|+|c2|−l (Langmead SCS notes). An invalid
    // overlap (≤ 0 or > min length) falls back to plain concatenation.

    private static string RandDna(Random rng, int len)
    {
        const string bases = "ACGT";
        var chars = new char[len];
        for (int i = 0; i < len; i++) chars[i] = bases[rng.Next(4)];
        return new string(chars);
    }

    /// <summary>Two random contigs plus an arbitrary (possibly invalid) overlap length 0..15.</summary>
    private static Arbitrary<(string c1, string c2, int overlap)> MergeArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            return (RandDna(rng, 3 + rng.Next(8)), RandDna(rng, 3 + rng.Next(8)), rng.Next(0, 16));
        }).ToArbitrary();

    /// <summary>A pair sharing a genuine overlap g: c1 = prefix+g, c2 = g+suffix, overlap = |g|.</summary>
    private static Arbitrary<(string c1, string c2, int overlap)> GenuineOverlapArbitrary() =>
        Gen.Choose(0, int.MaxValue).Select(seed =>
        {
            var rng = new Random(seed);
            int l = 1 + rng.Next(5);
            string g = RandDna(rng, l);
            string c1 = RandDna(rng, rng.Next(6)) + g;
            string c2 = g + RandDna(rng, rng.Next(6));
            return (c1, c2, l);
        }).ToArbitrary();

    /// <summary>
    /// INV-1 (P): The merged length is |c1|+|c2|−l for a valid overlap (else |c1|+|c2|), and is never
    /// shorter than the longer input contig.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Merge_Length_IsCorrectAndAtLeastLongerContig()
    {
        return Prop.ForAll(MergeArbitrary(), t =>
        {
            var (c1, c2, l) = t;
            string merged = SequenceAssembler.MergeContigs(c1, c2, l);
            bool validOverlap = l > 0 && l <= Math.Min(c1.Length, c2.Length);
            int expected = validOverlap ? c1.Length + c2.Length - l : c1.Length + c2.Length;
            return (merged.Length == expected && merged.Length >= Math.Max(c1.Length, c2.Length))
                .Label($"merged length {merged.Length} ≠ expected {expected}");
        });
    }

    /// <summary>
    /// INV-2 (P, superstring): for a genuine suffix/prefix overlap the merged contig contains both
    /// inputs — c1 is its prefix and c2 is its suffix — with the shared region collapsed once.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Merge_GenuineOverlap_YieldsSuperstring()
    {
        return Prop.ForAll(GenuineOverlapArbitrary(), t =>
        {
            var (c1, c2, l) = t;
            string merged = SequenceAssembler.MergeContigs(c1, c2, l);
            return (merged.StartsWith(c1, StringComparison.Ordinal)
                    && merged.EndsWith(c2, StringComparison.Ordinal)
                    && merged.Length == c1.Length + c2.Length - l)
                .Label($"merge of overlap {l}: '{merged}' is not a superstring of '{c1}' and '{c2}'");
        });
    }

    /// <summary>
    /// INV-3 (D): Merging is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Merge_IsDeterministic()
    {
        return Prop.ForAll(MergeArbitrary(), t =>
        {
            var (c1, c2, l) = t;
            return (SequenceAssembler.MergeContigs(c1, c2, l) == SequenceAssembler.MergeContigs(c1, c2, l))
                .Label("MergeContigs must be deterministic");
        });
    }

    /// <summary>
    /// INV-4 (boundary): non-positive and over-long overlaps fall back to concatenation; the
    /// canonical Langmead trace BAA + AAB at overlap 2 collapses to BAAB.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Merge_BoundaryAndGoldenCases()
    {
        Assert.Multiple(() =>
        {
            Assert.That(SequenceAssembler.MergeContigs("BAA", "AAB", 2), Is.EqualTo("BAAB"), "overlap collapsed once");
            Assert.That(SequenceAssembler.MergeContigs("ACGT", "TGCA", 0), Is.EqualTo("ACGTTGCA"), "no overlap → concat");
            Assert.That(SequenceAssembler.MergeContigs("ACGT", "TGCA", -1), Is.EqualTo("ACGTTGCA"), "negative overlap → concat");
            Assert.That(SequenceAssembler.MergeContigs("AC", "GT", 5), Is.EqualTo("ACGT"), "overlap > shorter contig → concat");
        });
    }

    #endregion

    #region ASSEMBLY-OLC-001: R: every read ends up in a contig; P: reported overlaps ≥ minOverlap; P: clean tiling reconstructs the genome; D: deterministic

    // AssembleOLC follows Overlap-Layout-Consensus (Compeau et al. 2011; Langmead OLC notes): find
    // suffix/prefix overlaps ≥ minOverlap, greedily chain reads by best overlap, emit the merged
    // superstring of each chain.

    private static SequenceAssembler.AssemblyResult Olc(string[] reads, int minOverlap) =>
        SequenceAssembler.AssembleOLC(reads,
            new SequenceAssembler.AssemblyParameters(MinOverlap: minOverlap, MinIdentity: 0.9, MinContigLength: 1));

    /// <summary>
    /// INV-1 (P, reconstruction): a clean overlapping tiling of a non-repetitive genome is
    /// reassembled into the original sequence as a single contig.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Olc_CleanTiling_ReconstructsGenome()
    {
        const string genome = "TTGACCAGTCGAATGCCTAGGCATTACGGT"; // 30 bp, non-repetitive
        // Windows of length 12, step 6 → consecutive reads overlap by 6.
        var reads = new[]
        {
            genome.Substring(0, 12),
            genome.Substring(6, 12),
            genome.Substring(12, 12),
            genome.Substring(18, 12),
        };
        var contigs = Olc(reads, minOverlap: 5).Contigs;

        Assert.That(contigs, Has.Count.EqualTo(1));
        Assert.That(contigs[0], Is.EqualTo(genome), "clean tiling must reconstruct the genome");
    }

    /// <summary>
    /// INV-2 (P): Every reported overlap is at least minOverlap and at most the shorter read length
    /// (a suffix/prefix overlap cannot exceed either read).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Olc_ReportedOverlaps_MeetMinOverlap()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            const int minOverlap = 3;
            var overlaps = SequenceAssembler.FindAllOverlaps(reads, minOverlap, 0.9);
            bool ok = overlaps.All(o =>
                o.OverlapLength >= minOverlap &&
                o.OverlapLength <= Math.Min(reads[o.ReadIndex1].Length, reads[o.ReadIndex2].Length));
            return ok.Label("an overlap was below minOverlap or longer than a read");
        });
    }

    /// <summary>
    /// INV-3 (R): With reads present and no length filter, the assembly produces at least one contig
    /// and every contig has positive length (each read is placed in some contig).
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Olc_ProducesAtLeastOneContig()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
        {
            var contigs = Olc(reads, minOverlap: 3).Contigs;
            return (contigs.Count >= 1 && contigs.All(c => c.Length >= 1))
                .Label($"expected ≥1 non-empty contig, got {contigs.Count}");
        });
    }

    /// <summary>
    /// INV-4 (D): OLC assembly is deterministic.
    /// </summary>
    [FsCheck.NUnit.Property]
    public Property Olc_IsDeterministic()
    {
        return Prop.ForAll(AlignedReadsArbitrary(), reads =>
            Olc(reads, 3).Contigs.SequenceEqual(Olc(reads, 3).Contigs)
                .Label("AssembleOLC must be deterministic"));
    }

    /// <summary>
    /// INV-5 (boundary): empty input yields an empty assembly.
    /// </summary>
    [Test]
    [Category("Property")]
    public void Olc_EmptyInput_YieldsNoContigs()
    {
        Assert.That(Olc(Array.Empty<string>(), 3).Contigs, Is.Empty);
    }

    #endregion
}
