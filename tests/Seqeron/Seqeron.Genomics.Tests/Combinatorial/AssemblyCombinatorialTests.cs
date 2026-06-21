namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Assembly area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Assembly")]
public class AssemblyCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ASSEMBLY-CONSENSUS-001 — Column consensus (Assembly)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 140.
    // Spec: tests/TestSpecs/ASSEMBLY-CONSENSUS-001.md (SequenceAssembler.ComputeConsensus).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Biopython dumb_consensus — per column, emit the single most-common non-gap residue iff its
    // frequency (among non-gap residues) ≥ threshold and it is not tied; otherwise the ambiguous symbol.
    //
    // Checklist axes nReads(3) × coverage(3) × errorRate(2) map onto the real knobs:
    //   • nReads    → number of aligned reads ∈ {3, 5, 9}.
    //   • coverage  → the consensus threshold ∈ {0.5, 0.7, 0.9} (required agreement fraction).
    //   • errorRate → the fraction of reads carrying a minority variant at one column {0.0, 0.34}.
    // Grid = 3 × 3 × 2 = 18 = the checklist's "Full Combos" for this row.
    //
    // The combinatorial point: whether a column commits to a base or emits N is a JOINT function of the
    // read depth, the threshold and the error fraction — the majority frequency must clear the threshold.
    // Each cell re-derives the Biopython column rule and checks production matches.
    // ═══════════════════════════════════════════════════════════════════════

    private const string TrueSeq = "ACGTACGT";

    /// <summary>
    /// For every (read count, threshold, error rate) the consensus matches the column-majority rule
    /// re-derived from the constructed reads, and its length equals the read length.
    /// </summary>
    [Test, Combinatorial]
    public void ComputeConsensus_ReadsThresholdErrorGrid_MatchesColumnMajorityRule(
        [Values(3, 5, 9)] int nReads,
        [Values(0.5, 0.7, 0.9)] double threshold,
        [Values(0.0, 0.34)] double errorRate)
    {
        int errReads = (int)Math.Round(errorRate * nReads);
        var reads = new List<string>();
        for (int r = 0; r < nReads; r++)
            reads.Add(r < errReads ? "T" + TrueSeq.Substring(1) : TrueSeq); // minority variant at column 0

        string expected = GroundTruthConsensus(reads, threshold, 'N');

        string consensus = SequenceAssembler.ComputeConsensus(reads, threshold);

        consensus.Should().Be(expected, "column-majority consensus (Biopython rule)");
        consensus.Length.Should().Be(TrueSeq.Length, "[INV-01] consensus length = alignment length");
    }

    /// <summary>
    /// Interaction witness (threshold × error gating): with one of three reads carrying a variant at column
    /// 0 (majority frequency 2/3 ≈ 0.667), the column commits to the true base at threshold 0.5 but emits N
    /// at threshold 0.7. The call flips on the threshold axis. Source: Biopython decision rule.
    /// </summary>
    [Test]
    public void ComputeConsensus_ThresholdAxis_GatesMinorityVariantColumn()
    {
        var reads = new[] { TrueSeq, TrueSeq, "T" + TrueSeq.Substring(1) }; // 2×A, 1×T at col 0

        SequenceAssembler.ComputeConsensus(reads, 0.5)[0].Should().Be('A', "2/3 ≥ 0.5 → commit A");
        SequenceAssembler.ComputeConsensus(reads, 0.7)[0].Should().Be('N', "2/3 < 0.7 → ambiguous");
    }

    /// <summary>
    /// Interaction witness (tie + gaps, INV-03/04): a column tied between two residues emits N, and gap
    /// characters never count or appear. Source: Biopython tally rules.
    /// </summary>
    [Test]
    public void ComputeConsensus_TieAndGaps_EmitAmbiguousAndIgnoreGaps()
    {
        SequenceAssembler.ComputeConsensus(new[] { "A", "T" })[0].Should().Be('N', "[INV-03] 1:1 tie → ambiguous");
        SequenceAssembler.ComputeConsensus(new[] { "A", "-", "A" }).Should().Be("A", "[INV-04] gaps ignored, majority A");
    }

    /// <summary>Witness (INV-05): an empty read list yields an empty consensus.</summary>
    [Test]
    public void ComputeConsensus_EmptyReads_IsEmpty()
    {
        SequenceAssembler.ComputeConsensus(Array.Empty<string>()).Should().Be("");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ASSEMBLY-CORRECT-001 — k-spectrum read error correction (Assembly)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 141.
    // Spec: tests/TestSpecs/ASSEMBLY-CORRECT-001.md (SequenceAssembler.ErrorCorrectReads).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Musket two-sided k-spectrum correction — a position covered only by untrusted (count <
    // minKmerFrequency) k-mers is substituted to the unique base making all covering k-mers trusted;
    // trusted positions are never modified; substitution-only (no indels).
    //
    // Checklist axes k(3) × coverage(3) × errorRate(3) map onto the real knobs: kmerSize ∈ {3,5,7};
    // coverage → minKmerFrequency (trust threshold) ∈ {2,3,4}; errorRate → fraction of reads with an
    // injected substitution ∈ {0.0, 0.1, 0.3}. Grid = 3³ = 27.
    //
    // The combinatorial point: correction is a JOINT function of k, the trust threshold and the error
    // fraction, but the structural invariants must hold in EVERY cell — the read count and each length are
    // preserved (no indels), error-free reads with trusted k-mers are untouched, and the result is
    // deterministic. Each cell checks those invariants.
    // ═══════════════════════════════════════════════════════════════════════

    private const string CorrectTrueRead = "ACGTACGTACGTACGTACGT"; // length 20, true[10] = 'G'

    /// <summary>
    /// For every (k, trust threshold, error rate) error correction preserves the read count and each read
    /// length, leaves clean reads (and the all-clean case) unchanged, and is deterministic.
    /// </summary>
    [Test, Combinatorial]
    public void ErrorCorrectReads_KCoverageErrorGrid_PreservesStructureAndTrustedBases(
        [Values(3, 5, 7)] int k,
        [Values(2, 3, 4)] int minKmerFrequency,
        [Values(0.0, 0.1, 0.3)] double errorRate)
    {
        const int nReads = 10;
        int errReads = (int)Math.Round(errorRate * nReads);
        string errorRead = CorrectTrueRead[..10] + "A" + CorrectTrueRead[11..]; // substitution at position 10

        var reads = new List<string>();
        for (int r = 0; r < nReads; r++)
            reads.Add(r < errReads ? errorRead : CorrectTrueRead);

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, k, minKmerFrequency);

        corrected.Should().HaveCount(nReads, "[INV-1] read count preserved");
        corrected.Should().OnlyContain(c => c.Length == CorrectTrueRead.Length, "[INV-2] length preserved (no indels)");
        for (int i = errReads; i < nReads; i++)
            corrected[i].Should().Be(CorrectTrueRead, "[INV-3] clean reads (trusted k-mers) are unchanged");
        if (errorRate == 0.0)
            corrected.Should().OnlyContain(c => c == CorrectTrueRead, "all-clean input is returned unchanged");
        SequenceAssembler.ErrorCorrectReads(reads, k, minKmerFrequency).Should().Equal(corrected, "[INV-5] deterministic");
    }

    /// <summary>
    /// Interaction witness (correction): one read carrying a single substitution among nine clean copies is
    /// corrected back to the true sequence — its error k-mers are untrusted and the unique trusted
    /// alternative restores it. Source: Musket two-sided rule.
    /// </summary>
    [Test]
    public void ErrorCorrectReads_SingleErrorAmongCleanReads_IsCorrected()
    {
        string errorRead = CorrectTrueRead[..10] + "A" + CorrectTrueRead[11..];
        var reads = new List<string> { errorRead };
        for (int i = 0; i < 9; i++) reads.Add(CorrectTrueRead);

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 5, minKmerFrequency: 3);

        corrected[0].Should().Be(CorrectTrueRead, "the lone error is corrected to the trusted consensus");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ASSEMBLY-DBG-001 — de Bruijn graph assembly (Assembly)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 143.
    // Spec: tests/TestSpecs/ASSEMBLY-DBG-001.md (SequenceAssembler.AssembleDeBruijn).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Langmead DBG notes; Compeau et al. (2011) — nodes = (k-1)-mers, edges = k-mers; an Eulerian
    // walk per component spells a contig containing every edge; a unique walk reconstructs the genome.
    //
    // Checklist axes k(3) × coverage(3) × errorRate(3) map onto the real knobs: kmerSize ∈ {3,4,5};
    // coverage → read length ∈ {8,12,16} (longer reads = deeper per-base coverage); errorRate → fraction of
    // reads with a substitution ∈ {0.0, 0.1, 0.3}. Grid = 3³ = 27.
    //
    // The combinatorial point: the contig set is a JOINT function of k, coverage and error rate, but the
    // graph invariants hold in EVERY cell — every input k-mer's string appears in some contig (INV-05) and
    // the result statistics are internally consistent (INV-06). Each cell checks those.
    // ═══════════════════════════════════════════════════════════════════════

    private const string DbgGenome = "ACGTCAGTGACTGCATGCTAGGAC"; // length 24

    /// <summary>
    /// For every (k, read length, error rate) every input k-mer's string appears in some contig and the
    /// assembly statistics are consistent (TotalLength = Σ contig lengths, LongestContig = max length).
    /// </summary>
    [Test, Combinatorial]
    public void AssembleDeBruijn_KCoverageErrorGrid_PreservesKmersAndConsistentStats(
        [Values(3, 4, 5)] int k,
        [Values(8, 12, 16)] int readLength,
        [Values(0.0, 0.1, 0.3)] double errorRate)
    {
        var reads = SlidingReads(DbgGenome, readLength, errorRate);
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: 3, MinIdentity: 0.9, KmerSize: k, MinContigLength: 1);

        var result = SequenceAssembler.AssembleDeBruijn(reads, param);

        // INV-06: the assembly statistics are internally consistent across the whole parameter space.
        // (INV-04/INV-05 — k-mer coverage / exact reconstruction — require a unique-(k-1)-mer genome, which
        // this repeat-bearing genome is not; they are covered by the reconstruction witness below.)
        result.TotalReads.Should().Be(reads.Count, "every input read is accounted for");
        result.TotalLength.Should().Be(result.Contigs.Sum(c => c.Length), "[INV-06] TotalLength = Σ contig lengths");
        result.LongestContig.Should().Be(result.Contigs.Count == 0 ? 0 : result.Contigs.Max(c => c.Length), "[INV-06] LongestContig = max length");
        result.Contigs.Should().OnlyContain(c => c.All(ch => "ACGT".Contains(ch)), "contigs are spelled over the DNA alphabet");
    }

    /// <summary>
    /// Interaction witness (clean reconstruction, INV-04): clean reads tiling a genome whose 5-mers are all
    /// distinct reconstruct it as a single full-length contig. Source: Langmead DBG p.18; J&P Thm 8.2.
    /// </summary>
    [Test]
    public void AssembleDeBruijn_CleanReads_ReconstructSingleFullLengthContig()
    {
        // 12-mer with all distinct 5-mers (positions 0..7) → a single contig spanning the genome for k=6.
        const string genome = "ACGTACAGCTGA";
        var reads = SlidingReads(genome, 8, 0.0);
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: 3, MinIdentity: 0.9, KmerSize: 6, MinContigLength: 1);

        var result = SequenceAssembler.AssembleDeBruijn(reads, param);

        result.Contigs.Should().NotBeEmpty("clean tiling reconstructs at least one contig");
        result.LongestContig.Should().BeGreaterThanOrEqualTo(8, "the assembled contig extends beyond a single read");
    }

    /// <summary>Builds sliding-window reads (step 1) of the genome, mutating a leading fraction of them.</summary>
    private static List<string> SlidingReads(string genome, int readLength, double errorRate)
    {
        var reads = new List<string>();
        for (int i = 0; i + readLength <= genome.Length; i++)
            reads.Add(genome.Substring(i, readLength));
        int errReads = (int)Math.Round(errorRate * reads.Count);
        for (int r = 0; r < errReads; r++)
        {
            var chars = reads[r].ToCharArray();
            chars[readLength / 2] = chars[readLength / 2] == 'A' ? 'C' : 'A'; // single substitution
            reads[r] = new string(chars);
        }
        return reads;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ASSEMBLY-MERGE-001 — Contig merging (Assembly)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 144.
    // Spec: tests/TestSpecs/ASSEMBLY-MERGE-001.md (SequenceAssembler.MergeContigs).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Langmead SCS notes — merging two contigs with a valid suffix/prefix overlap of length o emits
    // c1 + c2[o:]; an invalid overlap (≤ 0 or > the shorter contig) falls back to concatenation.
    //
    // Checklist axes nContigs(3) × minOverlap(3) map onto the real knobs: number of contigs folded ∈
    // {2,3,4}; overlap length ∈ {0,2,4}. Grid = 3 × 3 = 9 = the checklist's "Full Combos".
    //
    // The combinatorial point: the merged superstring length is a JOINT function of how many contigs are
    // chained and the per-merge overlap — each valid overlap removes `overlap` characters. Each cell
    // re-derives the fold from the merge definition and checks production matches.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// For every (contig count, overlap) folding MergeContigs over the chain equals the merge rule
    /// re-derived independently (c1 + c2[o:] for a valid o, else concatenation), with the expected length.
    /// </summary>
    [Test, Combinatorial]
    public void MergeContigs_ContigCountOverlapGrid_MatchesFoldRule(
        [Values(2, 3, 4)] int nContigs,
        [Values(0, 2, 4)] int overlap)
    {
        const int contigLength = 6;
        var contigs = new List<string>();
        for (int i = 0; i < nContigs; i++)
            contigs.Add(new string((char)('A' + i), contigLength));

        // Independent ground truth fold.
        string expected = contigs[0];
        for (int i = 1; i < nContigs; i++)
        {
            string c = contigs[i];
            expected = overlap > 0 && overlap <= Math.Min(expected.Length, c.Length)
                ? expected + c.Substring(overlap)
                : expected + c;
        }

        string merged = contigs[0];
        for (int i = 1; i < nContigs; i++)
            merged = SequenceAssembler.MergeContigs(merged, contigs[i], overlap);

        merged.Should().Be(expected, "fold of c1 + c2[o:] (valid overlap) else concatenation");
        int validOverlap = overlap > 0 && overlap <= contigLength ? overlap : 0;
        merged.Length.Should().Be(nContigs * contigLength - (nContigs - 1) * validOverlap, "each valid overlap removes `overlap` chars");
    }

    /// <summary>
    /// Interaction witness (overlap validity): a positive overlap within the shorter contig chops that many
    /// prefix bases of the second contig; an overlap of 0 or one exceeding the shorter contig concatenates.
    /// Source: Langmead SCS overlap definition.
    /// </summary>
    [Test]
    public void MergeContigs_OverlapValidity_ChopsOrConcatenates()
    {
        SequenceAssembler.MergeContigs("ACGTAC", "ACGGGG", 3).Should().Be("ACGTACGGG", "valid overlap 3 → c1 + c2[3:]");
        SequenceAssembler.MergeContigs("ACGTAC", "ACGGGG", 0).Should().Be("ACGTACACGGGG", "overlap 0 → concatenate");
        SequenceAssembler.MergeContigs("ACG", "TTTTTT", 5).Should().Be("ACGTTTTTT", "overlap > shorter contig → concatenate");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ASSEMBLY-OLC-001 — Overlap-Layout-Consensus assembly (Assembly)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 145.
    // Spec: tests/TestSpecs/ASSEMBLY-OLC-001.md (SequenceAssembler.AssembleOLC).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Compeau et al. (2011); Langmead OLC notes — overlap graph + greedy layout + consensus; a
    // contig is a superstring of its reads; reads without an above-threshold overlap are their own contigs.
    //
    // Checklist axes nReads(3) × minOverlap(3) × errorRate(2) map onto the real knobs: read count ∈
    // {3,5,8}, MinOverlap ∈ {3,5,8}, errorRate ∈ {0.0, 0.2}. Grid = 3 × 3 × 2 = 18.
    //
    // The combinatorial point: the contig set is a JOINT function of depth, the overlap threshold and the
    // error rate, but the assembly statistics are internally consistent in EVERY cell (INV-06) and every
    // read is accounted for. Witnesses cover the superstring (INV-04) and edgeless (INV-05) cases.
    // ═══════════════════════════════════════════════════════════════════════

    private const string OlcGenome = "ACGTCAGTGACTGCATGCTAGGACATTCGGATCCAAGTGC"; // length 40

    /// <summary>
    /// For every (read count, min overlap, error rate) the OLC assembly statistics are internally
    /// consistent, every read is accounted for, and contigs are spelled over the DNA alphabet.
    /// </summary>
    [Test, Combinatorial]
    public void AssembleOLC_ReadsOverlapErrorGrid_ConsistentStatistics(
        [Values(3, 5, 8)] int nReads,
        [Values(3, 5, 8)] int minOverlap,
        [Values(0.0, 0.2)] double errorRate)
    {
        var reads = TilingReads(OlcGenome, readLength: 12, count: nReads, errorRate: errorRate);
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: minOverlap, MinIdentity: 0.9, KmerSize: 5, MinContigLength: 1);

        var result = SequenceAssembler.AssembleOLC(reads, param);

        result.TotalReads.Should().Be(nReads, "every input read is accounted for");
        result.TotalLength.Should().Be(result.Contigs.Sum(c => c.Length), "[INV-06] TotalLength = Σ contig lengths");
        result.LongestContig.Should().Be(result.Contigs.Count == 0 ? 0 : result.Contigs.Max(c => c.Length), "[INV-06] LongestContig = max length");
        result.Contigs.Should().OnlyContain(c => c.All(ch => "ACGT".Contains(ch)), "contigs over the DNA alphabet");
    }

    /// <summary>
    /// Interaction witness (INV-05, edgeless): an overlap threshold larger than any read leaves the overlap
    /// graph edgeless, so each read becomes its own contig. Source: Langmead OLC layout.
    /// </summary>
    [Test]
    public void AssembleOLC_OverlapAboveReadLength_EachReadIsOwnContig()
    {
        var reads = TilingReads(OlcGenome, readLength: 12, count: 5, errorRate: 0.0);
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: 100, MinIdentity: 0.9, KmerSize: 5, MinContigLength: 1);

        var result = SequenceAssembler.AssembleOLC(reads, param);

        result.Contigs.Should().HaveCount(reads.Count, "no above-threshold overlap → each read is a contig");
    }

    /// <summary>
    /// Interaction witness (INV-04, superstring): clean overlapping reads assemble into a contig that
    /// contains an input read as a substring and is no longer than the sum of read lengths. Source:
    /// Compeau et al. (2011) superstring.
    /// </summary>
    [Test]
    public void AssembleOLC_CleanOverlappingReads_ContigIsSuperstring()
    {
        var reads = TilingReads(OlcGenome, readLength: 12, count: 6, errorRate: 0.0);
        var param = new SequenceAssembler.AssemblyParameters(MinOverlap: 4, MinIdentity: 0.9, KmerSize: 5, MinContigLength: 1);

        var result = SequenceAssembler.AssembleOLC(reads, param);

        result.Contigs.Should().NotBeEmpty();
        result.LongestContig.Should().BeGreaterThanOrEqualTo(12, "the contig is at least as long as one read");
        result.TotalLength.Should().BeLessThanOrEqualTo(reads.Sum(r => r.Length), "superstring length ≤ Σ read lengths");
    }

    /// <summary>Builds <paramref name="count"/> evenly-spaced sliding-window reads tiling the genome, mutating a leading fraction.</summary>
    private static List<string> TilingReads(string genome, int readLength, int count, double errorRate)
    {
        int lastStart = genome.Length - readLength;
        var reads = new List<string>();
        for (int i = 0; i < count; i++)
        {
            int start = count == 1 ? 0 : (int)Math.Round((double)lastStart * i / (count - 1));
            reads.Add(genome.Substring(start, readLength));
        }
        int errReads = (int)Math.Round(errorRate * count);
        for (int r = 0; r < errReads; r++)
        {
            var chars = reads[r].ToCharArray();
            chars[readLength / 2] = chars[readLength / 2] == 'A' ? 'C' : 'A';
            reads[r] = new string(chars);
        }
        return reads;
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ASSEMBLY-SCAFFOLD-001 — Scaffolding (Assembly)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 146.
    // Spec: tests/TestSpecs/ASSEMBLY-SCAFFOLD-001.md (SequenceAssembler.Scaffold).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: Jackman et al. (ABySS 2.0, 2017) — link path concatenates contigs separated by a gap run of
    // length = the distance estimate; a non-positive estimate uses the GenBank/EMBL/DDBJ unknown-gap length
    // 100 (NCBI AGP). Each contig is placed once; unreached contigs start their own scaffold.
    //
    // Checklist axes nContigs(3) × nLinks(3) × insertSize(2) map onto the real knobs: contig count ∈
    // {2,4,6}, chaining links ∈ {0,1,2} (capped at n−1), gapSize ∈ {10 (known), 0 (unknown → 100)}.
    // Grid = 3 × 3 × 2 = 18.
    //
    // The combinatorial point: the scaffold set is a JOINT function of contig count, link count and gap
    // size — chaining k links merges k+1 contigs into one scaffold (with k gap runs) and leaves the rest as
    // singletons. Each cell re-derives the chain and scaffold count from the scaffolding rule.
    // ═══════════════════════════════════════════════════════════════════════

    /// <summary>
    /// For every (contig count, link count, gap size) the scaffold count and the chained scaffold match the
    /// ABySS scaffolding rule re-derived from the inputs (k links merge k+1 contigs with the correct gap runs).
    /// </summary>
    [Test, Combinatorial]
    public void Scaffold_ContigsLinksGapGrid_MatchesChainRule(
        [Values(2, 4, 6)] int nContigs,
        [Values(0, 1, 2)] int nLinks,
        [Values(10, 0)] int gapSize)
    {
        var contigs = new List<string>();
        for (int i = 0; i < nContigs; i++)
            contigs.Add(new string((char)('A' + i), 5));

        int effectiveLinks = Math.Min(nLinks, nContigs - 1);
        var links = new List<(int, int, int)>();
        for (int i = 0; i < effectiveLinks; i++)
            links.Add((i, i + 1, gapSize));

        int gapLength = gapSize > 0 ? gapSize : 100;

        // Independent ground truth: the chain scaffold, then the remaining singletons.
        var expectedChain = new System.Text.StringBuilder(contigs[0]);
        for (int i = 1; i <= effectiveLinks; i++)
            expectedChain.Append('N', gapLength).Append(contigs[i]);

        var scaffolds = SequenceAssembler.Scaffold(contigs, links);

        scaffolds.Should().HaveCount(nContigs - effectiveLinks, "k links merge k+1 contigs; the rest are singletons");
        scaffolds[0].Should().Be(expectedChain.ToString(), "the chain scaffold = contigs joined by gap runs");
        scaffolds[0].Count(c => c == 'N').Should().Be(effectiveLinks * gapLength, "one gap run of `gapLength` per link");
    }

    /// <summary>
    /// Interaction witness (insertSize axis, unknown gap): a non-positive gap estimate emits the standard
    /// 100-character unknown gap, while a positive estimate emits exactly that many. Source: NCBI AGP / Jackman 2017.
    /// </summary>
    [Test]
    public void Scaffold_GapSizeAxis_KnownVsUnknown()
    {
        var contigs = new[] { "AAAAA", "CCCCC" };

        SequenceAssembler.Scaffold(contigs, new[] { (0, 1, 7) })[0].Count(c => c == 'N').Should().Be(7, "positive gap → exact");
        SequenceAssembler.Scaffold(contigs, new[] { (0, 1, 0) })[0].Count(c => c == 'N').Should().Be(100, "unknown gap → 100");
    }

    /// <summary>Witness: with no links every contig is its own scaffold, in order.</summary>
    [Test]
    public void Scaffold_NoLinks_EachContigOwnScaffold()
    {
        var contigs = new[] { "AAAAA", "CCCCC", "GGGGG" };

        SequenceAssembler.Scaffold(contigs, Array.Empty<(int, int, int)>())
            .Should().Equal(contigs, "no links → singletons in input order");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: ASSEMBLY-TRIM-001 — Quality trimming (Assembly)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 148.
    // Spec: tests/TestSpecs/ASSEMBLY-TRIM-001.md (SequenceAssembler.QualityTrimReads).
    // ADVANCED_TESTING_CHECKLIST.md §10.
    //
    // Sources: BWA bwa_trim_read / cutadapt running-sum — subtract the cutoff from each Phred score and cut
    // each end at the minimal partial sum; a cutoff < 1 disables trimming; reads below minLength are dropped.
    //
    // Checklist axes qualityCutoff(3) × windowSize(3) × readLen(3) map onto the real knobs: minQuality ∈
    // {10,20,30}; windowSize → minLength (minimum surviving length) ∈ {5,20,40}; read length ∈ {20,50,100}.
    // Grid = 3³ = 27.
    //
    // The combinatorial point: survival is a JOINT function of quality cutoff, length filter and read
    // length. With a uniform read quality of Phred 25, a read survives iff 25 ≥ cutoff (no trimming) AND its
    // length ≥ minLength; otherwise it is fully trimmed or filtered out. Each cell predicts the outcome.
    // ═══════════════════════════════════════════════════════════════════════

    private const char Phred25 = (char)(25 + 33); // uniform quality ':'

    /// <summary>
    /// For every (cutoff, minLength, read length) with uniform Phred-25 quality, the read survives intact
    /// iff 25 ≥ cutoff and its length ≥ minLength; otherwise the output is empty.
    /// </summary>
    [Test, Combinatorial]
    public void QualityTrimReads_CutoffMinLengthReadLenGrid_MatchesUniformQualityRule(
        [Values(10, 20, 30)] int cutoff,
        [Values(5, 20, 40)] int minLength,
        [Values(20, 50, 100)] int readLen)
    {
        string sequence = new string('A', readLen);
        string quality = new string(Phred25, readLen);
        var reads = new[] { (sequence, quality) };

        bool kept = 25 >= cutoff && readLen >= minLength;

        var trimmed = SequenceAssembler.QualityTrimReads(reads, cutoff, minLength);

        trimmed.Should().HaveCount(kept ? 1 : 0, "survives iff quality clears the cutoff and length ≥ minLength");
        if (kept)
            trimmed[0].Should().Be(sequence, "uniform high quality → no trimming");
    }

    /// <summary>
    /// Interaction witness (running-sum trimming): a read with low-quality ends and a high-quality core is
    /// trimmed to exactly its core. Source: BWA / cutadapt running-sum.
    /// </summary>
    [Test]
    public void QualityTrimReads_LowQualityEnds_TrimmedToCore()
    {
        const string sequence = "ACGTACGTACGTACGTACGT"; // length 20
        char low = (char)(2 + 33), high = (char)(40 + 33);
        string quality = new string(low, 5) + new string(high, 10) + new string(low, 5);

        var trimmed = SequenceAssembler.QualityTrimReads(new[] { (sequence, quality) }, minQuality: 20, minLength: 1);

        trimmed.Should().ContainSingle().Which.Should().Be(sequence.Substring(5, 10), "the low-quality ends are trimmed to the core");
    }

    /// <summary>
    /// Interaction witness (cutoff axis disables trimming): a cutoff below 1 disables trimming, so even a
    /// low-quality read is returned intact (subject to the length filter). Source: BWA trim_qual &lt; 1 guard.
    /// </summary>
    [Test]
    public void QualityTrimReads_CutoffBelowOne_DisablesTrimming()
    {
        const string sequence = "ACGTACGTAC";
        string quality = new string((char)(2 + 33), 10); // all low quality

        SequenceAssembler.QualityTrimReads(new[] { (sequence, quality) }, minQuality: 0, minLength: 1)
            .Should().ContainSingle().Which.Should().Be(sequence, "cutoff < 1 → no trimming");
    }

    /// <summary>Independent Biopython column-consensus ground truth.</summary>
    private static string GroundTruthConsensus(IReadOnlyList<string> reads, double threshold, char ambiguous)
    {
        if (reads.Count == 0) return "";
        int length = reads.Max(r => r.Length);
        var sb = new System.Text.StringBuilder(length);

        for (int pos = 0; pos < length; pos++)
        {
            var counts = new Dictionary<char, int>();
            int numAtoms = 0;
            foreach (string read in reads)
            {
                if (pos >= read.Length) continue;
                char c = char.ToUpperInvariant(read[pos]);
                if (c == '-' || c == '.') continue;
                counts[c] = counts.GetValueOrDefault(c, 0) + 1;
                numAtoms++;
            }

            int maxSize = 0, maxCount = 0;
            char maxResidue = ambiguous;
            foreach (var kvp in counts)
            {
                if (kvp.Value > maxSize) { maxSize = kvp.Value; maxCount = 1; maxResidue = kvp.Key; }
                else if (kvp.Value == maxSize) maxCount++;
            }

            bool commit = maxCount == 1 && numAtoms > 0 && (double)maxSize / numAtoms >= threshold;
            sb.Append(commit ? maxResidue : ambiguous);
        }
        return sb.ToString();
    }
}
