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
