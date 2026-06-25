namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Complexity area (SequenceComplexity,
/// Seqeron.Genomics.Analysis).
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of combinatorial testing.
/// Each grid cell carries a real business assertion; small grids use the exhaustive
/// <c>[Combinatorial]</c> product (a strict superset of pairwise).
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Complexity")]
public class ComplexityCombinatorialTests
{
    public enum SeqKind { Homopolymer, Repetitive, Diverse }

    private static string MakeSeq(SeqKind kind, int length) => kind switch
    {
        SeqKind.Homopolymer => new string('A', length),
        SeqKind.Repetitive => string.Concat(Enumerable.Repeat("AT", (length + 1) / 2))[..length],
        _ => string.Concat(Enumerable.Repeat("ACGT", (length + 3) / 4))[..length],
    };

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-COMPLEX-DUST-001 — DUST low-complexity score (Complexity)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 229.
    // Spec: tests/TestSpecs/SEQ-COMPLEX-DUST-001.md (canonical CalculateDustScore). ADVANCED §10.
    // Dimensions: seqType(3) × window(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Morgulis 2006 DUST; Li 2025): the DUST score is Σ_t c_t(c_t−1)/2 / (L−w+1) over the
    // counts c_t of the overlapping w-mers — high for repetitive/low-complexity sequence, low for
    // diverse sequence.
    //
    // Axis mapping (documented): seqType → composition (homopolymer/repetitive/diverse); window → the
    // DUST word size. The combinatorial point: the score equals an independent Σ c(c−1)/2/(L−w+1)
    // recomputation at every (seqType, word size) cell, and is non-negative.
    // ═══════════════════════════════════════════════════════════════════════

    private static double DustGroundTruth(string seq, int wordSize)
    {
        if (seq.Length < wordSize) return 0;
        var counts = new Dictionary<string, int>();
        int wordCount = seq.Length - wordSize + 1;
        for (int i = 0; i < wordCount; i++)
        {
            string w = seq.Substring(i, wordSize);
            counts[w] = counts.GetValueOrDefault(w) + 1;
        }
        double sum = counts.Values.Sum(c => (double)c * (c - 1) / 2.0);
        return sum / wordCount;
    }

    [Test, Combinatorial]
    public void ComplexDust_MatchesFormula_AcrossSeqTypeAndWordSize(
        [Values(SeqKind.Homopolymer, SeqKind.Repetitive, SeqKind.Diverse)] SeqKind seqType,
        [Values(2, 3, 4)] int wordSize)
    {
        string seq = MakeSeq(seqType, 60);

        double dust = SequenceComplexity.CalculateDustScore(seq, wordSize);
        dust.Should().BeApproximately(DustGroundTruth(seq, wordSize), 1e-9, "DUST = Σ c(c−1)/2 / (L−w+1)");
        dust.Should().BeGreaterThanOrEqualTo(0, "the DUST score is non-negative");
    }

    /// <summary>
    /// Interaction witness — DUST ranks complexity: a homopolymer (all words identical) scores far
    /// higher than a diverse sequence.
    /// </summary>
    [Test]
    public void ComplexDust_HomopolymerHigherThanDiverse()
    {
        SequenceComplexity.CalculateDustScore(MakeSeq(SeqKind.Homopolymer, 60), 3)
            .Should().BeGreaterThan(SequenceComplexity.CalculateDustScore(MakeSeq(SeqKind.Diverse, 60), 3),
                "a repetitive sequence is more low-complexity (higher DUST)");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-COMPLEX-KMER-001 — k-mer entropy complexity (Complexity)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 230.
    // Spec: tests/TestSpecs/SEQ-COMPLEX-KMER-001.md (canonical CalculateKmerEntropy). ADVANCED §10.
    // Dimensions: seqType(3) × k(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Shannon 1948): k-mer entropy = −Σ p·log2 p over the observed k-mer frequencies, in
    // [0, log2(#distinct)]; 0 for a homopolymer (one k-mer), maximal for a diverse sequence.
    //
    // Axis mapping (documented): seqType → composition; k → the k-mer length. The combinatorial point:
    // the entropy equals an independent −Σ p·log2 p recomputation at every cell, and a homopolymer
    // has entropy 0.
    // ═══════════════════════════════════════════════════════════════════════

    private static double KmerEntropyGroundTruth(string seq, int k)
    {
        if (seq.Length < k) return 0;
        var counts = new Dictionary<string, int>();
        int total = 0;
        for (int i = 0; i + k <= seq.Length; i++) { counts[seq.Substring(i, k)] = counts.GetValueOrDefault(seq.Substring(i, k)) + 1; total++; }
        double e = 0;
        foreach (int c in counts.Values) { double p = (double)c / total; e -= p * Math.Log2(p); }
        return e;
    }

    [Test, Combinatorial]
    public void ComplexKmer_MatchesFormula_AcrossSeqTypeAndK(
        [Values(SeqKind.Homopolymer, SeqKind.Repetitive, SeqKind.Diverse)] SeqKind seqType,
        [Values(1, 2, 3)] int k)
    {
        string seq = MakeSeq(seqType, 60);

        double entropy = SequenceComplexity.CalculateKmerEntropy(seq, k);
        entropy.Should().BeApproximately(KmerEntropyGroundTruth(seq, k), 1e-9, "k-mer entropy = −Σ p·log2 p");
        entropy.Should().BeGreaterThanOrEqualTo(-1e-12);

        if (seqType == SeqKind.Homopolymer)
            entropy.Should().BeApproximately(0.0, 1e-12, "a homopolymer has a single k-mer ⇒ entropy 0");
    }

    /// <summary>
    /// Interaction witness — k-mer entropy ranks complexity: diverse &gt; repetitive &gt; homopolymer.
    /// </summary>
    [Test]
    public void ComplexKmer_EntropyRanksComplexity()
    {
        double homo = SequenceComplexity.CalculateKmerEntropy(MakeSeq(SeqKind.Homopolymer, 60), 2);
        double rep = SequenceComplexity.CalculateKmerEntropy(MakeSeq(SeqKind.Repetitive, 60), 2);
        double div = SequenceComplexity.CalculateKmerEntropy(MakeSeq(SeqKind.Diverse, 60), 2);

        homo.Should().Be(0.0);
        rep.Should().BeGreaterThan(homo);
        div.Should().BeGreaterThan(rep, "a more diverse sequence has higher k-mer entropy");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: SEQ-COMPLEX-WINDOW-001 — Windowed complexity profile (Complexity)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 231.
    // Spec: tests/TestSpecs/SEQ-COMPLEX-WINDOW-001.md (canonical CalculateWindowedComplexity). ADVANCED §10.
    // Dimensions: seqType(3) × window(3) × step(2). Grid 3×3×2 = 18 (full, exhaustive).
    //
    // Model (Troyanskaya 2002 complexity profile): a sliding window emits, per fully-contained window,
    // its Shannon entropy (bits) and linguistic complexity with the window's coordinates.
    //
    // The combinatorial point: across composition, window size and step, the number of points equals
    // the tiling count ⌊(L−w)/step⌋+1, each window's coordinates are consistent, and the reported
    // Shannon entropy equals the independent per-window CalculateShannonEntropy and lies in [0, log2 4].
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void ComplexWindow_TilingAndPerWindowEntropy_AcrossSeqTypeWindowStep(
        [Values(SeqKind.Homopolymer, SeqKind.Repetitive, SeqKind.Diverse)] SeqKind seqType,
        [Values(10, 20, 30)] int windowSize,
        [Values(5, 10)] int step)
    {
        string seq = MakeSeq(seqType, 60);
        var points = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(seq), windowSize, step).ToList();

        int expected = (seq.Length - windowSize) / step + 1;
        points.Should().HaveCount(expected, "windows = ⌊(L−w)/step⌋ + 1");

        foreach (var p in points)
        {
            p.WindowEnd.Should().Be(p.WindowStart + windowSize - 1, "window spans windowSize");
            p.Position.Should().Be(p.WindowStart + windowSize / 2, "position is the window midpoint");
            p.ShannonEntropy.Should().BeInRange(0.0, 2.0 + 1e-9, "DNA entropy ≤ log2 4 = 2 bits");
            p.ShannonEntropy.Should().BeApproximately(
                SequenceComplexity.CalculateShannonEntropy(seq.Substring(p.WindowStart, windowSize)), 1e-9,
                "the reported entropy is the per-window Shannon entropy");
        }
    }

    /// <summary>
    /// Interaction witness — a diverse window carries more Shannon entropy than a homopolymer window
    /// at the same coordinates.
    /// </summary>
    [Test]
    public void ComplexWindow_DiverseHasHigherEntropyThanHomopolymer()
    {
        double homo = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(MakeSeq(SeqKind.Homopolymer, 40)), 20, 10).First().ShannonEntropy;
        double div = SequenceComplexity.CalculateWindowedComplexity(new DnaSequence(MakeSeq(SeqKind.Diverse, 40)), 20, 10).First().ShannonEntropy;
        homo.Should().Be(0.0, "a homopolymer window has zero entropy");
        div.Should().BeGreaterThan(homo, "a diverse window has positive entropy");
    }
}
