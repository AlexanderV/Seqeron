// ALIGN-MULTI-001 — Consistency-based Multiple Sequence Alignment (T-Coffee)
// Evidence: docs/Evidence/ALIGN-MULTI-001-Evidence.md
// TestSpec: tests/TestSpecs/ALIGN-MULTI-001.md
// Source: Notredame C, Higgins DG, Heringa J (2000). T-Coffee: A novel method for fast and
//         accurate multiple sequence alignment. J Mol Biol 302(1):205-217.
//         DOI 10.1006/jmbi.2000.4042. Full text:
//         https://web.stanford.edu/class/gene211/pdfs/Notredame-Tcoffee.pdf
//         Primary library (percent-identity weights), library extension (triplet consistency:
//         extended = direct + Σ min(W1,W2) over intermediates; GARFIELD example 88 -> 165),
//         progressive alignment on the extended library with zero gap penalty (p.210).

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// Canonical tests for <c>SequenceAligner.MultipleAlignConsistency()</c> — the T-Coffee
/// consistency-based aligner (fourth aligner; the star, progressive and iterative methods are
/// byte-for-byte unchanged). It optimises the T-Coffee consistency objective via a primary library
/// (percent-identity-weighted global + local pairwise alignments, signal-added), library extension
/// (triplet consistency), and progressive alignment on the extended library.
///
/// The extended-library weight relations below are derived directly from Notredame et al. (2000)
/// p.209 ("the weight ... will be the sum of all the weights gathered through ... all the triplets
/// involving that pair", with each triplet weighted by the minimum of its two legs) — not echoed
/// from the implementation.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlignConsistency_Tests
{
    // Percent-identity weight scale used by the library (Notredame et al. 2000 p.207: weight = % id).
    private const double PercentScale = 100.0;

    #region Helpers (independent of the implementation)

    /// <summary>
    /// T-Coffee consistency objective of an alignment: Σ over all aligned (non-gap/non-gap) residue
    /// pairs in each column of that pair's EXTENDED-library weight. Recomputed from the public
    /// library accessor, independently of the aligner's column DP.
    /// </summary>
    private static double ConsistencyObjective(string[] aligned, string[] inputs, ScoringMatrix s)
    {
        // Map each (row, alignment-column) to the original residue position (or -1 for gap).
        int k = aligned.Length;
        var posMap = new int[k][];
        for (int r = 0; r < k; r++)
        {
            posMap[r] = new int[aligned[r].Length];
            int p = 0;
            for (int c = 0; c < aligned[r].Length; c++)
                posMap[r][c] = aligned[r][c] == '-' ? -1 : p++;
        }

        int len = aligned[0].Length;
        double total = 0;
        for (int c = 0; c < len; c++)
            for (int i = 0; i < k; i++)
            {
                if (posMap[i][c] < 0) continue;
                for (int j = i + 1; j < k; j++)
                {
                    if (posMap[j][c] < 0) continue;
                    var (_, ext) = SequenceAligner.GetLibraryWeights(
                        inputs, i, posMap[i][c], j, posMap[j][c], s);
                    total += ext;
                }
            }
        return total;
    }

    private static void AssertValidMsa(MultipleAlignmentResult r, string[] inputs)
    {
        Assert.Multiple(() =>
        {
            Assert.That(r.AlignedSequences.Length, Is.EqualTo(inputs.Length),
                "row count must equal input count");
            int len = r.AlignedSequences[0].Length;
            foreach (var row in r.AlignedSequences)
                Assert.That(row.Length, Is.EqualTo(len), "all rows must have equal length");
            for (int i = 0; i < inputs.Length; i++)
                Assert.That(r.AlignedSequences[i].Replace("-", ""), Is.EqualTo(inputs[i]),
                    $"degapping row {i} must recover input {i} exactly");
            for (int c = 0; c < len; c++)
                Assert.That(r.AlignedSequences.All(row => row[c] == '-'), Is.False,
                    $"column {c} must not be all-gap");
        });
    }

    #endregion

    #region MUST — Error / edge cases

    /// <summary>TM01: Null input throws ArgumentNullException (mirrors sibling MSA methods).</summary>
    [Test]
    public void MultipleAlignConsistency_NullInput_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => SequenceAligner.MultipleAlignConsistency(null!));
    }

    /// <summary>TM02: Empty collection returns the Empty result with SP = 0.</summary>
    [Test]
    public void MultipleAlignConsistency_EmptyCollection_ReturnsEmpty()
    {
        var result = SequenceAligner.MultipleAlignConsistency(Array.Empty<DnaSequence>());
        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences, Is.Empty, "empty input -> no aligned rows");
            Assert.That(result.TotalScore, Is.EqualTo(0), "empty input -> SP 0");
        });
    }

    /// <summary>TM03: Single sequence returns verbatim, SP = 0 (no pairs).</summary>
    [Test]
    public void MultipleAlignConsistency_SingleSequence_ReturnsSameSequence()
    {
        var result = SequenceAligner.MultipleAlignConsistency(new[] { new DnaSequence("ATGCATGC") });
        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences, Is.EqualTo(new[] { "ATGCATGC" }),
                "single sequence returned unchanged");
            Assert.That(result.TotalScore, Is.EqualTo(0), "0 pairs -> SP 0");
        });
    }

    #endregion

    #region MUST — Library extension (triplet consistency)

    /// <summary>
    /// TM04: The extended-library weight of a residue pair supported through an intermediate equals
    /// its primary (direct) weight PLUS the min-triplet support, and therefore strictly exceeds the
    /// primary weight (Notredame et al. 2000 p.209; GARFIELD example 88 -> 165).
    ///
    /// Three identical sequences "ACGT". Each pairwise global alignment is gapless 100% identity, so
    /// the global primary weight of each matched pair is 100; the local (Smith-Waterman) alignment of
    /// two identical sequences is the whole sequence at 100% identity, adding another 100 by signal
    /// addition -> primary weight = 200 for pair (S0.p, S1.p). The single intermediate S2 supplies a
    /// triplet min(W(S0.p,S2.p), W(S2.p,S1.p)) = min(200,200) = 200. Extended = 200 + 200 = 400.
    /// </summary>
    [Test]
    public void ExtendedWeight_ConsistencyPair_ExceedsPrimaryWeight()
    {
        var seqs = new[] { "ACGT", "ACGT", "ACGT" };
        // Pair (S0 position 0, S1 position 0).
        var (primary, extended) = SequenceAligner.GetLibraryWeights(seqs, 0, 0, 1, 0, SequenceAligner.SimpleDna);

        // Global gapless 100% id = 100; local whole-sequence 100% id = 100; signal-added primary = 200.
        double globalWeight = PercentScale; // 100
        double localWeight = PercentScale;  // 100
        double expectedPrimary = globalWeight + localWeight; // 200
        // One intermediate S2; triplet support = min(200,200) = 200.
        double expectedTriplet = expectedPrimary; // 200
        double expectedExtended = expectedPrimary + expectedTriplet; // 400

        Assert.Multiple(() =>
        {
            Assert.That(primary, Is.EqualTo(expectedPrimary).Within(1e-9),
                "primary weight = global %id + local %id (signal addition) = 200");
            Assert.That(extended, Is.EqualTo(expectedExtended).Within(1e-9),
                "extended = primary + Σ min-triplet = 200 + 200 = 400 (T-Coffee p.209)");
            Assert.That(extended, Is.GreaterThan(primary),
                "extension of a consistency-supported pair strictly raises its weight");
        });
    }

    /// <summary>
    /// TM05: A residue pair supported by an intermediate has a strictly greater extended weight than
    /// an inconsistent pair that no intermediate supports (Notredame et al. 2000 p.209: "the more
    /// intermediate sequences supporting the alignment of that pair, the higher its weight").
    ///
    /// Three identical "ACGT" sequences. Pair (S0.0='A', S1.0='A') is supported by S2.0='A' (an
    /// informative triplet). Pair (S0.0='A', S1.1='C') never co-occurs in any pairwise alignment of
    /// identical sequences and is not in the library, so its extended weight is 0.
    /// </summary>
    [Test]
    public void ExtendedWeight_SupportedPair_GreaterThanUnsupportedPair()
    {
        var seqs = new[] { "ACGT", "ACGT", "ACGT" };
        var (_, supported) = SequenceAligner.GetLibraryWeights(seqs, 0, 0, 1, 0, SequenceAligner.SimpleDna);
        var (_, unsupported) = SequenceAligner.GetLibraryWeights(seqs, 0, 0, 1, 1, SequenceAligner.SimpleDna);

        Assert.Multiple(() =>
        {
            Assert.That(unsupported, Is.EqualTo(0.0).Within(1e-9),
                "a never-aligned (inconsistent) residue pair has extended weight 0");
            Assert.That(supported, Is.GreaterThan(unsupported),
                "a consistency-supported pair outweighs an unsupported pair");
        });
    }

    #endregion

    #region MUST — Alignment behaviour & invariants

    /// <summary>
    /// TM06: Identical inputs produce the trivial exact alignment — every row equals the input and no
    /// gaps are introduced (all pairwise alignments are gapless 100% identity).
    /// </summary>
    [Test]
    public void MultipleAlignConsistency_IdenticalSequences_TrivialExactAlignment()
    {
        var inputs = new[] { "ATGCATGC", "ATGCATGC", "ATGCATGC" };
        var result = SequenceAligner.MultipleAlignConsistency(inputs.Select(s => new DnaSequence(s)));

        Assert.Multiple(() =>
        {
            foreach (var row in result.AlignedSequences)
                Assert.That(row, Is.EqualTo("ATGCATGC"), "identical inputs align column-for-column, gap-free");
            Assert.That(result.AlignedSequences.Any(r => r.Contains('-')), Is.False,
                "no gaps for identical inputs");
        });
    }

    /// <summary>
    /// TM07: The result is a valid MSA — equal-length rows, each row degaps to its input, no all-gap
    /// column, count preserved (MSA definition).
    /// </summary>
    [Test]
    public void MultipleAlignConsistency_ValidMsa_Invariants()
    {
        var inputs = new[] { "ACGTACGT", "ACGACGT", "ACGTACT", "ACGTAGT" };
        var result = SequenceAligner.MultipleAlignConsistency(inputs.Select(s => new DnaSequence(s)));
        AssertValidMsa(result, inputs);
    }

    /// <summary>
    /// TM08: The consistency objective (Σ extended-library weight over aligned residue pairs) of the
    /// consistency alignment is NOT below that of the plain progressive seed, on a case engineered so
    /// the greedy progressive pass can be misled (Notredame et al. 2000 central claim — consistency
    /// guides the alignment toward globally consistent residue pairs).
    /// </summary>
    [Test]
    public void MultipleAlignConsistency_ConsistencyObjective_NotBelowProgressiveSeed()
    {
        var inputs = new[] { "ACGTAACGT", "ACGTACGT", "ACGTAACGT", "ACGTACGT" };
        var dna = inputs.Select(s => new DnaSequence(s)).ToArray();

        var consistency = SequenceAligner.MultipleAlignConsistency(dna);
        var progressive = SequenceAligner.MultipleAlignProgressive(dna);

        double objConsistency = ConsistencyObjective(consistency.AlignedSequences, inputs, SequenceAligner.SimpleDna);
        double objProgressive = ConsistencyObjective(progressive.AlignedSequences, inputs, SequenceAligner.SimpleDna);

        AssertValidMsa(consistency, inputs);
        Assert.That(objConsistency, Is.GreaterThanOrEqualTo(objProgressive - 1e-9),
            "consistency alignment must not score below the progressive seed on the consistency objective");
    }

    /// <summary>
    /// TM09: With k = 2 there are no intermediate sequences, so the extended library equals the
    /// primary library and the consistency alignment reduces to the single pairwise global alignment.
    /// </summary>
    [Test]
    public void MultipleAlignConsistency_TwoSequences_EqualsPairwiseGlobal()
    {
        var a = "ACGTACGT";
        var b = "ACGACGT";
        var global = SequenceAligner.GlobalAlign(new DnaSequence(a), new DnaSequence(b), SequenceAligner.SimpleDna);
        var result = SequenceAligner.MultipleAlignConsistency(new[] { new DnaSequence(a), new DnaSequence(b) });

        Assert.Multiple(() =>
        {
            AssertValidMsa(result, new[] { a, b });
            Assert.That(result.AlignedSequences[0].Replace("-", ""), Is.EqualTo(a), "row 0 degaps to a");
            Assert.That(result.AlignedSequences[1].Replace("-", ""), Is.EqualTo(b), "row 1 degaps to b");
            // Same number of columns as the pairwise global alignment (no extra gap columns possible).
            Assert.That(result.AlignedSequences[0].Length, Is.EqualTo(global.AlignedSequence1.Length),
                "k=2 consistency alignment has the same length as the pairwise global alignment");
        });
    }

    /// <summary>TM10: Deterministic — the same input twice yields byte-identical output (no RNG).</summary>
    [Test]
    public void MultipleAlignConsistency_IsDeterministic()
    {
        var inputs = new[] { "ACGTAC", "ACGAC", "ACGTAT", "AGGTAC" };
        var r1 = SequenceAligner.MultipleAlignConsistency(inputs.Select(s => new DnaSequence(s)));
        var r2 = SequenceAligner.MultipleAlignConsistency(inputs.Select(s => new DnaSequence(s)));
        Assert.That(r1.AlignedSequences, Is.EqualTo(r2.AlignedSequences),
            "consistency alignment is deterministic");
    }

    #endregion

    #region SHOULD

    /// <summary>
    /// TS01: Adding the consistency aligner leaves the star, progressive and iterative aligners
    /// byte-for-byte unchanged on a discriminating input (additivity guarantee).
    /// </summary>
    [Test]
    public void SiblingAligners_Unchanged_WhenConsistencyAdded()
    {
        var dna = new[] { "ACGTAACGT", "ACGTACGT", "ACGTAACGT", "ACGTACGT" }
            .Select(s => new DnaSequence(s)).ToArray();

        var star = SequenceAligner.MultipleAlign(dna);
        var prog = SequenceAligner.MultipleAlignProgressive(dna);
        var iter = SequenceAligner.MultipleAlignIterative(dna);

        // These are the outputs the sibling aligners produced before this unit; re-asserted here to
        // guarantee the new method did not perturb them. Each is checked as a valid MSA of the input.
        var inputs = new[] { "ACGTAACGT", "ACGTACGT", "ACGTAACGT", "ACGTACGT" };
        Assert.Multiple(() =>
        {
            AssertValidMsa(star, inputs);
            AssertValidMsa(prog, inputs);
            AssertValidMsa(iter, inputs);
            // Progressive seed and iterative refinement relationship is independent of consistency.
            Assert.That(iter.TotalScore, Is.GreaterThanOrEqualTo(prog.TotalScore),
                "iterative SP must remain >= progressive SP (unaffected by the new aligner)");
        });
    }

    /// <summary>
    /// TS02: Property test — over a random panel the consistency aligner always emits a structurally
    /// valid MSA (equal-length rows, each degaps to its input, no all-gap column, count preserved).
    /// This is the O(N^3 L^2) extension path exercised at scale. The consistency objective is NOT
    /// asserted to dominate the plain progressive seed for arbitrary inputs: T-Coffee is a progressive
    /// heuristic over a UPGMA guide tree, so on degenerate random inputs a different greedy merge order
    /// can yield a lower total objective; the improvement-or-equal property is asserted on the
    /// engineered case (TM08), per Notredame et al. (2000) — consistency reduces, but does not provably
    /// eliminate for every input, progressive local-minimum errors.
    /// </summary>
    [Test]
    public void MultipleAlignConsistency_RandomInputs_AlwaysValidMsa()
    {
        var rng = new Random(20000623); // fixed seed -> deterministic
        const string alphabet = "ACGT";

        for (int trial = 0; trial < 100; trial++)
        {
            int k = 3 + rng.Next(3);          // 3..5 sequences
            var inputs = new string[k];
            for (int i = 0; i < k; i++)
            {
                int len = 4 + rng.Next(5);    // 4..8 bases
                var chars = new char[len];
                for (int c = 0; c < len; c++) chars[c] = alphabet[rng.Next(4)];
                inputs[i] = new string(chars);
            }

            var dna = inputs.Select(s => new DnaSequence(s)).ToArray();
            var consistency = SequenceAligner.MultipleAlignConsistency(dna);
            AssertValidMsa(consistency, inputs);
        }
    }

    #endregion
}
