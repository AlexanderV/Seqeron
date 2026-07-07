// ALIGN-MULTI-001 — Iterative Refinement of Progressive Multiple Sequence Alignment
// Evidence: docs/Evidence/ALIGN-MULTI-001-Evidence.md
// TestSpec: tests/TestSpecs/ALIGN-MULTI-001.md
// Source: Edgar RC (2004). MUSCLE: multiple sequence alignment with high accuracy and
//         high throughput. Nucleic Acids Res 32(5):1792-1797. Stage 3 (tree-dependent
//         restricted partitioning), steps 3.1-3.4.
//         Barton GJ, Sternberg MJ (1987). J Mol Biol 198(2):327-337 (iterative refinement).
//         Wikipedia "Multiple sequence alignment" (sum-of-pairs score definition).

namespace Seqeron.Genomics.Tests.Unit.Alignment;

/// <summary>
/// Canonical tests for <c>SequenceAligner.MultipleAlignIterative()</c> — iterative refinement of a
/// progressive MSA via MUSCLE-style tree-dependent restricted partitioning (Edgar 2004, Stage 3).
/// Removes the single-pass "once a gap, always a gap" limitation of MultipleAlignProgressive.
///
/// Algorithm under test (third aligner; star MultipleAlign and progressive MultipleAlignProgressive
/// are byte-for-byte unchanged):
///   1. Build the progressive (Feng-Doolittle / UPGMA) seed + keep its guide tree.
///   2. For each internal guide-tree edge (leaves-first, deterministic), split the alignment into
///      two leaf groups, realign the two sub-profiles with the existing profile-profile NW.
///   3. Accept the re-alignment only if the sum-of-pairs (SP) score does not decrease.
///   4. Iterate full passes until no edge improves (convergence) or maxIterations is reached.
///
/// Sources (retrieved 2026-06-23):
///   - Edgar RC (2004) NAR 32(5):1792-1797 — Stage 3 steps 3.1-3.4:
///     "An edge is chosen from TREE2 ... divided into two subtrees by deleting the edge. The
///     profile of the multiple alignment in each subtree is computed. A new multiple alignment is
///     produced by re-aligning the two profiles. If the SP score is improved, the new alignment is
///     kept ... repeated until convergence or until a user-defined limit is reached."
///     https://academic.oup.com/nar/article/32/5/1792/2380623
///   - Barton &amp; Sternberg (1987) J Mol Biol 198:327-337 — iterative refinement of an existing
///     alignment. https://pubmed.ncbi.nlm.nih.gov/3430611/
///   - Wikipedia: Multiple sequence alignment — SP score = "sum of all of the pairs of characters
///     at each position in the alignment". https://en.wikipedia.org/wiki/Multiple_sequence_alignment
///
/// Expected SP scores and aligned columns below are hand-derived from the SP definition (not echoed
/// from the implementation); see per-test comments for the derivations.
/// </summary>
[TestFixture]
[Category("Alignment")]
[Category("ALIGN-MULTI-001")]
public class SequenceAligner_MultipleAlignIterative_Tests
{
    #region MUST — Error / edge cases (null, empty, 1, 2 sequences)

    /// <summary>M01: Null input throws ArgumentNullException (.NET convention, mirrors sibling MSA).</summary>
    [Test]
    public void MultipleAlignIterative_NullInput_ThrowsArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() =>
            SequenceAligner.MultipleAlignIterative(null!));
    }

    /// <summary>M02: maxIterations &lt; 1 throws ArgumentOutOfRangeException (positive cap required).</summary>
    [Test]
    public void MultipleAlignIterative_NonPositiveMaxIterations_Throws()
    {
        var input = new[] { new DnaSequence("ACGT"), new DnaSequence("ACGT") };
        Assert.Throws<ArgumentOutOfRangeException>(() =>
            SequenceAligner.MultipleAlignIterative(input, null, 0));
    }

    /// <summary>M03: Empty collection returns the Empty result.</summary>
    [Test]
    public void MultipleAlignIterative_EmptyCollection_ReturnsEmpty()
    {
        var result = SequenceAligner.MultipleAlignIterative(Array.Empty<DnaSequence>());

        Assert.Multiple(() =>
        {
            Assert.That(result, Is.EqualTo(MultipleAlignmentResult.Empty));
            Assert.That(result.AlignedSequences, Is.Empty);
            Assert.That(result.Consensus, Is.Empty);
            Assert.That(result.TotalScore, Is.EqualTo(0));
        });
    }

    /// <summary>M04: Single sequence is returned verbatim (no pairs, no gaps, SP = 0).</summary>
    [Test]
    public void MultipleAlignIterative_SingleSequence_ReturnsSameSequence()
    {
        var result = SequenceAligner.MultipleAlignIterative(new[] { new DnaSequence("ATGCATGC") });

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(1));
            Assert.That(result.AlignedSequences[0], Is.EqualTo("ATGCATGC"));
            Assert.That(result.Consensus, Is.EqualTo("ATGCATGC"));
            Assert.That(result.TotalScore, Is.EqualTo(0));
        });
    }

    /// <summary>
    /// M05: Two sequences — no internal guide-tree edge exists, so refinement is a no-op and the
    /// result equals the progressive seed exactly.
    /// Hand-derived NW("ACGT","ACT") with SimpleDna: "ACGT" / "AC-T"; SP = +3 (A,C,T) + GapExtend(−1)
    /// for the one gap column = 2.
    /// </summary>
    [Test]
    public void MultipleAlignIterative_TwoSequences_EqualsProgressiveSeed()
    {
        var input = new[] { new DnaSequence("ACGT"), new DnaSequence("ACT") };

        var iter = SequenceAligner.MultipleAlignIterative(input);
        var prog = SequenceAligner.MultipleAlignProgressive(input);

        Assert.Multiple(() =>
        {
            Assert.That(iter.AlignedSequences, Is.EqualTo(new[] { "ACGT", "AC-T" }));
            Assert.That(iter.TotalScore, Is.EqualTo(2));
            Assert.That(iter.AlignedSequences, Is.EqualTo(prog.AlignedSequences),
                "With no internal edge, refinement must equal the progressive seed");
        });
    }

    #endregion

    #region MUST — Iterative refinement corrects an early-gap-placement error (the headline fix)

    /// <summary>
    /// M06: HEADLINE CORRECTION CASE. Progressive MSA misplaces an internal gap; iterative
    /// refinement relocates it and raises the SP score by exactly +2. This is the "once a gap,
    /// always a gap" error being corrected.
    ///
    /// Input ["CGA","GAGAT","CGC","GAC"], SimpleDna (match +1, mismatch −1, gap −1).
    ///
    /// Progressive seed (hand-derived; the early profile merge fixes a gap before the final base):
    ///   -CG-A
    ///   GAGAT
    ///   -CG-C
    ///   --GAC
    ///   SP by column (6 pairs each): col0 {-,G,-,-}=−3; col1 {C,A,C,-}=−4; col2 {G,G,G,G}=+6;
    ///   col3 {-,A,-,A}=−3; col4 {A,T,C,C}=−4  ⇒  SP = −8.
    ///
    /// Iterative refinement re-splits at a guide-tree edge and realigns, relocating the internal
    /// gap of rows 0 and 2 to the terminal position:
    ///   -CGA-
    ///   GAGAT
    ///   -CGC-
    ///   --GAC
    ///   SP by column: col0 {-,G,-,-}=−3; col1 {C,A,C,-}=−4; col2 {G,G,G,G}=+6;
    ///   col3 {A,A,C,A}=0; col4 {-,T,-,C}=−5  ⇒  SP = −6.
    ///
    /// The gain (+2) comes from col3 −3 → 0 and col4 −4 → −5 (net +3 − 1 = +2). Both alignments
    /// degap to the same inputs, so this is a pure gap RELOCATION the single-pass method cannot make.
    /// </summary>
    [Test]
    public void MultipleAlignIterative_CorrectsEarlyGapPlacement_RaisesSpScore()
    {
        var input = new[] { "CGA", "GAGAT", "CGC", "GAC" }
            .Select(s => new DnaSequence(s)).ToArray();

        var prog = SequenceAligner.MultipleAlignProgressive(input);
        var iter = SequenceAligner.MultipleAlignIterative(input);

        Assert.Multiple(() =>
        {
            // Progressive seed is the suboptimal arrangement (hand-derived SP = −8).
            Assert.That(prog.AlignedSequences,
                Is.EqualTo(new[] { "-CG-A", "GAGAT", "-CG-C", "--GAC" }),
                "Progressive seed places the gap internally (the error to be corrected)");
            Assert.That(prog.TotalScore, Is.EqualTo(-8));

            // Refinement relocates the gap to the terminal column (hand-derived SP = −6).
            Assert.That(iter.AlignedSequences,
                Is.EqualTo(new[] { "-CGA-", "GAGAT", "-CGC-", "--GAC" }),
                "Refinement relocates the internal gap to the terminal position");
            Assert.That(iter.TotalScore, Is.EqualTo(-6),
                "Hand-derived corrected SP score");

            // The fix strictly improves SP by +2.
            Assert.That(iter.TotalScore - prog.TotalScore, Is.EqualTo(2),
                "Iterative refinement strictly improves SP by exactly 2 over the seed");

            // Integrity: both arrangements degap to the inputs (a relocation, not a content change).
            for (int i = 0; i < input.Length; i++)
                Assert.That(iter.AlignedSequences[i].Replace("-", ""),
                    Is.EqualTo(input[i].Sequence), $"Row {i} degaps to its input");
        });
    }

    #endregion

    #region MUST — Monotonicity: refined SP is never below the progressive seed

    /// <summary>
    /// M07: SP MONOTONICITY (Edgar 2004 step 3.4: "If the SP score is improved, the new alignment
    /// is kept, otherwise it is discarded"). Over a panel of inputs the iterative SP is always
    /// ≥ the progressive seed's SP — never worse. Each case is checked against the seed computed
    /// by the existing progressive aligner (no echoing of the iterative output).
    /// </summary>
    [Test]
    public void MultipleAlignIterative_SpScore_NeverBelowProgressiveSeed()
    {
        string[][] panel =
        {
            new[] { "CGA", "GAGAT", "CGC", "GAC" },          // headline win
            new[] { "ACGT", "ACGT", "AGT" },                 // already optimal seed
            new[] { "GGAC", "TGT", "CTGGT", "TGG" },         // win
            new[] { "ACTAG", "TCAT", "TAAGG", "CTG" },       // win
            new[] { "AAGA", "GAT", "TAAT", "TCTGA" },        // win
            new[] { "ATGCATGC", "ATGCATGC", "ATGCATGC" },    // identical -> no change
            new[] { "AAGAA", "AACAA", "GGTGG", "GGTGG" },    // two clusters
        };

        foreach (var caseSeqs in panel)
        {
            var seqs = caseSeqs.Select(s => new DnaSequence(s)).ToArray();
            int progSp = SequenceAligner.MultipleAlignProgressive(seqs).TotalScore;
            int iterSp = SequenceAligner.MultipleAlignIterative(seqs).TotalScore;

            Assert.That(iterSp, Is.GreaterThanOrEqualTo(progSp),
                $"Iterative SP ({iterSp}) must be ≥ progressive seed SP ({progSp}) for " +
                $"[{string.Join(",", caseSeqs)}]");
        }
    }

    #endregion

    #region MUST — Alignment integrity invariants (length, degap, no all-gap column, count)

    /// <summary>
    /// M08: Output rows are equal length; each row degaps to its input; row count is preserved;
    /// no column is entirely gaps (Wikipedia MSA invariants), on a refined alignment.
    /// </summary>
    [Test]
    public void MultipleAlignIterative_Invariants_Hold()
    {
        var originals = new[] { "CGA", "GAGAT", "CGC", "GAC", "ACGGA" };
        var result = SequenceAligner.MultipleAlignIterative(
            originals.Select(s => new DnaSequence(s)));

        int len = result.AlignedSequences[0].Length;
        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences.Length, Is.EqualTo(originals.Length),
                "Row count preserved");
            Assert.That(result.AlignedSequences.All(s => s.Length == len), Is.True,
                "All rows equal length");
            Assert.That(len, Is.GreaterThanOrEqualTo(originals.Max(s => s.Length)),
                "Column count ≥ max input length");
            for (int i = 0; i < originals.Length; i++)
                Assert.That(result.AlignedSequences[i].Replace("-", ""), Is.EqualTo(originals[i]),
                    $"Removing gaps recovers original sequence {i}");
            Assert.That(result.Consensus.Length, Is.EqualTo(len),
                "Consensus length = aligned length");
            for (int col = 0; col < len; col++)
            {
                bool allGaps = result.AlignedSequences.All(s => s[col] == '-');
                Assert.That(allGaps, Is.False, $"Column {col} must not be entirely gaps");
            }
        });
    }

    #endregion

    #region MUST — Idempotence: refining an already-optimal alignment changes nothing

    /// <summary>
    /// M09: IDEMPOTENCE. Refining identical sequences (the alignment is already optimal: gap-free,
    /// every column a perfect match) leaves it unchanged. SP = C(3,2)=3 pairs × 8 matched columns
    /// × Match(1) = 24; no gaps introduced.
    /// </summary>
    [Test]
    public void MultipleAlignIterative_AlreadyOptimal_ChangesNothing()
    {
        var input = new[]
        {
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGCATGC"),
            new DnaSequence("ATGCATGC")
        };

        var result = SequenceAligner.MultipleAlignIterative(input);

        Assert.Multiple(() =>
        {
            Assert.That(result.AlignedSequences,
                Is.EqualTo(new[] { "ATGCATGC", "ATGCATGC", "ATGCATGC" }));
            Assert.That(result.AlignedSequences.Any(s => s.Contains('-')), Is.False,
                "An already-optimal gap-free alignment must not gain gaps under refinement");
            Assert.That(result.TotalScore, Is.EqualTo(24));
        });
    }

    /// <summary>
    /// M10: IDEMPOTENCE (fixed point). Running the refiner on the already-refined result (its rows
    /// degapped back to inputs) reproduces the same SP score — the refined alignment is a fixed
    /// point of the procedure.
    /// </summary>
    [Test]
    public void MultipleAlignIterative_IsFixedPoint()
    {
        var input = new[] { "CGA", "GAGAT", "CGC", "GAC" }
            .Select(s => new DnaSequence(s)).ToArray();

        var first = SequenceAligner.MultipleAlignIterative(input);
        // Degap the refined rows back to raw inputs and refine again.
        var degapped = first.AlignedSequences
            .Select(s => new DnaSequence(s.Replace("-", ""))).ToArray();
        var second = SequenceAligner.MultipleAlignIterative(degapped);

        Assert.That(second.TotalScore, Is.EqualTo(first.TotalScore),
            "Refining the refined (degapped) result reproduces the same SP score (fixed point)");
    }

    #endregion

    #region MUST — Convergence within the cap

    /// <summary>
    /// M11: CONVERGENCE. The headline correction converges within ONE pass: capping maxIterations
    /// at 1 yields the same fully-refined SP score as the default cap, so no further passes change
    /// the result (Edgar 2004: "repeated until convergence or until a user-defined limit").
    /// </summary>
    [Test]
    public void MultipleAlignIterative_ConvergesWithinCap()
    {
        var input = new[] { "CGA", "GAGAT", "CGC", "GAC" }
            .Select(s => new DnaSequence(s)).ToArray();

        var oneIter = SequenceAligner.MultipleAlignIterative(input, null, maxIterations: 1);
        var manyIter = SequenceAligner.MultipleAlignIterative(input, null, maxIterations: 100);

        Assert.Multiple(() =>
        {
            Assert.That(oneIter.TotalScore, Is.EqualTo(-6),
                "One pass already reaches the converged SP score");
            Assert.That(manyIter.TotalScore, Is.EqualTo(oneIter.TotalScore),
                "Additional passes do not change the converged result");
            Assert.That(manyIter.AlignedSequences, Is.EqualTo(oneIter.AlignedSequences),
                "Converged alignment is identical regardless of remaining iteration budget");
        });
    }

    #endregion

    #region MUST — Determinism (no RNG; fixed iteration order)

    /// <summary>
    /// M12: DETERMINISM. Repeated runs on the same input produce byte-identical rows and the same
    /// SP score (no RNG; fixed edge order).
    /// </summary>
    [Test]
    public void MultipleAlignIterative_IsDeterministic()
    {
        var input = new[] { "GGAC", "TGT", "CTGGT", "TGG" }
            .Select(s => new DnaSequence(s)).ToArray();

        var r1 = SequenceAligner.MultipleAlignIterative(input);
        var r2 = SequenceAligner.MultipleAlignIterative(input);

        Assert.Multiple(() =>
        {
            Assert.That(r2.AlignedSequences, Is.EqualTo(r1.AlignedSequences));
            Assert.That(r2.TotalScore, Is.EqualTo(r1.TotalScore));
            Assert.That(r2.Consensus, Is.EqualTo(r1.Consensus));
        });
    }

    #endregion

    #region SHOULD — Custom scoring honored; sibling aligners unchanged

    /// <summary>
    /// S01: Custom scoring matrix is honored. Two identical 4-mers under BlastDna (match +2) score
    /// double the SimpleDna (match +1) result. No internal edge to refine, so this also confirms
    /// the seed pass-through carries the supplied matrix.
    /// SimpleDna: 1 pair × 4 matches × 1 = 4. BlastDna: 1 pair × 4 matches × 2 = 8.
    /// </summary>
    [Test]
    public void MultipleAlignIterative_WithCustomScoring_UsesProvidedMatrix()
    {
        var input = new[] { new DnaSequence("ATGC"), new DnaSequence("ATGC") };

        var simple = SequenceAligner.MultipleAlignIterative(input);
        var blast = SequenceAligner.MultipleAlignIterative(input, SequenceAligner.BlastDna);

        Assert.Multiple(() =>
        {
            Assert.That(simple.TotalScore, Is.EqualTo(4));
            Assert.That(blast.TotalScore, Is.EqualTo(8));
        });
    }

    /// <summary>
    /// S02: Adding the iterative aligner does not change the star or progressive aligners. On the
    /// discriminating two-cluster input the star MSA must still emit its documented length-6 gapped
    /// result and the progressive its gap-free length-5 result (byte-for-byte unchanged).
    /// </summary>
    [Test]
    public void SiblingAligners_Unchanged_OnDiscriminatingInput()
    {
        var input = new[]
        {
            new DnaSequence("AAGAA"),
            new DnaSequence("AACAA"),
            new DnaSequence("GGTGG"),
            new DnaSequence("GGTGG")
        };

        var star = SequenceAligner.MultipleAlign(input);
        var prog = SequenceAligner.MultipleAlignProgressive(input);

        Assert.Multiple(() =>
        {
            Assert.That(star.AlignedSequences,
                Is.EqualTo(new[] { "AAG-AA", "-AACAA", "-GGTGG", "-GGTGG" }),
                "Star MSA output must remain unchanged by the additive iterative aligner");
            Assert.That(prog.AlignedSequences,
                Is.EqualTo(new[] { "AAGAA", "AACAA", "GGTGG", "GGTGG" }),
                "Progressive MSA output must remain unchanged by the additive iterative aligner");
        });
    }

    #endregion

    #region SHOULD — Property-based: random panel never regresses, always degaps (O(n^2)+ invariant)

    /// <summary>
    /// S03: PROPERTY (O(k·L²)+ algorithm). Over many random inputs, the iterative SP is always
    /// ≥ the progressive seed SP (never a regression) AND every row degaps to its input and all
    /// rows are equal length. Deterministic fixed seed.
    /// </summary>
    [Test]
    public void MultipleAlignIterative_RandomInputs_NeverRegressAndStayValid()
    {
        var rng = new Random(20240623);
        char[] alpha = { 'A', 'C', 'G', 'T' };

        for (int trial = 0; trial < 500; trial++)
        {
            int k = rng.Next(3, 6);
            var raw = new string[k];
            for (int i = 0; i < k; i++)
            {
                int len = rng.Next(3, 9);
                var c = new char[len];
                for (int j = 0; j < len; j++) c[j] = alpha[rng.Next(4)];
                raw[i] = new string(c);
            }
            var seqs = raw.Select(s => new DnaSequence(s)).ToArray();

            int progSp = SequenceAligner.MultipleAlignProgressive(seqs).TotalScore;
            var iter = SequenceAligner.MultipleAlignIterative(seqs);

            Assert.That(iter.TotalScore, Is.GreaterThanOrEqualTo(progSp),
                $"Regression on trial {trial}: iter {iter.TotalScore} < prog {progSp} " +
                $"for [{string.Join(",", raw)}]");

            int len0 = iter.AlignedSequences[0].Length;
            for (int i = 0; i < k; i++)
            {
                Assert.That(iter.AlignedSequences[i].Length, Is.EqualTo(len0),
                    $"Unequal row length on trial {trial}");
                Assert.That(iter.AlignedSequences[i].Replace("-", ""), Is.EqualTo(raw[i]),
                    $"Row {i} does not degap to its input on trial {trial}");
            }
        }
    }

    #endregion
}
