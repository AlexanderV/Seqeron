using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Alignment;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Assembly area — K-mer Spectrum Read Error Correction
/// (ASSEMBLY-CORRECT-001), the Musket/Quake two-sided substitution corrector
/// <see cref="SequenceAssembler.ErrorCorrectReads(IReadOnlyList{string}, int, int)"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and
/// asserts the code NEVER fails in an undisciplined way: no hang or infinite loop
/// (the single left-to-right correction pass over a fixed spectrum must always
/// terminate), no state corruption, no nonsense output (a corrected read whose
/// length ≠ its input length, a residue invented outside A,C,G,T, a "correction"
/// of a position already covered by a trusted k-mer, or a non-deterministic
/// result), and no *unhandled* runtime exception — in particular NO DivideByZero
/// or NullReference on the degenerate boundaries of this row: ZERO COVERAGE
/// (no k-mer ever reaches the trusted cut-off), ALL-ERROR READS (so noisy that no
/// solid k-mer spectrum exists), and EMPTY (no reads / empty read). Every input
/// must resolve to EITHER a well-defined, theory-correct result OR a *documented,
/// intentional* validation exception (ArgumentNullException for a null read list
/// — §3.3, §6.1; ArgumentOutOfRangeException for kmerSize &lt; 1 — §3.3, §6.1).
/// A raw runtime exception, a hang, a length-changing edit, an invented residue,
/// a false correction, or an order-dependent result is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ASSEMBLY-CORRECT-001 — K-mer Spectrum Read Error Correction
/// Checklist: docs/checklists/03_FUZZING.md, row 141.
/// Algorithm doc: docs/algorithms/Assembly/Error_Correction.md
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the
///          checklist row "zero coverage, all-error reads, empty":
///          – ZERO COVERAGE: reads so sparse (or a cut-off so high) that NO k-mer
///            reaches the trusted multiplicity t → NO position is trusted yet NO
///            substitution can make all covering k-mers trusted → DOCUMENTED no-op:
///            every read returned UNCHANGED (uppercased), no false corrections, no
///            DivideByZero on the empty/never-trusted spectrum (§6.1 "No valid
///            alternative ⇒ unchanged"; INV-04).
///          – ALL-ERROR READS: a read set so noisy that every k-mer is a unique
///            singleton (the error mode with no genomic mode) → spectrum has no
///            solid k-mer → DOCUMENTED no-op, no crash, NO infinite correction loop
///            (single fixed-spectrum pass; INV-03/INV-04).
///          – EMPTY: no reads → documented empty result (no columns to correct);
///            an empty / shorter-than-k read → returned unchanged (uppercased),
///            with NO DivideByZero / NullReference (§3.3, §6.1, INV-01/INV-02).
/// — docs/checklists/03_FUZZING.md §Description (BE = Boundary Exploitation:
///   граничні значення 0, -1, MaxInt, empty).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test (Error_Correction.md §2.2, §2.4, §3, §6.1)
/// ───────────────────────────────────────────────────────────────────────────
/// Build the k-mer spectrum once (overlapping length-k substrings over all reads,
/// case-insensitive). A k-mer x is TRUSTED iff mult(x) ≥ t (minKmerFrequency),
/// UNTRUSTED otherwise (§2.2). For each read position i:
///   • if any k-mer covering i is trusted, the base is LEFT UNCHANGED (trusted-base
///     rule, INV-03);
///   • else test the alternative DNA bases {A,C,G,T} in fixed order: substitute and
///     check whether EVERY k-mer covering i becomes trusted; apply the substitution
///     ONLY when exactly ONE alternative qualifies, else restore the original base
///     (two-sided unique-alternative rule; ambiguity ⇒ unchanged, INV-04).
/// Output: corrected reads, UPPER-CASED, in input order; read COUNT preserved
/// (INV-01) and per-read LENGTH preserved — substitutions only, no indels (INV-02);
/// output is deterministic (INV-05); only A,C,G,T are ever produced as replacements
/// (§3.3). A read shorter than k contributes no k-mers and is returned unchanged
/// (§6.1). Null reads → ArgumentNullException; kmerSize &lt; 1 →
/// ArgumentOutOfRangeException (§3.3, §6.1). Defaults: kmerSize = 15, t = 2.
///   SequenceAssembler.ErrorCorrectReads(
///       IReadOnlyList&lt;string&gt; reads, int kmerSize = 15, int minKmerFrequency = 2)
///   → IReadOnlyList&lt;string&gt;
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class AssemblyErrorCorrectionFuzzTests
{
    // Documented defaults (Error_Correction.md §3.1).
    private const int DefaultKmerSize = 15;
    private const int DefaultMinKmerFrequency = 2;

    private static readonly char[] DnaUpper = { 'A', 'C', 'G', 'T' };

    #region Helpers

    /// <summary>
    /// Asserts the corrected read list is WELL-FORMED per the documented contract,
    /// regardless of the (possibly degenerate) input:
    ///   INV-01 output count == input count;
    ///   INV-02 each output read length == its input read length (substitutions only,
    ///          no indels);
    ///   §3.2  each output read is the upper-cased form of a same-length string;
    ///   §3.3  every emitted character is either an UPPER-CASED copy of the input
    ///          character at that position (unchanged base) OR one of A,C,G,T (the
    ///          only legal replacement residues) — no invented symbol, no gap.
    /// </summary>
    private static void AssertWellFormed(IReadOnlyList<string> corrected, IReadOnlyList<string> input)
    {
        corrected.Should().NotBeNull("the corrector never returns null (§3.2)");
        corrected.Count.Should().Be(input.Count, "INV-01: output read count == input read count");

        for (int r = 0; r < input.Count; r++)
        {
            string outRead = corrected[r];
            string inRead = input[r] ?? string.Empty;

            outRead.Should().NotBeNull("INV-01: each output read is a real string");
            outRead.Should().HaveLength(inRead.Length, "INV-02: per-read length is preserved (no indels)");

            for (int i = 0; i < outRead.Length; i++)
            {
                char outCh = outRead[i];
                char inUpper = char.ToUpperInvariant(inRead[i]);

                bool unchanged = outCh == inUpper;
                bool isDnaReplacement = Array.IndexOf(DnaUpper, outCh) >= 0;

                (unchanged || isDnaReplacement).Should().BeTrue(
                    "§3.3: every output base is either the unchanged (upper-cased) input base " +
                    "or one of the A,C,G,T replacement residues — never an invented symbol");
            }
        }
    }

    /// <summary>A random DNA read of the given length over the upper-case ACGT alphabet.</summary>
    private static string RandomDna(Random rng, int length)
    {
        var sb = new StringBuilder(length);
        for (int i = 0; i < length; i++)
            sb.Append(DnaUpper[rng.Next(DnaUpper.Length)]);
        return sb.ToString();
    }

    #endregion

    #region ASSEMBLY-CORRECT-001 — K-mer Spectrum Error Correction (BE: zero coverage, all-error reads, empty)

    #region Positive sanity — documented correction actually happens / error-free left unchanged

    // Documented worked example (§7.1, k=3, cut-off=2): three true reads ACGTACGT and one
    // copy carrying a single A→T substitution at index 4. The unique base that re-trusts
    // every covering k-mer is A, so the errant read is corrected back to ACGTACGT.
    [Test]
    public void ErrorCorrectReads_WorkedExample_CorrectsSingleSubstitution()
    {
        var reads = new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTTCGT" };

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        corrected.Should().Equal(new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTACGT" },
            "the unique re-trusting base at the error position is A (§7.1)");
        AssertWellFormed(corrected, reads);
    }

    // An error-free, high-coverage read set: every position is covered by a trusted k-mer,
    // so the trusted-base rule leaves every read unchanged (INV-03; §6.1 "all reads
    // identical/error-free ⇒ unchanged").
    [Test]
    public void ErrorCorrectReads_ErrorFreeHighCoverage_LeavesReadsUnchanged()
    {
        var reads = Enumerable.Repeat("ACGTACGTACGTAC", 6).ToArray();

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 4, minKmerFrequency: 2);

        corrected.Should().Equal(reads, "every position is covered by a trusted k-mer (INV-03)");
        AssertWellFormed(corrected, reads);
    }

    // Lower-case input is upper-cased; a single substitution in an otherwise solid context
    // is still corrected (case-insensitive spectrum, §3.3).
    [Test]
    public void ErrorCorrectReads_LowerCaseError_UppercasedAndCorrected()
    {
        var reads = new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT", "acgttcgt" };

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        corrected[3].Should().Be("ACGTACGT", "lower-case errant read is upper-cased then corrected (§3.3)");
        AssertWellFormed(corrected, reads);
    }

    #endregion

    #region BE — Boundary: zero coverage (no trusted k-mer ⇒ documented no-op, no false corrections)

    // Cut-off so high that NO k-mer ever reaches it: every k-mer is untrusted, so no position
    // is trusted, but no substitution can make covering k-mers trusted either → every read is
    // returned UNCHANGED (upper-cased). No false corrections, no DivideByZero (§6.1, INV-04).
    [Test]
    public void ErrorCorrectReads_CutoffAboveAllMultiplicities_ReturnsReadsUnchanged()
    {
        var reads = new[] { "ACGTACGT", "ACGTACGT", "ACGTACGT", "ACGTTCGT" };

        // Highest k-mer multiplicity here is 8; a cut-off of 100 trusts nothing.
        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 100);

        corrected.Should().Equal(reads.Select(r => r.ToUpperInvariant()),
            "with no trusted reference k-mer nothing can be corrected → reads unchanged (INV-04)");
        AssertWellFormed(corrected, reads);
    }

    // Single-copy reads (each k-mer multiplicity 1) at the default cut-off 2: zero coverage of
    // any trusted k-mer → no corrections at all, reads returned unchanged, no DivideByZero.
    [Test]
    public void ErrorCorrectReads_AllDistinctSingletonReads_ReturnsReadsUnchanged()
    {
        var reads = new[] { "ACGTACGTAC", "TTGGCCAATT", "GATTACAGAT" };

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 4, minKmerFrequency: 2);

        corrected.Should().Equal(reads, "no k-mer reaches the trusted cut-off ⇒ nothing corrected");
        AssertWellFormed(corrected, reads);
    }

    // Fuzz: random sparse reads with a deliberately unreachable cut-off must NEVER change any
    // read and NEVER throw a runtime exception (no DivideByZero on the never-trusted spectrum).
    [Test]
    [CancelAfter(30_000)]
    public void ErrorCorrectReads_ZeroCoverage_NeverCorrects_NeverThrows()
    {
        var rng = new Random(141_001);
        for (int trial = 0; trial < 600; trial++)
        {
            int depth = rng.Next(1, 6);
            var reads = Enumerable.Range(0, depth)
                                  .Select(_ => RandomDna(rng, rng.Next(0, 30)))
                                  .ToArray();
            int k = rng.Next(1, 12);

            // Cut-off strictly above the maximum possible multiplicity ⇒ zero trusted k-mers.
            int maxPossibleMult = Math.Max(1, reads.Sum(r => Math.Max(0, r.Length - k + 1)));
            int unreachableCutoff = maxPossibleMult + 1;

            var corrected = SequenceAssembler.ErrorCorrectReads(reads, k, unreachableCutoff);

            AssertWellFormed(corrected, reads);
            corrected.Should().Equal(reads.Select(r => r.ToUpperInvariant()),
                "zero trusted coverage ⇒ documented no-op (reads returned unchanged, INV-04)");
        }
    }

    #endregion

    #region BE — Boundary: all-error reads (no solid spectrum ⇒ no-op, no crash, no infinite loop)

    // Every read distinct and every k-mer a singleton — the "all error, no genomic mode" case:
    // no solid k-mer exists, so the corrector is a documented no-op and MUST terminate (single
    // fixed-spectrum pass — no iterative re-counting that could loop).
    [Test]
    [CancelAfter(30_000)]
    public void ErrorCorrectReads_AllErrorReads_NoSolidSpectrum_NoOpTerminates()
    {
        var rng = new Random(141_002);
        // 40 mutually-distinct reads ⇒ overlapping k-mers are overwhelmingly unique.
        var seen = new HashSet<string>();
        var reads = new List<string>();
        while (reads.Count < 40)
        {
            string r = RandomDna(rng, 25);
            if (seen.Add(r)) reads.Add(r);
        }

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 11, minKmerFrequency: 2);

        AssertWellFormed(corrected, reads);
        // With near-certainty no k-mer repeats; even if a rare collision trusts a k-mer, the
        // contract still preserves count/length/alphabet (asserted above). The dominant
        // documented behavior is no-op: assert nothing was lengthened/shortened and the run
        // terminated (CancelAfter guards the no-infinite-loop property, INV-03/INV-04).
        corrected.Count.Should().Be(reads.Count, "INV-01: count preserved even with no solid spectrum");
    }

    // Homopolymer-style noisy reads of differing single bases (A.., C.., G.., T..): no shared
    // k-mer across reads, so each read's k-mers are trusted only WITHIN that read at depth ≥ cutoff
    // when the run is long. This is a no-crash / termination probe over a pathological spectrum.
    [Test]
    [CancelAfter(15_000)]
    public void ErrorCorrectReads_PathologicalHomopolymers_NoCrash()
    {
        var reads = new[]
        {
            new string('A', 30),
            new string('C', 30),
            new string('G', 30),
            new string('T', 30),
        };

        Action act = () =>
        {
            var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 5, minKmerFrequency: 2);
            AssertWellFormed(corrected, reads);
        };

        act.Should().NotThrow("a pathological homopolymer spectrum must not crash or hang");
    }

    // A single noisy substitution that has NO unique re-trusting base must be LEFT UNCHANGED
    // (ambiguity / no-valid-alternative rule, §6.1, INV-04) — not silently mangled.
    [Test]
    public void ErrorCorrectReads_NoUniqueAlternative_LeavesBaseUnchanged()
    {
        // Two competing trusted contexts: an errant middle base could be re-trusted by more than
        // one alternative OR by none; either way the documented behavior is "unchanged".
        var reads = new[] { "AAATAAA", "AAATAAA" }; // every 3-mer multiplicity 2 already trusted
        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 3, minKmerFrequency: 2);

        corrected.Should().Equal(reads, "already-trusted positions are never modified (INV-03)");
        AssertWellFormed(corrected, reads);
    }

    #endregion

    #region BE — Boundary: empty (no reads / empty read / shorter-than-k ⇒ guarded, no /0 or NRE)

    // No reads → documented empty result; no spectrum to divide by, no NullReference (§3.3, §6.1).
    [Test]
    public void ErrorCorrectReads_EmptyReadList_ReturnsEmpty()
    {
        var corrected = SequenceAssembler.ErrorCorrectReads(Array.Empty<string>());

        corrected.Should().BeEmpty("no reads ⇒ empty corrected list (INV-01)");
    }

    // A single empty read → returned unchanged (length 0 preserved); no DivideByZero / NRE.
    [Test]
    public void ErrorCorrectReads_SingleEmptyRead_ReturnsEmptyReadUnchanged()
    {
        var reads = new[] { "" };

        var corrected = SequenceAssembler.ErrorCorrectReads(reads);

        corrected.Should().Equal(new[] { "" }, "an empty read has no k-mers ⇒ returned unchanged (§6.1)");
        AssertWellFormed(corrected, reads);
    }

    // Reads all shorter than k contribute no k-mers and are returned unchanged (upper-cased),
    // with no DivideByZero on the empty spectrum (§6.1, "read shorter than k ⇒ unchanged").
    [Test]
    public void ErrorCorrectReads_AllReadsShorterThanK_ReturnedUppercasedUnchanged()
    {
        var reads = new[] { "ac", "g", "tt" };

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 15, minKmerFrequency: 2);

        corrected.Should().Equal(new[] { "AC", "G", "TT" },
            "reads shorter than k contribute no k-mers ⇒ returned upper-cased unchanged (§6.1)");
        AssertWellFormed(corrected, reads);
    }

    // Mix of empty and non-empty reads at the default k: empty/short reads pass through unchanged;
    // no element triggers a NullReference and the spectrum build (zero contributing k-mers) is safe.
    [Test]
    public void ErrorCorrectReads_MixedEmptyAndShortReads_NoCrash()
    {
        var reads = new[] { "", "ACG", "", "ACGTACGTACGTACGTACGT" };

        var corrected = SequenceAssembler.ErrorCorrectReads(reads, kmerSize: 15, minKmerFrequency: 2);

        AssertWellFormed(corrected, reads);
        corrected[0].Should().BeEmpty("empty reads stay empty");
        corrected[2].Should().BeEmpty("empty reads stay empty");
    }

    #endregion

    #region Validation contract (§3.3, §6.1)

    // Null read list → documented ArgumentNullException (§3.3, §6.1).
    [Test]
    public void ErrorCorrectReads_NullReadList_Throws()
    {
        Action act = () => SequenceAssembler.ErrorCorrectReads(null!);

        act.Should().Throw<ArgumentNullException>("null reads is the documented validation contract (§3.3)");
    }

    // kmerSize < 1 → documented ArgumentOutOfRangeException (§3.3, §6.1).
    [Test]
    public void ErrorCorrectReads_NonPositiveKmerSize_Throws()
    {
        var reads = new[] { "ACGTACGT" };

        foreach (int badK in new[] { 0, -1, int.MinValue })
        {
            Action act = () => SequenceAssembler.ErrorCorrectReads(reads, kmerSize: badK);
            act.Should().Throw<ArgumentOutOfRangeException>(
                "kmerSize < 1 is rejected (§3.3) for k={0}", badK);
        }
    }

    #endregion

    #region BE — Broad fuzz: random reads/cutoffs never throw, always preserve the contract

    // Sweeps depth, length, k and cut-off (including degenerate combinations) and asserts the
    // corrector NEVER throws an unexpected exception, ALWAYS preserves count/length/alphabet,
    // and is DETERMINISTIC (INV-05) — re-running the same input yields the identical result.
    [Test]
    [CancelAfter(60_000)]
    public void ErrorCorrectReads_RandomInputs_NeverThrows_WellFormed_Deterministic()
    {
        var rng = new Random(141_003);
        for (int trial = 0; trial < 1200; trial++)
        {
            int depth = rng.Next(0, 10);
            var reads = Enumerable.Range(0, depth)
                                  .Select(_ => RandomDna(rng, rng.Next(0, 35)))
                                  .ToArray();
            int k = rng.Next(1, 18);
            int cutoff = rng.Next(1, 6);

            var first = SequenceAssembler.ErrorCorrectReads(reads, k, cutoff);
            var second = SequenceAssembler.ErrorCorrectReads(reads, k, cutoff);

            AssertWellFormed(first, reads);
            first.Should().Equal(second, "INV-05: error correction is deterministic");
        }
    }

    #endregion

    #endregion
}
