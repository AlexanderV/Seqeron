using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Metagenomics;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Metagenomics TETRA tetranucleotide z-score signature
/// (Teeling et al. 2004, BMC Bioinformatics 5:163; Environ Microbiol 6(9):938–947).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, malformed and boundary inputs to a unit and asserts the
/// code NEVER fails in an undisciplined way: no hang, no state corruption, no
/// *unhandled* runtime exception (DivideByZeroException, IndexOutOfRange, …) and
/// — critically for a numeric signature — no *out-of-contract output* (a NaN or
/// ±Inf z-score, a result that does not cover all 256 ACGT tetranucleotides, or a
/// non-zero z where the documented theory says it must be 0). Every input must
/// yield EITHER a well-defined, theory-correct result OR a documented validation
/// exception. — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-TETRA-001 — TETRA z-score signature (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 249. Strategies: BE, MC.
///   • BE = Boundary Exploitation — sequence &lt; 4 bp, empty sequence, single
///          base, all-same-base sequence, very long sequence.
///   • MC = Malformed Content — non-ACGT characters (IUPAC N/R/Y, digits, gaps,
///          unicode, lowercase), interleaved junk.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
///
/// Entry points under test:
///   • MetagenomicsAnalyzer.CalculateTetranucleotideZScores(string)  — 256-entry
///     ACGT tetranucleotide → Teeling z-score map.
///   • MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(string, string) —
///     Pearson r ∈ [−1, 1] of two signatures.
///   — src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs
///     (CalculateTetranucleotideZScores ~675, TetranucleotideZScore ~725,
///      ExtendWithReverseComplement ~759, CountOligonucleotides ~792,
///      TetranucleotideZScoreCorrelation ~707, PearsonCorrelation ~826);
///     docs/algorithms/Metagenomics/Genome_Binning.md.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The TETRA contract these fuzz tests pin (from the method XML docs + theory)
/// ───────────────────────────────────────────────────────────────────────────
/// For each of the 256 ACGT tetranucleotides n1n2n3n4 the algorithm computes a
/// maximal-order (2nd-order) Markov z-score over the sequence EXTENDED by its
/// reverse complement (strand-symmetric, Teeling):
///   E  = N(n1n2n3)·N(n2n3n4) / N(n2n3)                         (expected count)
///   var = E·[N(n2n3)−N(n1n2n3)]·[N(n2n3)−N(n2n3n4)] / N(n2n3)² (Schbath 1997)
///   z  = (N(n1n2n3n4) − E) / √var.
/// Documented degenerate fallbacks (TetranucleotideZScore, lines 736–749):
///   • N(n2n3) = 0           ⇒ z = 0   (zero-denominator guard);
///   • var ≤ 0               ⇒ z = 0   (no over/under-representation evidence).
/// These two guards are exactly what make EVERY returned z FINITE — there is no
/// path to a 0/0 or x/0 NaN/Inf. That finiteness is the HEADLINE invariant here.
///
/// Boundary / malformed handling fixed by the doc + source, pinned below:
///   • null / empty / single base / &lt; the bases needed for a 4-window →
///     ExtendWithReverseComplement returns "" (null/empty) or a too-short string,
///     so no tetranucleotide is ever counted ⇒ an all-zero 256-entry map. NOT a
///     crash, NOT an empty/partial map. (XML <returns>: "null, empty, or
///     single-base input produces an all-zero 256-entry map".)
///   • ALL-SAME-base sequence (e.g. "AAAA…") → extended strand is a run of A then
///     a run of T; only AAAA (and TTTT) ever occur, and for them prefix3 ==
///     middle2 so [N(n2n3)−N(n1n2n3)] = 0 ⇒ var = 0 ⇒ z = 0. So an all-same-base
///     input yields an ALL-ZERO signature with NO NaN — the canonical z=0 case.
///   • non-ACGT characters → ExtendWithReverseComplement drops them before
///     extension and CountOligonucleotides skips any window containing one
///     (Increment's IsAcgtOnly guard, line 808). A sequence made ENTIRELY of
///     non-ACGT symbols therefore behaves exactly like the empty case ⇒ all-zero.
///     Handled, never rejected, never crashed.
///   • very long sequence → the 256 int bins cannot overflow at realistic sizes
///     and the single-pass counter is O(n); it must COMPLETE under a time budget
///     with every z finite.
///
/// Correlation contract (TetranucleotideZScoreCorrelation + PearsonCorrelation):
///   • r ∈ [−1, 1] ALWAYS (Pearson bound), and FINITE — the `denom > 0 ? … : 0`
///     guard (line 842) maps a constant (all-zero) vector to r = 0, never 0/0.
///   • r(x, x) = 1 for a sequence with a NON-constant signature (self-correlation
///     of a non-degenerate vector). When the signature is all-zero (constant), the
///     documented value is r = 0, not 1 — pinned explicitly.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Anti-rubber-stamp positive sanity
/// ───────────────────────────────────────────────────────────────────────────
/// A real, compositionally-skewed sequence must produce a signature with at least
/// one non-zero (and finite) z-score, and its self-correlation must be exactly 1.
/// So a passing "no NaN" result cannot be a degenerate implementation that returns
/// an all-zero map (or r = 0) for every input.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism &amp; policy
/// ───────────────────────────────────────────────────────────────────────────
/// Every input is hand-built or drawn from a LOCALLY fixed-seed `new Random(seed)`
/// (no shared static Rng). CalculateTetranucleotideZScores / …Correlation are not
/// LimitationPolicy-guarded; the test assembly's module-initializer bootstrap sets
/// DefaultMode = Permissive regardless, so the real numeric path is exercised.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MetagenomicsTetraFuzzTests
{
    // The 256 ACGT tetranucleotides, the documented full key set of every signature.
    private static readonly string[] AllTetranucleotides =
        (from a in "ACGT"
         from b in "ACGT"
         from c in "ACGT"
         from d in "ACGT"
         select new string(new[] { a, b, c, d })).ToArray();

    private static void AssertFullFiniteSignature(IReadOnlyDictionary<string, double> z)
    {
        z.Should().HaveCount(256, "the signature must cover all 4^4 ACGT tetranucleotides");
        z.Keys.Should().BeEquivalentTo(AllTetranucleotides,
            "exactly the 256 ACGT tetranucleotides — no missing or extra key");
        z.Values.Should().OnlyContain(
            v => !double.IsNaN(v) && !double.IsInfinity(v),
            "EVERY z-score must be finite — the zero-denominator/zero-variance guards "
            + "forbid any NaN or ±Inf (HEADLINE invariant)");
    }

    #region META-TETRA-001 — TETRA z-score signature

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: NULL / EMPTY / SINGLE USABLE BASE (BE).
    // ExtendWithReverseComplement yields "" (null/empty) or a length-2 strand
    // (one usable base + its complement), so NO 4-window is ever counted ⇒ an
    // all-zero 256-entry map. KEY: a full, finite, all-zero map — never a crash,
    // never a partial map. NB: per the XML, only a SINGLE usable base is all-zero;
    // ≥2 usable ACGT bases already give a ≥4-nt extended strand (see next test).
    // — XML <returns>; ExtendWithReverseComplement / CountOligonucleotides.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void TetraZScores_NullEmptyOrSingleUsableBase_FullAllZeroSignature()
    {
        // null/empty → "" ; single usable ACGT base (incl. lowercase or junk-padded
        // to a single base) → length-2 extended strand → no tetranucleotide window.
        foreach (string? seq in new[] { null, "", "A", "C", "G", "T", "a", "N-A-N", "9G9" })
        {
            IReadOnlyDictionary<string, double> z = default!;
            Action act = () => z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq!);

            act.Should().NotThrow(
                $"a null/empty/single-usable-base sequence ({(seq is null ? "null" : $"\"{seq}\"")}) "
                + "is a defined boundary, not a crash");

            AssertFullFiniteSignature(z);
            z.Values.Should().OnlyContain(v => v == 0.0,
                "≤1 usable base ⇒ extended strand < 4 nt ⇒ no window counted ⇒ every z = 0");
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SHORTER-THAN-4 but ≥2 usable bases (BE) — finiteness boundary.
    // Per the XML, ≥2 usable ACGT bases already yield a ≥4-nt reverse-complement-
    // extended strand on which tetranucleotides CAN be counted. So these do NOT
    // give an all-zero map; the contract that must hold is the HEADLINE one: a
    // full 256-key map of FINITE z-scores (no NaN/Inf), with no crash.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void TetraZScores_TwoOrThreeUsableBases_FullFiniteSignature()
    {
        foreach (string seq in new[] { "AC", "ACG", "gt", "Ga", "tca", "A-C", "n c g" })
        {
            IReadOnlyDictionary<string, double> z = default!;
            Action act = () => z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);

            act.Should().NotThrow($"a short (\"{seq}\") sequence must not crash");
            AssertFullFiniteSignature(z);
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: ALL-SAME-BASE sequence (BE) — the zero-denominator NaN trap.
    // "AAAA…" extends with its reverse complement to a run of A then a run of T
    // ("AAAA…TTTT…"). Only a handful of tetranucleotides ever occur (AAAA, the A→T
    // junction words, TTTT); the OTHER ~250 never occur, so their N(n2n3) middle-
    // dimer count is 0 ⇒ the zero-denominator guard forces z = 0 (NO x/0 Inf, NO
    // 0/0 NaN — the HEADLINE invariant). The few that occur get a FINITE z. So the
    // contract: a full 256-key map of FINITE z-scores, and a tetranucleotide that
    // CANNOT occur in a homopolymer (e.g. one using a different base) is exactly 0.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void TetraZScores_AllSameBase_FiniteWithAbsentTetraAtZero()
    {
        foreach (char b in new[] { 'A', 'C', 'G', 'T' })
        {
            // a tetranucleotide that can NEVER occur in a {b}/{complement(b)} run.
            char other = b == 'G' ? 'A' : 'G';      // distinct base & complement
            string absent = new string(new[] { b, other, b, other });

            foreach (int len in new[] { 4, 10, 64, 500 })
            {
                string seq = new string(b, len);

                var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);

                AssertFullFiniteSignature(z);
                z[absent].Should().Be(0.0,
                    $"a tetranucleotide ('{absent}') absent from an all-'{b}' run (len {len}) "
                    + "has N(n2n3) = 0 ⇒ z forced to 0, NO NaN/Inf");
            }
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: NON-ACGT-ONLY content (MC).
    // Every character is dropped (extension) / every window skipped (counting) ⇒
    // identical to the empty case: a full, finite, all-zero signature. Handled,
    // not rejected. — ExtendWithReverseComplement filter; Increment IsAcgtOnly.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void TetraZScores_NonAcgtOnly_HandledAsAllZero()
    {
        foreach (string seq in new[]
                 {
                     "NNNNNNNNNNNN",       // IUPAC N
                     "NRYWSKMBDHVN",       // IUPAC ambiguity codes
                     "1234567890123456",   // digits
                     "----...----...",     // gaps / dots
                     "ΑΒΓΔαβγδ★☢ ΑΒΓΔ",     // unicode + spaces
                     "nnnnnnnnnnnn",       // lowercase n
                 })
        {
            IReadOnlyDictionary<string, double> z = default!;
            Action act = () => z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);

            act.Should().NotThrow($"non-ACGT content (\"{seq}\") is filtered, not a crash");

            AssertFullFiniteSignature(z);
            z.Values.Should().OnlyContain(v => v == 0.0,
                "with no usable ACGT base no tetranucleotide is counted ⇒ all z = 0");
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: case-insensitivity + interleaved junk (MC).
    // Lowercase ACGT must be upper-cased and counted; interleaved non-ACGT must be
    // dropped without breaking the window for the valid bases around it. The
    // signature must equal the one for the same sequence with the junk removed.
    // — XML: "case-insensitive; non-ACGT characters are skipped".
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void TetraZScores_LowercaseAndInterleavedJunk_MatchCleanedSequence()
    {
        // Note: non-ACGT is DROPPED (not a gap), so the cleaned sequence is the
        // ACGT subsequence in order. Build a clean sequence and a polluted copy.
        const string clean = "ACGTTGCAACGTTGCAACGTACGTGGCCAATTACGT";
        var polluted = new StringBuilder();
        var rng = new Random(20260626);
        foreach (char c in clean)
        {
            polluted.Append(char.ToLowerInvariant(c)); // lowercase the real base
            if (rng.Next(3) == 0) polluted.Append("NRY-9 ★"[rng.Next(7)]); // junk
        }

        var zClean = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(clean);
        var zPolluted = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(polluted.ToString());

        AssertFullFiniteSignature(zClean);
        AssertFullFiniteSignature(zPolluted);

        foreach (string tetra in AllTetranucleotides)
            zPolluted[tetra].Should().BeApproximately(zClean[tetra], 1e-12,
                "lowercasing + dropped junk must reproduce the cleaned-sequence signature");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random ACGT/junk batch (BE, MC) — finiteness everywhere.
    // A deterministic, locally-seeded generator produces sequences mixing ACGT,
    // lowercase, IUPAC codes, digits and gaps over a span of lengths (0 .. long).
    // EVERY produced signature must be a full 256-key map of FINITE z-scores —
    // the headline contract — for every malformed/boundary input.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void TetraZScores_RandomMalformedBatch_AlwaysFull256FiniteSignature()
    {
        var rng = new Random(20260626); // locally fixed seed — deterministic
        const string alphabet = "ACGTacgtNRYWSKnX-.0123456789 \t★Α";

        for (int i = 0; i < 400; i++)
        {
            int len = rng.Next(0, 40); // includes 0, < 4, and ≥ 4 lengths
            var chars = new char[len];
            for (int j = 0; j < len; j++)
                chars[j] = alphabet[rng.Next(alphabet.Length)];
            string seq = new string(chars);

            IReadOnlyDictionary<string, double> z = default!;
            Action act = () => z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);
            act.Should().NotThrow($"input #{i} (\"{seq}\") must not crash");

            AssertFullFiniteSignature(z);
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: VERY LONG sequence (BE) — no overflow / hang.
    // A long deterministic ACGT string must be counted in one O(n) pass with no
    // int-bin overflow and every z finite, within a time budget.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void TetraZScores_VeryLongSequence_CompletesWithFiniteSignature()
    {
        var rng = new Random(20260626);
        var sb = new StringBuilder(200_000);
        for (int i = 0; i < 200_000; i++)
            sb.Append("ACGT"[rng.Next(4)]);

        IReadOnlyDictionary<string, double> z = default!;
        Action act = () => z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(sb.ToString());
        act.Should().NotThrow("a 200k-bp sequence must not overflow the 256 int bins or hang");

        AssertFullFiniteSignature(z);
        z.Values.Should().Contain(v => v != 0.0,
            "a long random sequence has real over/under-representation ⇒ some non-zero z");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Positive sanity: a compositionally-skewed sequence must produce a
    // NON-degenerate, finite signature (at least one non-zero z). Guards against
    // an implementation that returns an all-zero map for everything (which would
    // pass every boundary test above).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void TetraZScores_SkewedSequence_HasNonZeroFiniteZScore()
    {
        // GC-rich, CpG-biased composition with strong tetranucleotide skew.
        const string skewed =
            "GCGCGCGCATATATATGCGCGCGCATATATATGGCCGGCCAATTAATTGCGCGCGC";

        var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(skewed);

        AssertFullFiniteSignature(z);
        z.Values.Should().Contain(v => v != 0.0,
            "a skewed composition must yield over/under-represented tetranucleotides");
    }

    #endregion

    #region META-TETRA-001 — z-score correlation

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: correlation involving a CONSTANT (all-zero) signature (BE).
    // null / empty / single-usable-base / non-ACGT-only inputs give an all-zero
    // (constant) 256-vector. Pearson of a constant vector hits the `denom > 0 ?
    // … : 0` guard ⇒ the documented r = 0, never a 0/0 NaN — regardless of the
    // other operand. Pinned for every pairing of these degenerate inputs with both
    // each other and a real signature. (NB: "AAAA" is NOT all-zero — its RC
    // extension breaks the homopolymer — so it is excluded here; see the next test.)
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void Correlation_WithConstantSignature_IsFiniteZeroNotNaN()
    {
        var allZero = new[] { "", "A", "G", "NNNN", "ΑΒΓΔ", "----", "123", null };
        const string real = "GCGCGCGCATATATATGGCCAATTACGTACGT"; // non-constant signature

        foreach (string? deg in allZero)
        {
            // sanity: each listed input really does have an all-zero signature.
            var zd = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(deg!);
            zd.Values.Should().OnlyContain(v => v == 0.0,
                "this input must have a constant (all-zero) signature for the r = 0 contract");

            foreach (string? other in new[] { deg, real, null })
            {
                double r1 = default, r2 = default;
                Action act = () =>
                {
                    r1 = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(deg!, other!);
                    r2 = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(other!, deg!);
                };
                act.Should().NotThrow("a constant signature must not crash the correlation");

                r1.Should().Be(0.0,
                    "Pearson with a zero-variance signature is the documented r = 0, never NaN");
                r2.Should().Be(0.0, "and symmetrically r = 0 with the operands swapped");
            }
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: self-correlation of a NON-degenerate signature (BE).
    // r(x, x) of a sequence with a non-constant signature is exactly 1 (Pearson
    // of a vector with itself). Pinned across skewed sequences.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void Correlation_SelfCorrelationOfRealSignature_IsExactlyOne()
    {
        foreach (string seq in new[]
                 {
                     "GCGCGCGCATATATATGCGCGCGCATATATATGGCCGGCCAATT",
                     "ACGTACGTTGCATGCAACGTACGTTTTTAAAACCCCGGGGACGT",
                     "AACCGGTTAACCGGTTACGTACGTGCATGCATAATTCCGGAACC",
                 })
        {
            // sanity: the signature is genuinely non-constant (so self-corr = 1
            // is a real test, not the r = 0 degenerate path).
            var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);
            z.Values.Should().Contain(v => v != 0.0, "the signature must be non-degenerate");

            double r = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, seq);
            r.Should().BeApproximately(1.0, 1e-9,
                "Pearson correlation of a non-constant signature with itself is 1");
        }
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: correlation always in [−1, 1] over a random batch (BE, MC).
    // For a deterministic, locally-seeded batch of arbitrary (incl. malformed)
    // sequence pairs the correlation must ALWAYS be finite and within [−1, 1] —
    // the Pearson bound — and symmetric: r(a, b) = r(b, a).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void Correlation_RandomPairs_AlwaysFiniteInRangeAndSymmetric()
    {
        var rng = new Random(20260626);
        const string alphabet = "ACGTacgtNRY-9 ★";

        string RandomSeq()
        {
            int len = rng.Next(0, 50);
            var chars = new char[len];
            for (int j = 0; j < len; j++)
                chars[j] = alphabet[rng.Next(alphabet.Length)];
            return new string(chars);
        }

        for (int i = 0; i < 200; i++)
        {
            string a = RandomSeq();
            string b = RandomSeq();

            double rab = default, rba = default;
            Action act = () =>
            {
                rab = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(a, b);
                rba = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(b, a);
            };
            act.Should().NotThrow($"pair #{i} must not crash the correlation");

            rab.Should().Match(v => !double.IsNaN(v) && !double.IsInfinity(v),
                "correlation is always finite");
            rab.Should().BeInRange(-1.0, 1.0, "Pearson r ∈ [−1, 1] always");
            rba.Should().BeApproximately(rab, 1e-12,
                "Pearson correlation is symmetric: r(a, b) = r(b, a)");
        }
    }

    #endregion
}
