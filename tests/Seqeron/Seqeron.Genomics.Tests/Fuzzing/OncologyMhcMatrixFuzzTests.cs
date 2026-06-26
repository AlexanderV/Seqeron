using System;
using System.Collections.Generic;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;
using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology position-specific scoring-matrix (PSSM) peptide–MHC binding predictors —
/// MHC-MATRIX-001 (ONCO-MHC-001, matrix-based / opt-in half of the unit).
/// The units under test are
/// <see cref="OncologyAnalyzer.PredictBindingHalfLifeBimas(string, OncologyAnalyzer.PmhcScoringMatrix)"/>
/// (BIMAS / Parker 1994 product rule → predicted half-time of dissociation),
/// <see cref="OncologyAnalyzer.PredictIc50Smm(string, OncologyAnalyzer.PmhcScoringMatrix)"/>
/// (SMM / Peters &amp; Sette 2005 additive log50k → <c>IC50 = 50000^(1−score)</c>), and the
/// <see cref="OncologyAnalyzer.PredictAndClassifySmm(string, OncologyAnalyzer.PmhcScoringMatrix)"/>
/// predict→classify wrapper, all in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs (the
/// <see cref="OncologyAnalyzer.PmhcScoringMatrix"/> record + the three predictors, ~lines 8330–8552).
///
/// This file is scoped STRICTLY to the MATRIX predictors. The threshold classifiers
/// (<see cref="OncologyAnalyzer.ClassifyBindingAffinity(double)"/>,
/// <see cref="OncologyAnalyzer.ClassifyBindingRank(double, OncologyAnalyzer.MhcClass)"/>,
/// <see cref="OncologyAnalyzer.IsValidPeptideLength(int, OncologyAnalyzer.MhcClass)"/>,
/// <see cref="OncologyAnalyzer.ClassifyMhcBinding(int, double, OncologyAnalyzer.MhcClass)"/>) are
/// covered by the sibling <c>OncologyMhcBindingFuzzTests</c> (checklist row 110), and the MHCflurry
/// neural predictor by <c>OncologyMhcNnFuzzTests</c>; this file does not re-test those.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed / boundary inputs and asserts the matrix predictors never fail in an
/// undisciplined way: no IndexOutOfRangeException when the peptide length ≠ the matrix row count, no
/// KeyNotFoundException leaking when a residue is absent from a position's row, no NaN / Infinity IC50
/// reaching downstream, no hang. Every input must resolve to EITHER a well-defined, THEORY-correct result
/// OR a *documented, intentional* validation exception (INV-07: <see cref="ArgumentException"/> for a
/// length mismatch / empty matrix, <see cref="ArgumentNullException"/> for a null peptide).
///
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzz strategies (docs/checklists/03_FUZZING.md §Description), row 246
/// ───────────────────────────────────────────────────────────────────────────
/// MC = Malformed Content, BE = Boundary Exploitation. Checklist targets for this row:
///   • "peptide length ≠ matrix length" (BE) — must be a DOCUMENTED <see cref="ArgumentException"/>
///     (INV-07), never an IndexOutOfRangeException from over-/under-indexing <c>matrix.Rows[i]</c>.
///   • "missing residue" (MC) — a residue not listed in a position's row is DEFINED behaviour, not a
///     crash: BIMAS contributes the neutral coefficient 1.0 (multiplicative identity), SMM contributes 0
///     (additive identity) (INV-06 / §2.2). Must NOT throw KeyNotFound.
///   • "empty matrix" (BE) — 0 position rows → DOCUMENTED <see cref="ArgumentException"/> (INV-07), never
///     a silent empty-product (1·const) / empty-sum (50000^(1−intercept)) result.
///   • "non-AA" (MC) — non-amino-acid characters in the peptide (digits, punctuation, whitespace, unicode)
///     are simply residues absent from every row → treated as the neutral contribution, never a crash.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (independently re-derived, NOT read off the code)
/// ───────────────────────────────────────────────────────────────────────────
/// MHC_Peptide_Binding_Classification.md (docs/algorithms/Oncology/MHC_Peptide_Binding_Classification.md):
///   • §2.2 / INV-06: BIMAS T½ = FinalConstant · ∏_i coeff[i][peptide[i]]; running score starts at 1.0;
///     an unlisted residue contributes the neutral coefficient 1.0. Higher T½ = stronger binder.
///   • §2.2 / INV-05: SMM score = intercept + Σ_i contrib[i][peptide[i]] (unlisted residue → 0);
///     IC50 = 50000^(1 − score), strictly decreasing in score, always finite &amp; &gt; 0
///     (score 0 → 50000 nM, score 1 → 1 nM, score 0.5 → √50000). Lower IC50 = stronger binder.
///   • §7.1 worked examples (HAND-DERIVED here, asserted by equality):
///       BIMAS: CONST=10, coeffs L=2.0, M=3.0, V=1.5 for "LMV" ⇒ T½ = 10·(2·3·1.5) = 90.
///       SMM: contributions for a peptide summing to score 1.0 (intercept 0) ⇒ IC50 = 50000^0 = 1 nM → Strong;
///            a peptide matching NO listed residue ⇒ score 0 ⇒ IC50 = 50000^1 = 50000 nM → NonBinder.
///   • INV-07 / §3.3: the predictors require a non-empty matrix and peptide.Length == matrix.Rows.Count,
///     else ArgumentException; a null peptide throws ArgumentNullException.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng. The assembly bootstrap
/// (_LimitationPolicyTestBootstrap) runs under Permissive, so the ONCO-MHC-001 LimitationPolicy guard in
/// the predictors is a no-op here (it is Moderate-gated and only throws under Strict regardless).
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyMhcMatrixFuzzTests
{
    // ── Hand-built matrix factories (built from the doc's DEFINITION, never from code-read values) ──────

    // One position row from RESIDUE→value pairs (residues stored upper-case, as the predictors upper-case
    // the peptide before lookup).
    private static IReadOnlyDictionary<char, double> Row(params (char Residue, double Value)[] entries)
    {
        var d = new Dictionary<char, double>(entries.Length);
        foreach (var (r, v) in entries)
        {
            d[char.ToUpperInvariant(r)] = v;
        }
        return d;
    }

    private static PmhcScoringMatrix Matrix(double finalConstant,
        params IReadOnlyDictionary<char, double>[] rows)
        => new PmhcScoringMatrix(rows, finalConstant);

    // BIMAS product oracle re-derived from §2.2 / INV-06, independent of the implementation's loop.
    private static double ExpectedBimas(string peptide, PmhcScoringMatrix matrix)
    {
        double score = 1.0;
        for (int i = 0; i < peptide.Length; i++)
        {
            char r = char.ToUpperInvariant(peptide[i]);
            score *= matrix.Rows[i].TryGetValue(r, out double c) ? c : 1.0; // neutral 1.0
        }
        return score * matrix.FinalConstant;
    }

    // SMM IC50 oracle re-derived from §2.2 / INV-05, independent of the implementation's loop.
    private static double ExpectedSmmIc50(string peptide, PmhcScoringMatrix matrix)
    {
        double score = matrix.FinalConstant; // intercept
        for (int i = 0; i < peptide.Length; i++)
        {
            char r = char.ToUpperInvariant(peptide[i]);
            score += matrix.Rows[i].TryGetValue(r, out double v) ? v : 0.0; // neutral 0.0
        }
        return Math.Pow(50000.0, 1.0 - score);
    }

    #region MHC-MATRIX-001 — Hand-derived BIMAS / SMM worked-example equality pins (§7.1)

    // §7.1 verbatim: "A BIMAS matrix with constant 10 and coefficients 2.0, 3.0, 1.5 for LMV gives
    // T½ = 10 · (2·3·1.5) = 90." Re-derived by hand (2·3·1.5 = 9; ·10 = 90), asserted by exact equality.
    [Test]
    public void Bimas_DocWorkedExample_LMV_EqualsNinety()
    {
        var matrix = Matrix(10.0,
            Row(('L', 2.0)),
            Row(('M', 3.0)),
            Row(('V', 1.5)));

        PredictBindingHalfLifeBimas("LMV", matrix).Should().Be(90.0);
    }

    // §2.2 / INV-06: an unlisted residue contributes the neutral coefficient 1.0. "LXV" (X absent from the
    // middle row, which lists only M) ⇒ T½ = 10 · (2 · 1.0 · 1.5) = 30 — hand-derived, exact equality.
    [Test]
    public void Bimas_MissingResidue_UsesNeutralCoefficientOne()
    {
        var matrix = Matrix(10.0,
            Row(('L', 2.0)),
            Row(('M', 3.0)),   // 'X' absent ⇒ neutral 1.0
            Row(('V', 1.5)));

        // 10 · (2 · 1.0 · 1.5) = 30.
        PredictBindingHalfLifeBimas("LXV", matrix).Should().Be(30.0);
    }

    // §7.1 / §2.2: SMM IC50 = 50000^(1 − score). A peptide whose contributions sum to score 1.0 with
    // intercept 0 ⇒ IC50 = 50000^0 = 1 nM → Strong. Hand-built so the three contributions sum to 1.0.
    [Test]
    public void Smm_ContributionsSumToOne_GiveIc50OneNm_Strong()
    {
        var matrix = Matrix(0.0, // intercept
            Row(('A', 0.5)),
            Row(('C', 0.3)),
            Row(('D', 0.2))); // 0.5 + 0.3 + 0.2 = 1.0

        PredictIc50Smm("ACD", matrix).Should().BeApproximately(1.0, 1e-9);

        var (ic50, strength) = PredictAndClassifySmm("ACD", matrix);
        ic50.Should().BeApproximately(1.0, 1e-9);
        strength.Should().Be(BindingStrength.Strong); // 1 nM < 50 nM
    }

    // §7.1 / §2.2: a peptide matching NONE of the listed residues scores 0 (intercept 0 + three neutral 0s)
    // ⇒ IC50 = 50000^(1−0) = 50000 nM → NonBinder. Hand-derived, exact-ish equality + category.
    [Test]
    public void Smm_NoMatchingResidue_GivesMaxIc50_NonBinder()
    {
        var matrix = Matrix(0.0,
            Row(('A', 0.5)),
            Row(('C', 0.3)),
            Row(('D', 0.2)));

        // "WWW" matches no row ⇒ score 0 ⇒ IC50 = 50000.
        PredictIc50Smm("WWW", matrix).Should().BeApproximately(50000.0, 1e-6);

        var (ic50, strength) = PredictAndClassifySmm("WWW", matrix);
        ic50.Should().BeApproximately(50000.0, 1e-6);
        strength.Should().Be(BindingStrength.NonBinder); // 50000 ≥ 500 nM
    }

    // INV-05 spot checks: score 0 → 50000, score 1 → 1, score 0.5 → √50000. Intercept carries the whole
    // score here (empty-contribution rows), so this isolates the 50000^(1−score) transform.
    [Test]
    public void Smm_TransformSpotChecks_ScoreToIc50()
    {
        // score = intercept (rows contribute 0 for the unlisted 'A').
        PredictIc50Smm("A", Matrix(0.0, Row(('Q', 9.9)))).Should().BeApproximately(50000.0, 1e-6); // 50000^1
        PredictIc50Smm("A", Matrix(1.0, Row(('Q', 9.9)))).Should().BeApproximately(1.0, 1e-9);      // 50000^0
        PredictIc50Smm("A", Matrix(0.5, Row(('Q', 9.9)))).Should().BeApproximately(Math.Sqrt(50000.0), 1e-6);
    }

    #endregion

    #region MHC-MATRIX-001 — BE: peptide length ≠ matrix length ⇒ documented ArgumentException (NOT IndexOOR)

    // INV-07 / §3.3: a peptide SHORTER than the matrix must throw a documented ArgumentException, never an
    // IndexOutOfRangeException (the loop runs to matrix.Rows.Count, which would over-index the peptide if
    // the guard were missing — but here the matrix-side check fires first regardless).
    [Test]
    public void Bimas_PeptideShorterThanMatrix_ThrowsArgumentException()
    {
        var matrix = Matrix(1.0, Row(('A', 2.0)), Row(('C', 2.0)), Row(('D', 2.0))); // 3 rows
        Action act = () => PredictBindingHalfLifeBimas("AC", matrix);                 // length 2
        act.Should().Throw<ArgumentException>();
        act.Should().NotThrow<IndexOutOfRangeException>("a length mismatch is a documented ArgumentException");
    }

    // INV-07: a peptide LONGER than the matrix must throw ArgumentException, never IndexOutOfRange from
    // indexing matrix.Rows[i] beyond its count.
    [Test]
    public void Bimas_PeptideLongerThanMatrix_ThrowsArgumentException()
    {
        var matrix = Matrix(1.0, Row(('A', 2.0)), Row(('C', 2.0))); // 2 rows
        Action act = () => PredictBindingHalfLifeBimas("ACDEF", matrix); // length 5
        act.Should().Throw<ArgumentException>();
    }

    [Test]
    public void Smm_LengthMismatch_BothDirections_ThrowArgumentException()
    {
        var matrix = Matrix(0.0, Row(('A', 0.5)), Row(('C', 0.5))); // 2 rows
        Action shorter = () => PredictIc50Smm("A", matrix);
        Action longer = () => PredictIc50Smm("ACDE", matrix);
        shorter.Should().Throw<ArgumentException>();
        longer.Should().Throw<ArgumentException>();
    }

    [Test]
    public void PredictAndClassifySmm_LengthMismatch_ThrowsArgumentException()
    {
        var matrix = Matrix(0.0, Row(('A', 0.5)), Row(('C', 0.5)));
        Action act = () => PredictAndClassifySmm("AAACCC", matrix);
        act.Should().Throw<ArgumentException>();
    }

    // BE fuzz: random length mismatches (peptide length deliberately != row count) over many seeds must
    // ALWAYS surface a documented ArgumentException and NEVER an IndexOutOfRange / KeyNotFound / any other
    // unhandled type.
    [Test]
    [CancelAfter(20_000)]
    public void Predictors_RandomLengthMismatch_AlwaysArgumentException_NeverIndexOOR()
    {
        var rng = new Random(246_0001);
        const string aa = "ACDEFGHIKLMNPQRSTVWY";
        for (int i = 0; i < 5_000; i++)
        {
            int rowCount = rng.Next(1, 16);
            var rows = new IReadOnlyDictionary<char, double>[rowCount];
            for (int r = 0; r < rowCount; r++)
            {
                rows[r] = Row((aa[rng.Next(aa.Length)], rng.NextDouble() * 2.0));
            }
            var matrix = new PmhcScoringMatrix(rows, rng.NextDouble());

            int pepLen;
            do { pepLen = rng.Next(0, 20); } while (pepLen == rowCount); // guarantee a mismatch
            var sb = new char[pepLen];
            for (int p = 0; p < pepLen; p++) sb[p] = aa[rng.Next(aa.Length)];
            string peptide = new string(sb);

            Action bimas = () => PredictBindingHalfLifeBimas(peptide, matrix);
            Action smm = () => PredictIc50Smm(peptide, matrix);
            bimas.Should().Throw<ArgumentException>($"len {pepLen} vs rows {rowCount} (BIMAS)");
            smm.Should().Throw<ArgumentException>($"len {pepLen} vs rows {rowCount} (SMM)");
        }
    }

    #endregion

    #region MHC-MATRIX-001 — BE: empty matrix (0 rows) ⇒ documented ArgumentException

    // INV-07 / §3.3: a matrix with 0 position rows is rejected with ArgumentException — NOT silently scored
    // as an empty product (1·const) for BIMAS or empty sum (50000^(1−intercept)) for SMM. Even an
    // empty-peptide + empty-matrix (where lengths "match" at 0) must still throw the no-rows exception.
    [Test]
    public void Bimas_EmptyMatrix_ThrowsArgumentException()
    {
        var empty = new PmhcScoringMatrix(Array.Empty<IReadOnlyDictionary<char, double>>(), 7.0);
        Action withPeptide = () => PredictBindingHalfLifeBimas("AC", empty);
        Action emptyPeptide = () => PredictBindingHalfLifeBimas("", empty); // lengths "match" at 0
        withPeptide.Should().Throw<ArgumentException>();
        emptyPeptide.Should().Throw<ArgumentException>("an empty matrix has no position rows (INV-07)");
    }

    [Test]
    public void Smm_EmptyMatrix_ThrowsArgumentException()
    {
        var empty = new PmhcScoringMatrix(Array.Empty<IReadOnlyDictionary<char, double>>(), 0.5);
        Action withPeptide = () => PredictIc50Smm("ACD", empty);
        Action emptyPeptide = () => PredictIc50Smm("", empty);
        withPeptide.Should().Throw<ArgumentException>();
        emptyPeptide.Should().Throw<ArgumentException>();
    }

    // A null Rows reference (degenerate default-constructed record) is also a "no rows" case ⇒ ArgumentException,
    // not a NullReferenceException leaking from matrix.Rows.Count.
    [Test]
    public void Predictors_NullRowsMatrix_ThrowArgumentException_NotNullReference()
    {
        var nullRows = new PmhcScoringMatrix(null!, 1.0);
        Action bimas = () => PredictBindingHalfLifeBimas("A", nullRows);
        Action smm = () => PredictIc50Smm("A", nullRows);
        // ArgumentException (the documented "no rows" guard), NOT a leaked NullReferenceException.
        bimas.Should().Throw<ArgumentException>();
        smm.Should().Throw<ArgumentException>();
        bimas.Should().NotThrow<NullReferenceException>();
        smm.Should().NotThrow<NullReferenceException>();
    }

    #endregion

    #region MHC-MATRIX-001 — Null peptide ⇒ ArgumentNullException (INV-07 / §3.3)

    [Test]
    public void Predictors_NullPeptide_ThrowArgumentNullException()
    {
        var matrix = Matrix(1.0, Row(('A', 2.0)));
        Action bimas = () => PredictBindingHalfLifeBimas(null!, matrix);
        Action smm = () => PredictIc50Smm(null!, matrix);
        Action chain = () => PredictAndClassifySmm(null!, matrix);
        bimas.Should().Throw<ArgumentNullException>();
        smm.Should().Throw<ArgumentNullException>();
        chain.Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region MHC-MATRIX-001 — MC: missing residue & non-amino-acid characters (defined neutral, no crash)

    // MC: non-amino-acid characters (digits, punctuation, whitespace, unicode) are residues absent from every
    // row ⇒ neutral contribution. BIMAS: a peptide of all non-AA chars scores 1.0·const (every coeff 1.0);
    // SMM: scores the intercept (every contribution 0). No KeyNotFound. Hand-derived equality.
    [Test]
    public void Bimas_AllNonAaPeptide_ScoresFinalConstant_NoThrow()
    {
        var matrix = Matrix(7.0, Row(('A', 2.0)), Row(('C', 3.0)), Row(('D', 4.0)));
        // "1#%" matches no listed residue at any position ⇒ ∏ neutral 1.0 ⇒ T½ = 7.0.
        double t = 0;
        Action act = () => t = PredictBindingHalfLifeBimas("1#%", matrix);
        act.Should().NotThrow();
        t.Should().Be(7.0);
        t.Should().Be(ExpectedBimas("1#%", matrix));
    }

    [Test]
    public void Smm_AllNonAaPeptide_ScoresIntercept_NoThrow()
    {
        var matrix = Matrix(0.4, Row(('A', 0.5)), Row(('C', 0.5))); // intercept 0.4
        // " @" matches nothing ⇒ score = 0.4 ⇒ IC50 = 50000^(1−0.4) = 50000^0.6.
        double ic50 = 0;
        Action act = () => ic50 = PredictIc50Smm(" @", matrix);
        act.Should().NotThrow();
        ic50.Should().BeApproximately(Math.Pow(50000.0, 0.6), 1e-6);
        ic50.Should().BeApproximately(ExpectedSmmIc50(" @", matrix), 1e-6);
        double.IsFinite(ic50).Should().BeTrue("the SMM transform must never produce NaN/Inf");
        ic50.Should().BeGreaterThan(0.0);
    }

    // Lower-casing: the predictors upper-case the peptide before lookup, so a lower-case AA still matches an
    // upper-case row key. Re-derived: "lmv" must equal the "LMV" worked example (T½ = 90).
    [Test]
    public void Bimas_LowerCasePeptide_MatchesUpperCaseRows()
    {
        var matrix = Matrix(10.0, Row(('L', 2.0)), Row(('M', 3.0)), Row(('V', 1.5)));
        PredictBindingHalfLifeBimas("lmv", matrix).Should().Be(90.0);
    }

    // MC fuzz: random peptides over a polluted alphabet (AAs + digits + punctuation + unicode), each scored
    // against a random matrix of matching length, must NEVER throw and must EXACTLY match the independent
    // oracle (BIMAS product, SMM IC50). This pins the missing-residue/non-AA neutral handling broadly and
    // guards against any NaN/Inf SMM output.
    [Test]
    [CancelAfter(30_000)]
    public void Predictors_RandomPollutedPeptides_MatchOracle_NoThrow_FiniteIc50()
    {
        var rng = new Random(246_0002);
        const string polluted = "ACDEFGHIKLMNPQRSTVWYacdefg0123456789#@%* \tЖλ";
        const string aa = "ACDEFGHIKLMNPQRSTVWY";

        for (int i = 0; i < 10_000; i++)
        {
            int len = rng.Next(1, 16);
            var rows = new IReadOnlyDictionary<char, double>[len];
            for (int r = 0; r < len; r++)
            {
                // sparse rows: list only a couple of residues so "missing residue" is exercised often.
                int listed = rng.Next(0, 4);
                var entries = new (char, double)[listed];
                for (int e = 0; e < listed; e++)
                {
                    entries[e] = (aa[rng.Next(aa.Length)], rng.NextDouble() * 2.0 - 0.5);
                }
                rows[r] = Row(entries);
            }

            // BIMAS matrix: positive-ish coefficients; SMM matrix: small contributions + small intercept.
            var bimasMatrix = new PmhcScoringMatrix(rows, rng.NextDouble() * 5.0);
            var smmMatrix = new PmhcScoringMatrix(rows, rng.NextDouble() - 0.5);

            var sb = new char[len];
            for (int p = 0; p < len; p++) sb[p] = polluted[rng.Next(polluted.Length)];
            string peptide = new string(sb);

            double t = 0, ic50 = 0;
            Action bimas = () => t = PredictBindingHalfLifeBimas(peptide, bimasMatrix);
            Action smm = () => ic50 = PredictIc50Smm(peptide, smmMatrix);
            bimas.Should().NotThrow($"peptide '{peptide}' (BIMAS)");
            smm.Should().NotThrow($"peptide '{peptide}' (SMM)");

            t.Should().BeApproximately(ExpectedBimas(peptide, bimasMatrix), 1e-9);
            ic50.Should().BeApproximately(ExpectedSmmIc50(peptide, smmMatrix), Math.Abs(ExpectedSmmIc50(peptide, smmMatrix)) * 1e-9 + 1e-9);
            double.IsFinite(ic50).Should().BeTrue($"SMM IC50 must be finite for '{peptide}'");
            ic50.Should().BeGreaterThan(0.0, $"SMM IC50 = 50000^(1−score) is always > 0 for '{peptide}'");
        }
    }

    #endregion

    #region MHC-MATRIX-001 — Theory invariants: directionality & SMM strict monotonicity

    // INV-06 directionality: a HIGHER BIMAS product ⇒ stronger predicted binder. Two matrices identical
    // except one position's coefficient (2.0 vs 4.0) ⇒ the larger coefficient yields the larger T½.
    [Test]
    public void Bimas_LargerCoefficient_GivesLargerHalfLife()
    {
        var weak = Matrix(1.0, Row(('A', 2.0)), Row(('C', 1.0)));
        var strong = Matrix(1.0, Row(('A', 4.0)), Row(('C', 1.0)));
        double tWeak = PredictBindingHalfLifeBimas("AC", weak);
        double tStrong = PredictBindingHalfLifeBimas("AC", strong);
        tStrong.Should().BeGreaterThan(tWeak, "higher BIMAS coefficient ⇒ higher T½ ⇒ stronger binder");
        tWeak.Should().Be(2.0);   // 1·(2·1)
        tStrong.Should().Be(4.0); // 1·(4·1)
    }

    // INV-05 strict monotonicity: IC50 = 50000^(1−score) is strictly DECREASING in score. Over random
    // score pairs, the larger total score must give the smaller (stronger) IC50.
    [Test]
    [CancelAfter(20_000)]
    public void Smm_Ic50_StrictlyDecreasingInScore()
    {
        var rng = new Random(246_0003);
        for (int i = 0; i < 20_000; i++)
        {
            double sLow = rng.NextDouble() * 1.5 - 0.25;  // ~[-0.25, 1.25)
            double sHigh = sLow + rng.NextDouble() * 1.0 + 1e-3; // strictly greater

            // Single-row matrix with intercept = score, peptide unlisted ⇒ total score == intercept.
            double ic50Low = PredictIc50Smm("A", Matrix(sLow, Row(('Q', 0.0))));
            double ic50High = PredictIc50Smm("A", Matrix(sHigh, Row(('Q', 0.0))));

            double.IsFinite(ic50Low).Should().BeTrue();
            double.IsFinite(ic50High).Should().BeTrue();
            ic50High.Should().BeLessThan(ic50Low,
                $"larger SMM score {sHigh} ⇒ smaller (stronger) IC50 than score {sLow}");
        }
    }

    #endregion
}
