using System;
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Oncology;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Fuzz tests for the Oncology MHCflurry pan-allele binding-<b>affinity</b> neural-network predictor —
/// MHC-NN-001. The units under test are
/// <see cref="MhcflurryAffinityPredictor.EncodePeptide(string)"/> (BLOSUM62 left_pad_centered_right_pad
/// peptide encoding, supported length window),
/// <see cref="MhcflurryAffinityPredictor.EncodePseudosequence(string)"/> (37×21 allele encoding),
/// <see cref="MhcflurryAffinityPredictor.GetPseudosequence(string)"/> (bundled allele→pseudosequence table),
/// <see cref="MhcflurryAffinityPredictor.ToIc50(double)"/> (the <c>IC50 = 50000^(1−x)</c> transform),
/// <see cref="MhcflurryAffinityPredictor.Network.ForwardRaw(ReadOnlySpan{double})"/> /
/// <see cref="MhcflurryAffinityPredictor.PredictIc50(IReadOnlyList{MhcflurryAffinityPredictor.Network}, string, string)"/>
/// (the per-network forward pass + the geometric-mean ensemble IC50), all in
/// src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/MhcflurryAffinityPredictor.cs.
///
/// This file is scoped STRICTLY to MHC-NN-001 (the MHCflurry NEURAL affinity predictor). It does NOT touch the
/// specification-driven threshold/matrix classifier ONCO-MHC-001 (row 110,
/// <see cref="OncologyAnalyzer.ClassifyBindingAffinity(double)"/> / <see cref="OncologyAnalyzer.ClassifyMhcBinding(int, double, OncologyAnalyzer.MhcClass)"/>),
/// which is fuzzed in OncologyMhcBindingFuzzTests.cs.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds malformed / boundary / degenerate inputs and asserts the predictor never fails in an
/// undisciplined way: no IndexOutOfRange from the peptide padding/encoding, no NaN/∞ affinity leaking out of
/// the sigmoid→to_ic50 chain, no negative IC50, no IC50 outside the documented nM range, no silent wrong
/// prediction on an unknown allele, no hang. Every input must resolve to EITHER a well-defined, theory-correct
/// result OR a *documented, intentional* validation exception:
///   • peptide length &lt; 5 or &gt; 15 (outside the left_pad_centered_right_pad window) → ArgumentOutOfRangeException;
///   • null peptide / null allele / null pseudosequence → ArgumentNullException;
///   • unknown allele name → KeyNotFoundException (NOT a silent wrong prediction);
///   • non-amino-acid residues (B/J/O/U/X, digits, unicode) → the BLOSUM62 X (unknown) column, never a throw;
///   • empty pseudosequence → a wrong-width network input → ArgumentException at the forward pass (not a crash);
///   • empty ensemble → ArgumentException.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test (independently derived from the source header + docs)
/// ───────────────────────────────────────────────────────────────────────────
/// MhcflurryAffinityPredictor.cs class/method XML docs + docs/algorithms/Oncology/MHC_Peptide_Binding_Classification.md §5.3:
///   • Supported peptide length window: PeptideMinLength=5 .. PeptideMaxLength=15 inclusive
///     (encodable_sequences.py "minimum length of 5"; models_class1_pan max_length=15). Outside ⇒
///     ArgumentOutOfRangeException (EncodePeptide XML doc / source lines 254–259).
///   • left_pad_centered_right_pad: 3 copies (left, centred, right) at 15 positions, each a 21-wide BLOSUM62
///     vector ⇒ PeptideFlatLength=945. Padding uses the X vector. (source lines 250–308.)
///   • Allele encoding: residues×21, position-major; a 37-residue pseudosequence ⇒ AlleleFlatLength=777.
///   • Non-canonical residues fall back to the X (unknown) column (allow_unsupported_amino_acids); lowercase is
///     folded case-insensitively to its canonical residue (BuildAminoAcidIndex, IndexOfResidue, lines 140–160).
///   • ToIc50(x) = 50000^(1−x); strictly decreasing; for x∈[0,1] the IC50 ∈ [1, 50000] nM (MaxIc50Nm=50000).
///   • Network output is a sigmoid ∈ (0,1) ⇒ a VALID (peptide, allele) gives a FINITE, POSITIVE IC50 in the
///     documented (1, 50000) nM range. The forward pass is DETERMINISTIC.
///   • A known strong binder (GILGFVFTL / HLA-A*02:01, flu M1) scores a much LOWER nM than a non-binder
///     (SIINFEKL / HLA-A*02:01, ovalbumin) for the same allele — the anchor in §M5 / MHCflurry oracle.
///   • GetPseudosequence: unknown allele ⇒ KeyNotFoundException; null ⇒ ArgumentNullException (lines 184–194).
///   • Empty ensemble ⇒ ArgumentException (lines 642–645).
///
/// The MHCflurry forward pass needs trained weights; this fixture loads the SAME embedded single-network weight
/// pack (the smallest models_class1_pan member) used by MhcflurryAffinityPredictor_PredictIc50_Tests.cs, so the
/// fuzz inputs are scored by the REAL trained network, not a stub. The predictor is composed of pure static
/// methods with no LimitationPolicy guard, so no Permissive bootstrap is required.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyMhcNnFuzzTests
{
    // The embedded single-network weight pack = the smallest models_class1_pan member (feedforward [512,512]),
    // shared with MhcflurryAffinityPredictor_PredictIc50_Tests. Loading it runs the full real forward pass.
    private const string SingleNetResource =
        "Seqeron.Genomics.Tests.TestData.Mhcflurry.mhcflurry_single_net.bin";

    // HLA-A*02:01 pseudosequence (37 residues) from the bundled MHCflurry table — the workhorse valid allele.
    private const string HlaA0201Pseudo = "YFAMYGEKVAHTHVDTLYGVRYDHYYTWAVLAYTWYA";

    // Documented IC50 bounds: output ∈ [0,1] ⇒ IC50 = 50000^(1−output) ∈ [1, 50000]. The sigmoid is open on
    // (0,1) so in practice the IC50 is strictly inside, but we assert the closed documented envelope to be safe.
    private const double MinIc50Nm = 1.0;
    private const double MaxIc50Nm = 50000.0;

    private static IReadOnlyList<MhcflurryAffinityPredictor.Network> LoadSingleNet()
    {
        var asm = typeof(OncologyMhcNnFuzzTests).Assembly;
        using Stream stream = asm.GetManifestResourceStream(SingleNetResource)
            ?? throw new InvalidOperationException($"Embedded resource '{SingleNetResource}' not found.");
        return MhcflurryAffinityPredictor.LoadWeightPack(stream);
    }

    // A predicted IC50 is in-contract iff it is finite, strictly positive, and inside the documented nM range.
    private static void AssertWellFormedIc50(double ic50, string because)
    {
        double.IsNaN(ic50).Should().BeFalse($"a predicted IC50 must never be NaN — {because}");
        double.IsInfinity(ic50).Should().BeFalse($"a predicted IC50 must never be ±∞ — {because}");
        ic50.Should().BeGreaterThan(0.0, $"a predicted IC50 is a positive concentration — {because}");
        ic50.Should().BeInRange(MinIc50Nm, MaxIc50Nm,
            $"a sigmoid-output IC50 must lie in the documented [{MinIc50Nm}, {MaxIc50Nm}] nM range — {because}");
    }

    // A small corpus of random valid (length-in-window, canonical-residue) peptides.
    private static string RandomCanonicalPeptide(Random rng, int len)
    {
        const string order = MhcflurryAffinityPredictor.AminoAcidOrder; // ACDEFGHIKLMNPQRSTVWYX
        var chars = new char[len];
        for (int i = 0; i < len; i++)
        {
            // pick from the 20 canonical residues (exclude the trailing X) so the input is a real peptide
            chars[i] = order[rng.Next(20)];
        }
        return new string(chars);
    }

    #region MHC-NN-001 — Positive sanity: valid (peptide, allele) ⇒ finite positive nM in range

    // A valid 9-mer on a real allele yields a finite, positive IC50 inside the documented [1, 50000] nM range,
    // and the forward pass is deterministic (the same input gives bit-identical output on repeat calls).
    [Test]
    public void PredictIc50_ValidPeptideAllele_FiniteInRangeAndDeterministic()
    {
        var nets = LoadSingleNet();

        double a = MhcflurryAffinityPredictor.PredictIc50(nets, "GILGFVFTL", "HLA-A*02:01");
        double b = MhcflurryAffinityPredictor.PredictIc50(nets, "GILGFVFTL", "HLA-A*02:01");

        AssertWellFormedIc50(a, "valid 9-mer / HLA-A*02:01");
        a.Should().Be(b, "the forward pass is a pure deterministic function of (peptide, allele)");
    }

    // Theory anchor (§M5 / MHCflurry oracle): a known HLA-A*02:01 strong binder (GILGFVFTL, flu M1) must score a
    // much LOWER nM than a non-binder (SIINFEKL, ovalbumin) for the SAME allele — not a silent equal/inverted
    // prediction. This pins that the network is doing real allele-aware scoring, not returning a constant.
    [Test]
    public void PredictIc50_KnownStrongBinder_ScoresLowerThanNonBinder_SameAllele()
    {
        var nets = LoadSingleNet();

        double strong = MhcflurryAffinityPredictor.PredictIc50(nets, "GILGFVFTL", "HLA-A*02:01");
        double nonBinder = MhcflurryAffinityPredictor.PredictIc50(nets, "SIINFEKL", "HLA-A*02:01");

        AssertWellFormedIc50(strong, "strong binder");
        AssertWellFormedIc50(nonBinder, "non-binder");
        strong.Should().BeLessThan(nonBinder,
            "the flu-M1 strong binder must have a lower (stronger) IC50 than the ovalbumin non-binder on HLA-A*02:01");
        (nonBinder / strong).Should().BeGreaterThan(10.0,
            "the affinities must differ by orders of magnitude, not be a coin flip");
    }

    // Fuzz: across all supported lengths (5..15) and many random canonical peptides on a real allele, EVERY
    // prediction is finite, positive and inside the documented nM range — no IndexOutOfRange from the
    // left_pad_centered_right_pad placement at any length, no NaN/∞ from the forward pass.
    [Test]
    [CancelAfter(60_000)]
    public void PredictIc50_RandomCanonicalPeptides_AllLengths_FiniteInRange()
    {
        var nets = LoadSingleNet();
        var rng = new Random(245_0001);

        for (int len = MhcflurryAffinityPredictor.PeptideMinLength; len <= MhcflurryAffinityPredictor.PeptideMaxLength; len++)
        {
            for (int i = 0; i < 200; i++)
            {
                string pep = RandomCanonicalPeptide(rng, len);
                double ic50 = MhcflurryAffinityPredictor.PredictIc50(nets, pep, "HLA-A*02:01");
                AssertWellFormedIc50(ic50, $"random {len}-mer '{pep}'");
            }
        }
    }

    #endregion

    #region MHC-NN-001 — BE: peptide length < 5 or > 15 (the encoding window boundary)

    // §EncodePeptide: the supported length window is [5, 15]. The inclusive endpoints encode to exactly the
    // 945-wide flat vector (no over/under-flow); just outside (4 and 16) throws the documented exception.
    [Test]
    public void EncodePeptide_WindowEndpoints_EncodeAndJustOutsideThrows()
    {
        MhcflurryAffinityPredictor.EncodePeptide(new string('A', 5)).Length
            .Should().Be(MhcflurryAffinityPredictor.PeptideFlatLength, "length 5 is the supported minimum");
        MhcflurryAffinityPredictor.EncodePeptide(new string('A', 15)).Length
            .Should().Be(MhcflurryAffinityPredictor.PeptideFlatLength, "length 15 is the supported maximum");

        Action below = () => MhcflurryAffinityPredictor.EncodePeptide(new string('A', 4));
        Action above = () => MhcflurryAffinityPredictor.EncodePeptide(new string('A', 16));
        below.Should().Throw<ArgumentOutOfRangeException>("length 4 is below the minimum");
        above.Should().Throw<ArgumentOutOfRangeException>("length 16 is above the maximum");
    }

    // Fuzz: every length outside [5, 15] — including 0 (empty), 1, the just-outside 4 and 16, and far-out 50,
    // 200 — throws ArgumentOutOfRangeException. Never an IndexOutOfRange or a silently-truncated encoding.
    [Test]
    public void EncodePeptide_OutOfWindowLengths_ThrowArgumentOutOfRange()
    {
        foreach (int len in new[] { 0, 1, 2, 3, 4, 16, 17, 20, 31, 50, 200 })
        {
            Action act = () => MhcflurryAffinityPredictor.EncodePeptide(new string('A', len));
            act.Should().Throw<ArgumentOutOfRangeException>($"length {len} is outside [5, 15]");
        }
    }

    // The same boundary holds end-to-end: PredictIc50 with an out-of-window peptide surfaces the encoder's
    // ArgumentOutOfRangeException — it must NOT reach the network with a wrong-width input and crash there.
    [Test]
    public void PredictIc50_OutOfWindowPeptide_ThrowsArgumentOutOfRange()
    {
        var nets = LoadSingleNet();
        foreach (string pep in new[] { "", "A", "AAAA", new string('A', 16), new string('A', 40) })
        {
            Action act = () => MhcflurryAffinityPredictor.PredictIc50(nets, pep, "HLA-A*02:01");
            act.Should().Throw<ArgumentOutOfRangeException>($"peptide of length {pep.Length} is out of the [5,15] window");
        }
    }

    // Null peptide is the documented ArgumentNullException (distinct from the length guard), at both the encoder
    // and the ensemble entry point.
    [Test]
    public void EncodePeptide_And_PredictIc50_NullPeptide_ThrowArgumentNull()
    {
        var nets = LoadSingleNet();
        ((Action)(() => MhcflurryAffinityPredictor.EncodePeptide(null!)))
            .Should().Throw<ArgumentNullException>();
        ((Action)(() => MhcflurryAffinityPredictor.PredictIc50(nets, null!, "HLA-A*02:01")))
            .Should().Throw<ArgumentNullException>();
    }

    #endregion

    #region MHC-NN-001 — MC: non-amino-acid residues fall back to the BLOSUM62 X column (never throw)

    // §IndexOfResidue / allow_unsupported_amino_acids: B, J, O, U, X, digits, unicode all encode as the X
    // (unknown) BLOSUM62 column — the X row is all-zero except its own channel (=1). No exception, no IndexOOR.
    [Test]
    public void EncodePeptide_NonAminoAcidResidues_EncodeAsXVector()
    {
        int xIdx = MhcflurryAffinityPredictor.AminoAcidOrder.IndexOf(MhcflurryAffinityPredictor.UnknownAminoAcid);

        // A 5-mer (min length) of pure non-AA junk: every position of block 0 must be the X vector.
        foreach (string junk in new[] { "BJOUX", "12345", "B*Z?#" })
        {
            double[] flat = MhcflurryAffinityPredictor.EncodePeptide(junk);
            flat.Length.Should().Be(MhcflurryAffinityPredictor.PeptideFlatLength);
            for (int p = 0; p < junk.Length; p++)
            {
                for (int c = 0; c < MhcflurryAffinityPredictor.EncodingWidth; c++)
                {
                    double expected = c == xIdx ? 1.0 : 0.0;
                    flat[p * MhcflurryAffinityPredictor.EncodingWidth + c].Should().Be(expected,
                        $"non-AA residue '{junk[p]}' at position {p} channel {c} must encode as the X column");
                }
            }
        }
    }

    // Unicode / control / non-ASCII (≥128) residues take the explicit `residue < 128 ? ... : X` branch — they
    // must encode as X, never index out of the 128-wide ASCII map.
    [Test]
    public void EncodePeptide_UnicodeAndNonAsciiResidues_FallBackToXNoThrow()
    {
        int xIdx = MhcflurryAffinityPredictor.AminoAcidOrder.IndexOf(MhcflurryAffinityPredictor.UnknownAminoAcid);
        // π, é, a CJK char, and a high surrogate-free BMP char — all ≥128.
        string pep = "π" + "é" + "中" + "Ω" + "Ж"; // 5 non-ASCII code units, length 5
        pep.Length.Should().Be(5, "each chosen char is a single UTF-16 code unit");

        double[] flat = MhcflurryAffinityPredictor.EncodePeptide(pep);
        for (int p = 0; p < pep.Length; p++)
        {
            flat[p * MhcflurryAffinityPredictor.EncodingWidth + xIdx].Should().Be(1.0,
                $"non-ASCII residue at position {p} must encode as the X column");
        }
    }

    // Case-insensitivity: a lowercase peptide encodes IDENTICALLY to its uppercase form (folded to canonical),
    // and a mixed-junk peptide still scores to a finite in-range IC50 through the full forward pass.
    [Test]
    public void EncodePeptide_Lowercase_FoldsToUppercase_And_JunkPeptideScoresInRange()
    {
        MhcflurryAffinityPredictor.EncodePeptide("gilgfvftl")
            .Should().Equal(MhcflurryAffinityPredictor.EncodePeptide("GILGFVFTL"),
                "lowercase residues fold case-insensitively to the canonical encoding");

        var nets = LoadSingleNet();
        // A peptide mixing canonical, lowercase, and unknown residues still produces a well-formed IC50.
        double ic50 = MhcflurryAffinityPredictor.PredictIc50(nets, "GxLBfVuTl", "HLA-A*02:01");
        AssertWellFormedIc50(ic50, "mixed canonical/lowercase/unknown peptide");
    }

    // Fuzz: random peptides over an ARBITRARY printable+junk alphabet (length in window) never throw and always
    // score to a finite, positive, in-range IC50 — the X fallback fully absorbs malformed content.
    [Test]
    [CancelAfter(60_000)]
    public void PredictIc50_RandomMalformedAlphabet_AlwaysFiniteInRange()
    {
        var nets = LoadSingleNet();
        var rng = new Random(245_0002);
        const string alphabet = "ACDEFGHIKLMNPQRSTVWYacdefBJOUXZ*?#0123456789 .-";

        for (int i = 0; i < 4_000; i++)
        {
            int len = rng.Next(MhcflurryAffinityPredictor.PeptideMinLength,
                               MhcflurryAffinityPredictor.PeptideMaxLength + 1);
            var chars = new char[len];
            for (int j = 0; j < len; j++)
            {
                chars[j] = alphabet[rng.Next(alphabet.Length)];
            }
            string pep = new(chars);

            double ic50 = MhcflurryAffinityPredictor.PredictIc50(nets, pep, "HLA-A*02:01");
            AssertWellFormedIc50(ic50, $"malformed peptide '{pep}'");
        }
    }

    #endregion

    #region MHC-NN-001 — MC/BE: unknown allele / null allele / empty pseudosequence

    // §GetPseudosequence: an allele absent from the bundled table throws KeyNotFoundException — a documented
    // reject, NOT a silent wrong prediction. Tested both directly and end-to-end through PredictIc50.
    [Test]
    public void UnknownAllele_ThrowsKeyNotFound_NotSilentPrediction()
    {
        var nets = LoadSingleNet();
        foreach (string bad in new[] { "HLA-Z*99:99", "NOPE", "HLA-A*02:01x", "", " " })
        {
            ((Action)(() => MhcflurryAffinityPredictor.GetPseudosequence(bad)))
                .Should().Throw<KeyNotFoundException>($"allele '{bad}' is not in the bundled table");
            ((Action)(() => MhcflurryAffinityPredictor.PredictIc50(nets, "GILGFVFTL", bad)))
                .Should().Throw<KeyNotFoundException>($"PredictIc50 must surface the unknown-allele reject for '{bad}'");
        }
    }

    // Null allele is the documented ArgumentNullException (distinct from the unknown-allele KeyNotFound).
    [Test]
    public void NullAllele_ThrowsArgumentNull()
    {
        var nets = LoadSingleNet();
        ((Action)(() => MhcflurryAffinityPredictor.GetPseudosequence(null!)))
            .Should().Throw<ArgumentNullException>();
        ((Action)(() => MhcflurryAffinityPredictor.PredictIc50(nets, "GILGFVFTL", null!)))
            .Should().Throw<ArgumentNullException>();
    }

    // §EncodePseudosequence: an empty pseudosequence is itself encodable (→ a zero-length vector, no throw), but
    // feeding it through the network yields a wrong-width (945 ≠ 1722) input ⇒ the documented ArgumentException
    // from the forward pass — never an IndexOutOfRange or a NaN affinity from a malformed allele.
    [Test]
    public void EmptyPseudosequence_EncodesEmpty_ButForwardPassThrowsWrongWidth()
    {
        MhcflurryAffinityPredictor.EncodePseudosequence("").Length
            .Should().Be(0, "an empty pseudosequence encodes to a zero-length vector");

        var nets = LoadSingleNet();
        Action act = () => MhcflurryAffinityPredictor.PredictIc50WithPseudosequence(nets, "GILGFVFTL", "");
        act.Should().Throw<ArgumentException>(
            "a 945-wide (peptide-only) input is the wrong width for the 1722-wide network");

        // Null pseudosequence is the documented ArgumentNullException.
        ((Action)(() => MhcflurryAffinityPredictor.EncodePseudosequence(null!)))
            .Should().Throw<ArgumentNullException>();
    }

    // A non-37-residue (but non-empty) pseudosequence of junk residues also yields a wrong-width input and the
    // same documented ArgumentException — the malformed allele never silently produces a bogus IC50.
    [Test]
    public void WrongLengthPseudosequence_ForwardPassThrowsWrongWidth()
    {
        var nets = LoadSingleNet();
        foreach (string ps in new[] { "AAAA", new string('X', 36), new string('Z', 38), "12345" })
        {
            Action act = () => MhcflurryAffinityPredictor.PredictIc50WithPseudosequence(nets, "GILGFVFTL", ps);
            act.Should().Throw<ArgumentException>($"a {ps.Length}-residue pseudosequence is the wrong network width");
        }
    }

    // A full-length (37) pseudosequence of UNKNOWN (X) residues is well-formed width-wise: it must score to a
    // finite, in-range IC50 (the all-X allele is a defined, non-crashing input), not a NaN.
    [Test]
    public void AllUnknownPseudosequence_CorrectLength_ScoresFiniteInRange()
    {
        var nets = LoadSingleNet();
        string allX = new('X', MhcflurryAffinityPredictor.PseudosequenceLength);
        double ic50 = MhcflurryAffinityPredictor.PredictIc50WithPseudosequence(nets, "GILGFVFTL", allX);
        AssertWellFormedIc50(ic50, "all-X 37-residue pseudosequence");
    }

    #endregion

    #region MHC-NN-001 — BE: ToIc50 transform & empty ensemble

    // §ToIc50: for ANY raw output in [0,1] the IC50 = 50000^(1−x) is finite, positive and in [1, 50000], and is
    // strictly decreasing. Fuzz the whole [0,1] domain plus the exact endpoints — no NaN, no negative, no leak.
    [Test]
    public void ToIc50_OverUnitInterval_FiniteDecreasingInRange()
    {
        MhcflurryAffinityPredictor.ToIc50(0.0).Should().BeApproximately(50000.0, 1e-3, "x=0 ⇒ the 50000 nM ceiling");
        MhcflurryAffinityPredictor.ToIc50(1.0).Should().BeApproximately(1.0, 1e-9, "x=1 ⇒ 1 nM floor");

        var rng = new Random(245_0003);
        double prev = double.PositiveInfinity;
        for (int i = 0; i <= 1000; i++)
        {
            double x = i / 1000.0; // monotone sweep 0..1
            double ic50 = MhcflurryAffinityPredictor.ToIc50(x);
            AssertWellFormedIc50(ic50, $"ToIc50({x})");
            ic50.Should().BeLessThanOrEqualTo(prev, "ToIc50 is strictly decreasing in the network output");
            prev = ic50;

            double rx = rng.NextDouble();
            AssertWellFormedIc50(MhcflurryAffinityPredictor.ToIc50(rx), $"ToIc50(random {rx})");
        }
    }

    // §PredictIc50WithPseudosequence: an empty network ensemble is the documented ArgumentException — never a
    // division-by-zero / NaN geometric mean (exp(0/0)).
    [Test]
    public void EmptyEnsemble_ThrowsArgumentException()
    {
        var empty = new List<MhcflurryAffinityPredictor.Network>();
        ((Action)(() => MhcflurryAffinityPredictor.PredictIc50(empty, "GILGFVFTL", "HLA-A*02:01")))
            .Should().Throw<ArgumentException>("an empty ensemble cannot produce a prediction");
        ((Action)(() => MhcflurryAffinityPredictor.PredictIc50WithPseudosequence(empty, "GILGFVFTL", HlaA0201Pseudo)))
            .Should().Throw<ArgumentException>();
    }

    // §ForwardRaw: a wrong-width raw input is rejected with ArgumentException (no out-of-bounds matmul). Sweep a
    // range of bad widths including 0 and the off-by-one neighbours of the true 1722.
    [Test]
    public void ForwardRaw_WrongWidthInputs_ThrowArgumentException()
    {
        var nets = LoadSingleNet();
        foreach (int width in new[] { 0, 1, 944, 945, 1721, 1723, 5000 })
        {
            Action act = () => nets[0].ForwardRaw(new double[width]);
            act.Should().Throw<ArgumentException>($"input width {width} ≠ {MhcflurryAffinityPredictor.NetworkInputLength}");
        }
    }

    #endregion
}
