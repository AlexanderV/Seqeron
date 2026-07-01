using System;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Annotation;
using static Seqeron.Genomics.Annotation.MiRnaAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for MIRNA-CLEAVAGE-001 — Drosha/Dicer cleavage-site prediction on a pri-/pre-miRNA
/// hairpin (<see cref="MiRnaAnalyzer.PredictDroshaDicerCleavage"/>): "given the basal ssRNA-dsRNA
/// junction, where do Drosha and Dicer cut, and what are the resulting 5p mature / 3p (miRNA*)
/// spans?".
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies (docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing")
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to a unit and asserts the code
/// NEVER fails in an undisciplined way: no hang, no unhandled runtime exception — and for a routine
/// that EXTRACTS FIXED-WIDTH WINDOWS by running a measuring "ruler" past a variable-length string,
/// the headline hazard is an IndexOutOfRange / negative-substring-length when the +11 bp Drosha
/// ruler or the +22 nt Dicer ruler runs off the end of a TOO-SHORT precursor, or when the basal
/// junction is itself out of range. No out-of-contract output either (a cleavage coordinate outside
/// [0,len), a negative mature length). Every input must resolve to EITHER a well-defined,
/// theory-correct result (including the null = "cannot predict" outcome) OR a documented validation
/// exception. — docs/ADVANCED_TESTING_CHECKLIST.md §8.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: MIRNA-CLEAVAGE-001 — Drosha/Dicer cleavage-site prediction
/// Checklist: docs/checklists/03_FUZZING.md, row 255.
/// Strategies (per §Description): BE = Boundary Exploitation, MC = Malformed Content.
/// Fuzz targets for row 255: "precursor too short, no stem, non-ACGU, empty seq".
/// Source doc: docs/algorithms/MiRNA/Pre_miRNA_Detection.md (§2.2.2 the cleavage measuring rules +
///   §7.1 the hsa-miR-21-5p worked example).
/// Source: src/.../Seqeron.Genomics.Annotation/MiRnaAnalyzer.cs
///   • DroshaDicerCleavage? PredictDroshaDicerCleavage(string sequence, int basalJunction) (~line 2250)
///   • readonly record struct DroshaDicerCleavage(BasalJunction, DroshaCut5Prime, DroshaCut3Prime,
///       MatureStart, MatureEnd, StarStart, StarEnd, MatureSequence, StarSequence, Sequence,
///       ThreePrimeOverhang, HasCnncMotif)                                              (~line 2196)
///
/// ───────────────────────────────────────────────────────────────────────────
/// The cleavage contract under test (independently RE-DERIVED from the doc, NOT read off the code)
/// ───────────────────────────────────────────────────────────────────────────
/// PredictDroshaDicerCleavage applies the PUBLISHED measuring ("ruler") rules verbatim — it does
/// NOT fold the hairpin and does NOT use a trained classifier (Pre_miRNA_Detection.md §2.2.2):
///   • Drosha basal-junction ruler (Han 2006 [9]): Drosha cleaves ~11 bp from the basal ssRNA-dsRNA
///     junction; the 5' cut is the 5' end of the 5p mature:
///         DroshaCut5' = basalJunction + 11.
///   • Dicer 5'-counting ruler (Park 2011 [10]): Dicer cleaves ~22 nt from the Drosha 5' end, fixing
///     the mature length at 22 nt:  mature = [DroshaCut5', DroshaCut5' + 21]  (22 nt inclusive).
///   • RNase III 2-nt 3' overhang (Lee 2003 [11]): each cut leaves a 2-nt 3' overhang; the 3p span's
///     Drosha-generated 3' end sits 2 nt 3' of the 5p mature 3' end (LINEAR-coordinate approximation),
///     with the same ~22-nt length:  StarEnd = MatureEnd + 2,  StarStart = StarEnd − 21.
///   • CNNC motif (Auyeung 2013 [12]): a C-N-N-C 16–18 nt 3' of the Drosha cut sets HasCnncMotif —
///     reported, NOT required.
///
/// Re-derived offsets (the hand-derived pins, computed here, NOT read from the implementation):
///   DroshaCut5' = j + 11
///   MatureStart = DroshaCut5',  MatureEnd = DroshaCut5' + 21  (length MatureEnd − MatureStart + 1 = 22)
///   StarEnd     = MatureEnd + 2 = j + 11 + 21 + 2 = j + 34
///   StarStart   = StarEnd − 21 = j + 13
/// So a non-null result requires StarEnd = j + 34 < len, i.e. len ≥ j + 35. Any shorter precursor —
/// or any j with j < 0 or j ≥ len, or whose +11/+22 ruler overruns — returns null, NEVER an
/// IndexOutOfRange / negative-substring.
///
/// Pinned invariants (theory, re-derived here):
///   INV-CV-1 (NULLABLE BOUNDARY — too short / out of range → null, never IndexOutOfRange):
///            null/empty → null; basalJunction &lt; 0 or ≥ len → null; a precursor too short for the
///            full +11 / +22 / +2 ruler (len &lt; j + 35) → null. NEVER a runtime exception, never a
///            negative Substring length.
///   INV-CV-2 (IN-CONTRACT COORDINATES, all non-null results): every reported coordinate ∈ [0, len);
///            MatureStart ≤ MatureEnd, StarStart ≤ StarEnd; the extracted MatureSequence /
///            StarSequence equal the indexed sub-spans of the normalised Sequence; ThreePrimeOverhang
///            is exactly 2.
///   INV-CV-3 (HAND-DERIVED RULER — the discriminative pin, NOT a no-op): for any valid (sequence, j)
///            the returned coordinates EQUAL the independently re-derived ruler offsets above
///            (DroshaCut5' = j + 11, mature length exactly 22, StarEnd = MatureEnd + 2, etc.). A test
///            that read the offsets off the code would be invalid; these are computed from the doc.
///   INV-CV-4 (miRBase WORKED EXAMPLE — Pre_miRNA_Detection.md §7.1): an 11-nt lower stem prefixed to
///            the miR-21 stem region with j = 0 places the +11 Drosha cut at the annotated 5p start
///            and reproduces miRBase hsa-miR-21-5p (MIMAT0000076, "UAGCUUAUCAGACUGAUGUUGA") exactly.
///   INV-CV-5 (NON-ACGU / no-stem content does NOT crash the ruler): the ruler is purely positional —
///            it does not test pairing — so DNA 'T' (→ 'U'), 'N', digits, punctuation and homopolymer
///            (no-stem) content are normalised/measured without a throw; the result is a well-formed
///            record or null, never a malformed coordinate.
///
/// If a source-derived assertion and the code disagree, the CODE is wrong (fixed minimally per the
/// doc): an IndexOutOfRange / negative Substring on a too-short precursor or an out-of-range junction,
/// or a coordinate outside [0,len), would each be a REAL bug to fix in the source — never softened
/// in the test.
///
/// LimitationPolicy: PredictDroshaDicerCleavage calls LimitationPolicy.Enforce("MIRNA-CLEAVAGE-001")
/// (the 3p/star span is a linear 2-nt-overhang approximation; Strict throws). The assembly
/// module-initializer (_LimitationPolicyTestBootstrap) already sets LimitationMode.Permissive
/// assembly-wide, so these calls return the approximate record rather than throwing.
///
/// All inputs are fixed / deterministically generated; the random helper uses a LOCALLY seeded
/// `new Random(seed)` (no shared static Rng), so every fuzz input is reproducible.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MiRnaCleavageFuzzTests
{
    #region Helpers

    /// <summary>Deterministic RNG — seed fixed locally so generated fuzz inputs are reproducible.</summary>
    private static string RandomRna(int length, int seed)
    {
        const string bases = "ACGU";
        var rng = new Random(seed);
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    // The published ruler constants, re-stated here from the doc (Han 2006 / Park 2011 / Lee 2003) —
    // NOT imported from the implementation, so the pins are an independent derivation.
    private const int DroshaBpFromJunction = 11; // Han 2006: ~11 bp from the basal junction
    private const int MatureLength = 22;          // Park 2011: Dicer 5'-counting fixes mature at ~22 nt
    private const int Overhang = 2;               // Lee 2003: RNase III 2-nt 3' overhang

    /// <summary>
    /// Asserts the universal cleavage-record invariants (INV-CV-2/3) by RE-DERIVING every coordinate
    /// from the documented ruler — every coordinate is in [0,len), spans are ordered, the mature is
    /// exactly 22 nt, the star end is the mature end + 2-nt overhang, and the extracted sub-sequences
    /// equal the indexed spans of the normalised sequence.
    /// </summary>
    private static void AssertCleavageMatchesRuler(DroshaDicerCleavage c, int basalJunction)
    {
        int len = c.Sequence.Length;

        // Hand-derived ruler offsets (the discriminative pin — computed from the doc, not the code).
        int expectedDrosha = basalJunction + DroshaBpFromJunction;          // j + 11
        int expectedMatureStart = expectedDrosha;
        int expectedMatureEnd = expectedMatureStart + MatureLength - 1;     // start + 21
        int expectedStarEnd = expectedMatureEnd + Overhang;                 // matureEnd + 2
        int expectedStarStart = expectedStarEnd - MatureLength + 1;         // starEnd - 21

        c.BasalJunction.Should().Be(basalJunction, "the record echoes the supplied basal junction (INV-CV-3)");
        c.DroshaCut5Prime.Should().Be(expectedDrosha,
            "Drosha cuts 11 bp from the basal junction: DroshaCut5' = j + 11 (Han 2006; INV-CV-3)");
        c.MatureStart.Should().Be(expectedMatureStart, "the 5p mature begins at the Drosha 5' cut (INV-CV-3)");
        c.MatureEnd.Should().Be(expectedMatureEnd,
            "Dicer 5'-counting fixes the mature at 22 nt: MatureEnd = DroshaCut5' + 21 (Park 2011; INV-CV-3)");
        c.StarEnd.Should().Be(expectedStarEnd,
            "the 3p Drosha 3' end sits 2 nt 3' of the mature 3' end: StarEnd = MatureEnd + 2 (Lee 2003; INV-CV-3)");
        c.StarStart.Should().Be(expectedStarStart,
            "the 3p span has the same 22-nt length: StarStart = StarEnd − 21 (INV-CV-3)");
        c.DroshaCut3Prime.Should().Be(expectedStarEnd, "DroshaCut3' is the 3p span's Drosha-generated 3' end (INV-CV-3)");

        // Mature length is EXACTLY 22 nt (the Dicer 5'-counting rule).
        (c.MatureEnd - c.MatureStart + 1).Should().Be(MatureLength,
            "the Dicer 5'-counting rule fixes the mature length at 22 nt (INV-CV-2)");
        (c.StarEnd - c.StarStart + 1).Should().Be(MatureLength, "the 3p span is also 22 nt (INV-CV-2)");
        c.ThreePrimeOverhang.Should().Be(Overhang, "every RNase III cut leaves a 2-nt 3' overhang (INV-CV-2)");

        // Every coordinate is a valid 0-based index inside the normalised sequence (no overrun).
        foreach (int coord in new[]
                 {
                     c.BasalJunction, c.DroshaCut5Prime, c.DroshaCut3Prime,
                     c.MatureStart, c.MatureEnd, c.StarStart, c.StarEnd,
                 })
            coord.Should().BeInRange(0, len - 1, "every cleavage coordinate is a 0-based index in [0,len) (INV-CV-2)");

        c.MatureStart.Should().BeLessThanOrEqualTo(c.MatureEnd, "the mature span is ordered (INV-CV-2)");
        c.StarStart.Should().BeLessThanOrEqualTo(c.StarEnd, "the star span is ordered (INV-CV-2)");

        // The extracted sub-sequences equal the indexed spans of the normalised sequence (no off-by-one).
        c.MatureSequence.Should().Be(
            c.Sequence.Substring(c.MatureStart, c.MatureEnd - c.MatureStart + 1),
            "MatureSequence equals the indexed mature span of the normalised sequence (INV-CV-2)");
        c.StarSequence.Should().Be(
            c.Sequence.Substring(c.StarStart, c.StarEnd - c.StarStart + 1),
            "StarSequence equals the indexed star span of the normalised sequence (INV-CV-2)");
        c.MatureSequence.Length.Should().Be(MatureLength, "the extracted mature is 22 nt (INV-CV-2)");
        c.Sequence.Should().NotContain("T", "the indexed Sequence is normalised RNA (T→U) (INV-CV-2)");
    }

    // hsa-miR-21-5p (MIMAT0000076) — the doc's §7.1 worked-example mature miRNA.
    private const string HsaMir21_5p = "UAGCUUAUCAGACUGAUGUUGA";

    // Pri-miRNA from Pre_miRNA_Detection.md §7.1: 11-nt lower stem + the miR-21 stem region, junction 0.
    // The +11 Drosha cut lands at the annotated 5p start, reproducing the miRBase mature exactly.
    private const string Mir21PriMiRna =
        "CCCCCCCCCCC" + "UAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGU";

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  MIRNA-CLEAVAGE-001 — Drosha/Dicer cleavage prediction : fuzz targets
    // ═══════════════════════════════════════════════════════════════════

    #region MIRNA-CLEAVAGE-001 — BE: empty / null sequence (the nullable boundary)

    /// <summary>
    /// BE — "empty seq": the empty and null sequence. PredictDroshaDicerCleavage short-circuits
    /// IsNullOrEmpty → null BEFORE the ruler runs, so it never indexes an empty string (INV-CV-1).
    /// The result is a well-defined null ("cannot predict"), never an exception.
    /// </summary>
    [Test]
    public void Predict_EmptyOrNull_ReturnsNullNeverThrows()
    {
        var actEmpty = () => PredictDroshaDicerCleavage("", 0);
        actEmpty.Should().NotThrow("an empty sequence is a documented null result, not an error (INV-CV-1)")
                .Subject.Should().BeNull("an empty sequence has no hairpin to cleave");

        var actNull = () => PredictDroshaDicerCleavage(null!, 0);
        actNull.Should().NotThrow("a null sequence is treated as empty, not an error (INV-CV-1)")
               .Subject.Should().BeNull("a null sequence yields no cleavage prediction");

        // A non-zero junction on an empty sequence is still the IsNullOrEmpty short-circuit.
        PredictDroshaDicerCleavage("", 5).Should().BeNull("empty short-circuits before the junction range check");
    }

    #endregion

    #region MIRNA-CLEAVAGE-001 — BE: precursor too short for the +11 / +22 / +2 ruler

    /// <summary>
    /// BE — "precursor too short": THE key index-overflow boundary. The ruler runs
    /// DroshaCut5' = j + 11, MatureEnd = DroshaCut5' + 21, StarEnd = MatureEnd + 2 = j + 34, so a
    /// non-null result needs len ≥ j + 35. A precursor SHORTER than that would, in a naive builder,
    /// drive Substring off the end (IndexOutOfRange) or to a negative length. The contract instead
    /// returns null. With j = 0 the threshold is len ≥ 35; we sweep EVERY short length 1..34 (all
    /// below the floor — including the off-by-one corner 34) and assert null + no throw, then the
    /// first valid length 35 as the counterpoint that the window opens exactly there.
    /// </summary>
    [Test]
    public void Predict_PrecursorTooShort_ReturnsNullNeverIndexOverflow()
    {
        const int j = 0;
        // Floor (re-derived): StarEnd = j + 34 must be < len ⇒ len ≥ j + 35 = 35.
        for (int len = 1; len <= 34; len++)
        {
            string tooShort = new string('A', len);
            var act = () => PredictDroshaDicerCleavage(tooShort, j);
            var result = act.Should().NotThrow(
                $"a {len}-nt precursor is too short for the +11/+22/+2 ruler — the guards must " +
                "prevent any Substring overflow (no IndexOutOfRange, no negative length)").Subject;
            result.Should().BeNull(
                $"a {len}-nt precursor (< j+35 = 35) cannot place the full Drosha+Dicer ruler (INV-CV-1)");
        }

        // The window opens at exactly len = 35: the first non-null result places StarEnd at index 34.
        var firstValid = PredictDroshaDicerCleavage(new string('G', 35), j);
        firstValid.Should().NotBeNull("len = j + 35 = 35 is the first length that fits the full ruler");
        firstValid!.Value.StarEnd.Should().Be(34, "the star 3' end sits at the last index of a 35-nt precursor");
        AssertCleavageMatchesRuler(firstValid.Value, j);
    }

    /// <summary>
    /// BE — too short relative to a NON-ZERO junction: the floor scales with j (len ≥ j + 35). Even a
    /// long sequence is "too short" when the junction is placed near its 3' end so the +34 ruler
    /// overruns. We pick a 60-nt sequence and sweep junctions straddling the floor j = len − 35 = 25:
    /// j ≥ 26 overruns (null), j = 25 is the last fitting junction (non-null at exactly the tail).
    /// </summary>
    [Test]
    public void Predict_JunctionNearThreePrimeEnd_RulerOverrun_ReturnsNull()
    {
        string seq = RandomRna(60, 12345); // 60 nt; floor junction = 60 − 35 = 25
        const int floorJunction = 60 - 35; // 25

        foreach (int j in new[] { floorJunction + 1, floorJunction + 5, 40, 59 })
        {
            var act = () => PredictDroshaDicerCleavage(seq, j);
            act.Should().NotThrow($"junction {j} overruns the 60-nt sequence but must not crash (INV-CV-1)")
               .Subject.Should().BeNull(
                   $"with j = {j} the +34 ruler (StarEnd = j+34 = {j + 34}) overruns the 60-nt sequence → null");
        }

        // The last fitting junction places StarEnd at the final index (59) — non-null, ruler-correct.
        var fitting = PredictDroshaDicerCleavage(seq, floorJunction);
        fitting.Should().NotBeNull("j = len − 35 = 25 is the last junction that fits the full ruler");
        fitting!.Value.StarEnd.Should().Be(59, "at the floor junction the star 3' end sits at the last index");
        AssertCleavageMatchesRuler(fitting.Value, floorJunction);
    }

    #endregion

    #region MIRNA-CLEAVAGE-001 — BE: basal junction out of range (negative / ≥ len)

    /// <summary>
    /// BE — basal junction out of range: a negative junction or one ≥ len has no defined basal cut.
    /// The guard `basalJunction &lt; 0 || basalJunction ≥ len` returns null BEFORE the ruler runs, so
    /// neither a negative index nor an over-the-end index ever reaches a Substring (INV-CV-1). We
    /// sweep the negative corners and the ≥ len corners on a comfortably long sequence and assert
    /// null + no throw — never an IndexOutOfRange / ArgumentOutOfRange off the junction.
    /// </summary>
    [Test]
    public void Predict_JunctionOutOfRange_ReturnsNullNeverThrows()
    {
        string seq = RandomRna(80, 2026); // long enough that the ruler itself would fit for in-range j

        foreach (int j in new[] { -1, -11, -100, int.MinValue })
        {
            var act = () => PredictDroshaDicerCleavage(seq, j);
            act.Should().NotThrow($"a negative junction ({j}) must not crash — it is a documented null (INV-CV-1)")
               .Subject.Should().BeNull($"a negative basal junction ({j}) has no defined basal cut → null");
        }

        foreach (int j in new[] { seq.Length, seq.Length + 1, seq.Length + 50, int.MaxValue })
        {
            var act = () => PredictDroshaDicerCleavage(seq, j);
            act.Should().NotThrow($"a junction past the end ({j}) must not crash — documented null (INV-CV-1)")
               .Subject.Should().BeNull($"a basal junction ≥ len ({j} ≥ {seq.Length}) is out of range → null");
        }

        // The boundary junction len − 1 is in range for the range check, but the +11 ruler immediately
        // overruns ⇒ still null (via the DroshaCut5' ≥ len guard), and still no throw.
        PredictDroshaDicerCleavage(seq, seq.Length - 1).Should().BeNull(
            "the last index is in range but its +11 Drosha cut overruns the sequence → null (INV-CV-1)");
    }

    #endregion

    #region MIRNA-CLEAVAGE-001 — MC: non-ACGU characters / no-stem homopolymer

    /// <summary>
    /// MC — "non-ACGU" + "no stem": the ruler is purely POSITIONAL — it never tests base pairing — so
    /// non-canonical content (DNA 'T' → 'U', 'N', digits, punctuation, whitespace) and homopolymer
    /// (no-stem) content are normalised and measured without a throw (INV-CV-5). We verify:
    ///   (a) a DNA-spelled precursor (T→U) yields the SAME coordinates and the same T→U-normalised
    ///       mature as its RNA spelling — the ruler does not depend on the alphabet, only on length;
    ///   (b) a homopolymer "no-stem" precursor of sufficient length still returns a well-formed,
    ///       ruler-correct record (the cleavage predictor measures, it does not validate the stem);
    ///   (c) arbitrary 'N'/digit/punctuation junk of sufficient length is measured into a well-formed
    ///       record (coordinates in range, mature = the indexed span) — never a crash, never a
    ///       malformed coordinate.
    /// </summary>
    [Test]
    public void Predict_NonAcguAndNoStem_NormalizedAndMeasured_NeverCrash()
    {
        const int j = 0;

        // (a) DNA spelling of the miR-21 pri-miRNA folds to the same ruler coordinates (T→U).
        string dnaPri = Mir21PriMiRna.Replace('U', 'T').ToLowerInvariant();
        var rna = PredictDroshaDicerCleavage(Mir21PriMiRna, j);
        var dna = PredictDroshaDicerCleavage(dnaPri, j);
        rna.Should().NotBeNull();
        dna.Should().NotBeNull("a DNA-spelled precursor is normalised (T→U), not rejected (INV-CV-5)");
        dna!.Value.DroshaCut5Prime.Should().Be(rna!.Value.DroshaCut5Prime, "the ruler is alphabet-independent");
        dna.Value.MatureSequence.Should().Be(rna.Value.MatureSequence,
            "T→U + uppercasing make the DNA spelling's mature identical to the RNA spelling (INV-CV-5)");
        AssertCleavageMatchesRuler(dna.Value, j);

        // (b) A homopolymer "no-stem" precursor: the positional ruler still produces a valid record.
        foreach (char b in new[] { 'A', 'C', 'G', 'U', 'T' }) // T normalises to U
        {
            string homo = new string(b, 50); // 50 ≥ j+35 → fits the ruler
            var act = () => PredictDroshaDicerCleavage(homo, j);
            var result = act.Should().NotThrow(
                $"a no-stem all-'{b}' homopolymer must not crash the positional ruler (INV-CV-5)").Subject;
            result.Should().NotBeNull($"the positional ruler fits a 50-nt all-'{b}' precursor (no stem test)");
            AssertCleavageMatchesRuler(result!.Value, j);
        }

        // (c) Arbitrary non-ACGU junk of sufficient length: measured into a well-formed record.
        foreach (string junk in new[]
                 {
                     new string('N', 50),                                   // ambiguity homopolymer
                     new string('N', 20) + "12345!@#$ \t%^&*" + new string('N', 20),
                     "GCGC" + new string('N', 40) + "GCGC",
                 })
        {
            var act = () => PredictDroshaDicerCleavage(junk, j);
            var result = act.Should().NotThrow(
                $"junk content (len {junk.Length}) must not crash the cleavage ruler (INV-CV-5)").Subject;
            // The ruler is positional, so a sufficiently long junk string still yields a valid record.
            if (result is not null)
                AssertCleavageMatchesRuler(result.Value, j);
        }
    }

    #endregion

    #region MIRNA-CLEAVAGE-001 — discriminative pin: the miRBase hsa-miR-21-5p worked example

    /// <summary>
    /// The discriminative pin (INV-CV-3/4): the predictor is NOT a no-op. Feeding the doc's §7.1
    /// pri-miRNA — an 11-nt lower stem prefixed to the miR-21 stem region with j = 0 — places the +11
    /// Drosha cut at the annotated 5p start and reproduces miRBase hsa-miR-21-5p (MIMAT0000076,
    /// "UAGCUUAUCAGACUGAUGUUGA") EXACTLY. We assert the hand-derived coordinates (DroshaCut5' = 11,
    /// mature 22 nt, StarEnd = MatureEnd + 2) and the literal miRBase mature sequence. A predictor
    /// that returned the wrong offset — or a test that read the offset off the code — would be invalid.
    /// </summary>
    [Test]
    public void Predict_Mir21PriMiRna_ReproducesMiRBaseMatureExactly()
    {
        const int j = 0;
        var result = PredictDroshaDicerCleavage(Mir21PriMiRna, j);
        result.Should().NotBeNull("the 69-nt miR-21 pri-miRNA fits the full Drosha+Dicer ruler");
        var c = result!.Value;

        // Hand-derived ruler coordinates (computed from the doc, NOT from the code).
        c.DroshaCut5Prime.Should().Be(11, "Drosha cuts 11 bp from the basal junction j = 0 (Han 2006; INV-CV-3)");
        c.MatureStart.Should().Be(11, "the 5p mature begins at the +11 Drosha cut (INV-CV-3)");
        c.MatureEnd.Should().Be(32, "Dicer 5'-counting fixes the mature at 22 nt: 11 + 21 = 32 (Park 2011; INV-CV-3)");
        c.StarEnd.Should().Be(34, "StarEnd = MatureEnd + 2-nt overhang = 32 + 2 (Lee 2003; INV-CV-3)");
        c.StarStart.Should().Be(13, "StarStart = StarEnd − 21 = 34 − 21 (INV-CV-3)");

        // The literal miRBase mature — the proof the ruler lands on the real biological boundary.
        c.MatureSequence.Should().Be(HsaMir21_5p,
            "the +11 Drosha cut + 22-nt Dicer ruler reproduces miRBase hsa-miR-21-5p exactly (INV-CV-4)");
        c.MatureSequence.Length.Should().Be(22, "the mature is exactly 22 nt (Dicer 5'-counting)");
        c.ThreePrimeOverhang.Should().Be(2, "each RNase III cut leaves a 2-nt 3' overhang (INV-CV-4)");

        // And the full universal record contract holds.
        AssertCleavageMatchesRuler(c, j);
    }

    /// <summary>
    /// Discriminative pin over MANY valid (sequence, junction) pairs: across fixed seeds, lengths and
    /// in-range junctions that fit the ruler (len ≥ j + 35), the predictor NEVER throws and ALWAYS
    /// returns a record whose coordinates equal the independently re-derived ruler offsets, all in
    /// [0,len), with a 22-nt mature and a 2-nt overhang. This proves the harness asserts against a
    /// predictor that really measures the ruler — not a no-op that always returns null.
    /// </summary>
    [Test]
    public void Predict_RandomValidPrecursors_AlwaysRulerCorrectAndCrashFree()
    {
        foreach (int seed in new[] { 1, 7, 42, 2026 })
        {
            foreach (int len in new[] { 35, 60, 90, 120 })
            {
                string raw = RandomRna(len, seed);
                // Junctions that fit the ruler: 0 .. len − 35 (the re-derived non-null window).
                int floorJunction = len - 35; // last junction whose +34 ruler still fits (StarEnd = j+34 < len)
                foreach (int j in new[] { 0, 1, floorJunction / 2, floorJunction })
                {
                    if (j < 0 || j > floorJunction) continue; // skip junctions that overrun this length
                    var act = () => PredictDroshaDicerCleavage(raw, j);
                    var result = act.Should().NotThrow(
                        $"a fitting (len {len}, j {j}, seed {seed}) precursor must not crash the ruler").Subject;
                    result.Should().NotBeNull(
                        $"len {len} ≥ j+35 = {j + 35} fits the full ruler (seed {seed})");
                    AssertCleavageMatchesRuler(result!.Value, j);
                }
            }
        }
    }

    #endregion
}
