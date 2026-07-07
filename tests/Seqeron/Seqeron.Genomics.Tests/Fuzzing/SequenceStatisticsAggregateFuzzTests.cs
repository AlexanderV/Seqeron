namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Statistics-area aggregate AMINO-ACID composition unit.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// and no *unhandled* runtime exception (IndexOutOfRangeException,
/// NullReferenceException, DivideByZeroException, OverflowException, …). Every
/// input must result in EITHER a well-defined, theory-correct value, OR a
/// *documented, intentional* validation exception. A raw runtime exception or a
/// hang on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: SEQ-STATS-001 — Sequence Composition Statistics (Statistics)
/// Checklist: docs/checklists/03_FUZZING.md, row 127.
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — empty, all-N, lowercase, mixed alphabet.
/// — docs/checklists/03_FUZZING.md §Description (strategy codes); row 127 targets.
///
/// ───────────────────────────────────────────────────────────────────────────
/// SCOPING vs row 121 (SEQ-COMPOSITION-001) and row 128 (SEQ-SUMMARY-001)
/// ───────────────────────────────────────────────────────────────────────────
/// SEQ-STATS-001 is documented as a CONSOLIDATED unit that SHARES the nucleotide
/// method CalculateNucleotideComposition with row 121 SEQ-COMPOSITION-001
/// (docs/algorithms/Sequence_Composition/Sequence_Composition.md §intro, lines
/// 21–23; CANONICAL_MAP.md "Consolidated / Same canonical fixture"). Row 121
/// (SequenceCompositionStatFuzzTests.cs) already fuzzes that method exhaustively
/// — counts, GcContent/AtContent, GcSkew/AtSkew. Row 128 SEQ-SUMMARY-001 owns the
/// nucleotide aggregator SummarizeNucleotideSequence.
///
/// This file therefore tests the DISTINCT aggregate-statistics surface that
/// SEQ-STATS-001's OWN algorithm doc lists as an entry point and that NO other
/// fuzz row covers: the PROTEIN composition aggregate record.
///   — docs/algorithms/Sequence_Composition/Sequence_Composition_Statistics.md
///     §5.1: "SequenceStatistics.CalculateAminoAcidComposition(string): protein
///     residue counts/ratios". §1 describes the unit as also exposing "a protein-
///     composition counterpart".
/// The per-amino-acid SCALAR metrics (MW row 124, pI row 125, hydrophobicity
/// row 123) test those scalars in isolation; this row asserts the AGGREGATE
/// AminoAcidComposition RECORD — its Counts map, Length, and the bundled
/// ChargedResidueRatio / AromaticResidueRatio / MW / pI / Hydrophobicity fields
/// — field by field, which is complementary to all of them.
///
/// ───────────────────────────────────────────────────────────────────────────
/// The aggregate amino-acid composition contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// API entry: SequenceStatistics.CalculateAminoAcidComposition(string)
///   (src/Seqeron/Algorithms/Seqeron.Genomics.Analysis/SequenceStatistics.cs
///    lines 99–139), returning an AminoAcidComposition record struct with fields
///    Length, Counts, MolecularWeight, IsoelectricPoint, Hydrophobicity,
///    ChargedResidueRatio, AromaticResidueRatio.
///
/// Documented behaviour (Sequence_Composition_Statistics.md §3.3, §5.1; source):
///   • null / empty input → degenerate composition: Length 0, EMPTY Counts,
///     MolecularWeight 0, IsoelectricPoint 7.0 (neutral sentinel — pI is undefined
///     for a zero-length protein; SequenceStatistics.cs lines 277–278, 103),
///     Hydrophobicity 0, ChargedResidueRatio 0, AromaticResidueRatio 0. No throw,
///     no DivideByZero.
///   • Counting is case-insensitive: input is upper-cased (ToUpperInvariant,
///     line 107) before classification.
///   • Counts holds every LETTER (char.IsLetter, line 109) keyed by its upper-case
///     form — including non-standard codes (B, J, O, U, X, Z, N, …). Non-letters
///     (digits, whitespace, punctuation, control chars) are NOT counted.
///   • Length = Σ Counts.Values (line 115) — the recognized-LETTER count, which can
///     be SMALLER than the input string length when junk is present.
///   • MolecularWeight / Hydrophobicity are computed over only the 20 STANDARD
///     residues (lines 116, 118); unknown letters contribute nothing.
///   • ChargedResidueRatio = (#D+#E+#K+#R+#H)/Length, 0 when Length = 0 (line 124).
///   • AromaticResidueRatio = (#F+#Y+#W)/Length, 0 when Length = 0 (line 129).
///
/// Invariants pinned (Sequence_Composition_Statistics.md §2.4 INV-03 analogue;
///   well-formedness of the aggregate record):
///   • Length = Σ Counts.Values, and Length ≥ 0.
///   • Every Counts key is an upper-case letter; every value ≥ 1.
///   • 0 ≤ ChargedResidueRatio ≤ 1 and 0 ≤ AromaticResidueRatio ≤ 1.
///   • When Length = 0: ratios 0, MW 0, Hydrophobicity 0, pI = 7.0 sentinel.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class SequenceStatisticsAggregateFuzzTests
{
    #region Helpers

    private const double Tolerance = 1e-9;

    /// <summary>The pI neutral sentinel returned for empty / no-ionizable input
    /// (SequenceStatistics.cs NeutralPhDefault = 7.0).</summary>
    private const double NeutralPi = 7.0;

    /// <summary>The 20 standard one-letter amino-acid codes.</summary>
    private const string StandardResidues = "ACDEFGHIKLMNPQRSTVWY";

    /// <summary>Generates a random string of arbitrary BMP code points (0x0000–0xFFFF),
    /// spanning control characters, the null byte, lone surrogate halves, unicode
    /// letters and digits — random-byte fuzz fodder for the classifier.</summary>
    private static string RandomBmpChars(Random rng, int length)
    {
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = (char)rng.Next(0x0000, 0x10000);
        return new string(chars);
    }

    /// <summary>
    /// Asserts the universal well-formedness contract that must hold for ANY input:
    /// the Length = Σ counts identity, well-keyed Counts, bounded ratios, and the
    /// zero-length guards (no DivideByZero on any bundled field).
    /// </summary>
    private static void AssertWellFormed(SequenceStatistics.AminoAcidComposition c)
    {
        c.Length.Should().BeGreaterThanOrEqualTo(0);

        // Length must equal the sum of all recognized-letter counts.
        c.Counts.Values.Sum().Should().Be(c.Length,
            "Length is defined as Σ Counts.Values (recognized letters only)");

        // Every key is an upper-case letter; every count is a positive occurrence.
        foreach (var (key, count) in c.Counts)
        {
            char.IsLetter(key).Should().BeTrue($"'{key}' must be a letter to be a Counts key");
            key.Should().Be(char.ToUpperInvariant(key), "keys are upper-cased before counting");
            count.Should().BeGreaterThan(0, "a present key has at least one occurrence");
        }

        // Ratios are well-defined fractions in [0,1].
        c.ChargedResidueRatio.Should().BeInRange(0.0, 1.0);
        c.AromaticResidueRatio.Should().BeInRange(0.0, 1.0);

        // pI is always within the standard pH window [0,14]. NOTE: the 7.0 sentinel
        // is returned ONLY for a null/EMPTY input STRING (IsNullOrEmpty short-circuit,
        // SequenceStatistics.cs line 297); a non-empty string with NO ionizable residue
        // still runs the bisection and converges to the termini-only pI (≈6.1), so
        // Length == 0 does NOT imply pI == 7.0 (Length counts letters, not all chars).
        c.IsoelectricPoint.Should().BeInRange(0.0, 14.0);
        c.MolecularWeight.Should().BeGreaterThanOrEqualTo(0.0);

        if (c.Length == 0)
        {
            // Zero-denominator guards: nothing divides by Length 0.
            c.ChargedResidueRatio.Should().Be(0.0, "no residues ⇒ charged ratio guarded to 0");
            c.AromaticResidueRatio.Should().Be(0.0, "no residues ⇒ aromatic ratio guarded to 0");
            c.MolecularWeight.Should().Be(0.0, "no residues ⇒ MW 0");
            c.Hydrophobicity.Should().Be(0.0, "no residues ⇒ GRAVY 0");
        }
    }

    #endregion

    // ═══════════════════════════════════════════════════════════════════
    //  SEQ-STATS-001 — aggregate amino-acid composition : fuzz targets (BE)
    // ═══════════════════════════════════════════════════════════════════

    #region Positive sanity — hand-computed exact result

    /// <summary>
    /// Positive baseline (not a boundary): a known protein must yield the documented
    /// aggregate fields EXACTLY. "KKRRDDEEHHFFYYWWAAAA" → 20 residues:
    /// charged D2 E2 K2 R2 H2 = 10/20 = 0.5; aromatic F2 Y2 W2 = 6/20 = 0.3.
    /// Length = 20 = Σ counts; every count is 2 except A which is 4. Confirms the
    /// suite asserts the BUSINESS contract, not just non-throwing.
    /// — Sequence_Composition_Statistics.md §5.1; SequenceStatistics.cs lines 120–129.
    /// </summary>
    [Test]
    public void AminoAcidComposition_KnownProtein_MatchesHandComputedFields()
    {
        var c = SequenceStatistics.CalculateAminoAcidComposition("KKRRDDEEHHFFYYWWAAAA");

        c.Length.Should().Be(20);
        c.Counts['A'].Should().Be(4);
        c.Counts['K'].Should().Be(2);
        c.Counts['R'].Should().Be(2);
        c.Counts['D'].Should().Be(2);
        c.Counts['E'].Should().Be(2);
        c.Counts['H'].Should().Be(2);
        c.Counts['F'].Should().Be(2);
        c.Counts['Y'].Should().Be(2);
        c.Counts['W'].Should().Be(2);

        // charged = D2+E2+K2+R2+H2 = 10; aromatic = F2+Y2+W2 = 6.
        c.ChargedResidueRatio.Should().BeApproximately(10.0 / 20.0, Tolerance);
        c.AromaticResidueRatio.Should().BeApproximately(6.0 / 20.0, Tolerance);

        c.MolecularWeight.Should().BeGreaterThan(0);
        c.IsoelectricPoint.Should().BeInRange(0.0, 14.0);

        AssertWellFormed(c);
    }

    #endregion

    #region BE — Boundary: empty / null

    /// <summary>
    /// BE: the empty string is the lower size boundary. Documented as a degenerate
    /// composition — empty Counts, Length 0, MW 0, pI 7.0 sentinel, ratios 0; NO
    /// DivideByZero on the Length-0 ratio denominators, no exception.
    /// — Sequence_Composition_Statistics.md §3.3; SequenceStatistics.cs lines 101–104.
    /// </summary>
    [Test]
    public void AminoAcidComposition_EmptyString_IsDegenerateAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.CalculateAminoAcidComposition(string.Empty);

        act.Should().NotThrow("the empty string is a defined boundary, not an error");

        var c = act();
        c.Length.Should().Be(0);
        c.Counts.Should().BeEmpty();
        c.MolecularWeight.Should().Be(0.0);
        c.IsoelectricPoint.Should().Be(NeutralPi);
        c.Hydrophobicity.Should().Be(0.0);
        c.ChargedResidueRatio.Should().Be(0.0);
        c.AromaticResidueRatio.Should().Be(0.0);
        AssertWellFormed(c);
    }

    /// <summary>
    /// BE: null is treated identically to empty (IsNullOrEmpty short-circuit,
    /// SequenceStatistics.cs line 101) — degenerate composition, no NullReference.
    /// </summary>
    [Test]
    public void AminoAcidComposition_Null_IsDegenerateAndDoesNotThrow()
    {
        var act = () => SequenceStatistics.CalculateAminoAcidComposition(null!);

        act.Should().NotThrow("null is documented as 'no sequence', not an error");

        var c = act();
        c.Length.Should().Be(0);
        c.Counts.Should().BeEmpty();
        c.IsoelectricPoint.Should().Be(NeutralPi);
        AssertWellFormed(c);
    }

    /// <summary>
    /// BE: a string of ONLY non-letters (digits, whitespace, punctuation) has no
    /// recognized residue, so the count surface is empty: Length 0, empty Counts,
    /// MW 0, ratios 0 — even though the input string is non-empty. Guards the
    /// "Length = recognized letters, not string length" contract at the boundary.
    /// Note: the pI is NOT the 7.0 empty sentinel here — the string is non-empty so
    /// CalculateIsoelectricPoint bisects and returns the termini-only pI (in [0,14],
    /// ≈6.1), since the 7.0 sentinel is reserved for a null/empty input STRING.
    /// — SequenceStatistics.cs line 109 (only char.IsLetter counted), lines 297–326.
    /// </summary>
    [Test]
    public void AminoAcidComposition_OnlyNonLetters_HasEmptyCountSurface()
    {
        var c = SequenceStatistics.CalculateAminoAcidComposition("12345  \t.,;-\0");

        c.Length.Should().Be(0, "no character is a letter");
        c.Counts.Should().BeEmpty();
        c.MolecularWeight.Should().Be(0.0);
        c.Hydrophobicity.Should().Be(0.0);
        c.ChargedResidueRatio.Should().Be(0.0);
        c.AromaticResidueRatio.Should().Be(0.0);
        c.IsoelectricPoint.Should().BeInRange(0.0, 14.0,
            "non-empty string ⇒ bisected termini-only pI, not the empty 7.0 sentinel");
        AssertWellFormed(c);
    }

    #endregion

    #region BE — Boundary: all-N

    /// <summary>
    /// BE: every character is 'N'. In the PROTEIN alphabet 'N' is asparagine — a
    /// standard residue — so all-N is a valid poly-asparagine: Counts = {'N': k},
    /// Length = k, MW &gt; 0 (N is in the weight table), but NO charged (D/E/K/R/H)
    /// and NO aromatic (F/Y/W) residues ⇒ both ratios are exactly 0 with NO
    /// DivideByZero. This contrasts with the NUCLEOTIDE meaning of N (ambiguous
    /// base) tested in row 121 — the aggregate amino-acid surface treats N as Asn.
    /// — SequenceStatistics.cs lines 107–124; AminoAcidWeights contains 'N'.
    /// </summary>
    [TestCase(1)]
    [TestCase(2)]
    [TestCase(50)]
    [TestCase(1000)]
    public void AminoAcidComposition_AllN_IsPolyAsparagine_RatiosZero(int length)
    {
        string seq = new string('N', length);

        var c = SequenceStatistics.CalculateAminoAcidComposition(seq);

        c.Length.Should().Be(length, "every N is a counted residue (asparagine)");
        c.Counts.Should().HaveCount(1);
        c.Counts['N'].Should().Be(length);
        c.MolecularWeight.Should().BeGreaterThan(0, "asparagine has a defined mass");
        c.ChargedResidueRatio.Should().Be(0.0, "N is not D/E/K/R/H");
        c.AromaticResidueRatio.Should().Be(0.0, "N is not F/Y/W");
        AssertWellFormed(c);
    }

    /// <summary>
    /// BE: lowercase 'n' must be classified identically to 'N' (ToUpperInvariant,
    /// line 107) — Counts key is the upper-case 'N'. Guards a case bug in counting.
    /// </summary>
    [Test]
    public void AminoAcidComposition_LowercaseN_CountedAsUppercaseN()
    {
        var c = SequenceStatistics.CalculateAminoAcidComposition("nnnn");

        c.Length.Should().Be(4, "lowercase n is upper-cased before counting");
        c.Counts.Should().ContainKey('N');
        c.Counts['N'].Should().Be(4);
        c.Counts.Should().NotContainKey('n', "keys are upper-cased");
        AssertWellFormed(c);
    }

    /// <summary>
    /// BE: a string of an UNKNOWN letter (e.g. 'X' = any residue, 'Z' = Glx) is still
    /// COUNTED (it is a letter) and contributes to Length and Counts, but contributes
    /// NO mass (not in the 20-residue weight table) and is neither charged nor
    /// aromatic ⇒ MW 0, ratios 0, but Length &gt; 0. This separates "recognized
    /// letter" (counted) from "standard residue" (weighed/scored).
    /// — SequenceStatistics.cs lines 109, 116 (MW skips unknown); §5.1.
    /// </summary>
    [TestCase('X')]
    [TestCase('Z')]
    [TestCase('B')]
    [TestCase('J')]
    [TestCase('O')]
    public void AminoAcidComposition_AllUnknownLetter_CountedButUnweighted(char letter)
    {
        string seq = new string(letter, 10);

        var c = SequenceStatistics.CalculateAminoAcidComposition(seq);

        c.Length.Should().Be(10, "an unknown letter is still a counted residue");
        c.Counts[char.ToUpperInvariant(letter)].Should().Be(10);
        c.MolecularWeight.Should().Be(0.0, $"'{letter}' is not a standard residue");
        c.ChargedResidueRatio.Should().Be(0.0);
        c.AromaticResidueRatio.Should().Be(0.0);
        AssertWellFormed(c);
    }

    #endregion

    #region BE — Boundary: lowercase / mixed case

    /// <summary>
    /// BE: counting is case-insensitive. A lowercase protein must produce the
    /// IDENTICAL aggregate (counts, Length, ratios, MW, pI, hydrophobicity) to its
    /// upper-case form (ToUpperInvariant, line 107). Guards against a case bug that
    /// would split 'a' and 'A' into separate Counts keys or miss standard residues.
    /// — Sequence_Composition_Statistics.md §6.1 (mixed case → identical result).
    /// </summary>
    [TestCase("mkvlwa")]
    [TestCase("MkVlWa")]
    [TestCase("kkrrddee")]
    [TestCase("ffyyww")]
    public void AminoAcidComposition_LowercaseOrMixedCase_EqualsUppercase(string seq)
    {
        var lower = SequenceStatistics.CalculateAminoAcidComposition(seq);
        var upper = SequenceStatistics.CalculateAminoAcidComposition(seq.ToUpperInvariant());

        lower.Length.Should().Be(upper.Length);
        lower.Counts.Should().BeEquivalentTo(upper.Counts, "case must not affect the count map");
        lower.MolecularWeight.Should().BeApproximately(upper.MolecularWeight, Tolerance);
        lower.IsoelectricPoint.Should().BeApproximately(upper.IsoelectricPoint, Tolerance);
        lower.Hydrophobicity.Should().BeApproximately(upper.Hydrophobicity, Tolerance);
        lower.ChargedResidueRatio.Should().BeApproximately(upper.ChargedResidueRatio, Tolerance);
        lower.AromaticResidueRatio.Should().BeApproximately(upper.AromaticResidueRatio, Tolerance);

        // Keys are upper-case even from lowercase input.
        lower.Counts.Keys.Should().OnlyContain(k => k == char.ToUpperInvariant(k));
        AssertWellFormed(lower);
    }

    /// <summary>
    /// BE: a single lowercase 'k' is a charged residue (ChargedResidueRatio 1.0),
    /// proving lowercase does not silently fail to be recognized as lysine.
    /// </summary>
    [Test]
    public void AminoAcidComposition_SingleLowercaseChargedResidue_RatioIsOne()
    {
        var c = SequenceStatistics.CalculateAminoAcidComposition("k");

        c.Length.Should().Be(1);
        c.Counts['K'].Should().Be(1);
        c.ChargedResidueRatio.Should().BeApproximately(1.0, Tolerance, "lysine is charged");
        AssertWellFormed(c);
    }

    #endregion

    #region BE — Boundary: mixed alphabet (DNA + protein + junk)

    /// <summary>
    /// BE: a mixed-alphabet input combining DNA bases, protein residues, degenerate
    /// codes and junk. EVERY LETTER is counted (case-folded) toward Length and
    /// Counts; non-letters are excluded. The ratios are computed over the
    /// recognized-letter Length, NOT the input string length, and stay well-defined.
    ///
    /// "K1R d-e F\tACGT" hand-trace (upper-cased letters only):
    ///   letters = K R D E F A C G T  → Length 9, each count 1.
    ///   non-letters '1',' ','-',' ','\t' → excluded.
    ///   charged = K+R+D+E = 4 ⇒ 4/9; aromatic = F = 1 ⇒ 1/9.
    /// Confirms field-by-field handling of a hostile mixed alphabet with no crash and
    /// no DivideByZero. — SequenceStatistics.cs lines 107–129; §3.3.
    /// </summary>
    [Test]
    public void AminoAcidComposition_MixedAlphabet_CountsLettersOnly_RatiosOverResidueLength()
    {
        const string seq = "K1R d-e F\tACGT";

        var c = SequenceStatistics.CalculateAminoAcidComposition(seq);

        // Letters present (upper-cased): K R D E F A C G T — each once.
        c.Length.Should().Be(9, "only the 9 letters are counted; digits/space/dash/tab excluded");
        foreach (char ch in "KRDEFACGT")
            c.Counts[ch].Should().Be(1, $"'{ch}' appears once");
        c.Counts.Values.Sum().Should().Be(9);

        // charged = K+R+D+E = 4 (A,C,G,T,F are not charged); aromatic = F = 1.
        c.ChargedResidueRatio.Should().BeApproximately(4.0 / 9.0, Tolerance);
        c.AromaticResidueRatio.Should().BeApproximately(1.0 / 9.0, Tolerance);

        AssertWellFormed(c);
    }

    /// <summary>
    /// BE: DNA-only input fed to the PROTEIN aggregate. A/C/G/T are all valid protein
    /// letters (Ala/Cys/Gly/Thr), so a DNA string is interpreted as a peptide:
    /// counted, weighed, Length = string length (all letters), ratios well-defined
    /// (none of A/C/G/T is charged ⇒ 0; none aromatic ⇒ 0). No crash. This documents
    /// the (mis)use boundary where a nucleotide string crosses into the protein path.
    /// </summary>
    [Test]
    public void AminoAcidComposition_DnaStringAsProtein_AllLettersCounted_RatiosZero()
    {
        var c = SequenceStatistics.CalculateAminoAcidComposition("ACGTACGTACGT");

        c.Length.Should().Be(12, "A/C/G/T are all protein letters");
        c.Counts['A'].Should().Be(3);
        c.Counts['C'].Should().Be(3);
        c.Counts['G'].Should().Be(3);
        c.Counts['T'].Should().Be(3);
        c.MolecularWeight.Should().BeGreaterThan(0, "Ala/Cys/Gly/Thr are standard residues");
        c.ChargedResidueRatio.Should().Be(0.0, "none of A/C/G/T is charged");
        c.AromaticResidueRatio.Should().Be(0.0, "none of A/C/G/T is aromatic");
        AssertWellFormed(c);
    }

    #endregion

    #region BE — Random / RB fuzz: never throws, always well-formed

    /// <summary>
    /// BE/RB: a large batch of arbitrary BMP strings (control chars, null byte, lone
    /// surrogate halves, unicode letters/digits) must NEVER throw and must ALWAYS
    /// produce a well-formed aggregate satisfying every invariant. This is the core
    /// fuzz guarantee: no DivideByZero on the ratio denominators, no IndexOutOfRange,
    /// no NullReference, no overflow on garbage. pI bisection must always terminate.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void AminoAcidComposition_RandomGarbageStrings_NeverThrow_AlwaysWellFormed()
    {
        var rng = new Random(20260620);

        for (int iteration = 0; iteration < 2000; iteration++)
        {
            int len = rng.Next(0, 200);
            string input = RandomBmpChars(rng, len);

            SequenceStatistics.AminoAcidComposition c = default;
            var act = () => c = SequenceStatistics.CalculateAminoAcidComposition(input);

            act.Should().NotThrow($"garbage input (len {len}) must never crash the aggregate");
            AssertWellFormed(c);
        }
    }

    /// <summary>
    /// BE: a randomly built STANDARD-residue protein must have its counts match an
    /// independent re-count, Length = Σ counts = string length (all letters), and the
    /// charged/aromatic ratios match a hand re-derivation over the same counts.
    /// Cross-checks the aggregate against a simple oracle over many shapes.
    /// </summary>
    [Test]
    [CancelAfter(30000)]
    public void AminoAcidComposition_RandomStandardProteins_FieldsMatchOracle()
    {
        var rng = new Random(424242);

        for (int iteration = 0; iteration < 1000; iteration++)
        {
            int len = rng.Next(1, 300);
            var chars = new char[len];
            for (int i = 0; i < len; i++)
                chars[i] = StandardResidues[rng.Next(StandardResidues.Length)];
            string seq = new string(chars);

            var c = SequenceStatistics.CalculateAminoAcidComposition(seq);

            c.Length.Should().Be(len, "all residues are standard letters");
            foreach (char aa in StandardResidues)
            {
                int expected = seq.Count(ch => ch == aa);
                if (expected == 0)
                    c.Counts.Should().NotContainKey(aa);
                else
                    c.Counts[aa].Should().Be(expected);
            }

            int charged = seq.Count(ch => "DEKRH".Contains(ch));
            int aromatic = seq.Count(ch => "FYW".Contains(ch));
            c.ChargedResidueRatio.Should().BeApproximately((double)charged / len, Tolerance);
            c.AromaticResidueRatio.Should().BeApproximately((double)aromatic / len, Tolerance);

            AssertWellFormed(c);
        }
    }

    #endregion
}
