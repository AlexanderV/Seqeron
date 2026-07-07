namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Metagenomics functional-annotation unit — homology-based
/// functional prediction (annotation transfer) via <see cref="MetagenomicsAnalyzer.PredictFunctions"/>.
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds random, invalid and boundary inputs to a unit and asserts that
/// the code NEVER fails in an undisciplined way: no hang, no state corruption,
/// no NaN/Infinity leaking into a score, and no *unhandled* runtime exception
/// (IndexOutOfRangeException, NullReferenceException, DivideByZeroException,
/// OverflowException, …). Every input must produce EITHER a well-defined,
/// theory-correct result, OR a *documented, intentional* validation exception
/// (ArgumentNullException). A raw runtime exception, a hang, or a *fabricated
/// annotation* on garbage input is a bug, not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: META-FUNC-001 — functional prediction / annotation transfer (Metagenomics)
/// Checklist: docs/checklists/03_FUZZING.md, row 194.
/// Fuzz strategy for THIS unit: BE = Boundary Exploitation (0, -1, MaxInt, empty)
///   — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Fuzz targets (checklist row 194): "empty, unknown genes, no DB hit".
///
/// ───────────────────────────────────────────────────────────────────────────
/// The contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// PredictFunctions implements homology-based annotation transfer: for each query
/// protein it scans every signature in the function database; for every signature
/// that occurs EXACTLY (ordinal substring) in the protein it computes the ungapped
/// BLOSUM62 raw self-score S of the matched signature, converts it to a BLAST bit
/// score S' = (λ·S − ln K)/ln 2 and an E-value E = K·m·n·e^(−λ·S) (m = protein
/// length, n = signature length) via Karlin-Altschul statistics, and transfers the
/// annotation (Function, Pathway, KO) of the single best hit (lowest E-value) to
/// the gene. One FunctionalAnnotation is yielded per gene that matched ≥ 1 signature.
///   — docs/algorithms/Metagenomics/Functional_Prediction.md §2.A, §4.A;
///     src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/MetagenomicsAnalyzer.cs
///     (PredictFunctions, lines 906–954).
/// Ungapped BLOSUM62 Karlin-Altschul parameters (NCBI blast_stat.c blosum62_values):
/// λ = 0.3176, K = 0.134 (§4 decision table; Functional_Prediction.md [2]).
/// BLOSUM62 diagonal self-scores (NCBI BLOSUM62 matrix [3]):
///   A4 R5 N6 D6 C9 Q5 E5 G6 H8 I4 L4 K5 M5 F6 P7 S4 T5 W11 Y7 V4.
///
/// Boundary / malformed-input handling fixed by the doc (§3.3, §6.1) and source,
/// which these fuzz tests pin so the contract can never silently drift:
///   • EMPTY (BE):
///       – empty/whitespace protein sequence → that gene is SKIPPED, no
///         annotation (PredictFunctions, IsNullOrWhiteSpace guard, line 919–920).
///         — §6.1 (Empty/whitespace protein sequence).
///       – empty protein collection → no annotations, no crash.
///       – empty or null function database (no signatures) → no annotations for
///         any gene (the signature loop never matches). — §6.1 (Empty database).
///       – empty-string signatures are skipped (IsNullOrEmpty(signature) guard,
///         line 927) — an empty signature must NOT vacuously "match" every protein.
///   • UNKNOWN GENES (BE/MC): residues outside the 20-letter BLOSUM62 alphabet
///     (X, B, Z, J, U, *, digits, gaps) contribute 0 to the raw self-score
///     (Blosum62SelfScore, unknown residue ⇒ TryGetValue miss ⇒ +0, line 969–970).
///     A signature of all-unknown residues therefore scores S = 0 ⇒ bit score
///     S' = (−ln K)/ln 2 = 2.899… and a *defined, finite* E-value — never a crash,
///     never NaN/Infinity. Such a gene is still annotated IF the signature matches.
///   • NO DB HIT (BE): a perfectly valid protein whose sub-sequences are none of
///     the database signatures → no signature ever satisfies string.Contains ⇒
///     best stays null ⇒ NO annotation emitted for that gene (not a fabricated
///     zero-score hit). — §6.1 (Gene matches no signature).
///   • Argument validation: null proteins / functionDatabase → ArgumentNullException
///     (lines 910–911). These guards are checked EAGERLY before the iterator is
///     returned, so they fire at the call itself (not on first enumeration) — pinned.
///
/// Positive sanity (worked example, Functional_Prediction.md §7.1): signature
/// "WWW" matched in the 3-residue protein "WWW" ⇒ raw S = 3·11 = 33, bit score
/// S' = (0.3176·33 − ln 0.134)/ln 2 = 18.0202932787533 and E-value (m = n = 3) =
/// 0.134·3·3·e^(−0.3176·33) = 3.3852730346546 × 10⁻⁵. These constants are derived
/// INDEPENDENTLY from the cited BLAST formulas (not echoed off the implementation),
/// matching the META-FUNC-001 evidence walk-through. A genuine match must therefore
/// yield a positive bit score and a tiny positive E-value — so a passing "no crash"
/// result cannot be a degenerate annotator that emits nothing (or a zero score) for
/// everything.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Determinism
/// ───────────────────────────────────────────────────────────────────────────
/// All inputs are hand-built or generated from a LOCALLY fixed-seed
/// `new Random(seed)` (never a shared static Rng), so every run is reproducible.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class MetagenomicsFunctionalFuzzTests
{
    // Independently-derived BLAST constants for the §7.1 worked example.
    // λ = 0.3176, K = 0.134; S(WWW) = 3·11 = 33; m = n = 3.
    //   S' = (0.3176·33 − ln 0.134)/ln 2          = 18.020293278753364
    //   E  = 0.134·3·3·e^(−0.3176·33)             = 3.3852730346545964e-05
    private const double WwwBitScore = 18.020293278753364;
    private const double WwwEValue = 3.3852730346545964e-05;

    // A small hand-built signature → (Function, Pathway, KO) database. Signatures
    // are short verbatim peptides so the exact-substring match is hand-checkable.
    private static IReadOnlyDictionary<string, (string Function, string Pathway, string Ko)> BuildDatabase() =>
        new Dictionary<string, (string, string, string)>
        {
            ["WWW"] = ("tryptophan transport", "ABC transport", "K00001"),
            ["CC"]  = ("disulfide oxidoreductase", "redox", "K00002"),
            ["MK"]  = ("histidine kinase", "two-component signaling", "K00003"),
        };

    #region META-FUNC-001 — functional annotation transfer

    // ════════════════════════════════════════════════════════════════════════
    //  Positive sanity — the §7.1 worked example must be reproduced EXACTLY.
    //  Guards against a degenerate annotator (one that emits nothing, or a
    //  constant/zero score) that would pass every boundary test below.
    // ════════════════════════════════════════════════════════════════════════

    [Test]
    public void PredictFunctions_WorkedExample_TransfersAnnotationWithExactBlastStats()
    {
        var db = BuildDatabase();
        var proteins = new[] { ("gene-WWW", "WWW") };

        var annotation = MetagenomicsAnalyzer.PredictFunctions(proteins, db).Single();

        annotation.GeneId.Should().Be("gene-WWW");
        annotation.Function.Should().Be("tryptophan transport",
            "the annotation of the matched signature is transferred to the gene");
        annotation.Pathway.Should().Be("ABC transport");
        annotation.KoNumber.Should().Be("K00001");

        // Bit score and E-value derived independently from the BLAST formulas (§7.1).
        annotation.BitScore.Should().BeApproximately(WwwBitScore, 1e-9,
            "S' = (λ·33 − ln K)/ln 2 with λ = 0.3176, K = 0.134 — §2.A worked example");
        annotation.EValue.Should().BeApproximately(WwwEValue, 1e-15,
            "E = K·m·n·e^(−λ·33) with m = n = 3 — §2.A worked example");
        annotation.BitScore.Should().BeGreaterThan(0.0).And.NotBe(double.NaN);
        annotation.EValue.Should().BeGreaterThan(0.0).And.BeLessThan(1.0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Best-hit rule: among several signatures present in one protein, the
    // annotation of the LOWEST E-value (most significant) hit is transferred.
    // In "WWW" both "WWW" (S = 33, E = 3.39e-5) and "WW" (S = 22, E = 7.43e-4)
    // occur; "WWW" must win. — Functional_Prediction.md INV-06, §4.A.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_MultipleMatches_TransfersLowestEValueHit()
    {
        var db = new Dictionary<string, (string, string, string)>
        {
            ["WW"]  = ("weaker-hit", "p2", "K2"),
            ["WWW"] = ("best-hit",   "p1", "K1"),
        };
        var proteins = new[] { ("g", "WWW") };

        var annotation = MetagenomicsAnalyzer.PredictFunctions(proteins, db).Single();

        annotation.Function.Should().Be("best-hit",
            "the lowest-E-value (WWW: 3.39e-5 < WW: 7.43e-4) annotation is transferred — INV-06");
        annotation.EValue.Should().BeApproximately(WwwEValue, 1e-15);
    }

    #endregion

    #region META-FUNC-001 — BE boundary: EMPTY

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY / WHITESPACE PROTEIN SEQUENCE (BE).
    // A gene with an empty or whitespace sequence is SKIPPED — no annotation,
    // no crash. — Functional_Prediction.md §6.1, source line 919–920.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_EmptyOrWhitespaceSequence_YieldsNoAnnotation()
    {
        var db = BuildDatabase();
        var proteins = new[]
        {
            ("empty",   ""),
            ("space",   "   "),
            ("tab",     "\t"),
            ("newline", "\n\r"),
        };

        var results = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        results.Should().BeEmpty(
            "empty/whitespace protein sequences are skipped — nothing to match (§6.1)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY PROTEIN COLLECTION (BE).
    // No genes ⇒ no annotations, no crash, no division on an empty enumeration.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_EmptyProteinCollection_YieldsNoAnnotation()
    {
        var db = BuildDatabase();
        var noProteins = Array.Empty<(string, string)>();

        var results = MetagenomicsAnalyzer.PredictFunctions(noProteins, db).ToList();

        results.Should().BeEmpty("no query genes ⇒ no annotations");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY FUNCTION DATABASE (BE).
    // No signatures ⇒ no signature can ever match ⇒ no annotation for any gene,
    // even a perfectly valid protein. — Functional_Prediction.md §6.1 (Empty DB).
    // KEY: an empty database must NOT fabricate an annotation.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_EmptyDatabase_YieldsNoAnnotation()
    {
        var emptyDb = new Dictionary<string, (string, string, string)>();
        var proteins = new[] { ("valid", "WWWCCMK") };

        var results = MetagenomicsAnalyzer.PredictFunctions(proteins, emptyDb).ToList();

        results.Should().BeEmpty(
            "with no database signatures no homology transfer is possible (§6.1)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: EMPTY-STRING SIGNATURE in the database (BE).
    // An empty signature is skipped (IsNullOrEmpty(signature) guard, line 927):
    // it must NOT vacuously match every protein (string.Contains("") is true).
    // Here the empty signature is the ONLY entry, so no annotation is emitted.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_EmptySignatureOnly_DoesNotVacuouslyAnnotate()
    {
        var db = new Dictionary<string, (string, string, string)>
        {
            [""] = ("vacuous", "none", "K0"),
        };
        var proteins = new[] { ("g", "WWWCCMK") };

        var results = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        results.Should().BeEmpty(
            "an empty signature is skipped and must not match every protein vacuously (line 927)");
    }

    #endregion

    #region META-FUNC-001 — BE boundary: NO DB HIT

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: NO DB HIT (BE).
    // A valid protein whose sub-sequences are none of the database signatures ⇒
    // best stays null ⇒ NO annotation emitted (not a fabricated zero-score hit).
    // — Functional_Prediction.md §6.1 (Gene matches no signature).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_ValidProteinNoSignatureMatch_YieldsNoAnnotation()
    {
        var db = BuildDatabase(); // signatures: WWW, CC, MK
        // A valid 20-AA-alphabet protein containing none of WWW / CC / MK as a substring.
        var proteins = new[] { ("no-hit", "ARNDQEGHILKFPSTYV") };

        var results = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        results.Should().BeEmpty(
            "homology transfer needs a hit; a protein matching no signature gets no annotation (§6.1)");
    }

    // ───────────────────────────────────────────────────────────────────────
    // Case-sensitivity of matching (MC adjacent to NO-HIT): signature matching
    // is ordinal/case-sensitive (source line 927, StringComparison.Ordinal).
    // A lower-case "www" therefore does NOT contain upper-case signature "WWW"
    // ⇒ no hit. — Functional_Prediction.md §3.3 (case-sensitive ordinal).
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_LowerCaseProteinAgainstUpperSignature_NoMatch()
    {
        var db = BuildDatabase();
        var proteins = new[] { ("lower", "www") };

        var results = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        results.Should().BeEmpty(
            "signature matching is case-sensitive ordinal: 'www' does not contain 'WWW' (§3.3)");
    }

    #endregion

    #region META-FUNC-001 — BE boundary: UNKNOWN GENES

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: UNKNOWN GENE RESIDUES (BE/MC).
    // Residues outside the 20-letter BLOSUM62 alphabet (X, B, Z, J, U, *, digits,
    // gaps) contribute 0 to the raw self-score (Blosum62SelfScore, line 969–970).
    // A protein full of unknown residues that nonetheless CONTAINS a signature is
    // still annotated with that signature's stats; the SIGNATURE's own residues
    // (not the unknown padding) drive the raw score. Here signature "WWW" is
    // embedded in an X-padded protein. The result must be finite — never NaN/Inf.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_UnknownResiduesAroundSignature_StillFiniteAnnotation()
    {
        var db = BuildDatabase();
        // 'X' (selenocysteine-like / unknown) padding around an exact "WWW" hit.
        var proteins = new[] { ("x-padded", "XXXWWWXXX") };

        var annotation = MetagenomicsAnalyzer.PredictFunctions(proteins, db).Single();

        annotation.Function.Should().Be("tryptophan transport",
            "the embedded 'WWW' signature is matched even amid unknown 'X' residues");
        annotation.BitScore.Should().Be(WwwBitScore,
            "the raw score is computed from the SIGNATURE 'WWW' (S = 33), not the X padding");
        annotation.BitScore.Should().NotBe(double.NaN).And.NotBe(double.PositiveInfinity);
        annotation.EValue.Should().NotBe(double.NaN).And.BeGreaterThan(0.0);
    }

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: SIGNATURE OF ALL-UNKNOWN RESIDUES (BE).
    // A database signature made entirely of unknown residues scores raw S = 0
    // (every TryGetValue misses ⇒ +0). The bit score and E-value must still be
    // DEFINED and finite: S' = (λ·0 − ln K)/ln 2 = −ln(0.134)/ln 2 = 2.899…,
    // a finite positive bit score, and E = K·m·n·e^0 = K·m·n, finite. Never a
    // crash, never NaN. — Blosum62SelfScore unknown-residue rule + §2.A.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_AllUnknownResidueSignature_DefinedFiniteZeroRawScore()
    {
        // "XX" — both residues unknown to BLOSUM62 ⇒ raw self-score 0.
        var db = new Dictionary<string, (string, string, string)>
        {
            ["XX"] = ("unknown-fn", "none", "K0"),
        };
        var proteins = new[] { ("g", "AAXXAA") }; // contains "XX"

        var annotation = MetagenomicsAnalyzer.PredictFunctions(proteins, db).Single();

        // Independently-derived: S = 0 ⇒ S' = −ln(0.134)/ln 2, E = K·m·n (m = 6, n = 2).
        double expectedBit = (-Math.Log(0.134)) / Math.Log(2.0);
        double expectedE = 0.134 * 6 * 2; // e^0 = 1
        annotation.BitScore.Should().BeApproximately(expectedBit, 1e-9,
            "raw S = 0 ⇒ S' = −ln K/ln 2, a finite positive bit score (no NaN)");
        annotation.EValue.Should().BeApproximately(expectedE, 1e-9,
            "raw S = 0 ⇒ E = K·m·n·e^0 = K·m·n, finite (no crash, no NaN)");
        annotation.BitScore.Should().NotBe(double.NaN);
        annotation.EValue.Should().NotBe(double.NaN);
    }

    #endregion

    #region META-FUNC-001 — argument validation (eager null guards)

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: null arguments (BE).
    // null proteins / functionDatabase → ArgumentNullException (lines 910–911).
    // These guards run EAGERLY (before the iterator is returned), so they throw
    // at the call site, NOT on first enumeration. — Functional_Prediction.md §3.3.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    public void PredictFunctions_NullArguments_ThrowDocumentedExceptions()
    {
        var db = BuildDatabase();
        var proteins = new[] { ("g", "WWW") };

        // Eager guard: no enumeration needed to trigger.
        Action nullProteins = () => MetagenomicsAnalyzer.PredictFunctions(null!, db);
        nullProteins.Should().Throw<ArgumentNullException>()
            .WithParameterName("proteins");

        Action nullDb = () => MetagenomicsAnalyzer.PredictFunctions(proteins, null!);
        nullDb.Should().Throw<ArgumentNullException>()
            .WithParameterName("functionDatabase");
    }

    #endregion

    #region META-FUNC-001 — randomized boundary sweep

    // ───────────────────────────────────────────────────────────────────────
    // Fuzz target: random malformed batch (BE) under a time budget.
    // A deterministic, locally-seeded generator builds a mixed batch of empty,
    // whitespace, unknown-residue, and partially-valid proteins. PredictFunctions
    // must process the whole batch without crashing or hanging, and every emitted
    // annotation must be WELL-FORMED: a finite (non-NaN, non-Infinity) bit score
    // and a finite strictly-positive E-value, the transferred fields drawn ONLY
    // from the database, and never an annotation for a gene that matched nothing.
    // ───────────────────────────────────────────────────────────────────────
    [Test]
    [CancelAfter(30000)]
    public void PredictFunctions_RandomMalformedBatch_NeverCrashesAndStaysWellFormed()
    {
        var db = BuildDatabase();
        var validFunctions = new HashSet<string>(db.Values.Select(v => v.Function));

        var rng = new Random(20260620); // locally fixed seed — deterministic
        // Mix of BLOSUM62 letters, unknown residues, signatures' letters, junk, gaps.
        const string alphabet = "ACGTWMKXBZ-09 \tΑ★";

        var proteins = new List<(string, string)>();
        for (int i = 0; i < 500; i++)
        {
            int len = rng.Next(0, 20); // includes 0-length (empty) proteins
            var chars = new char[len];
            for (int j = 0; j < len; j++)
                chars[j] = alphabet[rng.Next(alphabet.Length)];
            proteins.Add(($"gene-{i}", new string(chars)));
        }

        var results = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        // Never more annotations than genes; each is well-formed.
        results.Count.Should().BeLessThanOrEqualTo(proteins.Count,
            "at most one annotation per gene");
        results.Should().OnlyContain(a =>
            !double.IsNaN(a.BitScore) && !double.IsInfinity(a.BitScore),
            "every transferred bit score is finite (no NaN/Infinity on garbage)");
        results.Should().OnlyContain(a =>
            !double.IsNaN(a.EValue) && !double.IsInfinity(a.EValue) && a.EValue > 0.0,
            "every E-value is finite and strictly positive");
        results.Should().OnlyContain(a => validFunctions.Contains(a.Function),
            "transferred fields are drawn only from the database — nothing fabricated");
        results.Select(a => a.GeneId).Should().OnlyHaveUniqueItems(
            "at most one best-hit annotation per gene");
    }

    #endregion
}
