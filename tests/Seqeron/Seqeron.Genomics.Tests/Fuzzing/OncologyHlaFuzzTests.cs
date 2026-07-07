using static Seqeron.Genomics.Oncology.OncologyAnalyzer;

namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the Oncology HLA-nomenclature parsing + allele-specific HLA LOH unit — ONCO-HLA-001.
/// The units under test are
/// <see cref="OncologyAnalyzer.ParseHlaAllele(string)"/> (WHO-nomenclature grammar validator),
/// <see cref="OncologyAnalyzer.TryParseHlaAllele(string?, out OncologyAnalyzer.HlaAllele)"/> (non-throwing wrapper), and
/// <see cref="OncologyAnalyzer.DetectHlaLoh(OncologyAnalyzer.HlaAlleleCopyNumber)"/> (LOHHLA allele-specific LOH rule),
/// all in src/Seqeron/Algorithms/Seqeron.Genomics.Oncology/OncologyAnalyzer.cs.
///
/// This file is scoped STRICTLY to ONCO-HLA-001 (HLA allele nomenclature + allele-specific LOH classification).
/// It does NOT touch the MHC peptide–binding classifier (ONCO-MHC-001, row 110,
/// <see cref="OncologyAnalyzer.ClassifyMhcBinding"/>) nor the neoantigen peptide-window generator
/// (ONCO-NEO-001, row 109).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate / boundary / malformed inputs to a unit and asserts the code never fails in an
/// undisciplined way: no KeyNotFoundException / NullReferenceException / unhandled FormatException leaking from
/// a malformed allele string, no IndexOutOfRange when a field block is empty / truncated, no DivideByZero when
/// the locus has zero coverage (all copy numbers 0), and no FALSE allele-specific LOH on a homozygous locus.
/// Every input must resolve to EITHER a well-defined, theory-correct result OR a *documented, intentional*
/// exception (ArgumentNullException / ArgumentException / FormatException for parsing; ArgumentException for an
/// out-of-domain copy number / p value).
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ONCO-HLA-001 — HLA nomenclature parsing + allele-specific HLA LOH (Oncology)
/// Checklist: docs/checklists/03_FUZZING.md, row 117.
/// Fuzz strategy exercised for THIS unit:
///   • MC = Malformed Content (невалідний контент).
///     Targets (checklist row 117): "ambiguous allele, homozygous locus, no coverage".
/// — docs/checklists/03_FUZZING.md §Description (strategy codes).
/// Target mapping for this specification-/decision-rule-driven (no trained model) unit:
///   • "ambiguous allele" = an HLA allele string that is ambiguous / not fully resolved / malformed / not a
///       member of the nomenclature grammar (missing HLA- prefix, missing '*', non-numeric / empty field,
///       too few / too many fields, illegal suffix letter, ambiguity codes like P/G group designators).
///       ParseHlaAllele must throw a DOCUMENTED FormatException/ArgumentException — never a KeyNotFound /
///       NullReference / IndexOutOfRange — and TryParseHlaAllele must guard it to a quiet `false`.
///   • "homozygous locus" = the two homologous alleles at the locus are identical (same name) → there is no
///       heterozygosity to lose, so allele-specific LOH is undefined / not callable (doc §6.2). The locus must
///       NOT yield a FALSE allele-specific LOH call, and identical alleles must NOT trigger a DivideByZero
///       (the rule does NO division by an allele-difference term).
///   • "no coverage" = zero reads over the HLA locus → both estimated allele copy numbers are 0 → the result
///       must be a documented, guarded classification (DetectHlaLoh never divides by depth, so 0/0 cannot
///       arise; per ASM-01 two CN < 0.5 alleles are reported as homozygous loss, NOT allele-specific LOH).
///
/// ───────────────────────────────────────────────────────────────────────────
/// The documented contract under test
/// ───────────────────────────────────────────────────────────────────────────
/// HLA_Nomenclature_And_Allele_Specific_LOH.md (docs/algorithms/Oncology/HLA_Nomenclature_And_Allele_Specific_LOH.md):
///   • Grammar (§2.2): name = "HLA-" + gene + '*' + 2–4 colon-separated numeric fields + optional N/L/S/C/A/Q suffix.
///   • INV-01 (§2.4): a parsed allele has 2 ≤ field count ≤ 4 (validator rejects <2 or >4 fields).
///   • INV-02 (§2.4): the normalized `Name` round-trips the parsed gene/fields/suffix.
///   • INV-03 (§2.4): HLA LOH called ⇔ (exactly one allele CN &lt; 0.5) ∧ (imbalance p &lt; 0.01).
///   • INV-04 (§2.4): when LOH is called, the lost allele is the one with CN &lt; 0.5.
///   • INV-05 (§2.4): CN = 0.5 is retained, p = 0.01 is NOT significant (both comparisons strict).
///   • ASM-01 (§2.3): both alleles CN &lt; 0.5 → reported as homozygous loss (HlaLostAllele.Both), NOT
///       allele-specific LOH.
///   • §3.3: ParseHlaAllele throws ArgumentNullException(null) / ArgumentException(empty-whitespace) /
///       FormatException(grammar). TryParseHlaAllele returns false for those. DetectHlaLoh throws
///       ArgumentException for negative CN or p ∉ [0, 1].
///   • §6.1 edge cases: "HLA-A*02" → FormatException; "HLA-A*02:01:01:01:01" → FormatException;
///       "HLA-A*02:01X" → FormatException; "A*02:01" → FormatException; CN exactly 0.5 → not lost;
///       p exactly 0.01 → not significant; low CN but p ≥ 0.01 → no LOH.
///   • §6.2: a homozygous locus (two identical alleles) cannot be assessed for allele-specific loss.
///   • §7.1 worked example: ParseHlaAllele("HLA-A*24:02:01:02L") ⇒ Gene A, Fields [24,02,01,02], Suffix Low;
///       DetectHlaLoh(("HLA-A*02:01",1.8,"HLA-A*11:01",0.30,0.001)) ⇒ IsLoh, LostAllele Allele2.
///
/// All randomness is LOCALLY seeded (new Random(seed)); no shared static Rng.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public sealed class OncologyHlaFuzzTests
{
    // Documented LOHHLA thresholds (mirrored from the doc §2.2/§4.2 and the public consts), used to pin the
    // business contract in hand-built examples independently of the production constants.
    private const double LossCnThreshold = 0.5;   // strict `<` ⇒ allele lost
    private const double ImbalanceP = 0.01;       // strict `<` ⇒ allelic imbalance significant

    // ── Well-formed-result assertion helpers ─────────────────────────────────
    // A parsed HlaAllele is "well formed" iff its field count is in [2, 4], every field is a non-empty digit
    // group, the suffix is a defined enum member, and its normalized Name round-trips (INV-01, INV-02). This
    // guard stops a fuzz test rubber-stamping an out-of-contract parse.
    private static void AssertWellFormed(HlaAllele allele)
    {
        allele.Fields.Should().NotBeNull();
        allele.Fields.Count.Should().BeInRange(2, 4, "INV-01: 2–4 nomenclature fields");
        allele.Fields.Should().OnlyContain(f => f.Length > 0 && f.All(char.IsDigit),
            "each field is a non-empty digit group");
        allele.Suffix.Should().BeOneOf(
            HlaExpressionSuffix.None, HlaExpressionSuffix.Null, HlaExpressionSuffix.Low,
            HlaExpressionSuffix.Secreted, HlaExpressionSuffix.Cytoplasm, HlaExpressionSuffix.Aberrant,
            HlaExpressionSuffix.Questionable);
        allele.Gene.Should().NotBeNullOrEmpty();
        // INV-02: Name round-trips through a re-parse to the same components.
        HlaAllele reparsed = ParseHlaAllele(allele.Name);
        reparsed.Gene.Should().Be(allele.Gene);
        reparsed.Fields.Should().Equal(allele.Fields);
        reparsed.Suffix.Should().Be(allele.Suffix);
    }

    // A HlaLohResult is "well formed" iff LostAllele is a defined enum member and the IsLoh/LostAllele pairing is
    // self-consistent (LOH ⇒ exactly one of Allele1/Allele2; no LOH ⇒ None or Both).
    private static void AssertWellFormed(HlaLohResult r)
    {
        r.LostAllele.Should().BeOneOf(
            HlaLostAllele.None, HlaLostAllele.Allele1, HlaLostAllele.Allele2, HlaLostAllele.Both);
        if (r.IsLoh)
        {
            r.LostAllele.Should().BeOneOf(HlaLostAllele.Allele1, HlaLostAllele.Allele2);
            r.AllelicImbalanceSignificant.Should().BeTrue("LOH requires significant imbalance (INV-03)");
        }
        else
        {
            r.LostAllele.Should().BeOneOf(HlaLostAllele.None, HlaLostAllele.Both);
        }
    }

    // Independent reference oracle for the LOHHLA decision rule (§2.2/§2.4), built WITHOUT reusing the
    // implementation branch structure.
    private static (bool isLoh, HlaLostAllele lost) ExpectedLoh(double cn1, double cn2, double p)
    {
        bool sig = p < ImbalanceP;
        bool a1 = cn1 < LossCnThreshold;
        bool a2 = cn2 < LossCnThreshold;
        if (a1 && a2) return (false, HlaLostAllele.Both);
        if (sig && a1) return (true, HlaLostAllele.Allele1);
        if (sig && a2) return (true, HlaLostAllele.Allele2);
        return (false, HlaLostAllele.None);
    }

    #region ONCO-HLA-001 — Positive sanity (documented grammar, thresholds, worked examples)

    // §7.1 worked example: a 4-field name with a Low-expression suffix parses to all documented components.
    [Test]
    public void ParseHlaAllele_DocWorkedExample_FourFieldWithSuffix()
    {
        HlaAllele allele = ParseHlaAllele("HLA-A*24:02:01:02L");
        allele.Gene.Should().Be("A");
        allele.Fields.Should().Equal("24", "02", "01", "02");
        allele.Suffix.Should().Be(HlaExpressionSuffix.Low);
        allele.Name.Should().Be("HLA-A*24:02:01:02L");
        AssertWellFormed(allele);
    }

    // §7.1 worked example: clear allelic imbalance on a HETEROZYGOUS locus (one allele CN 0.30 < 0.5, p 0.001
    // significant) is called allele-specific HLA LOH, with allele 2 named as the lost homolog (INV-04).
    [Test]
    public void DetectHlaLoh_DocWorkedExample_ImbalancedHeterozygousLocus_IsLohAllele2()
    {
        HlaLohResult r = DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 1.8, "HLA-A*11:01", 0.30, 0.001));
        r.IsLoh.Should().BeTrue();
        r.LostAllele.Should().Be(HlaLostAllele.Allele2);
        r.AllelicImbalanceSignificant.Should().BeTrue();
        AssertWellFormed(r);
    }

    // POSITIVE contract: a BALANCED heterozygous locus (both alleles retained, CN ≈ 1) is NOT HLA-LOH even with
    // a tiny p — the loss-threshold guard (no allele CN < 0.5) prevents a false call (INV-03).
    [Test]
    public void DetectHlaLoh_BalancedHeterozygousLocus_IsNotLoh()
    {
        HlaLohResult r = DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-B*07:02", 1.05, "HLA-B*08:01", 0.98, 1e-6));
        r.IsLoh.Should().BeFalse("both homologs retained (CN ≥ 0.5) ⇒ no allele-specific loss");
        r.LostAllele.Should().Be(HlaLostAllele.None);
        AssertWellFormed(r);
    }

    // POSITIVE contract: a clearly imbalanced locus where ALLELE 1 is the depleted homolog is called LOH on
    // allele 1 (the other LOH branch, INV-04).
    [Test]
    public void DetectHlaLoh_Allele1Depleted_IsLohAllele1()
    {
        HlaLohResult r = DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-C*07:01", 0.12, "HLA-C*07:02", 1.7, 0.0005));
        r.IsLoh.Should().BeTrue();
        r.LostAllele.Should().Be(HlaLostAllele.Allele1);
        AssertWellFormed(r);
    }

    // §6.1 / INV-05: the loss threshold is strict `< 0.5` — CN exactly 0.5 is RETAINED, so no LOH.
    [Test]
    public void DetectHlaLoh_CopyNumberExactlyHalf_IsRetained()
    {
        HlaLohResult r = DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 1.5, "HLA-A*03:01", 0.5, 0.0001));
        r.IsLoh.Should().BeFalse("CN exactly 0.5 is not < 0.5 (strict, INV-05)");
        r.LostAllele.Should().Be(HlaLostAllele.None);
    }

    // §6.1 / INV-05: the imbalance threshold is strict `< 0.01` — p exactly 0.01 is NOT significant ⇒ no LOH,
    // even though one allele is below the loss threshold (over-calling guard).
    [Test]
    public void DetectHlaLoh_PValueExactlyThreshold_IsNotSignificantNoLoh()
    {
        HlaLohResult r = DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-A*02:01", 0.20, "HLA-A*11:01", 1.6, 0.01));
        r.IsLoh.Should().BeFalse("p exactly 0.01 is not < 0.01 (strict, INV-05)");
        r.AllelicImbalanceSignificant.Should().BeFalse();
        // The allele is depleted but the imbalance is not significant ⇒ retained call (None), not a false LOH.
        r.LostAllele.Should().Be(HlaLostAllele.None);
    }

    // §6.1: low CN but p ≥ 0.01 ⇒ no LOH (the over-calling guard rejects insignificant imbalance).
    [Test]
    public void DetectHlaLoh_LowCnButInsignificantP_NoLoh()
    {
        HlaLohResult r = DetectHlaLoh(
            new HlaAlleleCopyNumber("HLA-B*15:01", 0.10, "HLA-B*40:01", 2.0, 0.5));
        r.IsLoh.Should().BeFalse();
        r.LostAllele.Should().Be(HlaLostAllele.None);
    }

    #endregion

    #region ONCO-HLA-001 — MC: ambiguous / malformed allele (documented throw, never KeyNotFound / crash)

    // §6.1 / §3.3: each malformed or ambiguous allele string is rejected with a DOCUMENTED FormatException —
    // NEVER a KeyNotFoundException, NullReferenceException, IndexOutOfRangeException or unhandled crash.
    [TestCase("A*02:01")]                 // §6.1: missing 'HLA-' prefix
    [TestCase("HLA-A02:01")]              // missing '*' gene/field separator
    [TestCase("HLA-A*02")]               // §6.1: single field (below 2-field minimum, INV-01)
    [TestCase("HLA-A*02:01:01:01:01")]   // §6.1: five fields (above 4-field maximum, INV-01)
    [TestCase("HLA-A*02:01X")]           // §6.1: X ∉ {N,L,S,C,A,Q} illegal suffix
    [TestCase("HLA-A*02:0Z")]            // non-numeric field (Z is not a suffix letter)
    [TestCase("HLA-A*\uFF10\uFF12:\uFF10\uFF11")]          // non-ASCII (fullwidth) digits are not WHO nomenclature
    [TestCase("HLA-A*02:")]              // trailing empty field
    [TestCase("HLA-A*:01")]              // leading empty field
    [TestCase("HLA-A*02::01")]           // empty interior field
    [TestCase("HLA-A*")]                 // empty field block after '*'
    [TestCase("HLA-*02:01")]             // empty gene
    [TestCase("HLA-A*02:01P")]           // P-group ambiguity designator (not a single allele)
    [TestCase("HLA-A*02:01G")]           // G-group ambiguity designator (not a single allele)
    [TestCase("HLA-A*02:01:01:01N:02")]  // suffix mid-name (only trailing letter allowed)
    [TestCase("HLA-DRB1*15:01:01:01:01")]// five fields on a class II gene
    [TestCase("not an allele at all")]   // free-text garbage
    public void ParseHlaAllele_AmbiguousOrMalformed_ThrowsDocumentedException(string name)
    {
        Action act = () => ParseHlaAllele(name);
        act.Should().Throw<Exception>()
            .Which.Should().Match(ex => ex is FormatException || ex is ArgumentException,
                "malformed/ambiguous allele names are rejected with a documented Format/Argument exception, "
                + "never KeyNotFound / NullReference / IndexOutOfRange");
    }

    // §3.3: null is ArgumentNullException; empty / whitespace is ArgumentException — distinct documented guards.
    [Test]
    public void ParseHlaAllele_NullAndBlank_ThrowDocumentedArgumentExceptions()
    {
        Action nullAct = () => ParseHlaAllele(null!);
        nullAct.Should().Throw<ArgumentNullException>();

        foreach (string blank in new[] { "", "   ", "\t", "\n" })
        {
            Action blankAct = () => ParseHlaAllele(blank);
            blankAct.Should().Throw<ArgumentException>($"'{blank}' is empty/whitespace");
        }
    }

    // §3.3: TryParseHlaAllele is the documented non-throwing guard — every malformed / ambiguous / null input
    // resolves to a quiet `false` with `default` out, NEVER an escaping exception.
    [TestCase(null)]
    [TestCase("")]
    [TestCase("   ")]
    [TestCase("A*02:01")]
    [TestCase("HLA-A*02")]
    [TestCase("HLA-A*02:01X")]
    [TestCase("HLA-A*02:0Z")]
    [TestCase("HLA-A*02:01P")]
    [TestCase("HLA-A*02:01:01:01:01")]
    [TestCase("garbage")]
    public void TryParseHlaAllele_MalformedOrNull_ReturnsFalseNoThrow(string? name)
    {
        bool ok = true;
        HlaAllele allele = default;
        Action act = () => ok = TryParseHlaAllele(name, out allele);
        act.Should().NotThrow("the Try wrapper must swallow format/argument errors");
        ok.Should().BeFalse();
        allele.Should().Be(default(HlaAllele));
    }

    // INJ-adjacent MC: control / null bytes / unicode / very long junk inside the name must not crash the
    // validator — it must return a clean documented reject (false), never a KeyNotFound / overflow.
    [Test]
    public void TryParseHlaAllele_ControlBytesUnicodeAndLongJunk_RejectNoThrow()
    {
        // Each entry is structurally malformed (null byte inside a field, fullwidth non-ASCII digits, a
        // replacement-char field, no '*' separator, an all-colon empty-field run). Trimmable trailing
        // whitespace and an unvalidated gene token are intentionally NOT here — those are documented-accepted.
        string[] nasties =
        {
            "HLA-A*02\u0000:01", "HLA-A*\uFF10\uFF12:\uFF10\uFF11", "HLA-A*02:01\uFFFD",
            "HLA-" + new string('9', 5000), "HLA-A*" + new string(':', 100),
        };
        foreach (string nasty in nasties)
        {
            bool ok = true;
            Action act = () => ok = TryParseHlaAllele(nasty, out _);
            act.Should().NotThrow($"input '{nasty[..Math.Min(20, nasty.Length)]}…' must be guarded");
            ok.Should().BeFalse();
        }
    }

    // Random malformed-content fuzz: random gene/field/suffix permutations, many invalid, must each either
    // parse to a well-formed allele OR reject cleanly via the Try wrapper — never throw.
    [Test]
    [CancelAfter(20_000)]
    public void TryParseHlaAllele_RandomMalformedContent_NeverThrows_WellFormedWhenAccepted()
    {
        var rng = new Random(117_0001);
        string[] geneTokens = { "A", "B", "C", "DRB1", "DQB1", "", "1", "a", "ZZ", "Ω" };
        string[] fieldTokens = { "02", "01", "001", "0", "", "A", "1a", ":", "999", "07" };
        char[] suffixTokens = { '\0', 'N', 'L', 'S', 'C', 'A', 'Q', 'X', 'Z', '9' };

        for (int i = 0; i < 20_000; i++)
        {
            string gene = geneTokens[rng.Next(geneTokens.Length)];
            int nFields = rng.Next(0, 6);
            string fieldBlock = string.Join(":",
                Enumerable.Range(0, nFields).Select(_ => fieldTokens[rng.Next(fieldTokens.Length)]));
            char suffix = suffixTokens[rng.Next(suffixTokens.Length)];
            bool withPrefix = rng.Next(4) != 0; // mostly with prefix
            bool withStar = rng.Next(4) != 0;

            string name =
                (withPrefix ? "HLA-" : "")
                + gene
                + (withStar ? "*" : "")
                + fieldBlock
                + (suffix == '\0' ? "" : suffix.ToString());

            bool ok = false;
            HlaAllele allele = default;
            Action act = () => ok = TryParseHlaAllele(name, out allele);
            act.Should().NotThrow($"random content '{name}' must be guarded by the Try wrapper");
            if (ok)
            {
                AssertWellFormed(allele);
            }
        }
    }

    #endregion

    #region ONCO-HLA-001 — MC: homozygous locus (no false LOH, no DivideByZero on identical alleles)

    // §6.2: a HOMOZYGOUS locus has two identical alleles ⇒ no heterozygosity to lose. Modeled by identical
    // allele names with equal (balanced) copy numbers: the rule must NOT call allele-specific LOH and must NOT
    // DivideByZero on the identical homologs (the rule does NO division by an allele-difference term).
    [Test]
    public void DetectHlaLoh_HomozygousLocus_BalancedIdenticalAlleles_NoFalseLoh()
    {
        var input = new HlaAlleleCopyNumber("HLA-A*02:01", 1.0, "HLA-A*02:01", 1.0, 0.0001);
        HlaLohResult r = default;
        Action act = () => r = DetectHlaLoh(input);
        act.Should().NotThrow("identical homologs must not DivideByZero");
        r.IsLoh.Should().BeFalse("a homozygous locus cannot have allele-specific LOH (doc §6.2)");
        r.LostAllele.Should().Be(HlaLostAllele.None);
        AssertWellFormed(r);
    }

    // A homozygous locus that is entirely deleted (both identical homologs CN 0, e.g. homozygous deletion) must
    // be reported as homozygous loss (Both), NOT a false allele-specific LOH call (ASM-01) — and must not crash.
    [Test]
    public void DetectHlaLoh_HomozygousLocus_BothHomologsDeleted_IsHomozygousLossNotLoh()
    {
        var input = new HlaAlleleCopyNumber("HLA-B*07:02", 0.0, "HLA-B*07:02", 0.0, 0.0001);
        HlaLohResult r = DetectHlaLoh(input);
        r.IsLoh.Should().BeFalse("both homologs lost ⇒ homozygous loss, not allele-specific LOH (ASM-01)");
        r.LostAllele.Should().Be(HlaLostAllele.Both);
        AssertWellFormed(r);
    }

    // Even with a degenerate (depleted) allele-2 CN, identical allele NAMES still describe a homozygous locus.
    // The rule looks only at copy number (the names are informational, §3.2), so the classification is driven by
    // CN and must remain self-consistent and crash-free — never a name-keyed KeyNotFound on the equal names.
    [Test]
    public void DetectHlaLoh_HomozygousNames_AsymmetricCn_StaysSelfConsistentNoCrash()
    {
        var input = new HlaAlleleCopyNumber("HLA-A*02:01", 1.6, "HLA-A*02:01", 0.10, 0.0005);
        HlaLohResult r = default;
        Action act = () => r = DetectHlaLoh(input);
        act.Should().NotThrow();
        AssertWellFormed(r);
        // CN-driven contract: allele2 depleted + significant ⇒ matches the oracle (no name de-dup short-circuit).
        var (isLoh, lost) = ExpectedLoh(1.6, 0.10, 0.0005);
        r.IsLoh.Should().Be(isLoh);
        r.LostAllele.Should().Be(lost);
    }

    #endregion

    #region ONCO-HLA-001 — MC: no coverage (zero depth ⇒ no DivideByZero, documented guarded result)

    // "No coverage" = zero reads over the locus ⇒ both estimated copy numbers are 0. The rule must NOT
    // DivideByZero (it performs no division by depth) and must return a documented, guarded classification:
    // both CN < 0.5 ⇒ homozygous loss (Both), NOT allele-specific LOH.
    [Test]
    public void DetectHlaLoh_NoCoverage_ZeroCopyNumbers_HomozygousLossNoDivideByZero()
    {
        var input = new HlaAlleleCopyNumber("HLA-A*02:01", 0.0, "HLA-A*11:01", 0.0, 1.0);
        HlaLohResult r = default;
        Action act = () => r = DetectHlaLoh(input);
        act.Should().NotThrow("zero coverage must not DivideByZero");
        r.IsLoh.Should().BeFalse();
        r.LostAllele.Should().Be(HlaLostAllele.Both, "both CN 0 < 0.5 ⇒ homozygous loss (ASM-01)");
        AssertWellFormed(r);
    }

    // No coverage with a degenerate p value (p = 1, fully insignificant) — still no LOH, no crash; the
    // imbalance-significance flag must reflect the strict `< 0.01` rule (false here).
    [Test]
    public void DetectHlaLoh_NoCoverage_InsignificantP_NoLohFlagFalse()
    {
        HlaLohResult r = DetectHlaLoh(new HlaAlleleCopyNumber("HLA-C*07:02", 0.0, "HLA-C*06:02", 0.0, 1.0));
        r.AllelicImbalanceSignificant.Should().BeFalse();
        r.IsLoh.Should().BeFalse();
        r.LostAllele.Should().Be(HlaLostAllele.Both);
    }

    // §3.3: an out-of-domain copy number (negative, the "no-data sentinel" abuse) or an out-of-[0,1] p value is
    // a documented ArgumentException — never a silent NaN-leaking classification.
    [Test]
    public void DetectHlaLoh_OutOfDomainInputs_ThrowDocumentedArgumentException()
    {
        Action negCn = () => DetectHlaLoh(new HlaAlleleCopyNumber("HLA-A*02:01", -1.0, "HLA-A*11:01", 1.0, 0.001));
        negCn.Should().Throw<ArgumentException>("a negative copy number is out of domain");

        foreach (double badP in new[] { -0.0001, -1.0, 1.0001, 2.0, double.NaN, double.PositiveInfinity })
        {
            Action badPAct = () => DetectHlaLoh(new HlaAlleleCopyNumber("HLA-A*02:01", 1.0, "HLA-A*11:01", 0.2, badP));
            badPAct.Should().Throw<ArgumentException>($"p {badP} is outside [0, 1]");
        }
    }

    #endregion

    #region ONCO-HLA-001 — Invariants over broad random fuzz (oracle cross-check, no false LOH)

    // INV-03 / INV-04 / ASM-01: over many random valid (CN ≥ 0, p ∈ [0,1]) inputs — including zero-coverage,
    // homozygous, and balanced loci — DetectHlaLoh never throws, always returns a well-formed result, and
    // exactly matches an independent decision-rule oracle. This is the anti-rubber-stamp check: a false LOH on a
    // retained/homozygous locus would diverge from the oracle and fail.
    [Test]
    [CancelAfter(20_000)]
    public void DetectHlaLoh_RandomValidInputs_MatchOracle_NoFalseLoh()
    {
        var rng = new Random(117_0002);
        for (int i = 0; i < 50_000; i++)
        {
            // Mix in many degenerate edges: zero coverage, exactly-0.5 CN, exactly-0.01 p, identical CNs.
            double cn1 = PickCn(rng);
            double cn2 = PickCn(rng);
            double p = PickP(rng);

            var input = new HlaAlleleCopyNumber("HLA-A*02:01", cn1, "HLA-A*11:01", cn2, p);
            HlaLohResult r = default;
            Action act = () => r = DetectHlaLoh(input);
            act.Should().NotThrow($"cn1={cn1} cn2={cn2} p={p} is a valid domain input");
            AssertWellFormed(r);

            var (isLoh, lost) = ExpectedLoh(cn1, cn2, p);
            r.IsLoh.Should().Be(isLoh, $"cn1={cn1} cn2={cn2} p={p}");
            r.LostAllele.Should().Be(lost, $"cn1={cn1} cn2={cn2} p={p}");
            r.AllelicImbalanceSignificant.Should().Be(p < ImbalanceP);
        }
    }

    private static double PickCn(Random rng) => rng.Next(6) switch
    {
        0 => 0.0,                              // no coverage
        1 => 0.5,                              // exact retained boundary
        2 => rng.NextDouble() * 0.49,          // below loss threshold
        3 => 0.5 + rng.NextDouble() * 2.5,     // retained
        _ => rng.NextDouble() * 3.0,           // broad
    };

    private static double PickP(Random rng) => rng.Next(5) switch
    {
        0 => 0.01,                             // exact significance boundary
        1 => rng.NextDouble() * 0.01,          // significant
        2 => 1.0,                              // fully insignificant
        _ => rng.NextDouble(),                 // broad [0,1)
    };

    #endregion
}
