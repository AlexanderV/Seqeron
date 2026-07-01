// PROTMOTIF-PATTERN-001 — Protein Pattern Matching Methods
// Evidence: docs/Evidence/PROTMOTIF-PATTERN-001-Evidence.md
// TestSpec: tests/TestSpecs/PROTMOTIF-PATTERN-001.md
// Source: ExPASy/PROSITE pattern syntax (scanprosite_doc.html); de Castro et al. (2006) NAR 34:W362–W365;
//         Schneider TD, Stephens RM (1990). Nucleic Acids Res 18(20):6097–6100.
using System;
using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests.Unit.Analysis;

/// <summary>
/// PROTMOTIF-PATTERN-001: canonical test file for the four pattern-matching methods of
/// <see cref="ProteinMotifFinder"/>: FindMotifByPattern, ConvertPrositeToRegex,
/// FindMotifByProsite, FindDomains. Expected values are derived from PROSITE worked
/// examples (PS00001/05/16/17/29) and the information-content definition of
/// Schneider &amp; Stephens (1990) — never from running the implementation.
/// </summary>
[TestFixture]
[Category("PROTMOTIF-PATTERN-001")]
public class ProteinMotifFinder_FindMotifByPattern_Tests
{
    // Information content per position = log2(20/k); 20 = standard amino-acid alphabet.
    // Schneider & Stephens (1990), Nucleic Acids Res 18(20):6097-6100.
    private static readonly double Log2Of20 = Math.Log2(20.0); // ≈ 4.321928094887363 (k=1, fixed residue)
    private static readonly double Log2Of10 = Math.Log2(10.0); // ≈ 3.321928094887362 (k=2, two-residue class)
    private const double Tol = 1e-10;

    #region ConvertPrositeToRegex — PA-line grammar (M1–M5, INV-04)

    // M1 — PS00016 R-G-D (Cell attachment) -> RGD. https://prosite.expasy.org/PS00016
    [Test]
    public void ConvertPrositeToRegex_PS00016Literal_ProducesRGD()
    {
        Assert.That(ConvertPrositeToRegex("R-G-D"), Is.EqualTo("RGD"),
            "Literal residues map directly and '-' separators are dropped (PROSITE syntax)");
    }

    // M2 — PS00001 N-{P}-[ST]-{P} -> N[^P][ST][^P]. https://prosite.expasy.org/PS00001
    [Test]
    public void ConvertPrositeToRegex_PS00001ExclusionAndClass_ProducesNegatedAndClass()
    {
        Assert.That(ConvertPrositeToRegex("N-{P}-[ST]-{P}"), Is.EqualTo("N[^P][ST][^P]"),
            "{P} maps to [^P] and [ST] is preserved (PROSITE PS00001)");
    }

    // M3 — PS00017 [AG]-x(4)-G-K-[ST] -> [AG].{4}GK[ST]. https://prosite.expasy.org/PS00017
    [Test]
    public void ConvertPrositeToRegex_PS00017FixedRange_ProducesQuantifier()
    {
        Assert.That(ConvertPrositeToRegex("[AG]-x(4)-G-K-[ST]"), Is.EqualTo("[AG].{4}GK[ST]"),
            "x(4) maps to .{4} (PROSITE PS00017 P-loop)");
    }

    // M4 — PS00029 L-x(6)-L-x(6)-L-x(6)-L -> L.{6}L.{6}L.{6}L. https://prosite.expasy.org/PS00029
    [Test]
    public void ConvertPrositeToRegex_PS00029RepeatedRanges_ProducesFullRegex()
    {
        Assert.That(ConvertPrositeToRegex("L-x(6)-L-x(6)-L-x(6)-L"), Is.EqualTo("L.{6}L.{6}L.{6}L"),
            "Repeated x(6) segments each map to .{6} (PROSITE PS00029 leucine zipper)");
    }

    // M5 — Anchors and trailing period (PROSITE syntax: '<'->^, '>'->$, '.' terminates).
    [Test]
    public void ConvertPrositeToRegex_AnchorsAndTrailingPeriod_ProduceExpectedRegex()
    {
        Assert.Multiple(() =>
        {
            Assert.That(ConvertPrositeToRegex("<M-x-K>"), Is.EqualTo("^M.K$"),
                "'<' maps to ^ and '>' maps to $ (PROSITE terminal anchors)");
            Assert.That(ConvertPrositeToRegex("R-G-D."), Is.EqualTo("RGD"),
                "A trailing period terminates the pattern (PROSITE PA line)");
            Assert.That(ConvertPrositeToRegex("R-G-D.A-B-C"), Is.EqualTo("RGD"),
                "Characters after the terminating period are ignored");
        });
    }

    // M11 — PS00006 [ST]-x(2)-[DE] -> [ST].{2}[DE]: the x(n) gap form. https://prosite.expasy.org/PS00006
    [Test]
    public void ConvertPrositeToRegex_PS00006GapCount_ProducesRangeQuantifier()
    {
        Assert.That(ConvertPrositeToRegex("[ST]-x(2)-[DE]"), Is.EqualTo("[ST].{2}[DE]"),
            "x(2) maps to .{2} (PROSITE PS00006 casein-kinase-II site)");
    }

    // M12 — PS00028 ranges x(2,4)/x(3,5) -> .{2,4}/.{3,5}: the x(n,m) range form.
    // https://prosite.expasy.org/PS00028
    [Test]
    public void ConvertPrositeToRegex_PS00028Ranges_ProduceRangeQuantifiers()
    {
        Assert.That(ConvertPrositeToRegex("C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H"),
            Is.EqualTo("C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H"),
            "x(2,4) -> .{2,4} and x(3,5) -> .{3,5} (PROSITE PS00028 C2H2 zinc finger)");
    }

    // M13 — PS00004 [RK](2)-x-[ST] -> [RK]{2}.[ST]: a fixed count on a residue element (class).
    // PROSITE doc: "A(3) corresponds to A-A-A". https://prosite.expasy.org/PS00004
    [Test]
    public void ConvertPrositeToRegex_PS00004FixedCountOnClass_ProducesQuantifier()
    {
        Assert.That(ConvertPrositeToRegex("[RK](2)-x-[ST]"), Is.EqualTo("[RK]{2}.[ST]"),
            "[RK](2) -> [RK]{2} — a fixed count applies to the preceding element (PROSITE PS00004)");
    }

    // M14 — Fixed count on a single residue letter: A(3) -> A{3}.
    // PROSITE doc (scanprosite_doc / prosuser §IV.E): "A(3) corresponds to A-A-A".
    [Test]
    public void ConvertPrositeToRegex_FixedCountOnLetter_ProducesQuantifier()
    {
        Assert.That(ConvertPrositeToRegex("A(3)"), Is.EqualTo("A{3}"),
            "A(3) maps to A{3} per PROSITE 'A(3) corresponds to A-A-A'");
    }

    #endregion

    #region FindMotifByPattern — positions and IC scoring (M6–M8, INV-02, INV-03)

    // M6 — Literal RGD at 0-based position 3 in AAARGDAAA (PS00016).
    [Test]
    public void FindMotifByPattern_LiteralRGD_ReturnsExactPositionAndSubstring()
    {
        var matches = FindMotifByPattern("AAARGDAAA", "RGD", "RGD", "PS00016").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.EqualTo(1), "Exactly one RGD occurrence exists");
            Assert.That(matches[0].Start, Is.EqualTo(3), "RGD starts at 0-based position 3");
            Assert.That(matches[0].End, Is.EqualTo(5), "End is inclusive: 3 + len(3) - 1 = 5 (INV-02)");
            Assert.That(matches[0].Sequence, Is.EqualTo("RGD"), "Matched substring is RGD (INV-02)");
        });
    }

    // M7 — IC score for RGD = 3 fixed residues = 3*log2(20). Schneider & Stephens (1990).
    [Test]
    public void FindMotifByPattern_RGDScore_EqualsThreeTimesLog2Of20()
    {
        var match = FindMotifByPattern("AAARGDAAA", "RGD").Single();

        Assert.That(match.Score, Is.EqualTo(3.0 * Log2Of20).Within(Tol),
            "Three fixed residues contribute 3*log2(20/1) bits (INV-03; Schneider & Stephens 1990)");
    }

    // M8 — IC score for [ST].[RK] = log2(10) + 0 + log2(10). Two 2-residue classes + one wildcard.
    [Test]
    public void FindMotifByPattern_ClassPatternScore_EqualsTwoTimesLog2Of10()
    {
        // PS00005 PKC site [ST]-x-[RK] -> [ST].[RK]; "SAR" matches.
        var match = FindMotifByPattern("AASARAA", "[ST].[RK]").Single();

        Assert.That(match.Score, Is.EqualTo(2.0 * Log2Of10).Within(Tol),
            "Two 2-residue classes (log2(20/2) each) + wildcard (0) = 2*log2(10) bits (INV-03)");
    }

    // M8b — IC score for a negated class: PS00001 N[^P][ST][^P].
    // k = 1, 19, 2, 19 -> IC = log2(20) + log2(20/19) + log2(10) + log2(20/19). Schneider & Stephens (1990).
    [Test]
    public void FindMotifByPattern_NegatedClassScore_EqualsSumOfPerPositionIC()
    {
        // "NASA": N (k=1), A (in [^P], k=19), S (in [ST], k=2), A (in [^P], k=19).
        var match = FindMotifByPattern("NASA", "N[^P][ST][^P]").Single();

        double expected = Math.Log2(20.0 / 1) + Math.Log2(20.0 / 19)
                        + Math.Log2(20.0 / 2) + Math.Log2(20.0 / 19);
        Assert.That(match.Score, Is.EqualTo(expected).Within(Tol),
            "Negated class [^P] allows 20-1=19 residues -> log2(20/19) bits each (INV-03)");
    }

    #endregion

    #region FindMotifByProsite — end-to-end delegate (M9, INV-05)

    // M9 — PS00001 over AANASAAANGTAAAA: matches at 0-based starts {2,8} = {NASA, NGTA}.
    [Test]
    public void FindMotifByProsite_PS00001_FindsBothSequonsAtExactPositions()
    {
        var matches = FindMotifByProsite("AANASAAANGTAAAA", "N-{P}-[ST]-{P}", "N-glycosylation").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches.Select(m => m.Start), Is.EqualTo(new[] { 2, 8 }),
                "PS00001 N-{P}-[ST]-{P} matches at 0-based positions 2 and 8");
            Assert.That(matches[0].Sequence, Is.EqualTo("NASA"), "First sequon is NASA");
            Assert.That(matches[1].Sequence, Is.EqualTo("NGTA"), "Second sequon is NGTA");
        });
    }

    #endregion

    #region FindDomains — built-in signature delegate (M10)

    // M10 — P-loop / kinase ATP-binding signature [AG]-x(4)-G-K-[ST] (PS00017-equivalent).
    [Test]
    public void FindDomains_PLoopSignature_DetectsKinaseAtpBindingAtExactPosition()
    {
        // Synthetic: 3 prefix residues, then A-AAAA-G-K-S = [AG].{4}GK[ST] starting at index 3.
        const string sequence = "QQQACDEFGKSQQQ";
        var domains = FindDomains(sequence).ToList();

        var kinase = domains.Where(d => d.Accession == "PF00069").ToList();
        Assert.Multiple(() =>
        {
            Assert.That(kinase, Has.Count.EqualTo(1), "Exactly one P-loop/kinase ATP-binding match");
            Assert.That(kinase[0].Start, Is.EqualTo(3),
                "Walker-A [AG]x(4)GK[ST] starts at 0-based position 3 (the A)");
            Assert.That(kinase[0].End, Is.EqualTo(10),
                "Match spans 8 residues (A..S): End = 3 + 8 - 1 = 10");
            Assert.That(kinase[0].Name, Is.EqualTo("Protein Kinase ATP-binding"),
                "Domain name is the kinase ATP-binding signature (ProteinDomain.Name for PF00069)");
        });
    }

    #endregion

    #region Overlapping enumeration, E-value, rejection, case (S1–S4, INV-06, INV-07)

    // S1 — lookahead enables overlapping starts: "A.A" over "AAAA" -> starts 0 and 1.
    [Test]
    public void FindMotifByPattern_OverlappingPattern_EnumeratesAllStarts()
    {
        var starts = FindMotifByPattern("AAAA", "A.A").Select(m => m.Start).ToList();

        Assert.That(starts, Is.EqualTo(new[] { 0, 1 }),
            "Lookahead wrapper lists overlapping start positions 0 and 1");
    }

    // S2 — E-value = (N - L + 1) * 2^(-Score). For RGD in AAARGDAAA: (9-3+1)*2^(-3*log2(20)).
    [Test]
    public void FindMotifByPattern_EValue_EqualsExpectedRandomCount()
    {
        var match = FindMotifByPattern("AAARGDAAA", "RGD").Single();

        double expected = (9 - 3 + 1) * Math.Pow(2.0, -(3.0 * Log2Of20));
        Assert.That(match.EValue, Is.EqualTo(expected).Within(Tol),
            "E = (N-L+1)*2^(-Score) under uniform background (INV-07)");
    }

    // S3 — unsupported '*' (ScanProsite query extension) must be rejected, not silently dropped.
    [Test]
    public void ConvertPrositeToRegex_KleeneStar_ThrowsFormatException()
    {
        var ex = Assert.Throws<FormatException>(() => ConvertPrositeToRegex("<{C}*>"),
            "Kleene star '*' is not part of the PA-line grammar and must be rejected (INV-06)");
        Assert.That(ex!.Message, Does.Contain("*"), "Exception names the offending '*' construct");
    }

    // S3b — '?' and '+' are likewise not PA-line atoms; reject-don't-drop (INV-06).
    // PROSITE PA-line grammar (scanprosite_doc) has no '?'/'+' metacharacters.
    [Test]
    public void ConvertPrositeToRegex_OtherRegexMetacharacters_ThrowFormatException()
    {
        Assert.Multiple(() =>
        {
            Assert.Throws<FormatException>(() => ConvertPrositeToRegex("A?"),
                "'?' is not a PA-line atom and must be rejected (INV-06)");
            Assert.Throws<FormatException>(() => ConvertPrositeToRegex("A+"),
                "'+' is not a PA-line atom and must be rejected (INV-06)");
        });
    }

    // S4 — matching is case-insensitive: lower/mixed case give the same match as upper case.
    [Test]
    public void FindMotifByPattern_CaseInsensitive_ProducesSamePosition()
    {
        var upper = FindMotifByPattern("AAARGDAAA", "RGD").Single();
        var lower = FindMotifByPattern("aaargdaaa", "RGD").Single();
        var mixed = FindMotifByPattern("AaArGdAaA", "RGD").Single();

        Assert.Multiple(() =>
        {
            Assert.That(lower.Start, Is.EqualTo(upper.Start), "Lowercase yields same start (INV-05)");
            Assert.That(mixed.Start, Is.EqualTo(upper.Start), "Mixed case yields same start (INV-05)");
            Assert.That(lower.Sequence, Is.EqualTo("RGD"), "Matched substring is upper-cased");
        });
    }

    #endregion

    #region Null / empty / invalid inputs and substring invariant (C1, C2, INV-01, INV-02)

    // C1 — null/empty inputs return empty enumerations (no throw) across the methods.
    [Test]
    public void PatternMethods_NullOrEmptyInputs_ReturnEmptyWithoutThrowing()
    {
        Assert.Multiple(() =>
        {
            Assert.That(FindMotifByPattern(null!, "RGD"), Is.Empty, "Null sequence -> empty (INV-01)");
            Assert.That(FindMotifByPattern("", "RGD"), Is.Empty, "Empty sequence -> empty (INV-01)");
            Assert.That(FindMotifByPattern("AAARGDAAA", ""), Is.Empty, "Empty pattern -> empty (INV-01)");
            Assert.That(FindMotifByPattern("AAARGDAAA", null!), Is.Empty, "Null pattern -> empty (INV-01)");
            Assert.That(ConvertPrositeToRegex(""), Is.EqualTo(""), "Empty PROSITE -> empty regex");
            Assert.That(ConvertPrositeToRegex(null!), Is.EqualTo(""), "Null PROSITE -> empty regex");
            Assert.That(FindMotifByProsite("", "R-G-D"), Is.Empty, "Empty sequence -> empty (INV-01)");
            Assert.That(FindDomains(null!), Is.Empty, "Null sequence -> empty domains (INV-01)");
            Assert.That(FindDomains(""), Is.Empty, "Empty sequence -> empty domains (INV-01)");
        });
    }

    // INV-01 — an invalid .NET regex yields an empty enumeration rather than throwing.
    [Test]
    public void FindMotifByPattern_InvalidRegex_ReturnsEmpty()
    {
        Assert.That(FindMotifByPattern("AAARGDAAA", "[unclosed"), Is.Empty,
            "Malformed regex is caught and yields no matches (INV-01)");
    }

    // C2 — every reported match's Sequence equals the substring at [Start..End] (INV-02).
    [Test]
    public void FindMotifByPattern_AllMatches_SatisfySubstringInvariant()
    {
        const string sequence = "RGDxxRGDxxRGD";
        foreach (var m in FindMotifByPattern(sequence, "RGD"))
        {
            string expected = sequence.Substring(m.Start, m.End - m.Start + 1).ToUpperInvariant();
            Assert.That(m.Sequence, Is.EqualTo(expected),
                $"Match at {m.Start} must equal substring [Start..End] (INV-02)");
        }
    }

    #endregion
}
