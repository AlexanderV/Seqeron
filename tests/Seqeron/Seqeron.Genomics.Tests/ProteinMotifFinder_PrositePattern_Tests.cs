using System.Linq;
using NUnit.Framework;
using Seqeron.Genomics.Analysis;
using static Seqeron.Genomics.Analysis.ProteinMotifFinder;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// PROTMOTIF-PROSITE-001: PROSITE Pattern Matching — canonical test file.
/// Validates ConvertPrositeToRegex and FindMotifByProsite against the PROSITE
/// pattern specification (PROSITE User Manual §IV.E, https://prosite.expasy.org/prosuser.html).
/// Evidence: Hulo et al. (2007), De Castro et al. (2006), ScanProsite verified datasets.
/// </summary>
[TestFixture]
public class ProteinMotifFinder_PrositePattern_Tests
{
    #region PROSITE Pattern Constants — verified against official PROSITE entries

    // PROSITE PS00001: N-glycosylation site — https://prosite.expasy.org/PS00001
    private const string Ps00001Prosite = "N-{P}-[ST]-{P}";
    private const string Ps00001Regex = @"N[^P][ST][^P]";

    // PROSITE PS00004: cAMP/cGMP-dependent phosphorylation site — https://prosite.expasy.org/PS00004
    private const string Ps00004Prosite = "[RK](2)-x-[ST]";
    private const string Ps00004Regex = @"[RK]{2}.[ST]";

    // PROSITE PS00005: Protein kinase C phosphorylation site — https://prosite.expasy.org/PS00005
    private const string Ps00005Prosite = "[ST]-x-[RK]";
    private const string Ps00005Regex = @"[ST].[RK]";

    // PROSITE PS00008: N-myristoylation site — https://prosite.expasy.org/PS00008
    private const string Ps00008Prosite = "G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}";
    private const string Ps00008Regex = @"G[^EDRKHPFYW].{2}[STAGCN][^P]";

    // PROSITE PS00016: Cell attachment sequence (RGD) — https://prosite.expasy.org/PS00016
    private const string Ps00016Prosite = "R-G-D";
    private const string Ps00016Regex = @"RGD";

    // PROSITE PS00018: EF-hand calcium-binding domain — https://prosite.expasy.org/PS00018
    private const string Ps00018Prosite = "D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]";
    private const string Ps00018Regex = @"D[^W][DNS][^ILVFYW][DENSTG][DNQGHRK][^GP][LIVMC][DENQSTAGC].{2}[DE][LIVMFYW]";

    // PROSITE PS00028: Zinc finger C2H2 type — https://prosite.expasy.org/PS00028
    private const string Ps00028Prosite = "C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H";
    private const string Ps00028Regex = @"C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H";

    // PROSITE PS00029: Leucine zipper pattern — https://prosite.expasy.org/PS00029
    private const string Ps00029Prosite = "L-x(6)-L-x(6)-L-x(6)-L";
    private const string Ps00029Regex = @"L.{6}L.{6}L.{6}L";

    // PROSITE PS00267: Tachykinin family signature — https://prosite.expasy.org/PS00267
    // Pattern uses [G>] meaning "G or C-terminus" (PROSITE User Manual §IV.E)
    private const string Ps00267Prosite = "F-[IVFY]-G-[LM]-M-[G>].";
    private const string Ps00267Regex = @"F[IVFY]G[LM]M(?:G|$)";

    // PROSITE PS00539: Pyrokinins signature — https://prosite.expasy.org/PS00539
    // Pattern uses [G>] meaning "G or C-terminus" (PROSITE User Manual §IV.E)
    private const string Ps00539Prosite = "F-[GSTV]-P-R-L-[G>].";
    private const string Ps00539Regex = @"F[GSTV]PRL(?:G|$)";

    /// <summary>
    /// Human Transferrin (P02787, TRFE_HUMAN) — UniProt canonical sequence (698 aa).
    /// Used as real published dataset for N-glycosylation site matching.
    /// ScanProsite confirms exactly 2 PS00001 matches at positions 432–435 and 630–633 (1-based).
    /// Source: https://www.uniprot.org/uniprot/P02787
    /// </summary>
    private const string HumanTransferrinSequence =
        "MRLAVGALLVCAVLGLCLAVPDKTVRWCAVSEHEATKCQSFRDHMKSVIPSDGPSVACVK" + // 1-60
        "KASYLDCIRAIAANEADAVTLDAGLVYDAYLAPNNLKPVVAEFYGSKEDPQTFYYAVAVV" + // 61-120
        "KKDSGFQMNQLRGKKSCHTGLGRSAGWNIPIGLLYCDLPEPRKPLEKAVANFFSGSCAPC" + // 121-180
        "ADGTDFPQLCQLCPGCGCSTLNQYFGYSGAFKCLKDGAGDVAFVKHSTIFENLANKADRD" + // 181-240
        "QYELLCLDNTRKPVDEYKDCHLAQVPSHTVVARSMGGKEDLIWELLNQAQEHFGKDKSKE" + // 241-300
        "FQLFSSPHGKDLLFKDSAHGFLKVPPRMDAKMYLGYEYVTAIRNLREGTCPEAPTDECKP" + // 301-360
        "VKWCALSHHERLKCDEWSVNSVGKIECVSAETTEDCIAKIMNGEADAMSLDGGFVYIAGK" +  // 361-420
        "CGLVPVLAENYNKSDNCEDTPEAGYFAIAVVKKSASDLTWDNLKGKKSCHTAVGRTAGWN" + // 421-480
        "IPMGLLYNKINHCRFDEFFSEGCAPGSKKDSSLCKLCMGSGLNLCEPNNKEGYYGYTGAF" + // 481-540
        "RCLVEKGDVAFVKHQTVPQNTGGKNPDPWAKNLNEKDYELLCLDGTRKPVEEYANCHLAR" + // 541-600
        "APNHAVVTRKDKEACVHKILRQQQHLFGSNVTDCSGNFCLFRSETKDLLFRDDTVCLAKL" + // 601-660
        "HDRNTYEKYLGEEYVKAVGNLRKCSTSSLLEACTFRRP";                          // 661-698

    /// <summary>ScanProsite-verified position of first N-glycosylation site (1-based → 0-based).</summary>
    private const int TransferrinNGlyco1Start0Based = 431;

    /// <summary>ScanProsite-verified position of second N-glycosylation site (1-based → 0-based).</summary>
    private const int TransferrinNGlyco2Start0Based = 629;

    /// <summary>Expected matched sequence at first N-glycosylation site.</summary>
    private const string TransferrinNGlyco1Sequence = "NKSD";

    /// <summary>Expected matched sequence at second N-glycosylation site.</summary>
    private const string TransferrinNGlyco2Sequence = "NVTD";

    /// <summary>Number of N-glycosylation matches per ScanProsite scan of P02787 vs PS00001.</summary>
    private const int TransferrinExpectedNGlycoCount = 2;

    #endregion

    #region M1: Simple Literal Pattern — PROSITE PS00016

    [Test]
    public void ConvertPrositeToRegex_SimpleLiteralRGD_ProducesExactRegex()
    {
        // PS00016: R-G-D (Cell attachment sequence) — https://prosite.expasy.org/PS00016
        string result = ConvertPrositeToRegex(Ps00016Prosite);

        Assert.That(result, Is.EqualTo(Ps00016Regex),
            "R-G-D must convert to RGD per PROSITE pattern syntax: literal letters map directly");
    }

    #endregion

    #region M2: Any Amino Acid Marker — PROSITE User Manual

    [Test]
    public void ConvertPrositeToRegex_AnyAminoAcid_ProducesDot()
    {
        // PROSITE User Manual: "The symbol 'x' is used for a position where any amino acid is accepted."
        string result = ConvertPrositeToRegex("A-x-G");

        Assert.That(result, Is.EqualTo("A.G"),
            "'x' must convert to '.' (any character) per PROSITE User Manual §IV.E");
    }

    #endregion

    #region M3: Exact Repeat x(n) — PROSITE User Manual

    [Test]
    public void ConvertPrositeToRegex_ExactRepeatXN_ProducesQuantifier()
    {
        // PROSITE User Manual: "x(3) corresponds to x-x-x"
        string result = ConvertPrositeToRegex("x(3)-A");

        Assert.That(result, Is.EqualTo(".{3}A"),
            "x(3) must convert to .{3} per PROSITE User Manual");
    }

    #endregion

    #region M4: Range Repeat x(n,m) — PROSITE User Manual

    [Test]
    public void ConvertPrositeToRegex_RangeRepeatXNM_ProducesRangeQuantifier()
    {
        // PROSITE User Manual: "x(2,4) corresponds to x-x or x-x-x or x-x-x-x"
        string result = ConvertPrositeToRegex("A-x(2,4)-G");

        Assert.That(result, Is.EqualTo("A.{2,4}G"),
            "x(2,4) must convert to .{2,4} per PROSITE User Manual");
    }

    #endregion

    #region M5: Character Class [ABC] — PROSITE PS00005

    [Test]
    public void ConvertPrositeToRegex_CharacterClass_PreservesSquareBrackets()
    {
        // PROSITE User Manual: "[ALT] stands for Ala or Leu or Thr"
        // PS00005: [ST]-x-[RK] — https://prosite.expasy.org/PS00005
        string result = ConvertPrositeToRegex(Ps00005Prosite);

        Assert.That(result, Is.EqualTo(Ps00005Regex),
            "Character class [ST]-x-[RK] must convert to [ST].[RK] per PS00005");
    }

    #endregion

    #region M6: Exclusion Class {ABC} — PROSITE PS00001

    [Test]
    public void ConvertPrositeToRegex_ExclusionClass_ProducesNegatedCharClass()
    {
        // PROSITE User Manual: "{AM} stands for any amino acid except Ala and Met"
        // PS00001: N-{P}-[ST]-{P} — https://prosite.expasy.org/PS00001
        string result = ConvertPrositeToRegex(Ps00001Prosite);

        Assert.That(result, Is.EqualTo(Ps00001Regex),
            "Exclusion class {P} must convert to [^P]; full PS00001 pattern must match expected regex");
    }

    #endregion

    #region M7: N-Terminus Anchor — PROSITE User Manual

    [Test]
    public void ConvertPrositeToRegex_NTerminusAnchor_ProducesCaret()
    {
        // PROSITE User Manual: "When a pattern is restricted to the N-terminal, it starts with '<'"
        string result = ConvertPrositeToRegex("<M-x-K");

        Assert.That(result, Is.EqualTo("^M.K"),
            "'<' must convert to '^' (start anchor) per PROSITE User Manual");
    }

    #endregion

    #region M8: C-Terminus Anchor — PROSITE User Manual

    [Test]
    public void ConvertPrositeToRegex_CTerminusAnchor_ProducesDollar()
    {
        // PROSITE User Manual: "When a pattern is restricted to the C-terminal, it ends with '>'"
        string result = ConvertPrositeToRegex("A-x-G>");

        Assert.That(result, Is.EqualTo("A.G$"),
            "'>' must convert to '$' (end anchor) per PROSITE User Manual");
    }

    #endregion

    #region M9: Element Repetition A(n) — PROSITE PS00004

    [Test]
    public void ConvertPrositeToRegex_ElementRepetition_ProducesQuantifier()
    {
        // PROSITE User Manual: repetition with parentheses after any element
        // PS00004: [RK](2)-x-[ST] — https://prosite.expasy.org/PS00004
        string result = ConvertPrositeToRegex(Ps00004Prosite);

        Assert.That(result, Is.EqualTo(Ps00004Regex),
            "[RK](2) must produce [RK]{2} per PROSITE pattern repetition syntax (PS00004)");
    }

    #endregion

    #region M10: Complex PS00028 (Zinc Finger C2H2) — PROSITE PS00028

    [Test]
    public void ConvertPrositeToRegex_ComplexZincFinger_ProducesFullRegex()
    {
        // PS00028: C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H
        // https://prosite.expasy.org/PS00028
        string result = ConvertPrositeToRegex(Ps00028Prosite);

        Assert.That(result, Is.EqualTo(Ps00028Regex),
            "PS00028 zinc finger pattern must convert to correct regex with variable ranges");
    }

    #endregion

    #region M11: Complex PS00008 (N-myristoylation) — PROSITE PS00008

    [Test]
    public void ConvertPrositeToRegex_ComplexMyristoylation_ProducesFullRegex()
    {
        // PS00008: G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}
        // https://prosite.expasy.org/PS00008
        string result = ConvertPrositeToRegex(Ps00008Prosite);

        Assert.That(result, Is.EqualTo(Ps00008Regex),
            "PS00008 N-myristoylation pattern must convert including exclusion and character classes");
    }

    #endregion

    #region M12: Complex PS00018 (EF-hand) — PROSITE PS00018

    [Test]
    public void ConvertPrositeToRegex_ComplexEFHand_ProducesFullRegex()
    {
        // PS00018: D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW]
        // https://prosite.expasy.org/PS00018
        string result = ConvertPrositeToRegex(Ps00018Prosite);

        Assert.That(result, Is.EqualTo(Ps00018Regex),
            "PS00018 EF-hand pattern with 13 elements must convert correctly");
    }

    #endregion

    #region M13: Empty Input — Edge Case

    [Test]
    public void ConvertPrositeToRegex_EmptyString_ReturnsEmpty()
    {
        Assert.That(ConvertPrositeToRegex(""), Is.EqualTo(""),
            "Empty PROSITE pattern must return empty regex (INV-1)");
    }

    [Test]
    public void ConvertPrositeToRegex_Null_ReturnsEmpty()
    {
        Assert.That(ConvertPrositeToRegex(null!), Is.EqualTo(""),
            "Null PROSITE pattern must return empty regex");
    }

    #endregion

    #region M14: Trailing Period — PROSITE User Manual

    [Test]
    public void ConvertPrositeToRegex_TrailingPeriod_IsIgnored()
    {
        // PROSITE User Manual: "A period ends the pattern" in PA lines
        string result = ConvertPrositeToRegex("R-G-D.");

        Assert.That(result, Is.EqualTo("RGD"),
            "Trailing period in PROSITE data file format must be silently stripped");
    }

    #endregion

    #region M15: FindMotifByProsite N-glycosylation Match — PS00001

    [Test]
    public void FindMotifByProsite_NGlycosylation_MatchesAtCorrectPositions()
    {
        // Synthetic sequence with two N-glycosylation sites:
        // Position 2 (0-based): NAS → N(2), A(3)≠P, S(4)∈[ST], A(5)≠P → match "NASA"
        // Position 8 (0-based): NGT → N(8), G(9)≠P, T(10)∈[ST], A(11)≠P → match "NGTA"
        const string sequence = "AANASAAANGTAAAA";
        var matches = FindMotifByProsite(sequence, Ps00001Prosite, "N-glycosylation").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.EqualTo(2),
                "Must find exactly 2 N-glycosylation sites per pattern N-{P}-[ST]-{P}");
            Assert.That(matches[0].Start, Is.EqualTo(2),
                "First match (NASA) must start at 0-based position 2");
            Assert.That(matches[0].Sequence, Is.EqualTo("NASA"),
                "First match must capture the 4-character N-glycosylation sequon NASA");
            Assert.That(matches[1].Start, Is.EqualTo(8),
                "Second match (NGTA) must start at 0-based position 8");
            Assert.That(matches[1].Sequence, Is.EqualTo("NGTA"),
                "Second match must capture the 4-character N-glycosylation sequon NGTA");
        });
    }

    #endregion

    #region M16: FindMotifByProsite RGD Match — PS00016

    [Test]
    public void FindMotifByProsite_RGDPattern_MatchesAtExactPosition()
    {
        // PS00016: R-G-D (Cell attachment sequence) at 0-based position 3
        const string sequence = "AAARGDAAA";
        var matches = FindMotifByProsite(sequence, Ps00016Prosite, "RGD").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.EqualTo(1),
                "Must find exactly 1 RGD match in sequence");
            Assert.That(matches[0].Start, Is.EqualTo(3),
                "RGD must start at 0-based position 3");
            Assert.That(matches[0].End, Is.EqualTo(5),
                "RGD End must be 5 (inclusive 0-based, length 3)");
            Assert.That(matches[0].Sequence, Is.EqualTo("RGD"),
                "Matched sequence must be the 3-character RGD motif");
        });
    }

    #endregion

    #region M17: FindMotifByProsite No Match — Exclusion Blocks

    [Test]
    public void FindMotifByProsite_ExclusionBlocks_ReturnsNoMatches()
    {
        // N-{P}-[ST]-{P}: if second position is P, pattern must NOT match
        const string sequence = "AANPSAAA";
        var matches = FindMotifByProsite(sequence, Ps00001Prosite, "N-glycosylation").ToList();

        Assert.That(matches, Is.Empty,
            "Pattern must not match when excluded residue P appears at position {P}");
    }

    #endregion

    #region M18: FindMotifByProsite Empty Sequence — Edge Case

    [Test]
    public void FindMotifByProsite_EmptySequence_ReturnsNoMatches()
    {
        var matches = FindMotifByProsite("", Ps00016Prosite, "RGD").ToList();

        Assert.That(matches, Is.Empty,
            "Empty sequence must produce no matches (INV-3)");
    }

    #endregion

    #region M19: FindMotifByProsite Empty Pattern — Edge Case

    [Test]
    public void FindMotifByProsite_EmptyPattern_ReturnsNoMatches()
    {
        var matches = FindMotifByProsite("AAARGDAAA", "", "Custom").ToList();

        Assert.That(matches, Is.Empty,
            "Empty pattern must produce no matches (INV-4)");
    }

    #endregion

    #region M20: FindMotifByProsite Case Insensitivity — INV-5

    [Test]
    public void FindMotifByProsite_CaseInsensitive_ProducesSameResults()
    {
        // INV-5: FindMotifByProsite is case-insensitive
        const string upperSeq = "AAARGDAAA";
        const string lowerSeq = "aaargdaaa";
        const string mixedSeq = "AaArGdAaA";

        var upperMatches = FindMotifByProsite(upperSeq, Ps00016Prosite, "RGD").ToList();
        var lowerMatches = FindMotifByProsite(lowerSeq, Ps00016Prosite, "RGD").ToList();
        var mixedMatches = FindMotifByProsite(mixedSeq, Ps00016Prosite, "RGD").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(upperMatches, Has.Count.EqualTo(1),
                "Uppercase sequence must find 1 RGD match");
            Assert.That(lowerMatches, Has.Count.EqualTo(1),
                "Lowercase sequence must find 1 RGD match (case-insensitive)");
            Assert.That(mixedMatches, Has.Count.EqualTo(1),
                "Mixed-case sequence must find 1 RGD match (case-insensitive)");
            Assert.That(lowerMatches[0].Start, Is.EqualTo(upperMatches[0].Start),
                "Match position must be identical regardless of case");
        });
    }

    #endregion

    #region M21: FindMotifByProsite Multiple Matches — PS00005

    [Test]
    public void FindMotifByProsite_MultipleMatches_FindsAllPositions()
    {
        // PS00005: [ST]-x-[RK] — PKC phosphorylation site
        // Sequence designed with two non-overlapping sites:
        // Position 2: SAR (S∈[ST], A=any, R∈[RK])
        // Position 7: TGK (T∈[ST], G=any, K∈[RK])
        const string sequence = "AASARAATGKAAA";
        var matches = FindMotifByProsite(sequence, Ps00005Prosite, "PKC").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.GreaterThanOrEqualTo(2),
                "Must find at least 2 PKC phosphorylation sites");

            var first = matches.First(m => m.Start == 2);
            Assert.That(first.Sequence, Is.EqualTo("SAR"),
                "First PKC site at position 2 must match SAR");

            var second = matches.First(m => m.Start == 7);
            Assert.That(second.Sequence, Is.EqualTo("TGK"),
                "Second PKC site at position 7 must match TGK");
        });
    }

    #endregion

    #region M22: C-Terminus Inside Brackets [G>] — PROSITE PS00267, PS00539

    [Test]
    public void ConvertPrositeToRegex_PS00267_CTermInsideBrackets_ConvertsCorrectly()
    {
        // PROSITE User Manual §IV.E: "'>' can also occur inside square brackets for
        // the C-terminal element. 'F-[GSTV]-P-R-L-[G>]' means that either
        // 'F-[GSTV]-P-R-L-G' or 'F-[GSTV]-P-R-L>' are considered."
        // PS00267 (Tachykinin): F-[IVFY]-G-[LM]-M-[G>].
        string result = ConvertPrositeToRegex(Ps00267Prosite);
        Assert.That(result, Is.EqualTo(Ps00267Regex),
            "[G>] must produce (?:G|$) alternation — 'G or end-of-sequence'");
    }

    [Test]
    public void ConvertPrositeToRegex_PS00539_CTermInsideBrackets_ConvertsCorrectly()
    {
        // PS00539 (Pyrokinins): F-[GSTV]-P-R-L-[G>].
        string result = ConvertPrositeToRegex(Ps00539Prosite);
        Assert.That(result, Is.EqualTo(Ps00539Regex),
            "[G>] must produce (?:G|$) alternation — 'G or end-of-sequence'");
    }

    #endregion

    #region M23: FindMotifByProsite [G>] Match — PS00267

    [Test]
    public void FindMotifByProsite_Tachykinin_MatchesWithG()
    {
        // PS00267: F-[IVFY]-G-[LM]-M-[G>]
        // When the final position is G (not at C-terminus), the pattern matches normally.
        // Synthetic sequence embedding FVGLMG (matches [G>] via G branch)
        const string sequence = "AAAFVGLMGAAA";
        var matches = FindMotifByProsite(sequence, Ps00267Prosite, "Tachykinin").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.EqualTo(1),
                "Must find exactly 1 tachykinin match with terminal G");
            Assert.That(matches[0].Start, Is.EqualTo(3),
                "Tachykinin match must start at position 3 (F)");
            Assert.That(matches[0].Sequence, Is.EqualTo("FVGLMG"),
                "Matched sequence must be FVGLMG (G at final position)");
        });
    }

    [Test]
    public void FindMotifByProsite_Tachykinin_MatchesAtCTerminus()
    {
        // PS00267: F-[IVFY]-G-[LM]-M-[G>]
        // When the sequence ends right after position 5, the [G>] matches via C-terminus '>' branch.
        // Pattern matches 5 residues FVGLM at end of sequence ($ matches via [G>]).
        const string sequence = "AAAFVGLM";
        var matches = FindMotifByProsite(sequence, Ps00267Prosite, "Tachykinin").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.EqualTo(1),
                "Must find exactly 1 tachykinin match at C-terminus (> branch)");
            Assert.That(matches[0].Start, Is.EqualTo(3),
                "Tachykinin C-terminal match must start at position 3 (F)");
            Assert.That(matches[0].Sequence, Is.EqualTo("FVGLM"),
                "Matched sequence at C-terminus must be FVGLM (no G, ends at sequence boundary)");
        });
    }

    [Test]
    public void FindMotifByProsite_Tachykinin_NoMatchInMiddleWithoutG()
    {
        // PS00267: F-[IVFY]-G-[LM]-M-[G>]
        // In the middle of a sequence, [G>] requires G — the '>' branch ($) won't match mid-sequence.
        const string sequence = "AAAFVGLMAAAA";
        var matches = FindMotifByProsite(sequence, Ps00267Prosite, "Tachykinin").ToList();

        Assert.That(matches, Has.Count.EqualTo(0),
            "Must NOT match mid-sequence without G at final position — [G>] only matches G or C-terminus");
    }

    #endregion

    #region S1: Real Protein — Human Transferrin PS00001 (ScanProsite Verified)

    [Test]
    public void FindMotifByProsite_HumanTransferrin_FindsExactlyTwoNGlycoSites()
    {
        // ScanProsite scan of P02787 (TRFE_HUMAN) against PS00001 returns exactly 2 matches
        // at positions 432–435 and 630–633 (1-based).
        // Source: https://prosite.expasy.org/cgi-bin/prosite/scanprosite/PSScan.cgi?seq=P02787&sig=PS00001&output=json
        var matches = FindMotifByProsite(HumanTransferrinSequence, Ps00001Prosite, "N-glycosylation")
            .ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.EqualTo(TransferrinExpectedNGlycoCount),
                $"Human Transferrin must have exactly {TransferrinExpectedNGlycoCount} N-glycosylation sites per ScanProsite");

            // Verify each match starts with N and satisfies N-{P}-[ST]-{P}
            foreach (var match in matches)
            {
                Assert.That(match.Sequence[0], Is.EqualTo('N'),
                    $"N-glycosylation match at position {match.Start} must start with N");
                Assert.That(match.Sequence[1], Is.Not.EqualTo('P'),
                    $"N-glycosylation match at position {match.Start}: second residue must not be P");
                Assert.That(match.Sequence[2], Is.AnyOf('S', 'T'),
                    $"N-glycosylation match at position {match.Start}: third residue must be S or T");
                Assert.That(match.Sequence[3], Is.Not.EqualTo('P'),
                    $"N-glycosylation match at position {match.Start}: fourth residue must not be P");
            }

            // Verify exact ScanProsite-confirmed positions (1-based 432, 630 → 0-based 431, 629)
            Assert.That(matches[0].Start, Is.EqualTo(TransferrinNGlyco1Start0Based),
                "First N-glycosylation site must be at 0-based position 431 (ScanProsite: 432, 1-based)");
            Assert.That(matches[0].Sequence, Is.EqualTo(TransferrinNGlyco1Sequence),
                "First N-glycosylation matched sequence must be NKSD");
            Assert.That(matches[1].Start, Is.EqualTo(TransferrinNGlyco2Start0Based),
                "Second N-glycosylation site must be at 0-based position 629 (ScanProsite: 630, 1-based)");
            Assert.That(matches[1].Sequence, Is.EqualTo(TransferrinNGlyco2Sequence),
                "Second N-glycosylation matched sequence must be NVTD");
        });
    }

    #endregion

    #region S2: N-Terminal Anchored Pattern — PROSITE User Manual

    [Test]
    public void FindMotifByProsite_NTerminalAnchor_OnlyMatchesAtStart()
    {
        // PROSITE User Manual: "When a pattern is restricted to the N-terminal, it starts with '<'"
        // Pattern <M-x-K matches only if MxK is at the very beginning of the sequence
        const string matchingSeq = "MAKAAAA";    // MxK at position 0
        const string nonMatchingSeq = "AAMAKAA"; // MxK NOT at position 0

        var matchResults = FindMotifByProsite(matchingSeq, "<M-x-K", "N-term").ToList();
        var noMatchResults = FindMotifByProsite(nonMatchingSeq, "<M-x-K", "N-term").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matchResults, Has.Count.EqualTo(1),
                "N-terminal anchored pattern must match when MxK starts at position 0");
            Assert.That(matchResults[0].Start, Is.EqualTo(0),
                "N-terminal anchored match must start at position 0");
            Assert.That(noMatchResults, Is.Empty,
                "N-terminal anchored pattern must NOT match when MxK is not at the start");
        });
    }

    #endregion

    #region S3: C-Terminal Anchored Pattern — PROSITE User Manual

    [Test]
    public void FindMotifByProsite_CTerminalAnchor_OnlyMatchesAtEnd()
    {
        // PROSITE User Manual: "When a pattern is restricted to the C-terminal, it ends with '>'"
        const string matchingSeq = "AAAALRG";      // LRG at the end
        const string nonMatchingSeq = "AALRGAAA";  // LRG NOT at the end

        var matchResults = FindMotifByProsite(matchingSeq, "L-R-G>", "C-term").ToList();
        var noMatchResults = FindMotifByProsite(nonMatchingSeq, "L-R-G>", "C-term").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matchResults, Has.Count.EqualTo(1),
                "C-terminal anchored pattern must match when LRG ends the sequence");
            Assert.That(matchResults[0].Sequence, Is.EqualTo("LRG"),
                "C-terminal anchored match must capture LRG");
            Assert.That(noMatchResults, Is.Empty,
                "C-terminal anchored pattern must NOT match when LRG is not at the end");
        });
    }

    #endregion

    #region S4: Matched Sequence String — Implementation Contract

    [Test]
    public void FindMotifByProsite_MatchedSequence_IsExactSubstring()
    {
        // Verify that the Sequence field of each match is the exact substring from the input
        const string sequence = "AANASAAAA";
        var matches = FindMotifByProsite(sequence, Ps00001Prosite, "N-glycosylation").ToList();

        Assert.That(matches, Has.Count.EqualTo(1),
            "Must find exactly 1 N-glycosylation site");

        var match = matches[0];
        // The matched sequence must be extractable from the original sequence at the reported position
        string expected = sequence.Substring(match.Start, match.End - match.Start + 1).ToUpperInvariant();

        Assert.That(match.Sequence, Is.EqualTo(expected),
            "Matched sequence string must equal the exact substring at [Start..End] from the input");
    }

    #endregion

    #region S5: Leucine Zipper Pattern — PROSITE PS00029

    [Test]
    public void ConvertPrositeToRegex_LeucineZipperLongRepeat_ProducesCorrectRegex()
    {
        // PS00029: L-x(6)-L-x(6)-L-x(6)-L — https://prosite.expasy.org/PS00029
        string result = ConvertPrositeToRegex(Ps00029Prosite);

        Assert.That(result, Is.EqualTo(Ps00029Regex),
            "PS00029 leucine zipper with repeated x(6) segments must convert correctly");
    }

    #endregion

    #region C1: Zinc Finger Complex Pattern Match — PROSITE PS00028

    [Test]
    public void FindMotifByProsite_ZincFingerC2H2_MatchesSyntheticDomain()
    {
        // PS00028: C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H
        // Synthetic zinc finger domain: CxxCxxxLxxxxxxxxHxxxH (shortest variant, 23 residues)
        // C(0) AA(1-2) C(3) AAA(4-6) L(7) AAAAAAAA(8-15) H(16) AAA(17-19) H(20)
        const string zincFinger = "AAACAACAAALAAAAAAAAHAAAHAAAA";
        // The zinc finger starts at position 3: C-AA-C-AAA-L-AAAAAAAA-H-AAA-H
        var matches = FindMotifByProsite(zincFinger, Ps00028Prosite, "Zinc Finger").ToList();

        Assert.Multiple(() =>
        {
            Assert.That(matches, Has.Count.GreaterThanOrEqualTo(1),
                "Must find at least 1 zinc finger C2H2 domain match");
            Assert.That(matches[0].Start, Is.EqualTo(3),
                "Zinc finger match must start at position 3 (first C)");
            Assert.That(matches[0].Sequence, Does.StartWith("C"),
                "Zinc finger matched sequence must begin with first Cys");
            Assert.That(matches[0].Sequence, Does.EndWith("H"),
                "Zinc finger matched sequence must end with second His");
        });
    }

    #endregion

    #region Match Metadata Tests

    [Test]
    public void FindMotifByProsite_MatchProperties_AllFieldsPopulated()
    {
        // Every MotifMatch from FindMotifByProsite must have all fields populated
        const string sequence = "AAARGDAAA";
        var matches = FindMotifByProsite(sequence, Ps00016Prosite, "TestMotif").ToList();

        Assert.That(matches, Has.Count.EqualTo(1), "Must find 1 RGD match");

        var match = matches[0];
        Assert.Multiple(() =>
        {
            Assert.That(match.Start, Is.GreaterThanOrEqualTo(0),
                "Start must be non-negative");
            Assert.That(match.End, Is.GreaterThanOrEqualTo(match.Start),
                "End must be >= Start");
            Assert.That(match.Sequence, Is.Not.Empty,
                "Sequence must not be empty");
            Assert.That(match.MotifName, Is.EqualTo("TestMotif"),
                "MotifName must match the provided name");
            Assert.That(match.Pattern, Is.EqualTo(Ps00016Prosite),
                "Pattern must contain the original PROSITE pattern");
            Assert.That(match.Score, Is.GreaterThan(0),
                "Score must be positive for a non-trivial pattern");
            Assert.That(match.EValue, Is.GreaterThanOrEqualTo(0),
                "E-value must be non-negative");
        });
    }

    #endregion
}
