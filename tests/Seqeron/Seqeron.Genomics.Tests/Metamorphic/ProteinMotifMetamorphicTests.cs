using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics.Analysis;

namespace Seqeron.Genomics.Tests;

/// <summary>
/// Metamorphic tests for the ProteinMotif area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION, not from observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PROTMOTIF-FIND-001 — pattern-based protein motif finding (ProteinMotif).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 82.
///
/// API under test (ProteinMotifFinder.FindMotifByPattern):
///   Reports every (overlapping) occurrence of a regex motif in a protein, as (Start, End,
///   matched substring). A match is a purely local property of the residues it spans.
///
/// Relations (derived from local regex matching, NOT from output):
///   • SHIFT (prepend flank): prepending a non-matching flank of length k shifts every match
///          start by exactly k and preserves the matched substrings.
///   • MON  (broader pattern ⇒ ≥ matches): a pattern that generalises another matches a
///          superset of positions.
///   • INV  (flank change ⇒ same motif): editing residues outside a motif occurrence (without
///          creating/removing matches) leaves that occurrence detected at the same position.
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class ProteinMotifMetamorphicTests
{
    private static (int Start, int End, string Seq)[] Matches(string protein, string pattern) =>
        ProteinMotifFinder.FindMotifByPattern(protein, pattern).Select(m => (m.Start, m.End, m.Sequence)).ToArray();

    #region PROTMOTIF-FIND-001 SHIFT — prepending a flank shifts match starts

    [Test]
    [Description("SHIFT: prepending a non-matching flank of length k shifts every match start by exactly k and preserves the matched substrings.")]
    public void FindMotifByPattern_PrependFlank_ShiftsStarts()
    {
        const string protein = "AAACGGCDDD"; // contains one C-x(2)-C motif "CGGC" at index 3
        const string pattern = "C.{2}C";

        var baseMatches = Matches(protein, pattern);
        baseMatches.Should().NotBeEmpty(because: "the protein contains a C-x(2)-C motif");

        foreach (int k in new[] { 1, 3, 5 })
        {
            string flank = new string('W', k); // tryptophans — no cysteine, cannot match or extend the motif
            var shifted = Matches(flank + protein, pattern);

            shifted.Select(m => (m.Start, m.Seq)).Should().Equal(
                baseMatches.Select(m => (m.Start + k, m.Seq)),
                because: $"a non-matching {k}-residue flank shifts every match by {k} while preserving the matched residues");
        }
    }

    #endregion

    #region PROTMOTIF-FIND-001 MON — a broader pattern matches a superset

    [Test]
    [Description("MON: a generalised pattern matches a superset of the positions of the more specific pattern.")]
    public void FindMotifByPattern_BroaderPattern_MatchesSuperset()
    {
        const string protein = "CAACGGCAAS"; // C-x(2)-C twice, plus a C-x(2)-S

        var narrow = Matches(protein, "C.{2}C").Select(m => m.Start).ToHashSet();
        var broad = Matches(protein, "C.{2}[CS]").Select(m => m.Start).ToHashSet();

        broad.IsSupersetOf(narrow).Should().BeTrue(
            because: "C-x(2)-[CS] generalises C-x(2)-C, so it matches at least everywhere the narrower pattern does");
        broad.Count.Should().BeGreaterThan(narrow.Count,
            because: "the broader pattern additionally matches the C-x(2)-S occurrence");
    }

    #endregion

    #region PROTMOTIF-FIND-001 INV — flank edits don't change the detected motif

    [Test]
    [Description("INV: editing residues in the flanks (no cysteine, so no new/removed motif) leaves the motif occurrence detected at the same position with the same residues.")]
    public void FindMotifByPattern_FlankChange_SameMotif()
    {
        const string pattern = "C.{2}C";
        // Motif "CAAC" sits at a fixed index 5 behind a 5-residue prefix; flanks carry no cysteine.
        var baseline = Matches("WWWWW" + "CAAC" + "WWWWW", pattern);
        baseline.Should().ContainSingle().Which.Should().Be((5, 8, "CAAC"));

        foreach (var (prefix, suffix) in new[]
                 {
                     ("AAAAA", "GGGGG"),
                     ("DDDDD", "HHHHH"),
                     ("GGGGG", "AAAAA"),
                 })
        {
            Matches(prefix + "CAAC" + suffix, pattern).Should().Equal(baseline,
                because: "cysteine-free flank edits keep the lengths fixed and create no new motif, so the single occurrence is unchanged");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PROTMOTIF-PROSITE-001 — PROSITE-pattern motif finding (ProteinMotif).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 83.
    //
    // API under test (ProteinMotifFinder.FindMotifByProsite):
    //   Compiles a PROSITE pattern (e.g. "C-x-C") to a regex and reports occurrences. A more
    //   general pattern (wildcard where a specific residue stood) accepts a superset of strings.
    //
    // Relations (derived from PROSITE generalisation, NOT from output):
    //   • SUB  (specific ⊆ general): replacing a fixed residue by the wildcard x can only add
    //          matches, so the specific pattern's positions are a subset of the general one's.
    //   • INV  (non-matching flank ⇒ same detection): appending a flank that cannot match leaves
    //          the detected occurrences unchanged.
    // ───────────────────────────────────────────────────────────────────────────

    private static int[] PrositeStarts(string protein, string prosite) =>
        ProteinMotifFinder.FindMotifByProsite(protein, prosite).Select(m => m.Start).OrderBy(s => s).ToArray();

    #region PROTMOTIF-PROSITE-001 SUB — a specific pattern matches a subset of a general one

    [Test]
    [Description("SUB: generalising a fixed residue to the PROSITE wildcard x admits a superset, so the specific pattern's matches are a subset of the general pattern's.")]
    public void FindMotifByProsite_SpecificPattern_SubsetOfGeneral()
    {
        const string protein = "CGCAC"; // CGC at 0, CAC at 2

        var specific = PrositeStarts(protein, "C-G-C").ToHashSet();
        var general = PrositeStarts(protein, "C-x-C").ToHashSet();

        general.IsSupersetOf(specific).Should().BeTrue(
            because: "C-x-C generalises C-G-C, so it matches at least everywhere C-G-C does");
        general.Count.Should().BeGreaterThan(specific.Count,
            because: "the wildcard additionally matches the C-A-C occurrence");
    }

    #endregion

    #region PROTMOTIF-PROSITE-001 INV — a non-matching flank doesn't change detection

    [Test]
    [Description("INV: appending a flank that cannot match (no cysteine) leaves the detected PROSITE occurrences unchanged.")]
    public void FindMotifByProsite_NonMatchingFlank_SameDetection()
    {
        const string protein = "CGCAC";
        var baseline = PrositeStarts(protein, "C-x-C");

        foreach (var flank in new[] { "WWWWW", "AAAAA", "GHGHGH" }) // cysteine-free, cannot form C-x-C
            PrositeStarts(protein + flank, "C-x-C").Should().Equal(baseline,
                because: "a flank that cannot match adds no occurrences and (appended) shifts none");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PROTMOTIF-DOMAIN-001 — protein domain detection (ProteinMotif).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 84.
    //
    // API under test (ProteinMotifFinder.FindDomains):
    //   Scans a protein for known domain signatures (zinc finger, WD40, SH3, PDZ, kinase) via
    //   their regex patterns. A domain's confidence Score is the pattern's information content
    //   = Σ over conserved positions of log2(20/allowed) — set by the SIGNATURE, not by the
    //   matched substring's length.
    //
    // Relations (derived from local signature matching + IC scoring, NOT from output):
    //   • INV  (domain intact after non-domain insertion): inserting a non-domain segment
    //          elsewhere leaves the domain detected with the same matched residues.
    //   • MON  (more conserved signature ⇒ higher confidence): because confidence is the
    //          signature's information content, a signature that fixes more positions scores
    //          higher. (NOTE: this is "more conserved/longer SIGNATURE", not a longer matched
    //          substring — wildcard positions contribute zero information.)
    // ───────────────────────────────────────────────────────────────────────────

    // A C2H2 zinc-finger matching "C.{2,4}C.{3}[LIVMFYWC].{8}H.{3,5}H":
    //   C AA C AAA L AAAAAAAA H AAA H
    private const string ZincFingerProtein = "CAACAAALAAAAAAAAHAAAH";

    private static System.Collections.Generic.HashSet<string> ZincFingerMatches(string protein) =>
        ProteinMotifFinder.FindDomains(protein)
            .Where(d => d.Name == "Zinc Finger C2H2")
            .Select(d => protein.Substring(d.Start, d.End - d.Start + 1))
            .ToHashSet();

    #region PROTMOTIF-DOMAIN-001 INV — a non-domain insertion leaves the domain intact

    [Test]
    [Description("INV: inserting a non-domain segment up- or downstream leaves the zinc-finger domain detected with the same matched residues.")]
    public void FindDomains_NonDomainInsertion_DomainIntact()
    {
        var baseline = ZincFingerMatches(ZincFingerProtein);
        baseline.Should().NotBeEmpty(because: "the constructed sequence contains a C2H2 zinc finger");

        // 'P'/'G' blocks carry no C…C…H…H zinc-finger signature, so they are non-domain inserts.
        foreach (var (prefix, suffix) in new[] { ("PPPP", ""), ("", "GGGG"), ("PPPP", "GGGG") })
        {
            var matches = ZincFingerMatches(prefix + ZincFingerProtein + suffix);
            matches.IsSupersetOf(baseline).Should().BeTrue(
                because: "a non-domain insertion elsewhere cannot destroy or alter the zinc-finger occurrence");
        }
    }

    #endregion

    #region PROTMOTIF-DOMAIN-001 MON — a more conserved signature scores higher confidence

    [Test]
    [Description("MON: domain confidence is the signature's information content, so a signature fixing more positions scores higher than a less-specified one matching the same window.")]
    public void FindDomains_MoreConservedSignature_HigherConfidence()
    {
        // FindDomains scores each match by CalculateMotifScore (pattern information content);
        // exercise that scoring model directly on the same matching window.
        const string window = "CACAC";

        double lessSpecified = ProteinMotifFinder.FindMotifByPattern(window, "C.C", "domain").First(m => m.Start == 0).Score;
        double moreSpecified = ProteinMotifFinder.FindMotifByPattern(window, "C.C.C", "domain").First(m => m.Start == 0).Score;

        moreSpecified.Should().BeGreaterThan(lessSpecified,
            because: "fixing a third conserved cysteine adds information content, raising the domain confidence");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PROTMOTIF-CC-001 — coiled-coil prediction (ProteinMotif).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 163.
    //
    // API under test (ProteinMotifFinder.PredictCoiledCoils):
    //   Scores heptad a/d hydrophobic-core (I/L/V) occupancy in a sliding window and reports
    //   contiguous above-threshold regions (Lupas 1991; Mason & Arndt 2004).
    //
    // Relations (derived from the heptad-occupancy model, NOT from output):
    //   • INV  (deterministic): the prediction is a pure function of the sequence.
    //   • SHIFT (prepend flank shifts positions): a non-coiled-coil 5' flank moves the region's 3'
    //          boundary by exactly the flank length (the flank cannot extend the core past its end).
    // ───────────────────────────────────────────────────────────────────────────

    // Five perfect heptads (I/L/V at every a and d position) ⇒ one strong coiled-coil region.
    private static readonly string CoiledCoil = string.Concat(Enumerable.Repeat("LAALAAA", 5));

    #region PROTMOTIF-CC-001 INV — the prediction is deterministic

    [Test]
    [Description("INV: PredictCoiledCoils is a pure function, so repeated calls return the identical regions.")]
    public void CoiledCoils_SameSequence_SameRegions()
    {
        ProteinMotifFinder.PredictCoiledCoils(CoiledCoil).ToList()
            .Should().Equal(ProteinMotifFinder.PredictCoiledCoils(CoiledCoil).ToList(),
                because: "the heptad-occupancy prediction has no hidden state");
    }

    #endregion

    #region PROTMOTIF-CC-001 SHIFT — a prepended flank shifts the region's 3' boundary

    [Test]
    [Description("SHIFT: prepending a non-coiled-coil flank shifts the detected region's 3' end by exactly the flank length, keeping the same region count and peak score.")]
    public void CoiledCoils_PrependFlank_ShiftsEndBoundary()
    {
        var original = ProteinMotifFinder.PredictCoiledCoils(CoiledCoil).ToList();
        original.Should().ContainSingle(because: "five perfect heptads form one coiled-coil region");
        int originalEnd = original[0].End;
        double originalScore = original[0].Score;

        foreach (int flankLen in new[] { 5, 10 })
        {
            string flank = new string('P', flankLen); // proline: never a hydrophobic-core a/d residue
            var shifted = ProteinMotifFinder.PredictCoiledCoils(flank + CoiledCoil).ToList();

            shifted.Should().ContainSingle(because: "the proline flank adds no coiled-coil region of its own");
            shifted[0].End.Should().Be(originalEnd + flankLen,
                because: $"the {flankLen}-residue 5' flank shifts the 3' boundary by {flankLen} without extending it");
            shifted[0].Score.Should().BeApproximately(originalScore, 1e-12, because: "the core heptad occupancy is unchanged");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PROTMOTIF-COMMON-001 — common (PROSITE) motif scan (ProteinMotif).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 164.
    //
    // API under test (ProteinMotifFinder.FindCommonMotifs):
    //   Scans a protein against the dictionary of common PROSITE patterns and reports each match.
    //
    // Relations (derived from the per-occurrence reporting, NOT from output):
    //   • MON  (more occurrences ⇒ ≥ support): adding another occurrence of a motif yields a
    //          superset of matches (the per-checklist "more sharing → ≥ support" analog).
    //   • INV  (deterministic / order independent): the match set is a pure function of the
    //          sequence, independent of dictionary scan order.
    // ───────────────────────────────────────────────────────────────────────────

    #region PROTMOTIF-COMMON-001 MON — more motif occurrences yield more matches

    [Test]
    [Description("MON: a second occurrence of the RGD cell-attachment motif adds a match, so the doubled sequence reports more matches than the single one.")]
    public void CommonMotifs_MoreOccurrences_MoreMatches()
    {
        int single = ProteinMotifFinder.FindCommonMotifs("AARGDKKAA").Count();
        int doubled = ProteinMotifFinder.FindCommonMotifs("AARGDKKAARGDKKAA").Count();

        single.Should().BeGreaterThan(0, because: "the sequence contains one RGD motif");
        doubled.Should().BeGreaterThan(single, because: "a second RGD occurrence adds at least one more match");
    }

    #endregion

    #region PROTMOTIF-COMMON-001 INV — the scan is deterministic

    [Test]
    [Description("INV: FindCommonMotifs is a pure function, so repeated scans return the identical matches (independent of dictionary scan order).")]
    public void CommonMotifs_SameSequence_SameMatches()
    {
        const string seq = "AARGDKKAASARKAANFTAAA";
        ProteinMotifFinder.FindCommonMotifs(seq).Select(m => (m.Start, m.End)).ToList()
            .Should().Equal(ProteinMotifFinder.FindCommonMotifs(seq).Select(m => (m.Start, m.End)).ToList(),
                because: "the dictionary scan has no hidden state");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PROTMOTIF-LC-001 — low-complexity regions (SEG) (ProteinMotif).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 165.
    //
    // API under test (ProteinMotifFinder.FindLowComplexityRegions):
    //   SEG (Wootton & Federhen 1993): a window triggers when its entropy ≤ K1 (triggerComplexity);
    //   the segment extends over adjacent windows with entropy ≤ K2 (extensionComplexity).
    //
    // Relations (derived from the SEG threshold model, NOT from output):
    //   • MON  (more permissive trigger ⇒ superset): a region is low-complexity when its entropy is
    //          BELOW the complexity ceiling K1, so raising K1 (the permissive direction for this
    //          complexity-ceiling model) admits a superset of regions — the checklist's
    //          "lower threshold → superset" expressed in SEG's entropy-ceiling semantics.
    //   • SHIFT (prepend flank shifts regions): a high-complexity 5' flank relabels positions, so the
    //          regions reappear shifted by the flank length with unchanged complexity.
    // ───────────────────────────────────────────────────────────────────────────

    // A homopolymer run (entropy 0) and a ternary run (entropy log2 3 ≈ 1.585), separated by a
    // high-complexity spacer so they are distinct SEG segments.
    private const string LcHomo = "AAAAAAAAAAAA";
    private const string LcSpacer = "CDEFGHIKLMNP";
    private const string LcTernary = "DEFDEFDEFDEF";

    private static System.Collections.Generic.HashSet<int> CoveredPositions(
        System.Collections.Generic.IEnumerable<(int Start, int End, double Complexity)> regions)
    {
        var set = new System.Collections.Generic.HashSet<int>();
        foreach (var r in regions)
            for (int p = r.Start; p <= r.End; p++) set.Add(p);
        return set;
    }

    #region PROTMOTIF-LC-001 MON — raising the entropy ceiling admits a superset

    [Test]
    [Description("MON: low-complexity means entropy ≤ K1, so raising the trigger ceiling K1 (fixed K2) admits a superset of low-complexity coverage — the homopolymer triggers at any K1, the ternary run only once K1 reaches its entropy.")]
    public void LowComplexity_HigherTriggerCeiling_Superset()
    {
        string seq = LcHomo + LcSpacer + LcTernary;
        const int window = 12;
        const double k2 = 2.0; // extension ceiling above the ternary entropy (~1.585), below the spacer

        var strict = CoveredPositions(ProteinMotifFinder.FindLowComplexityRegions(seq, window, triggerComplexity: 0.5, extensionComplexity: k2));
        var permissive = CoveredPositions(ProteinMotifFinder.FindLowComplexityRegions(seq, window, triggerComplexity: 1.8, extensionComplexity: k2));

        strict.IsSubsetOf(permissive).Should().BeTrue(because: "every region triggered at K1=0.5 still triggers at K1=1.8");
        permissive.Count.Should().BeGreaterThan(strict.Count, because: "the ternary run (entropy ≈1.585) triggers only once K1 ≥ ~1.585");
    }

    #endregion

    #region PROTMOTIF-LC-001 SHIFT — a prepended flank shifts the regions

    [Test]
    [Description("SHIFT: a high-complexity 5' flank relabels positions without adding a low-complexity region, so each region reappears shifted by the flank length with unchanged complexity.")]
    public void LowComplexity_PrependFlank_ShiftsRegions()
    {
        string seq = "ACDEFGHIKLMN" + "QQQQQQQQQQQQ" + "RSTVWYACDEFG"; // isolated homopolymer LC region
        var original = ProteinMotifFinder.FindLowComplexityRegions(seq).ToList();
        original.Should().NotBeEmpty();

        foreach (var flank in new[] { "KLMNPQRSTVWY", "ACDEFGHIKLMNPQRS" }) // high-complexity, distinct residues
        {
            var shifted = ProteinMotifFinder.FindLowComplexityRegions(flank + seq).ToList();
            shifted.Select(r => (r.Start, r.End, r.Complexity))
                .Should().Equal(original.Select(r => (r.Start + flank.Length, r.End + flank.Length, r.Complexity)),
                    because: $"the {flank.Length}-residue high-complexity flank shifts every region by {flank.Length}");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: PROTMOTIF-PATTERN-001 — regex/PROSITE pattern matching (ProteinMotif).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 166.
    //
    // API under test (ProteinMotifFinder.FindMotifByPattern):
    //   Reports every (overlapping) match of a regex pattern in a protein.
    //
    // Relations (derived from regular-language matching, NOT from output):
    //   • SHIFT (prepend flank shifts matches): a flank with no match shifts every match start by the
    //          flank length.
    //   • SUB  (broader pattern ⇒ ≥ matches): if pattern P's language ⊇ pattern Q's at each position
    //          (P is a relaxation of Q), then P matches a superset of Q's start positions.
    // ───────────────────────────────────────────────────────────────────────────

    #region PROTMOTIF-PATTERN-001 SHIFT — a prepended flank shifts the matches

    [Test]
    [Description("SHIFT: prepending a flank that does not match the pattern shifts every match's start position by the flank length.")]
    public void Pattern_PrependFlank_ShiftsMatches()
    {
        const string seq = "AARGDKKRGDAA";
        const string pattern = "RGD";
        var original = ProteinMotifFinder.FindMotifByPattern(seq, pattern).Select(m => m.Start).ToList();
        original.Should().NotBeEmpty();

        foreach (var flank in new[] { "WW", "YPWYP" }) // contain no RGD and form none at the junction
        {
            var shifted = ProteinMotifFinder.FindMotifByPattern(flank + seq, pattern).Select(m => m.Start).ToList();
            shifted.Should().Equal(original.Select(s => s + flank.Length),
                because: $"a {flank.Length}-residue flank with no match relocates every match by {flank.Length}");
        }
    }

    #endregion

    #region PROTMOTIF-PATTERN-001 SUB — a broader pattern matches a superset

    [Test]
    [Description("SUB: 'RG.' is a relaxation of 'RGD' (it accepts any third residue), so every RGD match start is also an RG. match start — the broader pattern's matches are a superset.")]
    public void Pattern_BroaderPattern_Superset()
    {
        const string seq = "AARGDKKRGEKKRGDAA";
        var strict = ProteinMotifFinder.FindMotifByPattern(seq, "RGD").Select(m => m.Start).ToHashSet();
        var broad = ProteinMotifFinder.FindMotifByPattern(seq, "RG.").Select(m => m.Start).ToHashSet();

        strict.IsSubsetOf(broad).Should().BeTrue(because: "'RG.' accepts every string 'RGD' accepts");
        broad.Count.Should().BeGreaterThan(strict.Count, because: "'RG.' additionally matches RGE, which 'RGD' rejects");
    }

    #endregion
}
