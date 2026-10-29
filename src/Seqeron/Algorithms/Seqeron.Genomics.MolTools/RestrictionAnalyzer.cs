namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Analyzes restriction enzyme sites in DNA sequences.
/// Supports the 39 common restriction enzymes (<see cref="Enzymes"/>, the default set of every
/// member that enumerates enzymes), every REBASE enzyme with a defined cut resolved by name
/// (<see cref="RebaseEnzymes"/>, 754 enzymes; <see cref="EnzymeSet"/>), and caller-supplied <see cref="RestrictionEnzyme"/> definitions, which may be
/// non-palindromic and may cut outside their own recognition site (single-sided Type IIS cutters
/// such as BsaI, BsmBI, BbsI, SapI, FokI, MlyI).
/// </summary>
public static partial class RestrictionAnalyzer
{
    #region Built-in Enzyme Database

    private static readonly Dictionary<string, RestrictionEnzyme> _enzymes = new(StringComparer.OrdinalIgnoreCase)
    {
        // Type II restriction enzymes (most common in molecular biology)
        // Format: Name, Recognition sequence, Cut positions (forward, reverse from 5'), Organism,
        // REBASE methylation-site field. The last column is REBASE allenz.104 field <6> copied
        // verbatim (generated from the file, not typed; see RESTR-FIND-001-Evidence.md §5): the base
        // the enzyme's COGNATE methyltransferase modifies to protect its own site — NOT a statement
        // about whether Dam/Dcm/CpG methylation blocks this enzyme. "" = REBASE gives no site.

        // 6-cutter enzymes
        ["EcoRI"] = new("EcoRI", "GAATTC", 1, 5, "Escherichia coli", "3(6)"),
        ["BamHI"] = new("BamHI", "GGATCC", 1, 5, "Bacillus amyloliquefaciens", "5(4)"),
        ["HindIII"] = new("HindIII", "AAGCTT", 1, 5, "Haemophilus influenzae", "1(6)"),
        ["XhoI"] = new("XhoI", "CTCGAG", 1, 5, "Xanthomonas holcicola", "5(6)"),
        ["SalI"] = new("SalI", "GTCGAC", 1, 5, "Streptomyces albus", "5(6)"),
        ["PstI"] = new("PstI", "CTGCAG", 5, 1, "Providencia stuartii", "5(6)"),
        ["SphI"] = new("SphI", "GCATGC", 5, 1, "Streptomyces phaeochromogenes", ""),
        ["KpnI"] = new("KpnI", "GGTACC", 5, 1, "Klebsiella pneumoniae", "4(6)"),
        ["SacI"] = new("SacI", "GAGCTC", 5, 1, "Streptomyces achromogenes", "4(5)"),
        ["XbaI"] = new("XbaI", "TCTAGA", 1, 5, "Xanthomonas badrii", "6(6)"),
        ["NcoI"] = new("NcoI", "CCATGG", 1, 5, "Nocardia corallina", "?(4)"),
        ["NdeI"] = new("NdeI", "CATATG", 2, 4, "Neisseria denitrificans", "4(6)"),
        ["NheI"] = new("NheI", "GCTAGC", 1, 5, "Neisseria mucosa", "6(4)"),
        ["SpeI"] = new("SpeI", "ACTAGT", 1, 5, "Sphaerotilus natans", "1(6)"),
        ["AvrII"] = new("AvrII", "CCTAGG", 1, 5, "Anabaena variabilis", ""),
        ["ClaI"] = new("ClaI", "ATCGAT", 2, 4, "Caryophanon latum", "5(6)"),
        ["AgeI"] = new("AgeI", "ACCGGT", 1, 5, "Agrobacterium gelatinovorum", "3(5)"),

        // 8-cutter enzymes (rare cutters)
        ["NotI"] = new("NotI", "GCGGCCGC", 2, 6, "Nocardia otitidis-caviarum", "2(4)"),
        ["SfiI"] = new("SfiI", "GGCCNNNNNGGCC", 8, 5, "Streptomyces fimbriatus", "?(4)"), // Note: N = any base
        ["PacI"] = new("PacI", "TTAATTAA", 5, 3, "Pseudomonas alcaligenes", ""),
        ["AscI"] = new("AscI", "GGCGCGCC", 2, 6, "Arthrobacter species", "?(5)"),
        ["FseI"] = new("FseI", "GGCCGGCC", 6, 2, "Frankia species", "?(5)"),
        ["SwaI"] = new("SwaI", "ATTTAAAT", 4, 4, "Staphylococcus warneri", "?(6)"),

        // 4-cutter enzymes (frequent cutters)
        ["MspI"] = new("MspI", "CCGG", 1, 3, "Moraxella species", "1(5)"),
        ["HpaII"] = new("HpaII", "CCGG", 1, 3, "Haemophilus parainfluenzae", "2(5)"),
        ["TaqI"] = new("TaqI", "TCGA", 1, 3, "Thermus aquaticus", "4(6)"),
        ["AluI"] = new("AluI", "AGCT", 2, 2, "Arthrobacter luteus", "3(5)"), // Blunt
        ["RsaI"] = new("RsaI", "GTAC", 2, 2, "Rhodopseudomonas sphaeroides", "4(4)"), // Blunt
        ["HaeIII"] = new("HaeIII", "GGCC", 2, 2, "Haemophilus aegyptius", "3(5)"), // Blunt
        ["DpnI"] = new("DpnI", "GATC", 2, 2, "Diplococcus pneumoniae", ""), // Blunt, methylation-dependent
        ["Sau3AI"] = new("Sau3AI", "GATC", 0, 4, "Staphylococcus aureus", "4(5)"),
        ["MboI"] = new("MboI", "GATC", 0, 4, "Moraxella bovis", "2(6)"),

        // Common cloning enzymes
        ["EcoRV"] = new("EcoRV", "GATATC", 3, 3, "Escherichia coli", "2(6)"), // Blunt
        ["SmaI"] = new("SmaI", "CCCGGG", 3, 3, "Serratia marcescens", "2(4)"), // Blunt
        ["HincII"] = new("HincII", "GTYRAC", 3, 3, "Haemophilus influenzae", "5(6)"), // Degenerate, blunt
        ["ScaI"] = new("ScaI", "AGTACT", 3, 3, "Streptomyces caespitosus", "5(4)"), // Blunt
        ["StuI"] = new("StuI", "AGGCCT", 3, 3, "Streptomyces tubercidicus", ""), // Blunt
        ["BglII"] = new("BglII", "AGATCT", 1, 5, "Bacillus globigii", "5(4)"),
        ["ApaI"] = new("ApaI", "GGGCCC", 5, 1, "Acetobacter pasteurianus", "4(5)"),
    };

    /// <summary>
    /// The 39 common laboratory enzymes (<see cref="EnzymeSet.Common"/>) — the default set of every
    /// member that enumerates enzymes without being given a set or a list. For the full REBASE set
    /// (754 enzymes) see <see cref="RebaseEnzymes"/> / <see cref="GetEnzymes(EnzymeSet)"/>.
    /// </summary>
    public static IReadOnlyDictionary<string, RestrictionEnzyme> Enzymes => _enzymes;

    /// <summary>
    /// The largest distance, in bases, that a caller-supplied <see cref="RestrictionEnzyme"/> may leave between its top-strand
    /// and bottom-strand cuts (<c>|CutPositionForward - CutPositionReverse|</c>, the length of its overhang); a larger one is an
    /// <see cref="ArgumentException"/>. The overhang is materialised as a string of that length for every fragment end, so an
    /// unbounded stagger (<c>int.MaxValue</c> was accepted and exhausted memory) is not an input but a denial of service.
    /// The largest stagger among the 754 REBASE enzymes is 9 (TspRI, TscAI), so 64 leaves a sevenfold margin; the bound is a
    /// safety limit of this library, not a property of restriction enzymes.
    /// </summary>
    public const int MaxCutStagger = 64;

    /// <summary>
    /// The longest recognition sequence a caller-supplied <see cref="RestrictionEnzyme"/> may have; a longer one is an
    /// <see cref="ArgumentException"/>. Every site record carries the matched window and the enzyme, so a recognition sequence of
    /// <c>m</c> bases over a molecule of <c>n</c> gives output and memory of order <c>n * m</c> (an <c>N</c> pattern matches everywhere:
    /// a 200 KB pattern over a 2000-nt circle exhausted memory). The longest of the 754 REBASE sites is 15, so 64 leaves a fourfold
    /// margin; like <see cref="MaxCutStagger"/> this is a safety limit of this library, not a property of restriction enzymes.
    /// </summary>
    public const int MaxRecognitionLength = 64;

    /// <summary>
    /// The longest <see cref="RestrictionEnzyme.Name"/> of a caller-supplied enzyme (the name is repeated in every site and fragment
    /// record, so it multiplies the output the same way). The longest REBASE name is 9 characters.
    /// </summary>
    public const int MaxEnzymeNameLength = 64;

    /// <summary>
    /// Gets a restriction enzyme by name, resolved case-insensitively against the full REBASE set
    /// (<see cref="RebaseEnzymes"/>, which contains the 39 common enzymes of <see cref="Enzymes"/> as
    /// the same instances), or <see langword="null"/> when no enzyme of that name has a defined cut.
    /// </summary>
    /// <remarks>
    /// This is the name-resolution contract of every member that takes an enzyme <i>name</i>
    /// (<c>FindSites</c>, <c>Digest</c>, <c>PartialDigest</c>, <c>GetDigestSummary</c>,
    /// <c>CreateMap</c>, <c>AreCompatible</c>, <c>GetOverhangSequence</c>, the isoschizomer and
    /// compatibility queries), so e.g. <c>Digest(seq, "BsaI")</c> works although BsaI is not one of
    /// the 39 common enzymes. A REBASE enzyme whose cut position is not defined (e.g.
    /// <c>Dpi3084I</c>) is not resolvable: it cannot be digested.
    /// </remarks>
    /// <param name="name">Enzyme name (case-insensitive). Surrounding whitespace is ignored, as in
    /// <see cref="MethylaseSets.FromName"/>, and a hyphen is accepted wherever the table carries an
    /// underscore: REBASE spells the shipped enzyme <c>CviKI-1</c> with a hyphen, but the table is
    /// keyed by Biopython's class-name spelling <c>CviKI_1</c>, which cannot contain one (review
    /// 2026-09 B08 H13). Both spellings, and both in any casing, resolve to the same enzyme.</param>
    /// <exception cref="ArgumentNullException"><paramref name="name"/> is <see langword="null"/>.</exception>
    public static RestrictionEnzyme? GetEnzyme(string name)
    {
        ArgumentNullException.ThrowIfNull(name);

        string key = name.Trim();
        if (_enzymes.TryGetValue(key, out var enzyme) || RebaseEnzymes.TryGetValue(key, out enzyme))
            return enzyme;

        // The frozen table is keyed by Biopython's class names, which replace a REBASE hyphen with
        // an underscore. Accept the real REBASE spelling without touching the data file.
        if (key.Contains('-') && RebaseEnzymes.TryGetValue(key.Replace('-', '_'), out enzyme))
            return enzyme;

        return null;
    }

    /// <summary>
    /// Gets all enzymes whose recognition sequence has exactly the given length (in base pairs).
    /// </summary>
    /// <param name="length">Recognition-sequence length to match (Type II sites are 4–8 nt;
    /// see Wikipedia "Restriction enzyme").</param>
    /// <remarks>
    /// The length measured is the full recognition <i>string</i>, i.e. Biopython
    /// <c>Bio.Restriction.&lt;enzyme&gt;.size</c>. For an interrupted (divided) palindrome that is
    /// longer than the number of specified bases: <c>SfiI</c> (<c>GGCCNNNNNGGCC</c>) has
    /// <c>size == 13</c> although only 8 bases are specified, so it is returned by
    /// <c>GetEnzymesByCutLength(13)</c> and not by <c>GetEnzymesByCutLength(8)</c>. All 39
    /// built-in entries were audited against Biopython 1.88 <c>.size</c> (39/39 agreement).
    /// </remarks>
    public static IEnumerable<RestrictionEnzyme> GetEnzymesByCutLength(int length) =>
        GetEnzymesByCutLength(length, EnzymeSet.Common);

    /// <summary>
    /// Gets all enzymes whose recognition-sequence length falls within the inclusive range
    /// [<paramref name="minLength"/>, <paramref name="maxLength"/>] base pairs.
    /// </summary>
    /// <param name="minLength">Inclusive lower bound on recognition-sequence length (bp).</param>
    /// <param name="maxLength">Inclusive upper bound on recognition-sequence length (bp).</param>
    /// <remarks>
    /// Type II restriction recognition sites are 4–8 nucleotides in length
    /// (Wikipedia, "Restriction enzyme"). Both bounds are inclusive; an inverted range
    /// (<paramref name="minLength"/> &gt; <paramref name="maxLength"/>) yields no enzymes.
    /// </remarks>
    public static IEnumerable<RestrictionEnzyme> GetEnzymesByCutLength(int minLength, int maxLength) =>
        GetEnzymesByCutLength(minLength, maxLength, EnzymeSet.Common);

    /// <summary>
    /// Gets all enzymes that produce blunt ends — both strands are cut at the same position
    /// so each strand "terminates in a base pair" (Wikipedia, "Sticky and blunt ends").
    /// </summary>
    /// <remarks>
    /// Equivalent to Biopython <c>Bio.Restriction.&lt;enzyme&gt;.is_blunt()</c>; all 39 built-in
    /// entries were audited against it (39/39 agreement, 10 blunt cutters).
    /// </remarks>
    public static IEnumerable<RestrictionEnzyme> GetBluntCutters() =>
        GetBluntCutters(EnzymeSet.Common);

    /// <summary>
    /// Gets all enzymes that produce sticky (cohesive) ends — a staggered cut leaving a
    /// 5' or 3' single-stranded overhang (Wikipedia, "Sticky and blunt ends").
    /// </summary>
    /// <remarks>
    /// The complement of <see cref="GetBluntCutters"/>; the 5'/3' split is
    /// <see cref="RestrictionEnzyme.OverhangType"/>, audited against Biopython
    /// <c>is_5overhang()</c>/<c>is_3overhang()</c> for all 39 built-in entries (39/39 agreement,
    /// 21 with a 5' and 8 with a 3' overhang).
    /// </remarks>
    public static IEnumerable<RestrictionEnzyme> GetStickyCutters() =>
        GetStickyCutters(EnzymeSet.Common);

    #endregion

    #region Restriction Site Finding

    // ───────────────────────────────────────────────────────────────────────────────────────
    // Coordinate and feasibility conventions (reference: Biopython 1.88 Bio.Restriction,
    // Restriction.py — FormattedSeq.finditer L215-228, AbstractCut._drop L542-568,
    // Defined._modify L789-814 / Defined._rev_modify L817-825, Palindromic._search L967-978,
    // NonPalindromic._search L996-1019).
    //
    // • Every coordinate this class reports is 0-based and in FORWARD-strand (top-strand)
    //   space, for both strands. `Position` is the index of the 5'-most base of the
    //   recognition window on the forward strand. `CutPosition` is the number of forward-strand
    //   bases 5' of the cut, i.e. the cut lies between bases `CutPosition - 1` and
    //   `CutPosition`, on the strand named by `IsForwardStrand`.
    //   Biopython's `enzyme.search(seq)` returns the same cut as a 1-based top-strand position,
    //   so Biopython = this class's forward-strand `CutPosition` + 1.
    //
    // • For a forward-strand site at p the two cuts are `p + CutPositionForward` (top strand)
    //   and `p + CutPositionReverse` (bottom strand). For a site found on the reverse strand
    //   the enzyme is bound the other way round, so — with L the recognition length — the cuts
    //   are `p + L - CutPositionReverse` (top strand) and `p + L - CutPositionForward` (bottom
    //   strand); the reverse-strand record carries the latter, the cut on the strand it was
    //   read from. For a palindromic enzyme CutPositionForward + CutPositionReverse == L, so
    //   both formulas coincide and the output is identical to the historical one; they differ
    //   only for a custom NON-palindromic enzyme, where the historical `p + CutPositionReverse`
    //   was wrong (REBASE has exactly two such in-site cutters, CdiI CATC^G and SimI GG^GTC_N).
    //
    // • Distal (single-sided Type IIS) cutters need nothing extra: the offsets are signed and may
    //   lie outside [0, L] (BsaI GGTCTC(1/5) = offsets 7/11), and the formulas above are
    //   Biopython's own — forward `location + fst5`, reverse `location - fst3`, with
    //   CutPositionForward = fst5 and CutPositionReverse = size + fst3. For a non-palindromic
    //   site the reverse-strand occurrence therefore cuts on the OTHER side of the window
    //   (BsaI on GAGACC cuts 5' of it). Differential over all 153 single-sided distal enzymes of
    //   the Biopython 1.88 snapshot: 0 mismatches (review 2026-09 B08 F54).
    //
    // • Two-sided (Type IIB) cutters carry a second pair in RestrictionEnzyme.SecondCut
    //   (= scd5, size + scd3) and are reported as one record per cut (RestrictionSite.IsSecondCut),
    //   with Biopython TwoCuts._modify / _rev_modify (`location + scd5`, `location − scd3`); each cut
    //   is filtered on its own on a linear molecule, as `_drop` filters the cut list element-wise.
    //   Differential over all 25 two-cut enzymes of the snapshot: 0 mismatches (B08 F56).
    //
    // • Readings vs. windows (B08 F60, revised by H7/H45 → F101). This finder reports every strand
    //   READING: a window that matches the site on both strands yields a forward AND a reverse
    //   record. The two declared oracles disagree on what a digest should then do, and EMBOSS is
    //   followed:
    //     – Biopython's `search` yields one cut set per WINDOW, the forward reading first:
    //       `Palindromic._search` applies `_modify` only, and `NonPalindromic._search` takes the
    //       forward branch whenever the forward named group of the `(?=fwd)|(?=rev)` compsite
    //       matched. A reverse reading is therefore never reported for a window that also reads
    //       forward, even when it cuts other bonds.
    //     – EMBOSS 6.6.0 `restrict` reports every READING whose cut differs:
    //       `embPatRestrictScan` (nucleus/embpat.c) runs two independent scan loops — forward over
    //       `substr`/`binstr`, then `forward = ajFalse` over `revstr`/`binrev` — and afterwards
    //       keeps one hit per distinct `cut1`, preferring the forward hit when two coincide.
    //   A palindromic recognition window physically admits both binding orientations, so for a
    //   distal cutter it is cleaved on either side; Biopython's one-reading-per-window rule is an
    //   optimisation that happens to be exact only while the two readings cut the same bonds.
    //   Digest, CreateMap, GetDigestSummary and PartialDigest therefore dedupe by CUT
    //   (`ReverseReadingCutsAlike`), not by window.
    //   Scope, computed over the 754 REBASE enzymes with a defined cut: 564 have a palindromic
    //   recognition STRING, of which 9 have `fst5 + fst3 != 0` (AlfI, BdaI, BplI, FalI, HauII,
    //   NmeDI, RdeGBIII, SgeI, UcoMSI); of those 9 exactly 2 (HauII TGGCCA, SgeI CNNG) have a cut
    //   SET that differs under the reverse reading — the other 7 are two-sided (Type IIB) and
    //   their two cut pairs swap. Across the whole set 166 enzymes have a differing cut set, and
    //   exactly 5 of those admit a window matching on both strands: HauII, SgeI, and the
    //   degenerate non-palindromic sites MspJI CNNR, LpnPI CCDG, SgrTI CCDS. Only those 5 change
    //   under this rule. Differential vs EMBOSS 6.6.0 on those 5 enzymes: 44 of 170 rows diverged
    //   under the Biopython rule (always a cut missing, never an extra), 0 after (F101); the
    //   749 other enzymes keep the F59 agreement with Biopython (0 mismatches by name).
    //
    // • Feasibility (linear molecules): a cut can only fall on a phosphodiester bond that
    //   exists, so a site is reported only when BOTH its strand cuts lie in [1, length - 1].
    //   This is Biopython's `_drop` filter (`1 < cut_on_watson <= length` and the same for
    //   `cut_on_crick = cut_on_watson - ovhg`, with `ovhg = CutPositionForward -
    //   CutPositionReverse`) expressed in 0-based coordinates. For an in-site cutter it only
    //   ever excludes an enzyme that cuts at the very edge of its own recognition window — in
    //   the built-in table, Sau3AI and MboI (^GATC / GATC^); a distal cutter loses every site
    //   whose cut would fall beyond either end of the molecule (FokI needs 13 more bases after
    //   GGATG), exactly as Biopython drops it. One documented divergence: TscAI/TspRI CASTG
    //   (offsets 7/-2) are stored in Biopython 1.88 with `ovhg = 10`, contradicting their own
    //   fst5/fst3 and 9-nt `ovhgseq` NNCASTGNN, so Biopython also drops a site whose
    //   bottom-strand cut is the (existing) bond after base 0; this class keeps it.
    //
    // • Circular molecules: the search runs over `seq + seq[0 .. L-2]` so sites spanning the
    //   origin are found (Biopython appends `data[1:size]`), cut positions are reduced into
    //   [0, length) with a non-negative modulo (a distal offset may be negative or exceed the
    //   length), and no site is dropped — on a closed circle every bond exists. Biopython's
    //   `_drop` folds a circular cut back by ONE length only and does not re-sort, so on a
    //   circle shorter than the cut offset it can report positions outside [1, length].
    // ───────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Finds all restriction sites for a specific enzyme in a <b>linear</b> sequence.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="enzymeName">Name of the restriction enzyme (case-insensitive).</param>
    /// <returns>Collection of restriction sites found: every forward-strand hit first, in ascending
    /// <see cref="RestrictionSite.Position"/> order, then every reverse-strand hit. The reverse-strand
    /// hits come out in <b>descending</b> Position order on a linear molecule, because the scan runs
    /// 5'→3' along the reverse complement and <c>Position = n - i - len</c> falls as the scan index
    /// rises; on a circular molecule they are not monotone at all, since the position is additionally
    /// reduced modulo the length. Measured on <c>GAATTCAAGAATTCAAGAATTC</c> with EcoRI: positions
    /// 0, 8, 16 forward then 16, 8, 0 reverse. Sort by <c>Position</c> if a global order is needed
    /// (review 2026-09 B08 H31; <see cref="FindOffTargets"/> in <c>CrisprDesigner</c> behaves the
    /// same way).</returns>
    /// <remarks>
    /// Equivalent to <see cref="FindSites(DnaSequence, string, MoleculeTopology)"/> with
    /// <see cref="MoleculeTopology.Linear"/>; see that overload for the coordinate convention.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null, or
    /// <paramref name="enzymeName"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is empty or not a known enzyme (<see cref="GetEnzyme"/>).</exception>
    public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, string enzymeName)
        => FindSites(sequence, enzymeName, MoleculeTopology.Linear);

    /// <summary>
    /// Finds all restriction sites for a specific enzyme in a molecule of the given topology.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="enzymeName">Name of the restriction enzyme (case-insensitive).</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <remarks>
    /// <para>
    /// <see cref="RestrictionSite.Position"/> and <see cref="RestrictionSite.CutPosition"/> are
    /// 0-based forward-strand coordinates for <b>both</b> strands; <c>CutPosition</c> is the
    /// number of forward-strand bases 5' of the cut made on the strand named by
    /// <see cref="RestrictionSite.IsForwardStrand"/>. A palindromic enzyme therefore reports two
    /// records per physical site — a forward one carrying the top-strand cut and a reverse one
    /// carrying the bottom-strand cut — at the same <c>Position</c>. Overlapping sites are all
    /// reported (Biopython's site regexes are lookaheads, so it reports them too).
    /// </para>
    /// <para>
    /// <b>Linear</b> (the default): a site is reported only when both of its strand cuts fall on
    /// a phosphodiester bond that exists, i.e. both lie in <c>[1, Length - 1]</c>. This matches
    /// Biopython 1.88 <c>Bio.Restriction</c> (<c>AbstractCut._drop</c>), which drops a site when
    /// either cut would fall outside the molecule. It affects only enzymes cutting at the edge of
    /// their own recognition window: <c>Sau3AI</c>/<c>MboI</c> (<c>^GATC</c>) at the very 5' or
    /// 3' end. Example: <c>FindSites(new DnaSequence("GATCGATCGATC"), "Sau3AI")</c> reports the
    /// single surviving <i>window</i>, the one at position 4 — as two records, the forward one with
    /// cut 4 and the reverse one with cut 8, per the two-records-per-palindromic-site rule above.
    /// Its top-strand cut is 4, exactly
    /// <c>Bio.Restriction.Sau3AI.search(Seq("GATCGATCGATC")) == [5]</c> (1-based); the windows at 0
    /// and 8 are dropped because one of their cuts would fall outside the molecule
    /// (review 2026-09 B08 H31 — this used to read "the single site", which undercounted the
    /// records by one).
    /// </para>
    /// <para>
    /// <b>Circular</b>: the recognition window may wrap past the origin, so sites starting in the
    /// last <c>L-1</c> bases are found as well, and cut positions are reduced modulo the sequence
    /// length (a cut "after the last base" is the same bond as a cut at the origin). Nothing is
    /// dropped. This is Biopython's <c>search(..., linear=False)</c>, which searches
    /// <c>data + data[1:size]</c> and folds out-of-range cuts back by ±length. Example:
    /// <c>Bio.Restriction.EcoRI.search(Seq("TTCAAAAAAAGAA"), linear=False) == [12]</c>, i.e. the
    /// origin-spanning <c>GAATTC</c> at position 10 with a top-strand cut at 11.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, string enzymeName, MoleculeTopology topology)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        var enzyme = RequireEnzyme(enzymeName);
        return FindSitesCore(sequence.Sequence, enzyme, topology);
    }

    /// <summary>
    /// Finds all restriction sites for a specific enzyme in a raw <b>linear</b> sequence string.
    /// </summary>
    /// <remarks>
    /// The sequence is upper-cased before matching; characters that are not IUPAC nucleotide
    /// codes simply never match; an IUPAC ambiguity code in the sequence (e.g. <c>N</c>) is
    /// matched by a site <c>N</c> only, exactly as in Biopython (see
    /// <see cref="FindPatternOccurrences"/>). A null or empty sequence yields no sites — but an unknown
    /// enzyme name is still rejected, eagerly, before the sequence is looked at.
    /// </remarks>
    public static IEnumerable<RestrictionSite> FindSites(string sequence, string enzymeName)
        => FindSites(sequence, enzymeName, MoleculeTopology.Linear);

    /// <summary>
    /// Finds all restriction sites for a specific enzyme in a raw sequence string, for a molecule
    /// of the given topology.
    /// </summary>
    /// <remarks>
    /// The sequence is upper-cased before matching; characters that are not IUPAC nucleotide
    /// codes simply never match; an IUPAC ambiguity code in the sequence (e.g. <c>N</c>) is
    /// matched by a site <c>N</c> only, exactly as in Biopython —
    /// <c>FindSites("GGCCNNNNNGGCC", "SfiI")</c> reports the site at 0 (top-strand cut 8) where
    /// <c>SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]</c>. A null or empty sequence yields no sites — but an unknown
    /// enzyme name is still rejected, eagerly, before the sequence is looked at. See
    /// <see cref="FindSites(DnaSequence, string, MoleculeTopology)"/> for the coordinate and
    /// feasibility conventions.
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindSites(
        string sequence, string enzymeName, MoleculeTopology topology)
    {
        RequireTopology(topology);
        var enzyme = RequireEnzyme(enzymeName);

        // ASCII-only fold (review 2026-09 B08 F49, superseding F43): ToUpperInvariant would map U+017F 'ſ'
        // onto the IUPAC code 'S', so a non-IUPAC character would match.
        return string.IsNullOrEmpty(sequence)
            ? []
            : FindSitesCore(SequenceExtensions.ToUpperAscii(sequence), enzyme, topology);
    }

    /// <summary>
    /// Finds all restriction sites for a custom enzyme in a <b>linear</b> sequence.
    /// </summary>
    /// <exception cref="ArgumentException">
    /// <paramref name="enzyme"/> has an empty recognition sequence, or one containing a character
    /// that is not one of the 15 IUPAC nucleotide codes, a second cut pair that leaves a different stagger,
    /// two strand cuts more than <see cref="MaxCutStagger"/> bases apart, a recognition sequence longer than <see cref="MaxRecognitionLength"/>
    /// or a name longer than <see cref="MaxEnzymeNameLength"/>.
    /// </exception>
    public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, RestrictionEnzyme enzyme)
        => FindSites(sequence, enzyme, MoleculeTopology.Linear);

    /// <summary>
    /// Finds all restriction sites for a custom enzyme in a molecule of the given topology.
    /// </summary>
    /// <remarks>
    /// The recognition sequence is matched case-insensitively with IUPAC ambiguity semantics
    /// (<see cref="IupacHelper.MatchesIupac"/>) and is validated up front. Unlike the built-in
    /// table, a custom enzyme need not be palindromic: a non-palindromic enzyme yields distinct
    /// forward- and reverse-strand sites, and the reverse-strand record reports the cut the enzyme
    /// makes on the strand it read — see
    /// <see cref="FindSites(DnaSequence, string, MoleculeTopology)"/> for the full convention.
    /// </remarks>
    /// <exception cref="ArgumentException">
    /// <paramref name="enzyme"/> has an empty recognition sequence, or one containing a character
    /// that is not one of the 15 IUPAC nucleotide codes, a second cut pair that leaves a different stagger,
    /// two strand cuts more than <see cref="MaxCutStagger"/> bases apart, a recognition sequence longer than <see cref="MaxRecognitionLength"/>
    /// or a name longer than <see cref="MaxEnzymeNameLength"/>.
    /// </exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, RestrictionEnzyme enzyme, MoleculeTopology topology)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzyme);

        return FindSitesCore(sequence.Sequence, enzyme, topology);
    }

    /// <summary>
    /// Finds all restriction sites for multiple enzymes in a <b>linear</b> sequence.
    /// </summary>
    /// <remarks>The result is the concatenation of the per-enzyme searches, in the order given;
    /// a name repeated in <paramref name="enzymeNames"/> is searched (and reported) twice. The names
    /// are resolved eagerly: a null or empty array, and an unknown or empty name, are reported at the
    /// call (review 2026-09 B08 H24).</remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or
    /// <paramref name="enzymeNames"/> is null, or a name in it is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeNames"/> is empty, or a name in it
    /// is empty or not a REBASE enzyme with a defined cut.</exception>
    public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindSitesCore(sequence, enzymeNames, MoleculeTopology.Linear);
    }

    /// <summary>
    /// Finds all restriction sites for multiple enzymes in a molecule of the given topology.
    /// </summary>
    /// <remarks>The result is the concatenation of the per-enzyme searches, in the order given;
    /// a name repeated in <paramref name="enzymeNames"/> is searched (and reported) twice. The names
    /// are resolved eagerly, and an empty array is an error — the same contract as the methylation
    /// overload <see cref="FindSites(DnaSequence, MoleculeTopology, IReadOnlyList{Methylase}, string[])"/>,
    /// as <see cref="Digest(DnaSequence, MoleculeTopology, string[])"/>, <c>PartialDigest</c> and
    /// <c>GetDigestSummary</c>, and as the enzyme-record overloads. The one member where an empty
    /// array means something else is <see cref="CreateMap(DnaSequence, string[])"/>, whose parameter
    /// documents it as "every common enzyme" (review 2026-09 B08 H24).</remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or
    /// <paramref name="enzymeNames"/> is null, or a name in it is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeNames"/> is empty, or a name in it
    /// is empty or not a REBASE enzyme with a defined cut.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, MoleculeTopology topology, params string[] enzymeNames)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        return FindSitesCore(sequence, enzymeNames, topology);
    }

    /// <summary>
    /// Finds all restriction sites for several caller-supplied enzyme definitions — including
    /// enzymes that are not in the built-in table (Type IIS / Type IIB cutters) — in a molecule of
    /// the given topology.
    /// </summary>
    /// <remarks>The result is the concatenation of the per-enzyme searches
    /// (<see cref="FindSites(DnaSequence, RestrictionEnzyme, MoleculeTopology)"/>), in the order
    /// given; an enzyme listed twice is searched (and reported) twice, as with the name overload.
    /// The collection is validated eagerly.</remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or
    /// <paramref name="enzymes"/> is null, or <paramref name="enzymes"/> contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, or an enzyme's
    /// recognition sequence is empty or not IUPAC, its second cut pair is inconsistent, or the enzyme is larger than <see cref="MaxCutStagger"/> / <see cref="MaxRecognitionLength"/> / <see cref="MaxEnzymeNameLength"/> allow.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        var list = RequireEnzymeList(enzymes);
        return list.SelectMany(e => FindSitesCore(sequence.Sequence, e, topology));
    }

    /// <remarks>
    /// The names are resolved before anything is enumerated, so a misspelled name is reported at the
    /// call and not on first enumeration — the eager rule <see cref="RequireEnzyme"/> states and the
    /// raw-string and compatibility members already followed (review 2026-09 B08 H24).
    /// </remarks>
    private static IEnumerable<RestrictionSite> FindSitesCore(
        DnaSequence sequence, string[] enzymeNames, MoleculeTopology topology)
    {
        var enzymes = RequireEnzymeNames(enzymeNames, nameof(enzymeNames));
        return enzymes.SelectMany(e => FindSitesCore(sequence.Sequence, e, topology));
    }

    /// <summary>
    /// Refuses a <see cref="MoleculeTopology"/> value that is neither <see cref="MoleculeTopology.Linear"/> nor
    /// <see cref="MoleculeTopology.Circular"/>. Before B08 round 5 (A5R-5) such a value (reachable over MCP, where an
    /// enum deserializes from any number) meant linear in <c>FindSites</c> / <c>CreateMap</c> and circular in
    /// <c>Digest</c> / <c>GetDigestSummary</c> / <c>PartialDigest</c>, and was echoed back by the summary.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined value.</exception>
    private static void RequireTopology(MoleculeTopology topology)
    {
        if (!Enum.IsDefined(topology))
            throw new ArgumentOutOfRangeException(nameof(topology), topology,
                "The topology must be MoleculeTopology.Linear or MoleculeTopology.Circular.");
    }

    /// <summary>
    /// Resolves a non-empty array of enzyme names eagerly — the single contract for the
    /// <c>params string[]</c> shape: null array, empty array, null/empty name and unknown name are
    /// all reported at the call (review 2026-09 B08 H24). <see cref="CreateMap(DnaSequence, string[])"/>
    /// is the one deliberate exception: there an empty array means "every common enzyme", which is
    /// that parameter's documented meaning, so it does not come through here.
    /// </summary>
    private static List<RestrictionEnzyme> RequireEnzymeNames(string[] enzymeNames, string paramName)
    {
        if (enzymeNames is null)
            throw new ArgumentNullException(paramName);
        if (enzymeNames.Length == 0)
            throw new ArgumentException("At least one enzyme name is required", paramName);

        var resolved = new List<RestrictionEnzyme>(enzymeNames.Length);
        foreach (string name in enzymeNames)
            resolved.Add(RequireEnzyme(name, paramName));

        return resolved;
    }

    /// <summary>
    /// Finds all sites of the 39 common enzymes (<see cref="Enzymes"/>) in a <b>linear</b> sequence.
    /// For the full REBASE set use <see cref="FindAllSites(DnaSequence, MoleculeTopology, EnzymeSet)"/>.
    /// </summary>
    public static IEnumerable<RestrictionSite> FindAllSites(DnaSequence sequence)
        => FindAllSites(sequence, MoleculeTopology.Linear);

    /// <summary>
    /// Finds all sites of the 39 common enzymes (<see cref="Enzymes"/>) in a molecule of the given
    /// topology. For the full REBASE set use <see cref="FindAllSites(DnaSequence, MoleculeTopology, EnzymeSet)"/>.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<RestrictionSite> FindAllSites(DnaSequence sequence, MoleculeTopology topology)
        => FindAllSites(sequence, topology, EnzymeSet.Common);

    /// <summary>
    /// Enumerates every 0-based start index at which <paramref name="recognitionPattern"/> occurs
    /// in <paramref name="sequence"/>, under IUPAC ambiguity semantics, on the given strand only.
    /// </summary>
    /// <param name="sequence">Sequence to scan. Upper-cased before matching, so the scan is
    /// case-insensitive; a null or empty sequence yields no occurrences.</param>
    /// <param name="recognitionPattern">Recognition pattern, 1 or more of the 15 IUPAC nucleotide
    /// codes A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V (case-insensitive).</param>
    /// <param name="startIndex">First index to consider; occurrences starting before it are
    /// skipped. Defaults to 0.</param>
    /// <returns>The start indices in ascending order, lazily. Indices are 0-based and point at
    /// the 5'-most base of the match; <b>overlapping</b> occurrences are all reported.</returns>
    /// <remarks>
    /// <para>
    /// This is the raw pattern scan underlying <see cref="FindSites(DnaSequence, string)"/>,
    /// exposed for callers that hold a pattern rather than an enzyme and a sequence that
    /// <see cref="DnaSequence"/> would reject (e.g. one still carrying IUPAC ambiguity codes).
    /// It reports occurrences only — it applies no enzyme, no cut arithmetic, no reverse-strand
    /// scan and no linear-molecule cut-feasibility filter. To cover both strands, call it a
    /// second time with <see cref="DnaSequence.GetReverseComplementString"/> of the pattern.
    /// </para>
    /// <para>
    /// Matching is per position with Biopython 1.88 <c>Bio.Restriction</c>'s site-regex rule,
    /// composed from the canonical <see cref="IupacHelper"/>: a pattern <c>N</c> matches any of
    /// the 15 IUPAC codes in the sequence — an <c>N</c> included, so
    /// <c>FindPatternOccurrences("GGCCNNNNNGGCC", "GGCCNNNNNGGCC")</c> yields 0 exactly as
    /// <c>SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]</c>; every other pattern code matches only a
    /// concrete A/C/G/T base it expands to, so an ambiguity code in the sequence (<c>N</c>,
    /// <c>R</c>, ...) blocks it rather than widening it (<c>"GAANTC"</c> does not match
    /// <c>GAATTC</c>, nor <c>"GTYRAC"</c> match <c>GTYRAC</c>). A character outside the 15 IUPAC
    /// codes matches nothing.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentException"><paramref name="recognitionPattern"/> is null or empty,
    /// or contains a character that is not an IUPAC nucleotide code.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="startIndex"/> is negative.</exception>
    public static IEnumerable<int> FindPatternOccurrences(
        string? sequence, string recognitionPattern, int startIndex = 0)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(startIndex);
        string pattern = ValidatedPattern(recognitionPattern, nameof(recognitionPattern));

        // ASCII-only fold (review 2026-09 B08 F49, superseding F43): ToUpperInvariant would map U+017F 'ſ'
        // onto the IUPAC code 'S', so a non-IUPAC character would match.
        // A null or empty sequence yields no occurrences, as the <param> doc states and as the
        // sibling raw-string FindSites(string, string, MoleculeTopology) already did; the pattern is
        // still validated eagerly, before the sequence is looked at (review 2026-09 B08 H25).
        return string.IsNullOrEmpty(sequence)
            ? []
            : EnumeratePatternOccurrences(SequenceExtensions.ToUpperAscii(sequence), pattern, startIndex);
    }

    private static IEnumerable<int> EnumeratePatternOccurrences(string seq, string pattern, int startIndex)
    {
        // i <= Length - pattern.Length, not i + pattern.Length <= Length: the sum overflows for a start index near
        // int.MaxValue and the loop would then read out of range (B08 round 5 A5R-2).
        for (int i = startIndex; i <= seq.Length - pattern.Length; i++)
        {
            if (IupacPattern.MatchesAt(seq, i, pattern))
                yield return i;
        }
    }

    /// <summary>
    /// Resolves an enzyme name (<see cref="GetEnzyme"/>: the full REBASE set) or throws — the single name-resolution contract of this
    /// class, shared by the site finder, the digest, the map and the compatibility members.
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive; surrounding whitespace and a REBASE
    /// hyphen are tolerated, see <see cref="GetEnzyme"/>).</param>
    /// <param name="paramName">Name of the caller's parameter to blame, when it is not
    /// <c>enzymeName</c>.</param>
    /// <remarks>
    /// The two failures, which every name-taking public member documents with these exact types
    /// (review 2026-09 B08 H24, H19): a <see langword="null"/> name is
    /// <see cref="ArgumentNullException"/>, and an <i>empty</i> name or one that is not a REBASE enzyme
    /// with a defined cut is <see cref="ArgumentException"/> - the null-versus-empty split of
    /// <c>ArgumentException.ThrowIfNullOrEmpty</c>. Both are eager — no name-taking member defers them
    /// to first enumeration.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is empty or not a REBASE enzyme with
    /// a defined cut.</exception>
    private static RestrictionEnzyme RequireEnzyme(string enzymeName, string? paramName = null)
    {
        paramName ??= nameof(enzymeName);

        if (enzymeName is null)
            throw new ArgumentNullException(paramName);
        if (enzymeName.Length == 0)
            throw new ArgumentException("The enzyme name must not be empty.", paramName);

        return GetEnzyme(enzymeName)
            ?? throw new ArgumentException($"Unknown enzyme: {enzymeName}", paramName);
    }

    /// <summary>
    /// Validates a recognition pattern and returns it upper-cased, ready for matching.
    /// </summary>
    private static string ValidatedPattern(string? pattern, string paramName, string? enzymeName = null)
    {
        string what = enzymeName is null ? "Recognition pattern" : $"Enzyme '{enzymeName}' recognition sequence";

        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException($"{what} cannot be null or empty.", paramName);

        // ASCII-only fold (review 2026-09 B08 F49, superseding F43): ToUpperInvariant would map U+017F 'ſ'
        // onto the IUPAC code 'S', so a non-IUPAC character would match.
        string upper = SequenceExtensions.ToUpperAscii(pattern);
        for (int i = 0; i < upper.Length; i++)
        {
            if (!IupacHelper.IsNucleotideCode(upper[i]))
            {
                throw new ArgumentException(
                    $"{what} '{pattern}' contains '{pattern[i]}' at index {i}, which is not one of the " +
                    "15 IUPAC nucleotide codes (A, C, G, T, N, R, Y, S, W, K, M, B, D, H, V).",
                    paramName);
            }
        }

        return upper;
    }

    /// <summary>
    /// Dispatches a single-enzyme search; validates the recognition pattern eagerly so a malformed
    /// custom enzyme is reported at the call, not on first enumeration.
    /// </summary>
    private static IEnumerable<RestrictionSite> FindSitesCore(
        string seq, RestrictionEnzyme enzyme, MoleculeTopology topology, bool[]? blocked = null)
    {
        ValidateEnzymeLimits(enzyme, nameof(enzyme));
        string pattern = ValidatedPattern(enzyme.RecognitionSequence, nameof(enzyme), enzyme.Name);

        // A linear molecule shorter than the recognition sequence holds no window at all, but a
        // circular one does: the window reads past the origin and back into the molecule. The
        // length guard is therefore scoped to linear topology; on a circle only an empty molecule
        // has no window (review 2026-09 B08 H6 — the guard used to run before the dispatch, so
        // `FindSites("ACGT", "BsiWI", Circular)` answered no sites where Biopython cuts at 2).
        if (topology == MoleculeTopology.Circular)
            return seq.Length == 0 ? [] : FindSitesCircular(seq, enzyme, pattern, blocked);

        return seq.Length < pattern.Length ? [] : FindSitesLinear(seq, enzyme, pattern, blocked);
    }

    /// <summary>
    /// Rejects a caller-supplied enzyme that is larger than this library works with: a name longer than
    /// <see cref="MaxEnzymeNameLength"/>, a recognition sequence longer than <see cref="MaxRecognitionLength"/>, two strand cuts more
    /// than <see cref="MaxCutStagger"/> bases apart (the stagger is computed in <c>long</c>: <c>int</c> subtraction wrapped for offsets
    /// near the limits and reported a 1-nt overhang for a 4.3-billion-nt one), and then a <see cref="RestrictionEnzyme.SecondCut"/>
    /// pair whose stagger differs from the first pair's (Biopython stores one <c>ovhg</c> per enzyme; see the property's remarks).
    /// The size limits come first so that no work is spent on an enzyme that is refused anyway.
    /// </summary>
    private static void ValidateEnzymeLimits(RestrictionEnzyme enzyme, string paramName)
    {
        if (enzyme.Name is { Length: > MaxEnzymeNameLength })
        {
            throw new ArgumentException(
                $"The enzyme name is {enzyme.Name.Length} characters long; the longest supported name is {MaxEnzymeNameLength}.",
                paramName);
        }

        if (enzyme.RecognitionSequence is { Length: > MaxRecognitionLength })
        {
            throw new ArgumentException(
                $"Enzyme '{enzyme.Name}' has a recognition sequence of {enzyme.RecognitionSequence.Length} bases; " +
                $"the longest supported is {MaxRecognitionLength}.",
                paramName);
        }

        long firstStagger = (long)enzyme.CutPositionForward - enzyme.CutPositionReverse;
        if (Math.Abs(firstStagger) > MaxCutStagger)
        {
            throw new ArgumentException(
                $"Enzyme '{enzyme.Name}' cuts its two strands {Math.Abs(firstStagger)} bases apart " +
                $"({enzyme.CutPositionForward}, {enzyme.CutPositionReverse}); the largest supported stagger is {MaxCutStagger}.",
                paramName);
        }

        if (enzyme.SecondCut is { } second &&
            (long)second.Forward - second.Reverse != firstStagger)
        {
            throw new ArgumentException(
                $"Enzyme '{enzyme.Name}' second cut pair ({second.Forward}, {second.Reverse}) leaves a " +
                $"different stagger from its first pair ({enzyme.CutPositionForward}, " +
                $"{enzyme.CutPositionReverse}); both cuts of a two-sided enzyme must leave the same end.",
                paramName);
        }
    }

    /// <summary>
    /// The enzyme's cut pairs: the first, and the <see cref="RestrictionEnzyme.SecondCut"/> pair of a
    /// two-sided (Type IIB) cutter, each flagged with whether it is the second.
    /// </summary>
    private static IEnumerable<(int Forward, int Reverse, bool IsSecond)> CutPairs(RestrictionEnzyme enzyme)
    {
        yield return (enzyme.CutPositionForward, enzyme.CutPositionReverse, false);
        if (enzyme.SecondCut is { } second)
            yield return (second.Forward, second.Reverse, true);
    }

    private static IEnumerable<RestrictionSite> FindSitesLinear(
        string seq, RestrictionEnzyme enzyme, string pattern, bool[]? blocked)
    {
        int n = seq.Length;
        int len = pattern.Length;
        var pairs = CutPairs(enzyme).ToArray();

        for (int i = 0; i + len <= n; i++)
        {
            if (!IupacPattern.MatchesAt(seq, i, pattern) || MethylationMask.IsBlocked(blocked, i, len))
                continue;

            // Forward-strand site: top-strand cut i+fwd, bottom-strand cut i+rev — once per cut pair
            // (a Type IIB enzyme cuts on both sides; Biopython TwoCuts._modify). Each cut is kept or
            // dropped on its own, as Biopython's _drop filters the cut list element by element.
            foreach (var (fwd, rev, second) in pairs)
            {
                if (!CutsAreInternal(i + fwd, i + rev, n))
                    continue;

                yield return new RestrictionSite(
                    Position: i,
                    Enzyme: enzyme,
                    IsForwardStrand: true,
                    CutPosition: i + fwd,
                    RecognizedSequence: seq.Substring(i, len)) { IsSecondCut = second };
            }
        }

        // Check reverse strand (complement sequence)
        string revComp = IupacPattern.ReverseComplement(seq);
        for (int i = 0; i + len <= n; i++)
        {
            if (!IupacPattern.MatchesAt(revComp, i, pattern) || MethylationMask.IsBlocked(blocked, n - i - len, len))
                continue;

            // Convert back to forward strand coordinates
            int forwardPos = n - i - len;

            foreach (var (fwd, rev, second) in pairs)
            {
                // The enzyme reads the reverse strand, so its own 5' end sits at forward
                // coordinate forwardPos+len: top-strand cut forwardPos+len-rev, bottom-strand cut
                // forwardPos+len-fwd. For a palindromic enzyme fwd+rev == len, so the bottom-strand
                // cut equals the historical forwardPos+rev.
                int watsonCut = forwardPos + len - rev;
                int crickCut = forwardPos + len - fwd;

                if (!CutsAreInternal(watsonCut, crickCut, n))
                    continue;

                yield return new RestrictionSite(
                    Position: forwardPos,
                    Enzyme: enzyme,
                    IsForwardStrand: false,
                    CutPosition: crickCut,
                    RecognizedSequence: revComp.Substring(i, len)) { IsSecondCut = second };
            }
        }
    }

    private static IEnumerable<RestrictionSite> FindSitesCircular(
        string seq, RestrictionEnzyme enzyme, string pattern, bool[]? blocked)
    {
        int n = seq.Length;
        int len = pattern.Length;
        var pairs = CutPairs(enzyme).ToArray();

        // Biopython searches data + data[1:size] so a site may wrap past the origin; a start
        // index is therefore considered exactly once, for 0 <= i < n. A distal (Type IIS) cut
        // offset can be negative or exceed the length, so every cut is reduced with a
        // non-negative modulo into [0, n) (review 2026-09 B08 F54).
        string extended = WrapExtend(seq, n + len - 1);
        for (int i = 0; i < n; i++)
        {
            if (!IupacPattern.MatchesAt(extended, i, pattern) || MethylationMask.IsBlocked(blocked, i, len))
                continue;

            foreach (var (fwd, _, second) in pairs)
            {
                yield return new RestrictionSite(
                    Position: i,
                    Enzyme: enzyme,
                    IsForwardStrand: true,
                    CutPosition: Modulo((long)i + fwd, n),
                    RecognizedSequence: extended.Substring(i, len)) { IsSecondCut = second };
            }
        }

        string revComp = IupacPattern.ReverseComplement(seq);
        string revExtended = WrapExtend(revComp, n + len - 1);
        for (int i = 0; i < n; i++)
        {
            if (!IupacPattern.MatchesAt(revExtended, i, pattern) || MethylationMask.IsBlocked(blocked, Modulo(n - i - len, n), len))
                continue;

            int forwardPos = ((n - i - len) % n + n) % n;

            foreach (var (fwd, _, second) in pairs)
            {
                yield return new RestrictionSite(
                    Position: forwardPos,
                    Enzyme: enzyme,
                    IsForwardStrand: false,
                    CutPosition: Modulo((long)forwardPos + len - fwd, n),
                    RecognizedSequence: revExtended.Substring(i, len)) { IsSecondCut = second };
            }
        }
    }

    /// <summary>
    /// The molecule written out to <paramref name="totalLength"/> bases by repeating it, so that
    /// every start index in <c>[0, seq.Length)</c> has a full recognition-length window available
    /// past the origin.
    /// </summary>
    /// <remarks>
    /// For a molecule at least <c>len - 1</c> bases long this is exactly Biopython's circular
    /// extension <c>data + data[1:size]</c> (<c>FormattedSeq.finditer</c>, whose <c>data</c> carries
    /// a leading space, so the slice is the molecule's first <c>size - 1</c> bases). Biopython's
    /// slice silently truncates on a shorter molecule, which drops the windows that need more than
    /// one turn around the circle; repeating instead keeps all <c>n</c> start indices searchable
    /// (review 2026-09 B08 H6). The divergence is reported under
    /// <c>docs/Validation/review-2026-09/B08-restriction-slice.md</c> F100: on circular <c>ACGT</c>
    /// Biopython finds <c>BsiWI</c> (cut 2) but not <c>SnaBI</c>, whose window <c>TACGTA</c> at
    /// start 3 is physically present on the circle.
    /// </remarks>
    private static string WrapExtend(string seq, int totalLength)
    {
        if (seq.Length >= totalLength)
            return seq[..totalLength];

        return string.Create(totalLength, seq, static (span, molecule) =>
        {
            for (int written = 0; written < span.Length;)
            {
                int take = Math.Min(molecule.Length, span.Length - written);
                molecule.AsSpan(0, take).CopyTo(span[written..]);
                written += take;
            }
        });
    }

    /// <summary>
    /// The non-negative remainder of <paramref name="value"/> modulo <paramref name="n"/>. The value is a <c>long</c>:
    /// a caller-supplied cut offset of <c>int.MaxValue</c> added to a position overflows <c>int</c> and used to
    /// land on the wrong base of a circular molecule (B08 round 5 A5R-7).
    /// </summary>
    private static int Modulo(long value, int n) => (int)(((value % n) + n) % n);

    /// <summary>
    /// True when both strand cuts of a site fall on a phosphodiester bond that exists in a linear
    /// molecule of <paramref name="length"/> bases, i.e. both lie in <c>[1, length - 1]</c>
    /// (Biopython <c>AbstractCut._drop</c>, expressed in 0-based coordinates).
    /// </summary>
    private static bool CutsAreInternal(int watsonCut, int crickCut, int length) =>
        watsonCut >= 1 && watsonCut <= length - 1 &&
        crickCut >= 1 && crickCut <= length - 1;

    #endregion

    #region Digest Simulation

    // ───────────────────────────────────────────────────────────────────────────────────────
    // Digest conventions (reference: Biopython 1.88 Bio.Restriction, Restriction.py —
    // AbstractCut.catalyse/catalyze L571-610, RestrictionType.search L1180-1210,
    // AbstractCut.ovhgseq / is_5overhang / is_3overhang / is_blunt L480-540).
    //
    // • Fragment boundaries are TOP-STRAND (watson) cut positions only. Biopython's
    //   `enzyme.catalyse(seq)` slices the molecule at `enzyme.search(seq)`, which is the list of
    //   top-strand cut sites, so a digest is a partition of the top strand even for a staggered
    //   (sticky-end) cut. A forward-strand site record contributes its `CutPosition`, the
    //   top-strand cut. A reverse-strand record whose window also reads as the site on the forward
    //   strand — every reverse record of a palindromic enzyme, the 39 built-ins included — is the
    //   same window as a forward record and is dropped (one cut per window, the forward one:
    //   Biopython `Palindromic._search` searches one strand, `NonPalindromic._search`'s regex
    //   alternation takes the forward branch). Every other reverse-strand record is a distinct
    //   site of a non-palindromic enzyme (e.g. a Type IIS GAGACC for BsaI GGTCTC): it contributes
    //   its own top-strand cut `Position + L − CutPositionReverse` (Biopython `_rev_modify`,
    //   `location − fst3`), which lies on the OTHER side of the site (review 2026-09 B08 F55).
    //   The bottom-strand cut is not lost — it is reported as the fragment end's overhang (below).
    //
    // • Fragment end chemistry. A cut by an enzyme with forward cut `f` and reverse cut `r` has
    //   its top-strand cut at `watson` and its bottom-strand cut at `crick = watson + (r - f)` —
    //   for a site on either strand, and whether or not the cut lies inside the site (a distal
    //   cutter's window is bases of the flanking molecule, e.g. BsaI's 4 nt after GGTCTCN).
    //   The single-stranded extension spans the top-strand window
    //   `[min(watson, crick), max(watson, crick))`, of length `|r - f|`; `r > f` is a 5' overhang,
    //   `r < f` a 3' overhang, `r == f` a blunt end. The overhang is reported per fragment end,
    //   written 5'→3' on the strand that carries it: for a 5' overhang the downstream fragment's
    //   left end carries it on the top strand (the window verbatim) and the upstream fragment's
    //   right end carries the complementary extension on the bottom strand (the window's reverse
    //   complement); for a 3' overhang the two swap. Enzyme-level cross-check: for every one of
    //   the 39 built-in enzymes this window taken over the recognition sequence equals Biopython's
    //   `enzyme.ovhgseq` and the type equals `is_5overhang()`/`is_3overhang()`/`is_blunt()`
    //   (39/39, verified 2026-10-14 against Biopython 1.88) — e.g. EcoRI AATT 5', PstI TGCA 3',
    //   EcoRV "" blunt. Taken over the actual molecule rather than the pattern, a degenerate
    //   enzyme reports its real overhang (SfiI gives the three real bases, not Biopython's "NNN").
    //
    // • Fragment ORDER. Linear: 5'→3', the first fragment starting at position 0. Circular:
    //   ascending cut position, so the origin-spanning fragment is reported LAST (its
    //   `StartPosition` is the largest cut). Biopython emits the same multiset of fragments but
    //   puts the origin-spanning one FIRST (`fragments.append(d[r[-1]:] + d[1:r[0]])` before the
    //   rest); Seqeron's order is therefore Biopython's rotated by one, chosen so that
    //   `StartPosition` increases monotonically over the returned fragments in both topologies.
    // ───────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Simulates a restriction digest of a <b>linear</b> molecule and returns the resulting
    /// fragments. For circular-molecule (plasmid) support use the
    /// <see cref="Digest(DnaSequence, MoleculeTopology, string[])"/> overload.
    /// </summary>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="enzymeNames">Names of restriction enzymes to use.</param>
    /// <returns>Collection of DNA fragments after digestion, 5'→3'.</returns>
    /// <remarks>
    /// Equivalent to <see cref="Digest(DnaSequence, MoleculeTopology, string[])"/> with
    /// <see cref="MoleculeTopology.Linear"/>; see that overload for the full convention
    /// (top-strand boundaries, fragment order, and end overhangs).
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="enzymeNames"/> is null, or a name in it is null.</exception>
    public static IEnumerable<DigestFragment> Digest(DnaSequence sequence, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return DigestCore(sequence.Sequence, RequireEnzymeNames(enzymeNames, nameof(enzymeNames)));
    }

    /// <summary>
    /// Collects the distinct top-strand (forward-strand) cut positions for the requested enzymes,
    /// together with every enzyme's site at each position.
    /// </summary>
    /// <remarks>
    /// Each site is keyed by its top-strand cut, which is what Biopython's <c>catalyse</c> slices
    /// at: a forward-strand record's <see cref="RestrictionSite.CutPosition"/>, or a reverse-strand
    /// record's top-strand cut (see <c>TopStrandCut</c>). A reverse-strand record whose window also
    /// matches on the forward strand (the other record of a palindromic site) is skipped — it
    /// carries the bottom-strand cut of the same site and would double-count it as a second
    /// boundary; the reverse-strand sites of a non-palindromic enzyme are kept (F55). Several
    /// enzymes can cut at the same top-strand coordinate — isoschizomers (MspI/HpaII,
    /// MboI/Sau3AI) and nested sites (BamHI/BglII with MboI/Sau3AI, ClaI/TaqI, StuI/HaeIII,
    /// ScaI/RsaI, FseI/HaeIII; Biopython 1.88 <c>search</c> reports the same position for each
    /// pair; review 2026-09 B08 F48). All of them are kept, one site per enzyme (a name listed
    /// twice, in any casing, counts once), ordered by <see cref="StringComparer.Ordinal"/> on the
    /// table's enzyme name, so the result does not depend on the order of
    /// <paramref name="enzymeNames"/>. The first site of each list is the boundary's
    /// representative (see <see cref="DigestFragment.LeftEnzyme"/>).
    /// </remarks>
    private static (List<int> Cuts, Dictionary<int, List<RestrictionSite>> Sites) CollectForwardCuts(
        string seq, IEnumerable<RestrictionEnzyme> enzymes, MoleculeTopology topology, bool[]? blocked = null)
    {
        var cutPositions = new SortedSet<int>();
        var sitesByPosition = new Dictionary<int, List<RestrictionSite>>();

        foreach (var enzyme in enzymes)
        {
            string pattern = ValidatedPattern(enzyme.RecognitionSequence, nameof(enzyme), enzyme.Name);
            foreach (var site in FindSitesCore(seq, enzyme, topology, blocked))
            {
                // A reverse-strand record is dropped only when reading the window backwards
                // cuts the same bonds as reading it forwards, i.e. when it is the same cleavage
                // counted twice. That is EMBOSS 6.6.0's rule: embPatRestrictScan runs two
                // independent loops (forward over `substr`, then `forward = ajFalse` over
                // `revstr`, nucleus/embpat.c) and afterwards keeps one hit per distinct `cut1`,
                // preferring the forward one. Biopython drops one reading per *window* instead
                // (Palindromic._search applies `_modify` only; NonPalindromic._search takes the
                // forward branch whenever the forward group matched), which loses a physically
                // real cut for the five enzymes whose two readings cut differently and whose
                // window can match on both strands — HauII, SgeI, MspJI, LpnPI, SgrTI. EMBOSS is
                // followed here (review 2026-09 B08 H7; slice report F101).
                if (!site.IsForwardStrand
                    && IupacPattern.MatchesAtCircular(seq, site.Position, pattern)
                    && ReverseReadingCutsAlike(site.Enzyme))
                {
                    continue;
                }

                // FindSites already reports a circular cut reduced modulo the sequence length and
                // a linear site only when both strand cuts lie inside [1, Length-1]
                // (CutsAreInternal), so the top-strand cut needs no further normalisation.
                int watson = TopStrandCut(site, seq.Length, topology);
                cutPositions.Add(watson);
                if (!sitesByPosition.TryGetValue(watson, out var atCut))
                    sitesByPosition[watson] = atCut = [];
                if (!atCut.Exists(s => s.Enzyme.Name == site.Enzyme.Name))
                    atCut.Add(site);
            }
        }

        foreach (var atCut in sitesByPosition.Values)
            atCut.Sort((a, b) => string.CompareOrdinal(a.Enzyme.Name, b.Enzyme.Name));

        return (cutPositions.ToList(), sitesByPosition);
    }

    private static IEnumerable<DigestFragment> DigestCore(string seq, IEnumerable<RestrictionEnzyme> enzymes, bool[]? blocked = null)
    {
        var (cuts, sitesByPosition) = CollectForwardCuts(seq, enzymes, MoleculeTopology.Linear, blocked);
        int n = seq.Length;

        if (cuts.Count == 0)
        {
            // No cuts - return entire sequence as one fragment; both ends are the original
            // molecule's termini, so neither carries an enzyme-generated overhang.
            yield return new DigestFragment(
                Sequence: seq,
                StartPosition: 0,
                Length: n,
                LeftEnzyme: null,
                RightEnzyme: null,
                FragmentNumber: 1);
            yield break;
        }

        // Generate fragments. Every cut lies in [1, Length-1] (the linear cut-feasibility filter
        // CutsAreInternal), and the cut positions are distinct, so every interval
        // [boundary_i, boundary_i+1) is non-empty: k cuts always yield exactly k+1 fragments,
        // numbered 1..k+1 without a gap. FragmentNumber counts emitted fragments rather than
        // boundary indices so that the numbering stays contiguous even if a future cut source
        // were to admit a zero-length piece.
        var boundaries = new List<int>(cuts.Count + 2) { 0 };
        boundaries.AddRange(cuts);
        boundaries.Add(n);

        int fragmentNumber = 0;
        for (int i = 0; i < boundaries.Count - 1; i++)
        {
            int start = boundaries[i];
            int end = boundaries[i + 1];
            int length = end - start;

            if (length <= 0)
                continue;

            var leftSites = i == 0 ? null : sitesByPosition.GetValueOrDefault(start);
            var rightSites = i == boundaries.Count - 2 ? null : sitesByPosition.GetValueOrDefault(end);
            var leftSite = leftSites?[0];
            var rightSite = rightSites?[0];
            var leftEnd = DownstreamEnd(seq, leftSite, MoleculeTopology.Linear);
            var rightEnd = UpstreamEnd(seq, rightSite, MoleculeTopology.Linear);

            yield return new DigestFragment(
                Sequence: seq.Substring(start, length),
                StartPosition: start,
                Length: length,
                LeftEnzyme: leftSite?.Enzyme.Name,
                RightEnzyme: rightSite?.Enzyme.Name,
                FragmentNumber: ++fragmentNumber,
                LeftOverhangType: leftEnd.Type,
                LeftOverhang: leftEnd.Sequence,
                RightOverhangType: rightEnd.Type,
                RightOverhang: rightEnd.Sequence)
            {
                LeftEnzymes = EnzymeNames(leftSites),
                RightEnzymes = EnzymeNames(rightSites)
            };
        }
    }

    /// <summary>
    /// Simulates a restriction digest of a molecule with the given topology (linear or circular)
    /// and returns the resulting fragments.
    /// </summary>
    /// <remarks>
    /// <para>
    /// <b>Linear</b> (<see cref="MoleculeTopology.Linear"/>) is the default behaviour and is
    /// identical to <see cref="Digest(DnaSequence, string[])"/>: with <c>k</c> distinct
    /// top-strand cut positions a linear molecule yields <c>k+1</c> fragments; with zero
    /// cuts it yields the single whole linear molecule.
    /// </para>
    /// <para>
    /// <b>Circular</b> (<see cref="MoleculeTopology.Circular"/>) models a closed plasmid. The
    /// site search itself is circular, so a recognition site that spans the origin is found and
    /// cut (Biopython <c>search(..., linear=False)</c>). With <c>k</c> distinct top-strand cut
    /// positions a circular molecule yields exactly <c>k</c> fragments — "cutting a circular DNA
    /// molecule at k sites produces k fragments" (Univ. of Illinois MolBio "Restriction
    /// Digestion/Gel Electrophoresis": "If you cut a circle once, you get one linear
    /// fragment... cut it a second time to get 2 linear fragments"; Addgene Plasmids 101:
    /// a single cutter linearizes the plasmid into one fragment). Each fragment spans one cut
    /// to the next cut; the fragment that crosses the origin is the join of the
    /// last-cut→end piece and the start→first-cut piece, with length
    /// <c>(SequenceLength − lastCut) + firstCut</c>. A circular molecule with <b>zero</b> cut
    /// sites yields a single full-length (uncut) circular fragment — matching Biopython, whose
    /// <c>catalyse(..., linear=False)</c> returns the whole sequence as one fragment when the
    /// enzyme does not cut; a circular molecule with a <b>single</b> cut site is linearized into
    /// one full-length fragment.
    /// </para>
    /// <para>
    /// Fragments are ordered by ascending <see cref="DigestFragment.StartPosition"/>, so on a
    /// circle the origin-spanning fragment comes <b>last</b>. Biopython returns the same multiset
    /// with the origin-spanning fragment <b>first</b>; the two orders differ by a rotation of one.
    /// </para>
    /// <para>
    /// Each fragment end created by a cut also reports the enzyme's single-stranded extension —
    /// <see cref="DigestFragment.LeftOverhangType"/>/<see cref="DigestFragment.LeftOverhang"/> and
    /// <see cref="DigestFragment.RightOverhangType"/>/<see cref="DigestFragment.RightOverhang"/> —
    /// read 5'→3' on the strand that carries it, taken from the molecule itself. An end that is an
    /// original terminus of a linear molecule (the first fragment's left end, the last fragment's
    /// right end) reports <see langword="null"/>, as does the uncut circular fragment, which has
    /// no ends at all.
    /// </para>
    /// </remarks>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <param name="enzymeNames">Names of restriction enzymes to use.</param>
    /// <returns>Collection of DNA fragments after digestion.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="enzymeNames"/> is null, or a name in it is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<DigestFragment> Digest(
        DnaSequence sequence,
        MoleculeTopology topology,
        params string[] enzymeNames)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        var enzymes = RequireEnzymeNames(enzymeNames, nameof(enzymeNames));

        return topology == MoleculeTopology.Linear
            ? DigestCore(sequence.Sequence, enzymes)
            : DigestCircular(sequence.Sequence, enzymes);
    }

    /// <summary>
    /// Simulates a restriction digest with caller-supplied enzyme definitions — including enzymes
    /// that are not in the built-in table, such as Type IIS cutters whose cut lies outside the
    /// recognition site.
    /// </summary>
    /// <remarks>
    /// Same contract as <see cref="Digest(DnaSequence, MoleculeTopology, string[])"/> (top-strand
    /// boundaries, fragment order, per-end overhangs, several enzymes per boundary ordered by
    /// <see cref="StringComparer.Ordinal"/> on <see cref="RestrictionEnzyme.Name"/>), with the
    /// enzymes taken as given. An enzyme listed twice, or two enzymes with the same name, count once
    /// per boundary.
    /// </remarks>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <param name="enzymes">Enzyme definitions to digest with; at least one.</param>
    /// <returns>Collection of DNA fragments after digestion.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or
    /// <paramref name="enzymes"/> is null, or <paramref name="enzymes"/> contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, or an enzyme's
    /// recognition sequence is empty or not IUPAC.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static IEnumerable<DigestFragment> Digest(
        DnaSequence sequence,
        MoleculeTopology topology,
        IEnumerable<RestrictionEnzyme> enzymes)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        var list = RequireEnzymeList(enzymes);

        return topology == MoleculeTopology.Linear
            ? DigestCore(sequence.Sequence, list)
            : DigestCircular(sequence.Sequence, list);
    }

    /// <summary>
    /// Validates a caller-supplied enzyme collection: non-null, non-empty, no null element, every
    /// recognition sequence a valid IUPAC pattern. Returns a snapshot of it.
    /// </summary>
    private static IReadOnlyList<RestrictionEnzyme> RequireEnzymeList(IEnumerable<RestrictionEnzyme> enzymes)
    {
        ArgumentNullException.ThrowIfNull(enzymes);
        var list = enzymes.ToList();
        if (list.Count == 0)
            throw new ArgumentException("At least one enzyme is required", nameof(enzymes));
        foreach (var enzyme in list)
        {
            if (enzyme is null)
                throw new ArgumentNullException(nameof(enzymes), "The enzyme collection contains null.");
            ValidateEnzymeLimits(enzyme, nameof(enzymes));
            _ = ValidatedPattern(enzyme.RecognitionSequence, nameof(enzymes), enzyme.Name);
        }

        return list;
    }

    private static IEnumerable<DigestFragment> DigestCircular(string seq, IEnumerable<RestrictionEnzyme> enzymes, bool[]? blocked = null)
    {
        // Circular site search: a recognition window may wrap past the origin, so the cut it makes
        // must be seen too. Using the linear search here would silently leave an origin-spanning
        // site uncut and additionally drop an edge cutter's cuts through the linear-only
        // cut-feasibility filter, which does not apply to a closed circle.
        var (cuts, sitesByPosition) = CollectForwardCuts(seq, enzymes, MoleculeTopology.Circular, blocked);

        int n = seq.Length;

        if (cuts.Count == 0)
        {
            // No cuts on a closed circle: a single full-length, uncut circular fragment. It has no
            // ends, so no overhangs. Biopython's catalyse(..., linear=False) likewise returns the
            // whole sequence as one fragment for a non-cutter.
            yield return new DigestFragment(
                Sequence: seq,
                StartPosition: 0,
                Length: n,
                LeftEnzyme: null,
                RightEnzyme: null,
                FragmentNumber: 1);
            yield break;
        }

        // k cut sites on a circle → k fragments. Each fragment runs from one cut to the next
        // cut (cyclically). The fragment that wraps past the origin joins the last-cut→end
        // piece with the start→first-cut piece.
        for (int i = 0; i < cuts.Count; i++)
        {
            int start = cuts[i];
            int nextCut = cuts[(i + 1) % cuts.Count];

            string fragmentSeq;
            int length;
            if (i < cuts.Count - 1)
            {
                length = nextCut - start;
                fragmentSeq = seq.Substring(start, length);
            }
            else
            {
                // Origin-spanning fragment: last-cut → end, then start → first-cut. With a single
                // cut (nextCut == start) this is the whole linearized circle.
                length = (n - start) + nextCut;
                fragmentSeq = string.Concat(seq.AsSpan(start, n - start), seq.AsSpan(0, nextCut));
            }

            var leftSites = sitesByPosition.GetValueOrDefault(start);
            var rightSites = sitesByPosition.GetValueOrDefault(nextCut);
            var leftSite = leftSites?[0];
            var rightSite = rightSites?[0];
            var leftEnd = DownstreamEnd(seq, leftSite, MoleculeTopology.Circular);
            var rightEnd = UpstreamEnd(seq, rightSite, MoleculeTopology.Circular);

            yield return new DigestFragment(
                Sequence: fragmentSeq,
                StartPosition: start,
                Length: length,
                // On a circle every fragment end is a cut: both flanks carry an enzyme.
                LeftEnzyme: leftSite?.Enzyme.Name,
                RightEnzyme: rightSite?.Enzyme.Name,
                FragmentNumber: i + 1,
                LeftOverhangType: leftEnd.Type,
                LeftOverhang: leftEnd.Sequence,
                RightOverhangType: rightEnd.Type,
                RightOverhang: rightEnd.Sequence)
            {
                LeftEnzymes = EnzymeNames(leftSites),
                RightEnzymes = EnzymeNames(rightSites)
            };
        }
    }

    /// <summary>
    /// The top-strand (watson) cut of <paramref name="site"/>. A forward-strand record carries it
    /// as <see cref="RestrictionSite.CutPosition"/>; a reverse-strand record carries the
    /// bottom-strand cut <c>p + L − CutPositionForward</c>, and its top-strand cut is
    /// <c>p + L − CutPositionReverse</c>, i.e. the reported cut plus
    /// <c>CutPositionForward − CutPositionReverse</c> (Biopython <c>_rev_modify</c>:
    /// <c>location − fst3</c>). For the second cut of a two-sided (Type IIB) enzyme the offsets are
    /// its <see cref="RestrictionEnzyme.SecondCut"/> pair (<c>location − scd3</c>).
    /// </summary>
    private static int TopStrandCut(RestrictionSite site, int length, MoleculeTopology topology)
    {
        if (site.IsForwardStrand)
            return site.CutPosition;

        var (forward, reverse) = site.CutOffsets;
        long watson = (long)site.CutPosition + forward - reverse;
        return topology == MoleculeTopology.Circular ? Modulo(watson, length) : (int)watson;
    }

    /// <summary>
    /// Whether reading a window on the reverse strand cuts the same bonds as reading it on the
    /// forward strand, so a window that matches the recognition pattern both ways yields one
    /// cleavage, not two.
    /// </summary>
    /// <remarks>
    /// The forward reading of a window at <c>p</c> cuts at <c>(p + fwd, p + rev)</c> per cut pair;
    /// the reverse reading cuts at <c>(p + len - rev, p + len - fwd)</c>. The two sets are equal
    /// independently of <c>p</c>, so this is a property of the enzyme alone: a single-pair enzyme
    /// qualifies exactly when <c>fwd + rev == len</c> (equivalently Biopython's
    /// <c>fst5 + fst3 == 0</c>), and a two-sided (Type IIB) enzyme also qualifies when its two
    /// pairs are each other's mirror, which is how AlfI, BdaI, BplI, FalI, NmeDI, RdeGBIII and
    /// UcoMSI keep one cleavage per window although their individual pairs are asymmetric.
    /// Measured over the 754 REBASE enzymes with a defined cut: 166 have a cut set that differs
    /// under the reverse reading, and exactly 5 of those admit a window matching the pattern on
    /// both strands — HauII (TGGCCA), SgeI (CNNG), MspJI (CNNR), LpnPI (CCDG), SgrTI (CCDS).
    /// </remarks>
    private static bool ReverseReadingCutsAlike(RestrictionEnzyme enzyme)
    {
        // long: len - rev wraps in int for an offset within len of int.MinValue, and the wrapped value can equal fwd - the circular digest
        // and map then dropped a real cleavage (B08 A7-2, the arithmetic the other cut offsets already had since A5R-7).
        long len = enzyme.RecognitionSequence.Length;
        long fwd = enzyme.CutPositionForward;
        long rev = enzyme.CutPositionReverse;

        if (enzyme.SecondCut is not { } second)
            return fwd == len - rev;

        return (fwd == len - rev && second.Forward == len - second.Reverse)
            || (fwd == len - second.Reverse && second.Forward == len - rev);
    }

    /// <summary>
    /// The names of the enzymes cutting at one boundary, in the boundary's (ordinal) order; empty
    /// when the end is not a cut.
    /// </summary>
    private static IReadOnlyList<string> EnzymeNames(List<RestrictionSite>? sites)
        => sites is null ? [] : sites.ConvertAll(s => s.Enzyme.Name).AsReadOnly();

    /// <summary>
    /// The single-stranded extension at the <b>left</b> (5'-side) end of the fragment that starts
    /// at <paramref name="site"/>'s top-strand cut, written 5'→3' on the strand that carries it.
    /// </summary>
    private static (OverhangType? Type, string? Sequence) DownstreamEnd(
        string seq, RestrictionSite? site, MoleculeTopology topology)
    {
        if (site is null)
            return (null, null);

        var (type, window) = CutOverhangWindow(seq, site, topology);

        // 5' overhang: the extension is the top strand itself, so it reads as the window.
        // 3' overhang: the extension sits on the bottom strand, whose 5'→3' reading of the same
        // window is its reverse complement.
        return type switch
        {
            OverhangType.FivePrime => (type, window),
            OverhangType.ThreePrime => (type, IupacPattern.ReverseComplement(window)),
            _ => (type, string.Empty)
        };
    }

    /// <summary>
    /// The single-stranded extension at the <b>right</b> (3'-side) end of the fragment that ends
    /// at <paramref name="site"/>'s top-strand cut, written 5'→3' on the strand that carries it.
    /// </summary>
    private static (OverhangType? Type, string? Sequence) UpstreamEnd(
        string seq, RestrictionSite? site, MoleculeTopology topology)
    {
        if (site is null)
            return (null, null);

        var (type, window) = CutOverhangWindow(seq, site, topology);

        // Mirror image of DownstreamEnd: at the upstream side of the same cut a 5' overhang sits
        // on the bottom strand and a 3' overhang on the top strand.
        return type switch
        {
            OverhangType.FivePrime => (type, IupacPattern.ReverseComplement(window)),
            OverhangType.ThreePrime => (type, window),
            _ => (type, string.Empty)
        };
    }

    /// <summary>
    /// The overhang type of <paramref name="site"/>'s cut and the top-strand bases lying between
    /// its two strand cuts (empty for a blunt cutter).
    /// </summary>
    /// <remarks>
    /// Reads the enzyme's overhang window — the same <c>[min(cf, cr), max(cf, cr))</c> interval
    /// that <see cref="GetOverhangSequence(RestrictionEnzyme)"/> reads off the recognition pattern,
    /// via the shared <c>OverhangWindow</c> geometry — off the MOLECULE instead, translated
    /// so that the enzyme's forward cut sits on this site's reported cut position. For a built-in
    /// enzyme on a site that matches its pattern exactly the two therefore agree base for base;
    /// they differ only where the pattern is degenerate (<c>SfiI</c>: <c>"NNN"</c> from the pattern,
    /// the site's actual three bases from the molecule).
    /// </remarks>
    private static (OverhangType Type, string Window) CutOverhangWindow(
        string seq, RestrictionSite site, MoleculeTopology topology)
    {
        var enzyme = site.Enzyme;
        var type = enzyme.OverhangType;
        if (type == OverhangType.Blunt)
            return (type, string.Empty);

        var (forward, reverse) = site.CutOffsets;
        var (start, length) = OverhangWindow(forward, reverse);
        int n = seq.Length;
        // The window lies between the two strand cuts; its offset from the top-strand cut is the
        // same for both orientations of the site (bottom − top = reverse − forward either way), so
        // it is read relative to the top-strand cut. A two-sided enzyme's second cut uses its own
        // pair, whose stagger equals the first's (RestrictionEnzyme.SecondCut).
        int lo = TopStrandCut(site, n, topology) + (start - forward);

        if (topology == MoleculeTopology.Linear)
            return (type, seq.Substring(lo, length));

        // On a circle the window may run past either edge; read it modulo the length.
        var builder = new System.Text.StringBuilder(length);
        for (int k = 0; k < length; k++)
            builder.Append(seq[((lo + k) % n + n) % n]);

        return (type, builder.ToString());
    }

    /// <summary>
    /// Gets a summary of a <b>linear</b> digest with fragment sizes sorted descending.
    /// </summary>
    /// <remarks>
    /// Equivalent to <see cref="GetDigestSummary(DnaSequence, MoleculeTopology, string[])"/> with
    /// <see cref="MoleculeTopology.Linear"/>.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="enzymeNames"/> is null, or a name in it is null.</exception>
    public static DigestSummary GetDigestSummary(DnaSequence sequence, params string[] enzymeNames)
        => GetDigestSummary(sequence, MoleculeTopology.Linear, enzymeNames);

    /// <summary>
    /// Gets a summary of a digest of a molecule with the given topology, with fragment sizes
    /// sorted descending.
    /// </summary>
    /// <param name="sequence">DNA sequence to digest.</param>
    /// <param name="topology">Whether the molecule is linear or circular; see
    /// <see cref="Digest(DnaSequence, MoleculeTopology, string[])"/>. With <c>k</c> cut positions
    /// a linear molecule gives <c>k+1</c> fragments and a circular one <c>k</c>, so the two
    /// summaries differ.</param>
    /// <param name="enzymeNames">Names of restriction enzymes to use.</param>
    /// <remarks>
    /// <see cref="DigestSummary.FragmentSizes"/> are the fragment lengths in descending order (the
    /// gel-electrophoresis reading order); <see cref="DigestSummary.Topology"/> echoes
    /// <paramref name="topology"/>. The fragment lengths always sum to
    /// <paramref name="sequence"/>'s length in both topologies (Addgene Restriction Digest
    /// Protocol: "the sum of the individual fragments should equal the size of the original").
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="enzymeNames"/> is null, or a name in it is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static DigestSummary GetDigestSummary(
        DnaSequence sequence,
        MoleculeTopology topology,
        params string[] enzymeNames)
    {
        RequireTopology(topology);
        var fragments = Digest(sequence, topology, enzymeNames).ToList();
        return Summarize(fragments, enzymeNames.ToList(), topology);
    }

    /// <summary>
    /// Gets a summary of a digest with caller-supplied enzyme definitions — including enzymes that
    /// are not in the built-in table, such as Type IIS / Type IIB cutters.
    /// </summary>
    /// <remarks>
    /// Same contract as <see cref="GetDigestSummary(DnaSequence, MoleculeTopology, string[])"/>,
    /// over <see cref="Digest(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme})"/>;
    /// <see cref="DigestSummary.EnzymesUsed"/> lists the enzymes' <see cref="RestrictionEnzyme.Name"/>s
    /// in the order given.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or
    /// <paramref name="enzymes"/> is null, or <paramref name="enzymes"/> contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, or an enzyme
    /// definition is invalid (see <see cref="MaxCutStagger"/>, <see cref="MaxRecognitionLength"/> and <see cref="MaxEnzymeNameLength"/> for the size limits).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static DigestSummary GetDigestSummary(
        DnaSequence sequence,
        MoleculeTopology topology,
        IEnumerable<RestrictionEnzyme> enzymes)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        var list = RequireEnzymeList(enzymes);
        var fragments = Digest(sequence, topology, list).ToList();
        return Summarize(fragments, list.Select(e => e.Name).ToList(), topology);
    }

    private static DigestSummary Summarize(
        List<DigestFragment> fragments, List<string> enzymesUsed, MoleculeTopology topology)
    {
        var sizes = fragments.Select(f => f.Length).OrderByDescending(x => x).ToList();

        return new DigestSummary(
            TotalFragments: fragments.Count,
            FragmentSizes: sizes,
            LargestFragment: sizes.FirstOrDefault(),
            SmallestFragment: sizes.LastOrDefault(),
            AverageFragmentSize: sizes.Count > 0 ? sizes.Average() : 0,
            EnzymesUsed: enzymesUsed,
            Topology: topology);
    }

    #endregion

    #region Restriction Map

    /// <summary>
    /// Creates a restriction map of a <b>linear</b> molecule showing all enzyme sites.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="enzymeNames"/> is null, or a name in it is null.</exception>
    public static RestrictionMap CreateMap(DnaSequence sequence, params string[] enzymeNames)
        => CreateMap(sequence, MoleculeTopology.Linear, enzymeNames);

    /// <summary>
    /// Creates a restriction map of a molecule with the given topology, showing all enzyme sites.
    /// </summary>
    /// <param name="sequence">DNA sequence to map.</param>
    /// <param name="topology">Whether the molecule is linear or circular; see
    /// <see cref="FindSites(DnaSequence, string, MoleculeTopology)"/>.</param>
    /// <param name="enzymeNames">Enzymes to consider. Names are resolved case-insensitively and
    /// de-duplicated (any REBASE enzyme with a defined cut, <see cref="GetEnzyme"/>); an <b>empty</b>
    /// set means every enzyme of the common set (<see cref="Enzymes"/>); for the full REBASE set pass
    /// <see cref="RebaseEnzymes"/><c>.Values</c> to
    /// <see cref="CreateMap(DnaSequence, MoleculeTopology, IEnumerable{RestrictionEnzyme})"/>.</param>
    /// <remarks>
    /// <para>
    /// The <i>considered set</i> is what <see cref="RestrictionMap.NonCutters"/> is taken from:
    /// every considered enzyme with no site in <paramref name="sequence"/>, reported under its
    /// canonical table name. With an empty <paramref name="enzymeNames"/> the considered set is
    /// the whole built-in table, so <c>NonCutters</c> lists every built-in enzyme that does not
    /// cut (<c>SitesByEnzyme</c> ∪ <c>NonCutters</c> always covers the considered set exactly).
    /// </para>
    /// <para>
    /// <c>SitesByEnzyme</c> is keyed case-insensitively by canonical enzyme name and holds the
    /// positions of both strands' records. <c>TotalSites</c> and <c>UniqueCutters</c> count
    /// <b>cleavages</b>, by the same rule the digest's cut collection applies: every forward-strand
    /// record, plus every reverse-strand record whose backwards reading cuts other bonds. So a
    /// palindromic site with symmetric cuts counts once (which is what "forward-strand records only"
    /// amounted to for the 39 built-ins, and the only reading this doc used to give), a
    /// non-palindromic enzyme's reverse-strand sites count on their own, the two cut records of a
    /// two-sided (Type IIB) site count once, and a window that admits two different cleavages counts
    /// twice. Consequences for cloning, which is what <c>UniqueCutters</c> is for: a Type IIS
    /// enzyme is <b>not</b> a unique cutter merely because it has one forward-strand site —
    /// <c>CreateMap("A{10}GGTCTCA{10}GAGACCA{10}", "BsaI")</c> gives <c>TotalSites = 2</c> and an
    /// empty <c>UniqueCutters</c>, because the molecule carries GGTCTC on one strand and on the
    /// other; and an enzyme that binds one window but cleaves it twice (HauII, SgeI, MspJI, LpnPI,
    /// SgrTI) is not a unique cutter either. Verified against Biopython 1.88 over 14580 rows for the
    /// 749 enzymes where the two agree (review 2026-09 B08 H32, on the rule F101 settled).
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentException">A name in <paramref name="enzymeNames"/> is not a known
    /// enzyme (<see cref="GetEnzyme"/>).</exception>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or <paramref name="enzymeNames"/> is null, or a name in it is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static RestrictionMap CreateMap(
        DnaSequence sequence, MoleculeTopology topology, params string[] enzymeNames)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzymeNames);

        // Resolve the request to canonical, de-duplicated enzymes so that NonCutters can be
        // compared against the keys of SitesByEnzyme (which carry the table's own spelling) and
        // so a name repeated with different casing is not searched — and counted — twice.
        var considered = (enzymeNames.Length > 0
                ? enzymeNames.Select(n => RequireEnzyme(n))
                : _enzymes.Values)
            .DistinctBy(e => e.Name, StringComparer.OrdinalIgnoreCase)
            .ToList();

        return CreateMapCore(sequence.Sequence, topology, considered);
    }

    /// <summary>
    /// Creates a restriction map with caller-supplied enzyme definitions — including enzymes that
    /// are not in the built-in table, such as Type IIS / Type IIB cutters.
    /// </summary>
    /// <param name="sequence">DNA sequence to map.</param>
    /// <param name="topology">Whether the molecule is linear or circular.</param>
    /// <param name="enzymes">Enzyme definitions to consider; at least one. Two definitions with the
    /// same <see cref="RestrictionEnzyme.Name"/> (ignoring case) count once — the first is used.</param>
    /// <remarks>
    /// Same contract as <see cref="CreateMap(DnaSequence, MoleculeTopology, string[])"/>, with the
    /// considered set being <paramref name="enzymes"/> (an empty collection is an error here, not
    /// "every built-in"). <see cref="RestrictionMap.TotalSites"/> and
    /// <see cref="RestrictionMap.UniqueCutters"/> count cleavages, exactly as that overload
    /// documents: every forward-strand record, plus every reverse-strand record whose backwards
    /// reading cuts other bonds — so a non-palindromic enzyme's reverse-strand sites count, a
    /// palindromic site with symmetric cuts counts once, the two cuts of a two-sided (Type IIB) site
    /// count once, and a window that admits two different cleavages counts twice (B08 H32/F101).
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> or
    /// <paramref name="enzymes"/> is null, or <paramref name="enzymes"/> contains null.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymes"/> is empty, or an enzyme
    /// definition is invalid (see <see cref="MaxCutStagger"/>, <see cref="MaxRecognitionLength"/> and <see cref="MaxEnzymeNameLength"/> for the size limits).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="topology"/> is not a defined <see cref="MoleculeTopology"/> value.</exception>
    public static RestrictionMap CreateMap(
        DnaSequence sequence, MoleculeTopology topology, IEnumerable<RestrictionEnzyme> enzymes)
    {
        RequireTopology(topology);
        ArgumentNullException.ThrowIfNull(sequence);
        var considered = RequireEnzymeList(enzymes)
            .DistinctBy(e => e.Name, StringComparer.OrdinalIgnoreCase)
            .ToList();

        return CreateMapCore(sequence.Sequence, topology, considered);
    }

    private static RestrictionMap CreateMapCore(
        string seq, MoleculeTopology topology, List<RestrictionEnzyme> considered, bool[]? blocked = null)
    {
        var sites = considered.SelectMany(e => FindSitesCore(seq, e, topology, blocked)).ToList();

        var byEnzyme = sites
            .GroupBy(s => s.Enzyme.Name, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(
                g => g.Key,
                g => g.Select(s => s.Position).OrderBy(p => p).ToList(),
                StringComparer.OrdinalIgnoreCase);

        // Cleavages per enzyme, the same rule the digest's cut collection applies
        // (CollectForwardCuts): a forward-strand record always counts; a reverse-strand record
        // counts unless reading the window backwards cuts the same bonds as reading it forwards.
        // So a palindromic site with symmetric cuts (all 39 built-ins) counts once, a
        // non-palindromic enzyme's reverse-strand sites count, the two cut records of a two-sided
        // (Type IIB) site share (Position, strand) and count once, and a window that admits two
        // different cleavages counts twice — HauII, SgeI, MspJI, LpnPI, SgrTI, which is why
        // UniqueCutters excludes them where they bind once but cut twice (B08 H7/H32, F101).
        var physicalByEnzyme = sites
            .Where(s => s.IsForwardStrand
                        || !IupacPattern.MatchesAtCircular(seq, s.Position, SequenceExtensions.ToUpperAscii(s.Enzyme.RecognitionSequence))
                        || !ReverseReadingCutsAlike(s.Enzyme))
            .GroupBy(s => s.Enzyme.Name, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(
                g => g.Key,
                g => g.Select(s => (s.Position, s.IsForwardStrand)).Distinct().Count(),
                StringComparer.OrdinalIgnoreCase);
        var uniqueCutters = physicalByEnzyme.Where(kv => kv.Value == 1).Select(kv => kv.Key).ToList();
        var nonCutters = considered.Select(e => e.Name).Where(e => !byEnzyme.ContainsKey(e)).ToList();

        return new RestrictionMap(
            SequenceLength: seq.Length,
            Sites: sites,
            SitesByEnzyme: byEnzyme,
            TotalSites: physicalByEnzyme.Values.Sum(),
            UniqueCutters: uniqueCutters,
            NonCutters: nonCutters);
    }

    #endregion

    #region Compatibility Analysis

    // ───────────────────────────────────────────────────────────────────────────────────────
    // Enzyme-level end chemistry and the ligation-compatibility rule.
    //
    // REFERENCE PRIMITIVES (Biopython 1.88 Bio.Restriction): per enzyme, `ovhg` is the signed
    // stagger between the two strand cuts (negative = 5' overhang, positive = 3' overhang,
    // 0 = blunt) and `ovhgseq` is the overhang sequence — the TOP-strand bases lying between the
    // two cuts, written 5'→3'. For an enzyme whose site contains an ambiguity code inside that
    // window, `ovhgseq` carries the code (SfiI `GGCCNNNNNGGCC` 8/5 → `ovhgseq == "NNN"`), i.e. the
    // reference reports the overhang the enzyme is DEFINED to leave, not a particular site's bases.
    //
    // THIS CLASS: `GetOverhangSequence(enzyme)` is exactly that window, taken from the recognition
    // pattern — `RecognitionSequence[min(cf, cr) .. max(cf, cr))`. All 39 built-in entries were
    // compared base for base against `ovhgseq` this session: 39/39 identical.
    //
    // LIGATION RULE (derived from those primitives — see the Evidence file for why no published
    // compatible-ends table could be opened from this container):
    //   • One cut yields two ends whose single-stranded extensions are W (the downstream
    //     fragment's end, read on the top strand) and revcomp(W) (the upstream fragment's end,
    //     read on the bottom strand). They are complementary and antiparallel, which is exactly
    //     what T4 ligase needs, so they re-ligate — and an enzyme therefore OFFERS both.
    //   • Two ends anneal iff their extensions have the same length, the same polarity (both 5' or
    //     both 3' — two extensions of opposite polarity sit on the same strand of the joint and
    //     cannot base-pair) and one is the reverse complement of the other.
    //   • So enzymes E1 and E2 have compatible ends iff they share a polarity and
    //     W1 ∈ { W2, revcomp(W2) }. Blunt ends carry no extension and ligate to any blunt end.
    //   • An ambiguity code in a window stands for the set of bases it expands to, independently
    //     per position, so "W1 matches W2" means the two code strings INTERSECT position by
    //     position (`IupacHelper.MatchesIupac` composed over A/C/G/T) — i.e. some real pair of
    //     sites exists that ligates. A degenerate answer is therefore "possibly compatible".
    //
    // CROSS-CHECK against Biopython's own `%` operator (`RestrictionType.__mod__`, which tests
    // same overhang class then `ovhgseq` equality, with a regex expansion when one side is
    // `Ambiguous`):
    //   • 39 built-ins, all 1521 ordered pairs: 0 disagreements.
    //   • 389 REBASE enzymes whose `ovhgseq` is blunt or fully A/C/G/T, all 75 466 unordered
    //     pairs: 0 disagreements.
    //   • Beyond those, Biopython's `%` is not usable as an oracle and is NOT copied: it is
    //     asymmetric for 2386 of 283 881 REBASE pairs (`_mod1` dispatches on the left operand and
    //     only `Defined._mod2` falls back to the other side), and `Ambiguous._mod2` builds a broken
    //     regex when a code repeats in the overhang — it splits its own character class on the
    //     code, so `StyI % PscI` ("CWWG" vs "CATG") is False where the biology and this class say
    //     compatible. The rule above is symmetric by construction.
    //   • The reverse-complement branch is not reachable from the built-in table (every built-in
    //     window is its own reverse complement, a consequence of a palindromic site cut
    //     symmetrically), but it is real: 12 REBASE pairs are compatible only through it, e.g.
    //     `BsmI` ("CN", 3') with `BsrI` ("GN", 3') — BsmI's bottom-strand end can read `GG` and
    //     BsrI's `CC`. It is reachable here through the `RestrictionEnzyme` overloads.
    // ───────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// The single-stranded overhang a built-in enzyme leaves, written 5'→3' on the top strand —
    /// REBASE/Biopython <c>ovhgseq</c>.
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <returns><see cref="string.Empty"/> for a blunt cutter; otherwise the recognition-pattern
    /// bases between the two strand cuts, which may contain IUPAC ambiguity codes when the
    /// enzyme's own site does (<c>SfiI</c> → <c>"NNN"</c>).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is empty or not a known enzyme (<see cref="GetEnzyme"/>).</exception>
    public static string GetOverhangSequence(string enzymeName)
        => GetOverhangSequence(RequireEnzyme(enzymeName))!;

    /// <summary>
    /// The single-stranded overhang <paramref name="enzyme"/> leaves, written 5'→3' on the top
    /// strand — REBASE/Biopython <c>ovhgseq</c>.
    /// </summary>
    /// <param name="enzyme">Enzyme definition; may be a caller-supplied (non-built-in) enzyme.</param>
    /// <returns><see cref="string.Empty"/> for a blunt cutter; otherwise the bases between the two
    /// strand cuts, read off the recognition pattern extended by <c>N</c> on both sides. For an
    /// enzyme cutting inside its site this is the pattern slice (EcoRI <c>AATT</c>, SfiI
    /// <c>NNN</c>); a window reaching outside the site reads <c>N</c> there, because the bases are
    /// the flanking molecule's and can be anything: BsaI <c>GGTCTC(1/5)</c> → <c>NNNN</c>, TspRI
    /// <c>NNCASTGNN</c>, BaeI → <c>NNNNN</c>. This is Biopython 1.88's <c>ovhgseq</c>, which classes
    /// such an enzyme as <c>Ambiguous</c> — verified for all 178 distal or two-cut enzymes of its
    /// snapshot (review 2026-09 B08 F57). For a two-sided (Type IIB) enzyme it is the first cut
    /// pair's window; Biopython stores one <c>ovhgseq</c>, and both windows are identical for all
    /// 25 of its two-cut enzymes. Never <see langword="null"/> (the nullable return type is kept
    /// for source compatibility; before F57 a distal cutter returned <see langword="null"/>).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="enzyme"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzyme"/> has an empty or non-IUPAC recognition sequence, is larger than <see cref="MaxCutStagger"/> / <see cref="MaxRecognitionLength"/> / <see cref="MaxEnzymeNameLength"/> allow, or has a second cut pair that leaves a different stagger.</exception>
    public static string? GetOverhangSequence(RestrictionEnzyme enzyme)
    {
        ArgumentNullException.ThrowIfNull(enzyme);
        ValidateEnzymeLimits(enzyme, nameof(enzyme));
        // The same definition check as FindSites: an empty or non-IUPAC recognition sequence is an ArgumentException here too, and a
        // lower-case one is read as upper case (it used to come back in lower case, or as "NNNN" when empty).
        string site = ValidatedPattern(enzyme.RecognitionSequence, nameof(enzyme), enzyme.Name);

        var (start, length) = OverhangWindow(enzyme);
        if (length == 0)
            return string.Empty;

        // A base of the window outside the recognition pattern belongs to the flanking molecule
        // and is unknown, i.e. N. Biopython 1.88's stored ovhgseq (Restriction_Dictionary) equals
        // this padded window for all 178 distal / two-cut enzymes of its snapshot (F57).
        var window = new System.Text.StringBuilder(length);
        for (int k = start; k < start + length; k++)
            window.Append(k >= 0 && k < site.Length ? site[k] : 'N');

        return window.ToString();
    }

    /// <summary>
    /// The overhang window of <paramref name="enzyme"/> as an offset/length pair relative to the
    /// 5' end of its recognition sequence — the single definition of the window's geometry, shared
    /// by <see cref="GetOverhangSequence(RestrictionEnzyme)"/> (which reads it off the recognition
    /// pattern) and by <c>CutOverhangWindow</c> (which reads the same window off the molecule).
    /// </summary>
    /// <remarks>
    /// The window is <c>[min(cf, cr), max(cf, cr))</c>: the bases lying between the top-strand and
    /// bottom-strand cuts. Its length is <c>|ovhg|</c> and it is empty for a blunt cutter.
    /// </remarks>
    private static (int Start, int Length) OverhangWindow(RestrictionEnzyme enzyme) =>
        OverhangWindow(enzyme.CutPositionForward, enzyme.CutPositionReverse);

    /// <summary>The overhang window of one cut pair (see <see cref="OverhangWindow(RestrictionEnzyme)"/>).</summary>
    // The stagger is at most MaxCutStagger by the time a window is asked for (ValidateEnzymeLimits), so the narrowing is exact.
    private static (int Start, int Length) OverhangWindow(int forward, int reverse) =>
        (Math.Min(forward, reverse), (int)Math.Abs((long)forward - reverse));

    /// <summary>
    /// Enumerates every unordered pair of distinct common-set enzymes (<see cref="Enzymes"/>) whose ends can be ligated to
    /// each other.
    /// </summary>
    /// <returns>Tuples of the two enzyme names and the shared end: <c>"blunt"</c> for a blunt pair,
    /// otherwise the first enzyme's overhang sequence (see
    /// <see cref="GetOverhangSequence(string)"/>).</returns>
    /// <remarks>
    /// Only pairs of <i>different</i> enzymes are returned. Every enzyme is trivially compatible
    /// with itself, and 15 of the 39 built-ins (e.g. <c>EcoRI</c>, <c>PstI</c>, <c>SfiI</c>) have a
    /// unique overhang and so appear in no pair at all; use
    /// <see cref="FindCompatibleEnzymes(string)"/> to ask what a single enzyme can be ligated to.
    /// For the full REBASE set use <see cref="FindCompatibleEnzymes(EnzymeSet)"/>.
    /// </remarks>
    public static IEnumerable<(string Enzyme1, string Enzyme2, string CompatibleEnd)> FindCompatibleEnzymes()
        => FindCompatibleEnzymes(EnzymeSet.Common);

    /// <summary>
    /// Lists the common-set enzymes (<see cref="Enzymes"/>; for another set use
    /// <see cref="FindCompatibleEnzymes(string, EnzymeSet)"/>) whose ends can be ligated to an end left by
    /// <paramref name="enzymeName"/> — the enzyme itself included when it is a member of that set.
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <returns>Tuples of the compatible enzyme's name (as spelled in the table) and the end:
    /// <c>"blunt"</c> for a blunt pair, otherwise <b><paramref name="enzymeName"/>'s own</b> overhang
    /// window (<see cref="GetOverhangSequence(string)"/>), not a string both enzymes necessarily
    /// leave. The two coincide whenever the two windows intersect directly, but 12 REBASE pairs
    /// anneal only through the reverse complement and there the partner leaves the complementary
    /// window: all 12 are <c>GN</c> against <c>CN</c> —
    /// <c>FindCompatibleEnzymes("BsmI")</c> reports BsrI with <c>CN</c> (BsmI's end) and
    /// <c>FindCompatibleEnzymes("BsrI")</c> reports BsmI with <c>GN</c> (BsrI's end), and BsrI never
    /// leaves <c>CN</c>. Call <see cref="GetOverhangSequence(string)"/> on the partner for its own
    /// end (review 2026-09 B08 H30).</returns>
    /// <remarks>
    /// <para>
    /// This is the per-enzyme form of Biopython's <c>compatible_end()</c>, which has <b>three</b>
    /// implementations, one per end class, not one 5'-only filter:
    /// <c>Blunt.compatible_end</c> returns every blunt enzyme, <c>Ov5.compatible_end</c> every
    /// <c>x</c> with <c>x.is_5overhang() and x % cls</c>, and <c>Ov3.compatible_end</c> every
    /// <c>x</c> with <c>x.is_3overhang() and x % cls</c>. So the reference answers a non-empty list
    /// for blunt and 3'-overhang enzymes too, and includes the queried enzyme itself. Measured with
    /// Biopython 1.88: <c>EcoRV.compatible_end()</c> → 137 enzymes (every blunt cutter),
    /// <c>PstI.compatible_end()</c> → 13, <c>SacI.compatible_end()</c> → 9,
    /// <c>EcoRI.compatible_end()</c> → 51 (review 2026-09 B08 H33 — this doc used to claim the
    /// reference answers the empty list for every blunt and every 3'-overhang enzyme). Note also
    /// that the reference ignores its own <c>batch</c> argument and always iterates
    /// <c>AllEnzymes</c>; this overload honours its <c>EnzymeSet</c>.
    /// </para>
    /// <para>
    /// This overload applies no end-class filter of its own and is simply
    /// <see cref="AreCompatible(RestrictionEnzyme, RestrictionEnzyme)"/> over the chosen table. The
    /// queried enzyme appears in the result <b>only when it is itself a member of that table</b>:
    /// <c>FindCompatibleEnzymes("BsaI", EnzymeSet.Common)</c> cannot list BsaI, which is not one of
    /// the 39 common enzymes, although <c>GetEnzyme</c> resolves it (H33).
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is empty or not a known enzyme (<see cref="GetEnzyme"/>).</exception>
    public static IEnumerable<(string Enzyme, string CompatibleEnd)> FindCompatibleEnzymes(string enzymeName)
        => FindCompatibleEnzymes(enzymeName, EnzymeSet.Common);

    /// <summary>
    /// Checks whether two named enzymes (<see cref="GetEnzyme"/>) produce ligatable (compatible) ends.
    /// </summary>
    /// <param name="enzyme1Name">First enzyme name (case-insensitive).</param>
    /// <param name="enzyme2Name">Second enzyme name (case-insensitive).</param>
    /// <returns><see langword="true"/> when both are blunt cutters, or when both leave the same
    /// polarity of overhang and the overhang sequences can anneal (see the rule documented on this
    /// region).</returns>
    /// <remarks>
    /// An unknown enzyme name <b>throws</b>, as it does in <see cref="FindSites(DnaSequence, string)"/>
    /// and <see cref="CreateMap(DnaSequence, string[])"/>: a misspelled name is not an answer about
    /// compatibility, and returning <see langword="false"/> made the two indistinguishable. Use
    /// <see cref="TryAreCompatible"/> to probe a name that may not exist.
    /// </remarks>
    /// <exception cref="ArgumentNullException">Either name is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException">Either name is empty or not a known enzyme (<see cref="GetEnzyme"/>).</exception>
    public static bool AreCompatible(string enzyme1Name, string enzyme2Name)
    {
        var e1 = RequireEnzyme(enzyme1Name, nameof(enzyme1Name));
        var e2 = RequireEnzyme(enzyme2Name, nameof(enzyme2Name));

        return AreCompatible(e1, e2, out _);
    }

    /// <summary>
    /// Lenient form of <see cref="AreCompatible(string, string)"/> for callers that cannot
    /// guarantee both names resolve (<see cref="GetEnzyme"/>).
    /// </summary>
    /// <param name="enzyme1Name">First enzyme name (case-insensitive).</param>
    /// <param name="enzyme2Name">Second enzyme name (case-insensitive).</param>
    /// <param name="compatible">The compatibility verdict when both names resolved; otherwise
    /// <see langword="false"/>.</param>
    /// <returns><see langword="true"/> when both names were resolved, so
    /// <paramref name="compatible"/> is meaningful; <see langword="false"/> when either name is
    /// unknown.</returns>
    /// <exception cref="ArgumentNullException">Either name is <see langword="null"/> (an empty or unknown name answers <c>false</c>).</exception>
    public static bool TryAreCompatible(string enzyme1Name, string enzyme2Name, out bool compatible)
    {
        compatible = false;

        var e1 = GetEnzyme(enzyme1Name);
        var e2 = GetEnzyme(enzyme2Name);
        if (e1 is null || e2 is null)
            return false;

        compatible = AreCompatible(e1, e2, out _);
        return true;
    }

    /// <summary>
    /// Checks whether two enzyme definitions produce ligatable (compatible) ends. Accepts
    /// caller-supplied enzymes, not only built-ins.
    /// </summary>
    /// <param name="enzyme1">First enzyme definition.</param>
    /// <param name="enzyme2">Second enzyme definition.</param>
    /// <returns><see langword="true"/> when both are blunt cutters, or when both leave the same
    /// polarity of overhang and the overhang sequences can anneal (see the rule documented on this
    /// region). A cutter whose overhang window reaches outside its own site (Type IIS/IIB distal
    /// cut) takes part with its <see cref="GetOverhangSequence(RestrictionEnzyme)"/> window, whose
    /// outside bases are <c>N</c> — Biopython classes it <c>Ambiguous</c> (BsaI <c>ovhgseq</c>
    /// <c>NNNN</c>) and its <c>%</c> matches such an <c>N</c> against any base: BsaI is
    /// "possibly compatible" with every 4-nt 5'-overhang end (BsaI % EcoRI is True), never with a
    /// 3' or a blunt one (review 2026-09 B08 F57).</returns>
    /// <exception cref="ArgumentNullException">Either enzyme is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException">Either enzyme has an empty or non-IUPAC recognition sequence, is larger than <see cref="MaxCutStagger"/> / <see cref="MaxRecognitionLength"/> / <see cref="MaxEnzymeNameLength"/> allow, or has a second cut pair that leaves a different stagger.</exception>
    public static bool AreCompatible(RestrictionEnzyme enzyme1, RestrictionEnzyme enzyme2)
    {
        ArgumentNullException.ThrowIfNull(enzyme1);
        ArgumentNullException.ThrowIfNull(enzyme2);
        ValidateEnzymeLimits(enzyme1, nameof(enzyme1));
        ValidateEnzymeLimits(enzyme2, nameof(enzyme2));
        _ = ValidatedPattern(enzyme1.RecognitionSequence, nameof(enzyme1), enzyme1.Name);
        _ = ValidatedPattern(enzyme2.RecognitionSequence, nameof(enzyme2), enzyme2.Name);

        return AreCompatible(enzyme1, enzyme2, out _);
    }

    private static bool AreCompatible(RestrictionEnzyme e1, RestrictionEnzyme e2, out string? compatibleEnd)
    {
        compatibleEnd = null;

        // Blunt ends carry no single-stranded extension, so any blunt end ligates to any other.
        // Source: Wikipedia (Sticky and blunt ends) — "blunt ends are always compatible with each
        // other"; Biopython `Blunt._mod1` = `issubclass(other, Blunt)`.
        if (e1.IsBluntEnd && e2.IsBluntEnd)
        {
            compatibleEnd = "blunt";
            return true;
        }

        // Opposite polarity (one 5', one 3') cannot anneal even when the overhang strings are the
        // same letters: both extensions would end up on the same strand of the joint. Blunt against
        // sticky falls out here too. Biopython: `Ov5._mod1`/`Ov3._mod1` require the same class.
        if (e1.OverhangType != e2.OverhangType)
            return false;

        string? overhang1 = GetOverhangSequence(e1);
        string? overhang2 = GetOverhangSequence(e2);
        if (string.IsNullOrEmpty(overhang1) || string.IsNullOrEmpty(overhang2))
            return false;

        // Same polarity: the ends anneal iff one enzyme's overhang window can read as the other's
        // window or as its reverse complement (the other end of the same cut).
        if (!OverhangsAnneal(overhang1, overhang2))
            return false;

        compatibleEnd = overhang1;
        return true;
    }

    /// <summary>
    /// Whether two same-polarity overhang windows, each possibly carrying IUPAC ambiguity codes,
    /// describe at least one pair of ends that can anneal.
    /// </summary>
    private static bool OverhangsAnneal(string overhang1, string overhang2)
    {
        if (overhang1.Length != overhang2.Length)
            return false;

        return CodesIntersect(overhang1, overhang2)
            || CodesIntersect(overhang1, IupacPattern.ReverseComplement(overhang2));
    }

    /// <summary>
    /// Whether two equal-length IUPAC code strings have a common expansion — i.e. at every
    /// position the two codes admit a common base. Composed from the canonical
    /// <see cref="IupacHelper.MatchesIupac"/> so there is no second ambiguity table in this class.
    /// </summary>
    private static bool CodesIntersect(string codes1, string codes2)
    {
        for (int i = 0; i < codes1.Length; i++)
        {
            if (!CodesIntersect(codes1[i], codes2[i]))
                return false;
        }

        return true;
    }

    private static bool CodesIntersect(char code1, char code2)
    {
        foreach (char nucleotide in "ACGT")
        {
            if (IupacHelper.MatchesIupac(nucleotide, code1) && IupacHelper.MatchesIupac(nucleotide, code2))
                return true;
        }

        return false;
    }

    #endregion

    #region Isoschizomers

    /// <summary>
    /// Common-set enzymes (<see cref="Enzymes"/>; other sets:
    /// <see cref="FindIsoschizomers(string, EnzymeSet)"/>) that recognise the same sequence as <paramref name="enzymeName"/> but are
    /// not that enzyme — its isoschizomers (Biopython <c>isoschizomers()</c>; the union of
    /// <see cref="FindEquischizomers"/> and <see cref="FindNeoschizomers"/>).
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is empty or not a known enzyme (<see cref="GetEnzyme"/>).</exception>
    public static IEnumerable<RestrictionEnzyme> FindIsoschizomers(string enzymeName)
        => FindIsoschizomers(enzymeName, EnzymeSet.Common);

    /// <summary>
    /// Common-set isoschizomers of <paramref name="enzymeName"/> that also cut at the same positions,
    /// so they are interchangeable in a digest — equischizomers (Biopython
    /// <c>equischizomers()</c>: "same site, same position of restriction").
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is empty or not a known enzyme (<see cref="GetEnzyme"/>).</exception>
    public static IEnumerable<RestrictionEnzyme> FindEquischizomers(string enzymeName)
        => FindEquischizomers(enzymeName, EnzymeSet.Common);

    /// <summary>
    /// Common-set isoschizomers of <paramref name="enzymeName"/> that cut the same site at a
    /// different position, so they leave different ends — neoschizomers (Biopython
    /// <c>neoschizomers()</c>: "same site, different position of restriction").
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is empty or not a known enzyme (<see cref="GetEnzyme"/>).</exception>
    public static IEnumerable<RestrictionEnzyme> FindNeoschizomers(string enzymeName)
        => FindNeoschizomers(enzymeName, EnzymeSet.Common);

    /// <summary>
    /// The one same-recognition-sequence query the three isoschizomer accessors share. Comparison
    /// is on the recognition string as stored (both are upper-case in the table), and the queried
    /// enzyme itself is excluded, exactly as the reference's <c>del r[i]</c> does.
    /// </summary>
    private static IEnumerable<RestrictionEnzyme> SameSiteEnzymes(
        RestrictionEnzyme query, IReadOnlyDictionary<string, RestrictionEnzyme> table)
    {
        return table.Values.Where(e =>
            !ReferenceEquals(e, query) &&
            !string.Equals(e.Name, query.Name, StringComparison.OrdinalIgnoreCase) &&
            string.Equals(e.RecognitionSequence, query.RecognitionSequence, StringComparison.OrdinalIgnoreCase));
    }

    #endregion
}

/// <summary>
/// Represents a restriction enzyme with its properties.
/// </summary>
/// <param name="Name">Enzyme name.</param>
/// <param name="RecognitionSequence">Recognition sequence (IUPAC), one strand, 5'→3'.</param>
/// <param name="CutPositionForward">Top-strand cut offset: the number of bases from the 5' end of
/// the site (as written) to the top-strand cut. It may be negative (cut upstream of the site) or
/// exceed the site length (cut downstream of it), so single-sided Type IIS cutters are expressible:
/// BsaI GGTCTC(1/5) is <c>(GGTCTC, 7, 11)</c>, FokI GGATG(9/13) <c>(GGATG, 14, 18)</c>, MlyI
/// GAGTC(5/5) <c>(GAGTC, 10, 10)</c>. Equals Biopython <c>fst5</c>.</param>
/// <param name="CutPositionReverse">Bottom-strand cut offset, counted in top-strand bases from the
/// same 5' end of the site; equals Biopython <c>size + fst3</c>, and
/// <c>CutPositionForward − CutPositionReverse</c> is Biopython <c>ovhg</c> (negative = 5' overhang).
/// REBASE <c>(n/m)</c> after a site of length <c>L</c> maps to <c>(L + n, L + m)</c>. Enzymes that cut
/// on <b>both</b> sides of the site (Type IIB, e.g. BaeI, BcgI — Biopython <c>TwoCuts</c>) carry their
/// second cut pair in <see cref="SecondCut"/>.</param>
/// <param name="Organism">Source organism.</param>
/// <param name="RebaseMethylationField">The REBASE <c>allenz</c> field <c>&lt;6&gt;</c>
/// ("METHYLATION SITE") verbatim, e.g. <c>"3(6)"</c> for EcoRI; <c>""</c> (the default) when REBASE
/// records no site or the enzyme is user-defined. Parsed into
/// <see cref="CognateMethylationSites"/>; a malformed value throws <see cref="FormatException"/>
/// at construction.</param>
public sealed record RestrictionEnzyme(
    string Name,
    string RecognitionSequence,
    int CutPositionForward,
    int CutPositionReverse,
    string Organism,
    string RebaseMethylationField = "")
{
    /// <summary>
    /// The REBASE methylation-site field, validated at construction.
    /// </summary>
    public string RebaseMethylationField { get; init; } =
        ValidateRebaseMethylationField(RebaseMethylationField);

    /// <summary>
    /// The sites the enzyme's <b>cognate methyltransferase</b> modifies, parsed from
    /// <see cref="RebaseMethylationField"/> (REBASE <c>allenz.104</c> field <c>&lt;6&gt;</c>).
    /// Empty when REBASE records none.
    /// </summary>
    /// <remarks>
    /// REBASE (allenz.104 format header, "&lt;METHYLATION SITE&gt;"): "The site of methylation by
    /// the cognate methylase when known is indicated X(Y) … where X is the base within the
    /// recognition sequence that is modified. A negative number indicates the complementary
    /// strand, numbered from the 5' base of that strand, and Y is the specific type of methylation
    /// involved: (6) = N6-methyladenosine, (5) = 5-methylcytosine, (4) = N4-methylcytosine." In the
    /// data a <c>?</c> stands for X when only the modification type is known, and several sites
    /// are comma-separated.
    /// <para>
    /// This is the <i>protection</i> site of the restriction-modification system the enzyme comes
    /// from. It is <b>not</b> a statement about the enzyme's sensitivity to Dam, Dcm or CpG
    /// methylation by other methyltransferases; no such behaviour is modelled.
    /// </para>
    /// </remarks>
    public IReadOnlyList<CognateMethylationSite> CognateMethylationSites =>
        ParseRebaseMethylationField(RebaseMethylationField);

    private static string ValidateRebaseMethylationField(string field)
    {
        ArgumentNullException.ThrowIfNull(field);
        _ = ParseRebaseMethylationField(field);
        return field;
    }

    private static IReadOnlyList<CognateMethylationSite> ParseRebaseMethylationField(string field)
    {
        if (field.Length == 0)
            return Array.Empty<CognateMethylationSite>();

        var sites = new List<CognateMethylationSite>();
        foreach (var raw in field.Split(','))
        {
            var token = raw.Trim();
            int open = token.IndexOf('(');
            if (open <= 0 || token[^1] != ')')
                throw new FormatException($"Malformed REBASE methylation-site entry '{raw}' in '{field}'.");

            var positionText = token[..open];
            var typeText = token[(open + 1)..^1];

            var type = typeText switch
            {
                "6" => MethylatedBase.N6Methyladenine,
                "5" => MethylatedBase.C5Methylcytosine,
                "4" => MethylatedBase.N4Methylcytosine,
                _ => throw new FormatException(
                    $"Unknown REBASE methylation type '({typeText})' in '{field}'; expected (4), (5) or (6).")
            };

            if (positionText == "?")
            {
                sites.Add(new CognateMethylationSite(null, false, type));
                continue;
            }

            bool complementary = positionText[0] == '-';
            var digits = complementary ? positionText[1..] : positionText;
            if (digits.Length == 0 || !digits.All(char.IsAsciiDigit) ||
                !int.TryParse(digits, System.Globalization.NumberStyles.None,
                    System.Globalization.CultureInfo.InvariantCulture, out int position) ||
                position < 1)
                throw new FormatException($"Malformed REBASE methylation-site position '{positionText}' in '{field}'.");

            sites.Add(new CognateMethylationSite(position, complementary, type));
        }

        return sites;
    }

    /// <summary>
    /// The second cut pair of a two-sided (Type IIB) cutter, in the same convention as
    /// <see cref="CutPositionForward"/>/<see cref="CutPositionReverse"/>; <see langword="null"/>
    /// (the default) for an enzyme that cuts once per site.
    /// </summary>
    /// <remarks>
    /// Biopython 1.88 <c>TwoCuts</c>: <c>Forward = scd5</c>, <c>Reverse = size + scd3</c>; its
    /// <c>_modify</c> yields <c>location + fst5</c> and <c>location + scd5</c>, and
    /// <c>_rev_modify</c> <c>location − fst3</c> and <c>location − scd3</c>, so every site is cut
    /// twice and the recognition site is excised as a short fragment. Example: BaeI
    /// <c>(10/15)ACNNNNGTAYC(12/7)</c> is
    /// <c>new RestrictionEnzyme("BaeI", "ACNNNNGTAYC", -10, -15, "") { SecondCut = (23, 18) }</c>.
    /// Both pairs must leave the same stagger (<c>Forward − Reverse</c> equal to
    /// <c>CutPositionForward − CutPositionReverse</c>): Biopython stores a single <c>ovhg</c> per
    /// enzyme, and all 25 two-cut enzymes of its REBASE snapshot satisfy this, so
    /// <see cref="OverhangType"/> describes both cuts. A pair violating it is rejected with
    /// <see cref="ArgumentException"/> when the enzyme is used.
    /// <para>
    /// Excluded from JSON: System.Text.Json serialises properties only, and a value tuple's
    /// <c>Item1</c>/<c>Item2</c> are fields, so it would be written as <c>{}</c> (seen in the MCP
    /// <c>get_enzyme</c> output for BaeI before B08 WP12b). The pair is serialised through
    /// <see cref="SecondCutForward"/> and <see cref="SecondCutReverse"/> instead.
    /// </para>
    /// </remarks>
    [System.Text.Json.Serialization.JsonIgnore]
    public (int Forward, int Reverse)? SecondCut { get; init; }

    /// <summary>
    /// <see cref="SecondCut"/>'s top-strand offset (Biopython <c>scd5</c>), or <see langword="null"/>
    /// for a one-cut enzyme. Read-only view that carries the pair in JSON output (BaeI: 23).
    /// </summary>
    public int? SecondCutForward => SecondCut?.Forward;

    /// <summary>
    /// <see cref="SecondCut"/>'s bottom-strand offset (Biopython <c>size + scd3</c>), or
    /// <see langword="null"/> for a one-cut enzyme. Read-only view that carries the pair in JSON
    /// output (BaeI: 18).
    /// </summary>
    public int? SecondCutReverse => SecondCut?.Reverse;

    /// <summary>
    /// Gets whether this enzyme produces blunt ends.
    /// </summary>
    public bool IsBluntEnd => CutPositionForward == CutPositionReverse;

    /// <summary>
    /// Gets the type of overhang produced.
    /// </summary>
    public OverhangType OverhangType => CutPositionForward.CompareTo(CutPositionReverse) switch
    {
        0 => OverhangType.Blunt,
        < 0 => OverhangType.FivePrime,
        _ => OverhangType.ThreePrime
    };

    /// <summary>
    /// Gets the length of the recognition sequence.
    /// </summary>
    public int RecognitionLength => RecognitionSequence.Length;
}

/// <summary>
/// The modified base type in a REBASE methylation-site entry (allenz field <c>&lt;6&gt;</c>, the
/// number in parentheses).
/// </summary>
public enum MethylatedBase
{
    /// <summary>REBASE <c>(4)</c>: N4-methylcytosine (m4C).</summary>
    N4Methylcytosine = 4,
    /// <summary>REBASE <c>(5)</c>: 5-methylcytosine (m5C).</summary>
    C5Methylcytosine = 5,
    /// <summary>REBASE <c>(6)</c>: N6-methyladenine (m6A).</summary>
    N6Methyladenine = 6
}

/// <summary>
/// One base modified by an enzyme's cognate methyltransferase, as recorded by REBASE.
/// </summary>
/// <param name="Position">1-based position of the modified base within the recognition sequence,
/// counted from the 5' end of the strand named by <paramref name="IsComplementaryStrand"/>;
/// <see langword="null"/> when REBASE gives <c>?</c> (type known, base not determined).</param>
/// <param name="IsComplementaryStrand"><see langword="true"/> for a REBASE negative position: the
/// base is on the complementary strand, numbered from that strand's own 5' base.</param>
/// <param name="Type">The modification.</param>
public sealed record CognateMethylationSite(int? Position, bool IsComplementaryStrand, MethylatedBase Type);

/// <summary>
/// Topology of a DNA molecule being digested.
/// </summary>
/// <remarks>
/// Linear and circular molecules differ in fragment count: with <c>k</c> cut sites a linear
/// molecule yields <c>k+1</c> fragments whereas a circular molecule yields <c>k</c> fragments
/// (Univ. of Illinois MolBio "Restriction Digestion/Gel Electrophoresis"; Addgene Plasmids 101).
/// Only these two values are valid: any other number cast into the enum (an MCP caller can send one) is refused by
/// every <see cref="RestrictionAnalyzer"/> member that takes a topology with an <see cref="ArgumentOutOfRangeException"/>,
/// instead of meaning "linear" in one member and "circular" in another (B08 round 5 A5R-5).
/// </remarks>
public enum MoleculeTopology
{
    /// <summary>A linear molecule with two free ends (default digest behaviour).</summary>
    Linear,
    /// <summary>A closed circular molecule (e.g. a plasmid); fragments may span the origin.</summary>
    Circular
}

/// <summary>
/// Type of overhang produced by restriction enzyme.
/// </summary>
public enum OverhangType
{
    /// <summary>Blunt end (no overhang)</summary>
    Blunt,
    /// <summary>5' overhang (sticky end)</summary>
    FivePrime,
    /// <summary>3' overhang (sticky end)</summary>
    ThreePrime
}

/// <summary>
/// Represents a restriction site in a sequence. All coordinates are 0-based and in
/// <b>forward-strand</b> space, for both strands.
/// </summary>
/// <param name="Position">Index of the 5'-most base of the recognition window on the forward
/// strand.</param>
/// <param name="Enzyme">The enzyme definition that matched, exactly as supplied.</param>
/// <param name="IsForwardStrand"><see langword="true"/> when the recognition sequence was read on
/// the forward strand, <see langword="false"/> when it was read on the reverse strand. A
/// palindromic enzyme yields one record of each per physical site, at the same
/// <paramref name="Position"/>.</param>
/// <param name="CutPosition">Number of forward-strand bases 5' of the cut the enzyme makes on the
/// strand named by <paramref name="IsForwardStrand"/>; the cut lies between forward-strand bases
/// <c>CutPosition - 1</c> and <c>CutPosition</c>. Biopython
/// <c>Bio.Restriction.&lt;enzyme&gt;.search()</c> reports a forward-strand record's (top-strand)
/// cut as <c>CutPosition + 1</c> (1-based). A reverse-strand record carries the bottom-strand cut
/// <c>Position + L − CutPositionForward</c>; its top-strand cut, which Biopython reports for the
/// reverse-strand sites of a non-palindromic enzyme, is
/// <c>CutPosition + CutPositionForward − CutPositionReverse</c> (reduced modulo the length on a
/// circle). Either cut may lie outside the recognition window for a distal (Type IIS) cutter.</param>
/// <param name="RecognizedSequence">The matched window as read on the strand named by
/// <paramref name="IsForwardStrand"/> — so for a reverse-strand record it is the reverse
/// complement of the forward-strand window. Always
/// <c>Enzyme.RecognitionSequence.Length</c> characters long.</param>
/// <remarks>
/// A two-sided (Type IIB) enzyme — one with a <see cref="RestrictionEnzyme.SecondCut"/> — cuts every
/// site twice and is reported as <b>two</b> records per site and strand, one per cut pair, at the
/// same <paramref name="Position"/>; the record of the second pair has
/// <see cref="IsSecondCut"/> set and its <paramref name="CutPosition"/> follows the same rule with
/// the second pair's offsets. On a linear molecule each cut is kept or dropped on its own (Biopython
/// <c>_drop</c> filters cut by cut), so a site near an end may report only one of its two cuts.
/// </remarks>
public sealed record RestrictionSite(
    int Position,
    RestrictionEnzyme Enzyme,
    bool IsForwardStrand,
    int CutPosition,
    string RecognizedSequence)
{
    /// <summary>
    /// <see langword="true"/> when this record carries the cut of the enzyme's
    /// <see cref="RestrictionEnzyme.SecondCut"/> pair; always <see langword="false"/> for an enzyme
    /// that cuts once per site.
    /// </summary>
    public bool IsSecondCut { get; init; }

    /// <summary>
    /// The (forward, reverse) cut offsets this record's cut was made with: the enzyme's
    /// <see cref="RestrictionEnzyme.CutPositionForward"/>/<see cref="RestrictionEnzyme.CutPositionReverse"/>,
    /// or its <see cref="RestrictionEnzyme.SecondCut"/> pair when <see cref="IsSecondCut"/>.
    /// </summary>
    /// <remarks>
    /// Excluded from JSON (B08 WP12b): a value tuple's <c>Item1</c>/<c>Item2</c> are fields, which
    /// System.Text.Json does not write, so every MCP site record carried an empty
    /// <c>"cutOffsets":{}</c>. The value is derived — <see cref="Enzyme"/>'s offsets, or its
    /// <see cref="RestrictionEnzyme.SecondCutForward"/>/<see cref="RestrictionEnzyme.SecondCutReverse"/>
    /// when <see cref="IsSecondCut"/> — so nothing is lost.
    /// </remarks>
    [System.Text.Json.Serialization.JsonIgnore]
    public (int Forward, int Reverse) CutOffsets =>
        IsSecondCut && Enzyme.SecondCut is { } second
            ? second
            : (Enzyme.CutPositionForward, Enzyme.CutPositionReverse);
}

/// <summary>
/// Represents a DNA fragment after restriction digestion.
/// </summary>
/// <param name="Sequence">The fragment's <b>top-strand</b> bases, 5'→3'. A digest partitions the
/// top strand at the enzymes' top-strand (watson) cuts, exactly as Biopython 1.88
/// <c>Bio.Restriction</c> <c>enzyme.catalyse(seq)</c> slices at <c>enzyme.search(seq)</c>; the
/// bases of a sticky end therefore belong to the fragment downstream of the top-strand cut, and
/// the end's single-stranded extension is reported separately below.</param>
/// <param name="StartPosition">0-based index, in the original molecule, of the fragment's first
/// top-strand base. For the origin-spanning fragment of a circular digest this is the last cut
/// position, and the fragment continues past the origin.</param>
/// <param name="Length">Fragment length in bases. The lengths of all fragments of one digest sum
/// to the molecule's length in both topologies (Addgene Restriction Digest Protocol).</param>
/// <param name="LeftEnzyme">Name of the enzyme whose cut created the fragment's left (5'-side)
/// end, or <see langword="null"/> when that end is an original terminus of a linear molecule or
/// the fragment is an uncut circle. When several enzymes cut at the same top-strand coordinate
/// (e.g. the isoschizomers MspI/HpaII, or BamHI with MboI/Sau3AI) this is the <b>first</b> name of
/// <see cref="LeftEnzymes"/> — a deterministic, documented choice (ordinal order of the table
/// names), not the only enzyme; read <see cref="LeftEnzymes"/> for the complete set.</param>
/// <param name="RightEnzyme">Name of the enzyme whose cut created the fragment's right (3'-side)
/// end, under the same convention as <paramref name="LeftEnzyme"/>.</param>
/// <param name="FragmentNumber">1-based sequential number in the returned order (ascending
/// <paramref name="StartPosition"/>), contiguous over the digest.</param>
/// <param name="LeftOverhangType">End chemistry of the left end, as made by the
/// <paramref name="LeftEnzyme"/> cut alone — <see cref="OverhangType.Blunt"/>
/// for a blunt cutter, <see cref="OverhangType.FivePrime"/> or
/// <see cref="OverhangType.ThreePrime"/> for a staggered cut — or <see langword="null"/> when the
/// end was not created by a cut.</param>
/// <param name="LeftOverhang">The left end's single-stranded extension, written 5'→3' on the
/// strand that carries it and taken from the molecule itself: empty for a blunt end,
/// <see langword="null"/> when the end was not created by a cut. The extension spans the
/// <c>|CutPositionReverse − CutPositionForward|</c> bases between the cut's two strand cuts; over
/// the recognition sequence this window equals Biopython's <c>enzyme.ovhgseq</c> for every one of
/// the 39 built-in enzymes (e.g. EcoRI <c>AATT</c> 5', PstI <c>TGCA</c> 3').</param>
/// <param name="RightOverhangType">End chemistry of the right end, as
/// <paramref name="LeftOverhangType"/>.</param>
/// <param name="RightOverhang">The right end's single-stranded extension, as
/// <paramref name="LeftOverhang"/>. At one cut the upstream fragment's right overhang and the
/// downstream fragment's left overhang are reverse complements of each other — which is why they
/// re-ligate; for the palindromic built-in enzymes they are the same string.</param>
/// <remarks>
/// <b>Coinciding cuts (review 2026-09 B08 F48).</b> Every built-in collision but one joins enzymes
/// with identical ends (MspI/HpaII <c>CG</c> 5', MboI/Sau3AI/BamHI/BglII <c>GATC</c> 5',
/// ClaI/TaqI <c>CG</c> 5', StuI/HaeIII and ScaI/RsaI blunt — Biopython 1.88 <c>ovhgseq</c>), so the
/// reported end is the end of each of them. The exception is <b>FseI/HaeIII</b>: inside
/// <c>GGCCGGCC</c> both cut the top strand after base 6, but FseI leaves a 4-nt 3' overhang
/// (<c>CCGG</c>, Biopython <c>ovhg = 4</c>) while HaeIII cuts blunt. The overhang fields then report
/// the representative's end (FseI, ordinally first), i.e. the end of a digest by that enzyme
/// alone. The physical end of the combined digest also depends on bottom-strand nicks made at
/// other boundaries (here HaeIII's other cut, which coincides with FseI's bottom-strand cut) and is
/// not modelled; when <see cref="LeftEnzymes"/>/<see cref="RightEnzymes"/> name more than one
/// enzyme, check their individual ends via <see cref="RestrictionAnalyzer.GetEnzyme"/>(name).<see cref="RestrictionEnzyme.OverhangType"/> and <see cref="RestrictionAnalyzer.GetOverhangSequence(string)"/>.
/// </remarks>
public sealed record DigestFragment(
    string Sequence,
    int StartPosition,
    int Length,
    string? LeftEnzyme,
    string? RightEnzyme,
    int FragmentNumber,
    OverhangType? LeftOverhangType = null,
    string? LeftOverhang = null,
    OverhangType? RightOverhangType = null,
    string? RightOverhang = null)
{
    /// <summary>
    /// Every enzyme whose top-strand cut created the fragment's left end, ordered by
    /// <see cref="StringComparer.Ordinal"/> on the enzyme-table name and listed once each; empty
    /// when that end is not a cut. Its first element is <see cref="LeftEnzyme"/>. When not set
    /// explicitly it is <see cref="LeftEnzyme"/> alone (or empty when that is null).
    /// </summary>
    public IReadOnlyList<string> LeftEnzymes { get; init; } = LeftEnzyme is null ? [] : [LeftEnzyme];

    /// <summary>
    /// Every enzyme whose top-strand cut created the fragment's right end, under the convention
    /// of <see cref="LeftEnzymes"/>. Its first element is <see cref="RightEnzyme"/>.
    /// </summary>
    public IReadOnlyList<string> RightEnzymes { get; init; } = RightEnzyme is null ? [] : [RightEnzyme];

    /// <summary>
    /// Value equality over every member, comparing <see cref="LeftEnzymes"/> and
    /// <see cref="RightEnzymes"/> element-wise (a record would otherwise compare the lists by
    /// reference, so two identical digests would stop being equal).
    /// </summary>
    public bool Equals(DigestFragment? other)
        => other is not null
           && Sequence == other.Sequence
           && StartPosition == other.StartPosition
           && Length == other.Length
           && LeftEnzyme == other.LeftEnzyme
           && RightEnzyme == other.RightEnzyme
           && FragmentNumber == other.FragmentNumber
           && LeftOverhangType == other.LeftOverhangType
           && LeftOverhang == other.LeftOverhang
           && RightOverhangType == other.RightOverhangType
           && RightOverhang == other.RightOverhang
           && LeftEnzymes.SequenceEqual(other.LeftEnzymes)
           && RightEnzymes.SequenceEqual(other.RightEnzymes);

    /// <inheritdoc/>
    public override int GetHashCode()
        => HashCode.Combine(
            HashCode.Combine(Sequence, StartPosition, Length, LeftEnzyme, RightEnzyme, FragmentNumber),
            HashCode.Combine(LeftOverhangType, LeftOverhang, RightOverhangType, RightOverhang),
            LeftEnzymes.Count,
            RightEnzymes.Count);
}

/// <summary>
/// Summary of a restriction digest.
/// </summary>
/// <param name="TotalFragments">Number of fragments: <c>k+1</c> for a linear molecule with
/// <c>k</c> cut positions, <c>k</c> for a circular one.</param>
/// <param name="FragmentSizes">Fragment lengths in <b>descending</b> order (the
/// gel-electrophoresis reading order). They sum to the digested molecule's length.</param>
/// <param name="LargestFragment">Largest fragment length, or 0 when there are no fragments.</param>
/// <param name="SmallestFragment">Smallest fragment length, or 0 when there are no fragments.</param>
/// <param name="AverageFragmentSize">Mean fragment length, or 0 when there are no fragments.</param>
/// <param name="EnzymesUsed">The enzyme names as supplied by the caller, in order.</param>
/// <param name="Topology">The topology the digest was run with; <see cref="MoleculeTopology.Linear"/>
/// for the topology-free overload.</param>
public sealed record DigestSummary(
    int TotalFragments,
    IReadOnlyList<int> FragmentSizes,
    int LargestFragment,
    int SmallestFragment,
    double AverageFragmentSize,
    IReadOnlyList<string> EnzymesUsed,
    MoleculeTopology Topology = MoleculeTopology.Linear);

/// <summary>
/// A restriction map of a sequence.
/// </summary>
public sealed record RestrictionMap(
    int SequenceLength,
    IReadOnlyList<RestrictionSite> Sites,
    IReadOnlyDictionary<string, List<int>> SitesByEnzyme,
    int TotalSites,
    IReadOnlyList<string> UniqueCutters,
    IReadOnlyList<string> NonCutters);
