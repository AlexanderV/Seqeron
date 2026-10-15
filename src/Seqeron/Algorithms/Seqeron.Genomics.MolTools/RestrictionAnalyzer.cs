namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Analyzes restriction enzyme sites in DNA sequences.
/// Supports common restriction enzymes with palindromic recognition sequences.
/// </summary>
public static class RestrictionAnalyzer
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
    /// Gets all available restriction enzymes.
    /// </summary>
    public static IReadOnlyDictionary<string, RestrictionEnzyme> Enzymes => _enzymes;

    /// <summary>
    /// Gets a restriction enzyme by name, resolved case-insensitively, or <see langword="null"/>
    /// when the name is not in the built-in table.
    /// </summary>
    /// <param name="name">Enzyme name (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="name"/> is <see langword="null"/>.</exception>
    public static RestrictionEnzyme? GetEnzyme(string name)
    {
        ArgumentNullException.ThrowIfNull(name);
        return _enzymes.TryGetValue(name, out var enzyme) ? enzyme : null;
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
        _enzymes.Values.Where(e => e.RecognitionSequence.Length == length);

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
        _enzymes.Values.Where(e =>
            e.RecognitionSequence.Length >= minLength &&
            e.RecognitionSequence.Length <= maxLength);

    /// <summary>
    /// Gets all enzymes that produce blunt ends — both strands are cut at the same position
    /// so each strand "terminates in a base pair" (Wikipedia, "Sticky and blunt ends").
    /// </summary>
    /// <remarks>
    /// Equivalent to Biopython <c>Bio.Restriction.&lt;enzyme&gt;.is_blunt()</c>; all 39 built-in
    /// entries were audited against it (39/39 agreement, 10 blunt cutters).
    /// </remarks>
    public static IEnumerable<RestrictionEnzyme> GetBluntCutters() =>
        _enzymes.Values.Where(e => e.IsBluntEnd);

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
        _enzymes.Values.Where(e => !e.IsBluntEnd);

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
    // • Feasibility (linear molecules): a cut can only fall on a phosphodiester bond that
    //   exists, so a site is reported only when BOTH its strand cuts lie in [1, length - 1].
    //   This is Biopython's `_drop` filter (`1 < cut_on_watson <= length` and the same for
    //   `cut_on_crick = cut_on_watson - ovhg`, with `ovhg = CutPositionForward -
    //   CutPositionReverse`) expressed in 0-based coordinates. It only ever excludes an
    //   enzyme that cuts at the very edge of its own recognition window — in the built-in
    //   table, Sau3AI and MboI (^GATC / GATC^).
    //
    // • Circular molecules: the search runs over `seq + seq[0 .. L-2]` so sites spanning the
    //   origin are found (Biopython appends `data[1:size]`), cut positions are taken modulo
    //   the length, and no site is dropped — on a closed circle every bond exists.
    // ───────────────────────────────────────────────────────────────────────────────────────

    /// <summary>
    /// Finds all restriction sites for a specific enzyme in a <b>linear</b> sequence.
    /// </summary>
    /// <param name="sequence">DNA sequence to analyze.</param>
    /// <param name="enzymeName">Name of the restriction enzyme (case-insensitive).</param>
    /// <returns>Collection of restriction sites found: forward-strand hits first, then
    /// reverse-strand hits, each in ascending <see cref="RestrictionSite.Position"/> order.</returns>
    /// <remarks>
    /// Equivalent to <see cref="FindSites(DnaSequence, string, MoleculeTopology)"/> with
    /// <see cref="MoleculeTopology.Linear"/>; see that overload for the coordinate convention.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null, or
    /// <paramref name="enzymeName"/> is null or empty.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is not in the built-in table.</exception>
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
    /// single site at position 4 (cut 4), exactly as
    /// <c>Bio.Restriction.Sau3AI.search(Seq("GATCGATCGATC")) == [5]</c> (1-based).
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
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, string enzymeName, MoleculeTopology topology)
    {
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
    public static IEnumerable<RestrictionSite> FindSites(
        string sequence, string enzymeName, MoleculeTopology topology)
    {
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
    /// that is not one of the 15 IUPAC nucleotide codes.
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
    /// that is not one of the 15 IUPAC nucleotide codes.
    /// </exception>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, RestrictionEnzyme enzyme, MoleculeTopology topology)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzyme);

        return FindSitesCore(sequence.Sequence, enzyme, topology);
    }

    /// <summary>
    /// Finds all restriction sites for multiple enzymes in a <b>linear</b> sequence.
    /// </summary>
    /// <remarks>The result is the concatenation of the per-enzyme searches, in the order given;
    /// a name repeated in <paramref name="enzymeNames"/> is searched (and reported) twice.</remarks>
    public static IEnumerable<RestrictionSite> FindSites(DnaSequence sequence, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzymeNames);
        return FindSitesCore(sequence, enzymeNames, MoleculeTopology.Linear);
    }

    /// <summary>
    /// Finds all restriction sites for multiple enzymes in a molecule of the given topology.
    /// </summary>
    /// <remarks>The result is the concatenation of the per-enzyme searches, in the order given;
    /// a name repeated in <paramref name="enzymeNames"/> is searched (and reported) twice.</remarks>
    public static IEnumerable<RestrictionSite> FindSites(
        DnaSequence sequence, MoleculeTopology topology, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzymeNames);
        return FindSitesCore(sequence, enzymeNames, topology);
    }

    private static IEnumerable<RestrictionSite> FindSitesCore(
        DnaSequence sequence, string[] enzymeNames, MoleculeTopology topology)
    {
        foreach (var name in enzymeNames)
        {
            foreach (var site in FindSites(sequence, name, topology))
                yield return site;
        }
    }

    /// <summary>
    /// Finds all sites for all known enzymes in a <b>linear</b> sequence (comprehensive analysis).
    /// </summary>
    public static IEnumerable<RestrictionSite> FindAllSites(DnaSequence sequence)
        => FindAllSites(sequence, MoleculeTopology.Linear);

    /// <summary>
    /// Finds all sites for all known enzymes in a molecule of the given topology
    /// (comprehensive analysis).
    /// </summary>
    public static IEnumerable<RestrictionSite> FindAllSites(DnaSequence sequence, MoleculeTopology topology)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        return FindAllSitesCore(sequence, topology);
    }

    private static IEnumerable<RestrictionSite> FindAllSitesCore(DnaSequence sequence, MoleculeTopology topology)
    {
        foreach (var enzyme in _enzymes.Values)
        {
            foreach (var site in FindSitesCore(sequence.Sequence, enzyme, topology))
                yield return site;
        }
    }

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
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="recognitionPattern"/> is null or empty,
    /// or contains a character that is not an IUPAC nucleotide code.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="startIndex"/> is negative.</exception>
    public static IEnumerable<int> FindPatternOccurrences(
        string sequence, string recognitionPattern, int startIndex = 0)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentOutOfRangeException.ThrowIfNegative(startIndex);
        string pattern = ValidatedPattern(recognitionPattern, nameof(recognitionPattern));

        // ASCII-only fold (review 2026-09 B08 F49, superseding F43): ToUpperInvariant would map U+017F 'ſ'
        // onto the IUPAC code 'S', so a non-IUPAC character would match.
        return sequence.Length == 0
            ? []
            : EnumeratePatternOccurrences(SequenceExtensions.ToUpperAscii(sequence), pattern, startIndex);
    }

    private static IEnumerable<int> EnumeratePatternOccurrences(string seq, string pattern, int startIndex)
    {
        for (int i = startIndex; i + pattern.Length <= seq.Length; i++)
        {
            if (MatchesPatternAt(seq, i, pattern))
                yield return i;
        }
    }

    /// <summary>
    /// Resolves a built-in enzyme name or throws — the single name-resolution contract of this
    /// class, shared by the site finder, the digest, the map and the compatibility members.
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <param name="paramName">Name of the caller's parameter to blame, when it is not
    /// <c>enzymeName</c>.</param>
    private static RestrictionEnzyme RequireEnzyme(string enzymeName, string? paramName = null)
    {
        paramName ??= nameof(enzymeName);

        if (string.IsNullOrEmpty(enzymeName))
            throw new ArgumentNullException(paramName);

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
        string seq, RestrictionEnzyme enzyme, MoleculeTopology topology)
    {
        string pattern = ValidatedPattern(enzyme.RecognitionSequence, nameof(enzyme), enzyme.Name);

        if (seq.Length < pattern.Length)
            return [];

        return topology == MoleculeTopology.Circular
            ? FindSitesCircular(seq, enzyme, pattern)
            : FindSitesLinear(seq, enzyme, pattern);
    }

    private static IEnumerable<RestrictionSite> FindSitesLinear(
        string seq, RestrictionEnzyme enzyme, string pattern)
    {
        int n = seq.Length;
        int len = pattern.Length;
        int fwd = enzyme.CutPositionForward;
        int rev = enzyme.CutPositionReverse;

        for (int i = 0; i + len <= n; i++)
        {
            // Forward-strand site: top-strand cut i+fwd, bottom-strand cut i+rev.
            if (MatchesPatternAt(seq, i, pattern) && CutsAreInternal(i + fwd, i + rev, n))
            {
                yield return new RestrictionSite(
                    Position: i,
                    Enzyme: enzyme,
                    IsForwardStrand: true,
                    CutPosition: i + fwd,
                    RecognizedSequence: seq.Substring(i, len));
            }
        }

        // Check reverse strand (complement sequence)
        string revComp = DnaSequence.GetReverseComplementString(seq);
        for (int i = 0; i + len <= n; i++)
        {
            if (!MatchesPatternAt(revComp, i, pattern))
                continue;

            // Convert back to forward strand coordinates
            int forwardPos = n - i - len;

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
                RecognizedSequence: revComp.Substring(i, len));
        }
    }

    private static IEnumerable<RestrictionSite> FindSitesCircular(
        string seq, RestrictionEnzyme enzyme, string pattern)
    {
        int n = seq.Length;
        int len = pattern.Length;
        int fwd = enzyme.CutPositionForward;

        // Biopython searches data + data[1:size] so a site may wrap past the origin; a start
        // index is therefore considered exactly once, for 0 <= i < n.
        string extended = seq + seq[..(len - 1)];
        for (int i = 0; i < n; i++)
        {
            if (MatchesPatternAt(extended, i, pattern))
            {
                yield return new RestrictionSite(
                    Position: i,
                    Enzyme: enzyme,
                    IsForwardStrand: true,
                    CutPosition: (i + fwd) % n,
                    RecognizedSequence: extended.Substring(i, len));
            }
        }

        string revComp = DnaSequence.GetReverseComplementString(seq);
        string revExtended = revComp + revComp[..(len - 1)];
        for (int i = 0; i < n; i++)
        {
            if (!MatchesPatternAt(revExtended, i, pattern))
                continue;

            int forwardPos = ((n - i - len) % n + n) % n;

            yield return new RestrictionSite(
                Position: forwardPos,
                Enzyme: enzyme,
                IsForwardStrand: false,
                CutPosition: (forwardPos + len - fwd) % n,
                RecognizedSequence: revExtended.Substring(i, len));
        }
    }

    /// <summary>
    /// True when both strand cuts of a site fall on a phosphodiester bond that exists in a linear
    /// molecule of <paramref name="length"/> bases, i.e. both lie in <c>[1, length - 1]</c>
    /// (Biopython <c>AbstractCut._drop</c>, expressed in 0-based coordinates).
    /// </summary>
    private static bool CutsAreInternal(int watsonCut, int crickCut, int length) =>
        watsonCut >= 1 && watsonCut <= length - 1 &&
        crickCut >= 1 && crickCut <= length - 1;

    private static bool MatchesPatternAt(string seq, int position, string pattern)
    {
        for (int i = 0; i < pattern.Length; i++)
        {
            if (!SequenceBaseMatchesSiteCode(seq[position + i], pattern[i]))
                return false;
        }
        return true;
    }

    /// <summary>
    /// The per-position site-match rule of Biopython 1.88 <c>Bio.Restriction</c>, composed from
    /// the canonical <see cref="IupacHelper"/> predicates (no second ambiguity table here).
    /// </summary>
    /// <remarks>
    /// Biopython compiles a site's <c>N</c> to the regex <c>.</c> and every other code to the
    /// class of its A/C/G/T expansion (<c>SfiI.compsite</c> is <c>GGCC.....GGCC</c>,
    /// <c>HincII.compsite</c> is <c>GT[CT][AG]AC</c>), and the sequence keeps its own IUPAC
    /// letters. So a site <c>N</c> accepts <b>any</b> sequence code — an <c>N</c> included:
    /// <c>SfiI.search(Seq("GGCCNNNNNGGCC")) == [9]</c> — while a site A/C/G/T or a degenerate
    /// site code accepts only a concrete base it expands to: an ambiguity code in the sequence
    /// still blocks it (<c>EcoRI.search(Seq("AGAANTCA")) == []</c>,
    /// <c>HincII.search(Seq("AGTYRACA")) == []</c>). A sequence character outside the 15 IUPAC
    /// nucleotide codes matches nothing, not even a site <c>N</c> (Biopython would let a site
    /// <c>N</c> match any letter, <c>X</c> included — a deliberate narrowing).
    /// </remarks>
    private static bool SequenceBaseMatchesSiteCode(char sequenceBase, char siteCode) =>
        siteCode == 'N'
            ? IupacHelper.IsNucleotideCode(sequenceBase)
            : IupacHelper.MatchesIupac(sequenceBase, siteCode);

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
    //   (sticky-end) cut. This class therefore keeps only `IsForwardStrand` site records — whose
    //   `CutPosition` is the top-strand cut — and deliberately drops the reverse-strand record of
    //   the same physical site, which carries the BOTTOM-strand cut. For the 39 built-in enzymes,
    //   all palindromic, the two records are the two cuts of one site, so dropping the reverse one
    //   loses no site; the bottom-strand cut is not lost either — it is reported as the fragment
    //   end's overhang (below).
    //
    // • Fragment end chemistry. A cut by an enzyme with forward cut `f` and reverse cut `r` has
    //   its top-strand cut at `watson` and its bottom-strand cut at `crick = watson + (r - f)`.
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
    public static IEnumerable<DigestFragment> Digest(DnaSequence sequence, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (enzymeNames == null || enzymeNames.Length == 0)
            throw new ArgumentException("At least one enzyme is required", nameof(enzymeNames));
        return DigestCore(sequence, enzymeNames);
    }

    /// <summary>
    /// Collects the distinct top-strand (forward-strand) cut positions for the requested enzymes,
    /// together with every enzyme's site at each position.
    /// </summary>
    /// <remarks>
    /// Only <see cref="RestrictionSite.IsForwardStrand"/> records are kept: their
    /// <see cref="RestrictionSite.CutPosition"/> is the top-strand cut, which is what Biopython's
    /// <c>catalyse</c> slices at. The reverse-strand record of the same palindromic site carries
    /// the bottom-strand cut and would double-count the site as a second boundary. Several
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
        DnaSequence sequence, string[] enzymeNames, MoleculeTopology topology)
    {
        var cutPositions = new SortedSet<int>();
        var sitesByPosition = new Dictionary<int, List<RestrictionSite>>();

        foreach (var name in enzymeNames)
        {
            foreach (var site in FindSites(sequence, name, topology))
            {
                if (!site.IsForwardStrand)
                    continue;

                // FindSites already reports a circular cut reduced modulo the sequence length
                // (FindSitesCircular yields (i + fwd) % n) and a linear cut inside [1, Length-1]
                // (CutsAreInternal), so no further normalisation is needed here.
                cutPositions.Add(site.CutPosition);
                if (!sitesByPosition.TryGetValue(site.CutPosition, out var atCut))
                    sitesByPosition[site.CutPosition] = atCut = [];
                if (!atCut.Exists(s => s.Enzyme.Name == site.Enzyme.Name))
                    atCut.Add(site);
            }
        }

        foreach (var atCut in sitesByPosition.Values)
            atCut.Sort((a, b) => string.CompareOrdinal(a.Enzyme.Name, b.Enzyme.Name));

        return (cutPositions.ToList(), sitesByPosition);
    }

    private static IEnumerable<DigestFragment> DigestCore(DnaSequence sequence, string[] enzymeNames)
    {
        var (cuts, sitesByPosition) = CollectForwardCuts(sequence, enzymeNames, MoleculeTopology.Linear);
        string seq = sequence.Sequence;
        int n = sequence.Length;

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
    public static IEnumerable<DigestFragment> Digest(
        DnaSequence sequence,
        MoleculeTopology topology,
        params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (enzymeNames == null || enzymeNames.Length == 0)
            throw new ArgumentException("At least one enzyme is required", nameof(enzymeNames));

        if (topology == MoleculeTopology.Linear)
            return Digest(sequence, enzymeNames);

        return DigestCircular(sequence, enzymeNames);
    }

    private static IEnumerable<DigestFragment> DigestCircular(DnaSequence sequence, string[] enzymeNames)
    {
        // Circular site search: a recognition window may wrap past the origin, so the cut it makes
        // must be seen too. Using the linear search here would silently leave an origin-spanning
        // site uncut and additionally drop an edge cutter's cuts through the linear-only
        // cut-feasibility filter, which does not apply to a closed circle.
        var (cuts, sitesByPosition) = CollectForwardCuts(sequence, enzymeNames, MoleculeTopology.Circular);

        int n = sequence.Length;
        string seq = sequence.Sequence;

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
            OverhangType.ThreePrime => (type, DnaSequence.GetReverseComplementString(window)),
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
            OverhangType.FivePrime => (type, DnaSequence.GetReverseComplementString(window)),
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
    /// via the shared <see cref="OverhangWindow"/> geometry — off the MOLECULE instead, translated
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

        var (patternStart, length) = OverhangWindow(enzyme);
        int lo = site.CutPosition + (patternStart - enzyme.CutPositionForward);
        int n = seq.Length;

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
    public static DigestSummary GetDigestSummary(
        DnaSequence sequence,
        MoleculeTopology topology,
        params string[] enzymeNames)
    {
        var fragments = Digest(sequence, topology, enzymeNames).ToList();
        var sizes = fragments.Select(f => f.Length).OrderByDescending(x => x).ToList();

        return new DigestSummary(
            TotalFragments: fragments.Count,
            FragmentSizes: sizes,
            LargestFragment: sizes.FirstOrDefault(),
            SmallestFragment: sizes.LastOrDefault(),
            AverageFragmentSize: sizes.Count > 0 ? sizes.Average() : 0,
            EnzymesUsed: enzymeNames.ToList(),
            Topology: topology);
    }

    #endregion

    #region Restriction Map

    /// <summary>
    /// Creates a restriction map of a <b>linear</b> molecule showing all enzyme sites.
    /// </summary>
    public static RestrictionMap CreateMap(DnaSequence sequence, params string[] enzymeNames)
        => CreateMap(sequence, MoleculeTopology.Linear, enzymeNames);

    /// <summary>
    /// Creates a restriction map of a molecule with the given topology, showing all enzyme sites.
    /// </summary>
    /// <param name="sequence">DNA sequence to map.</param>
    /// <param name="topology">Whether the molecule is linear or circular; see
    /// <see cref="FindSites(DnaSequence, string, MoleculeTopology)"/>.</param>
    /// <param name="enzymeNames">Enzymes to consider. Names are resolved case-insensitively and
    /// de-duplicated; an <b>empty</b> set means every built-in enzyme.</param>
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
    /// positions of both strands' records; <c>TotalSites</c> and <c>UniqueCutters</c> count
    /// forward-strand records only, since a palindromic site yields one of each.
    /// </para>
    /// </remarks>
    /// <exception cref="ArgumentException">A name in <paramref name="enzymeNames"/> is not in the
    /// built-in table.</exception>
    public static RestrictionMap CreateMap(
        DnaSequence sequence, MoleculeTopology topology, params string[] enzymeNames)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        ArgumentNullException.ThrowIfNull(enzymeNames);

        // Resolve the request to canonical, de-duplicated enzyme names so that NonCutters can be
        // compared against the keys of SitesByEnzyme (which carry the table's own spelling) and
        // so a name repeated with different casing is not searched — and counted — twice.
        var considered = (enzymeNames.Length > 0
                ? enzymeNames.Select(n => RequireEnzyme(n))
                : _enzymes.Values)
            .Select(e => e.Name)
            .Distinct(StringComparer.OrdinalIgnoreCase)
            .ToList();

        var sites = FindSitesCore(sequence, considered.ToArray(), topology).ToList();

        var byEnzyme = sites
            .GroupBy(s => s.Enzyme.Name, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(
                g => g.Key,
                g => g.Select(s => s.Position).OrderBy(p => p).ToList(),
                StringComparer.OrdinalIgnoreCase);

        // UniqueCutters: enzymes with exactly one forward-strand site.
        // SitesByEnzyme includes both strands; for palindromic enzymes a single site
        // produces two entries (forward + reverse) at the same position.
        // Filter to forward-strand sites for unique cutter identification.
        var forwardByEnzyme = sites
            .Where(s => s.IsForwardStrand)
            .GroupBy(s => s.Enzyme.Name, StringComparer.OrdinalIgnoreCase)
            .ToDictionary(
                g => g.Key,
                g => g.Select(s => s.Position).Distinct().ToList(),
                StringComparer.OrdinalIgnoreCase);
        var uniqueCutters = forwardByEnzyme.Where(kv => kv.Value.Count == 1).Select(kv => kv.Key).ToList();
        var nonCutters = considered.Where(e => !byEnzyme.ContainsKey(e)).ToList();

        return new RestrictionMap(
            SequenceLength: sequence.Length,
            Sites: sites,
            SitesByEnzyme: byEnzyme,
            TotalSites: sites.Count(s => s.IsForwardStrand), // Count only forward strand
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
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is not a built-in enzyme.</exception>
    public static string GetOverhangSequence(string enzymeName)
        => GetOverhangSequence(RequireEnzyme(enzymeName))
           ?? throw new ArgumentException(
               $"Enzyme '{enzymeName}' cuts outside its own recognition sequence, so its overhang " +
               "is not determined by the recognition pattern.", nameof(enzymeName));

    /// <summary>
    /// The single-stranded overhang <paramref name="enzyme"/> leaves, written 5'→3' on the top
    /// strand — REBASE/Biopython <c>ovhgseq</c>.
    /// </summary>
    /// <param name="enzyme">Enzyme definition; may be a caller-supplied (non-built-in) enzyme.</param>
    /// <returns><see cref="string.Empty"/> for a blunt cutter; the recognition-pattern bases
    /// between the two strand cuts for a staggered cutter; or <see langword="null"/> when the
    /// overhang window is not contained in the recognition sequence, which happens only for a
    /// cutter that cuts outside its own site (a Type IIS-style distal cut). In that case the
    /// overhang depends on the flanking molecule and no enzyme-level answer exists.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="enzyme"/> is <see langword="null"/>.</exception>
    public static string? GetOverhangSequence(RestrictionEnzyme enzyme)
    {
        ArgumentNullException.ThrowIfNull(enzyme);

        var (start, length) = OverhangWindow(enzyme);
        if (length == 0)
            return string.Empty;
        if (start < 0 || start + length > enzyme.RecognitionSequence.Length)
            return null;

        return enzyme.RecognitionSequence.Substring(start, length);
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
        (Math.Min(enzyme.CutPositionForward, enzyme.CutPositionReverse),
         Math.Abs(enzyme.CutPositionForward - enzyme.CutPositionReverse));

    /// <summary>
    /// Enumerates every unordered pair of distinct built-in enzymes whose ends can be ligated to
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
    /// </remarks>
    public static IEnumerable<(string Enzyme1, string Enzyme2, string CompatibleEnd)> FindCompatibleEnzymes()
    {
        var enzymes = _enzymes.Values.ToList();

        for (int i = 0; i < enzymes.Count; i++)
        {
            for (int j = i + 1; j < enzymes.Count; j++)
            {
                var e1 = enzymes[i];
                var e2 = enzymes[j];

                if (AreCompatible(e1, e2, out string? compatibleEnd))
                {
                    yield return (e1.Name, e2.Name, compatibleEnd!);
                }
            }
        }
    }

    /// <summary>
    /// Lists the built-in enzymes whose ends can be ligated to an end left by
    /// <paramref name="enzymeName"/>, the enzyme itself included.
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <returns>Tuples of the compatible enzyme's name (as spelled in the table) and the shared
    /// end — <c>"blunt"</c> or the overhang sequence.</returns>
    /// <remarks>
    /// This is the per-enzyme form of Biopython's <c>compatible_end()</c>. The reference's own
    /// implementation filters its candidates with <c>x.is_5overhang()</c>, so it answers the empty
    /// list for every blunt and every 3'-overhang enzyme; this overload has no such restriction and
    /// is simply <see cref="AreCompatible(RestrictionEnzyme, RestrictionEnzyme)"/> over the table.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is not a built-in enzyme.</exception>
    public static IEnumerable<(string Enzyme, string CompatibleEnd)> FindCompatibleEnzymes(string enzymeName)
    {
        var query = RequireEnzyme(enzymeName);
        return Enumerate();

        IEnumerable<(string, string)> Enumerate()
        {
            foreach (var candidate in _enzymes.Values)
            {
                if (AreCompatible(query, candidate, out string? compatibleEnd))
                    yield return (candidate.Name, compatibleEnd!);
            }
        }
    }

    /// <summary>
    /// Checks whether two named built-in enzymes produce ligatable (compatible) ends.
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
    /// <exception cref="ArgumentException">Either name is not a built-in enzyme.</exception>
    public static bool AreCompatible(string enzyme1Name, string enzyme2Name)
    {
        var e1 = RequireEnzyme(enzyme1Name, nameof(enzyme1Name));
        var e2 = RequireEnzyme(enzyme2Name, nameof(enzyme2Name));

        return AreCompatible(e1, e2, out _);
    }

    /// <summary>
    /// Lenient form of <see cref="AreCompatible(string, string)"/> for callers that cannot
    /// guarantee both names are in the built-in table.
    /// </summary>
    /// <param name="enzyme1Name">First enzyme name (case-insensitive).</param>
    /// <param name="enzyme2Name">Second enzyme name (case-insensitive).</param>
    /// <param name="compatible">The compatibility verdict when both names resolved; otherwise
    /// <see langword="false"/>.</param>
    /// <returns><see langword="true"/> when both names were resolved, so
    /// <paramref name="compatible"/> is meaningful; <see langword="false"/> when either name is
    /// unknown.</returns>
    /// <exception cref="ArgumentNullException">Either name is <see langword="null"/>.</exception>
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
    /// region). <see langword="false"/> when either enzyme cuts outside its own recognition
    /// sequence, since there is then no enzyme-level overhang to compare — Biopython's
    /// <c>NotDefined</c> ends are likewise "compatible with nobody".</returns>
    /// <exception cref="ArgumentNullException">Either enzyme is <see langword="null"/>.</exception>
    public static bool AreCompatible(RestrictionEnzyme enzyme1, RestrictionEnzyme enzyme2)
    {
        ArgumentNullException.ThrowIfNull(enzyme1);
        ArgumentNullException.ThrowIfNull(enzyme2);

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
            || CodesIntersect(overhang1, DnaSequence.GetReverseComplementString(overhang2));
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
    /// Built-in enzymes that recognise the same sequence as <paramref name="enzymeName"/> but are
    /// not that enzyme — its isoschizomers (Biopython <c>isoschizomers()</c>; the union of
    /// <see cref="FindEquischizomers"/> and <see cref="FindNeoschizomers"/>).
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is not a built-in enzyme.</exception>
    public static IEnumerable<RestrictionEnzyme> FindIsoschizomers(string enzymeName)
        => SameSiteEnzymes(RequireEnzyme(enzymeName));

    /// <summary>
    /// Built-in isoschizomers of <paramref name="enzymeName"/> that also cut at the same positions,
    /// so they are interchangeable in a digest — equischizomers (Biopython
    /// <c>equischizomers()</c>: "same site, same position of restriction").
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is not a built-in enzyme.</exception>
    public static IEnumerable<RestrictionEnzyme> FindEquischizomers(string enzymeName)
    {
        var query = RequireEnzyme(enzymeName);
        return SameSiteEnzymes(query).Where(e =>
            e.CutPositionForward == query.CutPositionForward &&
            e.CutPositionReverse == query.CutPositionReverse);
    }

    /// <summary>
    /// Built-in isoschizomers of <paramref name="enzymeName"/> that cut the same site at a
    /// different position, so they leave different ends — neoschizomers (Biopython
    /// <c>neoschizomers()</c>: "same site, different position of restriction").
    /// </summary>
    /// <param name="enzymeName">Enzyme name (case-insensitive).</param>
    /// <exception cref="ArgumentNullException"><paramref name="enzymeName"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentException"><paramref name="enzymeName"/> is not a built-in enzyme.</exception>
    public static IEnumerable<RestrictionEnzyme> FindNeoschizomers(string enzymeName)
    {
        var query = RequireEnzyme(enzymeName);
        return SameSiteEnzymes(query).Where(e =>
            e.CutPositionForward != query.CutPositionForward ||
            e.CutPositionReverse != query.CutPositionReverse);
    }

    /// <summary>
    /// The one same-recognition-sequence query the three isoschizomer accessors share. Comparison
    /// is on the recognition string as stored (both are upper-case in the table), and the queried
    /// enzyme itself is excluded, exactly as the reference's <c>del r[i]</c> does.
    /// </summary>
    private static IEnumerable<RestrictionEnzyme> SameSiteEnzymes(RestrictionEnzyme query)
    {
        return _enzymes.Values.Where(e =>
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
/// <param name="CutPositionForward">Forward-strand cut offset from the 5' end of the site.</param>
/// <param name="CutPositionReverse">Reverse-strand cut offset, in forward-strand coordinates.</param>
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
/// <c>Bio.Restriction.&lt;enzyme&gt;.search()</c> reports the same top-strand cut as
/// <c>CutPosition + 1</c> (1-based).</param>
/// <param name="RecognizedSequence">The matched window as read on the strand named by
/// <paramref name="IsForwardStrand"/> — so for a reverse-strand record it is the reverse
/// complement of the forward-strand window. Always
/// <c>Enzyme.RecognitionSequence.Length</c> characters long.</param>
public sealed record RestrictionSite(
    int Position,
    RestrictionEnzyme Enzyme,
    bool IsForwardStrand,
    int CutPosition,
    string RecognizedSequence);

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
