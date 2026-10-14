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
        // Format: Name, Recognition sequence, Cut positions (forward, reverse from 5')

        // 6-cutter enzymes
        ["EcoRI"] = new("EcoRI", "GAATTC", 1, 5, "Escherichia coli"),
        ["BamHI"] = new("BamHI", "GGATCC", 1, 5, "Bacillus amyloliquefaciens"),
        ["HindIII"] = new("HindIII", "AAGCTT", 1, 5, "Haemophilus influenzae"),
        ["XhoI"] = new("XhoI", "CTCGAG", 1, 5, "Xanthomonas holcicola"),
        ["SalI"] = new("SalI", "GTCGAC", 1, 5, "Streptomyces albus"),
        ["PstI"] = new("PstI", "CTGCAG", 5, 1, "Providencia stuartii"),
        ["SphI"] = new("SphI", "GCATGC", 5, 1, "Streptomyces phaeochromogenes"),
        ["KpnI"] = new("KpnI", "GGTACC", 5, 1, "Klebsiella pneumoniae"),
        ["SacI"] = new("SacI", "GAGCTC", 5, 1, "Streptomyces achromogenes"),
        ["XbaI"] = new("XbaI", "TCTAGA", 1, 5, "Xanthomonas badrii"),
        ["NcoI"] = new("NcoI", "CCATGG", 1, 5, "Nocardia corallina"),
        ["NdeI"] = new("NdeI", "CATATG", 2, 4, "Neisseria denitrificans"),
        ["NheI"] = new("NheI", "GCTAGC", 1, 5, "Neisseria mucosa"),
        ["SpeI"] = new("SpeI", "ACTAGT", 1, 5, "Sphaerotilus natans"),
        ["AvrII"] = new("AvrII", "CCTAGG", 1, 5, "Anabaena variabilis"),
        ["ClaI"] = new("ClaI", "ATCGAT", 2, 4, "Caryophanon latum"),
        ["AgeI"] = new("AgeI", "ACCGGT", 1, 5, "Agrobacterium gelatinovorum"),

        // 8-cutter enzymes (rare cutters)
        ["NotI"] = new("NotI", "GCGGCCGC", 2, 6, "Nocardia otitidis-caviarum"),
        ["SfiI"] = new("SfiI", "GGCCNNNNNGGCC", 8, 5, "Streptomyces fimbriatus"), // Note: N = any base
        ["PacI"] = new("PacI", "TTAATTAA", 5, 3, "Pseudomonas alcaligenes"),
        ["AscI"] = new("AscI", "GGCGCGCC", 2, 6, "Arthrobacter species"),
        ["FseI"] = new("FseI", "GGCCGGCC", 6, 2, "Frankia species"),
        ["SwaI"] = new("SwaI", "ATTTAAAT", 4, 4, "Staphylococcus warneri"),

        // 4-cutter enzymes (frequent cutters)
        ["MspI"] = new("MspI", "CCGG", 1, 3, "Moraxella species"),
        ["HpaII"] = new("HpaII", "CCGG", 1, 3, "Haemophilus parainfluenzae"),
        ["TaqI"] = new("TaqI", "TCGA", 1, 3, "Thermus aquaticus"),
        ["AluI"] = new("AluI", "AGCT", 2, 2, "Arthrobacter luteus"), // Blunt
        ["RsaI"] = new("RsaI", "GTAC", 2, 2, "Rhodopseudomonas sphaeroides"), // Blunt
        ["HaeIII"] = new("HaeIII", "GGCC", 2, 2, "Haemophilus aegyptius"), // Blunt
        ["DpnI"] = new("DpnI", "GATC", 2, 2, "Diplococcus pneumoniae"), // Blunt, methylation-dependent
        ["Sau3AI"] = new("Sau3AI", "GATC", 0, 4, "Staphylococcus aureus"),
        ["MboI"] = new("MboI", "GATC", 0, 4, "Moraxella bovis"),

        // Common cloning enzymes
        ["EcoRV"] = new("EcoRV", "GATATC", 3, 3, "Escherichia coli"), // Blunt
        ["SmaI"] = new("SmaI", "CCCGGG", 3, 3, "Serratia marcescens"), // Blunt
        ["HincII"] = new("HincII", "GTYRAC", 3, 3, "Haemophilus influenzae"), // Degenerate, blunt
        ["ScaI"] = new("ScaI", "AGTACT", 3, 3, "Streptomyces caespitosus"), // Blunt
        ["StuI"] = new("StuI", "AGGCCT", 3, 3, "Streptomyces tubercidicus"), // Blunt
        ["BglII"] = new("BglII", "AGATCT", 1, 5, "Bacillus globigii"),
        ["ApaI"] = new("ApaI", "GGGCCC", 5, 1, "Acetobacter pasteurianus"),
    };

    /// <summary>
    /// Gets all available restriction enzymes.
    /// </summary>
    public static IReadOnlyDictionary<string, RestrictionEnzyme> Enzymes => _enzymes;

    /// <summary>
    /// Gets a restriction enzyme by name.
    /// </summary>
    public static RestrictionEnzyme? GetEnzyme(string name) =>
        _enzymes.TryGetValue(name, out var enzyme) ? enzyme : null;

    /// <summary>
    /// Gets all enzymes whose recognition sequence has exactly the given length (in base pairs).
    /// </summary>
    /// <param name="length">Recognition-sequence length to match (Type II sites are 4–8 nt;
    /// see Wikipedia "Restriction enzyme").</param>
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
    public static IEnumerable<RestrictionEnzyme> GetBluntCutters() =>
        _enzymes.Values.Where(e => e.IsBluntEnd);

    /// <summary>
    /// Gets all enzymes that produce sticky (cohesive) ends — a staggered cut leaving a
    /// 5' or 3' single-stranded overhang (Wikipedia, "Sticky and blunt ends").
    /// </summary>
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
    /// codes simply never match. A null or empty sequence yields no sites — but an unknown
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
    /// codes simply never match. A null or empty sequence yields no sites — but an unknown
    /// enzyme name is still rejected, eagerly, before the sequence is looked at. See
    /// <see cref="FindSites(DnaSequence, string, MoleculeTopology)"/> for the coordinate and
    /// feasibility conventions.
    /// </remarks>
    public static IEnumerable<RestrictionSite> FindSites(
        string sequence, string enzymeName, MoleculeTopology topology)
    {
        var enzyme = RequireEnzyme(enzymeName);

        return string.IsNullOrEmpty(sequence)
            ? []
            : FindSitesCore(sequence.ToUpperInvariant(), enzyme, topology);
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
    /// Matching is per position via the canonical <see cref="IupacHelper.MatchesIupac"/>:
    /// ambiguity is interpreted in the <b>pattern</b> only. A sequence position holding anything
    /// other than A, C, G or T — an ambiguity code such as <c>N</c> or <c>R</c> included — matches
    /// no pattern code at all, so it blocks a match rather than widening it.
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

        return sequence.Length == 0
            ? []
            : EnumeratePatternOccurrences(sequence.ToUpperInvariant(), pattern, startIndex);
    }

    private static IEnumerable<int> EnumeratePatternOccurrences(string seq, string pattern, int startIndex)
    {
        for (int i = startIndex; i + pattern.Length <= seq.Length; i++)
        {
            if (MatchesPatternAt(seq, i, pattern))
                yield return i;
        }
    }

    private static RestrictionEnzyme RequireEnzyme(string enzymeName)
    {
        if (string.IsNullOrEmpty(enzymeName))
            throw new ArgumentNullException(nameof(enzymeName));

        return GetEnzyme(enzymeName)
            ?? throw new ArgumentException($"Unknown enzyme: {enzymeName}", nameof(enzymeName));
    }

    /// <summary>
    /// Validates a recognition pattern and returns it upper-cased, ready for matching.
    /// </summary>
    private static string ValidatedPattern(string? pattern, string paramName, string? enzymeName = null)
    {
        string what = enzymeName is null ? "Recognition pattern" : $"Enzyme '{enzymeName}' recognition sequence";

        if (string.IsNullOrEmpty(pattern))
            throw new ArgumentException($"{what} cannot be null or empty.", paramName);

        string upper = pattern.ToUpperInvariant();
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
            if (!IupacHelper.MatchesIupac(seq[position + i], pattern[i]))
                return false;
        }
        return true;
    }

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
    /// together with one representative site per position.
    /// </summary>
    /// <remarks>
    /// Only <see cref="RestrictionSite.IsForwardStrand"/> records are kept: their
    /// <see cref="RestrictionSite.CutPosition"/> is the top-strand cut, which is what Biopython's
    /// <c>catalyse</c> slices at. The reverse-strand record of the same palindromic site carries
    /// the bottom-strand cut and would double-count the site as a second boundary. When several
    /// enzymes cut at the same coordinate only one representative is retained (the last seen), so
    /// the fragment-boundary enzyme name and overhang are that enzyme's.
    /// </remarks>
    private static (List<int> Cuts, Dictionary<int, RestrictionSite> Sites) CollectForwardCuts(
        DnaSequence sequence, string[] enzymeNames, MoleculeTopology topology)
    {
        var cutPositions = new SortedSet<int>();
        var sitesByPosition = new Dictionary<int, RestrictionSite>();

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
                sitesByPosition[site.CutPosition] = site;
            }
        }

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

            var leftSite = i == 0 ? null : sitesByPosition.GetValueOrDefault(start);
            var rightSite = i == boundaries.Count - 2 ? null : sitesByPosition.GetValueOrDefault(end);
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
                RightOverhang: rightEnd.Sequence);
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

            var leftSite = sitesByPosition.GetValueOrDefault(start);
            var rightSite = sitesByPosition.GetValueOrDefault(nextCut);
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
                RightOverhang: rightEnd.Sequence);
        }
    }

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
    private static (OverhangType Type, string Window) CutOverhangWindow(
        string seq, RestrictionSite site, MoleculeTopology topology)
    {
        var enzyme = site.Enzyme;
        var type = enzyme.OverhangType;
        if (type == OverhangType.Blunt)
            return (type, string.Empty);

        int watson = site.CutPosition;
        int crick = watson + (enzyme.CutPositionReverse - enzyme.CutPositionForward);
        int lo = Math.Min(watson, crick);
        int length = Math.Abs(crick - watson);
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
                ? enzymeNames.Select(RequireEnzyme)
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

    /// <summary>
    /// Finds enzymes that produce compatible (ligatable) ends.
    /// </summary>
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
    /// Checks if two enzymes produce compatible sticky ends.
    /// </summary>
    public static bool AreCompatible(string enzyme1Name, string enzyme2Name)
    {
        var e1 = GetEnzyme(enzyme1Name);
        var e2 = GetEnzyme(enzyme2Name);

        if (e1 == null || e2 == null) return false;

        return AreCompatible(e1, e2, out _);
    }

    private static bool AreCompatible(RestrictionEnzyme e1, RestrictionEnzyme e2, out string? compatibleEnd)
    {
        compatibleEnd = null;

        // Blunt ends are always compatible with other blunt ends
        // Source: Wikipedia (Sticky and blunt ends) — "blunt ends are always compatible with each other"
        if (e1.IsBluntEnd && e2.IsBluntEnd)
        {
            compatibleEnd = "blunt";
            return true;
        }

        // For sticky ends, both overhang type AND sequence must match.
        // Source: Wikipedia (Sticky and blunt ends) — "overhangs have to be complementary
        // in order for the ligase to work". A 5' overhang cannot pair with a 3' overhang,
        // even if the overhang sequence is palindromic (same string), because both
        // single-stranded extensions would be on the same strand and cannot base-pair.
        if (e1.OverhangType != e2.OverhangType)
            return false;

        string overhang1 = GetOverhang(e1);
        string overhang2 = GetOverhang(e2);

        if (overhang1 == overhang2 && !string.IsNullOrEmpty(overhang1))
        {
            compatibleEnd = overhang1;
            return true;
        }

        return false;
    }

    private static string GetOverhang(RestrictionEnzyme enzyme)
    {
        if (enzyme.IsBluntEnd) return "";

        int cutDiff = enzyme.CutPositionForward - enzyme.CutPositionReverse;

        // This is a simplified calculation - in reality would need actual sequence
        return enzyme.OverhangType switch
        {
            OverhangType.FivePrime => enzyme.RecognitionSequence.Substring(
                enzyme.CutPositionForward,
                Math.Abs(cutDiff)),
            OverhangType.ThreePrime => enzyme.RecognitionSequence.Substring(
                enzyme.CutPositionReverse,
                Math.Abs(cutDiff)),
            _ => ""
        };
    }

    #endregion
}

/// <summary>
/// Represents a restriction enzyme with its properties.
/// </summary>
public sealed record RestrictionEnzyme(
    string Name,
    string RecognitionSequence,
    int CutPositionForward,
    int CutPositionReverse,
    string Organism)
{
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
/// the fragment is an uncut circle. When several enzymes cut at the same coordinate only one
/// representative name is reported.</param>
/// <param name="RightEnzyme">Name of the enzyme whose cut created the fragment's right (3'-side)
/// end, under the same convention as <paramref name="LeftEnzyme"/>.</param>
/// <param name="FragmentNumber">1-based sequential number in the returned order (ascending
/// <paramref name="StartPosition"/>), contiguous over the digest.</param>
/// <param name="LeftOverhangType">End chemistry of the left end — <see cref="OverhangType.Blunt"/>
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
    string? RightOverhang = null);

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
