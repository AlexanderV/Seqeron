using System.Globalization;
using Seqeron.Genomics.Infrastructure;
using Seqeron.Genomics.Phylogenetics;

namespace Seqeron.Genomics.Annotation;

// TargetScan 7.0 conservation pipeline from a multi-species 3'UTR alignment:
//   targetscan_70.pl            — seed-match sites per species and their cross-species groups
//   targetscan_70_BL_bins.pl    — the UTR's conservation bin (median per-position branch length)
//   targetscan_70_BL_PCT.pl     — branch length (bin tree) and PCT per site
// Sources opened 2026-11-25: www.targetscan.org/vert_70/vert_70_data_download/targetscan_70.zip
// (targetscan_70.pl 7.0, sample input and targetscan_70_output.txt; SHA-256 6e23c72a…d749) and
// targetscan_70_BL_PCT.zip (SHA-256 5ecd0e90…5885). Friedman et al. (2009) Genome Res 19:92.
public static partial class MiRnaAnalyzer
{
    /// <summary>
    /// One seed-match site of a miRNA family in one species of a 3'UTR alignment, as a row of
    /// TargetScan 7.0's <c>targetscan_70.pl</c> output.
    /// </summary>
    /// <param name="Species">The species identifier (key of the alignment dictionary).</param>
    /// <param name="MsaStart">1-based first alignment column of the site (counting gaps; TargetScan <c>MSA_start</c>).</param>
    /// <param name="MsaEnd">1-based last alignment column of the site (TargetScan <c>MSA_end</c>).</param>
    /// <param name="UtrStart">1-based position of the site's first nucleotide in the ungapped UTR (<c>UTR_start</c>).</param>
    /// <param name="UtrEnd"><c>UtrStart</c> + site length − 1 (<c>UTR_end</c>; 8 nt for an 8mer, 7 for a 7mer, 6 for a 6mer).</param>
    /// <param name="SiteType">8mer, 7mer-m8, 7mer-A1 or 6mer.</param>
    /// <param name="OrfOverlap">TargetScan's <c>ORF_overlap</c> flag (lower-case UTR letters mark ORF overlap; see <see cref="FindTargetScanSiteGroups"/>). TargetScan gives such a site no PCT.</param>
    public readonly record struct TargetScanAlignedSite(
        string Species,
        int MsaStart,
        int MsaEnd,
        int UtrStart,
        int UtrEnd,
        TargetSiteType SiteType,
        bool OrfOverlap);

    /// <summary>
    /// A group of seed-match sites of one miRNA family that overlap across species in a 3'UTR alignment
    /// (TargetScan 7.0 <c>targetscan_70.pl</c>: <c>groupSitesThisGeneThisMiRNA</c> /
    /// <c>summarizePrintGroupsThisGeneThisMiRNA</c>).
    /// </summary>
    public sealed class TargetScanSiteGroup
    {
        private readonly IReadOnlyDictionary<TargetSiteType, IReadOnlyList<string>> _speciesBySiteType;

        internal TargetScanSiteGroup(
            int groupNumber,
            IReadOnlyList<TargetScanAlignedSite> sites,
            IReadOnlyList<TargetSiteType> siteTypes,
            IReadOnlyList<string> species,
            IReadOnlyDictionary<TargetSiteType, IReadOnlyList<string>> speciesBySiteType)
        {
            GroupNumber = groupNumber;
            Sites = sites;
            SiteTypes = siteTypes;
            Species = species;
            _speciesBySiteType = speciesBySiteType;
            GroupType = string.Join("+", siteTypes.Select(TargetScanSiteTypeName));
        }

        /// <summary>Group number, 1-based in the order TargetScan assigns them for this UTR and miRNA (the script continues one counter over all genes and families).</summary>
        public int GroupNumber { get; }

        /// <summary>The sites of the group, in TargetScan's output order (ordinal order of "species::MSA_start::MSA_end").</summary>
        public IReadOnlyList<TargetScanAlignedSite> Sites { get; }

        /// <summary>The distinct site types in the group, in TargetScan's type-id order (7mer-A1, 7mer-m8, 8mer, 6mer).</summary>
        public IReadOnlyList<TargetSiteType> SiteTypes { get; }

        /// <summary>TargetScan's <c>Group_type</c> column: the site-type names joined with "+" (e.g. "7mer-1a+8mer-1a+6mer").</summary>
        public string GroupType { get; }

        /// <summary>TargetScan's <c>Species_in_this_group</c>: the distinct species of the group, sorted as the script sorts them (numeric value, then ordinal).</summary>
        public IReadOnlyList<string> Species { get; }

        /// <summary>
        /// The species of the group that have a site of <paramref name="siteType"/> or of a site type
        /// containing it (an 8mer contains a 7mer-m8, a 7mer-A1 and a 6mer; a 7mer contains a 6mer) —
        /// TargetScan's <c>%GROUP_NUM_PLUS_TYPE_2_SPECIES_LIST</c>, printed as
        /// <c>Species_in_this_group_with_this_site_type</c> when the group is heterogeneous and equal to
        /// <see cref="Species"/> when it is not. This is the species list <c>targetscan_70_BL_PCT.pl</c>
        /// computes the branch length of for a site of that type (pass it as
        /// <see cref="PctConservation.SpeciesWithSite"/>). Empty when no site of the group has the type.
        /// </summary>
        /// <param name="siteType">The site type.</param>
        public IReadOnlyList<string> GetSpeciesWithSiteType(TargetSiteType siteType) =>
            _speciesBySiteType.TryGetValue(siteType, out var list) ? list : Array.Empty<string>();
    }

    /// <summary>
    /// Finds the seed-match sites of a miRNA family in every species of a multi-species 3'UTR
    /// alignment and groups the sites that overlap across species, as TargetScan 7.0's
    /// <c>targetscan_70.pl</c> does. The result gives, per group, each species' site (alignment and
    /// ungapped UTR coordinates, site type, ORF-overlap flag), the group type and the per-site-type
    /// species lists that TargetScan's branch-length / PCT step uses.
    /// </summary>
    /// <remarks>
    /// <para>Port of <c>targetscan_70.pl</c> 7.0 (Whitehead Institute; sites per Lewis et al. 2005,
    /// grouping as in TargetScan 7.0). Per species:</para>
    /// <list type="number">
    /// <item>Input normalisation: T/t → U (<c>s/T/U/gi</c>; a lower-case t becomes an upper-case U), then the
    /// first 14 A/C/G/T/U letters (any case) become N (<c>$BEG_UTR_MASK_LENGTH</c>: "miRNAs don't target
    /// right next to CDS").</item>
    /// <item>Seed matches (<c>get_seeds</c>): 7mer-m8 = reverse complement of the seed region (miRNA nt 2–8),
    /// 6mer = its last 6 nt, 7mer-A1 = 6mer + A, 8mer = 7mer-m8 + A, each matched on the upper-cased aligned
    /// sequence with any number of '-' gaps between consecutive nucleotides (<c>makeSeedMatchRegex</c>), at every
    /// start position (overlapping matches are found).</item>
    /// <item>Shorter matches inside longer ones are dropped (<c>findRemoveMatchSubsets</c>): an 8mer drops the
    /// 7mer-m8, 7mer-A1 and 6mer it contains, a 7mer the 6mer it contains.</item>
    /// <item>Coordinates: MSA start/end are the 1-based first/last alignment columns; UTR start = number of
    /// non-gap characters up to and including the first column (the 14 masked letters count), UTR end =
    /// start + site length − 1 (<c>getUTRcoords</c>).</item>
    /// </list>
    /// <para>Grouping (<c>groupSitesThisGeneThisMiRNA</c>): sites of different species are compared all
    /// against all in ordinal order of the key "species::MSA_start::MSA_end"; a pair overlaps when it shares
    /// the start or the end column, or overlaps by ≥ 2 alignment columns (<c>$REQUIRED_OVERLAP = 2</c>), or
    /// one lies inside the other. An overlapping pair joins the group of whichever site already has one, or
    /// opens a new group; two sites that already belong to different groups are NOT merged (the script's
    /// <c>groupThisPair</c>; reproduced because it defines TargetScan's groups). Remaining sites become
    /// single-site groups in key order. Same-species sites are never paired directly.</para>
    /// <para>Species lists: <see cref="TargetScanSiteGroup.Species"/> and
    /// <see cref="TargetScanSiteGroup.GetSpeciesWithSiteType"/> are made non-redundant and sorted as
    /// <c>makeListNonRedundant</c> does (ordinal, then stable by perl numeric value).</para>
    /// <para>ORF overlap (TargetScan 7.0: "Look for sites in ORF-overlapping regions (lowercase UTR) but flag
    /// as such"): <c>getMatches</c> compares <c>substr($alignment, $start, $len)</c> with its upper-cased copy
    /// using the 1-based <c>$start</c> as a 0-based offset, i.e. the window one column to the right of the
    /// site. That window defines TargetScan's <c>ORF_overlap</c> column and its PCT "NA", so it is reproduced.</para>
    /// <para>Not copied (reference quirks that do not define output for valid input): the script reads
    /// species "0" as missing; with <c>$FIND_SITES_ALL_SPECIES = 0</c> (TargetScan 7) only the species in
    /// which the family is annotated are searched — pass them as <paramref name="speciesWithMiRna"/>
    /// (<c>null</c> searches every species, the script's <c>$FIND_SITES_ALL_SPECIES = 1</c>). For the seed
    /// UUUUUUU the 7mer-m8 and 7mer-A1 matches occupy the same columns and the perl keeps whichever its
    /// hash order visits last; this method keeps the 7mer-m8.</para>
    /// </remarks>
    /// <param name="miRnaSeedOrMature">The 7-nt seed region (miRNA nt 2–8) or the mature miRNA (≥ 8 nt; nt 2–8 are used). RNA or DNA alphabet, any case; only A/C/G/U/T.</param>
    /// <param name="alignedUtrs">Species identifier → aligned 3'UTR (gaps '-'; upper case, lower case = ORF overlap). Rows may differ in length.</param>
    /// <param name="speciesWithMiRna">Species in which the miRNA family is annotated (only these are searched, as TargetScan 7); <c>null</c> = all species of the alignment.</param>
    /// <returns>The groups in TargetScan order, numbered from 1.</returns>
    /// <exception cref="ArgumentNullException">An argument, a species key or an aligned UTR is null.</exception>
    /// <exception cref="ArgumentException">The seed is not 7+ nt of A/C/G/U/T, or a species key is empty.</exception>
    public static IReadOnlyList<TargetScanSiteGroup> FindTargetScanSiteGroups(
        string miRnaSeedOrMature,
        IReadOnlyDictionary<string, string> alignedUtrs,
        IReadOnlyCollection<string>? speciesWithMiRna = null)
    {
        string seed = TargetScanSeedRegion(miRnaSeedOrMature);
        ArgumentNullException.ThrowIfNull(alignedUtrs);
        HashSet<string>? annotated = speciesWithMiRna is null ? null : new HashSet<string>(speciesWithMiRna, StringComparer.Ordinal);

        // get_seeds: types 1 (7mer-1a), 2 (7mer-m8), 3 (8mer-1a), 6 (6mer) are searched in TargetScan 7.
        string rseed2 = GetReverseComplement(seed);
        string rseed6 = rseed2[1..];
        var patterns = new (TargetSiteType Type, string Match)[]
        {
            (TargetSiteType.Seed7merA1, rseed6 + "A"),
            (TargetSiteType.Seed7merM8, rseed2),
            (TargetSiteType.Seed8mer, rseed2 + "A"),
            (TargetSiteType.Seed6mer, rseed6),
        };

        var found = new List<FoundSite>();
        foreach (var species in alignedUtrs.Keys.OrderBy(k => k, StringComparer.Ordinal))
        {
            if (species is null) throw new ArgumentNullException(nameof(alignedUtrs), "A species key is null.");
            if (species.Length == 0) throw new ArgumentException("A species key is empty.", nameof(alignedUtrs));
            string utr = alignedUtrs[species] ?? throw new ArgumentNullException(nameof(alignedUtrs), $"The aligned UTR of species '{species}' is null.");
            if (annotated is not null && !annotated.Contains(species))
                continue;
            found.AddRange(FindTargetScanSitesInSpecies(species, NormaliseTargetScanUtr(utr), patterns));
        }

        return GroupTargetScanSites(found);
    }

    // Per-species site with the 1-based alignment columns of its nucleotides (for subset removal).
    private sealed record FoundSite(string Species, int MsaStart, int MsaEnd, TargetSiteType Type, int[] Columns, int UtrStart, bool OrfOverlap)
    {
        public string Key { get; } = $"{Species}::{MsaStart}::{MsaEnd}";
    }

    private static string TargetScanSeedRegion(string miRnaSeedOrMature)
    {
        ArgumentNullException.ThrowIfNull(miRnaSeedOrMature);
        if (miRnaSeedOrMature.Length < 7)
            throw new ArgumentException("Expected the 7-nt seed region or a mature miRNA of ≥ 8 nt.", nameof(miRnaSeedOrMature));
        string seed = SequenceCase.ToUpperRna(miRnaSeedOrMature.Length == 7 ? miRnaSeedOrMature : miRnaSeedOrMature.Substring(1, 7));
        foreach (char c in seed)
            if (c is not ('A' or 'C' or 'G' or 'U'))
                throw new ArgumentException($"Seed region '{seed}' contains a character other than A/C/G/U/T.", nameof(miRnaSeedOrMature));
        return seed;
    }

    // readUTRs: s/T/U/gi, then 14 × s/[ACGTU]/N/i (the first 14 nucleotide letters, any case).
    private static string NormaliseTargetScanUtr(string utr)
    {
        var chars = utr.ToCharArray();
        int masked = 0;
        for (int i = 0; i < chars.Length; i++)
        {
            if (chars[i] is 'T' or 't')
                chars[i] = 'U';
            if (masked < 14 && chars[i] is 'A' or 'a' or 'C' or 'c' or 'G' or 'g' or 'U' or 'u')
            {
                chars[i] = 'N';
                masked++;
            }
        }
        return new string(chars);
    }

    private static List<FoundSite> FindTargetScanSitesInSpecies(
        string species, string utr, (TargetSiteType Type, string Match)[] patterns)
    {
        string upper = utr.ToUpperInvariant();
        var sites = new Dictionary<(int Start, int End), FoundSite>();
        foreach (var (type, match) in patterns)
        {
            for (int i = 0; i < upper.Length; i++)
            {
                int[]? cols = MatchWithGaps(upper, i, match);
                if (cols is null)
                    continue;
                int start = i + 1, end = cols[^1];
                // getMatches: substr($alignment, $start, $len) with the 1-based $start (one column right).
                bool orf = false;
                for (int k = start; k < Math.Min(utr.Length, start + (end - start + 1)); k++)
                    if (char.IsAsciiLetterLower(utr[k])) { orf = true; break; }
                int utrStart = 0;
                for (int k = 0; k < start; k++)
                    if (utr[k] != '-') utrStart++;
                var site = new FoundSite(species, start, end, type, cols, utrStart, orf);
                // Same columns for two types only for the seed UUUUUUU (7mer-m8 = 7mer-1a = AAAAAAA): keep 7mer-m8.
                if (!sites.TryGetValue((start, end), out var prior) || type == TargetSiteType.Seed7merM8 || prior.Type != TargetSiteType.Seed7merM8)
                    sites[(start, end)] = site;
            }
        }

        // findRemoveMatchSubsets (order-independent: every dropped subset of a dropped site is also a
        // subset of the site that dropped it).
        var drop = new HashSet<(int, int)>();
        foreach (var s in sites.Values)
        {
            int[] c = s.Columns;
            switch (s.Type)
            {
                case TargetSiteType.Seed7merA1: drop.Add((c[0], c[^2])); break;                 // 6mer, same start
                case TargetSiteType.Seed7merM8: drop.Add((c[1], c[^1])); break;                 // 6mer, same end
                case TargetSiteType.Seed8mer:
                    drop.Add((c[0], c[^2]));                                                    // 7mer-m8
                    drop.Add((c[1], c[^1]));                                                    // 7mer-1a
                    drop.Add((c[1], c[^2]));                                                    // 6mer
                    break;
            }
        }
        return sites.Where(kv => !drop.Contains(kv.Key)).Select(kv => kv.Value).ToList();
    }

    // The regex  m1 -{0,} m2 -{0,} … mk  anchored at index i: the 1-based columns of the k letters, or null.
    private static int[]? MatchWithGaps(string upper, int i, string match)
    {
        if (upper[i] != match[0])
            return null;
        var cols = new int[match.Length];
        cols[0] = i + 1;
        int j = i;
        for (int k = 1; k < match.Length; k++)
        {
            j++;
            while (j < upper.Length && upper[j] == '-') j++;
            if (j >= upper.Length || upper[j] != match[k])
                return null;
            cols[k] = j + 1;
        }
        return cols;
    }

    private static IReadOnlyList<TargetScanSiteGroup> GroupTargetScanSites(List<FoundSite> found)
    {
        var sorted = found.OrderBy(s => s.Key, StringComparer.Ordinal).ToList();
        var groupOf = new Dictionary<FoundSite, int>(ReferenceEqualityComparer.Instance);
        int groupNum = 0;

        // groupSitesThisGeneThisMiRNA: all vs. all, both orders, different species only.
        foreach (var s1 in sorted)
        {
            foreach (var s2 in sorted)
            {
                if (string.Equals(s1.Species, s2.Species, StringComparison.Ordinal))
                    continue;
                bool pair;
                if (s1.MsaStart == s2.MsaStart || s1.MsaEnd == s2.MsaEnd)
                    pair = true;
                else if (s1.MsaStart > s2.MsaStart && s1.MsaStart <= s2.MsaEnd)
                    pair = s2.MsaEnd - s1.MsaStart + 1 >= TargetScanRequiredOverlap;
                else if (s1.MsaEnd >= s2.MsaStart && s1.MsaEnd < s2.MsaEnd)
                    pair = s1.MsaEnd - s2.MsaStart + 1 >= TargetScanRequiredOverlap;
                else
                    pair = (s1.MsaStart > s2.MsaStart && s1.MsaEnd < s2.MsaEnd) || (s2.MsaStart > s1.MsaStart && s2.MsaEnd < s1.MsaEnd);
                if (!pair)
                    continue;

                // groupThisPair: only a site without a group is assigned; groups are never merged.
                bool has1 = groupOf.TryGetValue(s1, out int g1), has2 = groupOf.TryGetValue(s2, out int g2);
                if (has1 && has2) continue;
                if (has1) groupOf[s2] = g1;
                else if (has2) groupOf[s1] = g2;
                else { groupNum++; groupOf[s1] = groupNum; groupOf[s2] = groupNum; }
            }
        }
        foreach (var s in sorted)
            if (!groupOf.ContainsKey(s))
                groupOf[s] = ++groupNum;

        var groups = new List<TargetScanSiteGroup>(groupNum);
        foreach (var members in sorted.GroupBy(s => groupOf[s]).OrderBy(g => g.Key))
        {
            var list = members.ToList();
            var types = list.Select(s => s.Type).Distinct().OrderBy(TargetScanSiteTypeId).ToList();
            var perType = new Dictionary<TargetSiteType, List<string>>();
            void Add(TargetSiteType t, string sp)
            {
                if (!perType.TryGetValue(t, out var l))
                {
                    l = new List<string>();
                    perType[t] = l;
                }
                l.Add(sp);
            }
            foreach (var s in list)
            {
                Add(s.Type, s.Species);
                // "If a wide site is present, its subset sites are also present"
                switch (s.Type)
                {
                    case TargetSiteType.Seed7merA1 or TargetSiteType.Seed7merM8: Add(TargetSiteType.Seed6mer, s.Species); break;
                    case TargetSiteType.Seed8mer:
                        Add(TargetSiteType.Seed7merA1, s.Species);
                        Add(TargetSiteType.Seed7merM8, s.Species);
                        Add(TargetSiteType.Seed6mer, s.Species);
                        break;
                }
            }
            groups.Add(new TargetScanSiteGroup(
                members.Key,
                list.Select(s => new TargetScanAlignedSite(s.Species, s.MsaStart, s.MsaEnd, s.UtrStart,
                    s.UtrStart + TargetScanSiteLength(s.Type) - 1, s.Type, s.OrfOverlap)).ToList(),
                types,
                PerlNonRedundantSorted(list.Select(s => s.Species)),
                perType.ToDictionary(kv => kv.Key, kv => PerlNonRedundantSorted(kv.Value))));
        }
        return groups;
    }

    // $REQUIRED_OVERLAP in targetscan_70.pl.
    private const int TargetScanRequiredOverlap = 2;

    // getSiteTypeKeys: %SITE_ID_2_SITE_TYPE ids and names, %SITE_ID_2_LENGTH.
    private static int TargetScanSiteTypeId(TargetSiteType t) => t switch
    {
        TargetSiteType.Seed7merA1 => 1,
        TargetSiteType.Seed7merM8 => 2,
        TargetSiteType.Seed8mer => 3,
        TargetSiteType.Seed6mer => 6,
        _ => throw new ArgumentOutOfRangeException(nameof(t), t, "Not a TargetScan site type.")
    };

    private static string TargetScanSiteTypeName(TargetSiteType t) => t switch
    {
        TargetSiteType.Seed7merA1 => "7mer-1a",
        TargetSiteType.Seed7merM8 => "7mer-m8",
        TargetSiteType.Seed8mer => "8mer-1a",
        TargetSiteType.Seed6mer => "6mer",
        _ => throw new ArgumentOutOfRangeException(nameof(t), t, "Not a TargetScan site type.")
    };

    private static int TargetScanSiteLength(TargetSiteType t) => t switch
    {
        TargetSiteType.Seed8mer => 8,
        TargetSiteType.Seed6mer => 6,
        _ => 7
    };

    // makeListNonRedundant: unique, sort (string), then sort { $a <=> $b } — perl's sort is a stable merge
    // sort, so the result is ordered by numeric value, ties in string order.
    private static IReadOnlyList<string> PerlNonRedundantSorted(IEnumerable<string> values) =>
        values.Distinct(StringComparer.Ordinal)
              .OrderBy(PerlNumericValue)
              .ThenBy(v => v, StringComparer.Ordinal)
              .ToList();

    // Perl's numeric value of a string: the leading decimal number (after leading whitespace), else 0.
    private static double PerlNumericValue(string s)
    {
        int i = 0;
        while (i < s.Length && char.IsWhiteSpace(s[i])) i++;
        int start = i;
        if (i < s.Length && s[i] is '+' or '-') i++;
        int digits = 0;
        while (i < s.Length && char.IsAsciiDigit(s[i]))
        {
            i++;
            digits++;
        }
        if (i < s.Length && s[i] == '.')
        {
            i++;
            while (i < s.Length && char.IsAsciiDigit(s[i]))
            {
                i++;
                digits++;
            }
        }
        if (digits == 0) return 0.0;
        if (i < s.Length && s[i] is 'e' or 'E')
        {
            int j = i + 1;
            if (j < s.Length && s[j] is '+' or '-') j++;
            if (j < s.Length && char.IsAsciiDigit(s[j]))
            {
                while (j < s.Length && char.IsAsciiDigit(s[j])) j++;
                i = j;
            }
        }
        return double.Parse(s.AsSpan(start, i - start), NumberStyles.Float, CultureInfo.InvariantCulture);
    }

    // ── targetscan_70_BL_bins.pl: UTR conservation bin ─────────────────────────────────────

    /// <summary>
    /// A 3'UTR's TargetScan conservation bin (<c>targetscan_70_BL_bins.pl</c>).
    /// </summary>
    /// <param name="MedianBranchLength">Median over the reference species' aligned nucleotides of the branch length of the species sharing that nucleotide (<c>NaN</c> when the reference row has no nucleotide; TargetScan then prints an empty value).</param>
    /// <param name="Bin">The bin, 1–10: the largest k whose threshold the median exceeds (<see cref="TargetScanUtrBinThresholds"/>); 1 when it exceeds none.</param>
    public readonly record struct TargetScanUtrBin(double MedianBranchLength, int Bin);

    /// <summary>The 10 lower bounds of TargetScan 7.0's UTR conservation bins (<c>targetscan_70_BL_bins.pl</c> <c>@BL_thresholds</c>, for the bundled generic tree): 0, 1.21207417, …, 5.78196252.</summary>
    public static IReadOnlyList<double> TargetScanUtrBinThresholds => Array.AsReadOnly(TargetScanUtrBinThresholdValues);

    /// <summary>
    /// Assigns a 3'UTR alignment to one of TargetScan 7.0's ten conservation bins, as
    /// <c>targetscan_70_BL_bins.pl</c> does: for every alignment column where the reference species has a
    /// non-gap character, the branch length of the species whose character equals it (case-insensitive;
    /// the reference included) is computed on the generic tree (the total length of the minimal subtree
    /// connecting them, unrounded; one species ⇒ 0); the UTR's value is the median over those columns
    /// (Statistics::Lite <c>median</c> = <see cref="StatisticsHelper.Median"/>) and the bin is the largest
    /// k with median &gt; <see cref="TargetScanUtrBinThresholds"/>[k − 1].
    /// </summary>
    /// <remarks>
    /// The script reads the raw alignment (no T→U, no 5' masking). Species that are not leaves of the tree
    /// are ignored. Columns beyond the end of a shorter row count as "different". Not copied: when a gene
    /// lacks the reference row, the script silently reuses the previous gene's reference row; here the
    /// reference is required. The bin selects the tree for <see cref="ComputeBranchLengthScore"/>
    /// (<see cref="GetTargetScanBinTree"/>).
    /// </remarks>
    /// <param name="alignedUtrs">Species identifier → aligned 3'UTR (same alignment as for <see cref="FindTargetScanSiteGroups"/>).</param>
    /// <param name="referenceSpecies">The reference species (<c>$refGenome</c>; TargetScanHuman: "9606").</param>
    /// <param name="tree">The tree; <c>null</c> ⇒ the bundled TargetScan 7.0 generic tree (<see cref="GetTargetScanGenericTree"/>). The thresholds belong to that tree.</param>
    /// <exception cref="ArgumentNullException">An argument or an aligned UTR is null.</exception>
    /// <exception cref="ArgumentException">The reference species is missing from the alignment or is not a leaf of the tree.</exception>
    public static TargetScanUtrBin AssignTargetScanUtrBin(
        IReadOnlyDictionary<string, string> alignedUtrs,
        string referenceSpecies = TargetScanBinReferenceSpecies,
        PhylogeneticAnalyzer.PhyloNode? tree = null)
    {
        ArgumentNullException.ThrowIfNull(alignedUtrs);
        ArgumentNullException.ThrowIfNull(referenceSpecies);
        if (!alignedUtrs.TryGetValue(referenceSpecies, out string? reference))
            throw new ArgumentException($"The alignment has no row for the reference species '{referenceSpecies}'.", nameof(referenceSpecies));
        if (reference is null)
            throw new ArgumentNullException(nameof(alignedUtrs), "The reference species' aligned UTR is null.");
        tree ??= TargetScanGenericTree.Value;
        if (!PhylogeneticAnalyzer.GetLeaves(tree).Any(leaf => leaf.Name == referenceSpecies))
            throw new ArgumentException($"The reference species '{referenceSpecies}' is not a leaf of the tree.", nameof(referenceSpecies));

        var rows = new List<(string Species, string Seq)>(alignedUtrs.Count);
        foreach (var kv in alignedUtrs)
        {
            if (kv.Value is null)
                throw new ArgumentNullException(nameof(alignedUtrs), $"The aligned UTR of species '{kv.Key}' is null.");
            rows.Add((kv.Key, kv.Value));
        }
        var cache = new Dictionary<string, double>(StringComparer.Ordinal);
        var lengths = new List<double>(reference.Length);
        for (int i = 0; i < reference.Length; i++)
        {
            char r = reference[i];
            if (r == '-')
                continue;
            char rl = char.ToLowerInvariant(r);
            var same = rows.Where(x => i < x.Seq.Length && char.ToLowerInvariant(x.Seq[i]) == rl)
                           .Select(x => x.Species).OrderBy(s => s, StringComparer.Ordinal).ToList();
            string key = string.Join(" ", same);
            if (!cache.TryGetValue(key, out double bl))
                cache[key] = bl = same.Count < 2 ? 0.0 : ConnectingSubtreeLength(tree, same);
            lengths.Add(bl);
        }

        if (lengths.Count == 0)
            return new TargetScanUtrBin(double.NaN, 1);
        double median = StatisticsHelper.Median(lengths);
        int bin = 1;
        for (int k = 0; k < TargetScanUtrBinThresholdValues.Length; k++)
            if (median > TargetScanUtrBinThresholdValues[k])
                bin = k + 1;
        return new TargetScanUtrBin(median, bin);
    }

    // ── Bundled trees ──────────────────────────────────────────────────────────────────────

    private static readonly Lazy<PhylogeneticAnalyzer.PhyloNode> TargetScanGenericTree =
        new(() => PhylogeneticAnalyzer.ParseNewick(TargetScanTreeGenericNewick));

    // Parsed on first use (the lambda reads the generated array at run time, so the order of static
    // initialisation across the partial files does not matter).
    private static readonly Lazy<PhylogeneticAnalyzer.PhyloNode[]> TargetScanBinTrees =
        new(() => Array.ConvertAll(TargetScanBinTreeNewicks, PhylogeneticAnalyzer.ParseNewick));

    /// <summary>
    /// The bundled TargetScan 7.0 tree for a 3'UTR conservation bin (<c>PCT_parameters/Tree.bin_01..10.txt</c>,
    /// 84 species as NCBI taxonomy ids) as single-line Newick text — the tree <c>targetscan_70_BL_PCT.pl</c>
    /// computes a site's branch length on.
    /// </summary>
    /// <param name="bin">The bin, 1–10 (<see cref="AssignTargetScanUtrBin"/>).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="bin"/> is not 1–10.</exception>
    public static string GetTargetScanBinTreeNewick(int bin)
    {
        if (bin is < 1 or > 10)
            throw new ArgumentOutOfRangeException(nameof(bin), bin, "TargetScan bins are 1–10.");
        return TargetScanBinTreeNewicks[bin - 1];
    }

    /// <summary>The bundled TargetScan 7.0 generic tree (<c>PCT_parameters/Tree.generic.txt</c>, used for the UTR bins) as single-line Newick text.</summary>
    public static string TargetScanGenericTreeNewick => TargetScanTreeGenericNewick;

    /// <summary>A new parsed copy of the bundled TargetScan 7.0 tree for <paramref name="bin"/> (see <see cref="GetTargetScanBinTreeNewick"/>).</summary>
    /// <param name="bin">The bin, 1–10.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="bin"/> is not 1–10.</exception>
    public static PhylogeneticAnalyzer.PhyloNode GetTargetScanBinTree(int bin) =>
        PhylogeneticAnalyzer.ParseNewick(GetTargetScanBinTreeNewick(bin));

    /// <summary>A new parsed copy of the bundled TargetScan 7.0 generic tree (see <see cref="TargetScanGenericTreeNewick"/>).</summary>
    public static PhylogeneticAnalyzer.PhyloNode GetTargetScanGenericTree() =>
        PhylogeneticAnalyzer.ParseNewick(TargetScanTreeGenericNewick);

    // ── targetscan_70_BL_PCT.pl: branch length, PCT and conservation per site ──────────────

    /// <summary>
    /// One site of <see cref="ScoreTargetScanConservation"/>: the aligned site with TargetScan's branch
    /// length, PCT and conservation call (a row of <c>targetscan_70_output.BL_PCT.txt</c>).
    /// </summary>
    /// <param name="Site">The aligned site.</param>
    /// <param name="GroupNumber">Its group (<see cref="TargetScanSiteGroup.GroupNumber"/>).</param>
    /// <param name="SpeciesWithSite">The species whose branch length is computed (<see cref="TargetScanSiteGroup.GetSpeciesWithSiteType"/> for the site's type).</param>
    /// <param name="BranchLength">Branch-length score on the bin tree, 4 decimals (<see cref="ComputeBranchLengthScore"/>).</param>
    /// <param name="Pct">PCT, or <c>null</c> for TargetScan's "NA" (seed not one of the 111 PCT families, or ORF-overlapping site).</param>
    /// <param name="Conserved">TargetScan's "Conserved" column (<see cref="IsConservedByBranchLength"/>).</param>
    public readonly record struct TargetScanSiteConservation(
        TargetScanAlignedSite Site,
        int GroupNumber,
        IReadOnlyList<string> SpeciesWithSite,
        double BranchLength,
        double? Pct,
        bool Conserved);

    /// <summary>Result of <see cref="ScoreTargetScanConservation"/>.</summary>
    /// <param name="UtrBin">The UTR's bin (computed, or the caller's bin with a <c>NaN</c> median).</param>
    /// <param name="Groups">The site groups (<see cref="FindTargetScanSiteGroups"/>).</param>
    /// <param name="Sites">Every site with branch length, PCT and conservation, in group order.</param>
    public sealed record TargetScanConservationResult(
        TargetScanUtrBin UtrBin,
        IReadOnlyList<TargetScanSiteGroup> Groups,
        IReadOnlyList<TargetScanSiteConservation> Sites);

    /// <summary>
    /// End-to-end TargetScan 7.0 conservation of a miRNA family's sites in a multi-species 3'UTR
    /// alignment: site groups (<see cref="FindTargetScanSiteGroups"/>, <c>targetscan_70.pl</c>) → UTR bin
    /// (<see cref="AssignTargetScanUtrBin"/>, <c>targetscan_70_BL_bins.pl</c>) → per-site branch length on
    /// the bin's tree, PCT with the bundled sigmoid parameters and the conservation call
    /// (<c>targetscan_70_BL_PCT.pl</c>).
    /// </summary>
    /// <remarks>
    /// Per site (<c>processPreviousGroup</c>): the branch length is that of the group's species with the
    /// site's type (<see cref="TargetScanSiteGroup.GetSpeciesWithSiteType"/>; one species ⇒ 0) on
    /// <c>Tree.bin_NN</c> (<see cref="ComputeBranchLengthScore"/>); PCT is "NA" when the seed is not one of
    /// the 111 PCT families (<see cref="IsTargetScanPctSeedFamily"/>) or the site overlaps the ORF, else
    /// <c>getPCT</c>: 0 for branch length 0 or a 6mer, otherwise the bundled sigmoid
    /// (<see cref="PctFromBranchLength(double, PctSigmoidParameters, TargetSiteType)"/>); conserved per
    /// <see cref="IsConservedByBranchLength"/>. In the script the alternative "best of three branch lengths"
    /// branch for 8mers never runs (it tests <c>$siteType ne "8mer"</c> while the type is "8mer-1a"), so
    /// each site uses its own type's species list, as here. Every row of
    /// <c>targetscan_70_output.BL_PCT.txt</c> of the kit sample is reproduced from the raw alignment.
    /// </remarks>
    /// <param name="miRnaSeedOrMature">The 7-nt seed region or the mature miRNA (see <see cref="FindTargetScanSiteGroups"/>).</param>
    /// <param name="alignedUtrs">Species identifier → aligned 3'UTR.</param>
    /// <param name="speciesWithMiRna">Species in which the family is annotated (searched); <c>null</c> = all.</param>
    /// <param name="referenceSpecies">Reference species for the bin (default "9606").</param>
    /// <param name="utrBin">A bin to use instead of computing it (TargetScan reads bins from a file; 1–10).</param>
    /// <param name="thresholds">Conservation thresholds; <c>null</c> ⇒ <see cref="BranchLengthConservationThresholds.TargetScanHuman7"/>.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="utrBin"/> is not 1–10.</exception>
    public static TargetScanConservationResult ScoreTargetScanConservation(
        string miRnaSeedOrMature,
        IReadOnlyDictionary<string, string> alignedUtrs,
        IReadOnlyCollection<string>? speciesWithMiRna = null,
        string referenceSpecies = TargetScanBinReferenceSpecies,
        int? utrBin = null,
        BranchLengthConservationThresholds? thresholds = null)
    {
        if (utrBin is < 1 or > 10)
            throw new ArgumentOutOfRangeException(nameof(utrBin), utrBin, "TargetScan bins are 1–10.");
        string seed = TargetScanSeedRegion(miRnaSeedOrMature);
        var groups = FindTargetScanSiteGroups(seed, alignedUtrs, speciesWithMiRna);
        var bin = utrBin is int b ? new TargetScanUtrBin(double.NaN, b) : AssignTargetScanUtrBin(alignedUtrs, referenceSpecies);
        var tree = TargetScanBinTrees.Value[bin.Bin - 1];
        bool family = IsTargetScanPctSeedFamily(seed);

        var sites = new List<TargetScanSiteConservation>();
        foreach (var g in groups)
        {
            foreach (var s in g.Sites)
            {
                var species = g.GetSpeciesWithSiteType(s.SiteType);
                double bl = species.Count < 2 ? 0.0 : ComputeBranchLengthScore(tree, species);
                double? pct = null;
                if (family && !s.OrfOverlap)
                {
                    if (s.SiteType == TargetSiteType.Seed6mer || bl == 0.0)
                        pct = 0.0;
                    else
                    {
                        TryGetTargetScanPctParameters(seed, s.SiteType, out var p);
                        pct = PctFromBranchLength(bl, p, s.SiteType);
                    }
                }
                sites.Add(new TargetScanSiteConservation(s, g.GroupNumber, species, bl, pct,
                    IsConservedByBranchLength(bl, s.SiteType, thresholds)));
            }
        }
        return new TargetScanConservationResult(bin, groups, sites);
    }
}
