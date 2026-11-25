using System.Text;
using Seqeron.Genomics.Analysis;
using Seqeron.Genomics.Phylogenetics;

namespace Seqeron.Genomics.Annotation;

/// <summary>
/// Provides algorithms for microRNA analysis, target prediction, and seed matching.
/// </summary>
public static partial class MiRnaAnalyzer
{
    #region Records and Types

    /// <summary>
    /// Represents a microRNA.
    /// </summary>
    public readonly record struct MiRna(
        string Name,
        string Sequence,
        string SeedSequence,
        int SeedStart,
        int SeedEnd);

    /// <summary>
    /// Represents a predicted miRNA target site.
    /// </summary>
    /// <remarks>
    /// <see cref="Score"/> is the predicted repression: <c>max(0, −ContextScore)</c> where
    /// <see cref="ContextScore"/> is the TargetScan 5.0 context score of Grimson et al. (2007)
    /// (see <see cref="MiRnaAnalyzer.ScoreTargetSiteGrimson(string, string, int, int, TargetSiteType)"/>);
    /// 0 for sites TargetScan 5.0 has no coefficients for (6mer, offset 6mer, centered) and for sites too
    /// close to the stop codon. <see cref="ContextScore"/> is <c>null</c> in those cases.
    /// </remarks>
    public readonly record struct TargetSite(
        int Start,
        int End,
        string TargetSequence,
        string MiRnaName,
        TargetSiteType Type,
        int SeedMatchLength,
        double Score,
        double FreeEnergy,
        string Alignment)
    {
        /// <summary>The TargetScan 5.0 / Grimson (2007) context score as printed (log2 units; negative = repressed); <c>null</c> when undefined.</summary>
        public double? ContextScore { get; init; }
    }

    /// <summary>
    /// Types of miRNA target sites.
    /// </summary>
    public enum TargetSiteType
    {
        Seed8mer,      // 8mer: perfect seed + A at position 1
        Seed7merM8,    // 7mer-m8: perfect seed match positions 2-8
        Seed7merA1,    // 7mer-A1: positions 2-7 + A at position 1
        Seed6mer,      // 6mer: positions 2-7
        Offset6mer,    // offset 6mer: positions 3-8
        /// <summary>
        /// Never emitted. 3'-supplementary pairing is not a site class: Grimson et al. (2007) treat it as an
        /// additive term (the 3' pairing score) of a seed site — see
        /// <see cref="MiRnaAnalyzer.ScoreTargetSiteGrimson(string, string, int, int, TargetSiteType)"/>.
        /// Kept for backward compatibility of the public enum.
        /// </summary>
        Supplementary,
        /// <summary>Centered site (Shin et al. 2010): ≥ 11 contiguous Watson–Crick pairs at miRNA nt 4–15; emitted by <see cref="MiRnaAnalyzer.FindCenteredSites"/>, not by <see cref="MiRnaAnalyzer.FindTargetSites"/>.</summary>
        Centered
    }

    /// <summary>
    /// Represents a potential precursor miRNA (pre-miRNA) hairpin.
    /// </summary>
    public readonly record struct PreMiRna(
        int Start,
        int End,
        string Sequence,
        string MatureSequence,
        string StarSequence,
        string Structure,
        double FreeEnergy);

    /// <summary>
    /// miRNA-mRNA duplex with alignment details.
    /// </summary>
    public readonly record struct MiRnaDuplex(
        string MiRnaSequence,
        string TargetSequence,
        string AlignmentString,
        int Matches,
        int Mismatches,
        int GUWobbles,
        int Gaps,
        double FreeEnergy);

    /// <summary>
    /// Minimum-free-energy miRNA–target duplex (ViennaRNA <c>RNAduplex</c> model): intermolecular base
    /// pairs only, with stacks, bulges and internal loops between them.
    /// </summary>
    /// <param name="MiRnaSequence">Normalised miRNA (upper case, T→U).</param>
    /// <param name="TargetSequence">Normalised target (upper case, T→U).</param>
    /// <param name="MiRnaStart">0-based first miRNA position of the duplex (outermost pair).</param>
    /// <param name="MiRnaEnd">0-based last miRNA position of the duplex (innermost pair).</param>
    /// <param name="TargetStart">0-based first target position of the duplex (pairs <paramref name="MiRnaEnd"/>).</param>
    /// <param name="TargetEnd">0-based last target position of the duplex (pairs <paramref name="MiRnaStart"/>).</param>
    /// <param name="Structure">RNAduplex notation: miRNA[MiRnaStart..MiRnaEnd] as '('/'.', then '&amp;', then
    /// target[TargetStart..TargetEnd] as ')'/'.'.</param>
    /// <param name="BasePairs">(miRNA index, target index) pairs, 0-based, outermost first.</param>
    /// <param name="FreeEnergy">Turner 2004 ΔG°37 of the duplex (kcal/mol).</param>
    public readonly record struct MiRnaDuplexStructure(
        string MiRnaSequence,
        string TargetSequence,
        int MiRnaStart,
        int MiRnaEnd,
        int TargetStart,
        int TargetEnd,
        string Structure,
        IReadOnlyList<(int MiRnaIndex, int TargetIndex)> BasePairs,
        double FreeEnergy);

    /// <summary>
    /// Result of comparing two miRNA seed regions.
    /// </summary>
    public readonly record struct SeedComparison(
        int Matches,
        int Mismatches,
        bool IsSameFamily);

    #endregion

    #region Seed Matching

    /// <summary>
    /// Extracts the seed region from a miRNA sequence: nucleotides 2–8 of the mature miRNA
    /// (7 nt, the TargetScan "Seed+m8" family-defining region; Bartel 2009), upper-cased.
    /// Returns "" for null or sequences shorter than 8 nt. The alphabet is preserved (a DNA
    /// input keeps T); <see cref="CreateMiRna"/> performs the T→U normalisation, and
    /// <see cref="CompareSeedRegions"/> / <see cref="GroupBySeedFamily"/> treat T ≡ U.
    /// </summary>
    public static string GetSeedSequence(string miRnaSequence)
    {
        if (string.IsNullOrEmpty(miRnaSequence) || miRnaSequence.Length < 8)
            return "";

        return miRnaSequence.Substring(1, 7).ToUpperInvariant();
    }

    /// <summary>
    /// Creates a MiRna record from a sequence. The sequence is upper-cased and DNA T is
    /// converted to RNA U (the same normalisation TargetScan applies to its seed input:
    /// <c>targetscan_70.pl</c> <c>s/T/U/gi; uc()</c>) before the nt 2–8 seed is extracted.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    public static MiRna CreateMiRna(string name, string sequence)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        string upper = sequence.ToUpperInvariant().Replace('T', 'U');
        string seed = GetSeedSequence(upper);

        return new MiRna(
            Name: name,
            Sequence: upper,
            SeedSequence: seed,
            SeedStart: 1,
            SeedEnd: 7);
    }

    /// <summary>
    /// Compares the seed regions of two miRNAs, returning the number of matches,
    /// mismatches (Hamming distance), and whether they belong to the same seed family.
    /// </summary>
    /// <remarks>
    /// A miRNA family is the set of miRNAs sharing the same sequence at nucleotides 2–8
    /// (Bartel 2009, Cell 136:215; TargetScan <c>miR_Family_Info</c> "Seed+m8"). Seeds are
    /// compared after the TargetScan normalisation (upper-case, T→U; <c>targetscan_70.pl</c>
    /// lines 229–232), so a DNA- or lower-case-encoded seed is the same seed as its RNA form.
    /// Mismatches = Hamming distance over the common prefix (canonical
    /// <see cref="SequenceExtensions.HammingDistance"/>) plus the length difference. An empty
    /// (undefined) seed on either side yields a zeroed, non-family comparison.
    /// </remarks>
    public static SeedComparison CompareSeedRegions(MiRna mirna1, MiRna mirna2)
    {
        string seed1 = NormalizeSeed(mirna1.SeedSequence);
        string seed2 = NormalizeSeed(mirna2.SeedSequence);

        if (seed1.Length == 0 || seed2.Length == 0)
            return new SeedComparison(Matches: 0, Mismatches: 0, IsSameFamily: false);

        int length = Math.Min(seed1.Length, seed2.Length);
        int hamming = seed1.AsSpan(0, length).HammingDistance(seed2.AsSpan(0, length));
        int matches = length - hamming;

        // Account for length differences (if seeds have different lengths)
        int mismatches = hamming + Math.Abs(seed1.Length - seed2.Length);

        bool isSameFamily = mismatches == 0;

        return new SeedComparison(Matches: matches, Mismatches: mismatches, IsSameFamily: isSameFamily);
    }

    /// <summary>
    /// TargetScan seed normalisation (<c>targetscan_70.pl</c>: <c>s/T/U/gi; uc()</c>):
    /// upper-case, DNA T → RNA U. Null → empty.
    /// </summary>
    private static string NormalizeSeed(string? seed) =>
        string.IsNullOrEmpty(seed) ? "" : seed.ToUpperInvariant().Replace('T', 'U');

    /// <summary>
    /// Finds all potential target sites for a miRNA in an mRNA sequence.
    /// Scans for the 6mer core (RC of miRNA positions 2-7), then extends to classify
    /// site types per Bartel (2009) and TargetScan conventions.
    /// </summary>
    /// <remarks>
    /// Each site's duplex (<see cref="TargetSite.TargetSequence"/>, <see cref="TargetSite.Alignment"/>,
    /// <see cref="TargetSite.FreeEnergy"/>) is the full miRNA aligned antiparallel in the seed-match
    /// register, extending UPSTREAM of the site (see <c>CreateTargetSite</c>).
    /// <see cref="TargetSite.Score"/> is the predicted repression <c>max(0, −context score)</c> of the
    /// fitted TargetScan 5.0 model (Grimson et al. 2007: site type + 3' pairing + local AU + position;
    /// <see cref="ScoreTargetSiteGrimson(string, string, int, int, TargetSiteType)"/>); the signed
    /// score is <see cref="TargetSite.ContextScore"/>. 6mer and offset-6mer sites have no TargetScan 5.0
    /// coefficients and score 0, as do sites starting within the first 15 nt of the 3'UTR. For TargetScan
    /// context++ use <see cref="ScoreTargetSiteContextPlusPlus"/>. With the default
    /// <paramref name="minScore"/> 0 every site is returned; the offset-6mer suppression (an offset 6mer
    /// overlapping a 6mer-core site is dropped) does not depend on <paramref name="minScore"/>.
    /// </remarks>
    public static IEnumerable<TargetSite> FindTargetSites(
        string mRnaSequence,
        MiRna miRna,
        double minScore = 0.0)
    {
        if (string.IsNullOrEmpty(mRnaSequence) || string.IsNullOrEmpty(miRna.Sequence))
            yield break;

        string mrna = mRnaSequence.ToUpperInvariant().Replace('T', 'U');
        string mirna = miRna.Sequence;

        string seedRC = GetReverseComplement(miRna.SeedSequence);
        if (seedRC.Length < 7)
            yield break;

        // seedRC layout (7 chars): [RC of pos8] [RC of pos7] ... [RC of pos2]
        // 6mer core = seedRC[1..7] = RC of miRNA positions 2-7 (6 chars)
        // offset 6mer pattern = seedRC[0..6] = RC of miRNA positions 3-8 (6 chars)
        string sixmerCore = seedRC.Substring(1, 6);
        string offset6Pat = seedRC.Substring(0, 6);

        // Track positions covered by 6mer-core-based sites to suppress overlapping offset 6mer
        var coveredPositions = new HashSet<int>();

        // Pass 1: Find 6mer core matches → classify as 8mer/7mer-m8/7mer-A1/6mer
        for (int i = 0; i <= mrna.Length - 6; i++)
        {
            if (mrna.Substring(i, 6) != sixmerCore)
                continue;

            // Check for position 8 match upstream: mrna[i-1] == seedRC[0]
            bool hasPos8 = i > 0 && mrna[i - 1] == seedRC[0];
            // Check for A opposite miRNA position 1 downstream: mrna[i+6] == 'A'
            bool hasA1 = i + 6 < mrna.Length && mrna[i + 6] == 'A';

            TargetSiteType type;
            int siteStart;
            int siteLength;
            int seedMatchLen;

            if (hasPos8 && hasA1)
            {
                // 8mer: match to positions 2-8 + A opposite position 1
                type = TargetSiteType.Seed8mer;
                siteStart = i - 1;
                siteLength = 8;
                seedMatchLen = 8;
            }
            else if (hasPos8)
            {
                // 7mer-m8: match to positions 2-8
                type = TargetSiteType.Seed7merM8;
                siteStart = i - 1;
                siteLength = 7;
                seedMatchLen = 7;
            }
            else if (hasA1)
            {
                // 7mer-A1: match to positions 2-7 + A opposite position 1
                type = TargetSiteType.Seed7merA1;
                siteStart = i;
                siteLength = 7;
                seedMatchLen = 7;
            }
            else
            {
                // 6mer: match to positions 2-7 only
                type = TargetSiteType.Seed6mer;
                siteStart = i;
                siteLength = 6;
                seedMatchLen = 6;
            }

            // The mRNA base opposite miRNA nt 1 is the one just 3' of the 6mer core (i + 6).
            // Suppression of overlapping offset 6mers must not depend on the score threshold:
            // every 6mer-core site claims its positions whether or not it is emitted.
            for (int j = siteStart; j < siteStart + siteLength; j++)
                coveredPositions.Add(j);

            var site = CreateTargetSite(mrna, siteStart, siteLength, i + 6, mirna, miRna.Name, type, seedMatchLen);
            if (site.Score >= minScore)
                yield return site;
        }

        // Pass 2: Find offset 6mer matches (positions 3-8) not overlapping higher-priority sites
        for (int i = 0; i <= mrna.Length - 6; i++)
        {
            if (mrna.Substring(i, 6) != offset6Pat)
                continue;

            // Skip if this position overlaps with a site already found
            bool overlaps = false;
            for (int j = i; j < i + 6; j++)
            {
                if (coveredPositions.Contains(j))
                {
                    overlaps = true;
                    break;
                }
            }
            if (overlaps) continue;

            // Also skip if full seedRC matches here (would have been found in pass 1 at i+1)
            if (i + 7 <= mrna.Length && mrna[i + 6] == seedRC[6])
                continue;

            // Offset 6mer pairs miRNA nt 3-8 at mrna[i..i+5]; nt 1 is therefore opposite i + 7.
            var site = CreateTargetSite(mrna, i, 6, i + 7, mirna, miRna.Name, TargetSiteType.Offset6mer, 6);
            if (site.Score >= minScore)
            {
                yield return site;
            }
        }
    }

    /// <summary>
    /// Builds a <see cref="TargetSite"/>, aligning the full miRNA antiparallel to the mRNA in the
    /// register fixed by the seed match. miRNA nt k pairs with mRNA index
    /// <paramref name="nt1Index"/> − (k − 1), so the miRNA 3' end pairs UPSTREAM (5') of the seed
    /// match (Bartel 2009 Fig. 1; TargetScan <c>extractSubseqForAlignment</c> takes the pairing
    /// subsequence from utrStart − 16 up to the site end). <c>TargetSequence</c> is the mRNA segment
    /// opposite the miRNA: mrna[max(0, nt1Index − L + 1) .. nt1Index] (truncated at the mRNA ends).
    /// When nt1Index lies past the mRNA 3' end (a 6mer / offset 6mer flush with the end), the missing
    /// partner(s) are padded with 'N' for the duplex only, so miRNA nt 1 keeps its register and is
    /// reported as unpaired.
    /// </summary>
    private static TargetSite CreateTargetSite(string mrna, int pos, int length, int nt1Index, string mirna, string mirnaName, TargetSiteType type, int seedMatchLength)
    {
        int windowStart = Math.Max(0, nt1Index - mirna.Length + 1);
        int windowEnd = Math.Min(nt1Index, mrna.Length - 1);
        string targetSeq = mrna.Substring(windowStart, windowEnd - windowStart + 1);
        int missing3Prime = nt1Index - windowEnd;
        string duplexTarget = missing3Prime > 0 ? targetSeq + new string('N', missing3Prime) : targetSeq;

        var duplex = AlignMiRnaToTarget(mirna, duplexTarget);
        double? context = ScoreTargetSiteGrimson(mrna, mirna, pos, pos + length - 1, type)?.ContextScore;
        double score = context is double c ? Math.Max(0.0, -c) : 0.0;

        return new TargetSite(
            Start: pos,
            End: pos + length - 1,
            TargetSequence: targetSeq,
            MiRnaName: mirnaName,
            Type: type,
            SeedMatchLength: seedMatchLength,
            Score: score,
            FreeEnergy: duplex.FreeEnergy,
            Alignment: duplex.AlignmentString)
        {
            ContextScore = context
        };
    }

    #endregion

    #region Sequence Utilities

    /// <summary>
    /// Gets the RNA reverse complement of a sequence (Biopython <c>reverse_complement_rna</c>):
    /// the sequence is reversed and each base is complemented by the canonical
    /// <see cref="SequenceExtensions.GetRnaComplementBase(char)"/> — A↔U, G↔C, DNA T→A, IUPAC ambiguity
    /// codes complemented (R↔Y, K↔M, B↔V, D↔H; S, W, N unchanged), case-insensitive with upper-case
    /// output for recognised codes; any other character (gap, '*', …) passes through unchanged.
    /// Null or empty input yields "".
    /// </summary>
    public static string GetReverseComplement(string rnaSequence)
    {
        if (string.IsNullOrEmpty(rnaSequence))
            return "";

        return string.Create(rnaSequence.Length, rnaSequence, static (dest, src) =>
        {
            for (int i = 0; i < src.Length; i++)
                dest[i] = SequenceExtensions.GetRnaComplementBase(src[src.Length - 1 - i]);
        });
    }

    /// <summary>
    /// Checks if two bases can pair (Watson-Crick A·U / G·C or the G·U wobble).
    /// Inputs are case-insensitive and DNA T is treated as RNA U. Delegates to the canonical
    /// <see cref="RnaSecondaryStructure.CanPair(char, char)"/>.
    /// </summary>
    public static bool CanPair(char base1, char base2) => RnaSecondaryStructure.CanPair(base1, base2);

    /// <summary>
    /// Checks if a pair is a G·U wobble.
    /// Inputs are case-insensitive and DNA T is treated as RNA U. Delegates to the canonical
    /// <see cref="RnaSecondaryStructure.GetBasePairType(char, char)"/>.
    /// </summary>
    public static bool IsWobblePair(char base1, char base2) =>
        RnaSecondaryStructure.GetBasePairType(base1, base2) == RnaSecondaryStructure.BasePairType.Wobble;

    #endregion

    #region Alignment and Scoring

    /// <summary>
    /// Aligns a miRNA to a target sequence as an ungapped antiparallel duplex: miRNA position i
    /// (5'→3') faces target[target.Length − 1 − i], over min(len) positions. Each position is a
    /// Watson-Crick match ('|'), a G·U wobble (':') or a mismatch (' '); inputs are case-insensitive
    /// and DNA T is read as U. <see cref="MiRnaDuplex.FreeEnergy"/> is the Turner 2004
    /// nearest-neighbour ΔG°37 (kcal/mol) of that duplex with every Watson-Crick / G·U position
    /// paired (intermolecular initiation, stacks, internal loops, terminal AU/GU, exterior dangling
    /// ends / terminal mismatches; + the 0.43 symmetry term when both strands are identical), as ViennaRNA
    /// <c>eval_structure</c> evaluates the same <c>miRNA&amp;target</c> structure (Vienna rounds the
    /// parameters to 0.1 kcal/mol and omits the symmetry term); 0 when nothing pairs. Null/empty input
    /// ⇒ empty duplex.
    /// </summary>
    public static MiRnaDuplex AlignMiRnaToTarget(string miRnaSequence, string targetSequence)
    {
        if (string.IsNullOrEmpty(miRnaSequence) || string.IsNullOrEmpty(targetSequence))
        {
            return new MiRnaDuplex("", "", "", 0, 0, 0, 0, 0);
        }

        string mirna = miRnaSequence.ToUpperInvariant().Replace('T', 'U');
        string target = targetSequence.ToUpperInvariant().Replace('T', 'U');

        // miRNA pairs antiparallel to the target: miRNA 5'→3' index i pairs with the target
        // base read 3'→5', i.e. target[target.Length-1-i].
        // Source: Lewis BP et al. (2005) Cell 120:15–20 — targets are the reverse complement
        // of the miRNA seed; Watson-Crick (A-U, C-G) + G:U wobble (Crick 1966).
        var alignmentChars = new char[Math.Min(mirna.Length, target.Length)];
        int matches = 0, mismatches = 0, wobbles = 0;

        for (int i = 0; i < alignmentChars.Length; i++)
        {
            char m = mirna[i];
            int targetIdx = target.Length - 1 - i;
            if (targetIdx < 0) break;

            char t = target[targetIdx];

            if (CanPair(m, t))
            {
                if (IsWobblePair(m, t))
                {
                    alignmentChars[i] = ':';
                    wobbles++;
                }
                else
                {
                    alignmentChars[i] = '|';
                    matches++;
                }
            }
            else
            {
                alignmentChars[i] = ' ';
                mismatches++;
            }
        }

        double freeEnergy = CalculateDuplexEnergy(mirna, target, alignmentChars);

        return new MiRnaDuplex(
            MiRnaSequence: mirna,
            TargetSequence: target,
            AlignmentString: new string(alignmentChars),
            Matches: matches,
            Mismatches: mismatches,
            GUWobbles: wobbles,
            Gaps: 0,
            FreeEnergy: freeEnergy);
    }

    // NNDB Turner 2004 intermolecular initiation, ΔG°37 = 4.09 kcal/mol, applied once for bimolecular
    // duplex formation (rna.urmc.rochester.edu/NNDB/rna_2004/rna_2004_watson_crick_helices.html §1.1;
    // Xia et al. 1998). ViennaRNA 2.x rna_turner2004.par rounds it to 4.10 (DuplexInit37 = 410).
    private const double DuplexIntermolecularInitiation = 4.09;

    // NNDB Turner 2004 symmetry correction, ΔG°37 = +0.43 kcal/mol, for a duplex of two identical
    // (self-complementary) strands — same NNDB page, worked example "Self complementary duplex".
    private const double DuplexSymmetryCorrection = 0.43;

    /// <summary>
    /// Turner 2004 nearest-neighbour free energy (kcal/mol, 37 °C) of the reported ungapped duplex:
    /// every in-register Watson-Crick / G·U position of <paramref name="alignment"/> is a base pair
    /// (miRNA i · target[L−1−i]), evaluated by <see cref="EvaluateDuplexEnergy"/>.
    /// </summary>
    private static double CalculateDuplexEnergy(string mirna, string target, char[] alignment)
    {
        var pairs = new List<(int MiRnaIndex, int TargetIndex)>(alignment.Length);
        for (int i = 0; i < alignment.Length; i++)
            if (alignment[i] == '|' || alignment[i] == ':')
                pairs.Add((i, target.Length - 1 - i));
        return EvaluateDuplexEnergy(mirna, target, pairs);
    }

    /// <summary>
    /// Turner 2004 nearest-neighbour ΔG°37 (kcal/mol) of an intermolecular miRNA·target duplex given by
    /// its base pairs (miRNA index ascending, target index descending; normalised upper-case RNA), as
    /// ViennaRNA <c>eval_structure</c> evaluates the same <c>miRNA&amp;target</c> structure:
    /// ΔG = intermolecular initiation (4.09) [+ symmetry 0.43 for a symmetric duplex of two identical
    /// strands] + Σ helix stacks (incl. the GGUC/CUGG special stack) + Σ bulge / internal loops between
    /// helices + terminal AU/GU penalty on each exterior helix end + exterior-loop stacking at each of the
    /// two duplex ends (terminal mismatch when the end is extended on both strands, otherwise the single
    /// 5′/3′ dangling end) — NNDB Turner 2004 helix, bulge, internal-loop, dangling-end and
    /// terminal-mismatch pages. Every term comes from the canonical <see cref="RnaSecondaryStructure"/>
    /// Turner 2004 kernels. No pair ⇒ no duplex ⇒ 0. Non-ACGU neighbours (e.g. the 'N' padding past an
    /// mRNA end) are treated as absent.
    /// </summary>
    private static double EvaluateDuplexEnergy(string mirna, string target, IReadOnlyList<(int MiRnaIndex, int TargetIndex)> pairs)
    {
        if (pairs.Count == 0)
            return 0.0;

        // Concatenated complex s = miRNA + target; target index t sits at s[lm + t].
        string s = mirna + target;
        int lm = mirna.Length;

        double energy = DuplexIntermolecularInitiation + SymmetryTerm(mirna, target, pairs);

        // Helices = maximal runs of stacked pairs (miRNA +1, target −1), joined by bulge / internal loops.
        int runStart = 0;
        for (int k = 1; k <= pairs.Count; k++)
        {
            if (k < pairs.Count && pairs[k].MiRnaIndex == pairs[k - 1].MiRnaIndex + 1
                && pairs[k].TargetIndex == pairs[k - 1].TargetIndex - 1)
                continue;

            var helix = new List<RnaSecondaryStructure.BasePair>(k - runStart);
            for (int r = runStart; r < k; r++)
                helix.Add(DuplexBasePair(s, lm, pairs[r]));

            // Stacks + GGUC special + terminal AU/GU on both helix ends (canonical). The AU/GU end
            // penalty belongs only to helix ends facing the exterior loop; an end closing a bulge
            // (n > 1) or internal loop is charged by that loop's model, and a 1-nt bulge continues the
            // helix (NNDB), so the penalty is removed for loop-facing ends.
            energy += RnaSecondaryStructure.CalculateStemEnergy(s, helix);
            if (runStart > 0 && Turner2004Parameters.IsTerminalPenaltyPair(helix[0].Base1, helix[0].Base2))
                energy -= Turner2004Parameters.TerminalAuGuPenalty;
            if (k < pairs.Count && Turner2004Parameters.IsTerminalPenaltyPair(helix[^1].Base1, helix[^1].Base2))
                energy -= Turner2004Parameters.TerminalAuGuPenalty;

            if (k < pairs.Count)
                energy += DuplexLoopEnergy(s, lm, pairs[k - 1], pairs[k]);

            runStart = k;
        }

        var (m0, t0) = pairs[0];
        var (ml, tl) = pairs[^1];
        energy += OuterEndStacking(mirna, target, m0, t0) + InnerEndStacking(mirna, target, ml, tl);

        return Math.Round(energy, 2);
    }

    private static RnaSecondaryStructure.BasePair DuplexBasePair(string s, int lm, (int MiRnaIndex, int TargetIndex) p)
    {
        char b1 = s[p.MiRnaIndex], b2 = s[lm + p.TargetIndex];
        return new RnaSecondaryStructure.BasePair(p.MiRnaIndex, lm + p.TargetIndex, b1, b2,
            RnaSecondaryStructure.GetBasePairType(b1, b2)!.Value);
    }

    // NNDB symmetry term: only a duplex of two identical strands whose pair set maps onto itself under
    // strand exchange (miRNA a · target b ⇔ miRNA b · target a) is self-complementary.
    private static double SymmetryTerm(string mirna, string target, IReadOnlyList<(int MiRnaIndex, int TargetIndex)> pairs)
    {
        if (!string.Equals(mirna, target, StringComparison.Ordinal))
            return 0.0;
        var set = new HashSet<(int, int)>(pairs);
        foreach (var (a, b) in pairs)
            if (!set.Contains((b, a)))
                return 0.0;
        return DuplexSymmetryCorrection;
    }

    // Loop closed by the outer pair (miRNA k · target l) and the inner pair (miRNA i · target j),
    // k < i, l > j: stack, bulge or internal loop, scored by the canonical Turner 2004 kernels.
    private static double DuplexLoopEnergy(string s, int lm, (int MiRnaIndex, int TargetIndex) outer, (int MiRnaIndex, int TargetIndex) inner)
    {
        int k = outer.MiRnaIndex, ko = lm + outer.TargetIndex;
        int i = inner.MiRnaIndex, io = lm + inner.TargetIndex;
        int n1 = i - k - 1, n2 = ko - io - 1;

        if (n1 == 0 && n2 == 0)
            return Turner2004Parameters.StackingEnergies.TryGetValue($"{s[k]}{s[i]}/{s[ko]}{s[io]}", out double st) ? st : 0.0;

        if (n1 == 0 || n2 == 0)
        {
            int size = n1 + n2;
            char bulged = n1 > 0 ? s[k + 1] : s[ko - 1];
            // NNDB bulge degeneracy for a 1-nt bulge: identical neighbours on the bulged strand give
            // equivalent states (same counting as the canonical RnaSecondaryStructure MFE recursion).
            int states = 1;
            if (size == 1)
            {
                if (n1 == 1)
                {
                    if (s[k] == s[k + 1]) states++;
                    if (s[i] == s[k + 1]) states++;
                }
                else
                {
                    if (s[ko] == s[ko - 1]) states++;
                    if (s[io] == s[ko - 1]) states++;
                }
            }
            // The bulged nucleotide lies on the miRNA strand (n1 > 0: between outer pair base s[k] and inner
            // pair base s[i], the kernel's 5' strand) or on the target strand (n2 > 0: between s[io] and s[ko]).
            // The strand decides NNDB's special-C term ("a bulged C adjacent to at least one paired C"), so it is
            // passed exactly as the canonical MFE recursion does; without it the kernel applies the bonus only
            // when a paired C flanks the bulge on both strands.
            return RnaSecondaryStructure.CalculateBulgeLoopEnergy(size, bulged, s[k], s[ko], s[i], s[io], states,
                bulgeOn5PrimeStrand: n1 > 0);
        }

        return RnaSecondaryStructure.CalculateInternalLoopEnergy(
            n1, n2, s[k], s[ko], s[i], s[io], s[k + 1], s[ko - 1], s[i - 1], s[io + 1]);
    }

    private static bool IsRnaNucleotide(char c) => c is 'A' or 'C' or 'G' or 'U';

    // Exterior stacking on the duplex end holding the miRNA 5' side: pair (miRNA m · target t), seen
    // from the loop as (target t, miRNA m); 3' neighbour of the target base = target[t+1], 5' neighbour
    // of the miRNA base = miRNA[m−1].
    private static double OuterEndStacking(string mirna, string target, int m, int t) =>
        ExteriorEndStacking(target[t], mirna[m],
            m > 0 && IsRnaNucleotide(mirna[m - 1]) ? mirna[m - 1] : null,
            t + 1 < target.Length && IsRnaNucleotide(target[t + 1]) ? target[t + 1] : null);

    // Exterior stacking on the duplex end holding the strand break: pair (miRNA m · target t);
    // 3' neighbour of the miRNA base = miRNA[m+1], 5' neighbour of the target base = target[t−1].
    private static double InnerEndStacking(string mirna, string target, int m, int t) =>
        ExteriorEndStacking(mirna[m], target[t],
            t > 0 && IsRnaNucleotide(target[t - 1]) ? target[t - 1] : null,
            m + 1 < mirna.Length && IsRnaNucleotide(mirna[m + 1]) ? mirna[m + 1] : null);

    // Exterior-loop stacking on a helix end whose pair, read 5'→3' from the loop side, is (x, y):
    // x is the paired base whose 3' neighbour is <paramref name="threePrime"/>, y the base whose 5'
    // neighbour is <paramref name="fivePrime"/>; plus the terminal AU/GU penalty of the end pair.
    // NNDB Turner 2004 dangling-end page: "if a helix end is extended on both the 5' and 3' strands,
    // then a terminal mismatch exists (not the sum of 5' and 3' dangling ends)"; otherwise the single
    // dangling end applies.
    private static double ExteriorEndStacking(char x, char y, char? fivePrime, char? threePrime)
    {
        // NNDB dangle tables are keyed (top base, dangling base, bottom base) with the top base 5'→3':
        // a 3' dangle sits 3' of the top base x, a 5' dangle sits 5' of the bottom base y.
        if (fivePrime is char f && threePrime is char t)
            return RnaSecondaryStructure.GetTerminalMismatchEnergy(x, y, t, f);
        if (threePrime is char t3)
            return RnaSecondaryStructure.GetDanglingEndEnergy(x, y, t3, is3Prime: true);
        if (fivePrime is char f5)
            return RnaSecondaryStructure.GetDanglingEndEnergy(x, y, f5, is3Prime: false);
        return 0.0;
    }

    /// <summary>
    /// Minimum-free-energy intermolecular duplex of a miRNA and a target (the ViennaRNA
    /// <c>RNAduplex</c> / <c>duplexfold</c> model; Lorenz et al. 2011): only miRNA·target pairs are
    /// allowed (no intramolecular structure), and consecutive pairs are joined by stacks, bulges or
    /// internal loops of at most <paramref name="maxLoopSize"/> unpaired nt in total. Recursion
    /// (Vienna <c>duplexfold</c>): C(i,j) = init + exterior(i,j) or min over the next-outer pair (k,l) of
    /// C(k,l) + loop(k,l;i,j); the duplex energy is min C(i,j) + exterior on the inner end. All loop,
    /// stacking, dangling-end / terminal-mismatch and AU/GU terms are the NNDB Turner 2004 values from the
    /// canonical <see cref="RnaSecondaryStructure"/> kernels (incl. the GGUC/CUGG special stack and the
    /// 1-nt bulge degeneracy, which ViennaRNA omits); the reported energy is that of the returned
    /// structure (see <see cref="AlignMiRnaToTarget"/> for the ungapped fixed-register duplex).
    /// </summary>
    /// <param name="miRnaSequence">miRNA 5'→3' (RNA or DNA; case-insensitive; T read as U).</param>
    /// <param name="targetSequence">Target 5'→3' (RNA or DNA; case-insensitive; T read as U).</param>
    /// <param name="maxLoopSize">Largest bulge / internal loop (total unpaired nt); default 30 = ViennaRNA MAXLOOP.</param>
    /// <returns>The MFE duplex, or null when the inputs are null/empty or no base pair is possible.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoopSize"/> is negative.</exception>
    public static MiRnaDuplexStructure? FoldMiRnaDuplex(string miRnaSequence, string targetSequence, int maxLoopSize = 30)
    {
        ArgumentOutOfRangeException.ThrowIfNegative(maxLoopSize);
        if (string.IsNullOrEmpty(miRnaSequence) || string.IsNullOrEmpty(targetSequence))
            return null;

        string m = miRnaSequence.ToUpperInvariant().Replace('T', 'U');
        string t = targetSequence.ToUpperInvariant().Replace('T', 'U');
        int n1 = m.Length, n2 = t.Length;
        string s = m + t;

        var c = new double[n1, n2];
        var pred = new int[n1, n2];          // k·n2 + l of the next-outer pair; −1 = outermost pair
        var predGguc = new bool[n1, n2];     // predecessor reached through the GGUC/CUGG 3-stack
        double emin = double.PositiveInfinity;
        int iMin = -1, jMin = -1;

        for (int i = 0; i < n1; i++)
        {
            for (int j = n2 - 1; j >= 0; j--)
            {
                c[i, j] = double.PositiveInfinity;
                if (!RnaSecondaryStructure.CanPair(m[i], t[j]))
                    continue;

                double best = DuplexIntermolecularInitiation + OuterEndStacking(m, t, i, j) + AuGuEnd(m[i], t[j]);
                int bestPred = -1;
                bool bestGguc = false;

                for (int k = i - 1; k >= 0 && i - k - 1 <= maxLoopSize; k--)
                {
                    for (int l = j + 1; l < n2; l++)
                    {
                        if ((i - k - 1) + (l - j - 1) > maxLoopSize)
                            break;
                        if (double.IsPositiveInfinity(c[k, l]))
                            continue;
                        double e = c[k, l] + DuplexLoopEnergy(s, n1, (k, l), (i, j));
                        if (e < best)
                        {
                            best = e;
                            bestPred = k * n2 + l;
                            bestGguc = false;
                        }
                    }
                }

                // NNDB GU note: 5'GGUC/3'CUGG replaces its three stacks by −4.12 (canonical stem kernel).
                if (i >= 3 && j + 3 < n2 && !double.IsPositiveInfinity(c[i - 3, j + 3])
                    && string.CompareOrdinal(m, i - 3, "GGUC", 0, 4) == 0
                    && string.CompareOrdinal(t, j, "GGUC", 0, 4) == 0)
                {
                    var four = new List<RnaSecondaryStructure.BasePair>(4);
                    for (int r = 3; r >= 0; r--)
                        four.Add(DuplexBasePair(s, n1, (i - r, j + r)));
                    double e = c[i - 3, j + 3] + RnaSecondaryStructure.CalculateStemEnergy(s, four);
                    if (e < best)
                    {
                        best = e;
                        bestPred = (i - 3) * n2 + (j + 3);
                        bestGguc = true;
                    }
                }

                c[i, j] = best;
                pred[i, j] = bestPred;
                predGguc[i, j] = bestGguc;

                double total = best + InnerEndStacking(m, t, i, j) + AuGuEnd(m[i], t[j]);
                if (total < emin)
                {
                    emin = total;
                    iMin = i;
                    jMin = j;
                }
            }
        }

        if (iMin < 0)
            return null;

        // Backtrack from the innermost pair to the outermost one.
        var pairs = new List<(int MiRnaIndex, int TargetIndex)>();
        int ci = iMin, cj = jMin;
        while (true)
        {
            pairs.Add((ci, cj));
            int p = pred[ci, cj];
            if (p < 0)
                break;
            int pk = p / n2, pl = p % n2;
            if (predGguc[ci, cj])
            {
                pairs.Add((ci - 1, cj + 1));
                pairs.Add((ci - 2, cj + 2));
            }
            ci = pk;
            cj = pl;
        }
        pairs.Reverse();

        int mStart = pairs[0].MiRnaIndex, mEnd = pairs[^1].MiRnaIndex;
        int tStart = pairs[^1].TargetIndex, tEnd = pairs[0].TargetIndex;
        var sb = new StringBuilder(mEnd - mStart + tEnd - tStart + 3);
        var mPaired = new HashSet<int>(pairs.Select(q => q.MiRnaIndex));
        var tPaired = new HashSet<int>(pairs.Select(q => q.TargetIndex));
        for (int x = mStart; x <= mEnd; x++) sb.Append(mPaired.Contains(x) ? '(' : '.');
        sb.Append('&');
        for (int x = tStart; x <= tEnd; x++) sb.Append(tPaired.Contains(x) ? ')' : '.');

        return new MiRnaDuplexStructure(
            MiRnaSequence: m,
            TargetSequence: t,
            MiRnaStart: mStart,
            MiRnaEnd: mEnd,
            TargetStart: tStart,
            TargetEnd: tEnd,
            Structure: sb.ToString(),
            BasePairs: pairs,
            FreeEnergy: EvaluateDuplexEnergy(m, t, pairs));
    }

    private static double AuGuEnd(char b1, char b2) =>
        Turner2004Parameters.IsTerminalPenaltyPair(b1, b2) ? Turner2004Parameters.TerminalAuGuPenalty : 0.0;

    #endregion

    #region TargetScan context++ scoring (Agarwal et al. 2015) — opt-in

    /// <summary>
    /// Breakdown of the locally-computable TargetScan context++ score (Agarwal et al. 2015)
    /// for a single seed-matched target site. The context++ score is the SUM of per-feature
    /// contributions of a multiple-linear-regression model fit SEPARATELY per site type
    /// (8mer / 7mer-m8 / 7mer-A1 / 6mer). Features computable from the miRNA and the 3'UTR, plus
    /// TA_3UTR / SPS (bundled TargetScan 7.0 seed-region table), are always realised; Len_ORF / ORF8m
    /// need the ORF, PCT a conservation input, and an exact SA the RNAplfold value. Features without
    /// data are reported as <see cref="OmittedFeatures"/>; <see cref="ContextScorePartial"/> is then a
    /// PARTIAL context++ score, not the published headline CS. See algorithm doc §5.3.
    /// </summary>
    /// <param name="SiteType">The seed-match site type (must be one of the four canonical types).</param>
    /// <param name="Intercept">Site-type intercept term (raw coefficient; Agarwal_2015_parameters.txt).</param>
    /// <param name="LocalAuContribution">coeff(Local_AU) × scaled(local-AU fraction).</param>
    /// <param name="SRna1Contribution">Sum of sRNA1A/C/G indicator contributions (0 if miRNA nt1 = U).</param>
    /// <param name="SRna8Contribution">Sum of sRNA8A/C/G indicator contributions (0 if miRNA nt8 = U).</param>
    /// <param name="Site8Contribution">Sum of Site8A/C/G contributions (only 7mer-A1 / 6mer; else 0).</param>
    /// <param name="SaContribution">coeff(SA) × scaled(log10 of the unpaired probability of the 14-nt window centred on the match to miRNA nt 7–8). The probability is <see cref="ContextPlusPlusInputs.PlfoldUnpairedProbability"/> when supplied (exactly TargetScan's value), otherwise an estimate from the Turner-2004 McCaskill partition function on one 80-nt window (not RNAplfold's 40-nt-span window average — it can differ from TargetScan by several log10 units). 0 (and listed as omitted) when the window does not fit the 3'UTR or the probability is 0.</param>
    /// <param name="ThreePrimePairingContribution">coeff(3P_score) × scaled(3' supplementary-pairing raw score). Computed from miRNA + 3'UTR.</param>
    /// <param name="MinDistContribution">coeff(Min_dist) × scaled(log10 distance to nearest 3'UTR end). Computed from the 3'UTR.</param>
    /// <param name="Len3UtrContribution">coeff(Len_3UTR) × scaled(log10 3'UTR length). Computed from the 3'UTR.</param>
    /// <param name="Off6mContribution">coeff(Off6m) × offset-6mer count in the 3'UTR (used raw). Computed from miRNA + 3'UTR.</param>
    /// <param name="SpsContribution">coeff(SPS) × scaled(SPS). SPS is <see cref="ContextPlusPlusInputs.Sps"/> when supplied, otherwise the bundled TargetScan 7.0 value for the miRNA's seed region (nt 2–8) and site type (<see cref="TryGetTargetScanTaSps"/>). 0 (and listed as omitted) only when neither exists (seed region not 7 nt over ACGU).</param>
    /// <param name="TaContribution">coeff(TA_3UTR) × scaled(TA). TA is <see cref="ContextPlusPlusInputs.Ta"/> when supplied, otherwise the bundled TargetScan 7.0 value for the seed region (<see cref="TryGetTargetScanTaSps"/>). 0 (and listed as omitted) only when neither exists.</param>
    /// <param name="LenOrfContribution">coeff(Len_ORF) × scaled(log10 ORF length). The length is <see cref="ContextPlusPlusInputs.OrfLength"/>, or the length of <see cref="ContextPlusPlusInputs.OrfSequence"/> without gaps. 0 (and listed as omitted) when neither is supplied.</param>
    /// <param name="Orf8mContribution">coeff(ORF8m) × ORF-8mer count (used raw). The count is <see cref="ContextPlusPlusInputs.Orf8mCount"/>, or the TargetScan count in <see cref="ContextPlusPlusInputs.OrfSequence"/>. 0 (and listed as omitted) when neither is supplied.</param>
    /// <param name="PctContribution">coeff(PCT) × scaled(PCT value). PCT is computed from the Friedman et al. (2009) branch-length score (Bls) via the published sigmoid; 0 (and listed as omitted) unless a <see cref="ContextPlusPlusInputs.Conservation"/> input is supplied.</param>
    /// <param name="BranchLengthScore">The Friedman et al. (2009) branch-length score (Bls, 4 decimals as TargetScan reports it) computed from the supplied conservation tree, or 0 when no conservation input was supplied.</param>
    /// <param name="Pct">The probability of conserved targeting (PCT) computed from <paramref name="BranchLengthScore"/> as TargetScan's <c>getPCT</c> does (sigmoid, 4 decimals, truncated at 0; 0 when Bls = 0 or for a 6mer site), or 0 when no conservation input was supplied.</param>
    /// <param name="ContextScorePartial">Unrounded, uncapped sum of all realised contributions above. See <see cref="ContextScore"/> for the value TargetScan reports.</param>
    /// <param name="OmittedFeatures">context++ features NOT included for this call (honest residual; depends on which optional inputs were supplied).</param>
    public readonly record struct ContextPlusPlusScore(
        TargetSiteType SiteType,
        double Intercept,
        double LocalAuContribution,
        double SRna1Contribution,
        double SRna8Contribution,
        double Site8Contribution,
        double SaContribution,
        double ThreePrimePairingContribution,
        double MinDistContribution,
        double Len3UtrContribution,
        double Off6mContribution,
        double SpsContribution,
        double TaContribution,
        double LenOrfContribution,
        double Orf8mContribution,
        double PctContribution,
        double BranchLengthScore,
        double Pct,
        double ContextScorePartial,
        IReadOnlyList<string> OmittedFeatures)
    {
        /// <summary>
        /// The context++ score as TargetScan 7 reports it (<c>targetscan_70_context_scores.pl</c>,
        /// <c>readTargets</c>): every feature contribution rounded to 3 decimals
        /// (<c>getAgarwalContribution</c>: <c>sprintf("%.3f")</c>), summed with the intercept, the sum
        /// rounded to 3 decimals, then capped at the per-site-type ceiling <c>$maxContextScore</c>
        /// (8mer −0.03, 7mer-m8 −0.02, 7mer-A1 −0.01, 6mer 0). <c>null</c> when the site starts within
        /// the first 14 nt of the 3'UTR (1-based start &lt; <c>$MIN_DIST_TO_CDS</c> = 15): TargetScan
        /// reports such a site as <c>too_close</c> and gives it no score. Over the same contributions as
        /// <see cref="ContextScorePartial"/> (omitted features count 0, as TargetScan's NA inputs do).
        /// </summary>
        public double? ContextScore { get; init; }

        /// <summary>
        /// <c>true</c> when a conservation input was supplied without caller sigmoid parameters and
        /// the miRNA's seed region (nt 2–8) is not one of TargetScan 7's 111 PCT seed families:
        /// <c>targetscan_70_BL_PCT.pl</c> <c>getPCT</c> then reports PCT "NA" and
        /// <c>targetscan_70_context_scores.pl</c> uses PCT contribution 0. <see cref="Pct"/> is 0 and
        /// the feature is not listed as omitted (this is TargetScan's complete score for the site).
        /// </summary>
        public bool PctIsNa { get; init; }
    }

    /// <summary>
    /// Optional context++ feature inputs that the library cannot derive from the miRNA and the
    /// 3'UTR sequence alone. Each contribution is computed with the verbatim Agarwal et al. (2015)
    /// coefficient and scaling. All are <c>null</c> by default. SPS and TA then come from the bundled
    /// TargetScan 7.0 seed-region table (as <c>targetscan_70_context_scores.pl</c> reads them); the
    /// other null inputs leave their feature as an honest residual (reported in
    /// <see cref="ContextPlusPlusScore.OmittedFeatures"/>, contribution 0).
    /// </summary>
    /// <param name="Sps">Predicted seed-pairing stability (kcal/mol; min-max scaled per site type). Overrides the bundled TargetScan 7.0 value (<c>TA_SPS_by_seed_region.txt</c>, Garcia et al. 2011 Supp. Data 5; see <see cref="TryGetTargetScanTaSps"/>). <c>null</c> ⇒ the bundled value for the miRNA's seed region (nt 2–8) and site type.</param>
    /// <param name="Ta">Target-site abundance (log10 site count; min-max scaled). Overrides the bundled TargetScan 7.0 value (see <see cref="TryGetTargetScanTaSps"/>; <see cref="ComputeTa3Utr"/> computes it from your own 3'UTR set). <c>null</c> ⇒ the bundled value for the seed region.</param>
    /// <param name="OrfLength">ORF (CDS) length in nucleotides (log10-transformed, then min-max scaled). Overrides the length derived from <paramref name="OrfSequence"/>.</param>
    /// <param name="Orf8mCount">Count of cognate 8mer sites in the ORF (used raw, no scaling). Overrides the count derived from <paramref name="OrfSequence"/>.</param>
    /// <param name="PlfoldUnpairedProbability">Optional exact RNAplfold accessibility value for the site — the number TargetScan reads in <c>getSA_contribution</c>: row <c>utrStart + 7</c> (7mer-A1: <c>utrStart + 6</c>; 1-based seed-match start), column <c>l=14</c> of the <c>_lunp</c> file written by <c>RNAplfold -L 40 -W 80 -u 20</c> on the whole 3'UTR. When supplied (&gt; 0) it replaces the library's own SA estimate, so SA equals TargetScan's. <c>null</c> ⇒ the library estimates it (one 80-nt McCaskill window, no 40-nt span limit, no window averaging — an approximation of RNAplfold, see the algorithm doc).</param>
    /// <param name="Conservation">Multi-species conservation input (a phylogenetic tree + the set of species in which the site is conserved; the per-site-type PCT sigmoid parameters default to the bundled TargetScan 7.0 tables). When supplied, the library computes the Friedman et al. (2009) branch-length score (Bls) and the resulting PCT value, then applies the bundled Agarwal et al. (2015) PCT coefficient. <c>null</c> ⇒ PCT stays an honest residual.</param>
    /// <param name="OrfSequence">The transcript's ORF (CDS) sequence (RNA or DNA, any case; aligned rows with <c>-</c>/<c>.</c> gaps are accepted). Len_ORF and ORF8m are then derived as TargetScan's <c>targetscan_count_8mers.pl</c> does: gaps removed, length = remaining characters; ORF8m = number of non-overlapping, case-insensitive occurrences of reverse_complement(seed region nt 2–8) + "A" (the 8mer site) on the given strand (see <see cref="CountTargetScanOrf8merSites"/>). <see cref="OrfLength"/> / <see cref="Orf8mCount"/>, when supplied, take precedence. <c>null</c> ⇒ both features stay residual unless supplied directly.</param>
    public readonly record struct ContextPlusPlusInputs(
        double? Sps = null,
        double? Ta = null,
        double? OrfLength = null,
        int? Orf8mCount = null,
        PctConservation? Conservation = null,
        double? PlfoldUnpairedProbability = null,
        string? OrfSequence = null);

    // ── PCT (probability of conserved targeting) — Friedman et al. (2009) Genome Res 19:92 ──
    // Bls definition (Methods, retrieved verbatim this session): "The conservation of a given
    // sequence (e.g., an 8mer miR-1 site in a particular 3'UTR) was then assessed by summing the
    // total branch length in the phylogenetic tree connecting the subset of species having the
    // sequence perfectly aligned." (doi:10.1101/gr.082701.108; PMC2612969)
    // PCT definition (Methods): "the P_CT was defined as E[(S − B)/S]" — a signal-to-background
    // correction (targetscan.org/docs/pct.html: PCT ≈ (S/B − 1)/(S/B), or near zero for S/B < 1).
    // TargetScan parameterises the Bls→PCT relationship as a logistic (sigmoid) curve
    // (targetscan_70_BL_PCT.pl, calculatePCTthisBL, retrieved verbatim this session):
    //   $pct = $b0 + ( $b1 / (1 + $eConstant ** ( (0 - $b2) * $BL + $b3))); if ($pct < 0) { $pct = 0; }
    // i.e. PCT(Bls) = b0 + b1 / (1 + e^(−b2·Bls + b3)), printed with sprintf("%.4f"), truncated at 0;
    // getPCT returns "0.0" for Bls = 0 and for 6mer sites, and getBranchLength returns Bls with
    // sprintf("%.4f"). The per-miRNA-family b0..b3 are distributed in TargetScan's BL/PCT kit
    // (targetscan_70_BL_PCT.zip → PCT_parameters/{8mer,7mer_m8,7mer_1a}_PCT_parameters.txt, 111
    // seed families each; reuse with citation, TargetScan FAQ). They are bundled verbatim in
    // MiRnaAnalyzer.TargetScanPctParameters.cs (generated from the kit files, SHA-256 in its header)
    // and looked up by TryGetTargetScanPctParameters; a caller-supplied row still overrides.

    /// <summary>
    /// The four logistic (sigmoid) parameters that map a branch-length score (Bls) to a PCT value
    /// for one miRNA seed family + site type, per <c>targetscan_70_BL_PCT.pl</c>
    /// (<c>calculatePCTthisBL</c>): <c>PCT(Bls) = B0 + B1 / (1 + e^(−B2·Bls + B3))</c>, rounded to
    /// 4 decimals and truncated at 0.
    /// </summary>
    /// <remarks>
    /// TargetScan 7 distributes the fitted values in its BL/PCT kit
    /// (<c>targetscan_70_BL_PCT.zip</c>, <c>PCT_parameters/8mer_PCT_parameters.txt</c>,
    /// <c>7mer_m8_PCT_parameters.txt</c>, <c>7mer_1a_PCT_parameters.txt</c>; one row per seed
    /// region, columns 2–5 = b0..b3). There is no 6mer table (6mer PCT is always 0). The library
    /// bundles the three tables (111 seed families each, values verbatim; Friedman et al. 2009 /
    /// Agarwal et al. 2015, cite TargetScan when used): see
    /// <see cref="TryGetTargetScanPctParameters"/>. A caller may still supply its own row.
    /// </remarks>
    /// <param name="B0">Offset term (b0).</param>
    /// <param name="B1">Scale term (b1).</param>
    /// <param name="B2">Bls slope inside the logistic (b2).</param>
    /// <param name="B3">Logistic offset (b3).</param>
    public readonly record struct PctSigmoidParameters(double B0, double B1, double B2, double B3);

    /// <summary>
    /// Multi-species conservation evidence for one target site: the phylogenetic tree relating the
    /// aligned species, the set of species in which the seed match is conserved (perfectly aligned),
    /// and optionally the per-site-type PCT sigmoid parameters. The library derives the Friedman
    /// et al. (2009) branch-length score (Bls) from <see cref="Tree"/> + <see cref="SpeciesWithSite"/>
    /// and the PCT value from <see cref="SigmoidParameters"/>, or — when that is <c>null</c> — from
    /// the bundled TargetScan 7.0 row for the miRNA's seed (nt 2–8) and the site type
    /// (<see cref="TryGetTargetScanPctParameters"/>), as <c>targetscan_70_BL_PCT.pl</c> does.
    /// </summary>
    /// <param name="Tree">A phylogenetic tree (<see cref="PhylogeneticAnalyzer.PhyloNode"/>; e.g. parsed from Newick) whose LEAF names are species identifiers and whose <c>BranchLength</c> values are the edge lengths.</param>
    /// <param name="SpeciesWithSite">The species (leaf names) in which the same site type is aligned at this site, the reference species included (TargetScan: <c>Species_in_this_group_with_this_site_type</c>, or the group's species when all have one type). Bls = total branch length of the minimal subtree connecting these species (Friedman 2009), see <see cref="ComputeBranchLengthScore"/>.</param>
    /// <param name="SigmoidParameters">The b0..b3 row for this seed family and this site type (see <see cref="PctSigmoidParameters"/>); overrides the bundled table. <c>null</c> (default) ⇒ the bundled TargetScan 7.0 row for the miRNA's seed region; a seed that is not one of TargetScan's 111 PCT families has no PCT ("NA", contribution 0, as TargetScan). Ignored for a 6mer site (PCT 0).</param>
    public readonly record struct PctConservation(
        PhylogeneticAnalyzer.PhyloNode Tree,
        IReadOnlyCollection<string> SpeciesWithSite,
        PctSigmoidParameters? SigmoidParameters = null);

    // ── Agarwal et al. (2015) context++ fitted coefficients ──────────────────────────────
    // Source (retrieved verbatim this session): TargetScan distribution parameter file
    //   Agarwal_2015_parameters.txt
    //   https://raw.githubusercontent.com/nsoranzo/targetscan/main/Agarwal_2015_parameters.txt
    // Computation/scaling logic (retrieved verbatim this session): targetscan_70_context_scores.pl
    //   getAgarwalContribution(), getLocalAU_contribution(), get_sRNA1_8_contributions(),
    //   getSite8_contribution()
    //   https://raw.githubusercontent.com/nsoranzo/targetscan/main/targetscan_70_context_scores.pl
    // Peer-reviewed model: Agarwal V, Bell GW, Nam JW, Bartel DP (2015) eLife 4:e05005,
    //   doi:10.7554/eLife.05005 — "one [regression] model for each of the four site types".
    // Parameter-file column order is 8mer / 7mer-m8 / 7mer-A1 / 6mer.

    // Intercept row (raw coeff, no scaling): siteType2siteTypeContribution in the perl.
    private const double CtxIntercept8mer = -0.589;
    private const double CtxIntercept7merM8 = -0.224;
    private const double CtxIntercept7merA1 = -0.195;
    private const double CtxIntercept6mer = -0.079;

    // Local_AU row: coeff, then min/max for min-max scaling of the local-AU fraction.
    private const double CtxLocalAuCoeff8mer = -0.254, CtxLocalAuMin8mer = 0.308, CtxLocalAuMax8mer = 0.814;
    private const double CtxLocalAuCoeff7merM8 = -0.177, CtxLocalAuMin7merM8 = 0.277, CtxLocalAuMax7merM8 = 0.782;
    private const double CtxLocalAuCoeff7merA1 = -0.075, CtxLocalAuMin7merA1 = 0.342, CtxLocalAuMax7merA1 = 0.801;
    private const double CtxLocalAuCoeff6mer = -0.040, CtxLocalAuMin6mer = 0.295, CtxLocalAuMax6mer = 0.772;

    // sRNA position-1 identity indicators (binary, used raw): coeff per site type.
    private const double CtxSRna1A8mer = -0.018, CtxSRna1A7merM8 = 0.010, CtxSRna1A7merA1 = -0.025, CtxSRna1A6mer = -0.002;
    private const double CtxSRna1C8mer = -0.021, CtxSRna1C7merM8 = 0.014, CtxSRna1C7merA1 = -0.021, CtxSRna1C6mer = 0.004;
    private const double CtxSRna1G8mer = 0.060, CtxSRna1G7merM8 = 0.062, CtxSRna1G7merA1 = 0.030, CtxSRna1G6mer = 0.018;

    // sRNA position-8 identity indicators (binary, used raw): coeff per site type.
    private const double CtxSRna8A8mer = 0.022, CtxSRna8A7merM8 = 0.004, CtxSRna8A7merA1 = -0.049, CtxSRna8A6mer = -0.015;
    private const double CtxSRna8C8mer = 0.012, CtxSRna8C7merM8 = -0.031, CtxSRna8C7merA1 = 0.033, CtxSRna8C6mer = 0.016;
    private const double CtxSRna8G8mer = 0.015, CtxSRna8G7merM8 = -0.008, CtxSRna8G7merA1 = -0.017, CtxSRna8G6mer = 0.006;

    // Site8 identity indicators (binary, used raw): defined ONLY for 7mer-A1 and 6mer
    // (the perl computes them only for siteType 1 = 7mer-A1 and 4 = 6mer).
    private const double CtxSite8A7merA1 = 0.000, CtxSite8A6mer = -0.002;
    private const double CtxSite8C7merA1 = 0.036, CtxSite8C6mer = 0.015;
    private const double CtxSite8G7merA1 = 0.015, CtxSite8G6mer = 0.012;

    // 3P_score row (3' supplementary pairing): coeff, then min/max for min-max scaling of the raw
    // pairing score. Agarwal_2015_parameters.txt: min=1, max=3.5 for all four site types.
    private const double Ctx3PCoeff8mer = -0.040, Ctx3PCoeff7merM8 = -0.055, Ctx3PCoeff7merA1 = -0.060, Ctx3PCoeff6mer = -0.024;
    private const double Ctx3PMin = 1.0, Ctx3PMax = 3.5;

    // Min_dist row: coeff, then min/max for min-max scaling of log10(distance to nearest 3'UTR end).
    private const double CtxMinDistCoeff8mer = 0.118, CtxMinDistCoeff7merM8 = 0.056, CtxMinDistCoeff7merA1 = 0.045, CtxMinDistCoeff6mer = 0.036;
    private const double CtxMinDistMin8mer = 1.415, CtxMinDistMax8mer = 3.113;
    private const double CtxMinDistMin7merM8 = 1.491, CtxMinDistMax7merM8 = 3.096;
    private const double CtxMinDistMin7merA1 = 1.431, CtxMinDistMax7merA1 = 3.117;
    private const double CtxMinDistMin6mer = 1.477, CtxMinDistMax6mer = 3.106;

    // Len_3UTR row: coeff, then min/max for min-max scaling of log10(3'UTR length).
    private const double CtxLen3UtrCoeff8mer = 0.310, CtxLen3UtrCoeff7merM8 = 0.154, CtxLen3UtrCoeff7merA1 = 0.129, CtxLen3UtrCoeff6mer = 0.045;
    private const double CtxLen3UtrMin8mer = 2.392, CtxLen3UtrMax8mer = 3.637;
    private const double CtxLen3UtrMin7merM8 = 2.409, CtxLen3UtrMax7merM8 = 3.615;
    private const double CtxLen3UtrMin7merA1 = 2.413, CtxLen3UtrMax7merA1 = 3.630;
    private const double CtxLen3UtrMin6mer = 2.405, CtxLen3UtrMax6mer = 3.620;

    // Off6m row (offset-6mer count): coeff only — used raw (NOT in the min-max regex of getAgarwalContribution).
    private const double CtxOff6mCoeff8mer = -0.020, CtxOff6mCoeff7merM8 = -0.011, CtxOff6mCoeff7merA1 = -0.020, CtxOff6mCoeff6mer = -0.010;

    // SPS row (seed-pairing stability; bundled TargetScan table value or caller value): coeff, then min/max for min-max scaling.
    private const double CtxSpsCoeff8mer = 0.210, CtxSpsCoeff7merM8 = 0.135, CtxSpsCoeff7merA1 = 0.095, CtxSpsCoeff6mer = 0.035;
    private const double CtxSpsMin8mer = -11.13, CtxSpsMax8mer = -5.52;
    private const double CtxSpsMin7merM8 = -11.13, CtxSpsMax7merM8 = -5.49;
    private const double CtxSpsMin7merA1 = -8.41, CtxSpsMax7merA1 = -3.33;
    private const double CtxSpsMin6mer = -8.57, CtxSpsMax6mer = -3.33;

    // TA_3UTR row (target-site abundance; bundled TargetScan table value or caller value): coeff, then min/max for min-max scaling.
    private const double CtxTaCoeff8mer = 0.222, CtxTaCoeff7merM8 = 0.139, CtxTaCoeff7merA1 = 0.117, CtxTaCoeff6mer = 0.058;
    private const double CtxTaMin8mer = 3.113, CtxTaMax8mer = 3.865;
    private const double CtxTaMin7merM8 = 3.067, CtxTaMax7merM8 = 3.887;
    private const double CtxTaMin7merA1 = 3.145, CtxTaMax7merA1 = 3.887;
    private const double CtxTaMin6mer = 3.113, CtxTaMax6mer = 3.887;

    // Len_ORF row (ORF length; from OrfSequence or caller value): coeff, then min/max for min-max scaling of log10(ORF length).
    private const double CtxLenOrfCoeff8mer = 0.205, CtxLenOrfCoeff7merM8 = 0.100, CtxLenOrfCoeff7merA1 = 0.063, CtxLenOrfCoeff6mer = 0.029;
    private const double CtxLenOrfMin8mer = 2.788, CtxLenOrfMax8mer = 3.753;
    private const double CtxLenOrfMin7merM8 = 2.773, CtxLenOrfMax7merM8 = 3.729;
    private const double CtxLenOrfMin7merA1 = 2.773, CtxLenOrfMax7merA1 = 3.730;
    private const double CtxLenOrfMin6mer = 2.775, CtxLenOrfMax6mer = 3.731;

    // ORF8m row (ORF 8mer count; from OrfSequence or caller value): coeff only — used raw (NOT min-max scaled).
    private const double CtxOrf8mCoeff8mer = -0.118, CtxOrf8mCoeff7merM8 = -0.044, CtxOrf8mCoeff7merA1 = -0.058, CtxOrf8mCoeff6mer = -0.060;

    // PCT row (probability of conserved targeting; Friedman et al. 2009): coeff, then min/max for
    // min-max scaling of the PCT value. PCT enters getAgarwalContribution exactly like the other
    // scaled features (targetscan_70_context_scores.pl getPCT_contribution → getAgarwalContribution,
    // PCT matches the min-max regex branch). Verbatim from Agarwal_2015_parameters.txt (PCT row),
    // column order 8mer / 7mer-m8 / 7mer-A1 / 6mer (coeff×4, min×4, max×4):
    //   PCT  -0.103  -0.048  -0.048  0.005  0  0  0  0  0.816  0.364  0.449  0.193
    private const double CtxPctCoeff8mer = -0.103, CtxPctMin8mer = 0.0, CtxPctMax8mer = 0.816;
    private const double CtxPctCoeff7merM8 = -0.048, CtxPctMin7merM8 = 0.0, CtxPctMax7merM8 = 0.364;
    private const double CtxPctCoeff7merA1 = -0.048, CtxPctMin7merA1 = 0.0, CtxPctMax7merA1 = 0.449;
    private const double CtxPctCoeff6mer = 0.005, CtxPctMin6mer = 0.0, CtxPctMax6mer = 0.193;

    // SA row (structural accessibility): coeff, then min/max for min-max scaling of log10(plfold).
    // Verbatim from Agarwal_2015_parameters.txt (SA row), column order 8mer / 7mer-m8 / 7mer-A1 / 6mer:
    //   SA  -0.115  -0.134  -0.077  -0.028  -4.356  -5.218  -4.23  -5.082  -0.661  -0.725  -0.588  -0.666
    private const double CtxSaCoeff8mer = -0.115, CtxSaMin8mer = -4.356, CtxSaMax8mer = -0.661;
    private const double CtxSaCoeff7merM8 = -0.134, CtxSaMin7merM8 = -5.218, CtxSaMax7merM8 = -0.725;
    private const double CtxSaCoeff7merA1 = -0.077, CtxSaMin7merA1 = -4.23, CtxSaMax7merA1 = -0.588;
    private const double CtxSaCoeff6mer = -0.028, CtxSaMin6mer = -5.082, CtxSaMax6mer = -0.666;

    // SA structural-accessibility window parameters (targetscan_70_context_scores.pl).
    //   getSA_contribution reads column index 13 of the RNAplfold _lunp row 7 nt downstream of the
    //   seed-match start → the 14th unpaired-probability column → the probability that the 14-nt
    //   stretch ENDING at that row is unpaired (RNAplfold man page; Bernhart et al. 2006). The
    //   eLife paper (Agarwal 2015, Fig 4A): "the log10 value of the unpaired probability for a
    //   14-nt region centered on the match to miRNA nucleotides 7 and 8".
    private const int SaUnpairedWindowLength = 14;   // RNAplfold -u column 14 used by getSA
    private const int SaRowOffsetDownstream = 7;      // grep -A7 → the row 7 nt 3' of utrStart
    // RNAplfold local-fold parameters used by runRNAplfold_all_UTRs: "RNAplfold -L 40 -W 80 -u 20".
    private const int SaPlfoldWindowSize = 80;        // -W (sliding-window averaging size)

    // Number of flanking nucleotides used for the local-AU feature (getLocalAU_contribution
    // extracts 30 nt up- and downstream of the site).
    private const int LocalAuFlankLength = 30;

    // 3' supplementary-pairing alignment constants (targetscan_70_context_scores.pl):
    //   $DESIRED_UTR_ALIGNMENT_LENGTH = 23; site-type-specific UTR start offset is 16 nt
    //   upstream of the seed match.
    private const int DesiredUtrAlignmentLength = 23;
    private const int ThreePrimeUtrStartOffset = 16;

    /// <summary>
    /// Computes the TargetScan 7 context++ score (Agarwal et al. 2015) for one seed-matched target
    /// site, as an OPT-IN alternative to the heuristic <c>Score</c> on <see cref="TargetSite"/>.
    /// Port of <c>targetscan_70_context_scores.pl</c> with the coefficients of
    /// <c>Agarwal_2015_parameters.txt</c>, for a single 3'UTR isoform (AIR = 1).
    /// <para>
    /// Computed from the miRNA and the 3'UTR: site-type intercept, Local_AU, sRNA1 / sRNA8 identity,
    /// Site8 identity (7mer-A1 / 6mer), 3' supplementary pairing (3P_score), Min_dist, Len_3UTR and
    /// Off6m. Structural accessibility (SA) uses a caller-supplied RNAplfold value
    /// (<see cref="ContextPlusPlusInputs.PlfoldUnpairedProbability"/>) when given; otherwise it is
    /// estimated with the Turner-2004 McCaskill partition function on one 80-nt window
    /// (<see cref="RnaSecondaryStructure.CalculateRegionUnpairedProbability"/>), which approximates
    /// RNAplfold (no 40-nt span limit, no averaging over windows).
    /// </para>
    /// <para>
    /// SPS and TA_3UTR come from the bundled TargetScan 7.0 table for the miRNA's seed region
    /// (<see cref="TryGetTargetScanTaSps"/>, as <c>targetscan_70_context_scores.pl</c> reads
    /// <c>TA_SPS_by_seed_region.txt</c>) unless the caller supplies them (<paramref name="inputs"/>;
    /// <see cref="ComputeTa3Utr"/> computes TA from a 3'UTR set). Len_ORF and ORF8m come from
    /// <see cref="ContextPlusPlusInputs.OrfSequence"/> (as <c>targetscan_count_8mers.pl</c>) or from
    /// caller numbers. PCT needs a tree + the species with the site (Friedman et al. 2009
    /// branch-length score, then the TargetScan sigmoid with the bundled TargetScan 7.0
    /// per-seed-family row, see <see cref="TryGetTargetScanPctParameters"/>, or a caller-supplied
    /// row). A feature without data contributes 0 and is listed in
    /// <see cref="ContextPlusPlusScore.OmittedFeatures"/>; the sum is then a partial score.
    /// </para>
    /// <para>
    /// <see cref="ContextPlusPlusScore.ContextScorePartial"/> is the unrounded sum;
    /// <see cref="ContextPlusPlusScore.ContextScore"/> applies TargetScan's reporting rules
    /// (3-decimal rounding, per-site-type ceiling, no score within 15 nt of the CDS).
    /// </para>
    /// </summary>
    /// <param name="mRnaSequence">The mRNA / 3'UTR sequence the site was found in (RNA or DNA; T→U; case-insensitive). Treated as the 3'UTR for the length / Min_dist / Off6m features.</param>
    /// <param name="miRna">The miRNA (its <c>Sequence</c> supplies nt1, nt8, the seed and the 3' supplementary region).</param>
    /// <param name="site">A target site whose <c>Type</c> is one of the four canonical seed-match types.</param>
    /// <param name="inputs">Optional feature inputs: ORF sequence or ORF length / ORF-8mer count, conservation (PCT), RNAplfold value (SA), and SPS / TA overrides of the bundled table.</param>
    /// <returns>The per-feature context++ contribution breakdown and their partial sum.</returns>
    /// <exception cref="ArgumentException">Thrown when <paramref name="site"/> is not one of the four canonical seed-match site types (8mer / 7mer-m8 / 7mer-A1 / 6mer).</exception>
    public static ContextPlusPlusScore ScoreTargetSiteContextPlusPlus(
        string mRnaSequence,
        MiRna miRna,
        TargetSite site,
        ContextPlusPlusInputs inputs = default)
    {
        if (site.Type is not (TargetSiteType.Seed8mer or TargetSiteType.Seed7merM8
            or TargetSiteType.Seed7merA1 or TargetSiteType.Seed6mer))
        {
            throw new ArgumentException(
                "context++ scoring is defined only for the four canonical seed-match site types " +
                "(8mer, 7mer-m8, 7mer-A1, 6mer) — Agarwal et al. (2015).", nameof(site));
        }

        string mrna = string.IsNullOrEmpty(mRnaSequence)
            ? ""
            : mRnaSequence.ToUpperInvariant().Replace('T', 'U');
        string mirna = string.IsNullOrEmpty(miRna.Sequence)
            ? ""
            : miRna.Sequence.ToUpperInvariant().Replace('T', 'U');

        TargetSiteType type = site.Type;

        double intercept = type switch
        {
            TargetSiteType.Seed8mer => CtxIntercept8mer,
            TargetSiteType.Seed7merM8 => CtxIntercept7merM8,
            TargetSiteType.Seed7merA1 => CtxIntercept7merA1,
            _ => CtxIntercept6mer
        };

        double localAu = LocalAuContribution(mrna, site.Start, site.End, type);
        double sRna1 = SRna1Contribution(mirna, type);
        double sRna8 = SRna8Contribution(mirna, type);
        double site8 = Site8Contribution(mrna, site.Start, type);
        double sa = SaContribution(mrna, site.Start, type, inputs.PlfoldUnpairedProbability, out bool saIncluded);
        double threeP = ThreePrimePairingContribution(mrna, mirna, site.Start, site.End, type);
        double minDist = MinDistContribution(mrna, site.Start, site.End, type);
        double len3Utr = Len3UtrContribution(mrna, type);
        double off6m = Off6mContribution(mrna, mirna, type);

        // readTargets: $seedRegion = substr($matureMiRNA, 1, 7) keys TA_SPS_by_seed_region.txt; TA and
        // SPS are used only when the row exists ("Set these NAs to 0"). A caller value overrides.
        string seedRegion = mirna.Length >= 8 ? mirna.Substring(1, 7) : "";
        bool haveTaSps = TryGetTargetScanTaSps(seedRegion, type, out double tableTa, out double tableSps);
        double? spsValue = inputs.Sps ?? (haveTaSps ? tableSps : null);
        double? taValue = inputs.Ta ?? (haveTaSps ? tableTa : null);

        // getORFlength_contribution / getORF8mer_contribution read the files written by
        // targetscan_count_8mers.pl from the ORF sequence; caller numbers override.
        double? orfLengthValue = inputs.OrfLength
            ?? (inputs.OrfSequence is string orfForLength ? TargetScanOrfLength(orfForLength) : null);
        int? orf8mValue = inputs.Orf8mCount
            ?? (inputs.OrfSequence is string orfForCount && seedRegion.Length == 7
                ? CountTargetScanOrf8merSites(orfForCount, seedRegion)
                : null);

        double sps = spsValue is double spsRaw ? SpsContribution(spsRaw, type) : 0.0;
        double ta = taValue is double taRaw ? TaContribution(taRaw, type) : 0.0;
        double lenOrf = orfLengthValue is double orfLen ? LenOrfContribution(orfLen, type) : 0.0;
        double orf8m = orf8mValue is int orf8mCount ? Orf8mContribution(orf8mCount, type) : 0.0;

        double bls = 0.0, pct = 0.0, pctContribution = 0.0;
        bool pctDefined = false, pctIsNa = false;
        if (inputs.Conservation is PctConservation conservation)
        {
            bls = ComputeBranchLengthScore(conservation.Tree, conservation.SpeciesWithSite);
            // getPCT: a negative or non-finite BL has no PCT ("NA"), and context++ then uses
            // PCT contribution 0 (readTargets: if ($pct eq "NA") { $PCT_contribution = 0; }).
            if (double.IsFinite(bls) && bls >= 0.0)
            {
                pctDefined = true;
                // Caller row overrides; otherwise the bundled TargetScan row for the seed (nt 2–8).
                // getPCT: a family without PCT parameters ($thisMiRNAhasPCT false) gets "NA", which
                // context++ scores as contribution 0 — TargetScan's own full score for that site.
                if (conservation.SigmoidParameters is PctSigmoidParameters callerRow)
                {
                    // getPCT: 6mer ⇒ "0.0"; the BL_PCT output's 4-dp PCT is what context++ reads.
                    pct = PctFromBranchLength(bls, callerRow, type);
                }
                else if (!IsTargetScanPctSeedFamily(seedRegion))
                {
                    pctIsNa = true;
                }
                else if (type != TargetSiteType.Seed6mer
                    && TryGetTargetScanPctParameters(seedRegion, type, out var tableRow))
                {
                    pct = PctFromBranchLength(bls, tableRow, type);
                }
                // else: 6mer of a PCT family ⇒ getPCT "0.0" (pct stays 0).
                pctContribution = PctContribution(pct, type);
            }
        }

        double partial = intercept + localAu + sRna1 + sRna8 + site8 + sa
            + threeP + minDist + len3Utr + off6m
            + sps + ta + lenOrf + orf8m + pctContribution;

        var omittedFeatures = BuildOmittedFeatures(inputs, saIncluded, pctDefined,
            spsValue.HasValue, taValue.HasValue, orfLengthValue.HasValue, orf8mValue.HasValue);

        // When any context++ feature is omitted (TA/SPS/PCT/ORF not supplied), the score is partial,
        // not the full Agarwal (2015) context++ score. Strict mode throws rather than returning it.
        if (omittedFeatures.Count > 0)
            Seqeron.Genomics.Core.LimitationPolicy.Enforce("MIRNA-TARGET-001");

        // TargetScan's reported score (readTargets): 3-dp contributions, 3-dp sum, per-type ceiling,
        // and no score for a site that starts within the first 14 nt of the UTR ("too_close").
        double? reported = site.Start + 1 < MinDistToCds
            ? null
            : ReportedContextScore(type, intercept, localAu, sRna1, sRna8, site8, sa, threeP, minDist,
                len3Utr, off6m, sps, ta, lenOrf, orf8m, pctContribution);

        return new ContextPlusPlusScore(
            SiteType: type,
            Intercept: intercept,
            LocalAuContribution: localAu,
            SRna1Contribution: sRna1,
            SRna8Contribution: sRna8,
            Site8Contribution: site8,
            SaContribution: sa,
            ThreePrimePairingContribution: threeP,
            MinDistContribution: minDist,
            Len3UtrContribution: len3Utr,
            Off6mContribution: off6m,
            SpsContribution: sps,
            TaContribution: ta,
            LenOrfContribution: lenOrf,
            Orf8mContribution: orf8m,
            PctContribution: pctContribution,
            BranchLengthScore: bls,
            Pct: pct,
            ContextScorePartial: partial,
            OmittedFeatures: omittedFeatures)
        { ContextScore = reported, PctIsNa = pctIsNa };
    }

    // targetscan_70_context_scores.pl: "too_close" when utrStart < $MIN_DIST_TO_CDS (MinDistToCds,
    // shared with AnalyzeTargetContext) and $maxContextScore{1..4} = -0.01 (7mer-1a), -0.02 (7mer-m8),
    // -0.03 (8mer-1a), 0 (6mer).

    private static double MaxContextScore(TargetSiteType type) => type switch
    {
        TargetSiteType.Seed8mer => -0.03,
        TargetSiteType.Seed7merM8 => -0.02,
        TargetSiteType.Seed7merA1 => -0.01,
        _ => 0.0
    };

    // readTargets: each contribution is the 3-dp string from getAgarwalContribution (the intercept
    // and the sRNA/Site8 coefficients are 3-dp numbers already), the total is sprintf("%.3f") of
    // their sum, and "if ($totalcontextScore > $maxContextScore{$siteType})" caps it.
    private static double ReportedContextScore(TargetSiteType type, params double[] contributions)
    {
        double sum = 0.0;
        foreach (double c in contributions)
            sum += RoundPerl3(c);
        double total = RoundPerl3(sum);
        double cap = MaxContextScore(type);
        return total > cap ? cap : total;
    }

    /// <summary>
    /// <c>sprintf("%.3f", x)</c> as perl (C printf) does it: the exact binary value rounded to 3
    /// decimals, exact ties to even. .NET Core 3.0+ "F3" formatting is IEEE-exact and gives the same
    /// digits (checked on 0.0125 → 0.013, 0.0625 → 0.062, 0.1875 → 0.188, 0.0005 → 0.001).
    /// </summary>
    internal static double RoundPerl3(double x) => RoundPerl(x, "F3");

    // sprintf("%.Nf", x) for the fixed-point format given ("F3", "F4"); see RoundPerl3.
    private static double RoundPerl(double x, string format) =>
        double.Parse(x.ToString(format, System.Globalization.CultureInfo.InvariantCulture),
            System.Globalization.CultureInfo.InvariantCulture);

    // ── TA_3UTR: target-site abundance computed from a 3'UTR set ───────────────────────────
    // Garcia et al. (2011) Nat Struct Mol Biol 18:1139 (Online Methods, retrieved verbatim this
    // session): "TA in the human transcriptome was calculated as the number of non-overlapping
    // 3′UTR 8mer, 7mer-m8, and 7mer-A1 sites in the reference mRNAs."  Agarwal et al. (2015) eLife
    // 4:e05005 (Table 1) names the context++ feature TA_3UTR = "Number of sites in all annotated
    // 3′ UTRs (Arvey et al., 2010; Garcia et al., 2011)".  TargetScan stores ONE TA value per 7-nt
    // seed region in TA_SPS_by_seed_region.txt (column 4; values ≈ 1.6–4.4, retrieved verbatim this
    // session) and feeds it AS-IS to getAgarwalContribution (targetscan_70_context_scores.pl) — the
    // value is therefore log10(count) (the regression min/max are ≈ 3.1–3.9; Garcia 2011 notes "each
    // additional CG dinucleotide imparted an additional log10 reduction in TA").  So:
    //     TA = log10( N ),  N = total non-overlapping 8mer+7mer-m8+7mer-A1 sites of the seed
    //                           summed over every 3'UTR in the supplied set.
    // The site definitions are the canonical TargetScan ones (Bartel 2009; Lewis 2005), identical to
    // FindTargetSites: an 8mer/7mer-m8/7mer-A1 is anchored on a 6mer-core match (reverse complement
    // of miRNA positions 2-7) that additionally has a position-8 match (m8) and/or an A opposite
    // position 1 (A1).  The three types are mutually exclusive per anchor (each qualifying anchor is
    // one site), and counting is left-to-right greedy so that physically overlapping core anchors
    // (possible only for self-similar / periodic cores) are counted once — honouring Garcia's
    // explicit "non-overlapping" requirement (see CountSeedSitesInUtr).

    /// <summary>
    /// Computes the TargetScan context++ feature <c>TA_3UTR</c> (transcriptome-wide target-site
    /// abundance) for a miRNA seed from a supplied set of 3'UTR sequences, per Garcia et al. (2011)
    /// / Agarwal et al. (2015): <c>TA = log10(N)</c>, where <c>N</c> is the total number of
    /// non-overlapping canonical <b>8mer</b>, <b>7mer-m8</b>, and <b>7mer-A1</b> sites of the seed
    /// across all the 3'UTRs. This converts TA from an opaque caller-supplied number into a value
    /// computed from a 3'UTR set per the published definition; the result can be fed straight into
    /// <see cref="ScoreTargetSiteContextPlusPlus"/> via <see cref="ContextPlusPlusInputs.Ta"/>
    /// (which applies the bundled Agarwal TA coefficient and min-max scaling).
    /// </summary>
    /// <remarks>
    /// Site definitions match <see cref="FindTargetSites"/> and the TargetScan precompute: each
    /// occurrence of the 6mer core (reverse complement of miRNA positions 2-7) that also has a
    /// position-8 match (7mer-m8 / 8mer) and/or an A opposite miRNA position 1 (7mer-A1 / 8mer)
    /// counts as one site. A bare 6mer (neither m8 nor A1) and offset-6mer matches are NOT counted
    /// (Garcia 2011 counts only the three high-confidence site types). Sites are counted
    /// <b>non-overlapping</b> (Garcia 2011): the scan is left-to-right greedy, so two qualifying
    /// core anchors whose 6mer-core footprints overlap (possible only for self-similar / periodic
    /// cores) are counted once.
    /// </remarks>
    /// <param name="miRna">The miRNA whose seed (positions 2-8) defines the site set. Its
    /// <c>SeedSequence</c> (or, if absent, positions 2-8 of its <c>Sequence</c>) is used.</param>
    /// <param name="threePrimeUtrs">The set of 3'UTR sequences (RNA or DNA; T→U; case-insensitive)
    /// that constitute the transcriptome over which abundance is counted.</param>
    /// <returns>The log10 of the total site count (TA_3UTR), or <c>0</c> when the total count is 0
    /// (log10(1) = 0 in TargetScan's convention; an empty/short seed also yields 0).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="threePrimeUtrs"/> is null.</exception>
    public static double ComputeTa3Utr(MiRna miRna, IEnumerable<string> threePrimeUtrs)
    {
        ArgumentNullException.ThrowIfNull(threePrimeUtrs);
        long total = CountSeedSites3Utr(miRna, threePrimeUtrs);
        // TargetScan stores log10(count). A zero count maps to 0 (= log10(1)); log10(0) is undefined,
        // and TargetScan never emits a seed with no sites, so the floor of 0 is the safe convention.
        return total > 0 ? Math.Log10(total) : 0.0;
    }

    /// <summary>
    /// Counts the total number of non-overlapping canonical <b>8mer</b>, <b>7mer-m8</b>, and
    /// <b>7mer-A1</b> sites of a miRNA seed across a set of 3'UTRs — the raw <c>N</c> underlying
    /// <see cref="ComputeTa3Utr"/> (before the log10), per Garcia et al. (2011).
    /// </summary>
    /// <param name="miRna">The miRNA whose seed defines the site set.</param>
    /// <param name="threePrimeUtrs">The 3'UTR sequences to count over.</param>
    /// <returns>The total site count (≥ 0).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="threePrimeUtrs"/> is null.</exception>
    public static long CountSeedSites3Utr(MiRna miRna, IEnumerable<string> threePrimeUtrs)
    {
        ArgumentNullException.ThrowIfNull(threePrimeUtrs);

        string seed = string.IsNullOrEmpty(miRna.SeedSequence)
            ? GetSeedSequence(miRna.Sequence ?? "")
            : miRna.SeedSequence;
        if (seed.Length < 7)
            return 0;

        // seedRC layout (7 chars): [RC of pos8][RC of pos7]…[RC of pos2]; 6mer core = seedRC[1..7].
        string seedRC = GetReverseComplement(seed);
        if (seedRC.Length < 7)
            return 0;
        char pos8Rc = seedRC[0];                  // base required opposite miRNA position 8 (upstream)
        string sixmerCore = seedRC.Substring(1, 6); // reverse complement of miRNA positions 2-7

        long total = 0;
        foreach (string? utr in threePrimeUtrs)
        {
            if (string.IsNullOrEmpty(utr))
                continue;
            total += CountSeedSitesInUtr(
                utr.ToUpperInvariant().Replace('T', 'U'), sixmerCore, pos8Rc);
        }
        return total;
    }

    // Counts 8mer + 7mer-m8 + 7mer-A1 sites in ONE 3'UTR. Scan for 6mer-core occurrences; a site
    // counts iff it additionally has the position-8 match (m8) or the A opposite position 1 (A1) —
    // i.e. it is an 8mer (both), a 7mer-m8 (m8 only) or a 7mer-A1 (A1 only). A bare 6mer (neither) is
    // NOT a counted site type (Garcia 2011 counts only the three).
    //
    // Garcia (2011) Online Methods: TA = number of *non-overlapping* sites. Distinct qualifying
    // core anchors can still physically overlap when the 6mer core is self-similar (e.g. a
    // low-complexity / periodic core against a homopolymeric UTR), so counting every anchor would
    // over-count relative to the published definition. We therefore scan left-to-right greedily:
    // when a site is counted at i, the next candidate core must start at i+6 or later (its 6mer-core
    // footprint must not overlap an already-counted site). On non-self-overlapping cores (real miRNA
    // seeds against typical 3'UTRs) this is identical to per-anchor counting; on periodic cores it
    // enforces the non-overlapping requirement.
    internal static int CountSeedSitesInUtr(string utr, string sixmerCore, char pos8Rc)
    {
        int count = 0;
        int i = 0;
        while (i + 6 <= utr.Length)
        {
            if (string.CompareOrdinal(utr, i, sixmerCore, 0, 6) != 0)
            {
                i++;
                continue;
            }
            bool hasPos8 = i > 0 && utr[i - 1] == pos8Rc;     // 7mer-m8 / 8mer
            bool hasA1 = i + 6 < utr.Length && utr[i + 6] == 'A'; // 7mer-A1 / 8mer
            if (hasPos8 || hasA1)
            {
                count++;
                i += 6; // skip past this site's 6mer-core footprint — sites are non-overlapping
            }
            else
            {
                i++;
            }
        }
        return count;
    }

    /// <summary>
    /// context++ features that remain residual for a given call. SA is omitted when the 14-nt
    /// RNAplfold window has no row at this site (too close to the 5' end) or its unpaired probability
    /// is 0 (perl: contribution 0). PCT is omitted unless a conservation input is supplied, or when
    /// its branch-length score is negative / not finite (TargetScan "NA"). SPS and TA are omitted only
    /// when neither the caller nor the bundled table gives a value (seed region not 7 nt over ACGU);
    /// ORF length and ORF-8mer count when neither a number nor an ORF sequence is supplied.
    /// </summary>
    private static IReadOnlyList<string> BuildOmittedFeatures(ContextPlusPlusInputs inputs, bool saIncluded, bool pctDefined,
        bool haveSps, bool haveTa, bool haveOrfLength, bool haveOrf8m)
    {
        var omitted = new List<string>();
        if (!saIncluded)
            omitted.Add("SA (structural accessibility — no 14-nt window at this site, or unpaired probability 0)");
        if (inputs.Conservation is null)
            omitted.Add("PCT (probability of conserved targeting — requires a multi-species alignment + tree)");
        else if (!pctDefined)
            omitted.Add("PCT (probability of conserved targeting — NA: the branch-length score is negative or not finite)");
        if (!haveSps) omitted.Add("SPS (predicted seed-pairing stability — seed region nt 2-8 is not in the TargetScan table; supply Sps)");
        if (!haveTa) omitted.Add("TA_3UTR (target-site abundance — seed region nt 2-8 is not in the TargetScan table; supply Ta)");
        if (!haveOrfLength) omitted.Add("Len_ORF (ORF length — requires the transcript ORF: supply OrfSequence or OrfLength)");
        if (!haveOrf8m) omitted.Add("ORF8m (ORF 8mer count — requires the transcript ORF: supply OrfSequence or Orf8mCount)");
        return omitted;
    }

    // Local_AU: faithful port of getLocalAU_contribution. The local-AU fraction is the
    // position-weighted A/U content of the 30 nt up- and downstream of the site; it is then
    // min-max scaled and multiplied by the site-type Local_AU coefficient. Upstream positions
    // are weighted from the base immediately 5' of the site (i=0) outward; downstream from the
    // base immediately 3' of the site. Weighting offsets are +1 vs +2 per site type, exactly as
    // in the perl (8mer/7mer-m8 favour the upstream by one rank; 8mer/7mer-A1 the downstream).
    internal static double LocalAuContribution(string mrna, int siteStart, int siteEnd, TargetSiteType type)
    {
        // Perl uses 1-based utrStart/utrEnd; here Start/End are 0-based inclusive site coordinates.
        // 8mer (Seed8mer) and 7mer-m8 weight upstream by 1/(i+1); 7mer-A1 and 6mer by 1/(i+2).
        // 8mer and 7mer-A1 weight downstream by 1/(i+2); 7mer-m8 and 6mer by 1/(i+1).
        int upOffset = (type is TargetSiteType.Seed8mer or TargetSiteType.Seed7merM8) ? 1 : 2;
        int downOffset = (type is TargetSiteType.Seed8mer or TargetSiteType.Seed7merA1) ? 2 : 1;
        double? weighted = LocalAuFraction(mrna, siteStart, siteEnd, LocalAuFlankLength, upOffset, downOffset);
        if (weighted is not double fraction)
            return 0.0;

        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxLocalAuCoeff8mer, CtxLocalAuMin8mer, CtxLocalAuMax8mer),
            TargetSiteType.Seed7merM8 => (CtxLocalAuCoeff7merM8, CtxLocalAuMin7merM8, CtxLocalAuMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxLocalAuCoeff7merA1, CtxLocalAuMin7merA1, CtxLocalAuMax7merA1),
            _ => (CtxLocalAuCoeff6mer, CtxLocalAuMin6mer, CtxLocalAuMax6mer)
        };

        // The perl rounds the fraction before scaling: $totalLocalAUscoreFraction =
        // sprintf("%.3f", $totalLocalAUscoreRaw / $localAUcontributionMaxRaw).
        double scaled = (RoundPerl3(fraction) - min) / (max - min);
        return coeff * scaled;
    }

    /// <summary>
    /// TargetScan/Grimson (2007) position-weighted local A/U fraction — the raw value of
    /// <c>getLocalAU_contribution</c> (<c>targetscan_70_context_scores.pl</c>): up to
    /// <paramref name="flank"/> nt immediately 5' of the site (walked outward from
    /// <paramref name="siteStart"/>−1, weight 1/(i+<paramref name="upOffset"/>)) and up to
    /// <paramref name="flank"/> nt immediately 3' (from <paramref name="siteEnd"/>+1, weight
    /// 1/(i+<paramref name="downOffset"/>)); fraction = Σ weights at A/U ÷ Σ all weights. The site
    /// itself is excluded. Returns null when there is no flanking nucleotide at all.
    /// Expects an upper-case RNA sequence.
    /// </summary>
    private static double? LocalAuFraction(string mrna, int siteStart, int siteEnd, int flank, int upOffset, int downOffset)
    {
        double scoreSum = 0.0;
        double maxRaw = 0.0;

        for (int i = 0; i < flank; i++)
        {
            int idx = siteStart - 1 - i;
            if (idx < 0) break;
            double weight = 1.0 / (i + upOffset);
            char b = mrna[idx];
            if (b == 'A' || b == 'U')
                scoreSum += weight;
            maxRaw += weight;
        }

        for (int i = 0; i < flank; i++)
        {
            int idx = siteEnd + 1 + i;
            if (idx >= mrna.Length) break;
            double weight = 1.0 / (i + downOffset);
            char b = mrna[idx];
            if (b == 'A' || b == 'U')
                scoreSum += weight;
            maxRaw += weight;
        }

        return maxRaw == 0.0 ? null : scoreSum / maxRaw;
    }

    // sRNA position-1 indicators: contributions are 0 when miRNA nt1 is U (perl: only computed
    // when sRNA1_nt ne "U"). Otherwise the indicator for the actual nt (A/C/G) is 1, others 0.
    internal static double SRna1Contribution(string mirna, TargetSiteType type)
    {
        if (mirna.Length < 1) return 0.0;
        char nt1 = mirna[0];
        if (nt1 == 'U') return 0.0;

        (double a, double c, double g) = type switch
        {
            TargetSiteType.Seed8mer => (CtxSRna1A8mer, CtxSRna1C8mer, CtxSRna1G8mer),
            TargetSiteType.Seed7merM8 => (CtxSRna1A7merM8, CtxSRna1C7merM8, CtxSRna1G7merM8),
            TargetSiteType.Seed7merA1 => (CtxSRna1A7merA1, CtxSRna1C7merA1, CtxSRna1G7merA1),
            _ => (CtxSRna1A6mer, CtxSRna1C6mer, CtxSRna1G6mer)
        };

        return nt1 switch { 'A' => a, 'C' => c, 'G' => g, _ => 0.0 };
    }

    // sRNA position-8 indicators: 0 when miRNA nt8 is U; else indicator for the actual A/C/G.
    internal static double SRna8Contribution(string mirna, TargetSiteType type)
    {
        if (mirna.Length < 8) return 0.0;
        char nt8 = mirna[7];
        if (nt8 == 'U') return 0.0;

        (double a, double c, double g) = type switch
        {
            TargetSiteType.Seed8mer => (CtxSRna8A8mer, CtxSRna8C8mer, CtxSRna8G8mer),
            TargetSiteType.Seed7merM8 => (CtxSRna8A7merM8, CtxSRna8C7merM8, CtxSRna8G7merM8),
            TargetSiteType.Seed7merA1 => (CtxSRna8A7merA1, CtxSRna8C7merA1, CtxSRna8G7merA1),
            _ => (CtxSRna8A6mer, CtxSRna8C6mer, CtxSRna8G6mer)
        };

        return nt8 switch { 'A' => a, 'C' => c, 'G' => g, _ => 0.0 };
    }

    // Site8 indicators: only defined for 7mer-A1 and 6mer (the perl computes them only for those
    // site types). The relevant target base is the nucleotide opposite miRNA position 8, i.e. the
    // base immediately 5' of the 6mer core. For 7mer-A1 / 6mer the 6mer core starts at site.Start,
    // so the position-8 base is mrna[site.Start - 1]. 0 when that base is U or out of range.
    internal static double Site8Contribution(string mrna, int siteStart, TargetSiteType type)
    {
        if (type is not (TargetSiteType.Seed7merA1 or TargetSiteType.Seed6mer))
            return 0.0;

        int idx = siteStart - 1;
        if (idx < 0 || idx >= mrna.Length) return 0.0;
        char sitePos8 = mrna[idx];
        if (sitePos8 == 'U') return 0.0;

        (double a, double c, double g) = type == TargetSiteType.Seed7merA1
            ? (CtxSite8A7merA1, CtxSite8C7merA1, CtxSite8G7merA1)
            : (CtxSite8A6mer, CtxSite8C6mer, CtxSite8G6mer);

        return sitePos8 switch { 'A' => a, 'C' => c, 'G' => g, _ => 0.0 };
    }

    // ── SA: structural accessibility (getSA_contribution) ─────────────────────────────────
    // getSA_contribution: TargetScan runs "RNAplfold -L 40 -W 80 -u 20" on the whole UTR
    // (runRNAplfold_all_UTRs) and reads `grep -A7 -P '^$utrStart\t' _lunp | tail -1 | cut -f 2-`,
    // index 13: the l = 14 column of the row 7 nt downstream of the seed-match start (utrStart − 1
    // for 7mer-A1), i.e. the probability that the 14 nt ENDING at that row are unpaired (RNAplfold
    // man page; Bernhart et al. 2006). When fewer than 7 rows follow, `tail -1` returns the UTR's
    // last row. A missing / zero / NA value ("! $plfold") gives contribution 0. Otherwise
    // log10(plfold) is min-max scaled by the site-type SA coefficient.
    // The value itself is either caller-supplied (an RNAplfold _lunp value → identical to TargetScan)
    // or estimated here by the Turner-2004 McCaskill partition function on ONE 80-nt window centred
    // on the 14-nt region (RnaSecondaryStructure.CalculateRegionUnpairedProbability). The estimate is
    // NOT RNAplfold: RNAplfold averages the window partition functions of every 80-nt window that
    // holds the region and limits base pairs to a 40-nt span; the canonical partition function has
    // no span limit, so that average cannot be computed here (request to RNA-STRUCT, see report).
    // Sources: targetscan_70_context_scores.pl getSA_contribution / runRNAplfold_all_UTRs;
    // RNAplfold man page + Bernhart et al. (2006) Bioinformatics 22:614; Agarwal et al. (2015)
    // eLife 4:e05005 Fig 4A ("log10 value of the unpaired probability for a 14-nt region centered
    // on the match to miRNA nucleotides 7 and 8").
    internal static double SaContribution(string mrna, int siteStart, TargetSiteType type, out bool included)
        => SaContribution(mrna, siteStart, type, null, out included);

    internal static double SaContribution(string mrna, int siteStart, TargetSiteType type, double? suppliedPlfold, out bool included)
    {
        included = false;

        double plfold;
        if (suppliedPlfold is double supplied)
        {
            plfold = supplied;
        }
        else
        {
            // 1-based UTR start of the seed match; perl decrements it for 7mer-A1 (siteType 1) so the
            // positioning matches 8mer / 7mer-m8 sites (the legacy "|| siteType == 5" is dead code).
            int utrStart1 = siteStart + 1;
            if (type == TargetSiteType.Seed7merA1)
                utrStart1--;

            // grep -A7 … | tail -1: row utrStart + 7, or the last row when the UTR ends sooner.
            int row1 = Math.Min(utrStart1 + SaRowOffsetDownstream, mrna.Length);
            int windowEnd0 = row1 - 1;
            int windowStart0 = windowEnd0 - SaUnpairedWindowLength + 1;
            if (utrStart1 < 1 || windowStart0 < 0)
                return 0.0; // the l = 14 column is NA on rows < 14 → plfold "NA" → 0 (SA omitted)

            plfold = LocalRegionUnpairedProbability(mrna, windowStart0, windowEnd0);
        }

        // perl: if (! $plfold || $plfold eq "NA") { return 0; }
        if (!(plfold > 0))
            return 0.0;

        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxSaCoeff8mer, CtxSaMin8mer, CtxSaMax8mer),
            TargetSiteType.Seed7merM8 => (CtxSaCoeff7merM8, CtxSaMin7merM8, CtxSaMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxSaCoeff7merA1, CtxSaMin7merA1, CtxSaMax7merA1),
            _ => (CtxSaCoeff6mer, CtxSaMin6mer, CtxSaMax6mer)
        };

        included = true;
        return ScaledContribution(Math.Log10(plfold), coeff, min, max);
    }

    /// <summary>
    /// Probability that the region [<paramref name="windowStart0"/>..<paramref name="windowEnd0"/>]
    /// (0-based, inclusive) is entirely unpaired, from the Turner-2004 McCaskill partition function
    /// (canonical <see cref="RnaSecondaryStructure.CalculateRegionUnpairedProbability"/>, Z_open/Z),
    /// folded over a local context of W = max(80, region length) nt centred on the region and clamped
    /// to the sequence — the RNAplfold <c>-W 80</c> local-folding window used by TargetScan
    /// (<c>runRNAplfold_all_UTRs</c>; Bernhart et al. 2006). Caller guarantees the region is in range.
    /// </summary>
    private static double LocalRegionUnpairedProbability(string mrna, int windowStart0, int windowEnd0)
    {
        int regionLength = windowEnd0 - windowStart0 + 1;
        int foldWindow = Math.Max(SaPlfoldWindowSize, regionLength);

        // Local fold context: up to W nt centred on the region, clamped to the sequence. RNAplfold
        // averages over length-W windows; folding this local context captures the local
        // accessibility of the region (a base can only pair within the local window).
        int contextStart = Math.Max(0, windowEnd0 - (foldWindow - regionLength) / 2 - regionLength + 1);
        contextStart = Math.Min(contextStart, windowStart0);
        int contextEnd = Math.Min(mrna.Length - 1, contextStart + foldWindow - 1);
        contextStart = Math.Max(0, contextEnd - foldWindow + 1);
        string context = mrna.Substring(contextStart, contextEnd - contextStart + 1);
        int localWindowEnd = windowEnd0 - contextStart;

        return RnaSecondaryStructure.CalculateRegionUnpairedProbability(context, localWindowEnd, regionLength);
    }

    // ── Generic min-max scaling (getAgarwalContribution) ──────────────────────────────────
    // Perl: scaledScore = (raw - min) / (max - min); contribution = coeff × scaledScore.
    // NOT clamped to [0,1] (scaled values < 0 or > 1 are kept, exactly as in the perl).
    internal static double ScaledContribution(double raw, double coeff, double min, double max)
        => max == min ? coeff * raw : coeff * ((raw - min) / (max - min));

    // ── 3' supplementary pairing (3P_score) ───────────────────────────────────────────────
    // Faithful port of get3primePairingContribution, fed by extractSubseqForAlignment +
    // modifySubseqForAlignment (targetscan_70_context_scores.pl). The raw score is the maximum,
    // over all single-gap offsets in both the UTR and the miRNA, of: 1.0 per contiguous-stretch
    // base pair in the 3' "seed-supplementary" window (offset-adjusted positions 4..7) and 0.5
    // elsewhere, summed only over runs of ≥2 consecutive pairs, minus max(0,(offset-2)/2). The
    // raw score is then min-max scaled (min=1, max=3.5) and multiplied by the 3P_score coeff.
    private static readonly int[] ThreePrimeUtrStart  = { 0, 8, 8, 8, 8 }; // index by perl siteType 1..4
    private static readonly int[] ThreePrimeMiRnaStart = { 0, 7, 8, 8, 8 };
    private static readonly int[] ThreePrimeOverhang   = { 0, 1, 0, 0, 0 };

    internal static double ThreePrimePairingContribution(string mrna, string mirna, int siteStart, int siteEnd, TargetSiteType type)
    {
        double raw = ThreePrimePairingRaw(mrna, mirna, siteStart, siteEnd, type);
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (Ctx3PCoeff8mer, Ctx3PMin, Ctx3PMax),
            TargetSiteType.Seed7merM8 => (Ctx3PCoeff7merM8, Ctx3PMin, Ctx3PMax),
            TargetSiteType.Seed7merA1 => (Ctx3PCoeff7merA1, Ctx3PMin, Ctx3PMax),
            _ => (Ctx3PCoeff6mer, Ctx3PMin, Ctx3PMax)
        };
        return ScaledContribution(raw, coeff, min, max);
    }

    // Perl site-type integers: 1 = 7mer-A1, 2 = 7mer-m8, 3 = 8mer, 4 = 6mer.
    private static int PerlSiteType(TargetSiteType type) => type switch
    {
        TargetSiteType.Seed7merA1 => 1,
        TargetSiteType.Seed7merM8 => 2,
        TargetSiteType.Seed8mer => 3,
        _ => 4
    };

    // Returns the raw 3P pairing score (the value getAgarwalContribution receives for "3P_score").
    internal static double ThreePrimePairingRaw(string mrna, string mirna, int siteStart, int siteEnd, TargetSiteType type)
    {
        if (mrna.Length == 0 || mirna.Length == 0) return 0.0;

        int perlType = PerlSiteType(type);
        // C# 0-based inclusive Start/End ↔ perl 1-based utrStart/utrEnd.
        int utrStart = siteStart + 1;
        int utrEnd = siteEnd + 1;

        string subseq = ExtractSubseqForAlignment(mrna, utrStart, utrEnd, perlType);
        (string finalUtr, string mirnaForAlign) = ModifySubseqForAlignment(mirna, subseq, perlType);
        return ThreePrimePairingScore(perlType, finalUtr, mirnaForAlign);
    }

    // Port of extractSubseqForAlignment (siteType < 5 path only; the four canonical types).
    private static string ExtractSubseqForAlignment(string utr, int utrStart, int utrEnd, int perlType)
    {
        int realStart = utrStart - ThreePrimeUtrStartOffset;
        if (realStart < 0) realStart = 0;

        // real_end = utrEnd + 1 for siteType 1, 2, 4; else utrEnd.
        int realEnd = (perlType == 1 || perlType == 2 || perlType == 4) ? utrEnd + 1 : utrEnd;
        if (realStart >= realEnd) realStart = 0;

        int length = realEnd - realStart;
        // substr in 0-based: perl coordinates are 1-based, so the 0-based offset is realStart-1...
        // but the perl applies substr($seq, $real_start, $length) with $real_start already a
        // sequence offset (it is utrStart-16 where utrStart is 1-based, so it doubles as 0-based
        // start of the window). We mirror the perl arithmetic exactly using 0-based realStart.
        int begin = realStart;
        if (begin < 0) begin = 0;
        int avail = Math.Max(0, utr.Length - begin);
        int take = Math.Min(length, avail);
        string sub = take > 0 ? utr.Substring(begin, take) : "";

        // Pad leading N's up to DESIRED_UTR_ALIGNMENT_LENGTH (perl loop adds 23-len+1 N's).
        if (DesiredUtrAlignmentLength > sub.Length)
        {
            int nCount = DesiredUtrAlignmentLength - sub.Length + 1;
            sub = new string('N', nCount) + sub;
        }

        // Site at the very end of the UTR: perl appends trailing N's and converts a leading "NN".
        if (utr.Length < (begin + length))
        {
            int lengthDiff = begin + length - utr.Length;
            var sb = new StringBuilder(sub);
            for (int i = 0; i < lengthDiff; i++)
            {
                sb.Append('N');
                if (sb.Length >= 2 && sb[0] == 'N' && sb[1] == 'N')
                {
                    sb[0] = ' ';
                    sb[1] = ' ';
                }
            }
            sub = sb.ToString();
        }
        return sub;
    }

    // Port of modifySubseqForAlignment.
    private static (string finalUtr, string mirnaForAlign) ModifySubseqForAlignment(string matureMiRna, string subseq, int perlType)
    {
        int spacer1Length = matureMiRna.Length - DesiredUtrAlignmentLength;
        string spacer1 = spacer1Length > 0 ? new string(' ', spacer1Length) : "";

        int spacer2Length = spacer1Length < 0 ? -spacer1Length : 0;
        var spacer2 = new StringBuilder(new string(' ', spacer2Length));

        if (DesiredUtrAlignmentLength > subseq.Length)
        {
            int extra = DesiredUtrAlignmentLength - subseq.Length + 1;
            spacer2.Append(' ', extra);
        }

        // 7mer-A1 (1) and 8mer (3): chop one char off spacer2.
        if ((perlType == 1 || perlType == 3) && spacer2.Length > 0)
            spacer2.Length -= 1;

        string finalUtr = spacer1 + subseq;
        char[] rev = matureMiRna.ToCharArray();
        Array.Reverse(rev);
        string mirnaForAlign = spacer2.ToString() + new string(rev);
        return (finalUtr, mirnaForAlign);
    }

    // Port of get3primePairingContribution — returns the raw best score only.
    internal static double ThreePrimePairingScore(int perlType, string utrIn, string mirnaIn)
    {
        // Normalise: strip spaces/newlines, uppercase, T→U.
        string utr = NormaliseForPairing(utrIn);
        string mirna = NormaliseForPairing(mirnaIn);

        // Reverse both (perl: $utr = reverse $utr; $mirna = reverse $mirna).
        utr = Reverse(utr);
        mirna = Reverse(mirna);

        int utrStart = ThreePrimeUtrStart[perlType];
        int mirnaStart = ThreePrimeMiRnaStart[perlType];
        int overhang = ThreePrimeOverhang[perlType];

        string utrNum = utrStart < utr.Length ? utr.Substring(utrStart) : "";
        string mirnaNum = mirnaStart < mirna.Length ? mirna.Substring(mirnaStart) : "";

        int maxScoreLen = Math.Max(utrNum.Length, mirnaNum.Length);

        int[] UTR = ToBaseCodes(utrNum);
        int[] MIRNA = ToBaseCodes(mirnaNum);

        double bestScore = double.NegativeInfinity;

        for (int offset = 0; offset < maxScoreLen; offset++)
        {
            // "top" gap orientation: MIRNA shifted by offset.
            double topScore = PairingRunScore(UTR, MIRNA, offset, overhang, gapOnTop: true);
            topScore -= Math.Max(0.0, (offset - 2) / 2.0);
            if (topScore > bestScore) bestScore = topScore;

            // "bottom" gap orientation: UTR shifted by offset.
            double bottomScore = PairingRunScore(UTR, MIRNA, offset, overhang, gapOnTop: false);
            bottomScore -= Math.Max(0.0, (offset - 2) / 2.0);
            if (bottomScore > bestScore) bestScore = bottomScore;
        }

        return double.IsNegativeInfinity(bestScore) ? 0.0 : bestScore;
    }

    // One alignment offset: scores only runs of ≥2 consecutive base pairs, taking the best such run.
    internal static double PairingRunScore(int[] UTR, int[] MIRNA, int offset, int overhang, bool gapOnTop)
    {
        double score = 0.0;     // best committed run score for this offset
        double tempscore = 0.0; // current run score
        int prevmatch = 0;      // current run length

        int limit = gapOnTop
            ? Math.Min(MIRNA.Length - 1 - offset, UTR.Length - 1)
            : Math.Min(UTR.Length - 1 - offset, MIRNA.Length - 1);

        for (int i = 0; i <= limit; i++)
        {
            int u = gapOnTop ? UTR[i] : UTR[i + offset];
            int m = gapOnTop ? MIRNA[i + offset] : MIRNA[i];
            // Watson-Crick A:U/U:A → product 1×2 or 2×1 = 2; G:C/C:G → 3×4 or 4×3 = 12.
            bool isPair = (u * m == 2) || (u * m == 12);

            if (isPair)
            {
                int posCheck = gapOnTop ? (i + offset - overhang) : (i - overhang);
                if (prevmatch == 0) tempscore = 0.0;
                tempscore += (posCheck >= 4 && posCheck <= 7) ? 1.0 : 0.5;
                prevmatch++;
            }
            else
            {
                if (prevmatch >= 2 && tempscore > score) score = tempscore;
                tempscore = 0.0;
                prevmatch = 0;
            }
        }
        if (prevmatch >= 2 && tempscore > score) score = tempscore;
        return score;
    }

    private static string NormaliseForPairing(string s)
    {
        var sb = new StringBuilder(s.Length);
        foreach (char ch in s)
        {
            if (ch == ' ' || ch == '\n' || ch == '\r') continue;
            char u = char.ToUpperInvariant(ch);
            sb.Append(u == 'T' ? 'U' : u);
        }
        return sb.ToString();
    }

    private static string Reverse(string s)
    {
        char[] c = s.ToCharArray();
        Array.Reverse(c);
        return new string(c);
    }

    // tr/AUCGN/12345/ — any other character maps to 0 (cannot pair).
    private static int[] ToBaseCodes(string s)
    {
        var codes = new int[s.Length];
        for (int i = 0; i < s.Length; i++)
            codes[i] = s[i] switch { 'A' => 1, 'U' => 2, 'C' => 3, 'G' => 4, 'N' => 5, _ => 0 };
        return codes;
    }

    // ── Min_dist ───────────────────────────────────────────────────────────────────────────
    // Port of getMinDist_weighted_contribution for the single-isoform case (AIR_end = UTR length).
    // distTo5 = siteStart-1 (perl 1-based) ; distTo3 = UTRlen - siteEnd ; take the smaller, log10,
    // then min-max scale by the site-type Min_dist coefficients.
    internal static double MinDistContribution(string mrna, int siteStart, int siteEnd, TargetSiteType type)
    {
        if (mrna.Length == 0) return 0.0;
        // perl 1-based site coordinates and UTR length (= AIR_end).
        int perlStart = siteStart + 1;
        int perlEnd = siteEnd + 1;
        int airEnd = mrna.Length;
        if (airEnd < perlEnd) return 0.0; // UTR shorter than the site (defensive)

        int distTo5 = perlStart - 1;
        int distTo3 = airEnd - perlEnd;
        int nearest = Math.Min(distTo5, distTo3);

        double log10 = nearest > 0 ? Math.Log10(nearest) : 0.0;

        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxMinDistCoeff8mer, CtxMinDistMin8mer, CtxMinDistMax8mer),
            TargetSiteType.Seed7merM8 => (CtxMinDistCoeff7merM8, CtxMinDistMin7merM8, CtxMinDistMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxMinDistCoeff7merA1, CtxMinDistMin7merA1, CtxMinDistMax7merA1),
            _ => (CtxMinDistCoeff6mer, CtxMinDistMin6mer, CtxMinDistMax6mer)
        };
        return ScaledContribution(log10, coeff, min, max);
    }

    // ── Len_3UTR ──────────────────────────────────────────────────────────────────────────
    // Port of get_len3UTR_weighted_contribution (single-isoform: UTR length = AIR_end = mrna.Length).
    internal static double Len3UtrContribution(string mrna, TargetSiteType type)
    {
        int utrLength = mrna.Length;
        double log10 = utrLength > 0 ? Math.Log10(utrLength) : 0.0;

        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxLen3UtrCoeff8mer, CtxLen3UtrMin8mer, CtxLen3UtrMax8mer),
            TargetSiteType.Seed7merM8 => (CtxLen3UtrCoeff7merM8, CtxLen3UtrMin7merM8, CtxLen3UtrMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxLen3UtrCoeff7merA1, CtxLen3UtrMin7merA1, CtxLen3UtrMax7merA1),
            _ => (CtxLen3UtrCoeff6mer, CtxLen3UtrMin6mer, CtxLen3UtrMax6mer)
        };
        return ScaledContribution(log10, coeff, min, max);
    }

    // ── Off6m ─────────────────────────────────────────────────────────────────────────────
    // Port of getOffset6merSites + getOffset6mer_weighted_contribution (single-isoform). The
    // offset-6mer target pattern is the first 6 nt of reverse_complement(seed region 2-8); Off6m is
    // the count of (case-insensitive, non-overlapping) occurrences of that 6mer in the 3'UTR. Used RAW
    // (Off6m is not in the min-max regex of getAgarwalContribution).
    // Note: the public targetscan_70_context_scores.pl matches the pattern against
    // `$utr_seq{$transcriptID}` — a hash reference (the UTRs keyed by species), not the sequence — and
    // stores the site list in scalar context, so it reports Off6m = 1 for every site whatever the UTR
    // holds (reproduced: 4 CUACCU in a 295-nt UTR → perl Off6m −0.020 for an 8mer). This port counts
    // the occurrences, as Agarwal et al. (2015) define the feature ("number of offset-6mer sites in
    // the 3' UTR") and as the subroutine evidently intends.
    internal static double Off6mContribution(string mrna, string mirna, TargetSiteType type)
    {
        if (mrna.Length == 0 || mirna.Length < 8) return 0.0;

        string seedRegion = mirna.Substring(1, 7);          // miRNA nt 2-8
        string revComp = GetReverseComplement(seedRegion);  // 7 nt (canonical RNA reverse complement)
        string pattern = revComp.Substring(0, 6);           // first 6 nt

        // perl `while ($utr =~ /$siteToFind/gi)`: a /g match resumes after the previous match, so
        // occurrences are counted NON-overlapping (only differs for self-overlapping patterns such as
        // ACACAC or AAAAAA).
        int count = CountNonOverlapping(mrna, pattern);

        double coeff = type switch
        {
            TargetSiteType.Seed8mer => CtxOff6mCoeff8mer,
            TargetSiteType.Seed7merM8 => CtxOff6mCoeff7merM8,
            TargetSiteType.Seed7merA1 => CtxOff6mCoeff7merA1,
            _ => CtxOff6mCoeff6mer
        };
        return coeff * count; // used raw (no scaling)
    }

    // Number of non-overlapping occurrences of `pattern` in `text` (ordinal), scanning left to right
    // and resuming after each match — what a perl `while ($text =~ /$pattern/g)` loop counts.
    private static int CountNonOverlapping(string text, string pattern)
    {
        int count = 0;
        for (int i = 0; i + pattern.Length <= text.Length;)
        {
            if (string.CompareOrdinal(text, i, pattern, 0, pattern.Length) == 0)
            {
                count++;
                i += pattern.Length;
            }
            else
            {
                i++;
            }
        }
        return count;
    }

    // ── Len_ORF / ORF8m from an ORF sequence — TargetScan 7.0 targetscan_count_8mers.pl ──────
    // (TargetScan7_context_scores.zip, www.targetscan.org/vert_70/vert_70_data_download/, opened
    // 2026-11-25; README_CONTEXT_SCORES_70: "The ORF sequences can be aligned (with gaps) or not, but the
    // alignment is ignored; gaps are removed.") Per ORF row: `$sequence =~ s/-//g; s/\.//g; s/T/U/g;
    // s/t/u/g`; length = `length $sequence`; per family `$siteToFind = reverse_complement(seed region) .
    // "A"` (the 8mer: match to miRNA nt 2-8 + A opposite nt 1) and `while ($sequence =~ /$siteToFind/gi)
    // { $numThisSite++ }` — case-insensitive, non-overlapping, on the given (sense) strand only.
    // targetscan_70_context_scores.pl then uses log10(length) for Len_ORF (0 when the length is 0) and
    // the raw count for ORF8m (0 when the family has no row).

    // targetscan_count_8mers.pl: length of the ORF row after removing '-' and '.' gap characters.
    internal static int TargetScanOrfLength(string orfSequence)
    {
        ArgumentNullException.ThrowIfNull(orfSequence);
        int length = 0;
        foreach (char c in orfSequence)
            if (c is not ('-' or '.'))
                length++;
        return length;
    }

    /// <summary>
    /// Counts the 8mer sites of a miRNA seed region in an ORF (CDS) sequence exactly as TargetScan 7.0's
    /// <c>targetscan_count_8mers.pl</c> does for the context++ feature ORF8m (Agarwal et al. 2015):
    /// gap characters (<c>-</c>, <c>.</c>) are removed, T is read as U, and the site
    /// reverse_complement(seed region) + "A" is counted case-insensitively and non-overlapping
    /// (perl <c>/g</c> match), on the given strand only.
    /// </summary>
    /// <param name="orfSequence">The ORF sequence (RNA or DNA, any case, gaps allowed).</param>
    /// <param name="seedRegion">The 7-nt seed region (miRNA nt 2–8), RNA or DNA alphabet.</param>
    /// <returns>The number of 8mer sites (≥ 0).</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="seedRegion"/> is not 7 nt long.</exception>
    public static int CountTargetScanOrf8merSites(string orfSequence, string seedRegion)
    {
        ArgumentNullException.ThrowIfNull(orfSequence);
        ArgumentNullException.ThrowIfNull(seedRegion);
        if (seedRegion.Length != 7)
            throw new ArgumentException("The seed region (miRNA nt 2-8) must be 7 nt long.", nameof(seedRegion));

        string site = GetReverseComplement(seedRegion.ToUpperInvariant().Replace('T', 'U')) + "A";
        var orf = new System.Text.StringBuilder(orfSequence.Length);
        foreach (char c in orfSequence)
        {
            if (c is '-' or '.')
                continue;
            char u = char.ToUpperInvariant(c);
            orf.Append(u == 'T' ? 'U' : u);
        }
        return CountNonOverlapping(orf.ToString(), site);
    }

    // ── SPS / TA / Len_ORF / ORF8m contributions ─────────────────────────────────────────────
    internal static double SpsContribution(double sps, TargetSiteType type)
    {
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxSpsCoeff8mer, CtxSpsMin8mer, CtxSpsMax8mer),
            TargetSiteType.Seed7merM8 => (CtxSpsCoeff7merM8, CtxSpsMin7merM8, CtxSpsMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxSpsCoeff7merA1, CtxSpsMin7merA1, CtxSpsMax7merA1),
            _ => (CtxSpsCoeff6mer, CtxSpsMin6mer, CtxSpsMax6mer)
        };
        return ScaledContribution(sps, coeff, min, max);
    }

    internal static double TaContribution(double ta, TargetSiteType type)
    {
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxTaCoeff8mer, CtxTaMin8mer, CtxTaMax8mer),
            TargetSiteType.Seed7merM8 => (CtxTaCoeff7merM8, CtxTaMin7merM8, CtxTaMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxTaCoeff7merA1, CtxTaMin7merA1, CtxTaMax7merA1),
            _ => (CtxTaCoeff6mer, CtxTaMin6mer, CtxTaMax6mer)
        };
        return ScaledContribution(ta, coeff, min, max);
    }

    internal static double LenOrfContribution(double orfLength, TargetSiteType type)
    {
        double log10 = orfLength > 0 ? Math.Log10(orfLength) : 0.0;
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxLenOrfCoeff8mer, CtxLenOrfMin8mer, CtxLenOrfMax8mer),
            TargetSiteType.Seed7merM8 => (CtxLenOrfCoeff7merM8, CtxLenOrfMin7merM8, CtxLenOrfMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxLenOrfCoeff7merA1, CtxLenOrfMin7merA1, CtxLenOrfMax7merA1),
            _ => (CtxLenOrfCoeff6mer, CtxLenOrfMin6mer, CtxLenOrfMax6mer)
        };
        return ScaledContribution(log10, coeff, min, max);
    }

    internal static double Orf8mContribution(int count, TargetSiteType type)
    {
        double coeff = type switch
        {
            TargetSiteType.Seed8mer => CtxOrf8mCoeff8mer,
            TargetSiteType.Seed7merM8 => CtxOrf8mCoeff7merM8,
            TargetSiteType.Seed7merA1 => CtxOrf8mCoeff7merA1,
            _ => CtxOrf8mCoeff6mer
        };
        return coeff * count; // used raw (no scaling)
    }

    // ── PCT (probability of conserved targeting) — Friedman et al. (2009) ───────────────────

    /// <summary>
    /// Computes the Friedman et al. (2009) <b>branch-length score (Bls)</b> for a seed match: the
    /// total branch length of the minimal subtree of <paramref name="tree"/> that connects the
    /// <paramref name="speciesWithSite"/> (the species in which the site is perfectly aligned),
    /// reported to 4 decimals as TargetScan does.
    /// </summary>
    /// <remarks>
    /// <para>
    /// Per Friedman et al. (2009) Genome Res 19:92 (Methods): conservation "was then assessed by
    /// summing the total branch length in the phylogenetic tree connecting the subset of species
    /// having the sequence perfectly aligned." An edge belongs to that connecting subtree iff at
    /// least one species-with-site lies on each side of the edge (i.e. the edge is on a path
    /// between two such species). A single species (or none) yields Bls = 0. The result does not
    /// depend on where the tree is rooted.
    /// </para>
    /// <para>
    /// TargetScan 7 (<c>targetscan_70_BL_PCT.pl</c>, <c>getBranchLength</c>) computes the same
    /// connecting-subtree length by climbing from each species to the nearest already-included node
    /// and returns <c>sprintf("%.4f", $total_dist)</c>; that rounded value is what it reports and
    /// what it feeds to the PCT sigmoid and the conservation thresholds. This method returns the
    /// same 4-decimal value. TargetScan uses one tree per 3'UTR conservation bin
    /// (<c>PCT_parameters/Tree.bin_01..10.txt</c>); the caller passes the tree for the UTR's bin.
    /// </para>
    /// <para>
    /// Leaf matching uses the leaf <see cref="PhylogeneticAnalyzer.PhyloNode.Name"/> (ordinal); a
    /// species that is not a leaf of the tree is ignored, as in TargetScan ("Any other species will
    /// be ignored during BL calculations", README_BL_PCT_70). Duplicates count once.
    /// </para>
    /// </remarks>
    /// <param name="tree">A phylogenetic tree whose leaf names are species identifiers and whose <c>BranchLength</c> values are edge lengths (e.g. parsed via <see cref="PhylogeneticAnalyzer.ParseNewick"/>).</param>
    /// <param name="speciesWithSite">Species (leaf names) in which the site is conserved.</param>
    /// <returns>The branch-length score (≥ 0), rounded to 4 decimals.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="tree"/> or <paramref name="speciesWithSite"/> is null.</exception>
    public static double ComputeBranchLengthScore(
        PhylogeneticAnalyzer.PhyloNode tree,
        IReadOnlyCollection<string> speciesWithSite)
    {
        ArgumentNullException.ThrowIfNull(tree);
        ArgumentNullException.ThrowIfNull(speciesWithSite);

        return RoundPerl(ConnectingSubtreeLength(tree, speciesWithSite), "F4");
    }

    // Total branch length of the minimal subtree connecting the given leaves (unrounded; species that
    // are not leaves are ignored; fewer than two leaves ⇒ 0). Shared by ComputeBranchLengthScore and
    // the per-column branch lengths of AssignTargetScanUtrBin (targetscan_70_BL_bins.pl).
    private static double ConnectingSubtreeLength(PhylogeneticAnalyzer.PhyloNode tree, IEnumerable<string> speciesWithSite)
    {
        var present = new HashSet<string>(StringComparer.Ordinal);
        foreach (string species in speciesWithSite)
            if (species is not null)
                present.Add(species);
        if (present.Count == 0)
            return 0.0;

        double total = 0.0;
        int found = CountAndSumConnectingEdges(tree, present, ref total, totalPresent: -1);
        if (found < 2)
            return 0.0;
        // Second pass with the number of species that are actually leaves of the tree: an edge is
        // on a connecting path iff its lower side holds ≥ 1 and < found species-with-site.
        total = 0.0;
        CountAndSumConnectingEdges(tree, present, ref total, totalPresent: found);
        return total;
    }

    // Post-order walk returning the number of species-with-site leaves in the subtree of
    // <paramref name="node"/>. When totalPresent ≥ 0, each child edge whose subtree holds between 1
    // and totalPresent − 1 such leaves lies on a path between two of them and is added to total.
    // The root has no parent edge, so its BranchLength is never counted.
    private static int CountAndSumConnectingEdges(
        PhylogeneticAnalyzer.PhyloNode node,
        HashSet<string> present,
        ref double total,
        int totalPresent)
    {
        if (node.IsLeaf)
            return present.Contains(node.Name) ? 1 : 0;

        int count = 0;
        foreach (var child in node.Children)
        {
            int inChild = CountAndSumConnectingEdges(child, present, ref total, totalPresent);
            if (inChild > 0 && inChild < totalPresent)
                total += child.BranchLength;
            count += inChild;
        }
        return count;
    }

    /// <summary>
    /// Maps a branch-length score (Bls) to a PCT value with TargetScan's logistic
    /// (<c>targetscan_70_BL_PCT.pl</c>, <c>calculatePCTthisBL</c> via <c>getPCT</c>):
    /// <c>PCT(Bls) = B0 + B1 / (1 + e^(−B2·Bls + B3))</c>, rounded to 4 decimals and truncated at 0;
    /// Bls = 0 gives PCT 0.
    /// </summary>
    /// <remarks>
    /// <c>getPCT</c> calls the sigmoid only when Bls &gt; 0 and returns <c>"0.0"</c> when Bls = 0
    /// (a site seen in one species has no conserved targeting). <c>calculatePCTthisBL</c> then does
    /// <c>$pct = sprintf("%.4f", $pct); if ($pct &lt; 0) { $pct = "0.0"; }</c>. A rounded −0.0000
    /// is returned as 0. Bls is used as given; <see cref="ComputeBranchLengthScore"/> already
    /// returns TargetScan's 4-decimal value. For a site-type-aware call (6mer ⇒ 0) use
    /// <see cref="PctFromBranchLength(double, PctSigmoidParameters, TargetSiteType)"/>.
    /// </remarks>
    /// <param name="branchLengthScore">The Friedman et al. (2009) branch-length score (finite, ≥ 0).</param>
    /// <param name="p">The b0..b3 row for the seed family and site type (finite values).</param>
    /// <returns>PCT ≥ 0 (in practice ≤ 1), 4 decimals.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="branchLengthScore"/> is negative, NaN or infinite.</exception>
    /// <exception cref="ArgumentException">A sigmoid parameter is NaN or infinite.</exception>
    public static double PctFromBranchLength(double branchLengthScore, PctSigmoidParameters p)
    {
        if (!double.IsFinite(branchLengthScore) || branchLengthScore < 0.0)
            throw new ArgumentOutOfRangeException(nameof(branchLengthScore), branchLengthScore,
                "Branch-length score must be a finite value ≥ 0.");
        if (!double.IsFinite(p.B0) || !double.IsFinite(p.B1) || !double.IsFinite(p.B2) || !double.IsFinite(p.B3))
            throw new ArgumentException("PCT sigmoid parameters must be finite.", nameof(p));

        // getPCT: elsif ($branchLength == 0) { $pct = "0.0"; }
        if (branchLengthScore == 0.0)
            return 0.0;

        // calculatePCTthisBL: $pct = $b0 + ( $b1 / (1 + e ** ( (0 - $b2) * $BL + $b3)));
        double pct = p.B0 + p.B1 / (1.0 + Math.Exp(-p.B2 * branchLengthScore + p.B3));
        // $pct = sprintf ("%.4f", $pct); if ($pct < 0) { $pct = "0.0"; }
        pct = RoundPerl(pct, "F4");
        return pct <= 0.0 ? 0.0 : pct;
    }

    /// <summary>
    /// PCT for a site of the given type, as <c>getPCT</c> assigns it: 0 for a 6mer site (TargetScan
    /// fits no 6mer PCT; <c>($siteType ne "6mer") ? calculatePCTthisBL(…) : "0.0"</c>), otherwise
    /// <see cref="PctFromBranchLength(double, PctSigmoidParameters)"/>.
    /// </summary>
    /// <param name="branchLengthScore">The branch-length score (finite, ≥ 0).</param>
    /// <param name="p">The b0..b3 row for the seed family and <paramref name="siteType"/>.</param>
    /// <param name="siteType">The site type of this species' seed match.</param>
    /// <returns>PCT ≥ 0, 4 decimals.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="branchLengthScore"/> is negative, NaN or infinite.</exception>
    /// <exception cref="ArgumentException">A sigmoid parameter is NaN or infinite.</exception>
    public static double PctFromBranchLength(double branchLengthScore, PctSigmoidParameters p, TargetSiteType siteType)
    {
        double pct = PctFromBranchLength(branchLengthScore, p);
        return siteType == TargetSiteType.Seed6mer ? 0.0 : pct;
    }

    /// <summary>
    /// Looks up the bundled TargetScan 7.0 PCT sigmoid row (b0..b3) for a seed family and site type,
    /// as <c>targetscan_70_BL_PCT.pl</c> does (<c>readAllPCTdata</c> → <c>%familyPlusType2coeff</c>;
    /// <c>calculatePCTthisBL</c>).
    /// </summary>
    /// <remarks>
    /// <para>
    /// TargetScan keys the three parameter files (<c>PCT_parameters/8mer_PCT_parameters.txt</c>,
    /// <c>7mer_m8_PCT_parameters.txt</c>, <c>7mer_1a_PCT_parameters.txt</c>; 111 rows each) by the
    /// 7-nt seed region = miRNA nucleotides 2–8 in the U alphabet (<c>readMiRNAs</c>: column 2 of the
    /// miRNA family file, <c>s/T/U/gi</c>; the family name is mapped to this seed before the lookup).
    /// The key is matched exactly (ordinal, 7 nt). T is read as U as in the perl; lower case is
    /// upper-cased (the perl hash is case-sensitive, and every table key is upper case, so this only
    /// accepts lower-case input that the perl would miss). A miRNA family name (e.g. "let-7/98") is
    /// not a key: pass its seed (see <see cref="GetSeedSequence"/>).
    /// </para>
    /// <para>
    /// Only 8mer, 7mer-m8 and 7mer-A1 have rows; a 6mer (TargetScan PCT is always 0) and any other
    /// site type return <c>false</c>. Use <see cref="IsTargetScanPctSeedFamily"/> to ask whether the
    /// family has PCT at all (TargetScan reports "NA" otherwise). Values are the kit's decimal text
    /// verbatim (<c>MiRnaAnalyzer.TargetScanPctParameters.cs</c>; Friedman et al. 2009 Genome Res 19:92,
    /// Agarwal et al. 2015 eLife 4:e05005 — cite TargetScan when used).
    /// </para>
    /// </remarks>
    /// <param name="seedRegion">The 7-nt seed region (miRNA nt 2–8), RNA or DNA alphabet.</param>
    /// <param name="siteType">8mer, 7mer-m8 or 7mer-A1.</param>
    /// <param name="parameters">The b0..b3 row when found; <c>default</c> otherwise.</param>
    /// <returns><c>true</c> when the bundled table has a row for this seed and site type.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="seedRegion"/> is null.</exception>
    public static bool TryGetTargetScanPctParameters(
        string seedRegion, TargetSiteType siteType, out PctSigmoidParameters parameters)
    {
        ArgumentNullException.ThrowIfNull(seedRegion);
        parameters = default;
        if (!TargetScanPctParameterTable.TryGetValue(TargetScanPctKey(seedRegion), out var row))
            return false;
        switch (siteType)
        {
            case TargetSiteType.Seed8mer: parameters = row.Mer8; return true;
            case TargetSiteType.Seed7merM8: parameters = row.Mer7M8; return true;
            case TargetSiteType.Seed7merA1: parameters = row.Mer7A1; return true;
            default: return false;
        }
    }

    /// <summary>
    /// <c>true</c> when the seed region (miRNA nt 2–8) is one of the 111 seed families in TargetScan
    /// 7.0's PCT parameter tables (<c>targetscan_70_BL_PCT.pl</c> <c>$thisMiRNAhasPCT</c>). For any
    /// other family <c>getPCT</c> reports PCT "NA". Key normalisation as in
    /// <see cref="TryGetTargetScanPctParameters"/>.
    /// </summary>
    /// <param name="seedRegion">The 7-nt seed region, RNA or DNA alphabet.</param>
    /// <exception cref="ArgumentNullException"><paramref name="seedRegion"/> is null.</exception>
    public static bool IsTargetScanPctSeedFamily(string seedRegion)
    {
        ArgumentNullException.ThrowIfNull(seedRegion);
        return TargetScanPctParameterTable.ContainsKey(TargetScanPctKey(seedRegion));
    }

    // readMiRNAs: $mirSeedRegion =~ s/T/U/gi (upper-casing added, see TryGetTargetScanPctParameters).
    private static string TargetScanPctKey(string seedRegion) =>
        seedRegion.ToUpperInvariant().Replace('T', 'U');

    /// <summary>
    /// Looks up the bundled TargetScan 7.0 TA_3UTR and SPS values for a seed region and site type, as
    /// <c>targetscan_70_context_scores.pl</c> reads them from <c>TA_SPS_by_seed_region.txt</c>
    /// ("TA and SPS values (from Garcia et al., Supp Data 5)"; <c>read_TA_SPS</c>, <c>readTargets</c>).
    /// </summary>
    /// <remarks>
    /// <para>
    /// The file has one row for each of the 16 384 seed regions (miRNA nt 2–8, all 7-mers over ACGU):
    /// column 2 = SPS of 8mer and 7mer-m8 sites, column 3 = SPS of 7mer-A1 and 6mer sites
    /// (<c>$garcia{3|2}</c> / <c>$garcia{1|4}</c>), column 4 = TA (log10 site abundance, one value per
    /// seed). <c>readTargets</c> keys it with <c>substr($matureMiRNA, 1, 7)</c> (upper-cased since
    /// TargetScan 7.02). Here T is read as U and lower case is upper-cased; any other symbol (e.g. N)
    /// or a length other than 7 has no row (TargetScan then uses contribution 0).
    /// </para>
    /// <para>
    /// Values are the file's decimal text verbatim (<c>MiRnaAnalyzer.TargetScanTaSps.cs</c>, SHA-256 in
    /// its header; Garcia et al. 2011 Nat Struct Mol Biol 18:1139, Agarwal et al. 2015 eLife 4:e05005 —
    /// cite TargetScan when used). <see cref="ComputeTa3Utr"/> computes TA from your own 3'UTR set
    /// instead (Garcia's definition).
    /// </para>
    /// </remarks>
    /// <param name="seedRegion">The 7-nt seed region (miRNA nt 2–8), RNA or DNA alphabet.</param>
    /// <param name="siteType">8mer, 7mer-m8, 7mer-A1 or 6mer (selects the SPS column).</param>
    /// <param name="ta">TA_3UTR for the seed region when found; 0 otherwise.</param>
    /// <param name="sps">SPS (kcal/mol) for the seed region and site type when found; 0 otherwise.</param>
    /// <returns><c>true</c> when the table has a row for the seed region and the site type is one of the four canonical types.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="seedRegion"/> is null.</exception>
    public static bool TryGetTargetScanTaSps(string seedRegion, TargetSiteType siteType, out double ta, out double sps)
    {
        ArgumentNullException.ThrowIfNull(seedRegion);
        ta = 0.0;
        sps = 0.0;
        int index = TargetScanSeedRegionIndex(seedRegion);
        if (index < 0)
            return false;
        var table = TargetScanTaSpsTable.Value;
        switch (siteType)
        {
            case TargetSiteType.Seed8mer:
            case TargetSiteType.Seed7merM8:
                sps = table.Sps8merAnd7merM8[index];
                break;
            case TargetSiteType.Seed7merA1:
            case TargetSiteType.Seed6mer:
                sps = table.Sps7merA1And6mer[index];
                break;
            default:
                return false;
        }
        ta = table.Ta[index];
        return true;
    }

    // Row index in TA_SPS_by_seed_region.txt = base-4 number of the seed region (A=0, C=1, G=2, U=3;
    // nt 2 most significant): the file lists all 4^7 seeds in A<C<G<U order (checked by the generator).
    // -1 when the seed is not 7 nt over ACGU (T read as U, case-insensitive).
    private static int TargetScanSeedRegionIndex(string seedRegion)
    {
        if (seedRegion.Length != 7)
            return -1;
        int index = 0;
        foreach (char c in seedRegion)
        {
            int digit = char.ToUpperInvariant(c) switch
            {
                'A' => 0,
                'C' => 1,
                'G' => 2,
                'U' or 'T' => 3,
                _ => -1
            };
            if (digit < 0)
                return -1;
            index = index * 4 + digit;
        }
        return index;
    }

    // Parsed once on first use from the generated column texts (exact decimal text → double).
    private static readonly Lazy<(double[] Sps8merAnd7merM8, double[] Sps7merA1And6mer, double[] Ta)> TargetScanTaSpsTable =
        new(() => (ParseTargetScanColumn(TargetScanSps8merAnd7merM8Text),
                   ParseTargetScanColumn(TargetScanSps7merA1And6merText),
                   ParseTargetScanColumn(TargetScanTaText)));

    private static double[] ParseTargetScanColumn(string text)
    {
        string[] fields = text.Split(' ');
        if (fields.Length != TargetScanTaSpsSeedRegionCount)
            throw new InvalidOperationException("Bundled TargetScan TA/SPS column has the wrong row count.");
        var values = new double[fields.Length];
        for (int i = 0; i < fields.Length; i++)
            values[i] = double.Parse(fields[i], System.Globalization.NumberStyles.Float, System.Globalization.CultureInfo.InvariantCulture);
        return values;
    }

    /// <summary>
    /// Site-type-specific branch-length thresholds above which TargetScan calls a site conserved
    /// (<c>targetscan_70_BL_PCT.pl</c>, <c>%siteTypeToConsThreshold</c>, "If the BL is &gt;= these
    /// thresholds, it's conserved (for TS7 trees)"). They belong to the TargetScan 7 trees; 6mer
    /// sites are never conserved (threshold 100 in the script).
    /// </summary>
    /// <param name="Seed8mer">Threshold for 8mer sites.</param>
    /// <param name="Seed7merM8">Threshold for 7mer-m8 sites.</param>
    /// <param name="Seed7merA1">Threshold for 7mer-A1 sites.</param>
    public readonly record struct BranchLengthConservationThresholds(double Seed8mer, double Seed7merM8, double Seed7merA1)
    {
        /// <summary>TargetScanHuman 7: 8mer 1.8, 7mer-m8 2.8, 7mer-1a 3.6 (the script's active values).</summary>
        public static BranchLengthConservationThresholds TargetScanHuman7 => new(1.8, 2.8, 3.6);

        /// <summary>TargetScanMouse 7: 8mer 0.6, 7mer-m8 1.8, 7mer-1a 2.5 (the script's commented mouse block).</summary>
        public static BranchLengthConservationThresholds TargetScanMouse7 => new(0.6, 1.8, 2.5);
    }

    /// <summary>
    /// TargetScan's "Conserved" call for a site (<c>selectConservationForThisSite</c>): conserved
    /// iff Bls ≥ the threshold for its site type; a 6mer is never conserved.
    /// </summary>
    /// <param name="branchLengthScore">The branch-length score from <see cref="ComputeBranchLengthScore"/>.</param>
    /// <param name="siteType">The site type.</param>
    /// <param name="thresholds">Thresholds for the tree in use; <c>null</c> ⇒ <see cref="BranchLengthConservationThresholds.TargetScanHuman7"/>.</param>
    /// <returns><c>true</c> when the site is conserved.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="branchLengthScore"/> is negative, NaN or infinite.</exception>
    public static bool IsConservedByBranchLength(
        double branchLengthScore,
        TargetSiteType siteType,
        BranchLengthConservationThresholds? thresholds = null)
    {
        if (!double.IsFinite(branchLengthScore) || branchLengthScore < 0.0)
            throw new ArgumentOutOfRangeException(nameof(branchLengthScore), branchLengthScore,
                "Branch-length score must be a finite value ≥ 0.");
        var t = thresholds ?? BranchLengthConservationThresholds.TargetScanHuman7;
        return siteType switch
        {
            TargetSiteType.Seed8mer => branchLengthScore >= t.Seed8mer,
            TargetSiteType.Seed7merM8 => branchLengthScore >= t.Seed7merM8,
            TargetSiteType.Seed7merA1 => branchLengthScore >= t.Seed7merA1,
            _ => false
        };
    }

    // PCT enters getAgarwalContribution as a min-max-scaled feature × the site-type PCT coefficient
    // (targetscan_70_context_scores.pl getPCT_contribution → getAgarwalContribution; PCT matches the
    // min-max regex branch alongside TA_3UTR/SPS/Local_AU/3P_score/SA/Len_ORF/Len_3UTR/Min_dist).
    internal static double PctContribution(double pct, TargetSiteType type)
    {
        (double coeff, double min, double max) = type switch
        {
            TargetSiteType.Seed8mer => (CtxPctCoeff8mer, CtxPctMin8mer, CtxPctMax8mer),
            TargetSiteType.Seed7merM8 => (CtxPctCoeff7merM8, CtxPctMin7merM8, CtxPctMax7merM8),
            TargetSiteType.Seed7merA1 => (CtxPctCoeff7merA1, CtxPctMin7merA1, CtxPctMax7merA1),
            _ => (CtxPctCoeff6mer, CtxPctMin6mer, CtxPctMax6mer)
        };
        return ScaledContribution(pct, coeff, min, max);
    }

    #endregion

    #region Pre-miRNA Energy Parameters

    // Turner 2004 parameters live in the single source Seqeron.Genomics.Core.Turner2004Parameters. Neither the
    // miRNA–target duplex energy (CalculateDuplexEnergy) nor the pre-miRNA hairpin energy (AnalyzeHairpin) is
    // computed here: both are delegated to the canonical RnaSecondaryStructure Turner 2004 kernels.

    #endregion

    #region Pre-miRNA Prediction

    // Consecutive-pairing screen thresholds.
    // Stem: pre-miRNA stems carry a ~22-bp miRNA duplex; >= 18 consecutive bp required — Krol et al. (2004).
    private const int HeuristicMinStemBasePairs = 18;
    // Hairpin loop: NNDB Turner 2004 prohibits hairpin loops shorter than 3 nt; pre-miRNA terminal loops are
    // ~3-15 nt, up to ~25 nt — Bartel (2004), Cell 116:281.
    private const int HeuristicMinTerminalLoop = 3;
    private const int HeuristicMaxTerminalLoop = 25;

    /// <summary>
    /// Finds potential pre-miRNA hairpins. DEFAULT (<paramref name="useConsecutiveScreen"/> = <c>false</c>):
    /// every window of length [<paramref name="minHairpinLength"/>, <paramref name="maxHairpinLength"/>] is
    /// folded with the canonical Zuker–Stiegler MFE engine and accepted by the published criteria of
    /// <see cref="AssessHairpinByMfe"/> (single hairpin, ≥ 16 bp — Ambros 2003, loop 3–25 nt — Bartel 2004,
    /// MFEI ≥ <paramref name="minMfei"/> — Zhang 2006), i.e. exactly the windows of
    /// <see cref="FindPreMiRnaHairpinsByMfe"/>. Bulges and internal loops in the stem are therefore tolerated
    /// (e.g. hsa-mir-21 is detected). With <paramref name="useConsecutiveScreen"/> = <c>true</c> the previous
    /// consecutive-pairing screen is used instead: the window's two ends must pair UNINTERRUPTEDLY inward
    /// (Watson–Crick or G·U) for ≥ 18 bp with a 3–25 nt apical loop (the stem is extended as far as the
    /// mirrored bases pair), and the labels are blunt (see below).
    /// </summary>
    /// <remarks>
    /// <para>COST: the default folds every window — O(#windows · W³), ≈ 20 ms per 55–120-nt window on a
    /// desktop core — and the scan is parallelised over window starts (results keep the sequential order).
    /// Use <paramref name="useConsecutiveScreen"/> = <c>true</c> for long sequences when a fast screen is enough.</para>
    /// <para>Mature/star in the default mode: the window is read as a Drosha product (pre-miRNA), so the 5p
    /// product starts at the 5′ base of the outermost base pair of the MFE stem and is
    /// <paramref name="matureLength"/> nt long (Dicer 5′ counting, Park et al. 2011 — the rule of
    /// <see cref="PredictDroshaDicerCleavage(string, int, string?, DicerCountingRule)"/>; capped at the 5′ arm).
    /// The star is the opposite-arm product of the Dicer 2-nt-3′-overhang rule (miRDeep2 <c>find_star</c>,
    /// Friedländer et al. 2008; the same code as <c>PredictDroshaDicerCleavage</c>): it starts at
    /// partner(mature end − 2 − off) − off and ends at partner(mature start) + 2 (clipped to the window).
    /// A window carries no basal-junction knowledge, so no Drosha ruler is applied: mature/star are the
    /// arm-dominant ends of the stem, not an annotated guide/passenger assignment.</para>
    /// </remarks>
    /// <param name="sequence">RNA or DNA sequence (case-insensitive; T read as U). Null/empty ⇒ no candidates.</param>
    /// <param name="minHairpinLength">Minimum window length (default 55 nt; must be ≥ 1).</param>
    /// <param name="maxHairpinLength">Maximum window length (default 120 nt; a value below
    /// <paramref name="minHairpinLength"/> yields no candidates).</param>
    /// <param name="matureLength">Length of the 5p product (default 22 nt; must be ≥ 1; capped at the arm).</param>
    /// <param name="useConsecutiveScreen"><c>true</c> = previous consecutive-pairing screen (fast, no bulges);
    /// <c>false</c> (default) = MFE-fold based detection.</param>
    /// <param name="minMfei">MFEI threshold of the default mode (Zhang 2006, 0.85); ignored by the screen.</param>
    /// <returns>
    /// Default mode: <c>Structure</c>/<c>FreeEnergy</c> = canonical MFE dot-bracket/ΔG° of the window.
    /// Screen mode: <c>Structure</c>/<c>FreeEnergy</c> from the canonical
    /// <see cref="RnaSecondaryStructure.FindStemLoops"/> stem-loop (Turner 2004 stacks + AU/GU end penalties +
    /// full hairpin-loop model), mature = first <paramref name="matureLength"/> nt of the window and star = the
    /// last (blunt base-paired partner; no 2-nt overhang).
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="minHairpinLength"/> or
    /// <paramref name="matureLength"/> is &lt; 1.</exception>
    public static IEnumerable<PreMiRna> FindPreMiRnaHairpins(
        string sequence,
        int minHairpinLength = 55,
        int maxHairpinLength = 120,
        int matureLength = 22,
        bool useConsecutiveScreen = false,
        double minMfei = DefaultMinMfei)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minHairpinLength, 1);
        ArgumentOutOfRangeException.ThrowIfLessThan(matureLength, 1);
        return useConsecutiveScreen
            ? FindPreMiRnaHairpinsIterator(sequence, minHairpinLength, maxHairpinLength, matureLength)
            : FindPreMiRnaHairpinsByMfe(sequence, minHairpinLength, maxHairpinLength, minMfei)
                .Select(h => ToPreMiRna(h, matureLength));
    }

    /// <summary>
    /// Maps an MFE-assessed window onto <see cref="PreMiRna"/>: Structure/FreeEnergy of the fold, mature = 5p
    /// product from the 5′ base of the outermost pair (Dicer 5′ counting), star by miRDeep2 find_star
    /// (2-nt 3′ overhang) — the helpers of <see cref="PredictDroshaDicerCleavage(string, int, string?, DicerCountingRule)"/>.
    /// </summary>
    private static PreMiRna ToPreMiRna(PreMiRnaMfe h, int matureLength)
    {
        string db = h.DotBracket;
        int n = db.Length;
        int[] partner = BuildPairTable(db);
        int matureStart = db.IndexOf('(');                    // outermost pair (single hairpin ⇒ first '(')
        int lastOpen = h.TerminalLoopStart - 1;               // 5' arm ends at the last '('
        int matureEnd = Math.Min(matureStart + matureLength - 1, lastOpen);
        int starEnd = Math.Min(partner[matureStart] + RNaseIII3PrimeOverhang, n - 1);
        if (!TryFindOppositeArmStart(partner, matureStart, matureEnd, out int starStart))
            starStart = partner[matureEnd];                   // unreachable for ≥ 16 bp stems; blunt partner
        starStart = Math.Clamp(starStart, matureEnd + 1, starEnd);

        return new PreMiRna(
            Start: h.Start,
            End: h.End,
            Sequence: h.Sequence,
            MatureSequence: h.Sequence.Substring(matureStart, matureEnd - matureStart + 1),
            StarSequence: h.Sequence.Substring(starStart, starEnd - starStart + 1),
            Structure: db,
            FreeEnergy: h.FreeEnergy);
    }

    private static IEnumerable<PreMiRna> FindPreMiRnaHairpinsIterator(
        string sequence, int minHairpinLength, int maxHairpinLength, int matureLength)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minHairpinLength)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // Scan for potential hairpins
        for (int i = 0; i <= upper.Length - minHairpinLength; i++)
        {
            for (int len = minHairpinLength; len <= Math.Min(maxHairpinLength, upper.Length - i); len++)
            {
                string candidate = upper.Substring(i, len);

                // Check for hairpin structure
                var hairpin = AnalyzeHairpin(candidate, matureLength);
                if (hairpin != null)
                {
                    yield return new PreMiRna(
                        Start: i,
                        End: i + len - 1,
                        Sequence: candidate,
                        MatureSequence: hairpin.Value.Mature,
                        StarSequence: hairpin.Value.Star,
                        Structure: hairpin.Value.Structure,
                        FreeEnergy: hairpin.Value.Energy);
                }
            }
        }
    }

    private static (string Mature, string Star, string Structure, double Energy)? AnalyzeHairpin(string sequence, int matureLength)
    {
        int n = sequence.Length;

        // Uninterrupted pairing of the mirrored ends, inward, as far as the bases pair — but never
        // closing a hairpin loop shorter than the 3-nt NNDB minimum.
        int maxStem = (n - HeuristicMinTerminalLoop) / 2;
        int stemLength = 0;
        while (stemLength < maxStem && RnaSecondaryStructure.CanPair(sequence[stemLength], sequence[n - 1 - stemLength]))
            stemLength++;

        if (stemLength < HeuristicMinStemBasePairs)
            return null;

        int loopSize = n - 2 * stemLength; // >= 3 by construction of maxStem
        if (loopSize < HeuristicMinTerminalLoop || loopSize > HeuristicMaxTerminalLoop)
            return null;

        // Canonical stem-loop for exactly this stem/loop: FindStemLoops with the loop fixed at
        // [stemLength, n − stemLength − 1] extends the helix outward to the window ends, scores it with the
        // full Turner 2004 model and renders the dot-bracket with the canonical renderer.
        var stemLoop = RnaSecondaryStructure.FindStemLoops(
            sequence, minStemLength: stemLength, minLoopSize: loopSize, maxLoopSize: loopSize).Single();

        // 5'-arm "mature" label and the 3'-arm segment it pairs with (the "star" label).
        int armLength = Math.Min(matureLength, stemLength);
        string mature = sequence.Substring(0, armLength);
        string star = sequence.Substring(n - armLength, armLength);

        return (mature, star, stemLoop.DotBracketNotation, stemLoop.TotalFreeEnergy);
    }

    #endregion

    #region Pre-miRNA Prediction — MFE-structure-based (opt-in)

    /// <summary>
    /// A pre-miRNA hairpin candidate assessed from the REAL minimum-free-energy (MFE) secondary
    /// structure produced by the validated Zuker–Stiegler folder
    /// (<see cref="Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure"/>,
    /// RNA-STRUCT-001), rather than the consecutive-pairing heuristic of
    /// <see cref="FindPreMiRnaHairpins"/>.
    /// </summary>
    /// <param name="Start">0-based start of the candidate window in the input sequence.</param>
    /// <param name="End">0-based inclusive end of the candidate window in the input sequence.</param>
    /// <param name="Sequence">The folded candidate (upper-cased; T read as U).</param>
    /// <param name="DotBracket">The folded MFE dot-bracket structure of the candidate.</param>
    /// <param name="StemBasePairs">Number of base pairs in the single dominant hairpin stem.</param>
    /// <param name="TerminalLoopStart">0-based start (within the candidate) of the terminal loop.</param>
    /// <param name="TerminalLoopSize">Length of the single terminal (apical) loop, in nucleotides.</param>
    /// <param name="FreeEnergy">
    /// ΔG° (kcal/mol) of the folded structure — by construction equal to
    /// <see cref="Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMinimumFreeEnergy"/>
    /// for the same candidate (negative for a stable hairpin).
    /// </param>
    /// <param name="Amfe">
    /// Adjusted minimal folding free energy: AMFE = 100·|MFE| / length (kcal/mol per 100 nt).
    /// </param>
    /// <param name="Mfei">
    /// Minimal folding free energy index: MFEI = AMFE / (G+C)% (Zhang et al. 2006).
    /// </param>
    public readonly record struct PreMiRnaMfe(
        int Start,
        int End,
        string Sequence,
        string DotBracket,
        int StemBasePairs,
        int TerminalLoopStart,
        int TerminalLoopSize,
        double FreeEnergy,
        double Amfe,
        double Mfei);

    // Pre-miRNA structural thresholds.
    // Stem: a putative miRNA is embedded in one arm of a fold-back hairpin with >=16 complementary
    // bases to the opposite arm — Ambros et al. (2003), RNA 9:277-279.
    private const int MfeMinStemBasePairs = 16;
    // Terminal (apical) loop: pre-miRNA loops are ~3-15 nt, up to ~25 — Bartel (2004), Cell 116:281.
    private const int MfeMinTerminalLoop = 3;
    private const int MfeMaxTerminalLoop = 25;
    // Per-100-nt normalisation factor for AMFE — Zhang et al. (2006), Cell Mol Life Sci 63:246-254:
    // "AMFE means the MFE of a RNA sequence with 100 nt in length" = MFE/length * 100.
    private const double AmfePer100Nt = 100.0;
    // MFEI discriminative threshold: "MFEI of miRNA precursors is >= 0.85, remarkably higher than
    // other RNAs" — Zhang et al. (2006). Used as the default MFE-fold acceptance cutoff.
    private const double DefaultMinMfei = 0.85;

    /// <summary>
    /// Computes the minimal folding free energy index (MFEI) of Zhang et al. (2006),
    /// Cell Mol Life Sci 63:246-254: MFEI = AMFE / (G+C)% where AMFE = 100·|MFE| / length.
    /// MFE is taken as a magnitude (|ΔG°|) so that the published "MFEI > 0.85" criterion applies
    /// directly to the negative ΔG° returned by the Turner-2004 folder.
    /// </summary>
    /// <param name="freeEnergy">ΔG° (kcal/mol), as returned by the folder (typically negative).</param>
    /// <param name="length">Sequence length in nucleotides (must be &gt; 0).</param>
    /// <param name="gcPercent">G+C content as a percentage in (0,100].</param>
    /// <returns>MFEI (dimensionless). Returns 0 when length or GC% is non-positive.</returns>
    public static double CalculateMfeIndex(double freeEnergy, int length, double gcPercent)
    {
        if (length <= 0 || gcPercent <= 0) return 0.0;
        double amfe = AmfePer100Nt * Math.Abs(freeEnergy) / length;
        return amfe / gcPercent;
    }

    /// <summary>
    /// Finds pre-miRNA hairpins by folding each candidate window with the validated Zuker–Stiegler
    /// MFE engine (<see cref="Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure"/>,
    /// RNA-STRUCT-001) and deriving the hairpin features (single terminal loop, paired-stem count,
    /// ΔG°, AMFE, MFEI) from the ACTUAL MFE dot-bracket structure. This is the scan behind the default of
    /// <see cref="FindPreMiRnaHairpins"/> (which maps its windows onto <see cref="PreMiRna"/> with mature/star
    /// labels); the consecutive-pairing screen is its opt-in alternative (useConsecutiveScreen = true).
    /// </summary>
    /// <remarks>
    /// Acceptance criteria (all from primary sources):
    /// <list type="bullet">
    /// <item>The MFE structure is a SINGLE dominant hairpin: exactly one contiguous terminal loop
    /// (one maximal run of unpaired bases flanked by paired bases), with no branching.</item>
    /// <item>Stem base pairs ≥ <see cref="MfeMinStemBasePairs"/> (16) — Ambros et al. (2003).</item>
    /// <item>Terminal loop ∈ [<see cref="MfeMinTerminalLoop"/>, <see cref="MfeMaxTerminalLoop"/>]
    /// (3–25 nt) — Bartel (2004).</item>
    /// <item>MFEI ≥ <paramref name="minMfei"/> (default 0.85) — Zhang et al. (2006).</item>
    /// </list>
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence (RNA or DNA; T read as U); case-insensitive.</param>
    /// <param name="minHairpinLength">Minimum candidate window length (default 55 nt).</param>
    /// <param name="maxHairpinLength">Maximum candidate window length (default 120 nt).</param>
    /// <param name="minMfei">Minimum MFEI to accept a candidate (default 0.85, Zhang 2006).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    /// <returns>Accepted pre-miRNA hairpins, each with features from the real MFE structure.</returns>
    public static IEnumerable<PreMiRnaMfe> FindPreMiRnaHairpinsByMfe(
        string sequence,
        int minHairpinLength = 55,
        int maxHairpinLength = 120,
        double minMfei = DefaultMinMfei,
        int minLoopSize = 3)
    {
        ArgumentOutOfRangeException.ThrowIfLessThan(minHairpinLength, 1);
        return FindPreMiRnaHairpinsByMfeIterator(sequence, minHairpinLength, maxHairpinLength, minMfei, minLoopSize);
    }

    private static IEnumerable<PreMiRnaMfe> FindPreMiRnaHairpinsByMfeIterator(
        string sequence, int minHairpinLength, int maxHairpinLength, double minMfei, int minLoopSize)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minHairpinLength)
            yield break;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');

        // Window starts are independent folds: scan them in parallel, emit in the sequential (start, length) order.
        var perStart = Enumerable.Range(0, upper.Length - minHairpinLength + 1)
            .AsParallel().AsOrdered()
            .Select(i => ScanMfeWindowsAt(upper, i, minHairpinLength, maxHairpinLength, minMfei, minLoopSize));
        foreach (var windows in perStart)
            foreach (var w in windows)
                yield return w;
    }

    private static List<PreMiRnaMfe> ScanMfeWindowsAt(
        string upper, int i, int minHairpinLength, int maxHairpinLength, double minMfei, int minLoopSize)
    {
        var found = new List<PreMiRnaMfe>();
        for (int len = minHairpinLength; len <= Math.Min(maxHairpinLength, upper.Length - i); len++)
        {
            string candidate = upper.Substring(i, len);
            var assessed = AssessHairpinByMfe(candidate, minMfei, minLoopSize);
            if (assessed != null)
            {
                var a = assessed.Value;
                found.Add(new PreMiRnaMfe(
                    Start: i,
                    End: i + len - 1,
                    Sequence: candidate,
                    DotBracket: a.DotBracket,
                    StemBasePairs: a.StemBasePairs,
                    TerminalLoopStart: a.TerminalLoopStart,
                    TerminalLoopSize: a.TerminalLoopSize,
                    FreeEnergy: a.FreeEnergy,
                    Amfe: a.Amfe,
                    Mfei: a.Mfei));
            }
        }
        return found;
    }

    /// <summary>
    /// Folds a single candidate with the MFE engine and, if its MFE structure is a single dominant
    /// hairpin meeting the pre-miRNA structural thresholds, returns the derived features; otherwise
    /// returns <c>null</c>. Public so callers can assess one known sequence without the window scan.
    /// </summary>
    /// <param name="candidate">Candidate sequence (RNA or DNA; T read as U); case-insensitive.</param>
    /// <param name="minMfei">Minimum MFEI to accept (default 0.85, Zhang 2006).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    public static PreMiRnaMfe? AssessHairpinByMfe(
        string candidate,
        double minMfei = DefaultMinMfei,
        int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(candidate))
            return null;

        var mfe = Seqeron.Genomics.Analysis.RnaSecondaryStructure.CalculateMfeStructure(candidate, minLoopSize);
        string structure = mfe.DotBracket;
        int n = structure.Length;
        if (n == 0)
            return null;

        // A single dominant hairpin: the MFE structure must have exactly one terminal loop
        // (one maximal run of '.' bounded on the left by '(' and on the right by ')'), all pairs
        // nested with no branching. Detect this directly from the dot-bracket.
        if (!TryDescribeSingleHairpin(structure, out int stemBp, out int loopStart, out int loopSize))
            return null;

        if (stemBp < MfeMinStemBasePairs)
            return null;
        if (loopSize < MfeMinTerminalLoop || loopSize > MfeMaxTerminalLoop)
            return null;

        double gcPercent = GcPercent(mfe.Sequence);
        double amfe = AmfePer100Nt * Math.Abs(mfe.FreeEnergy) / n;
        double mfei = CalculateMfeIndex(mfe.FreeEnergy, n, gcPercent);
        if (mfei < minMfei)
            return null;

        return new PreMiRnaMfe(
            Start: 0,
            End: n - 1,
            Sequence: mfe.Sequence,
            DotBracket: structure,
            StemBasePairs: stemBp,
            TerminalLoopStart: loopStart,
            TerminalLoopSize: loopSize,
            FreeEnergy: mfe.FreeEnergy,
            Amfe: amfe,
            Mfei: mfei);
    }

    /// <summary>
    /// Returns true iff <paramref name="dotBracket"/> describes a SINGLE hairpin: a sequence of
    /// '(' (the 5' stem, possibly interrupted by internal/bulge '.'), exactly one apical run of
    /// '.' (the terminal loop), then a sequence of ')' (the 3' stem). No multiloop / branching:
    /// once a ')' is seen no further '(' may appear. Counts stem base pairs and the terminal loop.
    /// </summary>
    internal static bool TryDescribeSingleHairpin(
        string dotBracket, out int stemBasePairs, out int terminalLoopStart, out int terminalLoopSize)
    {
        stemBasePairs = 0;
        terminalLoopStart = -1;
        terminalLoopSize = 0;

        int opens = 0, closes = 0;
        bool seenClose = false;
        // The terminal loop is the run of dots immediately after the last '(' and before the
        // first ')'. Track the index just after the last '(' and the index of the first ')'.
        int lastOpenIndex = -1;
        int firstCloseIndex = -1;

        for (int k = 0; k < dotBracket.Length; k++)
        {
            char c = dotBracket[k];
            if (c == '(')
            {
                if (seenClose) return false; // a '(' after a ')' ⇒ branched / not a single hairpin
                opens++;
                lastOpenIndex = k;
            }
            else if (c == ')')
            {
                if (!seenClose) { seenClose = true; firstCloseIndex = k; }
                closes++;
            }
            else if (c != '.')
            {
                return false; // unexpected symbol
            }
        }

        if (opens == 0 || opens != closes) return false; // need a balanced, non-empty stem
        if (lastOpenIndex < 0 || firstCloseIndex < 0 || firstCloseIndex <= lastOpenIndex) return false;

        // Terminal loop = the dots strictly between the last '(' and the first ')'.
        terminalLoopStart = lastOpenIndex + 1;
        terminalLoopSize = firstCloseIndex - lastOpenIndex - 1;
        // Every position between them must be a dot (a single apical loop, no inner stem).
        for (int k = terminalLoopStart; k < firstCloseIndex; k++)
            if (dotBracket[k] != '.') return false;

        stemBasePairs = opens; // each '(' pairs with one ')'
        return true;
    }

    /// <summary>
    /// G+C content as a percentage in [0,100] — the canonical
    /// <see cref="SequenceExtensions.CalculateGcFractionFast"/> ((G+C)/(A+C+G+T+U), non-nucleotide symbols
    /// excluded from both counts, = Biopython <c>gc_fraction(…, "remove")</c>) scaled by 100.
    /// </summary>
    private static double GcPercent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0.0 : AmfePer100Nt * sequence.CalculateGcFractionFast();

    #endregion

    #region Drosha/Dicer cleavage-site prediction (Han 2006 / Park 2011 measuring rules) — opt-in

    /// <summary>
    /// Which end Dicer measures from when it places the apical (loop-side) cut of a pre-miRNA.
    /// </summary>
    public enum DicerCountingRule
    {
        /// <summary>
        /// 5' counting rule (Park et al. 2011, Nature 475:201): Dicer cleaves ~22 nt from the
        /// 5'-phosphorylated 5' end of the pre-miRNA, so the 5p product is 22 nt long. Default
        /// (pre-miRNAs made by Drosha carry a 5' monophosphate).
        /// </summary>
        FivePrimeCounting = 0,

        /// <summary>
        /// 3' counting rule (Park et al. 2011: Dicer "anchors not only the 3' end but also the 5' end"):
        /// Dicer cleaves ~22 nt from the 2-nt-overhang 3' end, so the 3p product is 22 nt long.
        /// </summary>
        ThreePrimeCounting = 1,
    }

    /// <summary>
    /// Predicted Drosha (basal) and Dicer (apical) cleavage coordinates of a folded pri-miRNA hairpin
    /// and the resulting 5p and 3p products. All coordinates are 0-based indices into
    /// <see cref="Sequence"/>; spans are inclusive.
    /// </summary>
    /// <remarks>
    /// Measuring rules (each constant cited at its declaration):
    /// <list type="bullet">
    /// <item>Drosha cleaves ~11 bp (base PAIRS of the folded stem, not nucleotides) from the basal
    /// ssRNA–dsRNA junction (Han et al. 2006; "an 11 bp basal stem was preferred", Auyeung et al. 2013).</item>
    /// <item>Dicer cleaves ~22 nt from the 5' end (5' counting) or from the 3' end (3' counting) of the
    /// pre-miRNA (Park et al. 2011).</item>
    /// <item>Each RNase III cut leaves a 2-nt 3' overhang; the opposite-arm end is the base-pairing
    /// partner + 2 (miRDeep2 <c>find_star</c>, Friedländer et al. 2008; the same rule Auyeung et al. 2013
    /// used to infer Microprocessor sites from RNAfold structures).</item>
    /// </list>
    /// The 5p product is <see cref="MatureStart"/>..<see cref="MatureEnd"/> and the 3p product is
    /// <see cref="StarStart"/>..<see cref="StarEnd"/>. Which arm is the dominant (guide) miRNA is not
    /// predicted; the names follow the 5p/3p arm, not strand selection.
    /// </remarks>
    /// <param name="BasalJunction">0-based index of the 5'-arm base of the basal junction pair (the first base pair of the stem).</param>
    /// <param name="DroshaCut5Prime">0-based 5' end of the pre-miRNA on the 5' arm = 5' base of the 12th stem pair (11 bp below it are removed).</param>
    /// <param name="DroshaCut3Prime">0-based 3' end of the pre-miRNA on the 3' arm = pairing partner of <see cref="DroshaCut5Prime"/> + 2 (2-nt 3' overhang).</param>
    /// <param name="MatureStart">0-based start of the 5p product (= <see cref="DroshaCut5Prime"/>).</param>
    /// <param name="MatureEnd">0-based inclusive end of the 5p product (Dicer cut on the 5' arm).</param>
    /// <param name="StarStart">0-based start of the 3p product (Dicer cut on the 3' arm).</param>
    /// <param name="StarEnd">0-based inclusive end of the 3p product (= <see cref="DroshaCut3Prime"/>).</param>
    /// <param name="MatureSequence">The predicted 5p product.</param>
    /// <param name="StarSequence">The predicted 3p product.</param>
    /// <param name="Sequence">The upper-cased (T→U) sequence the coordinates index into.</param>
    /// <param name="ThreePrimeOverhang">Length of the RNase III 3' overhang produced at each cut (2).</param>
    /// <param name="HasCnncMotif">True iff a CNNC starts 16–18 nt downstream of the 3p Drosha cut (positions +16..+18, Auyeung et al. 2013).</param>
    public readonly record struct DroshaDicerCleavage(
        int BasalJunction,
        int DroshaCut5Prime,
        int DroshaCut3Prime,
        int MatureStart,
        int MatureEnd,
        int StarStart,
        int StarEnd,
        string MatureSequence,
        string StarSequence,
        string Sequence,
        int ThreePrimeOverhang,
        bool HasCnncMotif)
    {
        /// <summary>Dot-bracket structure the stem was measured on (caller-supplied or the canonical MFE fold).</summary>
        public string Structure { get; init; } = "";

        /// <summary>The Dicer counting rule that placed the apical cut.</summary>
        public DicerCountingRule DicerRule { get; init; }

        /// <summary>True iff U and G sit at positions −14 and −13 upstream of the 5p Drosha cut (basal UG motif, Auyeung et al. 2013).</summary>
        public bool HasBasalUgMotif { get; init; }
    }

    // Drosha basal-junction "ruler": Drosha cleaves ~11 bp (≈ one helical turn) from the basal
    // stem-ssRNA junction — Han J et al. (2006), Cell 125:887-901, doi:10.1016/j.cell.2006.03.043:
    // "The cleavage site is determined mainly by the distance (approximately 11 bp) from the
    // stem-ssRNA junction." Counted in base pairs of the stem — Auyeung et al. (2013) Cell 152:844: "an
    // 11 bp basal stem was preferred" for all four selected pri-miRNAs; bulged nucleotides do not count.
    private const int DroshaCutBpFromBasalJunction = 11;

    // Dicer counting "ruler": ~22 nt from the 5' end (5' counting rule) — Park JE et al. (2011),
    // Nature 475:201-205, doi:10.1038/nature10198 — or from the 3' end (3' counting rule).
    private const int DicerCutNtFrom5PrimeEnd = 22;

    // RNase III staggered cut leaves a 2-nt 3' overhang (both Drosha and Dicer) — Lee Y et al. (2003);
    // Han et al. (2006).
    private const int RNaseIII3PrimeOverhang = 2;

    // Auyeung VC, Ulitsky I, McGeary SE, Bartel DP (2013), Cell 152:844-858. Positions are numbered from
    // the inferred Drosha sites: negative = upstream of the 5p cut, positive = downstream of the 3p cut.
    // CNNC: "enrichment was limited to a small range of positions 16–18 nt downstream of the site".
    private const int CnncMinNtDownstreamOfDroshaCut = 16;
    private const int CnncMaxNtDownstreamOfDroshaCut = 18;
    private const int CnncMotifLength = 4; // C-N-N-C
    // Basal UG: "a U at position −14 … immediately preceded the position [−13] that … displayed a
    // strong primary-sequence preference for a G".
    private const int BasalUgUPosition = 14;

    /// <summary>
    /// Predicts the Drosha and Dicer cleavage sites of a pri-miRNA hairpin with the published measuring
    /// rules applied to its folded stem, using the canonical MFE structure
    /// (<see cref="RnaSecondaryStructure.CalculateMfeStructure"/>) and the 5' counting rule.
    /// </summary>
    /// <param name="sequence">The pri-miRNA nucleotide sequence (RNA or DNA; T read as U; case-insensitive).</param>
    /// <param name="basalJunction">0-based index of the 5'-arm base of the basal junction pair, i.e. the first
    /// base pair of the stem above the basal ssRNA flank.</param>
    /// <returns>The prediction, or <c>null</c> (see <see cref="PredictDroshaDicerCleavage(string, int, string?, DicerCountingRule)"/>).</returns>
    public static DroshaDicerCleavage? PredictDroshaDicerCleavage(string sequence, int basalJunction) =>
        PredictDroshaDicerCleavage(sequence, basalJunction, null, DicerCountingRule.FivePrimeCounting);

    /// <summary>
    /// Predicts the Drosha (basal) and Dicer (apical) cleavage sites and the 5p/3p products of a
    /// pri-miRNA hairpin:
    /// <list type="number">
    /// <item>Drosha: walk up the unbranched stem from the basal junction pair counting BASE PAIRS; the
    /// 5' end of the pre-miRNA is the 5' base of pair 12 (11-bp lower stem, Han 2006); its 3' end is that
    /// base's partner + 2 (2-nt 3' overhang).</item>
    /// <item>Dicer: the product measured from the anchored end is 22 nt (Park 2011, <paramref name="dicerRule"/>);
    /// the opposite product ends are found with the 2-nt-overhang pairing rule (miRDeep2 <c>find_star</c>:
    /// start = partner(end − 2) with bulge compensation, end = partner(start) + 2).</item>
    /// <item>Auyeung (2013) basal UG (−14/−13) and downstream CNNC (+16..+18) flags. The apical UGU/GUG
    /// motif is not flagged: the paper defines it only as enrichment "in this region of the loop"
    /// (Fig. S7F window not given in the text; Fig. 7A says GUG at P22–P24, the mir-30a example is UGUG at
    /// P24–P27, and that UGUG is paired in the RNAfold MFE apical stem), so no reproducible rule exists.</item>
    /// </list>
    /// </summary>
    /// <param name="sequence">The pri-miRNA nucleotide sequence (RNA or DNA; T read as U; case-insensitive).</param>
    /// <param name="basalJunction">0-based index of the 5'-arm base of the basal junction pair.</param>
    /// <param name="dotBracket">Optional secondary structure (e.g. RNAfold output) of the same length; when
    /// null, the canonical MFE structure of <paramref name="sequence"/> is used.</param>
    /// <param name="dicerRule">Dicer counting rule (default 5' counting, Park 2011).</param>
    /// <returns>
    /// The prediction, or <c>null</c> when the sequence is null/empty, the junction is out of range or not
    /// the 5' base of a pair, the unbranched stem above it has fewer than 12 pairs, or a product would fall
    /// outside the sequence or not on opposite arms.
    /// </returns>
    /// <exception cref="ArgumentException"><paramref name="dotBracket"/> has a different length than
    /// <paramref name="sequence"/> or is not balanced.</exception>
    public static DroshaDicerCleavage? PredictDroshaDicerCleavage(
        string sequence, int basalJunction, string? dotBracket, DicerCountingRule dicerRule = DicerCountingRule.FivePrimeCounting)
    {
        if (string.IsNullOrEmpty(sequence))
            return null;
        if (dotBracket is not null && (dotBracket.Length != sequence.Length || !RnaSecondaryStructure.ValidateDotBracket(dotBracket)))
            throw new ArgumentException("Dot-bracket must be balanced and as long as the sequence.", nameof(dotBracket));
        if (basalJunction < 0 || basalJunction >= sequence.Length)
            return null;

        string upper = sequence.ToUpperInvariant().Replace('T', 'U');
        string structure = dotBracket ?? RnaSecondaryStructure.CalculateMfeStructure(upper).DotBracket;
        int n = upper.Length;
        int[] partner = BuildPairTable(structure);

        // Drosha: count base pairs up the unbranched stem from the junction pair (bulges skipped).
        if (partner[basalJunction] <= basalJunction)
            return null;
        var stem5 = WalkUnbranchedStem(partner, basalJunction);
        if (stem5.Count <= DroshaCutBpFromBasalJunction)
            return null;
        int droshaCut5Prime = stem5[DroshaCutBpFromBasalJunction];
        int droshaCut3Prime = partner[droshaCut5Prime] + RNaseIII3PrimeOverhang;
        if (droshaCut3Prime >= n)
            return null;

        int matureStart, matureEnd, starStart, starEnd;
        if (dicerRule == DicerCountingRule.FivePrimeCounting)
        {
            matureStart = droshaCut5Prime;
            matureEnd = matureStart + DicerCutNtFrom5PrimeEnd - 1;
            starEnd = droshaCut3Prime;
            if (matureEnd >= n || !TryFindOppositeArmStart(partner, matureStart, matureEnd, out starStart))
                return null;
        }
        else if (dicerRule == DicerCountingRule.ThreePrimeCounting)
        {
            starEnd = droshaCut3Prime;
            starStart = starEnd - DicerCutNtFrom5PrimeEnd + 1;
            matureStart = droshaCut5Prime;
            if (starStart < 0 || !TryFindOppositeArmEnd(partner, starStart, starEnd, out matureEnd))
                return null;
        }
        else
        {
            throw new ArgumentOutOfRangeException(nameof(dicerRule), dicerRule, "Unknown Dicer counting rule.");
        }

        // Products must sit on opposite arms of the hairpin (5p entirely 5' of the 3p).
        if (matureEnd < matureStart || starEnd < starStart || matureEnd >= starStart)
            return null;

        int u = droshaCut5Prime - BasalUgUPosition;
        bool basalUg = u >= 0 && upper[u] == 'U' && upper[u + 1] == 'G';

        return new DroshaDicerCleavage(
            BasalJunction: basalJunction,
            DroshaCut5Prime: droshaCut5Prime,
            DroshaCut3Prime: droshaCut3Prime,
            MatureStart: matureStart,
            MatureEnd: matureEnd,
            StarStart: starStart,
            StarEnd: starEnd,
            MatureSequence: upper.Substring(matureStart, matureEnd - matureStart + 1),
            StarSequence: upper.Substring(starStart, starEnd - starStart + 1),
            Sequence: upper,
            ThreePrimeOverhang: RNaseIII3PrimeOverhang,
            HasCnncMotif: HasCnncMotifDownstream(upper, droshaCut3Prime))
        {
            Structure = structure,
            DicerRule = dicerRule,
            HasBasalUgMotif = basalUg,
        };
    }

    /// <summary>Pair table from a balanced dot-bracket (canonical parser): partner index, or −1 if unpaired.</summary>
    private static int[] BuildPairTable(string dotBracket)
    {
        var partner = new int[dotBracket.Length];
        Array.Fill(partner, -1);
        foreach (var (i, j) in RnaSecondaryStructure.ParseDotBracket(dotBracket))
        {
            partner[i] = j;
            partner[j] = i;
        }
        return partner;
    }

    /// <summary>
    /// 5'-arm indices of the consecutive pairs of the unbranched helix chain that starts at pair
    /// (<paramref name="i0"/>, partner[i0]), from the base upward. The next pair is the first paired
    /// base 3' of the current 5' base, accepted only when its partner is the last paired base 5' of the
    /// current 3' base (interior loops and bulges are crossed; a multiloop branch ends the chain).
    /// </summary>
    private static List<int> WalkUnbranchedStem(int[] partner, int i0)
    {
        var stem = new List<int> { i0 };
        int i = i0, j = partner[i0];
        while (true)
        {
            int a = i + 1;
            while (a < j && partner[a] < 0) a++;
            if (a >= j) return stem; // hairpin loop
            int b = j - 1;
            while (b > i && partner[b] < 0) b--;
            if (partner[a] != b || b <= a) return stem; // multiloop branch
            stem.Add(a);
            i = a; j = b;
        }
    }

    // miRDeep2 find_star (Friedländer et al. 2008; rajewsky-lab/mirdeep2 src/miRDeep2_core_algorithm.pl):
    // star_beg = partner(mature_end − 2 − off) − off, where off skips unpaired (bulged) bases inward.
    private static bool TryFindOppositeArmStart(int[] partner, int start, int end, out int oppositeStart)
    {
        int anchor = end - RNaseIII3PrimeOverhang;
        int length = anchor - start + 1;
        for (int off = 0; off < length; off++)
        {
            if (partner[anchor - off] >= 0)
            {
                oppositeStart = partner[anchor - off] - off;
                return oppositeStart >= 0 && oppositeStart < partner.Length;
            }
        }
        oppositeStart = -1;
        return false;
    }

    // miRDeep2 find_star: star_end = partner(mature_beg + off) + off + 2.
    private static bool TryFindOppositeArmEnd(int[] partner, int start, int end, out int oppositeEnd)
    {
        int length = end - RNaseIII3PrimeOverhang - start + 1;
        for (int off = 0; off < length; off++)
        {
            if (partner[start + off] >= 0)
            {
                oppositeEnd = partner[start + off] + off + RNaseIII3PrimeOverhang;
                return oppositeEnd >= 0 && oppositeEnd < partner.Length;
            }
        }
        oppositeEnd = -1;
        return false;
    }

    /// <summary>
    /// True iff a CNNC motif (C, any, any, C) starts at positions +16..+18 downstream of the 3p Drosha
    /// cut (position +1 = the first nucleotide after the pre-miRNA 3' end; Auyeung et al. 2013).
    /// </summary>
    private static bool HasCnncMotifDownstream(string seq, int droshaCut3Prime)
    {
        for (int offset = CnncMinNtDownstreamOfDroshaCut; offset <= CnncMaxNtDownstreamOfDroshaCut; offset++)
        {
            int start = droshaCut3Prime + offset;
            if (start + CnncMotifLength > seq.Length)
                break;
            if (seq[start] == 'C' && seq[start + CnncMotifLength - 1] == 'C')
                return true;
        }
        return false;
    }

    #endregion

    #region Pre-miRNA natural-vs-background classifier (trained logistic regression) — opt-in

    /// <summary>
    /// Structure/sequence features of a candidate pre-miRNA hairpin, computed entirely from the
    /// REAL minimum-free-energy (MFE) secondary structure of the validated Zuker–Stiegler folder
    /// (<see cref="RnaSecondaryStructure.CalculateMfeStructure"/>, RNA-STRUCT-001) plus the base
    /// composition. These are the inputs to the trained natural-vs-background logistic-regression
    /// classifier (<see cref="ClassifyPreMiRna"/>). Each feature is published:
    /// <list type="bullet">
    /// <item><b>MFE</b> — minimum free energy ΔG° (kcal/mol) of the MFE fold (miPred "dG" family,
    /// Ng &amp; Mishra 2007, Bioinformatics 23(11):1321-1330). Bonnet et al. (2004) compare this MFE
    /// with the MFEs of dinucleotide shuffles; that randomisation test is
    /// <see cref="CalculateRandfoldPValue"/> (not a model input).</item>
    /// <item><b>AMFE / MFEI</b> — adjusted MFE = 100·|ΔG°|/length and MFEI = AMFE/(G+C)% of Zhang et
    /// al. (2006), Cell Mol Life Sci 63:246-254 (taken as magnitudes, so MFEI ≥ 0).</item>
    /// <item><b>PairedFraction</b> — fraction of NUCLEOTIDES paired in the MFE structure,
    /// 2·|bp|/L ∈ [0,1]; this is twice the normalised base-pairing propensity dP = |bp|/L ∈ [0,0.5] of
    /// Schultes et al. (1999) used by miPred (Ng &amp; Mishra 2007) and microPred (Batuwita &amp;
    /// Palade 2009), which is exposed as <see cref="PreMiRnaFeatures.BasePairingPropensity"/>.</item>
    /// <item><b>GcContent</b> — (G+C)/(A+C+G+U) of the canonical
    /// <see cref="SequenceExtensions.CalculateGcFractionFast"/>; <b>StemBasePairs</b> — |bp| of the MFE
    /// structure (all pairs, not only the terminal stem); <b>LoopSize</b> — longest enclosed unpaired run.</item>
    /// </list>
    /// </summary>
    /// <param name="FreeEnergy">ΔG° of the MFE structure (kcal/mol; negative for a stable hairpin).</param>
    /// <param name="Amfe">Adjusted MFE = 100·|ΔG°|/length (Zhang 2006).</param>
    /// <param name="Mfei">MFEI = AMFE/(G+C)% (Zhang 2006).</param>
    /// <param name="GcContent">G+C content as a fraction in [0,1].</param>
    /// <param name="PairedFraction">Fraction of bases paired in the MFE structure, in [0,1].</param>
    /// <param name="StemBasePairs">Number of base pairs in the MFE structure.</param>
    /// <param name="LoopSize">Largest single run of unpaired bases enclosed by the stem (apical loop).</param>
    /// <param name="Length">Sequence length in nucleotides.</param>
    public readonly record struct PreMiRnaFeatures(
        double FreeEnergy,
        double Amfe,
        double Mfei,
        double GcContent,
        double PairedFraction,
        int StemBasePairs,
        int LoopSize,
        int Length)
    {
        /// <summary>
        /// Normalised base-pairing propensity dP = |bp|/L ∈ [0, 0.5] (Schultes et al. 1999; miPred,
        /// Ng &amp; Mishra 2007): base pairs of the MFE structure per nucleotide. Equals
        /// <see cref="PairedFraction"/>/2. 0 for an empty record.
        /// </summary>
        public double BasePairingPropensity => Length > 0 ? (double)StemBasePairs / Length : 0.0;
    }

    /// <summary>
    /// Result of the trained natural-vs-background pre-miRNA classifier
    /// (<see cref="ClassifyPreMiRna"/>): the extracted features, the logistic probability that the
    /// candidate is a genuine (natural) pre-miRNA hairpin rather than a background/shuffled sequence,
    /// and the resulting boolean call at the chosen decision threshold.
    /// </summary>
    /// <param name="Features">The structure/sequence features fed to the model.</param>
    /// <param name="NaturalProbability">P(natural) ∈ [0,1] from the logistic-regression model.</param>
    /// <param name="IsNatural">True iff <see cref="NaturalProbability"/> ≥ the decision threshold.</param>
    public readonly record struct PreMiRnaClassification(
        PreMiRnaFeatures Features,
        double NaturalProbability,
        bool IsNatural);

    // ── Trained logistic-regression model ─────────────────────────────────────────────────
    // Model: P(natural) = sigmoid(b0 + Σ b_j · z_j), z_j = (x_j − mean_j)/std_j standardised by the
    // TRAINING-set mean and population (ddof 0) standard deviation. Fit: ridge-penalised maximum
    // likelihood, maximising (1/n)·Σ[y log p + (1−y) log(1−p)] − (λ/2)·Σ_j b_j² (λ = 1e-3, intercept
    // not penalised), solved by Newton–Raphson / IRLS — Hastie, Tibshirani & Friedman, "The Elements
    // of Statistical Learning" 2nd ed. §4.4.1 (eqs. 4.21–4.28) with the quadratic penalty of §4.4.4 /
    // §18.3; equivalent to scikit-learn LogisticRegression(C = 1/(λ·n_train)). The whole run is in
    // this file and is re-executed by the unit tests (TrainPreMiRnaClassifier must reproduce the
    // bundled constants), so the model is reproducible from the repository alone.
    //
    // TRAINING DATA (public domain): positives = the 13 human miRBase hairpin precursors in
    // PreMiRnaTrainingPositives (https://mirbase.org, Public Domain Dedication); negatives =
    // DinucleotideShuffle of each positive (Altschul & Erickson 1985 / uShuffle; the background
    // convention of Bonnet et al. 2004 / Xue et al. 2005 / Batuwita & Palade 2009), 4 per positive,
    // drawn from one Random(20060101) in the fixed order positive, its 4 shuffles, next positive …
    // → 65 examples. Split: first ⌊0.7·65⌋ = 45 examples by index = training set, last 20 = held-out.
    // Features: [FreeEnergy, Amfe, Mfei, GcContent, PairedFraction] of ExtractPreMiRnaFeatures, i.e.
    // of the canonical RNA-STRUCT-001 MFE folder AS IT IS NOW. A change of the folder or of the
    // shuffle changes the training features; the reproduction test then fails and the constants
    // below must be regenerated with TrainPreMiRnaClassifier (they are its "R" round-trip output).
    // StemBasePairs / LoopSize / Length are reported but not model inputs.
    private const double PreMiRnaLogRegBias = -4.278229873244336;
    private static readonly double[] PreMiRnaLogRegWeights =
        { -0.9994063824017212, 3.3045163837473233, 3.636520432606299, -0.8781055255888673, -0.260634668873017 };
    private static readonly double[] PreMiRnaFeatureMean =
        { -22.93488888888889, 29.25596144559518, 0.6119707188599345, 0.477134229110014, 0.6329947414221273 };
    private static readonly double[] PreMiRnaFeatureStd =
        { 8.568968205027105, 10.686847007183342, 0.20989563586718338, 0.061249038538481494, 0.10670275339810077 };

    /// <summary>
    /// The 13 human pre-miRNA hairpin precursors (miRBase, public domain) used as the positive class
    /// of the bundled classifier, in training order: hsa-mir-21, hsa-let-7a-1, hsa-let-7b,
    /// hsa-mir-106a, hsa-mir-31, hsa-let-7f-2, hsa-mir-16-1, hsa-mir-24-1, hsa-mir-93, hsa-mir-96,
    /// hsa-mir-98, hsa-mir-147a, hsa-mir-199a-2.
    /// </summary>
    internal static readonly string[] PreMiRnaTrainingPositives =
    {
        "UGUCGGGUAGCUUAUCAGACUGAUGUUGACUGUUGAAUCUCAUGGCAACACCAGUCGAUGGGCUGUCUGACA",
        "UGGGAUGAGGUAGUAGGUUGUAUAGUUUUAGGGUCACACCCACCACUGGGAGAUAACUAUACAAUCUACUGUCUUUCCUA",
        "CGGGGUGAGGUAGUAGGUUGUGUGGUUUCAGGGCAGUGAUGUUGCCCCUCGGAAGAUAACUAUACAACCUACUGCCUUCCCUG",
        "CCUUGGCCAUGUAAAAGUGCUUACAGUGCAGGUAGCUUUUUGAGAUCUACUGCAAUGUAAGCACUUCUUACAUUACCAUGG",
        "GGAGAGGAGGCAAGAUGCUGGCAUAGCUGUUGAACUGGGAACCUGCUAUGCCAACAUAUUGCCAUCUUUCC",
        "UGUGGGAUGAGGUAGUAGAUUGUAUAGUUUUAGGGUCAUACCCCAUCUUGGAGAUAACUAUACAGUCUACUGUCUUUCCCACG",
        "GUCAGCAGUGCCUUAGCAGCACGUAAAUAUUGGCGUUAAGAUUCUAAAAUUAUCUCCAGUAUUAACUGUGCUGCUGAAGUAAGGUUGAC",
        "CUCCGGUGCCUACUGAGCUGAUAUCAGUUCUCAUUUUACACACUGGCUCAGUUCAGCAGGAACAGGAG",
        "CUGGGGGCUCCAAAGUGCUGUUCGUGCAGGUAGUGUGAUUACCCAACCUACUGCUGAGCUAGCACUUCCCGAGCCCCCGG",
        "UGGCCGAUUUUGGCACUAGCACAUUUUUGCUUGUGUCUCUCCGCUCUGAGCAAUCAUGUGCAGUGCCAAUAUGGGAAA",
        "AGGAUUCUGCUCAUGCCAGGGUGAGGUAGUAAGUUGUAUUGUUGUGGGGUAGGGAUAUUAGGCCCCAAUUAGAAGAUAACUAUACAACUUACUACUUUCCCUGGUGUGUGGCAUAUUCA",
        "AAUCUAAAGACAACAUUUCUGCACACACACCAGACUAUGGAAGCCAGUGUGUGGAAAUGCUUCUGCUAGAUU",
        "AGGAAGCUUCUGGAGAUCCUGCUCCGUCGCCCCAGUGUUCAGACUACCUGUUCAGGACAAUGCCGUUGUACAGUAGUCUGCACAUUGGUUAGACUGGGCAAGGGAGAGCA",
    };

    internal const int PreMiRnaTrainingSeed = 20060101;
    internal const int PreMiRnaShufflesPerPositive = 4;
    internal const double PreMiRnaTrainingFraction = 0.70;
    internal const double PreMiRnaRidgeLambda = 1e-3;

    /// <summary>A fitted standardised logistic-regression model (bias, weights, training mean/std).</summary>
    internal sealed record PreMiRnaLogisticModel(double Bias, double[] Weights, double[] Mean, double[] Std);

    /// <summary>The five model inputs of a feature vector, in model order.</summary>
    internal static double[] PreMiRnaModelInputs(PreMiRnaFeatures f) =>
        new[] { f.FreeEnergy, f.Amfe, f.Mfei, f.GcContent, f.PairedFraction };

    /// <summary>
    /// Rebuilds the documented training/held-out data set: each positive of
    /// <see cref="PreMiRnaTrainingPositives"/> (label 1) followed by its
    /// <see cref="PreMiRnaShufflesPerPositive"/> dinucleotide shuffles (label 0), all shuffles drawn
    /// from one <c>Random(PreMiRnaTrainingSeed)</c>. Returns the model inputs, labels and sequences.
    /// </summary>
    internal static (double[][] X, int[] Y, string[] Sequences) BuildPreMiRnaTrainingSet()
    {
        var rng = new Random(PreMiRnaTrainingSeed);
        var x = new List<double[]>();
        var y = new List<int>();
        var seqs = new List<string>();
        void Add(string s, int label)
        {
            x.Add(PreMiRnaModelInputs(ExtractPreMiRnaFeatures(s)!.Value));
            y.Add(label);
            seqs.Add(s);
        }
        foreach (string pos in PreMiRnaTrainingPositives)
        {
            Add(pos, 1);
            for (int k = 0; k < PreMiRnaShufflesPerPositive; k++)
                Add(DinucleotideShuffle(pos, rng), 0);
        }
        return (x.ToArray(), y.ToArray(), seqs.ToArray());
    }

    /// <summary>Number of leading examples of the data set that form the training split (⌊0.7·n⌋).</summary>
    internal static int PreMiRnaTrainingCount(int n) => (int)Math.Floor(n * PreMiRnaTrainingFraction);

    /// <summary>
    /// Re-runs the documented training of the bundled classifier: builds the data set, standardises
    /// the training split by its mean and population standard deviation, and fits the ridge-penalised
    /// logistic regression by Newton–Raphson (<see cref="FitRidgeLogisticRegression"/>).
    /// </summary>
    internal static PreMiRnaLogisticModel TrainPreMiRnaClassifier()
    {
        var (x, y, _) = BuildPreMiRnaTrainingSet();
        int nTrain = PreMiRnaTrainingCount(x.Length);
        int d = x[0].Length;
        var mean = new double[d];
        var std = new double[d];
        for (int j = 0; j < d; j++)
        {
            var col = new double[nTrain];
            for (int i = 0; i < nTrain; i++) col[i] = x[i][j];
            mean[j] = col.Average();
            std[j] = Math.Sqrt(StatisticsHelper.PopulationVariance(col));
        }

        var z = new double[nTrain][];
        for (int i = 0; i < nTrain; i++)
        {
            z[i] = new double[d];
            for (int j = 0; j < d; j++)
                z[i][j] = std[j] == 0 ? 0.0 : (x[i][j] - mean[j]) / std[j];
        }

        var (bias, weights) = FitRidgeLogisticRegression(z, y.Take(nTrain).ToArray(), PreMiRnaRidgeLambda);
        return new PreMiRnaLogisticModel(bias, weights, mean, std);
    }

    /// <summary>
    /// Ridge-penalised logistic regression, maximising
    /// ℓ(b0, b) = (1/n)·Σ_i [y_i log p_i + (1−y_i) log(1−p_i)] − (λ/2)·‖b‖² (intercept not penalised),
    /// by Newton–Raphson / iteratively reweighted least squares (ESL §4.4.1) with step halving when a
    /// full step does not increase ℓ (ESL: "step size halving will guarantee convergence"). Stops when
    /// the largest parameter change is below 1e-13 (or after 200 iterations).
    /// </summary>
    internal static (double Bias, double[] Weights) FitRidgeLogisticRegression(double[][] x, int[] y, double lambda)
    {
        int n = x.Length, d = x[0].Length, m = d + 1;
        var theta = new double[m]; // theta[0] = bias

        double Objective(double[] t)
        {
            double ll = 0;
            for (int i = 0; i < n; i++)
            {
                double eta = t[0];
                for (int j = 0; j < d; j++) eta += t[j + 1] * x[i][j];
                ll -= y[i] == 1 ? Softplus(-eta) : Softplus(eta); // log p or log(1−p)
            }
            double pen = 0;
            for (int j = 1; j < m; j++) pen += t[j] * t[j];
            return ll / n - 0.5 * lambda * pen;
        }

        double current = Objective(theta);
        for (int iter = 0; iter < 200; iter++)
        {
            var grad = new double[m];
            var info = new double[m, m]; // negative Hessian
            for (int i = 0; i < n; i++)
            {
                double eta = theta[0];
                for (int j = 0; j < d; j++) eta += theta[j + 1] * x[i][j];
                double p = Sigmoid(eta), w = p * (1 - p), r = y[i] - p;
                for (int a = 0; a < m; a++)
                {
                    double xa = a == 0 ? 1.0 : x[i][a - 1];
                    grad[a] += r * xa / n;
                    for (int b = 0; b < m; b++)
                        info[a, b] += w * xa * (b == 0 ? 1.0 : x[i][b - 1]) / n;
                }
            }
            for (int j = 1; j < m; j++)
            {
                grad[j] -= lambda * theta[j];
                info[j, j] += lambda;
            }

            double[] step = SolveSymmetricPositiveDefinite(info, grad);
            double scale = 1.0;
            double[] trial = new double[m];
            double next = double.NegativeInfinity;
            for (int h = 0; h < 50; h++)
            {
                for (int a = 0; a < m; a++) trial[a] = theta[a] + scale * step[a];
                next = Objective(trial);
                if (next >= current) break;
                scale *= 0.5;
            }

            double change = 0;
            for (int a = 0; a < m; a++) change = Math.Max(change, Math.Abs(trial[a] - theta[a]));
            if (next < current) break; // no ascent possible: at the optimum to machine precision
            Array.Copy(trial, theta, m);
            current = next;
            if (change < 1e-13) break;
        }

        return (theta[0], theta.Skip(1).ToArray());
    }

    /// <summary>log(1 + e^x) without overflow.</summary>
    private static double Softplus(double x) => x > 0 ? x + Math.Log(1 + Math.Exp(-x)) : Math.Log(1 + Math.Exp(x));

    /// <summary>Solves A·x = b for a symmetric positive-definite A by Cholesky decomposition.</summary>
    private static double[] SolveSymmetricPositiveDefinite(double[,] a, double[] b)
    {
        int m = b.Length;
        var l = new double[m, m];
        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                double sum = a[i, j];
                for (int k = 0; k < j; k++) sum -= l[i, k] * l[j, k];
                l[i, j] = i == j ? Math.Sqrt(sum) : sum / l[j, j];
            }
        }
        var yv = new double[m];
        for (int i = 0; i < m; i++)
        {
            double sum = b[i];
            for (int k = 0; k < i; k++) sum -= l[i, k] * yv[k];
            yv[i] = sum / l[i, i];
        }
        var xv = new double[m];
        for (int i = m - 1; i >= 0; i--)
        {
            double sum = yv[i];
            for (int k = i + 1; k < m; k++) sum -= l[k, i] * xv[k];
            xv[i] = sum / l[i, i];
        }
        return xv;
    }

    // Decision threshold: 0.5 is the Bayes (minimum-error) cutoff of the fitted posterior under 0-1
    // loss (Hastie, Tibshirani & Friedman §2.4/§4.4). The posterior is calibrated to the training
    // class proportions (13 positives : 52 shuffles = 1 : 4), so 0.5 is the minimum-error cutoff for
    // that prior; callers with another prior or cost ratio should pass their own threshold.
    private const double DefaultNaturalThreshold = 0.5;

    /// <summary>
    /// Extracts the published structure/sequence features of a candidate pre-miRNA hairpin from its
    /// REAL minimum-free-energy secondary structure (folded once with the validated Zuker–Stiegler
    /// engine, RNA-STRUCT-001). The features are the inputs to <see cref="ClassifyPreMiRna"/>.
    /// </summary>
    /// <param name="sequence">Candidate nucleotide sequence (RNA or DNA; T read as U; case-insensitive).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    /// <returns>The extracted features, or <c>null</c> when the sequence is null/empty.</returns>
    public static PreMiRnaFeatures? ExtractPreMiRnaFeatures(string sequence, int minLoopSize = 3)
    {
        if (string.IsNullOrEmpty(sequence))
            return null;

        var mfe = RnaSecondaryStructure.CalculateMfeStructure(sequence, minLoopSize);
        string seq = mfe.Sequence;
        int n = seq.Length;
        if (n == 0)
            return null;

        int bp = mfe.BasePairs.Count;
        double pairedFraction = n > 0 ? (2.0 * bp) / n : 0.0; // each pair occupies two bases
        double gcFraction = GcPercent(seq) / AmfePer100Nt;    // GcPercent returns a percentage
        double amfe = AmfePer100Nt * Math.Abs(mfe.FreeEnergy) / n;
        double mfei = CalculateMfeIndex(mfe.FreeEnergy, n, GcPercent(seq));
        int loopSize = LargestEnclosedLoop(mfe.DotBracket);

        return new PreMiRnaFeatures(
            FreeEnergy: mfe.FreeEnergy,
            Amfe: amfe,
            Mfei: mfei,
            GcContent: gcFraction,
            PairedFraction: pairedFraction,
            StemBasePairs: bp,
            LoopSize: loopSize,
            Length: n);
    }

    /// <summary>
    /// Classifies a candidate pre-miRNA hairpin as a genuine (natural) pre-miRNA or a
    /// background/shuffled sequence using the bundled trained logistic-regression model over the
    /// published structure/sequence features (see <see cref="PreMiRnaFeatures"/>). This is an OPT-IN
    /// method — the existing detection methods (<see cref="FindPreMiRnaHairpins"/>,
    /// <see cref="FindPreMiRnaHairpinsByMfe"/>) and their defaults are unchanged. The model was
    /// trained from PUBLIC-DOMAIN miRBase precursors against dinucleotide-shuffled negatives
    /// (Bonnet et al. 2004 convention); NO trained GPL miRDeep2 code or weights are used.
    /// </summary>
    /// <param name="sequence">Candidate nucleotide sequence (RNA or DNA; T read as U; case-insensitive).</param>
    /// <param name="threshold">Decision threshold on P(natural) (default 0.5, the Bayes cutoff).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    /// <returns>The features, P(natural), and the boolean call; <c>null</c> for null/empty input.</returns>
    public static PreMiRnaClassification? ClassifyPreMiRna(
        string sequence,
        double threshold = DefaultNaturalThreshold,
        int minLoopSize = 3)
    {
        var features = ExtractPreMiRnaFeatures(sequence, minLoopSize);
        if (features is null)
            return null;

        double p = ScorePreMiRnaFeatures(features.Value);
        return new PreMiRnaClassification(features.Value, p, p >= threshold);
    }

    /// <summary>
    /// Evaluates the bundled trained logistic-regression model on a feature vector, returning
    /// P(natural) ∈ [0,1]. Exposed so callers can score features they extracted themselves.
    /// </summary>
    public static double ScorePreMiRnaFeatures(PreMiRnaFeatures features)
    {
        double[] raw = PreMiRnaModelInputs(features);

        double z = PreMiRnaLogRegBias;
        for (int j = 0; j < raw.Length; j++)
        {
            double std = PreMiRnaFeatureStd[j] == 0 ? 1.0 : PreMiRnaFeatureStd[j];
            double standardised = (raw[j] - PreMiRnaFeatureMean[j]) / std;
            z += PreMiRnaLogRegWeights[j] * standardised;
        }

        return Sigmoid(z);
    }

    /// <summary>
    /// Result of the randfold randomisation test (<see cref="CalculateRandfoldPValue"/>).
    /// </summary>
    /// <param name="FreeEnergy">ΔG° of the input's MFE fold (kcal/mol).</param>
    /// <param name="PValue">randfold probability (1 + #shuffles with ΔG° ≤ input ΔG°)/(N + 1).</param>
    /// <param name="Iterations">Number N of dinucleotide shuffles folded.</param>
    /// <param name="LowerOrEqualCount">Shuffles whose ΔG° is lower than or equal (|Δ| ≤ 1e-7) to the input's.</param>
    public readonly record struct RandfoldResult(double FreeEnergy, double PValue, int Iterations, int LowerOrEqualCount);

    /// <summary>
    /// The randfold randomisation test of Bonnet, Wuyts, Rouzé &amp; Van de Peer (2004),
    /// Bioinformatics 20(17):2911-2917, dinucleotide mode (<c>randfold -d</c>): the input is folded,
    /// <paramref name="iterations"/> dinucleotide shuffles (<see cref="DinucleotideShuffle"/>) are folded,
    /// and P = (1 + #{shuffle ΔG° ≤ input ΔG°}) / (N + 1), with |Δ| ≤ 1e-7 counted as equal — exactly the
    /// counting of <c>randfold.c</c> (<c>ct_eq = 1</c>; <c>proba = (ct_eq + ct_less)/(rand_num + 1)</c>).
    /// Folding uses the canonical Turner-2004 MFE engine (RNA-STRUCT-001) where randfold links ViennaRNA's
    /// <c>fold()</c>. Pre-miRNAs give small P (lower MFE than their shuffles).
    /// </summary>
    /// <param name="sequence">Candidate sequence (RNA or DNA; T read as U; case-insensitive).</param>
    /// <param name="iterations">Number of shuffles N (randfold's third argument, e.g. 999); 0 gives P = 1.</param>
    /// <param name="random">RNG for the shuffles (fixed seed ⇒ reproducible).</param>
    /// <param name="minLoopSize">Minimum hairpin loop size passed to the folder (NNDB minimum 3).</param>
    /// <returns>The test result, or <c>null</c> for null/empty input.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="random"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="iterations"/> is negative.</exception>
    public static RandfoldResult? CalculateRandfoldPValue(string sequence, int iterations, Random random, int minLoopSize = 3)
    {
        ArgumentNullException.ThrowIfNull(random);
        ArgumentOutOfRangeException.ThrowIfNegative(iterations);
        if (string.IsNullOrEmpty(sequence))
            return null;

        const double Epsilon = 0.0000001; // randfold.c: epsilon for comparing doubles
        double mfe = RnaSecondaryStructure.CalculateMfeStructure(sequence, minLoopSize).FreeEnergy;
        int lowerOrEqual = 0;
        for (int i = 0; i < iterations; i++)
        {
            double shuffled = RnaSecondaryStructure
                .CalculateMfeStructure(DinucleotideShuffle(sequence, random), minLoopSize).FreeEnergy;
            if (Math.Abs(shuffled - mfe) <= Epsilon || shuffled < mfe)
                lowerOrEqual++;
        }

        return new RandfoldResult(mfe, (1.0 + lowerOrEqual) / (iterations + 1), iterations, lowerOrEqual);
    }

    /// <summary>Logistic (sigmoid) link 1/(1+e^-z) — Hastie, Tibshirani &amp; Friedman §4.4.1.</summary>
    private static double Sigmoid(double z) => 1.0 / (1.0 + Math.Exp(-z));

    /// <summary>
    /// Largest single run of unpaired bases ('.') that is enclosed within the structure (i.e. lies
    /// after the first '(' and before the last ')') — the apical loop of a hairpin. Returns 0 when
    /// there is no enclosed unpaired run.
    /// </summary>
    internal static int LargestEnclosedLoop(string dotBracket)
    {
        int firstOpen = dotBracket.IndexOf('(');
        int lastClose = dotBracket.LastIndexOf(')');
        if (firstOpen < 0 || lastClose < 0 || lastClose <= firstOpen)
            return 0;

        int best = 0, run = 0;
        for (int k = firstOpen + 1; k < lastClose; k++)
        {
            if (dotBracket[k] == '.') { run++; if (run > best) best = run; }
            else run = 0;
        }
        return best;
    }

    /// <summary>
    /// Generates a dinucleotide-composition-preserving shuffle of <paramref name="sequence"/>: a
    /// uniformly random sequence among all sequences with the SAME first nucleotide, last nucleotide
    /// and exact dinucleotide (doublet) counts. Algorithm of Altschul &amp; Erickson (1985), Mol Biol
    /// Evol 2:526-538 (random Eulerian trail of the doublet multigraph: random "last-exit" edge per
    /// vertex forming an arborescence into the last vertex, the other edges of each vertex in random
    /// order, then the walk from the first vertex), with the last-exit arborescence drawn by Wilson's
    /// loop-erased random walk as in uShuffle (Jiang, Anderson, Gillespie &amp; Mayne 2008, BMC
    /// Bioinformatics 9:192; <c>ushuffle.c shuffle2</c>; Kandel et al. 1996). Wilson's walk samples the
    /// same distribution as Altschul–Erickson's rejection step (squid <c>StrDPShuffle</c>, used by
    /// randfold) but needs no retry, so every input yields a properly random shuffle.
    /// The random draws are consumed in exactly <c>ushuffle.c</c>'s order (vertices in order of
    /// first occurrence; <c>Next(n)</c> where uShuffle uses <c>random() % n</c>), so a recorded draw
    /// transcript reproduces uShuffle's output character for character. This is the pre-miRNA
    /// background convention of Bonnet et al. (2004) (randfold <c>-d</c>), Xue et al. (2005) and
    /// Batuwita &amp; Palade (2009).
    /// </summary>
    /// <param name="sequence">Input nucleotide sequence (case is normalised to upper, T→U; any other symbol is kept as a vertex of its own).</param>
    /// <param name="random">RNG (supply a fixed seed for reproducibility).</param>
    /// <returns>A shuffled sequence with identical dinucleotide counts; the input (normalised) for length &lt; 3.</returns>
    public static string DinucleotideShuffle(string sequence, Random random)
    {
        ArgumentNullException.ThrowIfNull(random);
        if (string.IsNullOrEmpty(sequence)) return sequence ?? "";

        string s = sequence.ToUpperInvariant().Replace('T', 'U');
        int n = s.Length;
        if (n < 3) return s; // uShuffle: k >= l is an exact copy; length 2 has a single doublet

        // Vertices = distinct symbols in order of first occurrence (uShuffle hash insertion order);
        // out-edge lists in sequence order. The root is the vertex of the last symbol.
        var vertexOf = new Dictionary<char, int>();
        var symbols = new List<char>();
        var vertexSeq = new int[n];
        for (int i = 0; i < n; i++)
        {
            if (!vertexOf.TryGetValue(s[i], out int v))
            {
                v = symbols.Count;
                vertexOf[s[i]] = v;
                symbols.Add(s[i]);
            }
            vertexSeq[i] = v;
        }

        int nv = symbols.Count;
        var edges = new List<int>[nv];
        for (int v = 0; v < nv; v++) edges[v] = new List<int>();
        for (int i = 0; i < n - 1; i++) edges[vertexSeq[i]].Add(vertexSeq[i + 1]);
        int root = vertexSeq[n - 1];

        // Wilson's algorithm: random arborescence of last-exit edges into the root.
        var inTree = new bool[nv];
        var next = new int[nv];
        inTree[root] = true;
        for (int i = 0; i < nv; i++)
        {
            int u = i;
            while (!inTree[u])
            {
                next[u] = random.Next(edges[u].Count);
                u = edges[u][next[u]];
            }
            u = i;
            while (!inTree[u])
            {
                inTree[u] = true;
                u = edges[u][next[u]];
            }
        }

        // Each non-root vertex: move its last-exit edge to the end, permute the others; root: permute all.
        for (int v = 0; v < nv; v++)
        {
            var list = edges[v];
            if (v != root)
            {
                int last = list.Count - 1;
                (list[last], list[next[v]]) = (list[next[v]], list[last]);
                FisherYatesShuffle(list, list.Count - 1, random);
            }
            else
            {
                FisherYatesShuffle(list, list.Count, random);
            }
        }

        // Walk the Eulerian trail from the first vertex.
        var pos = new int[nv];
        var sb = new StringBuilder(n);
        int cur = vertexSeq[0];
        sb.Append(symbols[cur]);
        while (pos[cur] < edges[cur].Count)
        {
            cur = edges[cur][pos[cur]++];
            sb.Append(symbols[cur]);
        }
        return sb.ToString();
    }

    // Fisher–Yates on the first `count` items, Durstenfeld order (i = count-1 … 1, j = Next(i+1)),
    // identical to uShuffle's permutei.
    private static void FisherYatesShuffle<T>(IList<T> list, int count, Random random)
    {
        for (int i = count - 1; i > 0; i--)
        {
            int j = random.Next(i + 1);
            (list[i], list[j]) = (list[j], list[i]);
        }
    }

    #endregion

    #region Target Context Analysis

    /// <summary>
    /// TargetScan <c>$MIN_DIST_TO_CDS</c>: a site whose 1-based 3'UTR start is &lt; 15 lies in the
    /// ribosome-occluded first 15 nt after the stop codon and receives no context score
    /// (<c>targetscan_70_context_scores.pl</c> "too_close"; Grimson et al. 2007).
    /// </summary>
    private const int MinDistToCds = 15;

    /// <summary>
    /// Annotates the 3'UTR context of a target site with the site-context determinants of
    /// Grimson et al. (2007, Mol Cell 27:91): AU-rich flanks, placement ≥ 15 nt downstream of the
    /// stop codon, and placement near the ends of the UTR.
    /// <paramref name="mRnaSequence"/> is interpreted as the 3'UTR (index 0 = first nt after the stop codon).
    /// </summary>
    /// <remarks>
    /// <list type="bullet">
    /// <item><b>AuContent</b> — TargetScan 5.0 position-weighted local A/U fraction of up to
    /// <paramref name="contextWindow"/> nt (TargetScan: 30) on each side of the site, the site itself
    /// excluded; weights are those of <paramref name="siteType"/> (see
    /// <see cref="ScoreTargetSiteGrimson(string, string, int, int, TargetSiteType)"/>). 0 when the site has no
    /// flanking nucleotide or <paramref name="siteType"/> has no TargetScan 5.0 weights.</item>
    /// <item><b>NearStart</b> — the site starts within the first 15 nt of the 3'UTR (1-based start
    /// &lt; 15, TargetScan <c>$MIN_DIST_TO_CDS</c>); such sites are cleared by the translating
    /// ribosome and are not effective (Grimson 2007).</item>
    /// <item><b>NearEnd</b> — descriptive flag: the site ends in the 3'-terminal 15% of the sequence.</item>
    /// <item><b>ContextScore</b> (log2 fold-change units, negative = repressed) — the UTR-context part of
    /// the TargetScan 5.0 context score: <c>local-AU contribution + position contribution</c>, i.e.
    /// round3(slope_AU·AuContent + intercept_AU) + round3(slope_pos·min(distance to the nearest UTR end, 1500)
    /// + intercept_pos) with the Grimson 2007 regression constants of <paramref name="siteType"/>. The
    /// site-type and 3'-pairing terms need the miRNA; the complete score is
    /// <see cref="ScoreTargetSiteGrimson(string, string, int, int, TargetSiteType)"/>. 0 when
    /// <b>NearStart</b> or when <paramref name="siteType"/> is not 8mer / 7mer-m8 / 7mer-A1.</item>
    /// </list>
    /// Empty input or an out-of-range site ([start, end] not within the sequence, or start &gt; end)
    /// yields all zeros / false. The coordinates are those of <see cref="TargetSite"/> (an 8mer includes its A1).
    /// </remarks>
    public static (double AuContent, bool NearStart, bool NearEnd, double ContextScore) AnalyzeTargetContext(
        string mRnaSequence,
        int targetStart,
        int targetEnd,
        int contextWindow = 30,
        TargetSiteType siteType = TargetSiteType.Seed7merM8)
    {
        if (string.IsNullOrEmpty(mRnaSequence) || targetStart < 0 || targetEnd >= mRnaSequence.Length || targetStart > targetEnd)
            return (0, false, false, 0);

        string mrna = NormalizeSeed(mRnaSequence);
        int n = mrna.Length;
        int perlType = GrimsonPerlType(siteType);

        double auContent = perlType == 0
            ? 0.0
            : GrimsonLocalAuFraction(mrna, targetStart, targetEnd, perlType, Math.Max(0, contextWindow));

        bool nearStart = targetStart + 1 < MinDistToCds;
        bool nearEnd = targetEnd > n * 0.85;

        if (nearStart || perlType == 0)
            return (auContent, nearStart, nearEnd, 0.0);

        double au = RoundPerl3(auContent * GrimsonAuSlope[perlType] + GrimsonAuIntercept[perlType]);
        double position = GrimsonPositionContribution(n, targetStart, targetEnd, perlType);
        return (auContent, false, nearEnd, au + position);
    }

    /// <summary>
    /// Structural accessibility of a target site: the equilibrium probability that every
    /// nucleotide of [<paramref name="siteStart"/>..<paramref name="siteEnd"/>] (0-based, inclusive)
    /// is unpaired, P = Z_open/Z from the Turner-2004 McCaskill partition function (canonical
    /// <see cref="RnaSecondaryStructure.CalculateRegionUnpairedProbability"/>), folded over a local
    /// RNAplfold-style context of W = max(80, site length) nt centred on the site (TargetScan runs
    /// <c>RNAplfold -W 80</c>; Bernhart et al. 2006).
    /// 1 = fully accessible, 0 = always paired.
    /// </summary>
    /// <remarks>
    /// Returns 0 for empty input or an out-of-range site (start &lt; 0, end ≥ length, start &gt; end).
    /// DNA input (T) and lower case are accepted. RNAplfold additionally limits base-pair span to
    /// L = 40 and averages over all W-windows covering the site; here the single centred window is
    /// folded without a span limit (same approach as the context++ SA feature).
    /// </remarks>
    public static double CalculateSiteAccessibility(string mRnaSequence, int siteStart, int siteEnd)
    {
        if (string.IsNullOrEmpty(mRnaSequence) || siteStart < 0 || siteEnd >= mRnaSequence.Length || siteStart > siteEnd)
            return 0;

        string mrna = mRnaSequence.ToUpperInvariant().Replace('T', 'U');
        return LocalRegionUnpairedProbability(mrna, siteStart, siteEnd);
    }

    #endregion

    #region miRNA Family Analysis

    /// <summary>
    /// Groups miRNAs by their seed sequence family.
    /// </summary>
    /// <remarks>
    /// Family = identical sequence at nucleotides 2–8 (Bartel 2009; TargetScan
    /// <c>miR_Family_Info</c> "Seed+m8"). The family key is the TargetScan-normalised seed
    /// (upper-case, T→U), so the grouping agrees with <see cref="CompareSeedRegions"/>.
    /// miRNAs without a defined seed (sequence shorter than 8 nt ⇒ empty seed) belong to no
    /// family and are omitted.
    /// </remarks>
    public static IEnumerable<(string SeedFamily, IReadOnlyList<MiRna> Members)> GroupBySeedFamily(IEnumerable<MiRna> miRnas)
    {
        ArgumentNullException.ThrowIfNull(miRnas);

        return miRnas
            .Select(m => (Key: NormalizeSeed(m.SeedSequence), MiRna: m))
            .Where(x => x.Key.Length > 0)
            .GroupBy(x => x.Key, x => x.MiRna)
            .Select(g => (g.Key, (IReadOnlyList<MiRna>)g.ToList()));
    }

    /// <summary>
    /// Finds miRNAs whose seed (nt 2–8) is within <paramref name="maxMismatches"/> of the query seed.
    /// </summary>
    /// <remarks>
    /// Distance is the seed mismatch count of <see cref="CompareSeedRegions"/> (Hamming distance
    /// plus length difference, after TargetScan normalisation). Entries with the same
    /// <see cref="MiRna.Name"/> as the query are skipped. A miRNA with an undefined (empty) seed
    /// is never similar to anything: if the query seed is empty, nothing is returned.
    /// </remarks>
    public static IEnumerable<MiRna> FindSimilarMiRnas(MiRna query, IEnumerable<MiRna> database, int maxMismatches = 1)
    {
        ArgumentNullException.ThrowIfNull(database);
        return FindSimilarMiRnasIterator(query, database, maxMismatches);
    }

    private static IEnumerable<MiRna> FindSimilarMiRnasIterator(MiRna query, IEnumerable<MiRna> database, int maxMismatches)
    {
        if (NormalizeSeed(query.SeedSequence).Length == 0)
            yield break;

        foreach (var mirna in database)
        {
            if (mirna.Name == query.Name)
                continue;

            if (NormalizeSeed(mirna.SeedSequence).Length == 0)
                continue;

            if (CompareSeedRegions(query, mirna).Mismatches <= maxMismatches)
                yield return mirna;
        }
    }

    #endregion

    #region Utility Methods

    /// <summary>
    /// Calculates the GC content of a sequence as a fraction in [0, 1].
    /// Delegates to the canonical <see cref="SequenceExtensions.CalculateGcFractionFast"/>
    /// ((G+C)/(A+C+G+T+U); other characters excluded from both counts).
    /// </summary>
    public static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcFractionFast();

    /// <summary>
    /// Enumerates the seed itself followed by every single-nucleotide substitution over
    /// A/C/G/U (1 + 3·L sequences for a length-L seed, in position-then-ACGU order).
    /// </summary>
    /// <remarks>
    /// The input is TargetScan-normalised first (upper-case, T→U) so that no variant duplicates
    /// the original. <paramref name="includeWobble"/> has no effect: every single substitution,
    /// including those whose target-site pairing would be a G:U wobble, is always enumerated.
    /// It is retained only for API compatibility.
    /// </remarks>
    /// <exception cref="ArgumentNullException"><paramref name="seedSequence"/> is null.</exception>
    public static IEnumerable<string> GenerateSeedVariants(string seedSequence, bool includeWobble = true)
    {
        ArgumentNullException.ThrowIfNull(seedSequence);
        return GenerateSeedVariantsIterator(NormalizeSeed(seedSequence));
    }

    private static IEnumerable<string> GenerateSeedVariantsIterator(string seedSequence)
    {
        yield return seedSequence;

        // Generate single-nucleotide variants at each position (all four bases)
        const string bases = "ACGU";

        for (int i = 0; i < seedSequence.Length; i++)
        {
            foreach (char b in bases)
            {
                if (b != seedSequence[i])
                {
                    char[] variant = seedSequence.ToCharArray();
                    variant[i] = b;
                    yield return new string(variant);
                }
            }
        }
    }

    #endregion
}
