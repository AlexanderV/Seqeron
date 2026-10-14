using System.Text;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Provides algorithms for codon optimization and sequence design for heterologous expression.
/// </summary>
public static class CodonOptimizer
{
    #region Records and Types

    /// <summary>
    /// Represents an organism's codon usage table.
    /// </summary>
    /// <param name="OrganismName">Organism name.</param>
    /// <param name="CodonFrequencies">
    /// Codon usage (RNA keys). For the built-in presets: relative frequency of each codon within
    /// its synonymous family (Kazusa "fraction" column).
    /// </param>
    /// <param name="CodonToAminoAcid">Codon → one-letter amino acid mapping.</param>
    public readonly record struct CodonUsageTable(
        string OrganismName,
        IReadOnlyDictionary<string, double> CodonFrequencies,
        IReadOnlyDictionary<string, string> CodonToAminoAcid)
    {
        /// <summary>
        /// Optional usage of each codon per thousand codons of the SAME table entry (Kazusa / CUTG
        /// "/1000" column, RNA keys), or <see langword="null"/> when the table carries only
        /// <see cref="CodonFrequencies"/>. Unlike the per-amino-acid fractions, these values weight
        /// each residue by its overall usage — the input the %MinMax reference implementation takes
        /// (Clarke &amp; Clark 2008; Clark-lab <c>calculateMinMax</c>), so
        /// <see cref="CalculateMinMaxProfile(string, CodonUsageTable, int)"/> uses them when present.
        /// Set for <see cref="Yeast"/> (Kazusa species=4932) and <see cref="Human"/> (Kazusa
        /// species=9606); <see cref="EColiK12"/> has none (see its remarks).
        /// </summary>
        public IReadOnlyDictionary<string, double>? PerThousand { get; init; }
    }

    /// <summary>
    /// Which codon-usage values <see cref="CalculateMinMaxProfile(string, CodonUsageTable, int, MinMaxUsageBasis)"/>
    /// feeds into the %MinMax sums.
    /// </summary>
    public enum MinMaxUsageBasis
    {
        /// <summary>
        /// <see cref="CodonUsageTable.PerThousand"/> when the table carries it, otherwise
        /// <see cref="CodonUsageTable.CodonFrequencies"/> (default; per-thousand usage is the input of
        /// the Clark-lab reference implementation).
        /// </summary>
        PerThousandWhenAvailable = 0,

        /// <summary>
        /// Always <see cref="CodonUsageTable.CodonFrequencies"/> (for the presets: per-amino-acid
        /// fractions, every residue weighted equally — the pre-2026-10 behaviour for Yeast/Human).
        /// </summary>
        CodonFrequencies = 1,
    }

    /// <summary>
    /// Result of codon optimization.
    /// </summary>
    public readonly record struct OptimizationResult(
        string OriginalSequence,
        string OptimizedSequence,
        string ProteinSequence,
        double OriginalCAI,
        double OptimizedCAI,
        double GcContentOriginal,
        double GcContentOptimized,
        int ChangedCodons,
        IReadOnlyList<(int Position, string Original, string Optimized)> Changes);

    /// <summary>
    /// One window of the %MinMax codon-usage profile of Clarke &amp; Clark (2008).
    /// </summary>
    /// <param name="WindowStartCodon">0-based codon index of the first codon in the window.</param>
    /// <param name="PercentMinMax">
    /// Signed %MinMax value for the window: positive values are %Max (a run of predominantly
    /// common codons), negative values are %Min (a run of predominantly rare codons), and 0
    /// means the window's codon usage equals the per-amino-acid average. A value of −100
    /// is a window encoded entirely with the rarest synonymous codons; +100 is one encoded
    /// entirely with the most common synonymous codons. (Clarke &amp; Clark 2008.)
    /// </param>
    public readonly record struct MinMaxWindow(int WindowStartCodon, double PercentMinMax);

    /// <summary>
    /// One rare-codon cluster (RCC) detected by the Sherlocc rule of Chartier et&#160;al. (2012).
    /// </summary>
    /// <param name="StartCodon">0-based codon index of the first codon in the cluster window.</param>
    /// <param name="EndCodon">0-based codon index of the last codon in the cluster window (inclusive).</param>
    /// <param name="RareCount">Number of rare ("pause") codons inside the window.</param>
    public readonly record struct RareCodonCluster(int StartCodon, int EndCodon, int RareCount);

    /// <summary>
    /// Parameters of the Sherlocc 1.0 stage-3 analysis
    /// (<see cref="FindConservedRareCodonClusters(IReadOnlyList{string}, CodonUsageTable, double, SherloccOptions?)"/>).
    /// Defaults are the constants of <c>sherlocc.pl</c> (github.com/mtthchrtr/sherlocc, commit b3b4d15).
    /// </summary>
    public sealed record SherloccOptions
    {
        /// <summary>Width (odd, in alignment columns) of the centred averaging window: columns i−3 … i+3 (l. 785, default 7).</summary>
        public int AveragingWindow { get; init; } = 7;

        /// <summary>Width (in alignment columns) of the cluster window i … i+6 (l. 838–842, default 7).</summary>
        public int ClusterWindow { get; init; } = 7;

        /// <summary>Minimum pause columns inside a cluster window (l. 848 <c>$win_pause_count &gt;= 4</c>, default 4).</summary>
        public int MinPauses { get; init; } = 4;

        /// <summary>
        /// Representativeness ("conservation") filter: a column is averaged only when at least
        /// <c>MinCoverage × sequences × AveragingWindow</c> window cells carry a usage value
        /// (l. 772 <c>$tresh=($row*(7))*(0.75)</c>, l. 799 <c>$freq_count &gt;= $tresh</c>; default 0.75 —
        /// README "windows occupied by alignment gaps … by 25% or more are excluded").
        /// </summary>
        public double MinCoverage { get; init; } = 0.75;

        /// <summary>
        /// <see langword="false"/> (default): usage values are <see cref="CodonUsageTable.PerThousand"/>
        /// (Sherlocc's CUTG per-thousand column, <c>files/cut</c>); a table without it is rejected.
        /// <see langword="true"/>: use <see cref="CodonUsageTable.CodonFrequencies"/> instead (the
        /// threshold must then be on that scale).
        /// </summary>
        public bool UseCodonFrequencies { get; init; }
    }

    /// <summary>
    /// One conserved rare-codon cluster of a codon alignment (Sherlocc 1.0 stage 3).
    /// </summary>
    /// <param name="StartColumn">0-based alignment column of the first pause of the cluster (Sherlocc prints it 1-based).</param>
    /// <param name="EndColumn">0-based alignment column of the last pause of the cluster (inclusive).</param>
    /// <param name="PauseCount">Number of pause columns in StartColumn…EndColumn.</param>
    /// <param name="MeanUsage">
    /// Sherlocc "Usage Frequency Average" (<c>$cluster_avrg</c>, l. 875–879): mean of the column
    /// averages over the span; a non-representative column inside the span has no average and
    /// contributes 0 (Perl undefined value), exactly as in the script.
    /// </param>
    /// <param name="MiddlePoint">
    /// Sherlocc "Middle Point" (l. 1137–1140), 1-based: <c>stop − (length / 2)</c> with the decimals dropped.
    /// </param>
    public readonly record struct ConservedRareCodonCluster(
        int StartColumn, int EndColumn, int PauseCount, double MeanUsage, int MiddlePoint);

    /// <summary>
    /// Result of <see cref="FindConservedRareCodonClusters(IReadOnlyList{string}, CodonUsageTable, double, SherloccOptions?)"/>.
    /// </summary>
    /// <param name="SequenceCount">Number of aligned sequences (rows).</param>
    /// <param name="AlignmentLength">Number of codon columns.</param>
    /// <param name="ColumnAverages">Per column: window-and-family average usage (<c>$avrg</c>), or <see langword="null"/> when the column is not representative ("NA").</param>
    /// <param name="ColumnUsageCounts">Per column: number of window cells with a usage value (<c>$freq_count</c>).</param>
    /// <param name="PauseColumns">0-based pause ("slow") columns: average ≤ threshold.</param>
    /// <param name="Clusters">Merged conserved rare-codon clusters in column order.</param>
    /// <param name="TotalClusterLength">Σ (EndColumn − StartColumn + 1) over the clusters.</param>
    /// <param name="ClusterFraction">TotalClusterLength / AlignmentLength (Sherlocc "Fraction of the pfam occupied by rare codon clusters").</param>
    /// <param name="MeanClusterUsage">Mean of the clusters' <see cref="ConservedRareCodonCluster.MeanUsage"/> (<c>clusters_freq_avrg</c>), or <see langword="null"/> without clusters.</param>
    public sealed record SherloccResult(
        int SequenceCount,
        int AlignmentLength,
        IReadOnlyList<double?> ColumnAverages,
        IReadOnlyList<int> ColumnUsageCounts,
        IReadOnlyList<int> PauseColumns,
        IReadOnlyList<ConservedRareCodonCluster> Clusters,
        int TotalClusterLength,
        double ClusterFraction,
        double? MeanClusterUsage);

    /// <summary>
    /// Optimization strategy options.
    /// </summary>
    public enum OptimizationStrategy
    {
        MaximizeCAI,           // Use most frequent codons
        BalancedOptimization,  // Balance CAI with other factors
        HarmonizeExpression,   // Match host codon usage distribution
        MinimizeSecondary,     // Avoid mRNA secondary structures
        AvoidRareCodeons       // Only replace rare codons
    }

    #endregion

    #region Genetic-code families

    // Codon → amino acid and the synonymous-codon families of an NCBI genetic code, RNA
    // spelling, in NCBI codon order (UUU, UUC, UUA, … GGG) — derived from the canonical
    // GeneticCode tables, never from a private copy (DUP_MAP §2). The NCBI order fixes the
    // tie-break of every "most frequent synonymous codon" choice below: when two synonymous
    // codons share the maximal table frequency, the one that comes first in NCBI order wins, so
    // optimisation is deterministic (Biopython CodonAdaptationIndex.optimize only warns on such
    // ties and picks one arbitrarily).
    //
    // A family is the set of codons GeneticCode.CodonTable assigns to one residue — DNA Chisel's
    // get_backtranslation_table(genetic_table) (Biopython CodonTable.forward_table grouped by
    // amino acid) and the grouping of the canonical CodonUsageAnalyzer core (RSCU/CAI/ENC). The
    // context-dependent ("dual-coding") stops of NCBI tables 27/28/31 therefore sit in the family
    // of the amino acid they encode, exactly as in Biopython/DNA Chisel and CodonUsageAnalyzer.
    private sealed class CodeFamilies
    {
        public CodeFamilies(GeneticCode code)
        {
            Code = code;
            CodonToAminoAcid = code.CodonTable.ToDictionary(kv => kv.Key, kv => kv.Value.ToString());
            Synonyms = code.CodonTable
                .GroupBy(kv => kv.Value)
                .ToDictionary(g => g.Key, g => g.Select(kv => kv.Key).ToArray());
        }

        public GeneticCode Code { get; }

        public IReadOnlyDictionary<string, string> CodonToAminoAcid { get; }

        public IReadOnlyDictionary<char, string[]> Synonyms { get; }
    }

    // One immutable family set per NCBI table (GeneticCode instances are per-table singletons).
    private static readonly IReadOnlyDictionary<int, CodeFamilies> FamiliesByTable =
        GeneticCode.SupportedTableNumbers.ToDictionary(n => n, n => new CodeFamilies(GeneticCode.GetByTableNumber(n)));

    private static readonly CodeFamilies StandardFamilies = FamiliesByTable[GeneticCode.Standard.TableNumber];

    private static readonly IReadOnlyDictionary<string, string> StandardCodonToAminoAcid =
        StandardFamilies.CodonToAminoAcid;

    // null → the NCBI Standard code (table 1), the historical behaviour of every API here.
    private static CodeFamilies FamiliesFor(GeneticCode? geneticCode) =>
        geneticCode is null ? StandardFamilies : FamiliesByTable[geneticCode.TableNumber];

    /// <summary>
    /// Translates one RNA triplet with the canonical genetic-code table of <paramref name="families"/>.
    /// Returns <see langword="false"/> (with <paramref name="aminoAcid"/> = 'X') for a triplet
    /// that is not a valid IUPAC codon, so callers can leave such a triplet untouched instead of
    /// throwing. IUPAC-ambiguous triplets translate per Biopython (GCN → A, UAR → *, NNN → X).
    /// </summary>
    private static bool TryTranslateCodon(string codon, CodeFamilies families, out char aminoAcid)
    {
        aminoAcid = 'X';
        if (codon.Length != 3)
            return false;

        foreach (char c in codon)
        {
            char dna = ToUpperAscii(c);
            if (dna == 'U') dna = 'T';
            if (!IupacHelper.IsNucleotideCode(dna))
                return false;
        }

        aminoAcid = families.Code.Translate(codon);
        return true;
    }

    // Synonymous codons of the amino acid encoded by <paramref name="codon"/> under the code of
    // <paramref name="families"/>, or an empty span when the triplet is unusable (non-IUPAC
    // symbols) or its amino acid is itself ambiguous (B, Z, J, X — no unique synonymous family).
    private static string[] SynonymsFor(string codon, CodeFamilies families)
    {
        if (!TryTranslateCodon(codon, families, out char aa))
            return Array.Empty<string>();
        return families.Synonyms.TryGetValue(aa, out var syn) ? syn : Array.Empty<string>();
    }

    #endregion

    #region Predefined Codon Usage Tables

    /// <summary>
    /// E. coli K12 codon usage frequencies (relative fraction per amino acid).
    /// Source: Kazusa Codon Usage Database, species=316407 (E. coli K-12 substr. W3110, 4332 CDS).
    /// URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407
    /// </summary>
    /// <remarks>
    /// <see cref="CodonUsageTable.PerThousand"/> is <see langword="null"/>: the per-thousand column
    /// of this Kazusa entry could not be obtained from a primary source (review 2026-10), and the
    /// EMBOSS CUTG146 <c>Eecoli.cut</c> / Clark-lab <c>EcolCUB.txt</c> tables are different E. coli
    /// entries whose fractions disagree with these at 9 / 43 codons. %MinMax on this preset therefore
    /// uses the per-amino-acid fractions (equal residue weighting).
    /// </remarks>
    public static readonly CodonUsageTable EColiK12 = new(
        "Escherichia coli K12",
        new Dictionary<string, double>
        {
            // Phenylalanine (F)
            { "UUU", 0.57 }, { "UUC", 0.43 },
            // Leucine (L)
            { "UUA", 0.13 }, { "UUG", 0.13 }, { "CUU", 0.10 }, { "CUC", 0.10 }, { "CUA", 0.04 }, { "CUG", 0.50 },
            // Isoleucine (I)
            { "AUU", 0.51 }, { "AUC", 0.42 }, { "AUA", 0.07 },
            // Methionine (M)
            { "AUG", 1.00 },
            // Valine (V)
            { "GUU", 0.26 }, { "GUC", 0.22 }, { "GUA", 0.15 }, { "GUG", 0.37 },
            // Serine (S)
            { "UCU", 0.15 }, { "UCC", 0.15 }, { "UCA", 0.12 }, { "UCG", 0.15 }, { "AGU", 0.15 }, { "AGC", 0.28 },
            // Proline (P)
            { "CCU", 0.16 }, { "CCC", 0.12 }, { "CCA", 0.19 }, { "CCG", 0.53 },
            // Threonine (T)
            { "ACU", 0.16 }, { "ACC", 0.44 }, { "ACA", 0.13 }, { "ACG", 0.27 },
            // Alanine (A)
            { "GCU", 0.16 }, { "GCC", 0.27 }, { "GCA", 0.21 }, { "GCG", 0.36 },
            // Tyrosine (Y)
            { "UAU", 0.57 }, { "UAC", 0.43 },
            // Stop (*)
            { "UAA", 0.64 }, { "UAG", 0.07 }, { "UGA", 0.29 },
            // Histidine (H)
            { "CAU", 0.57 }, { "CAC", 0.43 },
            // Glutamine (Q)
            { "CAA", 0.35 }, { "CAG", 0.65 },
            // Asparagine (N)
            { "AAU", 0.45 }, { "AAC", 0.55 },
            // Lysine (K)
            { "AAA", 0.76 }, { "AAG", 0.24 },
            // Aspartic acid (D)
            { "GAU", 0.63 }, { "GAC", 0.37 },
            // Glutamic acid (E)
            { "GAA", 0.69 }, { "GAG", 0.31 },
            // Cysteine (C)
            { "UGU", 0.44 }, { "UGC", 0.56 },
            // Tryptophan (W)
            { "UGG", 1.00 },
            // Arginine (R)
            { "CGU", 0.38 }, { "CGC", 0.40 }, { "CGA", 0.06 }, { "CGG", 0.10 }, { "AGA", 0.04 }, { "AGG", 0.02 },
            // Glycine (G)
            { "GGU", 0.34 }, { "GGC", 0.41 }, { "GGA", 0.11 }, { "GGG", 0.15 }
        },
        StandardCodonToAminoAcid);

    /// <summary>
    /// Saccharomyces cerevisiae (yeast) codon usage frequencies (relative fraction per amino acid).
    /// Source: Kazusa Codon Usage Database, species=4932.
    /// URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=4932
    /// </summary>
    /// <remarks>
    /// <see cref="CodonUsageTable.PerThousand"/>: the same Kazusa species=4932 entry, frequency per
    /// thousand codons as published (one decimal), verbatim from the Clark-lab CHARMING repository
    /// file <c>ScerCUB.txt</c> (github.com/wrightgs/CHARMING, commit 2c75d1b). Its per-amino-acid
    /// ratios round to these fractions at 62/64 codons; the stop family (UAA 1.1, UAG 0.5, UGA 0.7)
    /// differs by ≤ 0.013 only because its per-thousand values are rounded to 0.1.
    /// </remarks>
    public static readonly CodonUsageTable Yeast = new(
        "Saccharomyces cerevisiae",
        new Dictionary<string, double>
        {
            { "UUU", 0.59 }, { "UUC", 0.41 },
            { "UUA", 0.28 }, { "UUG", 0.29 }, { "CUU", 0.13 }, { "CUC", 0.06 }, { "CUA", 0.14 }, { "CUG", 0.11 },
            { "AUU", 0.46 }, { "AUC", 0.26 }, { "AUA", 0.27 },
            { "AUG", 1.00 },
            { "GUU", 0.39 }, { "GUC", 0.21 }, { "GUA", 0.21 }, { "GUG", 0.19 },
            { "UCU", 0.26 }, { "UCC", 0.16 }, { "UCA", 0.21 }, { "UCG", 0.10 }, { "AGU", 0.16 }, { "AGC", 0.11 },
            { "CCU", 0.31 }, { "CCC", 0.15 }, { "CCA", 0.42 }, { "CCG", 0.12 },
            { "ACU", 0.35 }, { "ACC", 0.22 }, { "ACA", 0.30 }, { "ACG", 0.14 },
            { "GCU", 0.38 }, { "GCC", 0.22 }, { "GCA", 0.29 }, { "GCG", 0.11 },
            { "UAU", 0.56 }, { "UAC", 0.44 },
            { "UAA", 0.47 }, { "UAG", 0.23 }, { "UGA", 0.30 },
            { "CAU", 0.64 }, { "CAC", 0.36 },
            { "CAA", 0.69 }, { "CAG", 0.31 },
            { "AAU", 0.59 }, { "AAC", 0.41 },
            { "AAA", 0.58 }, { "AAG", 0.42 },
            { "GAU", 0.65 }, { "GAC", 0.35 },
            { "GAA", 0.70 }, { "GAG", 0.30 },
            { "UGU", 0.63 }, { "UGC", 0.37 },
            { "UGG", 1.00 },
            { "CGU", 0.14 }, { "CGC", 0.06 }, { "CGA", 0.07 }, { "CGG", 0.04 }, { "AGA", 0.48 }, { "AGG", 0.21 },
            { "GGU", 0.47 }, { "GGC", 0.19 }, { "GGA", 0.22 }, { "GGG", 0.12 }
        },
        StandardCodonToAminoAcid)
    {
        PerThousand = YeastKazusaPerThousand(),
    };

    /// <summary>
    /// Human codon usage frequencies (relative fraction per amino acid).
    /// Source: Kazusa Codon Usage Database, species=9606.
    /// URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
    /// </summary>
    /// <remarks>
    /// <see cref="CodonUsageTable.PerThousand"/>: 1000 × count / 40 662 582 from the codon counts
    /// ("Number" column) of the same Kazusa species=9606 entry (93 487 CDS, 40 662 582 codons), as
    /// shipped in the PyPI package codon-optimizer 0.2.0
    /// (<c>reference_seq/Homo_sapiens_codon_frequency.csv</c>, citing kazusa.or.jp/codon). The
    /// per-amino-acid ratios of these counts round to the fractions below at all 64 codons.
    /// </remarks>
    public static readonly CodonUsageTable Human = new(
        "Homo sapiens",
        new Dictionary<string, double>
        {
            { "UUU", 0.46 }, { "UUC", 0.54 },
            { "UUA", 0.08 }, { "UUG", 0.13 }, { "CUU", 0.13 }, { "CUC", 0.20 }, { "CUA", 0.07 }, { "CUG", 0.40 },
            { "AUU", 0.36 }, { "AUC", 0.47 }, { "AUA", 0.17 },
            { "AUG", 1.00 },
            { "GUU", 0.18 }, { "GUC", 0.24 }, { "GUA", 0.12 }, { "GUG", 0.46 },
            { "UCU", 0.19 }, { "UCC", 0.22 }, { "UCA", 0.15 }, { "UCG", 0.05 }, { "AGU", 0.15 }, { "AGC", 0.24 },
            { "CCU", 0.29 }, { "CCC", 0.32 }, { "CCA", 0.28 }, { "CCG", 0.11 },
            { "ACU", 0.25 }, { "ACC", 0.36 }, { "ACA", 0.28 }, { "ACG", 0.11 },
            { "GCU", 0.27 }, { "GCC", 0.40 }, { "GCA", 0.23 }, { "GCG", 0.11 },
            { "UAU", 0.44 }, { "UAC", 0.56 },
            { "UAA", 0.30 }, { "UAG", 0.24 }, { "UGA", 0.47 },
            { "CAU", 0.42 }, { "CAC", 0.58 },
            { "CAA", 0.27 }, { "CAG", 0.73 },
            { "AAU", 0.47 }, { "AAC", 0.53 },
            { "AAA", 0.43 }, { "AAG", 0.57 },
            { "GAU", 0.46 }, { "GAC", 0.54 },
            { "GAA", 0.42 }, { "GAG", 0.58 },
            { "UGU", 0.46 }, { "UGC", 0.54 },
            { "UGG", 1.00 },
            { "CGU", 0.08 }, { "CGC", 0.18 }, { "CGA", 0.11 }, { "CGG", 0.20 }, { "AGA", 0.21 }, { "AGG", 0.21 },
            { "GGU", 0.16 }, { "GGC", 0.34 }, { "GGA", 0.25 }, { "GGG", 0.25 }
        },
        StandardCodonToAminoAcid)
    {
        PerThousand = HumanKazusaPerThousand(),
    };

    // Kazusa species=4932 frequency per thousand codons (one decimal, as published), verbatim from
    // CHARMING ScerCUB.txt (Clark lab; github.com/wrightgs/CHARMING). A method, not a static field,
    // so the preset initialiser above does not depend on static field order.
    private static IReadOnlyDictionary<string, double> YeastKazusaPerThousand() => new Dictionary<string, double>
    {
            { "UUU", 26.1 }, { "UUC", 18.4 }, { "UUA", 26.2 }, { "UUG", 27.2 }, { "UCU", 23.5 }, { "UCC", 14.2 }, { "UCA", 18.7 }, { "UCG", 8.6 },
            { "UAU", 18.8 }, { "UAC", 14.8 }, { "UAA", 1.1 }, { "UAG", 0.5 }, { "UGU", 8.1 }, { "UGC", 4.8 }, { "UGA", 0.7 }, { "UGG", 10.4 },
            { "CUU", 12.3 }, { "CUC", 5.4 }, { "CUA", 13.4 }, { "CUG", 10.5 }, { "CCU", 13.5 }, { "CCC", 6.8 }, { "CCA", 18.3 }, { "CCG", 5.3 },
            { "CAU", 13.6 }, { "CAC", 7.8 }, { "CAA", 27.3 }, { "CAG", 12.1 }, { "CGU", 6.4 }, { "CGC", 2.6 }, { "CGA", 3.0 }, { "CGG", 1.7 },
            { "AUU", 30.1 }, { "AUC", 17.2 }, { "AUA", 17.8 }, { "AUG", 20.9 }, { "ACU", 20.3 }, { "ACC", 12.7 }, { "ACA", 17.8 }, { "ACG", 8.0 },
            { "AAU", 35.7 }, { "AAC", 24.8 }, { "AAA", 41.9 }, { "AAG", 30.8 }, { "AGU", 14.2 }, { "AGC", 9.8 }, { "AGA", 21.3 }, { "AGG", 9.2 },
            { "GUU", 22.1 }, { "GUC", 11.8 }, { "GUA", 11.8 }, { "GUG", 10.8 }, { "GCU", 21.2 }, { "GCC", 12.6 }, { "GCA", 16.2 }, { "GCG", 6.2 },
            { "GAU", 37.6 }, { "GAC", 20.2 }, { "GAA", 45.6 }, { "GAG", 19.2 }, { "GGU", 23.9 }, { "GGC", 9.8 }, { "GGA", 10.9 }, { "GGG", 6.0 },
    };

    // Kazusa species=9606 codon counts ("Number" column; 40 662 582 codons), as shipped in PyPI
    // codon-optimizer 0.2.0 reference_seq/Homo_sapiens_codon_frequency.csv; per thousand =
    // 1000 * count / total.
    private static IReadOnlyDictionary<string, double> HumanKazusaPerThousand()
    {
        var counts = new Dictionary<string, long>
        {
                { "UUU", 714298 }, { "UUC", 824692 }, { "UUA", 311881 }, { "UUG", 525688 }, { "UCU", 618711 }, { "UCC", 718892 }, { "UCA", 496448 }, { "UCG", 179419 },
                { "UAU", 495699 }, { "UAC", 622407 }, { "UAA", 40285 }, { "UAG", 32109 }, { "UGU", 430311 }, { "UGC", 513028 }, { "UGA", 63237 }, { "UGG", 535595 },
                { "CUU", 536515 }, { "CUC", 796638 }, { "CUA", 290751 }, { "CUG", 1611801 }, { "CCU", 713233 }, { "CCC", 804620 }, { "CCA", 688038 }, { "CCG", 281570 },
                { "CAU", 441711 }, { "CAC", 613713 }, { "CAA", 501911 }, { "CAG", 1391973 }, { "CGU", 184609 }, { "CGC", 423516 }, { "CGA", 250760 }, { "CGG", 464485 },
                { "AUU", 650473 }, { "AUC", 846466 }, { "AUA", 304565 }, { "AUG", 896005 }, { "ACU", 533609 }, { "ACC", 768147 }, { "ACA", 614523 }, { "ACG", 246105 },
                { "AAU", 689701 }, { "AAC", 776603 }, { "AAA", 993621 }, { "AAG", 1295568 }, { "AGU", 493429 }, { "AGC", 791383 }, { "AGA", 494682 }, { "AGG", 486463 },
                { "GUU", 448607 }, { "GUC", 588138 }, { "GUA", 287712 }, { "GUG", 1143534 }, { "GCU", 750096 }, { "GCC", 1127679 }, { "GCA", 643471 }, { "GCG", 299495 },
                { "GAU", 885429 }, { "GAC", 1020595 }, { "GAA", 1177632 }, { "GAG", 1609975 }, { "GGU", 437126 }, { "GGC", 903565 }, { "GGA", 669873 }, { "GGG", 669768 },
        };
        long total = counts.Values.Sum();
        return counts.ToDictionary(kv => kv.Key, kv => 1000.0 * kv.Value / total);
    }

    #endregion

    #region Codon Optimization

    /// <summary>
    /// Optimizes a coding sequence for expression in a target organism by replacing codons with
    /// synonymous ones, preserving the encoded protein.
    /// </summary>
    /// <remarks>
    /// The input is upper-cased, read as RNA (T → U) and trimmed to complete codons. Stop codons
    /// are kept as they are; codons of single-codon amino acids (Met AUG, Trp UGG) have no
    /// synonym and never change; a triplet that is not a valid IUPAC codon, or whose amino acid
    /// is itself ambiguous (B, Z, J, X), is left untouched and contributes 'X' to
    /// <see cref="OptimizationResult.ProteinSequence"/>. Amino acids and synonymous families are
    /// taken from the canonical <see cref="GeneticCode"/> table selected by
    /// <paramref name="geneticCode"/> (default <see cref="GeneticCode.Standard"/>), so the
    /// protein is preserved under that code (<c>Translator.Translate(…, geneticCode)</c>).
    /// <para>Strategies:</para>
    /// <list type="bullet">
    /// <item><description><b>MaximizeCAI</b> — every codon becomes the most frequent synonymous
    /// codon of the target table ("one amino acid – one codon"; Puigbò et&#160;al. 2007 OPTIMIZER,
    /// DNA Chisel <c>use_best_codon</c>/<c>MaximizeCAI</c>, Biopython
    /// <c>CodonAdaptationIndex.optimize</c>). Ties are broken by NCBI codon order.</description></item>
    /// <item><description><b>AvoidRareCodeons</b> — only codons whose table frequency is strictly
    /// below <paramref name="rareCodonThreshold"/> are replaced, by the most frequent synonymous
    /// codon; codons at or above the threshold are kept.</description></item>
    /// <item><description><b>BalancedOptimization</b> (default) — as MaximizeCAI, then a
    /// GC-balancing pass moves the overall GC fraction into
    /// [<paramref name="gcTargetMin"/>, <paramref name="gcTargetMax"/>] with further synonymous
    /// swaps (DNA Chisel <c>EnforceGCContent(mini, maxi)</c> combined with codon
    /// optimisation).</description></item>
    /// <item><description><b>MinimizeSecondary</b> — BalancedOptimization (codon selection + GC
    /// pass) followed by the structure pass
    /// <see cref="ReduceSecondaryStructure(string, CodonUsageTable, int, double, double, SecondaryStructureScoring, double, int)"/>
    /// with its default window settings, scored by <paramref name="secondaryScoring"/>. The
    /// default is the pair-count heuristic, which costs ~0.06 s per kilobase; pass
    /// <see cref="SecondaryStructureScoring.Turner2004Mfe"/> for the sourced Turner 2004
    /// free-energy model (~9 s per kilobase).</description></item>
    /// <item><description><b>HarmonizeExpression</b> — the codon usage of the output matches the
    /// target table as closely as integer codon counts allow (DNA Chisel
    /// <c>match_codon_usage</c> / <c>MatchTargetCodonUsage</c>), deterministically: per amino
    /// acid the number of occurrences of each synonymous codon is the largest-remainder rounding
    /// of frequency × (number of residues of that amino acid), and positions that already carry
    /// an allotted codon keep it, so the edit count is minimal.</description></item>
    /// </list>
    /// </remarks>
    /// <param name="codingSequence">Coding sequence (DNA or RNA, any case); empty → empty result.</param>
    /// <param name="targetOrganism">Target codon-usage table (per-amino-acid relative fractions).</param>
    /// <param name="strategy">Codon-selection strategy (see remarks).</param>
    /// <param name="gcTargetMin">Lower bound of the target GC fraction (BalancedOptimization / MinimizeSecondary).</param>
    /// <param name="gcTargetMax">Upper bound of the target GC fraction (BalancedOptimization / MinimizeSecondary).</param>
    /// <param name="rareCodonThreshold">Frequency below which a codon counts as rare (strict &lt;).</param>
    /// <param name="secondaryScoring">Window scoring model of the MinimizeSecondary structure pass (default <see cref="SecondaryStructureScoring.PairCount"/>); ignored by the other strategies.</param>
    /// <param name="gcWindow">
    /// <see langword="null"/> (default): the GC pass constrains the <i>global</i> GC fraction only
    /// (DNA Chisel <c>EnforceGCContent(mini, maxi)</c>). A positive value: the GC pass constrains
    /// every sliding window of that many nucleotides to
    /// [<paramref name="gcTargetMin"/>, <paramref name="gcTargetMax"/>], i.e. DNA Chisel
    /// <c>EnforceGCContent(mini, maxi, window=gcWindow)</c>; see
    /// <see cref="BalanceGcContentWindowed"/> and <see cref="EvaluateGcContent"/>.
    /// </param>
    /// <param name="geneticCode">
    /// NCBI genetic code that defines the encoded protein and the synonymous-codon families
    /// (<see langword="null"/>, the default: <see cref="GeneticCode.Standard"/>, identical to the
    /// historical output). DNA Chisel <c>EnforceTranslation(genetic_table=…)</c>: under table 2
    /// (vertebrate mitochondrial) AGA/AGG are stops and are never used for Arg, UGA is Trp and AUA
    /// is Met; dual-coding codons of tables 27/28/31 belong to the family of the amino acid they
    /// encode (Biopython <c>forward_table</c>, DNA Chisel back-translation table).
    /// </param>
    public static OptimizationResult OptimizeSequence(
        string codingSequence,
        CodonUsageTable targetOrganism,
        OptimizationStrategy strategy = OptimizationStrategy.BalancedOptimization,
        double gcTargetMin = 0.40,
        double gcTargetMax = 0.60,
        double rareCodonThreshold = 0.15,
        SecondaryStructureScoring secondaryScoring = SecondaryStructureScoring.PairCount,
        int? gcWindow = null,
        GeneticCode? geneticCode = null)
    {
        if (gcWindow is <= 0)
            throw new ArgumentOutOfRangeException(nameof(gcWindow), gcWindow, "The GC window must be a positive number of nucleotides.");

        if (string.IsNullOrEmpty(codingSequence))
        {
            return new OptimizationResult("", "", "", 0, 0, 0, 0, 0, new List<(int, string, string)>());
        }

        string rna = ToUpperRna(codingSequence);

        if (rna.Length % 3 != 0)
        {
            // Trim to complete codons
            rna = rna.Substring(0, (rna.Length / 3) * 3);
        }

        var originalCodons = SplitIntoCodons(rna);
        var families = FamiliesFor(geneticCode);

        double originalCAI = CalculateCAI(rna, targetOrganism, geneticCode: geneticCode);
        var proteinBuilder = new StringBuilder();
        foreach (string codon in originalCodons)
            proteinBuilder.Append(TryTranslateCodon(codon, families, out char aa) ? aa : 'X');

        List<string> optimizedCodons = strategy == OptimizationStrategy.HarmonizeExpression
            ? MatchTargetCodonUsage(originalCodons, targetOrganism, families)
            : originalCodons
                .Select(c => SelectOptimalCodon(c, targetOrganism, strategy, rareCodonThreshold, families))
                .ToList();

        // GC balancing for the strategies that declare a GC target.
        if (strategy is OptimizationStrategy.BalancedOptimization or OptimizationStrategy.MinimizeSecondary)
        {
            if (gcWindow is int windowSize)
                BalanceGcContentWindowed(optimizedCodons, targetOrganism, gcTargetMin, gcTargetMax, windowSize, rareCodonThreshold, families);
            else
                BalanceGcContent(optimizedCodons, targetOrganism, gcTargetMin, gcTargetMax, rareCodonThreshold, families);
        }

        string optimizedSequence = string.Concat(optimizedCodons);

        // MinimizeSecondary additionally runs the dedicated structure pass on the GC-balanced
        // sequence, so the strategy actually does what its name says.
        if (strategy == OptimizationStrategy.MinimizeSecondary)
        {
            optimizedSequence = ReduceSecondaryStructure(
                optimizedSequence, targetOrganism, scoring: secondaryScoring, geneticCode: geneticCode);
            optimizedCodons = SplitIntoCodons(optimizedSequence);
        }

        var changes = new List<(int Position, string Original, string Optimized)>();
        for (int i = 0; i < originalCodons.Count; i++)
        {
            if (originalCodons[i] != optimizedCodons[i])
                changes.Add((i * 3, originalCodons[i], optimizedCodons[i]));
        }

        double optimizedCAI = CalculateCAI(optimizedSequence, targetOrganism, geneticCode: geneticCode);

        return new OptimizationResult(
            OriginalSequence: rna,
            OptimizedSequence: optimizedSequence,
            ProteinSequence: proteinBuilder.ToString(),
            OriginalCAI: originalCAI,
            OptimizedCAI: optimizedCAI,
            GcContentOriginal: CalculateGcContent(rna),
            GcContentOptimized: CalculateGcContent(optimizedSequence),
            ChangedCodons: changes.Count,
            Changes: changes);
    }

    // Most frequent synonymous codon of <paramref name="codon"/> in <paramref name="table"/>;
    // ties (and codons without usable synonyms) resolve to the first codon in NCBI order.
    private static string BestSynonymousCodon(string codon, CodonUsageTable table, CodeFamilies families)
    {
        var synonyms = SynonymsFor(codon, families);
        if (synonyms.Length == 0)
            return codon;

        string best = synonyms[0];
        double bestFrequency = table.CodonFrequencies.GetValueOrDefault(best, 0);
        for (int i = 1; i < synonyms.Length; i++)
        {
            double frequency = table.CodonFrequencies.GetValueOrDefault(synonyms[i], 0);
            if (frequency > bestFrequency)
            {
                best = synonyms[i];
                bestFrequency = frequency;
            }
        }

        return best;
    }

    private static string SelectOptimalCodon(
        string currentCodon, CodonUsageTable table, OptimizationStrategy strategy, double rareCodonThreshold, CodeFamilies families)
    {
        if (!TryTranslateCodon(currentCodon, families, out char aminoAcid) || aminoAcid == '*')
            return currentCodon; // stop codons and unusable triplets are preserved

        if (strategy == OptimizationStrategy.AvoidRareCodeons &&
            table.CodonFrequencies.GetValueOrDefault(currentCodon, 0) >= rareCodonThreshold)
        {
            return currentCodon; // not rare → untouched
        }

        return BestSynonymousCodon(currentCodon, table, families);
    }

    /// <summary>
    /// Rewrites the codons so that, per amino acid, the codon counts are the largest-remainder
    /// rounding of the target table's relative frequencies — the deterministic optimum of DNA
    /// Chisel's <c>MatchTargetCodonUsage</c> objective (score = −Σ_aa n_aa·Σ_codon |f_seq − f_table|),
    /// which that library reaches by randomised local search. Positions already carrying an
    /// allotted codon keep it, so only the surplus positions are edited. Stop codons and triplets
    /// without a usable synonymous family are left unchanged.
    /// </summary>
    private static List<string> MatchTargetCodonUsage(List<string> codons, CodonUsageTable table, CodeFamilies families)
    {
        var result = new List<string>(codons);

        // Group the editable positions by amino acid (NCBI-ordered families).
        var positionsByAminoAcid = new Dictionary<char, List<int>>();
        for (int i = 0; i < codons.Count; i++)
        {
            if (!TryTranslateCodon(codons[i], families, out char aa) || aa == '*' || !families.Synonyms.ContainsKey(aa))
                continue;
            if (!positionsByAminoAcid.TryGetValue(aa, out var positions))
                positionsByAminoAcid[aa] = positions = new List<int>();
            positions.Add(i);
        }

        foreach (var (aminoAcid, positions) in positionsByAminoAcid)
        {
            var family = families.Synonyms[aminoAcid];
            int n = positions.Count;

            double total = family.Sum(c => table.CodonFrequencies.GetValueOrDefault(c, 0));
            if (total <= 0)
                continue; // no usage data for this amino acid → leave the codons alone

            // Largest-remainder (Hamilton) allocation of n positions over the family.
            var quota = new double[family.Length];
            var allotted = new int[family.Length];
            int assigned = 0;
            for (int k = 0; k < family.Length; k++)
            {
                quota[k] = n * table.CodonFrequencies.GetValueOrDefault(family[k], 0) / total;
                allotted[k] = (int)Math.Floor(quota[k]);
                assigned += allotted[k];
            }

            foreach (int k in Enumerable.Range(0, family.Length)
                         .OrderByDescending(k => quota[k] - Math.Floor(quota[k]))
                         .ThenBy(k => k)
                         .Take(Math.Max(0, n - assigned)))
            {
                allotted[k]++;
            }

            // Keep positions that already carry an allotted codon (minimal edit), then fill the rest.
            var free = new List<int>();
            var remaining = (int[])allotted.Clone();
            foreach (int position in positions)
            {
                int k = Array.IndexOf(family, codons[position]);
                if (k >= 0 && remaining[k] > 0)
                    remaining[k]--;
                else
                    free.Add(position);
            }

            int next = 0;
            foreach (int position in free)
            {
                while (next < family.Length && remaining[next] == 0)
                    next++;
                if (next >= family.Length)
                    break;
                result[position] = family[next];
                remaining[next]--;
            }
        }

        return result;
    }

    /// <summary>
    /// Moves the overall GC fraction of <paramref name="codons"/> into
    /// [<paramref name="minGc"/>, <paramref name="maxGc"/>] with synonymous codon swaps (DNA
    /// Chisel <c>EnforceGCContent(mini, maxi)</c> resolved by synonymous mutations). Only swaps
    /// that actually move GC toward the target are applied, candidates that would overshoot the
    /// opposite bound are used only when no in-range candidate exists, and the pass stops as soon
    /// as the sequence is inside the window. Codons whose frequency in
    /// <paramref name="table"/> is below <paramref name="minCodonFrequency"/> are not used.
    /// </summary>
    private static void BalanceGcContent(
        List<string> codons,
        CodonUsageTable table,
        double minGc,
        double maxGc,
        double minCodonFrequency,
        CodeFamilies families)
    {
        int length = codons.Count * 3;
        if (length == 0)
            return;

        int gc = codons.Sum(CountGc);
        if (IsWithin(gc)) return;

        for (int i = 0; i < codons.Count; i++)
        {
            bool needMoreGc = gc < minGc * length;
            string current = codons[i];
            if (!TryTranslateCodon(current, families, out char aa) || aa == '*')
                continue;

            int currentGc = CountGc(current);
            string? best = null;
            int bestGc = currentGc;
            double bestFrequency = double.NegativeInfinity;
            bool bestInRange = false;

            foreach (string alternative in SynonymsFor(current, families))
            {
                if (alternative == current)
                    continue;
                double frequency = table.CodonFrequencies.GetValueOrDefault(alternative, 0);
                if (frequency < minCodonFrequency)
                    continue;

                int alternativeGc = CountGc(alternative);
                // The swap must move GC in the required direction.
                if (needMoreGc ? alternativeGc <= currentGc : alternativeGc >= currentGc)
                    continue;

                bool inRange = IsWithin(gc - currentGc + alternativeGc);
                // Prefer a candidate that lands inside the window; among those, the most frequent
                // codon (least CAI cost). Otherwise the one that moves GC furthest toward the
                // window, then the most frequent codon.
                bool better = (inRange, inRange ? 0 : Math.Abs(alternativeGc - currentGc), frequency)
                    .CompareTo((bestInRange, bestInRange ? 0 : Math.Abs(bestGc - currentGc), bestFrequency)) > 0;
                if (best is null || better)
                {
                    best = alternative;
                    bestGc = alternativeGc;
                    bestFrequency = frequency;
                    bestInRange = inRange;
                }
            }

            if (best is null)
                continue;

            codons[i] = best;
            gc += bestGc - currentGc;
            if (IsWithin(gc))
                return;
        }

        bool IsWithin(int gcCount) => gcCount >= minGc * length && gcCount <= maxGc * length;
    }

    /// <summary>
    /// A breach segment of a GC-content constraint, as DNA Chisel's
    /// <c>EnforceGCContent.evaluate</c> reports it: a half-open nucleotide interval
    /// [<paramref name="Start"/>, <paramref name="End"/>).
    /// </summary>
    /// <param name="Start">0-based start of the breach segment.</param>
    /// <param name="End">Exclusive end of the breach segment.</param>
    public sealed record GcBreachLocation(int Start, int End);

    /// <summary>
    /// Result of evaluating a GC-content constraint on a sequence
    /// (DNA Chisel <c>EnforceGCContent(mini, maxi, window).evaluate(problem)</c>).
    /// </summary>
    /// <param name="Score">
    /// Negated sum of the out-of-bound deviations over all evaluated windows
    /// (0 when every window satisfies the bounds), i.e. DNA Chisel's
    /// <c>SpecEvaluation.score</c>.
    /// </param>
    /// <param name="BreachWindowStarts">0-based start positions of the windows that breach the bounds.</param>
    /// <param name="BreachLocations">
    /// Breach segments after grouping neighbouring breaching windows, as DNA Chisel reports
    /// <c>SpecEvaluation.locations</c>.
    /// </param>
    public sealed record GcContentEvaluation(
        double Score,
        IReadOnlyList<int> BreachWindowStarts,
        IReadOnlyList<GcBreachLocation> BreachLocations);

    /// <summary>
    /// Evaluates a GC-content constraint on <paramref name="sequence"/> exactly as DNA Chisel's
    /// <c>EnforceGCContent(mini, maxi, window).evaluate(problem)</c> does (DNA Chisel 3.2.16,
    /// <c>builtin_specifications/EnforceGCContent.py</c> + <c>biotools/gc_content.py</c>).
    /// </summary>
    /// <remarks>
    /// With <paramref name="window"/> = <see langword="null"/> the global GC fraction
    /// (G+C over sequence length) is compared with the bounds and a breach is reported as the whole
    /// sequence. With a window size <c>w</c> the GC fraction of every sliding window
    /// <c>[i, i+w)</c>, <c>i = 0 … len-w</c>, is compared with the bounds (so there are
    /// <c>len - w + 1</c> windows; a window longer than the sequence yields no windows, score 0).
    /// Per window the breach extent is <c>max(0, mini - gc) + max(0, gc - maxi)</c>; both bounds are
    /// inclusive (a window exactly at a bound does not breach) and
    /// <see cref="GcContentEvaluation.Score"/> is the negated sum of those extents. Breaching
    /// windows are grouped the way DNA Chisel groups them (<c>group_nearby_segments</c> with
    /// <c>max_start_spread = max(1, locationsSpan)</c>: a window joins the current group while its
    /// start is less than <paramref name="locationsSpan"/> nucleotides after the group's first
    /// start), and each group is reported as [first start, last start + w).
    /// G and C are counted by the canonical counting primitive, so lower case and RNA (U) input is
    /// handled; any other symbol simply does not count as G/C, as in DNA Chisel.
    /// </remarks>
    /// <param name="sequence">Sequence to evaluate (DNA or RNA, any case).</param>
    /// <param name="mini">Lower bound of the GC fraction (inclusive).</param>
    /// <param name="maxi">Upper bound of the GC fraction (inclusive).</param>
    /// <param name="window">Sliding-window size in nucleotides, or <see langword="null"/> for global GC.</param>
    /// <param name="locationsSpan">
    /// Grouping resolution of the reported breach segments (DNA Chisel
    /// <c>EnforceGCContent.locations_span</c>, 50).
    /// </param>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is <see langword="null"/>.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="window"/> is not positive.</exception>
    public static GcContentEvaluation EvaluateGcContent(
        string sequence,
        double mini = 0.0,
        double maxi = 1.0,
        int? window = null,
        int locationsSpan = 50)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        if (window is <= 0)
            throw new ArgumentOutOfRangeException(nameof(window), window, "The GC window must be a positive number of nucleotides.");

        var noBreaches = new GcContentEvaluation(0, Array.Empty<int>(), Array.Empty<GcBreachLocation>());
        if (sequence.Length == 0)
            return noBreaches;

        if (window is not int w)
        {
            double globalGc = (double)sequence.AsSpan().CountGcAndValidNucleotides().GcCount / sequence.Length;
            double globalBreach = Math.Max(0, mini - globalGc) + Math.Max(0, globalGc - maxi);
            return globalBreach <= 0
                ? noBreaches
                : new GcContentEvaluation(-globalBreach, new[] { 0 }, new[] { new GcBreachLocation(0, sequence.Length) });
        }

        if (w > sequence.Length)
            return noBreaches; // no window fits: DNA Chisel evaluates an empty window array

        // Prefix sums of G/C indicators, as dnachisel's cumulative-sum window GC.
        var prefix = new int[sequence.Length + 1];
        for (int i = 0; i < sequence.Length; i++)
        {
            char c = char.ToUpperInvariant(sequence[i]);
            prefix[i + 1] = prefix[i] + (c is 'G' or 'C' ? 1 : 0);
        }

        double totalBreach = 0;
        var starts = new List<int>();
        for (int start = 0; start + w <= sequence.Length; start++)
        {
            double gc = (double)(prefix[start + w] - prefix[start]) / w;
            double breach = Math.Max(0, mini - gc) + Math.Max(0, gc - maxi);
            if (breach > 0)
            {
                totalBreach += breach;
                starts.Add(start);
            }
        }

        if (starts.Count == 0)
            return noBreaches;

        var locations = new List<GcBreachLocation>();
        if (starts.Count == 1)
        {
            locations.Add(new GcBreachLocation(starts[0], starts[0] + w));
        }
        else
        {
            int spread = Math.Max(1, locationsSpan);
            int groupStart = starts[0];
            int lastStart = starts[0];
            for (int k = 1; k < starts.Count; k++)
            {
                if (starts[k] - groupStart < spread)
                {
                    lastStart = starts[k];
                    continue;
                }

                locations.Add(new GcBreachLocation(groupStart, lastStart + w));
                groupStart = lastStart = starts[k];
            }

            locations.Add(new GcBreachLocation(groupStart, lastStart + w));
        }

        return new GcContentEvaluation(-totalBreach, starts, locations);
    }

    /// <summary>
    /// Brings every sliding window of <paramref name="window"/> nucleotides of
    /// <paramref name="codons"/> into [<paramref name="minGc"/>, <paramref name="maxGc"/>] with
    /// synonymous codon swaps — DNA Chisel's <c>EnforceGCContent(mini, maxi, window=…)</c>
    /// constraint resolved by synonymous mutations, with a deterministic greedy search instead of
    /// DNA Chisel's randomised constraint solver.
    /// </summary>
    /// <remarks>
    /// Search: repeatedly take the left-most window that breaches the bounds (breach extent and
    /// window geometry exactly as in <see cref="EvaluateGcContent"/>) and try every synonymous
    /// replacement of every codon overlapping that window whose frequency in
    /// <paramref name="table"/> is at least <paramref name="minCodonFrequency"/>. The replacement
    /// that lowers the total breach of the whole sequence the most is applied; ties are broken by
    /// the higher codon frequency, then the left-most codon position, then NCBI codon order, so the
    /// output is deterministic. When no replacement improves the total breach, that window is left
    /// alone (its bounds are unreachable with synonymous codons) and the search moves on to the
    /// next breaching window. Stop codons and triplets without a usable synonymous family are never
    /// changed, so the encoded protein is preserved.
    /// </remarks>
    private static void BalanceGcContentWindowed(
        List<string> codons,
        CodonUsageTable table,
        double minGc,
        double maxGc,
        int window,
        double minCodonFrequency,
        CodeFamilies families)
    {
        int length = codons.Count * 3;
        if (length == 0 || window <= 0 || window > length)
            return;

        var gcFlags = new int[length];
        for (int i = 0; i < codons.Count; i++)
        {
            string codon = codons[i];
            for (int k = 0; k < 3 && k < codon.Length; k++)
            {
                char c = char.ToUpperInvariant(codon[k]);
                gcFlags[i * 3 + k] = c is 'G' or 'C' ? 1 : 0;
            }
        }

        int windowCount = length - window + 1;

        double WindowBreach(int start)
        {
            int gc = 0;
            for (int p = start; p < start + window; p++)
                gc += gcFlags[p];
            double fraction = (double)gc / window;
            return Math.Max(0, minGc - fraction) + Math.Max(0, fraction - maxGc);
        }

        double TotalBreach()
        {
            double total = 0;
            for (int start = 0; start < windowCount; start++)
                total += WindowBreach(start);
            return total;
        }

        void WriteCodon(int index, string codon)
        {
            codons[index] = codon;
            for (int k = 0; k < 3; k++)
            {
                char c = char.ToUpperInvariant(codon[k]);
                gcFlags[index * 3 + k] = c is 'G' or 'C' ? 1 : 0;
            }
        }

        // Rounds of a left-to-right pass: a repair can re-break an already visited window, so the
        // pass is repeated while the total breach strictly decreases (at most MaxRounds rounds).
        const int MaxRounds = 8;
        double previousTotal = double.PositiveInfinity;
        for (int round = 0; round < MaxRounds; round++)
        {
            double totalBefore = TotalBreach();
            if (totalBefore <= 0 || totalBefore >= previousTotal - 1e-12)
                return;
            previousTotal = totalBefore;

            int frozen = 0;
            int guard = Math.Max(16, codons.Count * 16);
            while (frozen < windowCount && guard-- > 0)
            {
                int breaching = -1;
                for (int start = frozen; start < windowCount; start++)
                {
                    if (WindowBreach(start) > 0)
                    {
                        breaching = start;
                        break;
                    }
                }

                if (breaching < 0)
                    break;

                double windowBefore = WindowBreach(breaching);
                double totalCurrent = TotalBreach();
                int firstCodon = breaching / 3;
                int lastCodon = Math.Min(codons.Count - 1, (breaching + window - 1) / 3);

                int bestCodonIndex = -1;
                string? bestCodon = null;
                double bestWindow = windowBefore;
                double bestTotal = totalCurrent;
                double bestFrequency = double.NegativeInfinity;

                for (int index = firstCodon; index <= lastCodon; index++)
                {
                    string current = codons[index];
                    if (!TryTranslateCodon(current, families, out char aa) || aa == '*')
                        continue;

                    foreach (string alternative in SynonymsFor(current, families))
                    {
                        if (alternative == current)
                            continue;
                        double frequency = table.CodonFrequencies.GetValueOrDefault(alternative, 0);
                        if (frequency < minCodonFrequency)
                            continue;

                        WriteCodon(index, alternative);
                        double windowAfter = WindowBreach(breaching);
                        double totalAfter = TotalBreach();
                        WriteCodon(index, current);

                        // Only moves that improve the targeted window, or keep it and lower the
                        // total breach, are candidates (lexicographic objective).
                        bool improves = windowAfter < windowBefore - 1e-12 ||
                            (Math.Abs(windowAfter - windowBefore) <= 1e-12 && totalAfter < totalCurrent - 1e-12);
                        if (!improves)
                            continue;

                        bool better = bestCodon is null ||
                            windowAfter < bestWindow - 1e-12 ||
                            (Math.Abs(windowAfter - bestWindow) <= 1e-12 &&
                                (totalAfter < bestTotal - 1e-12 ||
                                    (Math.Abs(totalAfter - bestTotal) <= 1e-12 && frequency > bestFrequency)));
                        if (!better)
                            continue;

                        bestCodonIndex = index;
                        bestCodon = alternative;
                        bestWindow = windowAfter;
                        bestTotal = totalAfter;
                        bestFrequency = frequency;
                    }
                }

                if (bestCodon is null)
                {
                    // This window's bounds are unreachable with synonymous codons: leave it.
                    frozen = breaching + 1;
                    continue;
                }

                WriteCodon(bestCodonIndex, bestCodon);
            }
        }
    }

    // G+C count of a codon via the canonical counting primitive (SequenceExtensions).
    private static int CountGc(string codon) => codon.AsSpan().CountGcAndValidNucleotides().GcCount;

    #endregion

    #region CAI Calculation

    /// <summary>
    /// Calculates the Codon Adaptation Index (CAI) of Sharp &amp; Li (1987) against a codon-usage
    /// frequency table: <c>CAI = exp((1/L) Σ ln w_k)</c> with <c>w_ij = f_ij / max_j f_ij</c> over the
    /// synonymous codons of amino acid i. Delegates to the canonical
    /// <see cref="CodonUsageAnalyzer.CalculateCai(string, IReadOnlyDictionary{string, double}, GeneticCode)"/>
    /// core (CodonW <c>cai_out</c> conventions), under the genetic code <paramref name="geneticCode"/>
    /// (default Standard): the value equals
    /// <c>CodonUsageAnalyzer.CalculateCai(codingSequence, table, geneticCode)</c>, so stops and
    /// single-codon families are those of that code (under table 2, for instance, AGA/AGG are
    /// stops, the Met family is {AUA, AUG} and the Trp family {UGA, UGG}, so neither is excluded).
    /// </summary>
    /// <remarks>
    /// Stop codons are never scored. A relative adaptiveness below 0.0001 (codon absent from the
    /// table while a synonym is present) is replaced by 0.01 (CodonW; Bulmer 1988). An amino acid
    /// with no frequency data in <paramref name="table"/> is not scored. Triplets with symbols other
    /// than A/C/G/T/U (any case) are skipped without shifting the frame; a trailing partial codon is
    /// ignored. Returns 0 when no codon is scored.
    /// </remarks>
    /// <param name="codingSequence">Coding sequence (DNA or RNA; case-insensitive).</param>
    /// <param name="table">Reference codon usage table (frequencies keyed by RNA or DNA codon).</param>
    /// <param name="excludeSingleCodonAminoAcids">
    /// <see langword="true"/> (default): codons of single-codon amino acids (Met/AUG, Trp/UGG) are
    /// excluded, as Sharp &amp; Li (1987) prescribe ("codon families containing a single codon … should
    /// be excluded", quoted by Xia 2007, Evol. Bioinform. 3:53-58) and CodonW, seqinr and Biopython
    /// implement. <see langword="false"/>: they are scored with w = 1 (EMBOSS <c>cai</c> convention),
    /// which inflates CAI of Met/Trp-rich genes.
    /// </param>
    /// <param name="geneticCode">
    /// NCBI genetic code defining the stop codons and the synonymous families
    /// (<see langword="null"/>, the default: <see cref="GeneticCode.Standard"/>, identical to the
    /// historical value).
    /// </param>
    /// <exception cref="ArgumentException"><paramref name="table"/> has no frequency dictionary.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A frequency is negative or not finite.</exception>
    public static double CalculateCAI(
        string codingSequence,
        CodonUsageTable table,
        bool excludeSingleCodonAminoAcids = true,
        GeneticCode? geneticCode = null)
    {
        if (string.IsNullOrEmpty(codingSequence))
            return 0;
        if (table.CodonFrequencies is null)
            throw new ArgumentException("Codon usage table has no frequencies.", nameof(table));

        // The canonical core keys codons in DNA spelling.
        var reference = new Dictionary<string, double>(table.CodonFrequencies.Count);
        foreach (var (codon, frequency) in table.CodonFrequencies)
            reference[ToUpperAscii(codon).Replace('U', 'T')] = frequency;

        return CodonUsageAnalyzer.CalculateCai(
            codingSequence, reference, geneticCode ?? GeneticCode.Standard,
            excludeSingleCodonFamilies: excludeSingleCodonAminoAcids);
    }

    #endregion

    #region Sequence Modification

    /// <summary>
    /// Removes restriction-enzyme recognition sites from a coding sequence with synonymous codon
    /// substitutions, preserving the encoded protein.
    /// </summary>
    /// <remarks>
    /// Each site is matched with IUPAC ambiguity semantics (<see cref="IupacHelper.MatchesIupac"/>,
    /// the same matcher <see cref="RestrictionAnalyzer"/> uses), on <b>both strands</b>: a
    /// restriction enzyme cuts double-stranded DNA, so a non-palindromic site such as BsaI
    /// GGTCTC is also eliminated where its reverse complement GAGACC occurs (REBASE lists one
    /// strand only). For every occurrence, all synonymous replacements of the codons overlapping
    /// it are considered and the one that removes the occurrence with the <b>highest usage
    /// frequency</b> in <paramref name="table"/> is applied (DNA Chisel resolves an
    /// <c>AvoidPattern</c> constraint against the codon-optimisation objective the same way);
    /// ties resolve to the leftmost codon and NCBI codon order, so the result is deterministic.
    /// Exactly one codon is changed per removed occurrence. An occurrence that no single
    /// synonymous substitution can remove (e.g. inside a run of Met/Trp codons) is left in place
    /// and does not stop the other occurrences from being processed. The sequence is upper-cased,
    /// read as RNA (T → U) and trimmed to complete codons.
    /// </remarks>
    /// <param name="codingSequence">Coding sequence (DNA or RNA); empty → empty result.</param>
    /// <param name="restrictionSites">Recognition sequences (DNA or RNA, IUPAC codes allowed).</param>
    /// <param name="table">Codon-usage table of the expression host, used to pick replacements.</param>
    /// <param name="geneticCode">
    /// NCBI genetic code whose synonymous families the substitutions come from, so the protein is
    /// preserved under that code (<see langword="null"/>, the default:
    /// <see cref="GeneticCode.Standard"/>, identical to the historical output).
    /// </param>
    public static string RemoveRestrictionSites(
        string codingSequence,
        IEnumerable<string> restrictionSites,
        CodonUsageTable table,
        GeneticCode? geneticCode = null)
    {
        if (string.IsNullOrEmpty(codingSequence))
            return "";

        ArgumentNullException.ThrowIfNull(restrictionSites);

        string rna = ToUpperRna(codingSequence);
        var codons = SplitIntoCodons(rna);
        var families = FamiliesFor(geneticCode);

        // Match on DNA spelling: 'U' is not an IUPAC DNA code.
        var patterns = new List<string>();
        foreach (var site in restrictionSites)
        {
            if (string.IsNullOrEmpty(site))
                continue;
            string pattern = ToUpperAscii(site).Replace('U', 'T');
            if (!pattern.All(IupacHelper.IsNucleotideCode))
                throw new ArgumentException(
                    $"Restriction site '{site}' contains a character that is not an IUPAC nucleotide code.",
                    nameof(restrictionSites));

            AddPattern(patterns, pattern);
            AddPattern(patterns, DnaSequence.GetReverseComplementString(pattern));
        }

        foreach (string pattern in patterns)
        {
            int from = 0;
            // A substitution can create a new occurrence elsewhere, whose removal could in
            // principle undo the first one; cap the rewrites per pattern so the loop always
            // terminates (one rewrite per codon plus a small margin).
            int budget = codons.Count + 8;
            while (budget-- > 0)
            {
                string dna = string.Concat(codons).Replace('U', 'T');
                int position = FindPattern(dna, pattern, from);
                if (position < 0)
                    break;

                if (TryRemoveOccurrence(codons, position, pattern, table, families))
                    continue; // re-scan from the same offset: the edit may expose/shift matches

                from = position + 1; // unremovable occurrence — keep it and look further along
            }
        }

        return string.Concat(codons);

        static void AddPattern(List<string> patterns, string pattern)
        {
            if (!patterns.Contains(pattern))
                patterns.Add(pattern);
        }
    }

    // First index >= from where the IUPAC pattern matches (-1 when absent).
    private static int FindPattern(string dnaSequence, string pattern, int from)
    {
        for (int i = Math.Max(0, from); i + pattern.Length <= dnaSequence.Length; i++)
        {
            bool match = true;
            for (int k = 0; k < pattern.Length && match; k++)
                match = IupacHelper.MatchesIupac(dnaSequence[i + k], pattern[k]);
            if (match)
                return i;
        }

        return -1;
    }

    // Replaces the single synonymous codon that removes the occurrence at <paramref name="position"/>
    // with the smallest loss of codon usage; returns false when no synonymous substitution works.
    private static bool TryRemoveOccurrence(
        List<string> codons, int position, string pattern, CodonUsageTable table, CodeFamilies families)
    {
        int firstCodon = position / 3;
        int lastCodon = Math.Min((position + pattern.Length - 1) / 3, codons.Count - 1);

        int bestIndex = -1;
        string? bestCodon = null;
        double bestFrequency = double.NegativeInfinity;

        for (int i = firstCodon; i <= lastCodon; i++)
        {
            string original = codons[i];
            if (!TryTranslateCodon(original, families, out char aa) || aa == '*')
                continue;

            foreach (string alternative in SynonymsFor(original, families))
            {
                if (alternative == original)
                    continue;

                string dna = ConcatCodons(codons, 0, codons.Count, i, alternative).Replace('U', 'T');
                if (FindPattern(dna, pattern, position) == position)
                    continue;

                double frequency = table.CodonFrequencies.GetValueOrDefault(alternative, 0);
                if (frequency > bestFrequency)
                {
                    bestFrequency = frequency;
                    bestIndex = i;
                    bestCodon = alternative;
                }
            }
        }

        if (bestCodon is null)
            return false;

        codons[bestIndex] = bestCodon;
        return true;
    }

    /// <summary>
    /// How <see cref="ReduceSecondaryStructure(string, CodonUsageTable, int, double, double, SecondaryStructureScoring, double, int)"/>
    /// scores the structure of a window.
    /// </summary>
    public enum SecondaryStructureScoring
    {
        /// <summary>
        /// Fraction of position pairs (i, j), j ≥ i + 4, that can form a canonical pair
        /// (Watson-Crick or G·U wobble). A pairing-potential heuristic, not an energy model;
        /// the historical default.
        /// </summary>
        PairCount,

        /// <summary>
        /// Minimum free energy of the window in kcal/mol, from the canonical Turner 2004
        /// nearest-neighbour folding engine
        /// <see cref="RnaSecondaryStructure.CalculateMinimumFreeEnergy(string, int)"/>
        /// (Zuker-style DP, NNDB Turner 2004 parameters; the energy model of ViennaRNA RNAfold).
        /// Minimising structure means maximising the MFE (making it less negative).
        /// </summary>
        Turner2004Mfe
    }

    /// <summary>
    /// Reduces mRNA secondary structure by replacing codons inside structured windows
    /// with synonymous codons that lower the window's structure score.
    /// </summary>
    /// <remarks>
    /// Windows of <paramref name="windowSize"/> nucleotides (rounded up to whole codons) are
    /// walked left to right; a window whose score exceeds the mode's threshold is rewritten codon
    /// by codon, each codon taking the synonymous codon that minimises the score of the
    /// <b>current</b> window content (the baseline is re-evaluated after every accepted change),
    /// among codons whose usage frequency in <paramref name="table"/> is at least
    /// <paramref name="minCodonFrequency"/>. The sequence is upper-cased, read as RNA and trimmed
    /// to complete codons; a sequence shorter than <paramref name="windowSize"/> is returned
    /// normalised but otherwise unchanged. The translated protein is preserved in every mode:
    /// only synonymous codons of the genetic code <paramref name="geneticCode"/> (default
    /// <see cref="GeneticCode.Standard"/>) are used (a stop codon can only become another stop
    /// codon).
    /// <para><b>Scoring modes</b> (<paramref name="scoring"/>):</para>
    /// <list type="bullet">
    /// <item><description><see cref="SecondaryStructureScoring.PairCount"/> (default, unchanged
    /// behaviour) — the score is the fraction of position pairs (i, j), j ≥ i + 4, whose bases can
    /// form a canonical pair, Watson-Crick A·U / G·C <b>or the G·U wobble</b>, taken from the
    /// canonical <see cref="RnaSecondaryStructure.CanPair"/> (ViennaRNA default pair set), not a
    /// private Watson-Crick-only copy. A window is rewritten when that fraction is strictly above
    /// <paramref name="structureThreshold"/>. This is a pairing-potential heuristic: it ignores
    /// stacking, loop penalties and pair nesting.</description></item>
    /// <item><description><see cref="SecondaryStructureScoring.Turner2004Mfe"/> (opt-in) — the
    /// score is the window's minimum free energy from
    /// <see cref="RnaSecondaryStructure.CalculateMinimumFreeEnergy(string, int)"/> (Turner 2004
    /// nearest-neighbour parameters, the energy model of ViennaRNA RNAfold), and the search
    /// maximises it: a window is rewritten when its MFE is strictly below
    /// <paramref name="mfeThreshold"/> kcal/mol, and a synonymous codon is accepted only when it
    /// raises (makes less negative) the window MFE. <paramref name="structureThreshold"/> is
    /// ignored in this mode. Folding is O(w³) per evaluation, so window energies are memoised by
    /// window content and the pass stops once <paramref name="maxFoldEvaluations"/> distinct folds
    /// have been computed (deterministically, leaving the remaining windows as they are), which
    /// bounds the runtime of a long CDS. Note that DNA Chisel has no free-energy objective (its
    /// <c>AvoidHairpins</c> is the IDT reverse-complement stem rule), so the energy model is taken
    /// from ViennaRNA RNAfold / NNDB Turner 2004 through this repository's own folding
    /// engine.</description></item>
    /// </list>
    /// </remarks>
    /// <param name="codingSequence">Coding sequence (DNA or RNA).</param>
    /// <param name="table">Codon-usage table of the expression host.</param>
    /// <param name="windowSize">Window width in nucleotides (default 40).</param>
    /// <param name="structureThreshold">Pair-count score above which a window is rewritten (default 0.5; <see cref="SecondaryStructureScoring.PairCount"/> only).</param>
    /// <param name="minCodonFrequency">Minimum usage frequency of a replacement codon (default 0.1).</param>
    /// <param name="scoring">Window scoring model (default <see cref="SecondaryStructureScoring.PairCount"/>).</param>
    /// <param name="mfeThreshold">MFE in kcal/mol below which a window is rewritten (default 0, i.e. every window with a thermodynamically favourable structure; <see cref="SecondaryStructureScoring.Turner2004Mfe"/> only).</param>
    /// <param name="maxFoldEvaluations">Budget of distinct window folds in <see cref="SecondaryStructureScoring.Turner2004Mfe"/> mode (default 20000); the pass stops when it is exhausted.</param>
    /// <param name="geneticCode">
    /// NCBI genetic code whose synonymous families the rewrites use (<see langword="null"/>, the
    /// default: <see cref="GeneticCode.Standard"/>, identical to the historical output).
    /// </param>
    public static string ReduceSecondaryStructure(
        string codingSequence,
        CodonUsageTable table,
        int windowSize = 40,
        double structureThreshold = 0.5,
        double minCodonFrequency = 0.1,
        SecondaryStructureScoring scoring = SecondaryStructureScoring.PairCount,
        double mfeThreshold = 0.0,
        int maxFoldEvaluations = 20_000,
        GeneticCode? geneticCode = null)
    {
        if (string.IsNullOrEmpty(codingSequence))
            return codingSequence;

        string rna = ToUpperRna(codingSequence);
        if (rna.Length < windowSize)
            return rna;

        var codons = SplitIntoCodons(rna);
        var families = FamiliesFor(geneticCode);
        int windowCodons = windowSize / 3 + 1;

        bool turner = scoring == SecondaryStructureScoring.Turner2004Mfe;
        // Window energies are memoised by window content: the candidate windows of overlapping
        // positions repeat heavily, and a fold is orders of magnitude costlier than a lookup.
        var mfeCache = turner ? new Dictionary<string, double>(StringComparer.Ordinal) : null;
        int folds = 0;

        // Higher = more structured in both modes, so the same minimisation drives both:
        // the pair-count fraction, or the negated MFE (maximising MFE = minimising −MFE).
        double Score(string window)
        {
            if (!turner)
                return CalculateLocalStructure(window);

            if (!mfeCache!.TryGetValue(window, out double mfe))
            {
                mfe = RnaSecondaryStructure.CalculateMinimumFreeEnergy(window);
                mfeCache[window] = mfe;
                folds++;
            }

            return -mfe;
        }

        double threshold = turner ? -mfeThreshold : structureThreshold;

        for (int i = 0; i + windowCodons <= codons.Count; i++)
        {
            if (turner && folds >= maxFoldEvaluations)
                break; // fold budget exhausted — remaining windows are left unchanged

            if (Score(ConcatCodons(codons, i, windowCodons, -1, "")) <= threshold)
                continue;

            for (int j = i; j < i + windowCodons; j++)
            {
                string current = codons[j];
                double bestScore = Score(ConcatCodons(codons, i, windowCodons, -1, current));
                string bestAlternative = current;

                foreach (string alternative in SynonymsFor(current, families))
                {
                    if (alternative == current ||
                        table.CodonFrequencies.GetValueOrDefault(alternative, 0) < minCodonFrequency)
                    {
                        continue;
                    }

                    double score = Score(ConcatCodons(codons, i, windowCodons, j, alternative));

                    if (score < bestScore)
                    {
                        bestScore = score;
                        bestAlternative = alternative;
                    }
                }

                codons[j] = bestAlternative;
            }
        }

        return string.Concat(codons);
    }

    private static double CalculateLocalStructure(string sequence)
    {
        int complementaryPairs = 0;
        int n = sequence.Length;

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 4; j < n; j++)
            {
                if (RnaSecondaryStructure.CanPair(sequence[i], sequence[j]))
                    complementaryPairs++;
            }
        }

        double maxPairs = (n * (n - 4)) / 2.0;
        return maxPairs > 0 ? complementaryPairs / maxPairs : 0;
    }

    #endregion

    #region Analysis Functions

    /// <summary>
    /// Reports every in-frame codon whose usage frequency in <paramref name="table"/> is strictly
    /// below <paramref name="threshold"/> (a "rare" codon for the target organism).
    /// </summary>
    /// <remarks>
    /// The sequence is read in frame 0 (case-insensitive, DNA T or RNA U). Only the 64 unambiguous
    /// codons over {A,C,G,U} are screened: a triplet containing an IUPAC ambiguity code (N, R, Y, …)
    /// or any other symbol has no codon-usage frequency and is skipped without shifting the frame,
    /// and a trailing partial triplet is ignored — the same codon set as the canonical counter
    /// <see cref="CodonUsageAnalyzer.CountCodons(string)"/> (EMBOSS <c>ajCodSetTripletsS</c>).
    /// A valid codon that is absent from <paramref name="table"/> (never observed in the reference
    /// genes) has frequency 0 and is reported whenever <paramref name="threshold"/> &gt; 0.
    /// Stop codons are screened like any other row of the table. The amino acid is the
    /// translation under <paramref name="geneticCode"/> (default <see cref="GeneticCode.Standard"/>).
    /// </remarks>
    /// <returns>
    /// (0-based nucleotide position of the codon, RNA codon, one-letter amino acid or <c>*</c>,
    /// table frequency) in sequence order.
    /// </returns>
    /// <param name="geneticCode">
    /// NCBI genetic code used to translate each reported codon (<see langword="null"/>, the
    /// default: <see cref="GeneticCode.Standard"/>, identical to the historical output). The
    /// selection of rare codons itself is code-independent (it is a per-codon table lookup); the
    /// code only names the amino acid, so e.g. UGA is reported as <c>*</c> under table 1 and as
    /// <c>W</c> under table 2.
    /// </param>
    public static IEnumerable<(int Position, string Codon, string AminoAcid, double Frequency)> FindRareCodons(
        string codingSequence,
        CodonUsageTable table,
        double threshold = 0.15,
        GeneticCode? geneticCode = null)
        => FindRareCodons(codingSequence, table, threshold, inclusiveThreshold: false, geneticCode: geneticCode);

    /// <summary>
    /// <see cref="FindRareCodons(string, CodonUsageTable, double)"/> with a selectable threshold
    /// comparison: <paramref name="inclusiveThreshold"/> = <see langword="true"/> also reports codons
    /// whose frequency EQUALS <paramref name="threshold"/> (frequency ≤ threshold — the Sherlocc
    /// "slow position" comparison, <c>if($avrg &lt;= $avrg_lim)</c> in sherlocc.pl, README "equal or
    /// lower than the specified threshold"; Chartier et&#160;al. 2012); <see langword="false"/> is
    /// the default strict &lt; contract of the 3-argument overload.
    /// </summary>
    /// <param name="geneticCode">
    /// NCBI genetic code used to translate each reported codon (<see langword="null"/>, the
    /// default: <see cref="GeneticCode.Standard"/>, identical to the historical output). The
    /// selection of rare codons itself is code-independent (it is a per-codon table lookup); the
    /// code only names the amino acid, so e.g. UGA is reported as <c>*</c> under table 1 and as
    /// <c>W</c> under table 2.
    /// </param>
    public static IEnumerable<(int Position, string Codon, string AminoAcid, double Frequency)> FindRareCodons(
        string codingSequence,
        CodonUsageTable table,
        double threshold,
        bool inclusiveThreshold,
        GeneticCode? geneticCode = null)
    {
        var codons = SplitIntoScreenedCodons(codingSequence);
        var code = geneticCode ?? GeneticCode.Standard;

        for (int i = 0; i < codons.Length; i++)
        {
            string? codon = codons[i];
            if (codon is null)
                continue;

            double freq = table.CodonFrequencies.GetValueOrDefault(codon, 0);
            if (IsBelowThreshold(freq, threshold, inclusiveThreshold))
                yield return (i * 3, codon, code.Translate(codon).ToString(), freq);
        }
    }

    // Rare-codon comparison: strict (freq < threshold) or Sherlocc-inclusive (freq <= threshold;
    // sherlocc.pl: if($avrg <= $avrg_lim){ #We have a pause).
    private static bool IsBelowThreshold(double frequency, double threshold, bool inclusive) =>
        inclusive ? frequency <= threshold : frequency < threshold;

    // Default sliding-window width (in codons) for the %MinMax profile.
    // "%MinMax results are typically averaged over an 18-codon sliding window."
    // Clarke TF, Clark PL (2008) "Rare Codons Cluster", PLoS ONE 3(10):e3412.
    private const int DefaultMinMaxWindowCodons = 18;

    // Sherlocc rare-codon-cluster (RCC) detection parameters.
    // "a seven position-wide window ... containing at least four pause positions out of seven."
    // Chartier M, Gaudreault F, Najmanovich R (2012) Bioinformatics 28(11):1438-1445,
    // doi:10.1093/bioinformatics/bts149.
    private const int DefaultClusterWindowCodons = 7;
    private const int DefaultClusterMinRareCodons = 4;

    /// <summary>
    /// Computes the %MinMax codon-usage profile of Clarke &amp; Clark (2008) over a sliding window.
    /// </summary>
    /// <remarks>
    /// For each amino acid <c>i</c> with <c>n</c> synonymous codons, let <c>Xij</c> be the usage
    /// frequency of the codon actually used, <c>Xmax,i</c> / <c>Xmin,i</c> the usage frequencies of
    /// the most / least common synonymous codon, and <c>Xavg,i</c> the arithmetic mean of the
    /// synonymous codon frequencies. Over a window of <paramref name="windowSize"/> codons:
    /// if Σ Xij &gt; Σ Xavg,i the window yields %Max = Σ(Xij − Xavg,i) / Σ(Xmax,i − Xavg,i) × 100
    /// (returned as a positive value); if Σ Xij &lt; Σ Xavg,i it yields %Min =
    /// Σ(Xavg,i − Xij) / Σ(Xavg,i − Xmin,i) × 100 (returned as a negative value). Source:
    /// Clarke &amp; Clark (2008), PLoS ONE 3(10):e3412; reproduced term-for-term from the Clark
    /// lab reference code (<c>calculateMinMax</c> in CHARMING.py, Wright et&#160;al. 2022), whose
    /// synonymous families are those of the Standard code including the stop family
    /// {UAA, UAG, UGA}. Families here come from <paramref name="geneticCode"/> (default
    /// <see cref="GeneticCode.Standard"/>, the reference behaviour); codons of
    /// single-codon amino acids (Met, Trp) contribute 0 to both numerator and denominator, and an
    /// ambiguous triplet (non-ACGU symbol) occupies its window position but contributes nothing.
    /// <para>
    /// <b>Frequency scale.</b> The reference implementation takes a codon usage table in
    /// frequency per thousand codons (Kazusa "/1000" column), so the window sums weight each
    /// residue by its overall usage. This overload therefore uses
    /// <see cref="CodonUsageTable.PerThousand"/> when the table carries it (the <see cref="Yeast"/>
    /// and <see cref="Human"/> presets: Kazusa 4932 / 9606 per-thousand usage, so the presets
    /// reproduce the reference), and otherwise <c>CodonFrequencies</c> — a caller table whose
    /// <c>CodonFrequencies</c> hold per-thousand values also reproduces the reference exactly. With
    /// per-amino-acid relative fractions (<see cref="EColiK12"/>, or
    /// <see cref="MinMaxUsageBasis.CodonFrequencies"/>) every residue has equal weight, which gives
    /// the reference value only for windows of a single amino acid.
    /// </para>
    /// </remarks>
    /// <param name="codingSequence">DNA or RNA coding sequence (T is normalised to U).</param>
    /// <param name="table">
    /// Reference codon-usage table: per-thousand usage for reference-identical values, or
    /// per-amino-acid relative fractions (see remarks).
    /// </param>
    /// <param name="windowSize">Sliding-window width in codons (default 18, per Clarke &amp; Clark 2008).</param>
    /// <returns>
    /// One <see cref="MinMaxWindow"/> per window position (codon indices
    /// 0 .. codonCount − windowSize). Empty if the sequence has fewer than
    /// <paramref name="windowSize"/> complete codons.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">windowSize &lt; 1.</exception>
    /// <param name="geneticCode">
    /// NCBI genetic code defining the synonymous families of the %MinMax sums
    /// (<see langword="null"/>, the default: <see cref="GeneticCode.Standard"/>, identical to the
    /// historical profile).
    /// </param>
    public static IReadOnlyList<MinMaxWindow> CalculateMinMaxProfile(
        string codingSequence,
        CodonUsageTable table,
        int windowSize = DefaultMinMaxWindowCodons,
        GeneticCode? geneticCode = null)
        => CalculateMinMaxProfile(
            codingSequence, table, windowSize, MinMaxUsageBasis.PerThousandWhenAvailable, geneticCode);

    /// <summary>
    /// <see cref="CalculateMinMaxProfile(string, CodonUsageTable, int)"/> with an explicit choice of
    /// the usage values: <see cref="MinMaxUsageBasis.PerThousandWhenAvailable"/> (default of the
    /// 3-argument overload) or <see cref="MinMaxUsageBasis.CodonFrequencies"/> (always the
    /// per-amino-acid fractions / caller values in <c>CodonFrequencies</c>).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">windowSize &lt; 1, or <paramref name="basis"/> undefined.</exception>
    /// <param name="geneticCode">
    /// NCBI genetic code defining the synonymous families of the %MinMax sums
    /// (<see langword="null"/>, the default: <see cref="GeneticCode.Standard"/>, identical to the
    /// historical profile).
    /// </param>
    public static IReadOnlyList<MinMaxWindow> CalculateMinMaxProfile(
        string codingSequence,
        CodonUsageTable table,
        int windowSize,
        MinMaxUsageBasis basis,
        GeneticCode? geneticCode = null)
    {
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1 codon.");
        IReadOnlyDictionary<string, double> usage = basis switch
        {
            MinMaxUsageBasis.PerThousandWhenAvailable => table.PerThousand ?? table.CodonFrequencies,
            MinMaxUsageBasis.CodonFrequencies => table.CodonFrequencies,
            _ => throw new ArgumentOutOfRangeException(nameof(basis), basis, "Unknown %MinMax usage basis."),
        };

        var families = FamiliesFor(geneticCode);
        var profile = new List<MinMaxWindow>();
        var codons = SplitIntoScreenedCodons(codingSequence);
        if (codons.Length < windowSize)
            return profile;

        // Per-codon (Xij), per-family average (Xavg), max (Xmax) and min (Xmin) frequencies.
        var xij = new double[codons.Length];
        var xavg = new double[codons.Length];
        var xmax = new double[codons.Length];
        var xmin = new double[codons.Length];
        var familyStats = new Dictionary<char, (double Avg, double Max, double Min)>();
        for (int i = 0; i < codons.Length; i++)
        {
            string? codon = codons[i];
            if (codon is null)
                continue; // ambiguous triplet: all four terms stay 0 (no contribution).

            xij[i] = usage.GetValueOrDefault(codon, 0);
            char aa = families.Code.Translate(codon);
            if (!familyStats.TryGetValue(aa, out var stats))
            {
                double sum = 0, max = double.MinValue, min = double.MaxValue;
                int n = 0;
                foreach (string syn in families.Synonyms.GetValueOrDefault(aa, Array.Empty<string>()))
                {
                    double f = usage.GetValueOrDefault(syn, 0);
                    sum += f;
                    n++;
                    if (f > max) max = f;
                    if (f < min) min = f;
                }
                stats = n > 0 ? (sum / n, max, min) : (0.0, 0.0, 0.0);
                familyStats[aa] = stats;
            }

            (xavg[i], xmax[i], xmin[i]) = stats;
        }

        for (int start = 0; start + windowSize <= codons.Length; start++)
        {
            double sumXij = 0, sumXavg = 0, sumMaxDelta = 0, sumMinDelta = 0;
            for (int k = start; k < start + windowSize; k++)
            {
                sumXij += xij[k];
                sumXavg += xavg[k];
                sumMaxDelta += xmax[k] - xavg[k];
                sumMinDelta += xavg[k] - xmin[k];
            }

            double percent;
            if (sumXij > sumXavg)
            {
                // %Max — positive value.
                percent = sumMaxDelta > 0 ? (sumXij - sumXavg) / sumMaxDelta * 100.0 : 0.0;
            }
            else if (sumXij < sumXavg)
            {
                // %Min — returned as a negative value (rare-codon side).
                percent = sumMinDelta > 0 ? -((sumXavg - sumXij) / sumMinDelta * 100.0) : 0.0;
            }
            else
            {
                percent = 0.0;
            }

            profile.Add(new MinMaxWindow(start, percent));
        }

        return profile;
    }

    /// <summary>
    /// Detects rare-codon clusters (RCCs) using the Sherlocc rule of Chartier et&#160;al. (2012):
    /// a window of <paramref name="windowSize"/> codons is a cluster when it contains at least
    /// <paramref name="minRareCodons"/> rare ("pause") codons.
    /// </summary>
    /// <remarks>
    /// A codon is "rare"/"pause" when its usage frequency in <paramref name="table"/> is equal to
    /// or lower than <paramref name="rareThreshold"/> — the Sherlocc comparison: sherlocc.pl
    /// <c>if($avrg &lt;= $avrg_lim){ #We have a pause</c>, README "positions with a codon usage
    /// frequency average equal or lower than the specified threshold … considered as 'slow'
    /// positions" (github.com/mtthchrtr/sherlocc, commit b3b4d15). Until review 2026-10 the
    /// comparison was strict &lt;; that behaviour is available through the overload with
    /// <c>inclusiveThreshold: false</c>. The unambiguous-codon screen is the one of
    /// <see cref="FindRareCodons(string, CodonUsageTable, double)"/> (an ambiguous triplet is never
    /// a pause). Simplification kept: Sherlocc compares a centred 7-codon window AVERAGE of per-thousand
    /// usage (averaged over every sequence of a protein-family alignment) with the threshold, whereas
    /// this method compares each codon's own table value. Overlapping windows are merged into maximal clusters so a long
    /// rare run is reported once. This is opt-in; <see cref="FindRareCodons"/> (per-codon) is
    /// unchanged. The rule is <b>genetic-code independent</b>: a pause is a per-codon usage-table
    /// lookup, no codon is translated and no synonymous family is formed, so no
    /// <see cref="GeneticCode"/> parameter is offered (the same holds for
    /// <see cref="FindConservedRareCodonClusters(IReadOnlyList{string}, CodonUsageTable, double, SherloccOptions?)"/>,
    /// <see cref="CalculateCodonUsage(string)"/> and <see cref="CompareCodonUsage(string, string)"/>). Defaults reproduce the published Sherlocc rule "a seven position-wide window …
    /// containing at least four pause positions out of seven" (Chartier et&#160;al. 2012,
    /// doi:10.1093/bioinformatics/bts149).
    /// </remarks>
    /// <param name="codingSequence">DNA or RNA coding sequence (T is normalised to U).</param>
    /// <param name="table">Reference codon-usage table (per-amino-acid relative fractions).</param>
    /// <param name="rareThreshold">Per-codon rare-frequency cutoff (default 0.15, frequency ≤ cutoff, per Sherlocc).</param>
    /// <param name="windowSize">Cluster window width in codons (default 7, per Sherlocc).</param>
    /// <param name="minRareCodons">Minimum rare codons in a window to call a cluster (default 4, per Sherlocc).</param>
    /// <returns>Maximal, non-overlapping <see cref="RareCodonCluster"/> regions in codon-index order.</returns>
    /// <exception cref="ArgumentOutOfRangeException">windowSize &lt; 1 or minRareCodons &lt; 1.</exception>
    public static IReadOnlyList<RareCodonCluster> FindRareCodonClusters(
        string codingSequence,
        CodonUsageTable table,
        double rareThreshold = 0.15,
        int windowSize = DefaultClusterWindowCodons,
        int minRareCodons = DefaultClusterMinRareCodons)
        => FindRareCodonClusters(codingSequence, table, rareThreshold, windowSize, minRareCodons, inclusiveThreshold: true);

    /// <summary>
    /// <see cref="FindRareCodonClusters(string, CodonUsageTable, double, int, int)"/> with a selectable
    /// pause comparison: <paramref name="inclusiveThreshold"/> = <see langword="true"/> (the default of
    /// the 5-argument overload, Sherlocc) marks frequency ≤ threshold as a pause;
    /// <see langword="false"/> marks only frequency &lt; threshold (the pre-2026-10 behaviour, the
    /// same comparison as <see cref="FindRareCodons(string, CodonUsageTable, double)"/>).
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException">windowSize &lt; 1 or minRareCodons &lt; 1.</exception>
    public static IReadOnlyList<RareCodonCluster> FindRareCodonClusters(
        string codingSequence,
        CodonUsageTable table,
        double rareThreshold,
        int windowSize,
        int minRareCodons,
        bool inclusiveThreshold)
    {
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1 codon.");
        if (minRareCodons < 1)
            throw new ArgumentOutOfRangeException(nameof(minRareCodons), minRareCodons, "Minimum rare codons must be at least 1.");

        var clusters = new List<RareCodonCluster>();
        var codons = SplitIntoScreenedCodons(codingSequence);
        if (codons.Length < windowSize)
            return clusters;

        // Mark each codon as rare (pause) when its table frequency is at or below (Sherlocc) / strictly
        // below the threshold. An ambiguous triplet has no usage frequency and is never a pause position.
        var isRare = new bool[codons.Length];
        for (int i = 0; i < codons.Length; i++)
            isRare[i] = codons[i] is { } codon
                && IsBelowThreshold(table.CodonFrequencies.GetValueOrDefault(codon, 0), rareThreshold, inclusiveThreshold);

        int? mergedStart = null;
        int mergedEnd = -1;
        int windowRare = 0;
        for (int start = 0; start + windowSize <= codons.Length; start++)
        {
            if (start == 0)
            {
                for (int k = 0; k < windowSize; k++)
                    if (isRare[k]) windowRare++;
            }
            else
            {
                if (isRare[start - 1]) windowRare--;
                if (isRare[start + windowSize - 1]) windowRare++;
            }

            if (windowRare >= minRareCodons)
            {
                int end = start + windowSize - 1;
                if (mergedStart is null)
                {
                    mergedStart = start;
                    mergedEnd = end;
                }
                else if (start <= mergedEnd + 1)
                {
                    // Overlapping or adjacent qualifying window — extend the current cluster.
                    mergedEnd = end;
                }
                else
                {
                    clusters.Add(BuildCluster(mergedStart.Value, mergedEnd, isRare));
                    mergedStart = start;
                    mergedEnd = end;
                }
            }
        }

        if (mergedStart is not null)
            clusters.Add(BuildCluster(mergedStart.Value, mergedEnd, isRare));

        return clusters;
    }

    private static RareCodonCluster BuildCluster(int start, int end, bool[] isRare)
    {
        int count = 0;
        for (int k = start; k <= end; k++)
            if (isRare[k]) count++;
        return new RareCodonCluster(start, end, count);
    }

    /// <summary>
    /// Sherlocc 1.0 conserved rare-codon cluster analysis (Chartier, Gaudreault &amp; Najmanovich 2012;
    /// <c>sherlocc.pl</c> stage 3, github.com/mtthchrtr/sherlocc commit b3b4d15) of a codon alignment of
    /// homologous coding sequences, every sequence scored with the same usage table.
    /// </summary>
    /// <remarks>
    /// <para>Input convention: equal-length aligned nucleotide CDS (DNA or RNA, any case), one codon per
    /// 3 alignment characters; a codon column made only of <c>-</c>/<c>.</c> is a gap (Sherlocc <c>x</c>).
    /// A triplet with no usage value (partially gapped, ambiguous, or absent from the table) is also
    /// an <c>x</c> cell.</para>
    /// <para>Algorithm (verbatim from the script): (1) every cell gets its codon's per-thousand usage
    /// (stage 2, l. 494–522); (2) for every column i the values of columns i−3 … i+3 of ALL sequences
    /// are summed, skipping out-of-range and <c>x</c> cells (l. 780–798); (3) the column is
    /// representative only when the number of valued cells ≥ 0.75 × sequences × 7 (l. 772, 799) — the
    /// conservation/gap filter; (4) a representative column with average ≤ threshold is a pause
    /// (l. 803); (5) every window of columns i … i+6 holding ≥ 4 pauses yields a cluster from its first
    /// to its last pause (l. 838–854); (6) successive clusters are merged when the next one starts
    /// inside the current one, otherwise a new cluster starts — abutting clusters are NOT merged
    /// (l. 894–961); (7) each cluster's mean usage, middle point, the total length and the fraction of
    /// the alignment are those of the result table (l. 875–879, 1131–1155).</para>
    /// <para>Not ported: the FoldIndex/Uversky unstructured-region prediction (a separate output of the
    /// script) and stages 1–2 retrieval (UniProt/ENA web queries). The cluster P-values of the paper are
    /// not computed by <c>sherlocc.pl</c> and their model is not available (see the CODON-RARE-001
    /// Evidence), so no P-value is reported.</para>
    /// </remarks>
    /// <param name="alignedCodingSequences">Aligned in-frame coding sequences (gaps <c>-</c> or <c>.</c>).</param>
    /// <param name="table">Codon-usage table applied to every sequence (needs <see cref="CodonUsageTable.PerThousand"/> unless <see cref="SherloccOptions.UseCodonFrequencies"/>).</param>
    /// <param name="threshold">Sherlocc threshold (per thousand by default; README examples 14, 15, 16); pause = average ≤ threshold.</param>
    /// <param name="options">Window/filter parameters; <see langword="null"/> = Sherlocc defaults.</param>
    /// <exception cref="ArgumentNullException">A list, sequence or table is null.</exception>
    /// <exception cref="ArgumentException">Empty alignment, unequal lengths, length not a multiple of 3, or a table without per-thousand usage.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Threshold is NaN, or an option is out of range.</exception>
    public static SherloccResult FindConservedRareCodonClusters(
        IReadOnlyList<string> alignedCodingSequences,
        CodonUsageTable table,
        double threshold,
        SherloccOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignedCodingSequences);
        var tables = new CodonUsageTable[alignedCodingSequences.Count];
        Array.Fill(tables, table);
        return FindConservedRareCodonClusters(alignedCodingSequences, tables, threshold, options);
    }

    /// <summary>
    /// <see cref="FindConservedRareCodonClusters(IReadOnlyList{string}, CodonUsageTable, double, SherloccOptions?)"/>
    /// with one usage table per sequence — Sherlocc scores each family member with the codon usage of
    /// its own organism (stage 2 looks up <c>$cu_hash{$org_id}</c> per protein, l. 497).
    /// </summary>
    /// <param name="alignedCodingSequences">Aligned in-frame coding sequences (gaps <c>-</c> or <c>.</c>).</param>
    /// <param name="tablesPerSequence">One table per sequence, same order.</param>
    /// <param name="threshold">Sherlocc threshold; pause = average ≤ threshold.</param>
    /// <param name="options">Window/filter parameters; <see langword="null"/> = Sherlocc defaults.</param>
    public static SherloccResult FindConservedRareCodonClusters(
        IReadOnlyList<string> alignedCodingSequences,
        IReadOnlyList<CodonUsageTable> tablesPerSequence,
        double threshold,
        SherloccOptions? options = null)
    {
        ArgumentNullException.ThrowIfNull(alignedCodingSequences);
        ArgumentNullException.ThrowIfNull(tablesPerSequence);
        options ??= new SherloccOptions();
        if (double.IsNaN(threshold))
            throw new ArgumentOutOfRangeException(nameof(threshold), threshold, "Threshold must be a number.");
        if (options.AveragingWindow < 1 || options.AveragingWindow % 2 == 0)
            throw new ArgumentOutOfRangeException(nameof(options), options.AveragingWindow, "AveragingWindow must be a positive odd number of columns.");
        if (options.ClusterWindow < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.ClusterWindow, "ClusterWindow must be at least 1 column.");
        if (options.MinPauses < 1)
            throw new ArgumentOutOfRangeException(nameof(options), options.MinPauses, "MinPauses must be at least 1.");
        if (!(options.MinCoverage >= 0 && options.MinCoverage <= 1))
            throw new ArgumentOutOfRangeException(nameof(options), options.MinCoverage, "MinCoverage must be in [0, 1].");

        int rows = alignedCodingSequences.Count;
        if (rows == 0)
            throw new ArgumentException("The alignment must contain at least one sequence.", nameof(alignedCodingSequences));
        if (tablesPerSequence.Count != rows)
            throw new ArgumentException("One codon-usage table per sequence is required.", nameof(tablesPerSequence));

        int length = -1;
        for (int r = 0; r < rows; r++)
        {
            string? s = alignedCodingSequences[r];
            if (s is null)
                throw new ArgumentNullException(nameof(alignedCodingSequences), $"Sequence {r} is null.");
            if (s.Length % 3 != 0)
                throw new ArgumentException($"Sequence {r} length {s.Length} is not a multiple of 3.", nameof(alignedCodingSequences));
            if (length >= 0 && s.Length != length)
                throw new ArgumentException("All aligned sequences must have the same length.", nameof(alignedCodingSequences));
            length = s.Length;
            var t = tablesPerSequence[r];
            if (t.CodonFrequencies is null)
                throw new ArgumentNullException(nameof(tablesPerSequence), $"Table {r} has no codon frequencies.");
            if (!options.UseCodonFrequencies && t.PerThousand is null)
                throw new ArgumentException(
                    $"Table {r} ('{t.OrganismName}') has no per-thousand usage; Sherlocc uses CUTG per-thousand values " +
                    "(set SherloccOptions.UseCodonFrequencies to score CodonFrequencies instead).", nameof(tablesPerSequence));
        }

        int columns = length / 3;

        // Stage 2: one usage value per cell, NaN for an 'x' cell.
        var matrix = new double[rows, columns];
        for (int r = 0; r < rows; r++)
        {
            var usage = options.UseCodonFrequencies ? tablesPerSequence[r].CodonFrequencies : tablesPerSequence[r].PerThousand!;
            string s = alignedCodingSequences[r];
            for (int c = 0; c < columns; c++)
            {
                string codon = s.Substring(c * 3, 3).ToUpperInvariant().Replace('T', 'U');
                matrix[r, c] = usage.TryGetValue(codon, out double v) ? v : double.NaN;
            }
        }

        // Stage 3: window-and-family averages, representativeness filter, pauses.
        int half = options.AveragingWindow / 2;
        double minCells = rows * options.AveragingWindow * options.MinCoverage;
        var averages = new double?[columns];
        var counts = new int[columns];
        var pauses = new List<int>();
        var isPause = new bool[columns];
        for (int i = 0; i < columns; i++)
        {
            double sum = 0;
            int count = 0;
            for (int r = 0; r < rows; r++)
            {
                for (int c = i - half; c <= i + half; c++)
                {
                    if (c < 0 || c >= columns || double.IsNaN(matrix[r, c]))
                        continue;
                    sum += matrix[r, c];
                    count++;
                }
            }

            counts[i] = count;
            if (count > 0 && count >= minCells)
            {
                double avg = sum / count;
                averages[i] = avg;
                if (avg <= threshold)
                {
                    isPause[i] = true;
                    pauses.Add(i);
                }
            }
        }

        // Cluster windows i … i+W−1 with ≥ MinPauses pauses → (first pause, last pause).
        var merged = new List<(int Start, int Stop)>();
        for (int i = 0; i + options.ClusterWindow <= columns; i++)
        {
            int first = -1, last = -1, n = 0;
            for (int j = i; j < i + options.ClusterWindow; j++)
            {
                if (!isPause[j]) continue;
                if (first < 0) first = j;
                last = j;
                n++;
            }

            if (n < options.MinPauses)
                continue;
            if (merged.Count > 0 && first >= merged[^1].Start && first <= merged[^1].Stop)
                merged[^1] = (merged[^1].Start, last);
            else if (merged.Count == 0 || first > merged[^1].Stop)
                merged.Add((first, last));
        }

        var clusters = new List<ConservedRareCodonCluster>(merged.Count);
        int totalLength = 0;
        double meanSum = 0;
        foreach (var (start, stop) in merged)
        {
            double sum = 0;
            int pauseCount = 0;
            for (int k = start; k <= stop; k++)
            {
                sum += averages[k] ?? 0.0;
                if (isPause[k]) pauseCount++;
            }

            int spanLength = stop - start + 1;
            double mean = sum / spanLength;
            int middle = (int)((stop + 1) - spanLength / 2.0);
            clusters.Add(new ConservedRareCodonCluster(start, stop, pauseCount, mean, middle));
            totalLength += spanLength;
            meanSum += mean;
        }

        return new SherloccResult(
            rows,
            columns,
            averages,
            counts,
            pauses,
            clusters,
            totalLength,
            columns > 0 ? (double)totalLength / columns : 0.0,
            clusters.Count > 0 ? meanSum / clusters.Count : null);
    }

    /// <summary>
    /// Counts the in-frame (frame 0) codons of a coding sequence, returning RNA-spelled keys
    /// (<c>AUG</c>, <c>GCU</c>, …) mapped to raw counts (the "Number" column of an EMBOSS
    /// <c>cusp</c> / Kazusa codon usage table).
    /// </summary>
    /// <remarks>
    /// Input is case-insensitive and may be DNA (T) or RNA (U). Only the 64 unambiguous
    /// codons over {A,C,G,U} are counted: a triplet containing an IUPAC ambiguity code
    /// (N, R, Y, …) or any other non-nucleotide character is skipped without shifting the
    /// frame, and an incomplete trailing triplet is ignored — the contract documented for
    /// EMBOSS <c>ajCodSetTripletsS</c> ("Skips triplets with ambiguity codes and any
    /// incomplete triplet at the end"), and consistent with Biopython
    /// <c>CodonAdaptationIndex</c>, whose count table is closed over the 64 ACGT codons.
    /// Stop codons are counted like any other codon (they are rows of the cusp/Kazusa table).
    /// Delegates to the canonical counter <see cref="CodonUsageAnalyzer.CountCodons(string)"/>.
    /// Codons that do not occur are absent from the dictionary (no zero entries).
    /// </remarks>
    /// <param name="codingSequence">In-frame coding sequence (DNA or RNA); null/empty → empty result.</param>
    public static Dictionary<string, int> CalculateCodonUsage(string codingSequence)
    {
        var usage = new Dictionary<string, int>();

        if (string.IsNullOrEmpty(codingSequence))
            return usage;

        // The canonical counter normalises case and U/T itself and reports DNA spelling;
        // this API reports RNA spelling.
        foreach (var (codon, count) in CodonUsageAnalyzer.CountCodons(codingSequence))
            usage[codon.Replace('T', 'U')] = count;

        return usage;
    }

    /// <summary>
    /// Compares the codon usage of two sequences as the total-variation-distance similarity of
    /// their codon frequency distributions: <c>1 − ½·Σ_c |f₁(c) − f₂(c)|</c>, where
    /// <c>f_i(c) = count_i(c) / Σ count_i</c> over the codons counted by
    /// <see cref="CalculateCodonUsage(string)"/> (unambiguous, in-frame, complete codons only).
    /// Returns a value in [0, 1]; returns 0 when either sequence has no countable codon.
    /// </summary>
    public static double CompareCodonUsage(string sequence1, string sequence2)
    {
        var usage1 = CalculateCodonUsage(sequence1);
        var usage2 = CalculateCodonUsage(sequence2);

        var allCodons = usage1.Keys.Union(usage2.Keys).ToList();
        if (allCodons.Count == 0)
            return 0;

        int total1 = usage1.Values.Sum();
        int total2 = usage2.Values.Sum();

        if (total1 == 0 || total2 == 0)
            return 0;

        double l1Distance = 0;
        foreach (var codon in allCodons)
        {
            double freq1 = usage1.GetValueOrDefault(codon, 0) / (double)total1;
            double freq2 = usage2.GetValueOrDefault(codon, 0) / (double)total2;
            l1Distance += Math.Abs(freq1 - freq2);
        }

        return 1 - (l1Distance / 2);
    }

    #endregion

    #region Utility Methods

    // Frame-0 RNA-spelled codons for the per-position screens (rare codons, %MinMax, clusters):
    // index k = codon k; an ambiguous triplet is null (skipped, frame preserved). Delegates the
    // splitting and ACGT screen to the canonical CodonUsageAnalyzer core.
    private static string?[] SplitIntoScreenedCodons(string? codingSequence)
    {
        var codons = CodonUsageAnalyzer.SplitInFrameCodons(codingSequence);
        for (int k = 0; k < codons.Length; k++)
            codons[k] = codons[k]?.Replace('T', 'U');
        return codons;
    }

    // Upper-case RNA spelling (T read as U) used by every rewriting API of this class.
    private static string ToUpperRna(string sequence) => ToUpperAscii(sequence).Replace('T', 'U');

    // ASCII-only case folding (TRANS-PROT-001), as Core's SequenceExtensions.ToUpperAscii (internal to
    // Core): ToUpperInvariant folds U+017F 'ſ' → 'S', an IUPAC code, so "GGſ" would read as GGS → Gly.
    private static char ToUpperAscii(char c) => char.IsAsciiLetterLower(c) ? (char)(c - 32) : c;

    private static string ToUpperAscii(string s) => string.Create(s.Length, s, static (dest, src) =>
    {
        for (int i = 0; i < src.Length; i++)
            dest[i] = ToUpperAscii(src[i]);
    });

    // Every complete frame-0 triplet, including ambiguous ones (the rewriting APIs must keep
    // them in place to preserve the reading frame and the sequence length); a trailing partial
    // triplet is dropped.
    private static List<string> SplitIntoCodons(string sequence)
    {
        var codons = new List<string>(sequence.Length / 3);
        for (int i = 0; i + 2 < sequence.Length; i += 3)
        {
            codons.Add(sequence.Substring(i, 3));
        }
        return codons;
    }

    // Concatenates codons [start, start+count) with the codon at <paramref name="replaceAt"/>
    // swapped for <paramref name="replacement"/> (replaceAt &lt; 0 → no replacement), without
    // mutating the list.
    private static string ConcatCodons(List<string> codons, int start, int count, int replaceAt, string replacement)
    {
        var builder = new StringBuilder(count * 3);
        for (int k = start; k < start + count; k++)
            builder.Append(k == replaceAt ? replacement : codons[k]);
        return builder.ToString();
    }

    private static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcFractionFast();

    /// <summary>
    /// Creates a <see cref="CodonUsageTable"/> from codon frequencies, filling
    /// <see cref="CodonUsageTable.CodonToAminoAcid"/> with the canonical
    /// <see cref="GeneticCode"/> mapping (RNA spelling) of <paramref name="geneticCode"/>, so a
    /// caller-supplied table carries the same codon → amino-acid assignment as the built-in
    /// presets. Frequency keys are normalised to upper-case RNA spelling.
    /// </summary>
    /// <param name="organismName">Name stored in the table.</param>
    /// <param name="codonFrequencies">Per-amino-acid relative codon frequencies (DNA or RNA keys).</param>
    /// <param name="geneticCode">
    /// NCBI genetic code of the <see cref="CodonUsageTable.CodonToAminoAcid"/> mapping
    /// (<see langword="null"/>, the default: <see cref="GeneticCode.Standard"/>, identical to the
    /// historical table).
    /// </param>
    public static CodonUsageTable CreateCodonUsageTable(
        string organismName,
        IReadOnlyDictionary<string, double> codonFrequencies,
        GeneticCode? geneticCode = null)
    {
        ArgumentNullException.ThrowIfNull(codonFrequencies);

        var frequencies = new Dictionary<string, double>(codonFrequencies.Count);
        foreach (var (codon, frequency) in codonFrequencies)
            frequencies[ToUpperRna(codon)] = frequency;

        return new CodonUsageTable(organismName, frequencies, FamiliesFor(geneticCode).CodonToAminoAcid);
    }

    /// <summary>
    /// Builds a codon-usage table from a reference gene set (one concatenated in-frame coding
    /// sequence, or a single gene), as the relative frequency of each codon within its
    /// synonymous family of the genetic code <paramref name="geneticCode"/> (default Standard).
    /// </summary>
    /// <remarks>
    /// Codons are counted by the canonical counter (frame 0, complete unambiguous triplets only;
    /// see <see cref="CalculateCodonUsage(string)"/>). <b>A codon that does not occur in the
    /// reference set is given a count of 0.5</b>, "following the description in the original
    /// paper" — Sharp &amp; Li (1987) NAR 15:1281-1295, as implemented by Biopython
    /// <c>Bio.SeqUtils.CodonAdaptationIndex</c> (1.88). Without that pseudo-count an unobserved
    /// codon would be missing from the table and scored with the CodonW zero substitute
    /// (w = 0.01) by <see cref="CalculateCAI"/>, penalising codons about which the reference set
    /// simply carries no information. All 64 codons are therefore present in the result, and the
    /// relative adaptiveness w = f / max f derived from it equals Biopython's index exactly.
    /// Stop codons form their own family (as in Biopython). The frequencies of each family sum
    /// to 1. Stop codons of the chosen code form their own family; in tables 27/28/31 a
    /// dual-coding codon counts in the family of the amino acid it encodes (Biopython
    /// <c>forward_table</c>, DNA Chisel back-translation table).
    /// </remarks>
    /// <param name="referenceSequence">Reference coding sequence(s), in frame 0 (DNA or RNA).</param>
    /// <param name="organismName">Name stored in the resulting table.</param>
    /// <param name="geneticCode">
    /// NCBI genetic code defining the synonymous families the frequencies are normalised within,
    /// and the resulting table's codon → amino-acid mapping (<see langword="null"/>, the default:
    /// <see cref="GeneticCode.Standard"/>, identical to the historical table). DNA Chisel builds
    /// its families the same way from a Biopython codon table
    /// (<c>get_backtranslation_table(genetic_table)</c>).
    /// </param>
    public static CodonUsageTable CreateCodonTableFromSequence(
        string referenceSequence, string organismName, GeneticCode? geneticCode = null)
    {
        var families = FamiliesFor(geneticCode);
        var usage = CalculateCodonUsage(referenceSequence);
        var frequencies = new Dictionary<string, double>(64);

        foreach (var family in families.Synonyms.Values)
        {
            // Sharp & Li (1987) / Biopython: codons absent from the reference set count as 0.5.
            double total = 0;
            var counts = new double[family.Length];
            for (int k = 0; k < family.Length; k++)
            {
                int observed = usage.GetValueOrDefault(family[k], 0);
                counts[k] = observed == 0 ? 0.5 : observed;
                total += counts[k];
            }

            for (int k = 0; k < family.Length; k++)
                frequencies[family[k]] = counts[k] / total;
        }

        return new CodonUsageTable(organismName, frequencies, families.CodonToAminoAcid);
    }

    #endregion
}
