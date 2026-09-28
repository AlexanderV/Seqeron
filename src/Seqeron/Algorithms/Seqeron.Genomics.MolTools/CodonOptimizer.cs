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
    public readonly record struct CodonUsageTable(
        string OrganismName,
        IReadOnlyDictionary<string, double> CodonFrequencies,
        IReadOnlyDictionary<string, string> CodonToAminoAcid);

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

    #region Standard Genetic Code

    // Codon → amino acid and the synonymous-codon families of the NCBI Standard code (table 1),
    // RNA spelling, in NCBI codon order (UUU, UUC, UUA, … GGG) — derived from the canonical
    // GeneticCode.Standard, never from a private copy (DUP_MAP §2). The NCBI order fixes the
    // tie-break of every "most frequent synonymous codon" choice below: when two synonymous
    // codons share the maximal table frequency, the one that comes first in NCBI order wins, so
    // optimisation is deterministic (Biopython CodonAdaptationIndex.optimize only warns on such
    // ties and picks one arbitrarily).
    private static readonly IReadOnlyDictionary<string, string> StandardCodonToAminoAcid =
        GeneticCode.Standard.CodonTable.ToDictionary(kv => kv.Key, kv => kv.Value.ToString());

    private static readonly IReadOnlyDictionary<char, string[]> SynonymousCodons =
        GeneticCode.Standard.CodonTable
            .GroupBy(kv => kv.Value)
            .ToDictionary(g => g.Key, g => g.Select(kv => kv.Key).ToArray());

    /// <summary>
    /// Translates one RNA triplet with the canonical <see cref="GeneticCode.Standard"/> table.
    /// Returns <see langword="false"/> (with <paramref name="aminoAcid"/> = 'X') for a triplet
    /// that is not a valid IUPAC codon, so callers can leave such a triplet untouched instead of
    /// throwing. IUPAC-ambiguous triplets translate per Biopython (GCN → A, UAR → *, NNN → X).
    /// </summary>
    private static bool TryTranslateCodon(string codon, out char aminoAcid)
    {
        aminoAcid = 'X';
        if (codon.Length != 3)
            return false;

        foreach (char c in codon)
        {
            char dna = char.ToUpperInvariant(c);
            if (dna == 'U') dna = 'T';
            if (!IupacHelper.IsNucleotideCode(dna))
                return false;
        }

        aminoAcid = GeneticCode.Standard.Translate(codon);
        return true;
    }

    // Synonymous codons of the amino acid encoded by <paramref name="codon"/>, or an empty span
    // when the triplet is unusable (non-IUPAC symbols) or its amino acid is itself ambiguous
    // (B, Z, J, X — no unique synonymous family).
    private static string[] SynonymsFor(string codon)
    {
        if (!TryTranslateCodon(codon, out char aa))
            return Array.Empty<string>();
        return SynonymousCodons.TryGetValue(aa, out var syn) ? syn : Array.Empty<string>();
    }

    #endregion

    #region Predefined Codon Usage Tables

    /// <summary>
    /// E. coli K12 codon usage frequencies (relative fraction per amino acid).
    /// Source: Kazusa Codon Usage Database, species=316407 (E. coli K-12 substr. W3110, 4332 CDS).
    /// URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=316407
    /// </summary>
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
        StandardCodonToAminoAcid);

    /// <summary>
    /// Human codon usage frequencies (relative fraction per amino acid).
    /// Source: Kazusa Codon Usage Database, species=9606.
    /// URL: https://www.kazusa.or.jp/codon/cgi-bin/showcodon.cgi?species=9606
    /// </summary>
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
        StandardCodonToAminoAcid);

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
    /// <see cref="OptimizationResult.ProteinSequence"/>. Amino acids are taken from the canonical
    /// <see cref="GeneticCode.Standard"/> table.
    /// <para>Strategies:</para>
    /// <list type="bullet">
    /// <item><description><b>MaximizeCAI</b> — every codon becomes the most frequent synonymous
    /// codon of the target table ("one amino acid – one codon"; Puigbò et&#160;al. 2007 OPTIMIZER,
    /// DNA Chisel <c>use_best_codon</c>/<c>MaximizeCAI</c>, Biopython
    /// <c>CodonAdaptationIndex.optimize</c>). Ties are broken by NCBI codon order.</description></item>
    /// <item><description><b>AvoidRareCodeons</b> — only codons whose table frequency is strictly
    /// below <paramref name="rareCodonThreshold"/> are replaced, by the most frequent synonymous
    /// codon; codons at or above the threshold are kept.</description></item>
    /// <item><description><b>BalancedOptimization</b> (default) and <b>MinimizeSecondary</b> —
    /// as MaximizeCAI, then a GC-balancing pass moves the overall GC fraction into
    /// [<paramref name="gcTargetMin"/>, <paramref name="gcTargetMax"/>] with further synonymous
    /// swaps (DNA Chisel <c>EnforceGCContent(mini, maxi)</c> combined with codon optimisation).
    /// <c>MinimizeSecondary</c> shares this codon selection; the dedicated structure pass is
    /// <see cref="ReduceSecondaryStructure"/>.</description></item>
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
    public static OptimizationResult OptimizeSequence(
        string codingSequence,
        CodonUsageTable targetOrganism,
        OptimizationStrategy strategy = OptimizationStrategy.BalancedOptimization,
        double gcTargetMin = 0.40,
        double gcTargetMax = 0.60,
        double rareCodonThreshold = 0.15)
    {
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

        double originalCAI = CalculateCAI(rna, targetOrganism);
        var proteinBuilder = new StringBuilder();
        foreach (string codon in originalCodons)
            proteinBuilder.Append(TryTranslateCodon(codon, out char aa) ? aa : 'X');

        List<string> optimizedCodons = strategy == OptimizationStrategy.HarmonizeExpression
            ? MatchTargetCodonUsage(originalCodons, targetOrganism)
            : originalCodons
                .Select(c => SelectOptimalCodon(c, targetOrganism, strategy, rareCodonThreshold))
                .ToList();

        // GC balancing for the strategies that declare a GC target.
        if (strategy is OptimizationStrategy.BalancedOptimization or OptimizationStrategy.MinimizeSecondary)
            BalanceGcContent(optimizedCodons, targetOrganism, gcTargetMin, gcTargetMax, rareCodonThreshold);

        string optimizedSequence = string.Concat(optimizedCodons);

        var changes = new List<(int Position, string Original, string Optimized)>();
        for (int i = 0; i < originalCodons.Count; i++)
        {
            if (originalCodons[i] != optimizedCodons[i])
                changes.Add((i * 3, originalCodons[i], optimizedCodons[i]));
        }

        double optimizedCAI = CalculateCAI(optimizedSequence, targetOrganism);

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
    private static string BestSynonymousCodon(string codon, CodonUsageTable table)
    {
        var synonyms = SynonymsFor(codon);
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

    private static string SelectOptimalCodon(string currentCodon, CodonUsageTable table, OptimizationStrategy strategy, double rareCodonThreshold)
    {
        if (!TryTranslateCodon(currentCodon, out char aminoAcid) || aminoAcid == '*')
            return currentCodon; // stop codons and unusable triplets are preserved

        if (strategy == OptimizationStrategy.AvoidRareCodeons &&
            table.CodonFrequencies.GetValueOrDefault(currentCodon, 0) >= rareCodonThreshold)
        {
            return currentCodon; // not rare → untouched
        }

        return BestSynonymousCodon(currentCodon, table);
    }

    /// <summary>
    /// Rewrites the codons so that, per amino acid, the codon counts are the largest-remainder
    /// rounding of the target table's relative frequencies — the deterministic optimum of DNA
    /// Chisel's <c>MatchTargetCodonUsage</c> objective (score = −Σ_aa n_aa·Σ_codon |f_seq − f_table|),
    /// which that library reaches by randomised local search. Positions already carrying an
    /// allotted codon keep it, so only the surplus positions are edited. Stop codons and triplets
    /// without a usable synonymous family are left unchanged.
    /// </summary>
    private static List<string> MatchTargetCodonUsage(List<string> codons, CodonUsageTable table)
    {
        var result = new List<string>(codons);

        // Group the editable positions by amino acid (NCBI-ordered families).
        var positionsByAminoAcid = new Dictionary<char, List<int>>();
        for (int i = 0; i < codons.Count; i++)
        {
            if (!TryTranslateCodon(codons[i], out char aa) || aa == '*' || !SynonymousCodons.ContainsKey(aa))
                continue;
            if (!positionsByAminoAcid.TryGetValue(aa, out var positions))
                positionsByAminoAcid[aa] = positions = new List<int>();
            positions.Add(i);
        }

        foreach (var (aminoAcid, positions) in positionsByAminoAcid)
        {
            var family = SynonymousCodons[aminoAcid];
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
        double minCodonFrequency)
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
            if (!TryTranslateCodon(current, out char aa) || aa == '*')
                continue;

            int currentGc = CountGc(current);
            string? best = null;
            int bestGc = currentGc;
            double bestFrequency = double.NegativeInfinity;
            bool bestInRange = false;

            foreach (string alternative in SynonymsFor(current))
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

    // G+C count of a codon via the canonical counting primitive (SequenceExtensions).
    private static int CountGc(string codon) => codon.AsSpan().CountGcAndValidNucleotides().GcCount;

    #endregion

    #region CAI Calculation

    /// <summary>
    /// Calculates the Codon Adaptation Index (CAI) of Sharp &amp; Li (1987) against a codon-usage
    /// frequency table: <c>CAI = exp((1/L) Σ ln w_k)</c> with <c>w_ij = f_ij / max_j f_ij</c> over the
    /// synonymous codons of amino acid i. Delegates to the canonical
    /// <see cref="CodonUsageAnalyzer.CalculateCai(string, IReadOnlyDictionary{string, double}, GeneticCode)"/>
    /// core (CodonW <c>cai_out</c> conventions), under the Standard genetic code.
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
    /// <exception cref="ArgumentException"><paramref name="table"/> has no frequency dictionary.</exception>
    /// <exception cref="ArgumentOutOfRangeException">A frequency is negative or not finite.</exception>
    public static double CalculateCAI(string codingSequence, CodonUsageTable table, bool excludeSingleCodonAminoAcids = true)
    {
        if (string.IsNullOrEmpty(codingSequence))
            return 0;
        if (table.CodonFrequencies is null)
            throw new ArgumentException("Codon usage table has no frequencies.", nameof(table));

        // The canonical core keys codons in DNA spelling.
        var reference = new Dictionary<string, double>(table.CodonFrequencies.Count);
        foreach (var (codon, frequency) in table.CodonFrequencies)
            reference[codon.ToUpperInvariant().Replace('U', 'T')] = frequency;

        return CodonUsageAnalyzer.CalculateCai(
            codingSequence, reference, GeneticCode.Standard, excludeSingleCodonFamilies: excludeSingleCodonAminoAcids);
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
    public static string RemoveRestrictionSites(string codingSequence, IEnumerable<string> restrictionSites, CodonUsageTable table)
    {
        if (string.IsNullOrEmpty(codingSequence))
            return "";

        ArgumentNullException.ThrowIfNull(restrictionSites);

        string rna = ToUpperRna(codingSequence);
        var codons = SplitIntoCodons(rna);

        // Match on DNA spelling: 'U' is not an IUPAC DNA code.
        var patterns = new List<string>();
        foreach (var site in restrictionSites)
        {
            if (string.IsNullOrEmpty(site))
                continue;
            string pattern = site.ToUpperInvariant().Replace('U', 'T');
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

                if (TryRemoveOccurrence(codons, position, pattern, table))
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
    private static bool TryRemoveOccurrence(List<string> codons, int position, string pattern, CodonUsageTable table)
    {
        int firstCodon = position / 3;
        int lastCodon = Math.Min((position + pattern.Length - 1) / 3, codons.Count - 1);

        int bestIndex = -1;
        string? bestCodon = null;
        double bestFrequency = double.NegativeInfinity;

        for (int i = firstCodon; i <= lastCodon; i++)
        {
            string original = codons[i];
            if (!TryTranslateCodon(original, out char aa) || aa == '*')
                continue;

            foreach (string alternative in SynonymsFor(original))
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
    /// Reduces mRNA secondary structure by replacing codons inside self-complementary windows
    /// with synonymous codons that lower the window's self-complementarity.
    /// </summary>
    /// <remarks>
    /// The score of a window is the fraction of position pairs (i, j), j ≥ i + 4, whose bases can
    /// form a canonical pair — Watson-Crick A·U / G·C <b>or the G·U wobble</b>, taken from the
    /// canonical <see cref="RnaSecondaryStructure.CanPair"/> (ViennaRNA default pair set), not a
    /// private Watson-Crick-only copy. Windows scoring above
    /// <paramref name="structureThreshold"/> are rewritten codon by codon, each codon taking the
    /// synonymous codon that minimises the score of the <b>current</b> window content (the
    /// baseline is re-evaluated after every accepted change), among codons whose usage frequency
    /// in <paramref name="table"/> is at least <paramref name="minCodonFrequency"/>.
    /// The sequence is upper-cased, read as RNA and trimmed to complete codons; a sequence
    /// shorter than <paramref name="windowSize"/> is returned normalised but otherwise unchanged.
    /// <para>
    /// <b>Limitation.</b> The window score is a base-pair-count heuristic, not a thermodynamic
    /// folding model: it ignores stacking, loop penalties and pair nesting. A free-energy-guided
    /// search would evaluate <see cref="RnaSecondaryStructure.CalculateMinimumFreeEnergy"/>
    /// (Turner 2004) for every candidate codon, which is O(w³) per candidate (~6 ms per 43-nt
    /// window here, i.e. minutes per kilobase) and is therefore left to the caller.
    /// </para>
    /// </remarks>
    /// <param name="codingSequence">Coding sequence (DNA or RNA).</param>
    /// <param name="table">Codon-usage table of the expression host.</param>
    /// <param name="windowSize">Window width in nucleotides (default 40).</param>
    /// <param name="structureThreshold">Score above which a window is rewritten (default 0.5).</param>
    /// <param name="minCodonFrequency">Minimum usage frequency of a replacement codon (default 0.1).</param>
    public static string ReduceSecondaryStructure(
        string codingSequence,
        CodonUsageTable table,
        int windowSize = 40,
        double structureThreshold = 0.5,
        double minCodonFrequency = 0.1)
    {
        if (string.IsNullOrEmpty(codingSequence))
            return codingSequence;

        string rna = ToUpperRna(codingSequence);
        if (rna.Length < windowSize)
            return rna;

        var codons = SplitIntoCodons(rna);
        int windowCodons = windowSize / 3 + 1;

        for (int i = 0; i + windowCodons <= codons.Count; i++)
        {
            if (CalculateLocalStructure(ConcatCodons(codons, i, windowCodons, -1, "")) <= structureThreshold)
                continue;

            for (int j = i; j < i + windowCodons; j++)
            {
                string current = codons[j];
                double bestScore = CalculateLocalStructure(ConcatCodons(codons, i, windowCodons, -1, current));
                string bestAlternative = current;

                foreach (string alternative in SynonymsFor(current))
                {
                    if (alternative == current ||
                        table.CodonFrequencies.GetValueOrDefault(alternative, 0) < minCodonFrequency)
                    {
                        continue;
                    }

                    double score = CalculateLocalStructure(ConcatCodons(codons, i, windowCodons, j, alternative));

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
    /// Stop codons are screened like any other row of the table. The amino acid is the NCBI
    /// Standard-code translation (<see cref="GeneticCode.Standard"/>).
    /// </remarks>
    /// <returns>
    /// (0-based nucleotide position of the codon, RNA codon, one-letter amino acid or <c>*</c>,
    /// table frequency) in sequence order.
    /// </returns>
    public static IEnumerable<(int Position, string Codon, string AminoAcid, double Frequency)> FindRareCodons(
        string codingSequence,
        CodonUsageTable table,
        double threshold = 0.15)
    {
        var codons = SplitIntoScreenedCodons(codingSequence);

        for (int i = 0; i < codons.Length; i++)
        {
            string? codon = codons[i];
            if (codon is null)
                continue;

            double freq = table.CodonFrequencies.GetValueOrDefault(codon, 0);
            if (freq < threshold)
                yield return (i * 3, codon, GeneticCode.Standard.Translate(codon).ToString(), freq);
        }
    }

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
    /// {UAA, UAG, UGA}. Families here come from <see cref="GeneticCode.Standard"/>; codons of
    /// single-codon amino acids (Met, Trp) contribute 0 to both numerator and denominator, and an
    /// ambiguous triplet (non-ACGU symbol) occupies its window position but contributes nothing.
    /// <para>
    /// <b>Frequency scale.</b> The reference implementation takes a codon usage table in
    /// frequency per thousand codons (Kazusa "/1000" column), so the window sums weight each
    /// residue by its overall usage. Passing a <see cref="CodonUsageTable"/> whose
    /// <c>CodonFrequencies</c> hold per-thousand values reproduces the reference exactly. The
    /// built-in presets (<see cref="EColiK12"/>, <see cref="Yeast"/>, <see cref="Human"/>) hold
    /// per-amino-acid relative fractions; with them every residue has equal weight, which gives
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
    public static IReadOnlyList<MinMaxWindow> CalculateMinMaxProfile(
        string codingSequence,
        CodonUsageTable table,
        int windowSize = DefaultMinMaxWindowCodons)
    {
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1 codon.");

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

            xij[i] = table.CodonFrequencies.GetValueOrDefault(codon, 0);
            char aa = GeneticCode.Standard.Translate(codon);
            if (!familyStats.TryGetValue(aa, out var stats))
            {
                double sum = 0, max = double.MinValue, min = double.MaxValue;
                int n = 0;
                foreach (string syn in GeneticCode.Standard.GetCodonsForAminoAcid(aa))
                {
                    double f = table.CodonFrequencies.GetValueOrDefault(syn, 0);
                    sum += f;
                    n++;
                    if (f > max) max = f;
                    if (f < min) min = f;
                }
                stats = (sum / n, max, min);
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
    /// A codon is "rare"/"pause" when its usage frequency in <paramref name="table"/> is strictly
    /// below <paramref name="rareThreshold"/> — the same per-codon criterion (and the same
    /// unambiguous-codon screen: an ambiguous triplet is never a pause) as
    /// <see cref="FindRareCodons"/>. Overlapping windows are merged into maximal clusters so a long
    /// rare run is reported once. This is opt-in; <see cref="FindRareCodons"/> (per-codon) is
    /// unchanged. Defaults reproduce the published Sherlocc rule "a seven position-wide window …
    /// containing at least four pause positions out of seven" (Chartier et&#160;al. 2012,
    /// doi:10.1093/bioinformatics/bts149).
    /// </remarks>
    /// <param name="codingSequence">DNA or RNA coding sequence (T is normalised to U).</param>
    /// <param name="table">Reference codon-usage table (per-amino-acid relative fractions).</param>
    /// <param name="rareThreshold">Per-codon rare-frequency cutoff (default 0.15, strict &lt;).</param>
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
    {
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1 codon.");
        if (minRareCodons < 1)
            throw new ArgumentOutOfRangeException(nameof(minRareCodons), minRareCodons, "Minimum rare codons must be at least 1.");

        var clusters = new List<RareCodonCluster>();
        var codons = SplitIntoScreenedCodons(codingSequence);
        if (codons.Length < windowSize)
            return clusters;

        // Mark each codon as rare (pause) when its table frequency is strictly below the threshold.
        // An ambiguous triplet has no usage frequency and is never a pause position.
        var isRare = new bool[codons.Length];
        for (int i = 0; i < codons.Length; i++)
            isRare[i] = codons[i] is { } codon
                && table.CodonFrequencies.GetValueOrDefault(codon, 0) < rareThreshold;

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
    private static string ToUpperRna(string sequence) => sequence.ToUpperInvariant().Replace('T', 'U');

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
    /// <see cref="GeneticCode.Standard"/> mapping (RNA spelling) that the optimizer uses, so a
    /// caller-supplied table carries the same codon → amino-acid assignment as the built-in
    /// presets. Frequency keys are normalised to upper-case RNA spelling.
    /// </summary>
    /// <param name="organismName">Name stored in the table.</param>
    /// <param name="codonFrequencies">Per-amino-acid relative codon frequencies (DNA or RNA keys).</param>
    public static CodonUsageTable CreateCodonUsageTable(
        string organismName,
        IReadOnlyDictionary<string, double> codonFrequencies)
    {
        ArgumentNullException.ThrowIfNull(codonFrequencies);

        var frequencies = new Dictionary<string, double>(codonFrequencies.Count);
        foreach (var (codon, frequency) in codonFrequencies)
            frequencies[ToUpperRna(codon)] = frequency;

        return new CodonUsageTable(organismName, frequencies, StandardCodonToAminoAcid);
    }

    /// <summary>
    /// Builds a codon-usage table from a reference gene set (one concatenated in-frame coding
    /// sequence, or a single gene), as the relative frequency of each codon within its
    /// synonymous family of the Standard genetic code.
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
    /// to 1.
    /// </remarks>
    /// <param name="referenceSequence">Reference coding sequence(s), in frame 0 (DNA or RNA).</param>
    /// <param name="organismName">Name stored in the resulting table.</param>
    public static CodonUsageTable CreateCodonTableFromSequence(string referenceSequence, string organismName)
    {
        var usage = CalculateCodonUsage(referenceSequence);
        var frequencies = new Dictionary<string, double>(64);

        foreach (var family in SynonymousCodons.Values)
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

        return new CodonUsageTable(organismName, frequencies, StandardCodonToAminoAcid);
    }

    #endregion
}
