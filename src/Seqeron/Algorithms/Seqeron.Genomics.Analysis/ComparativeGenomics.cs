using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Provides comparative genomics algorithms for analyzing relationships between genomes.
/// Includes synteny detection, genome rearrangements, and ortholog identification.
/// </summary>
public static class ComparativeGenomics
{
    /// <summary>
    /// Represents a syntenic block between two genomes.
    /// </summary>
    public readonly record struct SyntenicBlock(
        string Genome1Id,
        int Start1,
        int End1,
        string Genome2Id,
        int Start2,
        int End2,
        bool IsInverted,
        int GeneCount,
        double Identity);

    /// <summary>
    /// Represents an orthologous gene pair.
    /// </summary>
    public readonly record struct OrthologPair(
        string Gene1Id,
        string Gene2Id,
        double Identity,
        double Coverage,
        int AlignmentLength);

    /// <summary>
    /// Represents a genome rearrangement event.
    /// </summary>
    public readonly record struct RearrangementEvent(
        RearrangementType Type,
        string GenomeId,
        int Position,
        int Length,
        string? TargetPosition = null);

    /// <summary>
    /// Types of genome rearrangements.
    /// </summary>
    public enum RearrangementType
    {
        Inversion,
        Translocation,
        Deletion,
        Insertion,
        Duplication,
        Transposition
    }

    /// <summary>
    /// Represents a gene for comparative analysis.
    /// </summary>
    public readonly record struct Gene(
        string Id,
        string GenomeId,
        int Start,
        int End,
        char Strand,
        string? Sequence = null);

    /// <summary>
    /// Result of comparative genome analysis.
    /// </summary>
    public readonly record struct ComparisonResult(
        IReadOnlyList<SyntenicBlock> SyntenicBlocks,
        IReadOnlyList<OrthologPair> Orthologs,
        IReadOnlyList<RearrangementEvent> Rearrangements,
        double OverallSynteny,
        int ConservedGenes,
        int GenomeSpecificGenes1,
        int GenomeSpecificGenes2);

    // --- MCScanX collinearity scoring constants (Wang et al. 2012, NAR 40(7):e49) ---
    // Score(v) = max(MatchScore(v), max_u(Score(u) + MatchScore(v) + GapPenalty * NumberofGaps(u,v)))

    /// <summary>Score awarded per anchored (collinear) gene pair. MCScanX default MatchScore = 50.</summary>
    private const int MatchScore = 50;

    /// <summary>Penalty applied per intervening gene between two anchors. MCScanX default GapPenalty = −1.</summary>
    private const int GapPenalty = -1;

    /// <summary>
    /// Maximum number of intervening genes permitted between consecutive anchors.
    /// MCScanX requires NumberofGaps(u,v) &lt; 25 (default MAX_GAPS = 25).
    /// </summary>
    private const int DefaultMaxGaps = 25;

    /// <summary>
    /// Minimum chain score for a collinear block to be reported.
    /// MCScanX reports non-overlapping chains with scores over 250.
    /// </summary>
    private const int MinChainScore = 250;

    /// <summary>
    /// Minimum number of collinear gene pairs (anchors) per reported block.
    /// MCScanX default = 5 (5 × MatchScore = 250).
    /// </summary>
    private const int DefaultMinAnchors = MinChainScore / MatchScore;

    /// <summary>
    /// Detects syntenic (collinear) blocks between two gene-ordered genomes from a set of
    /// orthologous anchors, following the MCScanX collinearity model: a dynamic-programming
    /// chaining of anchor pairs that rewards adjacency and penalizes intervening genes
    /// (Wang et al. 2012, <i>Nucleic Acids Research</i> 40(7):e49).
    /// </summary>
    /// <param name="genome1Genes">Genes of genome 1 in chromosomal order.</param>
    /// <param name="genome2Genes">Genes of genome 2 in chromosomal order.</param>
    /// <param name="orthologMap">Anchor map: genome-1 gene id → orthologous genome-2 gene id.</param>
    /// <param name="minAnchors">
    /// Minimum collinear gene pairs per block. MCScanX default 5 (score threshold 250).
    /// </param>
    /// <param name="maxGap">
    /// Maximum intervening genes between consecutive anchors. MCScanX default 25 (NumberofGaps &lt; 25).
    /// </param>
    /// <returns>The reported non-overlapping collinear blocks (forward and inverted).</returns>
    /// <exception cref="ArgumentNullException">Any required argument is null.</exception>
    public static IEnumerable<SyntenicBlock> FindSyntenicBlocks(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        IReadOnlyDictionary<string, string> orthologMap,
        int minAnchors = DefaultMinAnchors,
        int maxGap = DefaultMaxGaps)
    {
        ArgumentNullException.ThrowIfNull(genome1Genes);
        ArgumentNullException.ThrowIfNull(genome2Genes);
        ArgumentNullException.ThrowIfNull(orthologMap);

        if (genome1Genes.Count == 0 || genome2Genes.Count == 0)
            return Array.Empty<SyntenicBlock>();

        // Build genome-2 position index.
        var genome2Positions = new Dictionary<string, int>();
        for (int i = 0; i < genome2Genes.Count; i++)
            genome2Positions[genome2Genes[i].Id] = i;

        // Collect anchors as (pos1, pos2) ordered by genome-1 position.
        var anchors = new List<(int pos1, int pos2)>();
        for (int i = 0; i < genome1Genes.Count; i++)
        {
            if (orthologMap.TryGetValue(genome1Genes[i].Id, out string? ortholog) &&
                genome2Positions.TryGetValue(ortholog, out int pos2))
            {
                anchors.Add((i, pos2));
            }
        }

        if (anchors.Count < minAnchors)
            return Array.Empty<SyntenicBlock>();

        var chains = ChainCollinearAnchors(anchors, minAnchors, maxGap);
        return BuildBlocks(chains, genome1Genes, genome2Genes);
    }

    /// <summary>
    /// Greedily extracts non-overlapping collinear chains using the MCScanX scoring scheme.
    /// Anchors are sorted by genome-1 position; consecutive anchors extend a chain when they
    /// keep a consistent genome-2 direction and the number of intervening genes is &lt; maxGap.
    /// A chain is reported when its MCScanX score ≥ <see cref="MinChainScore"/> and it has ≥ minAnchors pairs.
    /// </summary>
    private static List<List<(int pos1, int pos2)>> ChainCollinearAnchors(
        List<(int pos1, int pos2)> anchors,
        int minAnchors,
        int maxGap)
    {
        var reported = new List<List<(int pos1, int pos2)>>();

        // Sort by genome-1 position (anchors are already appended in this order, but be explicit).
        anchors.Sort((a, b) => a.pos1.CompareTo(b.pos1));

        var current = new List<(int pos1, int pos2)>();
        int? direction = null; // +1 forward, -1 reverse.

        void Flush()
        {
            if (IsReportable(current, minAnchors))
                reported.Add(new List<(int pos1, int pos2)>(current));
            current.Clear();
            direction = null;
        }

        foreach (var anchor in anchors)
        {
            if (current.Count == 0)
            {
                current.Add(anchor);
                continue;
            }

            var last = current[^1];
            int delta2 = anchor.pos2 - last.pos2;
            // NumberofGaps = intervening genes between the two anchors in genome 2.
            int numberOfGaps = Math.Abs(delta2) - 1;
            int currentDir = delta2 > 0 ? 1 : -1;

            bool extends = delta2 != 0
                           && numberOfGaps < maxGap
                           && (direction == null || direction == currentDir);

            if (extends)
            {
                current.Add(anchor);
                direction = currentDir;
            }
            else
            {
                Flush();
                current.Add(anchor);
            }
        }

        Flush();
        return reported;
    }

    /// <summary>
    /// Computes the MCScanX chain score for an ordered run of anchors and tests the report rule:
    /// score ≥ <see cref="MinChainScore"/> AND anchor count ≥ minAnchors.
    /// </summary>
    private static bool IsReportable(List<(int pos1, int pos2)> chain, int minAnchors)
    {
        if (chain.Count < minAnchors)
            return false;

        // Score = n * MatchScore + GapPenalty * (total intervening genes).
        int totalGaps = 0;
        for (int i = 1; i < chain.Count; i++)
            totalGaps += Math.Abs(chain[i].pos2 - chain[i - 1].pos2) - 1;

        int score = chain.Count * MatchScore + GapPenalty * totalGaps;
        return score >= MinChainScore;
    }

    private static List<SyntenicBlock> BuildBlocks(
        List<List<(int pos1, int pos2)>> chains,
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes)
    {
        var blocks = new List<SyntenicBlock>(chains.Count);

        foreach (var chain in chains)
        {
            int start1 = genome1Genes[chain.Min(a => a.pos1)].Start;
            int end1 = genome1Genes[chain.Max(a => a.pos1)].End;
            int start2 = genome2Genes[chain.Min(a => a.pos2)].Start;
            int end2 = genome2Genes[chain.Max(a => a.pos2)].End;

            // Inverted iff genome-2 order decreases along the (genome-1-ordered) chain.
            bool isInverted = chain[^1].pos2 < chain[0].pos2;

            blocks.Add(new SyntenicBlock(
                Genome1Id: genome1Genes[0].GenomeId,
                Start1: Math.Min(start1, end1),
                End1: Math.Max(start1, end1),
                Genome2Id: genome2Genes[0].GenomeId,
                Start2: Math.Min(start2, end2),
                End2: Math.Max(start2, end2),
                IsInverted: isInverted,
                GeneCount: chain.Count,
                Identity: 1.0));
        }

        return blocks;
    }

    /// <summary>
    /// Renders syntenic blocks as a human-readable text summary, one line per block.
    /// Visualization helper for <see cref="FindSyntenicBlocks"/>.
    /// </summary>
    public static string VisualizeSynteny(IReadOnlyList<SyntenicBlock> blocks)
    {
        ArgumentNullException.ThrowIfNull(blocks);

        var sb = new StringBuilder();
        for (int i = 0; i < blocks.Count; i++)
        {
            var b = blocks[i];
            string orientation = b.IsInverted ? "reverse" : "forward";
            sb.Append(b.Genome1Id).Append(':').Append(b.Start1).Append('-').Append(b.End1)
              .Append(" <=> ")
              .Append(b.Genome2Id).Append(':').Append(b.Start2).Append('-').Append(b.End2)
              .Append(" [").Append(orientation).Append(", ")
              .Append(b.GeneCount).Append(" genes]");
            if (i < blocks.Count - 1)
                sb.Append('\n');
        }

        return sb.ToString();
    }

    // --- Reciprocal Best Hit (RBH) ortholog-detection constants ---
    // Moreno-Hagelsieb & Latimer (2008), Bioinformatics 24(3):319–324,
    // https://doi.org/10.1093/bioinformatics/btm585 :
    // "two genes residing in two different genomes are deemed orthologs if their protein
    //  products find each other as the best hit in the opposite genome."

    /// <summary>
    /// Minimum fraction of the shorter sequence that must be covered by the match.
    /// Moreno-Hagelsieb &amp; Latimer (2008) require "coverage of at least 50% of any of
    /// the protein sequences in the alignments" (default 0.5).
    /// </summary>
    private const double DefaultMinCoverage = 0.5;

    /// <summary>
    /// Minimum similarity score for a hit to qualify, mapping the significance gate
    /// (max E-value 1e-6 in Moreno-Hagelsieb &amp; Latimer 2008) onto the alignment-free
    /// similarity used here (default 0.3). See Evidence Assumption 1.
    /// </summary>
    private const double DefaultMinIdentity = 0.3;

    /// <summary>
    /// Identifies orthologous gene pairs as <b>reciprocal best hits</b> (RBH/BBH): two genes,
    /// one in each genome, are orthologs iff each gene's best qualifying hit in the other genome
    /// is the other gene. This is the symmetrical-best-hit criterion of Tatusov, Koonin &amp; Lipman
    /// (1997, <i>Science</i> 278:631–637) and the operational RBH definition of Moreno-Hagelsieb &amp;
    /// Latimer (2008, <i>Bioinformatics</i> 24:319–324). A one-directional best hit is NOT an ortholog.
    /// </summary>
    /// <param name="genome1Genes">Genes of genome 1 (each must carry a non-empty <see cref="Gene.Sequence"/>).</param>
    /// <param name="genome2Genes">Genes of genome 2.</param>
    /// <param name="minIdentity">Minimum similarity score for a qualifying hit (default 0.3).</param>
    /// <param name="minCoverage">Minimum coverage fraction for a qualifying hit (default 0.5).</param>
    /// <returns>The reciprocal best-hit ortholog pairs (an unordered matching: each gene appears at most once).</returns>
    /// <exception cref="ArgumentNullException">Either gene list is null.</exception>
    public static IEnumerable<OrthologPair> FindOrthologs(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        double minIdentity = DefaultMinIdentity,
        double minCoverage = DefaultMinCoverage)
        // Orthology by RBH is exactly the reciprocal-best-hit criterion; delegate to the
        // canonical RBH implementation so the two entry points cannot diverge.
        => FindReciprocalBestHits(genome1Genes, genome2Genes, minIdentity, minCoverage);

    /// <summary>
    /// Identifies within-genome paralog pairs (recent in-paralogs) as <b>mutual best hits inside the
    /// same genome</b>. Paralogs arise by gene duplication within an organism's history (Fitch 1970,
    /// <i>Syst. Zool.</i> 19:99–106); recent in-paralogs are within-species pairs reciprocally more
    /// similar to each other than to outside genes (Remm, Storm &amp; Sonnhammer 2001, <i>J. Mol. Biol.</i>
    /// 314:1041–1052). Each returned pair is an unordered pair of distinct genes that are mutual best
    /// hits within <paramref name="genes"/>.
    /// </summary>
    /// <param name="genes">Genes of a single genome (each must carry a non-empty <see cref="Gene.Sequence"/>).</param>
    /// <param name="minIdentity">Minimum similarity score for a qualifying hit (default 0.3).</param>
    /// <param name="minCoverage">Minimum coverage fraction for a qualifying hit (default 0.5).</param>
    /// <returns>Within-genome reciprocal best-hit paralog pairs (each unordered pair reported once).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="genes"/> is null.</exception>
    public static IEnumerable<OrthologPair> FindParalogs(
        IReadOnlyList<Gene> genes,
        double minIdentity = DefaultMinIdentity,
        double minCoverage = DefaultMinCoverage)
    {
        ArgumentNullException.ThrowIfNull(genes);

        var seqGenes = genes.Where(g => !string.IsNullOrEmpty(g.Sequence)).ToList();
        if (seqGenes.Count < 2)
            return Array.Empty<OrthologPair>();

        // Best within-genome hit for each gene (excluding itself).
        var bestHit = new Dictionary<string, (Gene hit, double identity, double coverage, int alignLen)>();
        foreach (var g in seqGenes)
        {
            var others = seqGenes.Where(o => o.Id != g.Id).ToList();
            var bh = FindBestHit(g, others, minIdentity, minCoverage);
            if (bh != null)
                bestHit[g.Id] = bh.Value;
        }

        // Mutual best hits, each unordered pair reported once.
        var paralogs = new List<OrthologPair>();
        var seen = new HashSet<string>();
        foreach (var g in seqGenes)
        {
            if (!bestHit.TryGetValue(g.Id, out var hit))
                continue;
            if (!bestHit.TryGetValue(hit.hit.Id, out var back) || back.hit.Id != g.Id)
                continue;

            string key = string.CompareOrdinal(g.Id, hit.hit.Id) < 0
                ? g.Id + "\0" + hit.hit.Id
                : hit.hit.Id + "\0" + g.Id;
            if (!seen.Add(key))
                continue;

            paralogs.Add(new OrthologPair(
                Gene1Id: g.Id,
                Gene2Id: hit.hit.Id,
                Identity: hit.identity,
                Coverage: hit.coverage,
                AlignmentLength: hit.alignLen));
        }

        return paralogs;
    }

    /// <summary>
    /// Returns the single best qualifying hit of <paramref name="query"/> among <paramref name="targets"/>.
    /// "Best" = maximum similarity score; ties are broken by larger coverage, then by ordinal gene id, so
    /// the best hit is unique and deterministic (Moreno-Hagelsieb &amp; Latimer 2008: sort hits by score,
    /// break ties deterministically). A hit qualifies only when identity ≥ minIdentity and coverage ≥ minCoverage.
    /// </summary>
    private static (Gene hit, double identity, double coverage, int alignLen)? FindBestHit(
        Gene query,
        IReadOnlyList<Gene> targets,
        double minIdentity,
        double minCoverage)
    {
        (Gene hit, double identity, double coverage, int alignLen)? best = null;

        foreach (var target in targets)
        {
            if (string.IsNullOrEmpty(target.Sequence))
                continue;

            var (identity, coverage, alignLen) = CalculateSequenceSimilarity(
                query.Sequence!, target.Sequence!);

            if (identity < minIdentity || coverage < minCoverage)
                continue;

            if (best == null ||
                identity > best.Value.identity ||
                (identity == best.Value.identity && coverage > best.Value.coverage) ||
                (identity == best.Value.identity && coverage == best.Value.coverage &&
                 string.CompareOrdinal(target.Id, best.Value.hit.Id) < 0))
            {
                best = (target, identity, coverage, alignLen);
            }
        }

        return best;
    }

    /// <summary>
    /// Finds <b>reciprocal best hits</b> (RBH, a.k.a. bidirectional/symmetrical best hits) between two
    /// genomes for ortholog identification. Per Moreno-Hagelsieb &amp; Latimer (2008,
    /// <i>Bioinformatics</i> 24:319–324): "two genes residing in two different genomes are deemed orthologs
    /// if their protein products find each other as the best hit in the opposite genome." This is the
    /// symmetrical-best-hit criterion underlying the COGs of Tatusov, Koonin &amp; Lipman (1997,
    /// <i>Science</i> 278:631–637). A one-directional best hit is NOT an ortholog.
    /// </summary>
    /// <remarks>
    /// A candidate hit qualifies only when it passes the significance gate (identity ≥ <paramref name="minIdentity"/>
    /// and coverage ≥ <paramref name="minCoverage"/>, mapping the 1×10⁻⁶ E-value and ≥50% coverage gates of
    /// Moreno-Hagelsieb &amp; Latimer 2008). The best hit is the qualifying candidate with the maximum similarity
    /// score; ties are broken deterministically (larger coverage, then ordinal gene id) so the returned matching
    /// is unique and order-independent. Each returned pair carries the actual identity, coverage, and alignment
    /// length of the hit. This is the dedicated RBH entry point; <see cref="FindOrthologs"/> applies the same
    /// criterion under the "orthologs" name.
    /// </remarks>
    /// <param name="genome1Genes">Genes of genome 1 (each must carry a non-empty <see cref="Gene.Sequence"/>).</param>
    /// <param name="genome2Genes">Genes of genome 2.</param>
    /// <param name="minIdentity">Minimum similarity score for a qualifying hit (default 0.3).</param>
    /// <param name="minCoverage">Minimum coverage fraction for a qualifying hit (default 0.5, per the ≥50% gate).</param>
    /// <returns>The reciprocal best-hit pairs (an unordered matching: each gene appears at most once).</returns>
    /// <exception cref="ArgumentNullException">Either gene list is null.</exception>
    public static IEnumerable<OrthologPair> FindReciprocalBestHits(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        double minIdentity = DefaultMinIdentity,
        double minCoverage = DefaultMinCoverage)
    {
        ArgumentNullException.ThrowIfNull(genome1Genes);
        ArgumentNullException.ThrowIfNull(genome2Genes);

        var seqGenes1 = genome1Genes.Where(g => !string.IsNullOrEmpty(g.Sequence)).ToList();
        var seqGenes2 = genome2Genes.Where(g => !string.IsNullOrEmpty(g.Sequence)).ToList();

        if (seqGenes1.Count == 0 || seqGenes2.Count == 0)
            return Array.Empty<OrthologPair>();

        // Best qualifying hit of each genome-1 gene into genome 2 (with full hit metrics),
        // and the best-hit id of each genome-2 gene back into genome 1.
        var best1To2 = new Dictionary<string, (Gene hit, double identity, double coverage, int alignLen)>();
        foreach (var g1 in seqGenes1)
        {
            var bh = FindBestHit(g1, seqGenes2, minIdentity, minCoverage);
            if (bh != null)
                best1To2[g1.Id] = bh.Value;
        }

        var best2To1Id = new Dictionary<string, string>();
        foreach (var g2 in seqGenes2)
        {
            var bh = FindBestHit(g2, seqGenes1, minIdentity, minCoverage);
            if (bh != null)
                best2To1Id[g2.Id] = bh.Value.hit.Id;
        }

        // Keep only reciprocal pairs: g1's best hit is g2 AND g2's best hit is g1.
        var rbh = new List<OrthologPair>();
        foreach (var g1 in seqGenes1)
        {
            if (best1To2.TryGetValue(g1.Id, out var hit) &&
                best2To1Id.TryGetValue(hit.hit.Id, out string? backId) &&
                backId == g1.Id)
            {
                rbh.Add(new OrthologPair(
                    Gene1Id: g1.Id,
                    Gene2Id: hit.hit.Id,
                    Identity: hit.identity,
                    Coverage: hit.coverage,
                    AlignmentLength: hit.alignLen));
            }
        }

        return rbh;
    }

    private static (double identity, double coverage, int alignLength) CalculateSequenceSimilarity(
        string seq1, string seq2)
    {
        // Simple k-mer based similarity (faster than full alignment)
        const int k = 5;

        if (seq1.Length < k || seq2.Length < k)
            return (0, 0, 0);

        var kmers1 = new HashSet<string>();
        for (int i = 0; i <= seq1.Length - k; i++)
            kmers1.Add(seq1.Substring(i, k).ToUpperInvariant());

        var kmers2 = new HashSet<string>();
        for (int i = 0; i <= seq2.Length - k; i++)
            kmers2.Add(seq2.Substring(i, k).ToUpperInvariant());

        int shared = kmers1.Intersect(kmers2).Count();
        int total = kmers1.Union(kmers2).Count();

        // identity = k-mer Jaccard similarity (alignment-free; Ondov et al. 2016). Used only as the
        // best-hit ranking score (Evidence Assumption 1).
        double identity = total > 0 ? (double)shared / total : 0;
        // coverage = fraction of the SHORTER sequence's k-mers that are shared; maps the
        // ">= 50% coverage of any of the protein sequences" gate of Moreno-Hagelsieb & Latimer (2008)
        // onto k-mer space. Identical sequences => coverage 1.0.
        int minKmerCount = Math.Min(kmers1.Count, kmers2.Count);
        double coverage = minKmerCount > 0 ? (double)shared / minKmerCount : 0;
        int alignLen = Math.Min(seq1.Length, seq2.Length);

        return (identity, coverage, alignLen);
    }

    // --- Signed-permutation breakpoint model (Bafna & Pevzner 1998; Tannier et al. 2009) ---
    // A genome's marker order is a signed permutation; the target (genome 2) is relabelled to the
    // identity. The permutation is extended with a left sentinel 0 and a right sentinel n+1
    // (Bafna & Pevzner 1998, SIAM J. Discrete Math. 11(2):224–240). A consecutive pair (x, y) of
    // the extended permutation is a BREAKPOINT iff it is not an identity adjacency, i.e. iff
    // y != x + 1 (this single signed test subsumes both the (x,y) and (-y,-x) clauses of the
    // breakpoint definition because a reversal negates the signs of the block it reverses —
    // Hunter College CompBio Lecture 16). Each breakpoint marks a rearrangement boundary;
    // breakpoint count is the breakpoint distance d_BP = n - (common adjacencies),
    // a lower bound d >= b/2 on the reversal distance (Tannier et al. 2009; Hunter Lecture 16).

    /// <summary>Left sentinel prepended to the extended signed permutation (π₀ = 0). Bafna &amp; Pevzner (1998).</summary>
    private const int LeftSentinel = 0;

    /// <summary>
    /// Detects genome rearrangements between two gene orders as <b>breakpoints</b> of the signed
    /// gene-order permutation. Orthologous markers (via <paramref name="orthologMap"/>) are read in
    /// genome-1 order and relabelled to genome 2's rank with a sign for relative strand; the
    /// permutation is extended with sentinels <c>0</c> and <c>n+1</c>. Every consecutive pair that is
    /// not an identity adjacency (<c>y ≠ x + 1</c>) is a breakpoint and is reported as one
    /// <see cref="RearrangementEvent"/>, classified by <see cref="ClassifyRearrangement"/>. This is the
    /// formally defined breakpoint model of Bafna &amp; Pevzner (1998, <i>SIAM J. Discrete Math.</i>
    /// 11(2):224–240) and Tannier, Zheng &amp; Sankoff (2009); the breakpoint count equals the
    /// breakpoint distance <c>d_BP = n − (common adjacencies)</c>.
    /// </summary>
    /// <param name="genome1Genes">Genes of genome 1 in chromosomal order.</param>
    /// <param name="genome2Genes">Genes of genome 2 in chromosomal order.</param>
    /// <param name="orthologMap">Anchor map: genome-1 gene id → orthologous genome-2 gene id.</param>
    /// <returns>One <see cref="RearrangementEvent"/> per breakpoint, in genome-1 order.</returns>
    /// <exception cref="ArgumentNullException">Any required argument is null.</exception>
    public static IEnumerable<RearrangementEvent> DetectRearrangements(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        IReadOnlyDictionary<string, string> orthologMap)
    {
        ArgumentNullException.ThrowIfNull(genome1Genes);
        ArgumentNullException.ThrowIfNull(genome2Genes);
        ArgumentNullException.ThrowIfNull(orthologMap);

        return DetectRearrangementsIterator(genome1Genes, genome2Genes, orthologMap);
    }

    private static IEnumerable<RearrangementEvent> DetectRearrangementsIterator(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        IReadOnlyDictionary<string, string> orthologMap)
    {
        // Genome-2 rank index (1-based so sentinel 0 cannot collide with a real marker).
        var genome2Rank = new Dictionary<string, int>(genome2Genes.Count);
        for (int i = 0; i < genome2Genes.Count; i++)
            genome2Rank[genome2Genes[i].Id] = i + 1;

        // Build the signed permutation: for each genome-1 ortholog (in order), its genome-2 rank,
        // signed by relative strand (+ if strands agree, − if they differ). Carry the genome-1
        // gene index so we can report a coordinate per breakpoint.
        var perm = new List<(int signedRank, int g1Index)>();
        for (int i = 0; i < genome1Genes.Count; i++)
        {
            var g1 = genome1Genes[i];
            if (orthologMap.TryGetValue(g1.Id, out string? ortholog) &&
                genome2Rank.TryGetValue(ortholog, out int rank))
            {
                int sign = g1.Strand == genome2Genes[rank - 1].Strand ? 1 : -1;
                perm.Add((sign * rank, i));
            }
        }

        // A permutation of fewer than 2 markers has no internal adjacency, hence no breakpoint.
        if (perm.Count < 2)
            yield break;

        int n = perm.Count;
        int rightSentinel = n + 1; // π_{n+1} = n + 1 (Bafna & Pevzner 1998).

        // Walk consecutive pairs of the extended permutation (0, perm..., n+1).
        // The two sentinel pairs at the ends are evaluated against the real markers' rank order,
        // not their genome-2 ranks, to remain a permutation of 0..n+1. We therefore compare the
        // *order ranks* of consecutive markers; a pair is a breakpoint iff its order-rank successor
        // is not exactly +1 with matching sign — equivalently iff |Δrank| != 1 or the signs flip.
        // To stay faithful to the signed-permutation definition we relabel the genome-2 ranks to
        // their order positions 1..n (the relative permutation), preserving sign.
        var orderOfRank = new Dictionary<int, int>(n);
        var sortedByAbs = perm.Select((p, idx) => (abs: Math.Abs(p.signedRank), idx))
                              .OrderBy(t => t.abs)
                              .ToList();
        for (int r = 0; r < sortedByAbs.Count; r++)
            orderOfRank[sortedByAbs[r].idx] = r + 1; // 1..n in genome-2 order

        // relabelled[i] = signed order position of marker i (sign preserved from genome-2 strand).
        var relabelled = new int[n];
        for (int i = 0; i < n; i++)
            relabelled[i] = Math.Sign(perm[i].signedRank) * orderOfRank[i];

        // Extended sequence of signed values: [0, relabelled..., n+1].
        // Evaluate each consecutive pair (x, y): breakpoint iff y != x + 1.
        int prev = LeftSentinel;
        for (int i = 0; i <= n; i++)
        {
            int curr = i < n ? relabelled[i] : rightSentinel;

            if (curr != prev + 1)
            {
                // This boundary lies just before marker i (or at the right sentinel).
                int g1Index = i < n ? perm[i].g1Index : perm[n - 1].g1Index;
                var anchorGene = genome1Genes[g1Index];

                yield return new RearrangementEvent(
                    Type: ClassifyBoundary(prev, curr),
                    GenomeId: genome1Genes[0].GenomeId,
                    Position: anchorGene.Start,
                    // Length encodes the signed pair so ClassifyRearrangement can re-derive the type:
                    // |Δ| of the signed values across the breakpoint.
                    Length: Math.Abs(curr - prev),
                    TargetPosition: $"{prev}->{curr}");
            }

            prev = curr;
        }
    }

    /// <summary>
    /// Classifies a single breakpoint boundary (the signed pair <paramref name="x"/>→<paramref name="y"/>
    /// of the extended permutation) into a <see cref="RearrangementType"/>: a sign reversal across the
    /// boundary (one value negative w.r.t. an otherwise-consecutive relation) indicates an
    /// <b>inversion</b> (a reversal negates the signs of the block it reverses — Hunter College CompBio
    /// Lecture 16); an orientation-preserving discontinuity (both values positive, non-consecutive)
    /// indicates a <b>transposition</b> (a block relocated to a new position preserving orientation —
    /// Bafna &amp; Pevzner 1998).
    /// </summary>
    private static RearrangementType ClassifyBoundary(int x, int y)
    {
        // A reversal produces an internal adjacency of the form (-(k+1), -k); at a reversal *boundary*
        // exactly one of the two flanking values carries a flipped sign relative to a forward chain.
        // Detect the sign flip: if the boundary separates values of opposite sign (and is not the
        // sentinel start where x=0), it was created by a reversal => Inversion. Otherwise the block
        // moved without re-orientation => Transposition.
        bool signFlip = (Math.Sign(x) != Math.Sign(y)) && x != LeftSentinel;
        bool negativeInvolved = x < 0 || y < 0;

        if (signFlip || negativeInvolved)
            return RearrangementType.Inversion;

        return RearrangementType.Transposition;
    }

    /// <summary>
    /// Classifies a detected <see cref="RearrangementEvent"/> by its breakpoint signature. The event's
    /// <see cref="RearrangementEvent.TargetPosition"/> records the signed pair <c>x-&gt;y</c> spanning the
    /// breakpoint; a sign reversal across the pair denotes an <b>inversion</b> (a reversal negates the
    /// signs of the reversed block — Hunter College CompBio Lecture 16), while an orientation-preserving
    /// discontinuity denotes a <b>transposition</b> (Bafna &amp; Pevzner 1998, <i>SIAM J. Discrete Math.</i>
    /// 11(2):224–240). Only these two operation classes are derivable from a single signed in-order
    /// permutation; translocation / deletion / insertion / duplication require chromosome identifiers or
    /// gene-set differences and are out of scope for this method (see Evidence Assumption 3).
    /// </summary>
    /// <param name="rearrangement">A breakpoint event produced by <see cref="DetectRearrangements"/>.</param>
    /// <returns>The rearrangement type implied by the breakpoint signature.</returns>
    public static RearrangementType ClassifyRearrangement(RearrangementEvent rearrangement)
    {
        // Re-derive (x, y) from TargetPosition "x->y" when available; otherwise trust the stored Type.
        if (rearrangement.TargetPosition is string tp)
        {
            int arrow = tp.IndexOf("->", StringComparison.Ordinal);
            if (arrow > 0 &&
                int.TryParse(tp.AsSpan(0, arrow), out int x) &&
                int.TryParse(tp.AsSpan(arrow + 2), out int y))
            {
                return ClassifyBoundary(x, y);
            }
        }

        return rearrangement.Type;
    }

    // --- CompareGenomes: pan-genome partition (Tettelin et al. 2005) + syntenic-gene fraction ---
    // Tettelin H, et al. (2005) PNAS 102(39):13950–13955, https://doi.org/10.1073/pnas.0506758102:
    // the species pan-genome = "a core genome containing genes present in all strains" plus
    // "a dispensable genome composed of genes absent from one or more strains and genes that are
    // unique to each strain." For two genomes: a gene that has an ortholog in the other genome is
    // CORE (conserved); a gene with no ortholog is DISPENSABLE (genome-specific). Shared genes are
    // found as reciprocal best hits (Moreno-Hagelsieb & Latimer 2008). OverallSynteny is the
    // "fraction of syntenic genes" (genes inside syntenic blocks ÷ smaller genome size), a standard
    // synteny-conservation metric, clamped to ≤ 1.

    /// <summary>Default RBH minimum similarity for the conserved-gene gate (inherited from FindReciprocalBestHits / COMPGEN-RBH-001).</summary>
    private const double DefaultCompareMinIdentity = DefaultMinIdentity;

    /// <summary>
    /// Default minimum collinear anchors per syntenic block used by CompareGenomes.
    /// Passed to FindSyntenicBlocks as minAnchors; note MCScanX still enforces the score ≥ 250
    /// (≥ 5 anchors) report rule, so OverallSynteny is non-zero only for ≥ 5 collinear orthologs.
    /// </summary>
    private const int DefaultCompareMinSyntenicBlockSize = 3;

    /// <summary>A fraction (e.g. OverallSynteny) cannot exceed 1.0.</summary>
    private const double MaxFraction = 1.0;

    /// <summary>
    /// Performs comprehensive two-genome comparison. Partitions genes into the <b>core (conserved)</b>
    /// set and each genome's <b>dispensable (genome-specific)</b> set following the pan-genome model of
    /// Tettelin et al. (2005, <i>PNAS</i> 102(39):13950–13955): a gene with an ortholog in the other
    /// genome is core; a gene with no ortholog is genome-specific. Shared genes are reciprocal best hits
    /// (Moreno-Hagelsieb &amp; Latimer 2008). <see cref="ComparisonResult.OverallSynteny"/> is the
    /// fraction of syntenic genes — genes inside syntenic blocks divided by the smaller genome's gene
    /// count, clamped to ≤ 1. Orthologs, syntenic blocks, and rearrangements are produced by the
    /// validated sub-methods (<see cref="FindReciprocalBestHits"/>, <see cref="FindSyntenicBlocks"/>,
    /// <see cref="DetectRearrangements"/>).
    /// </summary>
    /// <param name="genome1Genes">Genes of genome 1 in chromosomal order (each carrying a sequence for ortholog detection).</param>
    /// <param name="genome2Genes">Genes of genome 2 in chromosomal order.</param>
    /// <param name="minOrthologIdentity">Minimum RBH similarity for a conserved (shared) gene.</param>
    /// <param name="minSyntenicBlockSize">Minimum collinear anchors per syntenic block (see remarks on the MCScanX score threshold).</param>
    /// <returns>The comparison result: syntenic blocks, ortholog pairs, rearrangements, overall synteny, and the core/dispensable gene partition.</returns>
    /// <exception cref="ArgumentNullException">Either gene list is null.</exception>
    public static ComparisonResult CompareGenomes(
        IReadOnlyList<Gene> genome1Genes,
        IReadOnlyList<Gene> genome2Genes,
        double minOrthologIdentity = DefaultCompareMinIdentity,
        int minSyntenicBlockSize = DefaultCompareMinSyntenicBlockSize)
    {
        ArgumentNullException.ThrowIfNull(genome1Genes);
        ArgumentNullException.ThrowIfNull(genome2Genes);

        // Find orthologs (the shared/core genes) as reciprocal best hits.
        var orthologs = FindReciprocalBestHits(genome1Genes, genome2Genes, minOrthologIdentity).ToList();

        // Build ortholog map
        var orthologMap = orthologs.ToDictionary(o => o.Gene1Id, o => o.Gene2Id);

        // Find syntenic blocks
        var syntenicBlocks = FindSyntenicBlocks(
            genome1Genes, genome2Genes, orthologMap, minSyntenicBlockSize).ToList();

        // Detect rearrangements
        var rearrangements = DetectRearrangements(genome1Genes, genome2Genes, orthologMap).ToList();

        // Pan-genome partition (Tettelin et al. 2005): core = genes with an ortholog in the other
        // genome; dispensable/genome-specific = genes with no ortholog. The RBH matching maps each
        // gene at most once, so core + specific = genome size for each genome.
        var orthologGenes1 = new HashSet<string>(orthologs.Select(o => o.Gene1Id));
        var orthologGenes2 = new HashSet<string>(orthologs.Select(o => o.Gene2Id));

        int specific1 = genome1Genes.Count(g => !orthologGenes1.Contains(g.Id));
        int specific2 = genome2Genes.Count(g => !orthologGenes2.Contains(g.Id));

        // OverallSynteny = fraction of syntenic genes = genes inside syntenic blocks ÷ smaller genome.
        int smallerGenome = Math.Min(genome1Genes.Count, genome2Genes.Count);
        double synteny = syntenicBlocks.Count > 0 && smallerGenome > 0
            ? (double)syntenicBlocks.Sum(b => b.GeneCount) / smallerGenome
            : 0;

        return new ComparisonResult(
            SyntenicBlocks: syntenicBlocks,
            Orthologs: orthologs,
            Rearrangements: rearrangements,
            OverallSynteny: Math.Min(MaxFraction, synteny),
            ConservedGenes: orthologs.Count,
            GenomeSpecificGenes1: specific1,
            GenomeSpecificGenes2: specific2);
    }

    // --- Reversal distance: unsigned breakpoint lower bound (Bafna & Pevzner 1998) ---
    // The reversal distance is the minimum number of reversals transforming one gene order into
    // another. This method returns the BREAKPOINT LOWER BOUND on that distance for UNSIGNED
    // permutations, NOT the exact (signed) Hannenhalli–Pevzner distance:
    //   * A pair (π_i, π_{i+1}) of the extended permutation is a breakpoint iff π_{i+1} ≠ π_i + 1
    //     — unsigned form: |π_{i+1} − π_i| ≠ 1 (Bafna & Pevzner 1998, SIAM J. Discrete Math.
    //     11(2):224–240, §2; unsigned form per Hübotter 2020).
    //   * A single reversal removes at most two breakpoints, so b(π) ≤ 2t ⇒ d(π) ≥ b(π)/2
    //     (Hunter College CompBio Lecture 16).
    //   * The smallest integer satisfying d ≥ b/2 is ⌈b/2⌉, returned as (b + 1) / 2 below.

    // Breakpoint reduction cap: a reversal cuts two adjacencies, so it removes ≤ 2 breakpoints,
    // giving the divisor in d(π) ≥ b(π)/2 (Hunter College CompBio Lecture 16; Bafna & Pevzner 1998).
    private const int MaxBreakpointsRemovedPerReversal = 2;

    /// <summary>
    /// Calculates a lower bound on the reversal distance between two gene orders using the
    /// <b>unsigned breakpoint</b> model of Bafna &amp; Pevzner (1998, <i>SIAM J. Discrete Math.</i>
    /// 11(2):224–240, §2). Both inputs are treated as <b>unsigned</b> permutations of the same
    /// marker set; the result is <c>⌈b/2⌉</c> where <c>b</c> is the number of breakpoints of the
    /// extended relative permutation (a pair is a breakpoint iff the values are not consecutive
    /// integers). Because a reversal removes at most two breakpoints, this value is a guaranteed
    /// lower bound on the true reversal distance, NOT the exact signed Hannenhalli–Pevzner distance.
    /// </summary>
    /// <param name="permutation1">Source gene order (distinct integer marker ids).</param>
    /// <param name="permutation2">Target gene order over the same marker set.</param>
    /// <returns>Lower bound <c>⌈b/2⌉ ≥ 0</c> on the reversal distance.</returns>
    /// <exception cref="ArgumentException">The two orders have different lengths.</exception>
    public static int CalculateReversalDistance(
        IReadOnlyList<int> permutation1,
        IReadOnlyList<int> permutation2)
    {
        if (permutation1.Count != permutation2.Count)
            throw new ArgumentException("Permutations must have the same length");

        int n = permutation1.Count;
        // A permutation of fewer than two markers has no internal adjacency, hence no breakpoint.
        if (n <= 1) return 0;

        // Relabel to the relative permutation: target maps to the identity 0..n-1, so the
        // extended permutation is (-1, relative..., n) (Hunter Lecture 16: extended = (0, π, n+1)).
        var positionMap = new Dictionary<int, int>(n);
        for (int i = 0; i < n; i++)
            positionMap[permutation2[i]] = i;

        var relative = permutation1.Select(x => positionMap[x]).ToList();

        // Count breakpoints of the extended permutation. A pair is a breakpoint iff its two values
        // are not consecutive integers (unsigned definition, Bafna & Pevzner 1998 §2).
        int breakpoints = 0;

        // Left sentinel boundary: (-1, relative[0]) is a breakpoint iff relative[0] != 0.
        if (relative[0] != 0)
            breakpoints++;

        // Internal boundaries.
        for (int i = 0; i < n - 1; i++)
        {
            if (Math.Abs(relative[i + 1] - relative[i]) != 1)
                breakpoints++;
        }

        // Right sentinel boundary: (relative[n-1], n) is a breakpoint iff relative[n-1] != n-1.
        if (relative[n - 1] != n - 1)
            breakpoints++;

        // Lower bound d(π) ≥ b/2; smallest satisfying integer is ⌈b/2⌉ = (b + 1) / 2.
        return (breakpoints + 1) / MaxBreakpointsRemovedPerReversal;
    }

    // --- Conserved gene clusters: common intervals of permutations (Uno & Yagiura 1998/2000;
    //     Heber & Stoye 2001; common-interval definition per Bui-Xuan, Habib & Paul 2013) ---
    // Each genome is read as the sequence of ortholog-GROUP labels of its genes (in chromosomal
    // order). A CONSERVED GENE CLUSTER is a COMMON INTERVAL: "a set of integers that is an interval
    // of each Pk" (Bui-Xuan, Habib & Paul 2013, Def. 1, citing Uno & Yagiura) — i.e. a set of
    // ortholog groups that occupies a *contiguous window* (interval) in EVERY genome. An interval of
    // a permutation is the SET of ALL elements in some window [i,j] (Bui-Xuan et al. §2; Didier et
    // al. 2013 §2 for the sequence/duplicate case), so no foreign group may sit between members.
    // Intervals are defined only for i < j, so a cluster has size >= 2; the whole set is always a
    // (trivial) common interval. The strict model is gap-free; see Evidence Assumption 1 for maxGap.

    /// <summary>Smallest meaningful common interval has size 2 (interval defined only for i &lt; j; Bui-Xuan et al. 2013 §2).</summary>
    private const int MinCommonIntervalSize = 2;

    /// <summary>
    /// Finds conserved gene clusters across multiple genomes as <b>common intervals</b> of the
    /// ortholog-group permutations: a set of ortholog-group labels that occupies a contiguous window
    /// (interval) in <i>every</i> genome. This is the common-interval gene-cluster model of Uno &amp;
    /// Yagiura (2000, <i>Algorithmica</i> 26(2):290–309) and Heber &amp; Stoye (2001, CPM, LNCS 2089:207–218);
    /// formal definition per Bui-Xuan, Habib &amp; Paul (2013, arXiv:1304.5140, Def. 1): "a set of integers
    /// that is an interval of each Pₖ". Genes are mapped to their ortholog group via
    /// <paramref name="orthologGroups"/>; genes with no group are treated as window-breaking non-members.
    /// </summary>
    /// <param name="genomes">Genomes, each a gene list in chromosomal order.</param>
    /// <param name="orthologGroups">Map: gene id → ortholog-group id (the alphabet of group labels).</param>
    /// <param name="minClusterSize">Minimum number of distinct ortholog groups per reported cluster (≥ 2). Default 3.</param>
    /// <param name="maxGap">
    /// Retained for API/MCP backward compatibility. The validated behaviour is the strict (gap-free)
    /// common-interval model: clusters are contiguous windows in every genome. See Evidence Assumption 1.
    /// </param>
    /// <returns>The conserved clusters, each as the sorted list of its ortholog-group labels, in a deterministic order.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="genomes"/> or <paramref name="orthologGroups"/> is null.</exception>
    public static IEnumerable<IReadOnlyList<string>> FindConservedClusters(
        IReadOnlyList<IReadOnlyList<Gene>> genomes,
        IReadOnlyDictionary<string, string> orthologGroups,
        int minClusterSize = 3,
        int maxGap = 2)
    {
        ArgumentNullException.ThrowIfNull(genomes);
        ArgumentNullException.ThrowIfNull(orthologGroups);
        _ = maxGap; // strict gap-free model; parameter kept for signature/MCP compatibility (Evidence Assumption 1).

        // A common interval is a *family* notion (K ≥ 2 permutations): with fewer than two genomes
        // every window is trivially "common", so the conserved-cluster question is vacuous.
        if (genomes.Count < MinCommonIntervalSize)
            return Array.Empty<IReadOnlyList<string>>();

        int effectiveMin = Math.Max(minClusterSize, MinCommonIntervalSize);

        // Read each genome as its ordered sequence of ortholog-group labels.
        var groupSequences = new List<string[]>(genomes.Count);
        foreach (var genome in genomes)
        {
            var labels = new string[genome.Count];
            for (int i = 0; i < genome.Count; i++)
                labels[i] = orthologGroups.TryGetValue(genome[i].Id, out string? grp) ? grp : NoGroupSentinel;
            groupSequences.Add(labels);
        }

        // Candidate cluster sets = the group set of every contiguous window of the first genome that
        // contains no non-member sentinel. (Every common interval is, in particular, an interval of
        // genome 0, so enumerating genome-0 windows is complete.)
        var candidates = new HashSet<string>(StringComparer.Ordinal); // canonical key per distinct set
        var byKey = new Dictionary<string, SortedSet<string>>(StringComparer.Ordinal);

        string[] g0 = groupSequences[0];
        for (int start = 0; start < g0.Length; start++)
        {
            var set = new SortedSet<string>(StringComparer.Ordinal);
            for (int end = start; end < g0.Length; end++)
            {
                string label = g0[end];
                if (string.Equals(label, NoGroupSentinel, StringComparison.Ordinal))
                    break; // a window cannot contain a non-member group (interval = set of ALL window elements).

                set.Add(label);
                if (set.Count < effectiveMin)
                    continue;

                string key = string.Join("", set);
                if (candidates.Add(key))
                    byKey[key] = new SortedSet<string>(set, StringComparer.Ordinal);
            }
        }

        // Keep only candidate sets that are an interval in EVERY genome (common interval).
        var reported = new List<IReadOnlyList<string>>();
        foreach (var key in candidates)
        {
            var set = byKey[key];
            bool commonToAll = true;
            for (int g = 0; g < groupSequences.Count && commonToAll; g++)
                commonToAll = IsIntervalOf(groupSequences[g], set);

            if (commonToAll)
                reported.Add(set.ToList());
        }

        // Deterministic, order-independent output: by size, then lexicographically by joined labels.
        reported.Sort((a, b) =>
        {
            int c = a.Count.CompareTo(b.Count);
            if (c != 0) return c;
            return string.CompareOrdinal(string.Join("", a), string.Join("", b));
        });

        return reported;
    }

    /// <summary>Sentinel label for genes with no ortholog group; it can never be a cluster member and breaks windows.</summary>
    private const string NoGroupSentinel = "\0__NOGROUP__\0";

    /// <summary>
    /// Tests whether <paramref name="set"/> is an interval of the group sequence: there exists a
    /// contiguous window whose set of group labels equals exactly <paramref name="set"/> (Bui-Xuan,
    /// Habib &amp; Paul 2013 §2: an interval is the set of all elements of a window; Didier et al. 2013
    /// §2 for sequences with duplicates — any matching location suffices).
    /// </summary>
    private static bool IsIntervalOf(string[] sequence, SortedSet<string> set)
    {
        int target = set.Count;
        for (int start = 0; start < sequence.Length; start++)
        {
            if (!set.Contains(sequence[start]))
                continue;

            var window = new HashSet<string>(StringComparer.Ordinal);
            for (int end = start; end < sequence.Length; end++)
            {
                if (!set.Contains(sequence[end]))
                    break; // foreign element inside the window → not this location.

                window.Add(sequence[end]);
                if (window.Count == target)
                    return true; // window's set == set (no foreign elements, all members present).
            }
        }

        return false;
    }

    /// <summary>
    /// Calculates Average Nucleotide Identity (ANI) between two genomes.
    /// </summary>
    public static double CalculateANI(
        string genome1Sequence,
        string genome2Sequence,
        int fragmentSize = 1000,
        double minFragmentIdentity = 0.7)
    {
        if (string.IsNullOrEmpty(genome1Sequence) || string.IsNullOrEmpty(genome2Sequence))
            return 0;

        var identities = new List<double>();

        // Fragment genome1 and compare to genome2
        for (int i = 0; i <= genome1Sequence.Length - fragmentSize; i += fragmentSize / 2)
        {
            string fragment = genome1Sequence.Substring(i, fragmentSize);
            double bestIdentity = FindBestFragmentMatch(fragment, genome2Sequence);

            if (bestIdentity >= minFragmentIdentity)
            {
                identities.Add(bestIdentity);
            }
        }

        return identities.Count > 0 ? identities.Average() : 0;
    }

    private static double FindBestFragmentMatch(string fragment, string genome)
    {
        // Use SuffixTree for efficient longest common substring search
        var suffixTree = global::SuffixTree.SuffixTree.Build(genome.ToUpperInvariant());
        string lcs = suffixTree.LongestCommonSubstring(fragment.ToUpperInvariant());

        // Calculate identity based on LCS length relative to fragment length
        double identity = fragment.Length > 0 ? (double)lcs.Length / fragment.Length : 0;

        return Math.Min(identity, 1.0);
    }

    /// <summary>
    /// Generates a dot plot comparison between two sequences.
    /// Uses SuffixTree for efficient O(m+k) word matching.
    /// </summary>
    public static IEnumerable<(int x, int y)> GenerateDotPlot(
        string sequence1,
        string sequence2,
        int wordSize = 10,
        int stepSize = 1)
    {
        if (string.IsNullOrEmpty(sequence1) || string.IsNullOrEmpty(sequence2))
            yield break;

        // Build SuffixTree on sequence2 for efficient pattern matching
        var suffixTree = global::SuffixTree.SuffixTree.Build(sequence2.ToUpperInvariant());

        // Find matching words from sequence1 in sequence2
        for (int i = 0; i <= sequence1.Length - wordSize; i += stepSize)
        {
            string word = sequence1.Substring(i, wordSize).ToUpperInvariant();
            var positions = suffixTree.FindAllOccurrences(word);

            foreach (int j in positions)
            {
                yield return (i, j);
            }
        }
    }
}
