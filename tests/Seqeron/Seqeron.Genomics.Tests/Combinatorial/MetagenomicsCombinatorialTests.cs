namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Metagenomics area.
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of
/// combinatorial testing. Each grid cell carries a real business assertion;
/// small grids use the exhaustive <c>[Combinatorial]</c> product.
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Metagenomics")]
public class MetagenomicsCombinatorialTests
{
    private static string DiverseDna(int n, uint seed)
    {
        const string bases = "ACGT";
        var chars = new char[n];
        uint state = seed;
        for (int i = 0; i < n; i++)
        {
            state = state * 1664525u + 1013904223u;
            chars[i] = bases[(int)((state >> 16) & 3u)];
        }
        return new string(chars);
    }

    // root(1) → SpeciesA(2), SpeciesB(3)
    private static TaxonomyTree TwoSpecies() => new(new[]
    {
        new TaxonNode(1, "root", "root", 1),
        new TaxonNode(2, "SpeciesA", "species", 1),
        new TaxonNode(3, "SpeciesB", "species", 1),
    });

    private static readonly string GenomeA = DiverseDna(400, 0xA0A0u);
    private static readonly string GenomeB = DiverseDna(400, 0xB0B0u); // distinct → no shared k-mers

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-CLASS-001 — Taxonomic read classification (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 53.
    // Spec: tests/TestSpecs/META-CLASS-001.md (canonical ClassifyReads / BuildKmerDatabase).
    // Dimensions: kmerSize(3) × database(2) × readLen(3). Grid 3×2×3 = 18.
    //
    // Model (Kraken; Wood & Salzberg 2014): a read is classified by mapping its canonical k-mers
    // to taxa via the database and assigning the leaf of the maximum-scoring root-to-leaf path; a
    // read whose k-mers are absent from the database is left unclassified (root). A read drawn
    // wholly from one reference genome therefore classifies to that genome's taxon.
    //
    // The combinatorial point: k-mer length, the reference database content and read length
    // interact — a SpeciesA read is recovered for every k/length, while a SpeciesB read resolves
    // to SpeciesB only when the database includes SpeciesB and is otherwise unclassified.
    // ═══════════════════════════════════════════════════════════════════════

    public enum DatabaseKind { SpeciesAOnly, BothSpecies }

    [Test, Combinatorial]
    public void MetaClass_AssignsReadsByDatabaseContent(
        [Values(15, 21, 31)] int kmerSize,
        [Values(DatabaseKind.SpeciesAOnly, DatabaseKind.BothSpecies)] DatabaseKind database,
        [Values(50, 100, 200)] int readLen)
    {
        var taxonomy = TwoSpecies();
        var refs = database == DatabaseKind.BothSpecies
            ? new[] { (2, GenomeA), (3, GenomeB) }
            : new[] { (2, GenomeA) };
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(refs, taxonomy, kmerSize);

        string readA = GenomeA.Substring(50, readLen);
        var ca = MetagenomicsAnalyzer.ClassifyReads(new[] { ("rA", readA) }, db, taxonomy, kmerSize).Single();
        ca.TaxonId.Should().Be(2, "a read from genome A maps to SpeciesA");
        ca.TaxonName.Should().Be("SpeciesA");
        ca.TotalKmers.Should().Be(readLen - kmerSize + 1, "all-ACGT reads contribute every k-mer to Q");
        ca.MatchedKmers.Should().BeLessThanOrEqualTo(ca.TotalKmers);
        ca.Confidence.Should().BeInRange(0.0, 1.0);

        string readB = GenomeB.Substring(50, readLen);
        var cb = MetagenomicsAnalyzer.ClassifyReads(new[] { ("rB", readB) }, db, taxonomy, kmerSize).Single();
        cb.TaxonId.Should().Be(database == DatabaseKind.BothSpecies ? 3 : TaxonomyTree.RootId,
            "SpeciesB resolves only when the database contains it; otherwise unclassified");
    }

    /// <summary>
    /// Interaction witness: a read drawn wholly from one genome has every k-mer match its taxon,
    /// giving full confidence (1.0); a read with no database k-mers is unclassified with zero matches.
    /// </summary>
    [Test]
    public void MetaClass_FullMatchConfidence_AndUnclassified()
    {
        var taxonomy = TwoSpecies();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(new[] { (2, GenomeA) }, taxonomy, 21);

        var pure = MetagenomicsAnalyzer.ClassifyReads(new[] { ("a", GenomeA.Substring(10, 120)) }, db, taxonomy, 21).Single();
        pure.TaxonId.Should().Be(2);
        pure.MatchedKmers.Should().Be(pure.TotalKmers, "every k-mer of an A read is in the database");
        pure.Confidence.Should().BeApproximately(1.0, 1e-9);

        var alien = MetagenomicsAnalyzer.ClassifyReads(new[] { ("x", GenomeB.Substring(0, 120)) }, db, taxonomy, 21).Single();
        alien.TaxonId.Should().Be(TaxonomyTree.RootId, "no matching k-mers ⇒ unclassified");
        alien.MatchedKmers.Should().Be(0);
        alien.Confidence.Should().Be(0);
    }

    /// <summary>
    /// Interaction witness: ClassifyReads emits exactly one classification per input read, in order.
    /// </summary>
    [Test]
    public void MetaClass_OnePerRead_InOrder()
    {
        var taxonomy = TwoSpecies();
        var db = MetagenomicsAnalyzer.BuildKmerDatabase(new[] { (2, GenomeA), (3, GenomeB) }, taxonomy, 21);

        var reads = new[] { ("r0", GenomeA.Substring(0, 80)), ("r1", GenomeB.Substring(0, 80)), ("r2", GenomeA.Substring(100, 80)) };
        var results = MetagenomicsAnalyzer.ClassifyReads(reads, db, taxonomy, 21).ToList();

        results.Select(r => r.ReadId).Should().Equal("r0", "r1", "r2");
        results[0].TaxonId.Should().Be(2);
        results[1].TaxonId.Should().Be(3);
        results[2].TaxonId.Should().Be(2);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-PROF-001 — Taxonomic profile generation (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 54.
    // Spec: tests/TestSpecs/META-PROF-001.md (canonical GenerateTaxonomicProfile).
    // Dimensions: nReads(3) × nTaxa(3) × normalization(2). Grid 3×3×2 = 18.
    //
    // Model: a taxonomic profile turns per-read classifications into RELATIVE abundances
    // normalised over the CLASSIFIED reads (unclassified reads count toward TotalReads but not
    // the abundances), plus Shannon (−Σp·ln p) and Simpson (Σp²) indices over the species
    // proportions. The normalization axis varies whether unclassified "noise" reads are present.
    //
    // The combinatorial point: read count, taxon count and read composition interact — species
    // abundances always sum to 1 over classified reads and equal countₜ/classified, while
    // TotalReads/ClassifiedReads reflect the unclassified fraction; the diversity indices match
    // their closed-form definitions for every cell.
    // ═══════════════════════════════════════════════════════════════════════

    public enum ReadComposition { OnlyClassified, WithUnclassified }

    private static MetagenomicsAnalyzer.TaxonomicClassification Classified(string id, int t) =>
        new(id, 100 + t, $"Sp{t}", "species", 0, 1.0, 0, 0,
            "Bacteria", "Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", $"G{t}", $"Sp{t}");

    private static MetagenomicsAnalyzer.TaxonomicClassification Unclassified(string id) =>
        new(id, 1, "root", "root", 0, 0, 0, 0, "Unclassified", "", "", "", "", "", "");

    [Test, Combinatorial]
    public void MetaProf_RelativeAbundanceAndDiversity(
        [Values(20, 40, 80)] int nReads,
        [Values(2, 3, 5)] int nTaxa,
        [Values(ReadComposition.OnlyClassified, ReadComposition.WithUnclassified)] ReadComposition composition)
    {
        int unclassified = composition == ReadComposition.WithUnclassified ? nReads / 4 : 0;
        int classified = nReads - unclassified;
        var counts = new int[nTaxa];
        var list = new List<MetagenomicsAnalyzer.TaxonomicClassification>();
        for (int i = 0; i < classified; i++) { int t = i % nTaxa; counts[t]++; list.Add(Classified($"r{i}", t)); }
        for (int i = 0; i < unclassified; i++) list.Add(Unclassified($"u{i}"));

        var p = MetagenomicsAnalyzer.GenerateTaxonomicProfile(list);

        p.TotalReads.Should().Be(nReads);
        p.ClassifiedReads.Should().Be(classified);
        p.SpeciesAbundance.Should().HaveCount(nTaxa);
        p.SpeciesAbundance.Values.Sum().Should().BeApproximately(1.0, 1e-9, "abundances are relative to classified reads");

        for (int t = 0; t < nTaxa; t++)
            p.SpeciesAbundance[$"Sp{t}"].Should().BeApproximately((double)counts[t] / classified, 1e-9);

        var abund = p.SpeciesAbundance.Values.ToList();
        p.ShannonDiversity.Should().BeApproximately(-abund.Sum(a => a * Math.Log(a)), 1e-9);
        p.SimpsonDiversity.Should().BeApproximately(abund.Sum(a => a * a), 1e-9);
    }

    /// <summary>
    /// Interaction witness: an even community maximises Shannon diversity (ln S) and minimises
    /// Simpson dominance (1/S); a single-taxon community has Shannon 0 and Simpson 1.
    /// </summary>
    [Test]
    public void MetaProf_DiversityExtremes()
    {
        var even = Enumerable.Range(0, 4 * 5).Select(i => Classified($"r{i}", i % 4)).ToList();
        var pe = MetagenomicsAnalyzer.GenerateTaxonomicProfile(even);
        pe.ShannonDiversity.Should().BeApproximately(Math.Log(4), 1e-9, "four equal taxa ⇒ H = ln 4");
        pe.SimpsonDiversity.Should().BeApproximately(0.25, 1e-9, "even community ⇒ Σp² = 1/S");

        var mono = Enumerable.Range(0, 10).Select(i => Classified($"r{i}", 0)).ToList();
        var pm = MetagenomicsAnalyzer.GenerateTaxonomicProfile(mono);
        pm.ShannonDiversity.Should().BeApproximately(0.0, 1e-9, "one taxon ⇒ H = 0");
        pm.SimpsonDiversity.Should().BeApproximately(1.0, 1e-9, "one taxon ⇒ Σp² = 1");
    }

    /// <summary>
    /// Interaction witness: unclassified reads inflate TotalReads but not the abundances — the
    /// classified fraction and relative abundances are unchanged by adding noise.
    /// </summary>
    [Test]
    public void MetaProf_UnclassifiedReads_DoNotChangeRelativeAbundance()
    {
        var classifiedOnly = new[] { Classified("a", 0), Classified("b", 0), Classified("c", 1) };
        var withNoise = classifiedOnly.Append(Unclassified("n1")).Append(Unclassified("n2")).ToList();

        var p1 = MetagenomicsAnalyzer.GenerateTaxonomicProfile(classifiedOnly);
        var p2 = MetagenomicsAnalyzer.GenerateTaxonomicProfile(withNoise);

        p2.TotalReads.Should().Be(5);
        p2.ClassifiedReads.Should().Be(3);
        p2.SpeciesAbundance["Sp0"].Should().BeApproximately(p1.SpeciesAbundance["Sp0"], 1e-9, "abundance is over classified reads only");
        p2.SpeciesAbundance["Sp0"].Should().BeApproximately(2.0 / 3.0, 1e-9);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-ALPHA-001 — Alpha diversity indices (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 55.
    // Spec: tests/TestSpecs/META-ALPHA-001.md (canonical CalculateAlphaDiversity).
    // Dimensions: index(3: Shannon/Simpson/Chao1) × nSpecies(3) × evenness(3). Grid 3×3×3 = 27.
    //
    // Model: within-sample diversity — Shannon H = −Σp·ln p (Shannon-Wiener), Simpson D = Σp²
    // (Simpson 1949), and Chao1 = S_obs + f1²/(2f2) (or the bias-corrected form when f2=0)
    // estimating unseen richness from singletons/doubletons (Chao 1984).
    //
    // The combinatorial point: the index, the species count and the evenness of the abundance
    // distribution interact — each reported index equals its closed-form definition for every
    // community shape, and richness/evenness derivations stay consistent.
    // ═══════════════════════════════════════════════════════════════════════

    public enum AlphaIndex { Shannon, Simpson, Chao1 }
    public enum Evenness { Even, Moderate, Skewed }

    private static Dictionary<string, double> Community(int nSpecies, Evenness evenness)
    {
        var d = new Dictionary<string, double>();
        for (int i = 0; i < nSpecies; i++)
        {
            double count = evenness switch
            {
                Evenness.Even => 10,
                Evenness.Moderate => Math.Max(1, 10 - 2 * i),
                _ => i == 0 ? 50 : 1, // Skewed: one dominant, the rest singletons
            };
            d[$"sp{i}"] = count;
        }
        return d;
    }

    private static double IndepShannon(IReadOnlyCollection<double> c)
    {
        double sum = c.Sum();
        return -c.Where(v => v > 0).Sum(v => { double p = v / sum; return p * Math.Log(p); });
    }

    private static double IndepSimpson(IReadOnlyCollection<double> c)
    {
        double sum = c.Sum();
        return c.Where(v => v > 0).Sum(v => { double p = v / sum; return p * p; });
    }

    private static double IndepChao1(List<double> v)
    {
        int s = v.Count;
        int f1 = v.Count(x => Math.Abs(x - 1) < 1e-9), f2 = v.Count(x => Math.Abs(x - 2) < 1e-9);
        if (f1 == 0) return s;
        return f2 > 0 ? s + (double)(f1 * f1) / (2 * f2) : s + (double)(f1 * (f1 - 1)) / 2;
    }

    [Test, Combinatorial]
    public void MetaAlpha_IndicesMatchClosedForm(
        [Values(AlphaIndex.Shannon, AlphaIndex.Simpson, AlphaIndex.Chao1)] AlphaIndex index,
        [Values(3, 5, 8)] int nSpecies,
        [Values(Evenness.Even, Evenness.Moderate, Evenness.Skewed)] Evenness evenness)
    {
        var comm = Community(nSpecies, evenness);
        var a = MetagenomicsAnalyzer.CalculateAlphaDiversity(comm);
        var vals = comm.Values.ToList();

        a.ObservedSpecies.Should().Be(nSpecies);

        switch (index)
        {
            case AlphaIndex.Shannon:
                a.ShannonIndex.Should().BeApproximately(IndepShannon(vals), 1e-9);
                a.PielouEvenness.Should().BeApproximately(a.ShannonIndex / Math.Log(nSpecies), 1e-9);
                break;
            case AlphaIndex.Simpson:
                a.SimpsonIndex.Should().BeApproximately(IndepSimpson(vals), 1e-9);
                a.InverseSimpson.Should().BeApproximately(1 / IndepSimpson(vals), 1e-9);
                break;
            default:
                a.Chao1Estimate.Should().BeApproximately(IndepChao1(vals), 1e-9);
                a.Chao1Estimate.Should().BeGreaterThanOrEqualTo(nSpecies, "Chao1 never under-estimates observed richness");
                break;
        }
    }

    /// <summary>
    /// Interaction witness: an even community is maximally diverse — Pielou evenness 1.0 and
    /// Shannon ln(S) — whereas a skewed community is strictly less even.
    /// </summary>
    [Test]
    public void MetaAlpha_EvennessExtremes()
    {
        var even = MetagenomicsAnalyzer.CalculateAlphaDiversity(Community(5, Evenness.Even));
        even.PielouEvenness.Should().BeApproximately(1.0, 1e-9);
        even.ShannonIndex.Should().BeApproximately(Math.Log(5), 1e-9);

        var skewed = MetagenomicsAnalyzer.CalculateAlphaDiversity(Community(5, Evenness.Skewed));
        skewed.PielouEvenness.Should().BeLessThan(even.PielouEvenness, "a dominant taxon lowers evenness");
    }

    /// <summary>
    /// Interaction witness: Chao1 exceeds observed richness when singletons are present (unseen
    /// species are inferred) but equals it when every species is abundant.
    /// </summary>
    [Test]
    public void MetaAlpha_Chao1_ReflectsSingletons()
    {
        MetagenomicsAnalyzer.CalculateAlphaDiversity(Community(5, Evenness.Skewed))
            .Chao1Estimate.Should().BeGreaterThan(5, "singletons imply unseen richness");
        MetagenomicsAnalyzer.CalculateAlphaDiversity(Community(5, Evenness.Even))
            .Chao1Estimate.Should().Be(5, "no singletons ⇒ Chao1 = observed");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-BETA-001 — Beta diversity between samples (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 56.
    // Spec: tests/TestSpecs/META-BETA-001.md (canonical CalculateBetaDiversity).
    // Dimensions: metric(3: BC/Jaccard/UniFrac) × nSamples(3) × nSpecies(3). Grid 3×3×3 = 27.
    //
    // Model: between-sample dissimilarity — Bray-Curtis 1 − 2·Σmin(aᵢ,bᵢ)/Σ(aᵢ+bᵢ) (abundance)
    // and Jaccard 1 − shared/union (presence/absence), both in [0,1]. UniFrac requires a
    // phylogenetic tree and is not computed here (documented placeholder 0).
    //
    // The combinatorial point: the metric, sample count and species count interact — each
    // implemented dissimilarity equals its closed form for every consecutive sample pair (which
    // share all but their edge species), and the shared/unique counts are reproduced.
    // ═══════════════════════════════════════════════════════════════════════

    public enum BetaMetric { BrayCurtis, Jaccard, UniFrac }

    // Sliding species window: sample s covers species [s, s+nSpecies) so consecutive samples
    // share nSpecies−1 species and differ at one edge each.
    private static Dictionary<string, double> Sample(int s, int nSpecies) =>
        Enumerable.Range(0, nSpecies).ToDictionary(k => $"sp{s + k}", k => 10.0 + (s + k));

    private static double IndepBrayCurtis(Dictionary<string, double> a, Dictionary<string, double> b)
    {
        var union = a.Keys.Union(b.Keys);
        double sumMin = 0, sumTotal = 0;
        foreach (var sp in union)
        {
            double x = a.GetValueOrDefault(sp), y = b.GetValueOrDefault(sp);
            sumMin += Math.Min(x, y); sumTotal += x + y;
        }
        return sumTotal > 0 ? 1 - 2 * sumMin / sumTotal : 0;
    }

    private static (int Shared, int U1, int U2) Overlap(Dictionary<string, double> a, Dictionary<string, double> b)
    {
        int shared = a.Keys.Count(k => b.ContainsKey(k));
        return (shared, a.Count - shared, b.Count - shared);
    }

    [Test, Combinatorial]
    public void MetaBeta_DissimilaritiesMatchClosedForm(
        [Values(BetaMetric.BrayCurtis, BetaMetric.Jaccard, BetaMetric.UniFrac)] BetaMetric metric,
        [Values(2, 3, 4)] int nSamples,
        [Values(3, 5, 8)] int nSpecies)
    {
        var samples = Enumerable.Range(0, nSamples).Select(s => Sample(s, nSpecies)).ToArray();

        for (int i = 0; i + 1 < nSamples; i++)
        {
            var b = MetagenomicsAnalyzer.CalculateBetaDiversity($"S{i}", samples[i], $"S{i + 1}", samples[i + 1]);
            var (shared, u1, u2) = Overlap(samples[i], samples[i + 1]);

            switch (metric)
            {
                case BetaMetric.BrayCurtis:
                    b.BrayCurtis.Should().BeInRange(0.0, 1.0);
                    b.BrayCurtis.Should().BeApproximately(IndepBrayCurtis(samples[i], samples[i + 1]), 1e-9);
                    break;
                case BetaMetric.Jaccard:
                    b.JaccardDistance.Should().BeInRange(0.0, 1.0);
                    b.JaccardDistance.Should().BeApproximately(1.0 - (double)shared / (shared + u1 + u2), 1e-9);
                    b.SharedSpecies.Should().Be(shared);
                    b.UniqueToSample1.Should().Be(u1);
                    b.UniqueToSample2.Should().Be(u2);
                    break;
                default:
                    b.UniFracDistance.Should().Be(0, "UniFrac requires a phylogenetic tree and is not computed");
                    break;
            }
        }
    }

    /// <summary>
    /// Interaction witness: identical samples are maximally similar (BC = Jaccard = 0), while
    /// fully disjoint samples are maximally dissimilar (BC = Jaccard = 1).
    /// </summary>
    [Test]
    public void MetaBeta_IdenticalAndDisjointExtremes()
    {
        var s = Sample(0, 5);
        var same = MetagenomicsAnalyzer.CalculateBetaDiversity("a", s, "b", new Dictionary<string, double>(s));
        same.BrayCurtis.Should().BeApproximately(0.0, 1e-9);
        same.JaccardDistance.Should().BeApproximately(0.0, 1e-9);

        var disjoint = MetagenomicsAnalyzer.CalculateBetaDiversity("a", Sample(0, 3), "b", Sample(100, 3));
        disjoint.BrayCurtis.Should().BeApproximately(1.0, 1e-9);
        disjoint.JaccardDistance.Should().BeApproximately(1.0, 1e-9);
        disjoint.SharedSpecies.Should().Be(0);
    }

    /// <summary>
    /// Worked example: two samples sharing one of three species with differing abundances give
    /// Jaccard 1 − 1/3 and a Bray-Curtis matching the abundance formula.
    /// </summary>
    [Test]
    public void MetaBeta_WorkedExample()
    {
        var a = new Dictionary<string, double> { ["x"] = 4, ["y"] = 6 };
        var b = new Dictionary<string, double> { ["y"] = 2, ["z"] = 8 };
        var beta = MetagenomicsAnalyzer.CalculateBetaDiversity("a", a, "b", b);

        beta.SharedSpecies.Should().Be(1);   // only y
        beta.JaccardDistance.Should().BeApproximately(1.0 - 1.0 / 3.0, 1e-9);
        // sumMin = min(6,2)=2 ; sumTotal = (4)+(6+2)+(8) = 20 ; BC = 1 − 2·2/20 = 0.8
        beta.BrayCurtis.Should().BeApproximately(0.8, 1e-9);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-BIN-001 — Metagenome contig binning (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 57.
    // Spec: tests/TestSpecs/META-BIN-001.md (canonical BinContigs).
    // Dimensions: nContigs(3) × features(3: GC/tetra/coverage) × nBins(3). Grid 3×3×3 = 27.
    //
    // Model (MetaBAT-style binning): contigs are clustered (k-means) on composition+abundance
    // features — GC content, normalised coverage and tetranucleotide frequency — into genome
    // bins; bins below minBinSize are discarded. A bin's reported GC/coverage are the means of
    // its members; completeness is totalLength/expectedGenomeSize.
    //
    // The combinatorial point: contig count, the separating feature and the bin count interact;
    // across all of them the partition is disjoint, every bin clears minBinSize, and the bin
    // summaries equal the member averages. A separate witness verifies group recovery.
    // ═══════════════════════════════════════════════════════════════════════

    public enum BinFeature { Gc, Tetra, Coverage }

    private static double GcFraction(string s) => s.Count(c => c is 'G' or 'C') / (double)s.Length;

    private static List<(string ContigId, string Sequence, double Coverage)> ContigSet(int nGroups, int perGroup, BinFeature feature)
    {
        var tnfMotifs = new[] { "ACGT", "AGCT", "ATCG", "AGTC" }; // all 50% GC, distinct tetranucleotides
        var contigs = new List<(string, string, double)>();
        for (int g = 0; g < nGroups; g++)
            for (int i = 0; i < perGroup; i++)
            {
                string seq;
                double cov;
                switch (feature)
                {
                    case BinFeature.Gc:
                        int gcCount = (int)Math.Round(250 * (nGroups == 1 ? 0 : (double)g / (nGroups - 1)));
                        seq = new string('G', gcCount) + new string('A', 250 - gcCount); // GC = g/(K-1)
                        cov = 20; break;
                    case BinFeature.Coverage:
                        seq = DiverseDna(250, 0xC0u); cov = (g + 1) * 50; break;
                    default: // Tetra
                        seq = string.Concat(Enumerable.Repeat(tnfMotifs[g % tnfMotifs.Length], 64)); cov = 20; break;
                }
                contigs.Add(($"c{g}_{i}", seq, cov));
            }
        return contigs;
    }

    [Test, Combinatorial]
    public void MetaBin_PartitionIsConsistent_AcrossFeaturesAndBins(
        [Values(3, 5, 8)] int perGroup,
        [Values(BinFeature.Gc, BinFeature.Tetra, BinFeature.Coverage)] BinFeature feature,
        [Values(2, 3, 4)] int nBins)
    {
        var contigs = ContigSet(nBins, perGroup, feature);
        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: nBins, minBinSize: 1000, expectedGenomeSize: 100000).ToList();

        var binned = bins.SelectMany(b => b.ContigIds).ToList();
        binned.Should().OnlyHaveUniqueItems("bins are disjoint");
        binned.Count.Should().BeLessThanOrEqualTo(contigs.Count);

        foreach (var bin in bins)
        {
            bin.TotalLength.Should().BeGreaterThanOrEqualTo(1000, "bins below minBinSize are dropped");
            var ids = bin.ContigIds.ToHashSet();
            var members = contigs.Where(c => ids.Contains(c.ContigId)).ToList();

            bin.GcContent.Should().BeApproximately(members.Average(m => GcFraction(m.Sequence)), 1e-9, "bin GC is the member mean");
            bin.Coverage.Should().BeApproximately(members.Average(m => m.Coverage), 1e-9, "bin coverage is the member mean");
            bin.Completeness.Should().BeApproximately(Math.Min(bin.TotalLength / 100000.0 * 100, 100), 1e-9);
            bin.Contamination.Should().BeGreaterThanOrEqualTo(0);
        }
    }

    /// <summary>
    /// Interaction witness: two contig groups separated by coverage are recovered as two pure
    /// bins — every contig in a bin comes from the same source group.
    /// </summary>
    [Test]
    public void MetaBin_RecoversCoverageSeparatedGroups()
    {
        var contigs = ContigSet(nGroups: 2, perGroup: 5, BinFeature.Coverage);
        var bins = MetagenomicsAnalyzer.BinContigs(contigs, numBins: 2, minBinSize: 1000, expectedGenomeSize: 100000).ToList();

        bins.Should().HaveCount(2, "both genomes form a qualifying bin");
        foreach (var bin in bins)
        {
            var groups = bin.ContigIds.Select(id => id.Split('_')[0]).Distinct();
            groups.Should().ContainSingle("each bin holds contigs from a single source group");
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-FUNC-001 — Functional annotation by homology transfer (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 194.
    // Spec: tests/TestSpecs/META-FUNC-001.md (canonical PredictFunctions / FunctionalBitScore). ADVANCED §10.
    // Dimensions: dbSize(3) × nGenes(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Karlin-Altschul; BLAST): each protein is annotated by the best database signature it
    // contains, ranked by E-value = K·m·n·e^(−λ·S) (lowest wins); the bit score is
    // S' = (λ·S − ln K)/ln 2 over the ungapped BLOSUM62 self-score S.
    //
    // The combinatorial point: database size and the number of query genes interact — each gene's
    // annotation is its single contained signature, with E-value and bit score equal to the
    // Karlin-Altschul formulas (cross-checked via the exposed helpers), and decoy DB entries never hit.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly string[] FuncSignatures = { "WWWWWWWW", "YYYYYYYY", "FFFFFFFF" };
    private static readonly string[] FuncDecoys = { "MMMMMMMM", "PPPPPPPP", "EEEEEEEE", "QQQQQQQQ", "NNNNNNNN", "RRRRRRRR" };

    [Test, Combinatorial]
    public void MetaFunc_BestHitAnnotationAndScores_AcrossDbSizeAndGeneCount(
        [Values(3, 5, 7)] int dbSize,
        [Values(1, 2, 3)] int nGenes)
    {
        // Each gene contains exactly its own signature, padded with A.
        var proteins = Enumerable.Range(0, nGenes)
            .Select(i => ($"gene{i}", "AAAA" + FuncSignatures[i] + "AAAA"))
            .ToList();

        // DB = the nGenes real signatures + decoys to reach dbSize.
        var db = new Dictionary<string, (string Function, string Pathway, string Ko)>(StringComparer.Ordinal);
        for (int i = 0; i < nGenes; i++) db[FuncSignatures[i]] = ($"func{i}", $"path{i}", $"K0000{i}");
        int d = 0;
        while (db.Count < dbSize) { db[FuncDecoys[d]] = ($"decoy{d}", "-", "-"); d++; }

        var results = MetagenomicsAnalyzer.PredictFunctions(proteins, db).ToList();

        results.Should().HaveCount(nGenes, "every gene matches its embedded signature");
        foreach (var r in results)
        {
            int idx = int.Parse(r.GeneId.Substring(4));
            string sig = FuncSignatures[idx];
            r.Function.Should().Be($"func{idx}", "the best hit is the contained signature");
            int rawScore = MetagenomicsAnalyzer.Blosum62SelfScore(sig);
            r.BitScore.Should().BeApproximately(MetagenomicsAnalyzer.FunctionalBitScore(rawScore), 1e-9, "bit score = (λS−lnK)/ln2");
            r.EValue.Should().BeApproximately(MetagenomicsAnalyzer.ExpectedValue(rawScore, proteins[idx].Item2.Length, sig.Length), 1e-9,
                "E-value = K·m·n·e^(−λS)");
            r.EValue.Should().BeGreaterThan(0);
        }
    }

    /// <summary>
    /// Interaction witness — when a protein contains two signatures, the higher-raw-score one wins
    /// (lower E-value); bit score is monotone increasing in the raw score.
    /// </summary>
    [Test]
    public void MetaFunc_BestHitIsLowestEValue()
    {
        var db = new Dictionary<string, (string, string, string)>(StringComparer.Ordinal)
        {
            ["WWWWWWWW"] = ("tryptophan", "p", "k"),  // W BLOSUM diagonal 11 (higher self-score)
            ["CCCCCCCC"] = ("cysteine", "p", "k"),    // C BLOSUM diagonal 9 (lower self-score)
        };
        var protein = ("g", "WWWWWWWWAAAACCCCCCCC");

        MetagenomicsAnalyzer.PredictFunctions(new[] { protein }, db).Single().Function
            .Should().Be("tryptophan", "the higher-scoring signature gives the lower E-value");

        MetagenomicsAnalyzer.FunctionalBitScore(MetagenomicsAnalyzer.Blosum62SelfScore("WWWWWWWW"))
            .Should().BeGreaterThan(MetagenomicsAnalyzer.FunctionalBitScore(MetagenomicsAnalyzer.Blosum62SelfScore("CCCCCCCC")),
                "bit score increases with raw score");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-PATHWAY-001 — Pathway over-representation analysis (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 195.
    // Spec: tests/TestSpecs/META-PATHWAY-001.md (canonical FindPathwayEnrichment / HypergeometricUpperTail).
    // ADVANCED §10.
    // Dimensions: nGenes(3) × pathwaySet(2). Grid 3×2 = 6 (full, exhaustive ⊇ pairwise).
    //
    // Model (hypergeometric ORA): a pathway's enrichment p-value is the right-tail probability
    // P(X ≥ overlap) of drawing at least the observed overlap of query genes from the background,
    // P(X≥x) = Σ_{i≥x} C(M,i)C(N−M,n−i)/C(N,n).
    //
    // Axis mapping (documented): nGenes → query size; pathwaySet → which pathway database. The
    // combinatorial point: each reported pathway's p-value equals an independent hypergeometric
    // recomputation, results are sorted ascending by p-value, and the overlap/size bookkeeping is correct.
    // ═══════════════════════════════════════════════════════════════════════

    private static double Comb(int n, int k)
    {
        if (k < 0 || k > n) return 0;
        double r = 1;
        for (int i = 0; i < k; i++) r = r * (n - i) / (i + 1);
        return r;
    }

    private static double BruteHypergeometricUpperTail(int x, int bigN, int bigM, int n)
    {
        double p = 0;
        for (int i = x; i <= Math.Min(bigM, n); i++)
            p += Comb(bigM, i) * Comb(bigN - bigM, n - i) / Comb(bigN, n);
        return p;
    }

    [Test, Combinatorial]
    public void MetaPathway_EnrichmentMatchesHypergeometric_AcrossQuerySizeAndDb(
        [Values(2, 3, 4)] int nGenes,
        [Values(0, 1)] int pathwaySet)
    {
        var db = pathwaySet == 0
            ? new Dictionary<string, IReadOnlyCollection<string>>
            {
                ["P1"] = new[] { "g1", "g2", "g3", "g4" },
                ["P2"] = new[] { "g5", "g6", "g7", "g8" },
            }
            : new Dictionary<string, IReadOnlyCollection<string>>
            {
                ["A"] = new[] { "g1", "g2", "g3" },
                ["B"] = new[] { "g3", "g4", "g5", "g6" },
                ["C"] = new[] { "g7", "g8" },
            };
        var background = Enumerable.Range(1, 10).Select(i => $"g{i}").ToArray();
        var query = Enumerable.Range(1, nGenes).Select(i => $"g{i}").ToArray();

        var results = MetagenomicsAnalyzer.FindPathwayEnrichment(query, db, background);

        results.Select(r => r.PValue).Should().BeInAscendingOrder("pathways are sorted by p-value");
        foreach (var r in results)
        {
            r.PValue.Should().BeApproximately(
                MetagenomicsAnalyzer.HypergeometricUpperTail(r.Overlap, r.BackgroundSize, r.PathwaySize, r.QuerySize), 1e-9,
                "the reported p-value is the hypergeometric upper tail");
            r.PValue.Should().BeInRange(0.0, 1.0 + 1e-9);
        }
    }

    /// <summary>
    /// Interaction witness — HypergeometricUpperTail matches an independent C(·) summation, and a
    /// full overlap is more significant (smaller p) than a partial one.
    /// </summary>
    [Test]
    public void MetaPathway_HypergeometricUpperTailIsExact()
    {
        foreach (var (x, n, m, k) in new[] { (2, 10, 4, 4), (3, 20, 5, 6), (1, 8, 2, 3) })
            MetagenomicsAnalyzer.HypergeometricUpperTail(x, n, m, k)
                .Should().BeApproximately(BruteHypergeometricUpperTail(x, n, m, k), 1e-9);

        MetagenomicsAnalyzer.HypergeometricUpperTail(4, 20, 4, 4)
            .Should().BeLessThan(MetagenomicsAnalyzer.HypergeometricUpperTail(2, 20, 4, 4), "a larger overlap is more significant");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-RESIST-001 — Antibiotic-resistance gene detection (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 196.
    // Spec: tests/TestSpecs/META-RESIST-001.md (canonical FindAntibioticResistanceGenes). ADVANCED §10.
    // Dimensions: dbSize(3) × identity(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (ResFinder; Zankari 2012): a contig's best ungapped match to each reference gene gives a
    // percent identity (identical positions / alignment length) and reference coverage; a gene is
    // reported only when both clear their thresholds, with the single best-matching reference per contig.
    //
    // Axis mapping (documented): dbSize → reference-database size; identity → the identity threshold.
    // Engineered construct: a contig embedding a reference copy with interior mismatches ⇒ 85% identity,
    // full coverage. The combinatorial point: the resistance gene is reported exactly when the identity
    // threshold ≤ 0.85, independent of how many decoy references are present.
    // ═══════════════════════════════════════════════════════════════════════

    private const string ResistRef = "ACGTACGTACGTACGTACGT"; // 20 nt reference resistance gene

    [Test, Combinatorial]
    public void MetaResist_DetectsResistanceGene_AcrossDbSizeAndIdentity(
        [Values(1, 3, 5)] int dbSize,
        [Values(0.8, 0.9, 1.0)] double identity)
    {
        // Embed a copy with 3 interior mismatches (positions 5,10,15) ⇒ 17/20 = 0.85 identity, full coverage.
        var variant = ResistRef.ToCharArray();
        variant[5] = 'A'; variant[10] = 'A'; variant[15] = 'A';
        string contig = "TTTT" + new string(variant) + "TTTT";

        var refs = new List<(string, string, string, string)> { ("arr1", ResistRef, "ARR-1", "beta-lactam") };
        string[] decoys = { "GGGGCCCCGGGGCCCCGGGG", "TTAATTAATTAATTAATTAA", "CCCCGGGGCCCCGGGGCCCC", "AATTAATTAATTAATTAATT" };
        int d = 0;
        while (refs.Count < dbSize) { refs.Add(($"decoy{d}", decoys[d], $"DCY{d}", "none")); d++; }

        var hits = MetagenomicsAnalyzer.FindAntibioticResistanceGenes(
            new[] { ("contig1", contig) }, refs, identityThreshold: identity).ToList();

        bool expectHit = identity <= 0.85 + 1e-9;
        if (expectHit)
        {
            var hit = hits.Should().ContainSingle().Subject;
            hit.ResistanceGene.Should().Be("ARR-1");
            hit.AntibioticClass.Should().Be("beta-lactam");
            hit.PercentIdentity.Should().BeApproximately(0.85, 1e-9);
            hit.Coverage.Should().BeApproximately(1.0, 1e-9, "the full reference is covered");
        }
        else
        {
            hits.Should().BeEmpty("85% identity is below the threshold");
        }
    }

    /// <summary>
    /// Interaction witness — an exact reference copy is detected at every threshold, while an
    /// unrelated contig yields no hit.
    /// </summary>
    [Test]
    public void MetaResist_ExactCopyAlwaysDetected_UnrelatedNever()
    {
        var refs = new[] { ("arr1", ResistRef, "ARR-1", "beta-lactam") };

        MetagenomicsAnalyzer.FindAntibioticResistanceGenes(new[] { ("c", "GG" + ResistRef + "GG") }, refs, 1.0)
            .Should().ContainSingle(h => h.PercentIdentity >= 1.0 - 1e-9, "an exact copy is 100% identity");

        MetagenomicsAnalyzer.FindAntibioticResistanceGenes(new[] { ("c", new string('G', 30)) }, refs, 0.8)
            .Should().BeEmpty("an unrelated contig has no passing match");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-TAXA-001 — Differential-abundance testing (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 197.
    // Spec: tests/TestSpecs/META-TAXA-001.md (canonical FindSignificantTaxa / MannWhitneyU). ADVANCED §10.
    // Dimensions: nSamples(3) × nTaxa(3) × test(2). Grid 3×3×2 = 18 (full, exhaustive).
    //
    // Model (Mann & Whitney 1947): per taxon, a two-group rank-sum test with the asymptotic normal
    // approximation; the two-tailed p-value flags significance. The continuity correction subtracts
    // 0.5 from |U − m_U|, enlarging the p-value.
    //
    // Axis mapping (documented): nSamples → samples per group; nTaxa → number of taxa; test → the
    // continuity-correction flag. The combinatorial point: FindSignificantTaxa delegates per taxon to
    // MannWhitneyU (p-values agree), flags significance at p < threshold, and the U statistic is the
    // larger of U1/U2 — at every (nSamples, nTaxa, correction) cell.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void MetaTaxa_DelegatesToMannWhitneyPerTaxon_AcrossSamplesTaxaCorrection(
        [Values(3, 4, 5)] int nSamples,
        [Values(2, 3, 4)] int nTaxa,
        [Values(true, false)] bool continuityCorrection)
    {
        var profiles = new List<IReadOnlyDictionary<string, double>>();
        var groups = new List<int>();
        for (int s = 0; s < 2 * nSamples; s++)
        {
            int group = s < nSamples ? 1 : 2;
            var profile = new Dictionary<string, double>();
            for (int t = 0; t < nTaxa; t++)
            {
                // taxon0 is strongly differential; the rest overlap between groups.
                profile[$"taxon{t}"] = t == 0
                    ? (group == 1 ? s : 100.0 + s)
                    : (s % 3) + t * 0.1;
            }
            profiles.Add(profile);
            groups.Add(group);
        }

        var results = MetagenomicsAnalyzer.FindSignificantTaxa(profiles, groups, 0.05, continuityCorrection);

        results.Should().HaveCount(nTaxa, "one result per taxon");
        foreach (var st in results)
        {
            var g1 = Enumerable.Range(0, nSamples).Select(i => profiles[i].GetValueOrDefault(st.Taxon, 0.0)).ToList();
            var g2 = Enumerable.Range(nSamples, nSamples).Select(i => profiles[i].GetValueOrDefault(st.Taxon, 0.0)).ToList();
            var mw = MetagenomicsAnalyzer.MannWhitneyU(g1, g2, continuityCorrection);

            st.PValue.Should().BeApproximately(mw.PValue, 1e-12, "the taxon p-value is the Mann–Whitney p-value");
            st.U.Should().BeApproximately(Math.Max(mw.U1, mw.U2), 1e-12, "U is the larger of U1/U2");
            st.Significant.Should().Be(st.PValue < 0.05, "significance is p < threshold");
        }
    }

    /// <summary>
    /// Interaction witness — the Mann–Whitney U on a fully separated pair is exact, and the continuity
    /// correction enlarges the p-value (test axis).
    /// </summary>
    [Test]
    public void MetaTaxa_MannWhitneyExact_AndContinuityCorrectionRaisesP()
    {
        var g1 = new double[] { 1, 2, 3 };
        var g2 = new double[] { 4, 5, 6 };
        var mw = MetagenomicsAnalyzer.MannWhitneyU(g1, g2, useContinuityCorrection: false);
        mw.U1.Should().Be(0, "group1 ranks 1+2+3 ⇒ U1 = 6 − 6 = 0");
        mw.U2.Should().Be(9, "U2 = n1·n2 − U1 = 9");

        double pCc = MetagenomicsAnalyzer.MannWhitneyU(g1, g2, true).PValue;
        double pNoCc = MetagenomicsAnalyzer.MannWhitneyU(g1, g2, false).PValue;
        pCc.Should().BeGreaterThanOrEqualTo(pNoCc, "the continuity correction enlarges the p-value");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-CHECKM-001 — CheckM marker-gene bin quality (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 248.
    // Spec: tests/TestSpecs/META-CHECKM-001.md (MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts). ADVANCED §10.
    //
    // Sources: Parks et al. (2015) CheckM, Genome Res. 25:1043 (Eqs. 1–2).
    //
    // Model: completeness = 100·(1/|M|)·Σ_s present_s/|s|;  contamination = 100·(1/|M|)·Σ_s Σ_{g∈s}(N_g−1)/|s|
    //   over a domain's collocated marker sets M.
    //
    // Dimensions: domain(bac/ar) × markerCount(3) × copies(2). Grid 2×3×2 = 12 (exhaustive).
    //
    // The combinatorial point: domain selects the bundled marker-set family; markerCount (the number of
    // distinct markers present) drives completeness; copies (single vs duplicated) drives contamination.
    // Every cell reproduces BOTH CheckM equations from an independent recomputation over the bundled sets.
    // ═══════════════════════════════════════════════════════════════════════

    public enum CheckMDomain { Bacterial, Archaeal }

    private static (double Completeness, double Contamination) CheckMGroundTruth(
        IReadOnlyList<MetagenomicsAnalyzer.MarkerSet> sets, IReadOnlyDictionary<string, int> counts)
    {
        double compSum = 0.0, contSum = 0.0;
        foreach (var set in sets)
        {
            var ids = set.MarkerIds;
            int present = ids.Count(id => counts.TryGetValue(id, out int n) && n >= 1);
            double multi = ids.Sum(id => counts.TryGetValue(id, out int n) && n > 1 ? n - 1 : 0);
            compSum += (double)present / ids.Count;
            contSum += multi / ids.Count;
        }
        return (100.0 * compSum / sets.Count, 100.0 * contSum / sets.Count);
    }

    [Test, Combinatorial]
    public void CheckM_ReproducesCompletenessAndContamination_AcrossDomainMarkerCountAndCopies(
        [Values(CheckMDomain.Bacterial, CheckMDomain.Archaeal)] CheckMDomain domain,
        [Values(6, 12, 18)] int markerCount,
        [Values(1, 2)] int copies)
    {
        var sets = domain == CheckMDomain.Bacterial
            ? MetagenomicsAnalyzer.BundledBacterialMarkerSets()
            : MetagenomicsAnalyzer.BundledArchaealMarkerSets();

        var allMarkers = sets.SelectMany(s => s.MarkerIds).Distinct().ToList();
        int present = Math.Min(markerCount, allMarkers.Count);
        var counts = allMarkers.Take(present).ToDictionary(m => m, _ => copies);

        var q = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, counts);
        var (expComp, expCont) = CheckMGroundTruth(sets, counts);

        q.Completeness.Should().BeApproximately(expComp, 1e-9, "completeness reproduces CheckM Eq. 1 over the bundled sets");
        q.Contamination.Should().BeApproximately(expCont, 1e-9, "contamination reproduces CheckM Eq. 2 over the bundled sets");
        if (copies == 1)
            q.Contamination.Should().BeApproximately(0.0, 1e-9, "single-copy markers carry no contamination");
        else
            q.Contamination.Should().BeGreaterThan(0.0, "duplicated markers raise contamination");
    }

    /// <summary>
    /// Interaction witness: a complete single-copy set scores 100/0, and duplicating every marker raises
    /// contamination while leaving completeness at 100.
    /// </summary>
    [Test]
    public void CheckM_CompleteSingleCopyIs100And0_DuplicationRaisesContamination()
    {
        var sets = MetagenomicsAnalyzer.BundledBacterialMarkerSets();
        var allMarkers = sets.SelectMany(s => s.MarkerIds).Distinct().ToList();

        var single = allMarkers.ToDictionary(m => m, _ => 1);
        var dup = allMarkers.ToDictionary(m => m, _ => 2);

        var qSingle = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, single);
        var qDup = MetagenomicsAnalyzer.EstimateBinQualityFromMarkerCounts(sets, dup);

        qSingle.Completeness.Should().BeApproximately(100.0, 1e-9);
        qSingle.Contamination.Should().BeApproximately(0.0, 1e-9);
        qDup.Completeness.Should().BeApproximately(100.0, 1e-9, "every marker is still present");
        qDup.Contamination.Should().BeGreaterThan(qSingle.Contamination, "duplication raises contamination");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: META-TETRA-001 — Tetranucleotide z-score composition vector (Metagenomics)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 249.
    // Spec: tests/TestSpecs/META-TETRA-001.md
    //       (MetagenomicsAnalyzer.CalculateTetranucleotideZScores / TetranucleotideZScoreCorrelation). ADVANCED §10.
    //
    // Sources: Teeling et al. (2004) TETRA, BMC Bioinformatics 5:163 (strand-symmetric tetramer z-scores).
    //
    // Model: each tetramer's z-score is observed − expected over the Markov standard deviation; the
    //   estimator folds the sequence with its reverse complement, so the TETRA vector is STRAND-SYMMETRIC.
    //   The z-score correlation is the binning similarity metric (1 for identical composition).
    //
    // Dimensions: seqLen(3) × GCcontent(3) × strand(2). Grid 3×3×2 = 18 (exhaustive).
    //
    // The combinatorial point: across length and GC content the z-score vector is well-defined and
    // self-consistent (self-correlation 1). The strand axis reflects TETRA's strand symmetry — a
    // sequence correlates far more strongly with its OWN reverse complement than with a compositionally
    // different sequence (the basis of strand-independent binning).
    // ═══════════════════════════════════════════════════════════════════════

    private static char MetaComp(char b) => b switch { 'A' => 'T', 'T' => 'A', 'G' => 'C', 'C' => 'G', _ => b };
    private static string MetaRevComp(string s) => new(s.Reverse().Select(MetaComp).ToArray());

    // Non-periodic sequence of the requested length whose composition realises the GC level (LCG-driven,
    // all 4 bases present) — avoids the degenerate all-zero z-vector of a strictly periodic sequence.
    private static string TetraSeq(int len, double gc)
    {
        var chars = new char[len];
        uint state = (uint)(len * 2654435761u) ^ (uint)(gc * 1000.0);
        for (int i = 0; i < len; i++)
        {
            state = state * 1664525u + 1013904223u;
            double u = ((state >> 8) & 0xFFFF) / 65536.0;
            uint pick = (state >> 24) & 1u;
            chars[i] = u < gc ? (pick == 0 ? 'G' : 'C') : (pick == 0 ? 'A' : 'T');
        }
        return new string(chars);
    }

    public enum TetraStrand { Forward, Reverse }

    [Test, Combinatorial]
    public void Tetra_SelfConsistentAndStrandSymmetric_AcrossLengthGcAndStrand(
        [Values(120, 240, 480)] int seqLen,
        [Values(0.3, 0.5, 0.7)] double gc,
        [Values(TetraStrand.Forward, TetraStrand.Reverse)] TetraStrand strand)
    {
        string fwd = TetraSeq(seqLen, gc);
        string seq = strand == TetraStrand.Forward ? fwd : MetaRevComp(fwd);

        var z = MetagenomicsAnalyzer.CalculateTetranucleotideZScores(seq);
        z.Should().NotBeEmpty("a tetranucleotide z-score vector is defined");

        MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, seq)
            .Should().BeApproximately(1.0, 1e-9, "the z-score vector correlates perfectly with itself");

        // Strand symmetry: a sequence resembles its own reverse complement far more than a
        // compositionally different sequence of the same length.
        double selfRc = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, MetaRevComp(seq));
        double crossGc = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(seq, TetraSeq(seqLen, gc < 0.5 ? 0.8 : 0.2));
        selfRc.Should().BeGreaterThan(crossGc, "TETRA is strand-symmetric: revcomp is closer than a different-GC sequence");
    }

    /// <summary>
    /// Interaction witness (GC axis): two sequences of very different GC content have a TETRA z-score
    /// correlation below 1 — the metric discriminates composition (the basis of TETRA binning).
    /// </summary>
    [Test]
    public void Tetra_DifferentComposition_CorrelationBelowOne()
    {
        double corr = MetagenomicsAnalyzer.TetranucleotideZScoreCorrelation(
            TetraSeq(480, 0.3), TetraSeq(480, 0.75));
        corr.Should().BeLessThan(1.0, "compositionally distinct sequences are not perfectly correlated");
    }
}
