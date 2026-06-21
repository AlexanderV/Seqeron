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
}
