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
}
