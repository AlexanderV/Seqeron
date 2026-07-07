namespace Seqeron.Genomics.Tests.Fuzzing;

/// <summary>
/// Fuzz tests for the PanGenome area — phylogenetic marker selection (PANGEN-MARKER-001),
/// the panX (Ding 2018) / Roary (Page 2015) practice of keeping <em>single-copy core</em>
/// gene clusters and scoring each by its number of <em>parsimony-informative sites</em>
/// (PIS), the column criterion of Zvelebil &amp; Baum (2008). The unit under test is the
/// pair
///   PanGenomeAnalyzer.CountParsimonyInformativeSites(IReadOnlyList&lt;string&gt;)
///   PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, coreClusters, totalGenomes, maxMarkers)
/// (src/Seqeron/Algorithms/Seqeron.Genomics.Metagenomics/PanGenomeAnalyzer.cs,
///  CountParsimonyInformativeSites lines 934–976; SelectPhylogeneticMarkers lines 997–1043).
///
/// ───────────────────────────────────────────────────────────────────────────
/// What fuzzing verifies
/// ───────────────────────────────────────────────────────────────────────────
/// Fuzzing feeds degenerate, boundary and out-of-domain inputs to the PIS counter and
/// the marker selector and asserts the code NEVER fails in an undisciplined way: no hang
/// or infinite loop, no unhandled runtime exception, no corruption of the descending-PIS
/// ordering or the single-copy-core filter. The headline hazards for THIS unit — the ones
/// the BE checklist row ("empty core, more markers than core") points at — are:
///   • an EMPTY core (no candidate clusters) ⇒ the selector must return an EMPTY marker
///     sequence, never throw and never emit a phantom marker (§3.3, §6.1 "null/empty
///     inputs");
///   • MORE markers than core (`maxMarkers` ≫ the number of qualifying clusters, up to
///     int.MaxValue) ⇒ the cap is a `Take`, so the result is simply ALL qualifying
///     markers — never an over-allocation, an IndexOutOfRange, or a throw (INV-06,
///     result size ≤ maxMarkers is trivially satisfied);
///   • degenerate alignments to the PIS counter: &lt; 2 rows, empty rows, null,
///     unequal-length rows ⇒ PIS 0 (no countable common alignment), never an
///     IndexOutOfRange on a short row or a null-deref (§3.3, §6.1).
/// Every input must resolve to EITHER a well-defined, theory-correct result (PIS in
/// [0, alignment length] per INV-01; only single-copy core clusters with ≥ 1 PIS
/// selected per INV-04/INV-05; descending PIS, ≤ maxMarkers per INV-06) OR the
/// documented degenerate result (PIS 0 / empty marker sequence for null / empty /
/// sub-spec input). A raw exception, a hang, a PIS outside [0, length], a selected
/// non-single-copy-core or 0-PIS cluster, a mis-ordered or over-long result is a bug,
/// not a passing test.
/// — docs/ADVANCED_TESTING_CHECKLIST.md §8 "Fuzzing".
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: PANGEN-MARKER-001 — Phylogenetic marker selection (single-copy core + PIS)
/// Checklist: docs/checklists/03_FUZZING.md, row 193.
/// Algorithm doc: docs/algorithms/PanGenome/Phylogenetic_Marker_Selection.md
///                (Test Unit ID PANGEN-MARKER-001).
/// Fuzz strategy exercised for THIS unit:
///   • BE = Boundary Exploitation — the degenerate boundaries called out in the row
///          ("empty core, more markers than core"), interpreted against the contract:
///            – EMPTY core: coreClusters is empty (or every candidate is filtered out by
///              the single-copy-core / 0-PIS rules) → an EMPTY marker sequence, never a
///              throw (§3.3, §6.1).
///            – MORE markers than core: maxMarkers far exceeds the number of qualifying
///              clusters (up to int.MaxValue) → ALL qualifying markers returned, in
///              descending-PIS order; the Take cap never over-runs (INV-06).
///            – EMPTY / sub-spec PIS inputs: 0/1 rows, empty rows, null rows, unequal
///              lengths → PIS 0 (§3.3, §6.1, INV-01).
/// — docs/checklists/03_FUZZING.md §Description (strategy codes:
///    "BE = Boundary Exploitation (граничні значення: 0, -1, MaxInt, empty)").
///
/// ───────────────────────────────────────────────────────────────────────────
/// The marker-selection contract under test (Phylogenetic_Marker_Selection.md §2–§3)
/// ───────────────────────────────────────────────────────────────────────────
/// A column of an alignment is PARSIMONY-INFORMATIVE iff it has ≥ 2 distinct character
/// states and each of ≥ 2 of those states occurs in ≥ 2 rows (Zvelebil &amp; Baum 2008,
/// §2.2). Monomorphic columns (one state) and singleton columns (a variant in only one
/// row) contribute 0 (INV-02). A cluster is selected as a MARKER iff (a) it is
/// single-copy core — GenomeCount == totalGenomes AND GeneIds.Count == totalGenomes
/// (panX "all strains represented exactly once"; Roary paralog filtering) (INV-04); and
/// (b) its member alignment has ≥ 1 PIS (INV-05). Selected markers are returned ordered
/// by DESCENDING PIS, ties by ordinal cluster id, capped at maxMarkers (INV-06).
///   • INV-01: 0 ≤ PIS ≤ alignment length.
///   • INV-07: PIS is invariant to row order and to bijective state relabeling.
///
/// HAND-CHECKABLE worked example (Phylogenetic_Marker_Selection.md §7.1): the alignment
///   s1 = "AAAAA", s2 = "AAACA", s3 = "AACCG", s4 = "ACCTG"
/// has column 3 = (A,A,C,C) [two states ×2] and column 5 = (A,A,G,G) [two states ×2]
/// parsimony-informative; column 1 is monomorphic, column 2 has a singleton C, column 4
/// is four singletons → PIS = 2. Verified directly below.
/// </summary>
[TestFixture]
[Category("Fuzzing")]
public class PanGenomeMarkerFuzzTests
{
    #region Helpers

    private const int Seed = 193_0001; // local-only deterministic seed for the fuzz sweep

    private static (string GeneId, string Sequence) Gene(string id, string seq) => (id, seq);

    /// <summary>Builds a genome → list of (geneId, sequence) dictionary.</summary>
    private static Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> Genomes(
        params (string Genome, (string GeneId, string Sequence)[] Genes)[] entries)
    {
        var dict = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
        foreach (var (genome, genes) in entries)
            dict[genome] = genes.ToList();
        return dict;
    }

    /// <summary>
    /// Constructs a <see cref="PanGenomeAnalyzer.GeneCluster"/> with explicit gene ids
    /// and genome membership. GenomeCount is derived from <paramref name="genomeIds"/> so
    /// a caller can build single-copy-core, paralog, or below-core clusters by choosing
    /// the gene/genome lists explicitly (mirrors the canonical unit-test helper).
    /// </summary>
    private static PanGenomeAnalyzer.GeneCluster Cluster(
        string id, string[] geneIds, string[] genomeIds) =>
        new(id, geneIds, genomeIds, genomeIds.Length, 1.0, geneIds.Length > 0 ? geneIds[0] : "");

    /// <summary>
    /// Asserts a marker-selection result is WELL-FORMED per Phylogenetic_Marker_Selection.md
    /// §2–§3: every returned cluster is single-copy core (INV-04), every returned cluster
    /// has ≥ 1 PIS over its member alignment (INV-05), the result is ordered by descending
    /// PIS (INV-06), and the result size is ≤ <paramref name="maxMarkers"/> (INV-06). Used
    /// by the sweep so a leaked non-core cluster, a 0-PIS marker, a mis-ordered result, or
    /// an over-long result fails everywhere.
    /// </summary>
    private static void AssertWellFormedMarkers(
        IReadOnlyList<PanGenomeAnalyzer.GeneCluster> markers,
        IReadOnlyDictionary<string, string> geneSequences,
        int totalGenomes,
        int maxMarkers)
    {
        markers.Count.Should().BeLessThanOrEqualTo(maxMarkers,
            "INV-06: the result is capped at maxMarkers via Take.");

        int prevPis = int.MaxValue;
        foreach (var m in markers)
        {
            // INV-04: single-copy core — present in all genomes, exactly one gene each.
            m.GenomeCount.Should().Be(totalGenomes,
                "INV-04: a selected marker is core (present in every genome).");
            m.GeneIds.Count.Should().Be(totalGenomes,
                "INV-04: a selected marker is single-copy (exactly one gene per genome).");

            var members = m.GeneIds
                .Select(g => geneSequences.TryGetValue(g, out var s) ? s : string.Empty)
                .ToList();
            int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(members);

            // INV-05: only variable (≥ 1 PIS) clusters are kept.
            pis.Should().BeGreaterThanOrEqualTo(1,
                "INV-05: every selected marker has at least one parsimony-informative site.");

            // INV-01: PIS never exceeds the alignment length and is never negative.
            int len = members.Count == 0 ? 0 : (members[0]?.Length ?? 0);
            pis.Should().BeInRange(0, len, "INV-01: 0 ≤ PIS ≤ alignment length.");

            // INV-06: descending PIS ordering.
            pis.Should().BeLessThanOrEqualTo(prevPis,
                "INV-06: markers are ordered by descending PIS.");
            prevPis = pis;
        }
    }

    private static Dictionary<string, string> GeneIndex(
        IReadOnlyDictionary<string, IReadOnlyList<(string GeneId, string Sequence)>> genomes)
    {
        var idx = new Dictionary<string, string>();
        foreach (var genes in genomes.Values)
            foreach (var (id, seq) in genes)
                idx[id] = seq ?? string.Empty;
        return idx;
    }

    #endregion

    #region PANGEN-MARKER-001 — Phylogenetic marker selection

    #region Positive sanity — documented worked example & a hand-checkable selection

    // Phylogenetic_Marker_Selection.md §7.1 worked example for the PIS counter, made
    // hand-checkable: alignment s1=AAAAA, s2=AAACA, s3=AACCG, s4=ACCTG. Col 3 (A,A,C,C)
    // and col 5 (A,A,G,G) are parsimony-informative; cols 1 (mono), 2 (singleton C),
    // 4 (four singletons) are not ⇒ PIS = 2.
    [Test]
    [CancelAfter(20000)]
    public void Sanity_PisCounter_DocWorkedExample_ReturnsTwo()
    {
        var aln = new[] { "AAAAA", "AAACA", "AACCG", "ACCTG" };

        int pis = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln);

        pis.Should().Be(2,
            "§7.1/INV-03: cols 3 (A,A,C,C) and 5 (A,A,G,G) have two states each in ≥2 rows; "
            + "col 1 monomorphic, col 2 singleton C, col 4 four singletons (Zvelebil 2008).");
    }

    // A full marker selection that is hand-checkable end-to-end. Two genomes; two
    // single-copy core clusters whose 2-row member alignments are informative, plus one
    // cluster that is NOT single-copy core (a paralog) and must be dropped (INV-04).
    //   markerHi: members "AC", "AG" → col 2 = (C,G) two singletons ⇒ NOT informative;
    //             so build informative rows instead: "AACC","AAGG"-style needs ≥2 of each
    //             state in a column. With only 2 rows a column is informative iff both
    //             rows share the same single state pattern — impossible (each state in 1
    //             row). So with exactly totalGenomes=2 NO 2-row cluster is informative.
    // Therefore use totalGenomes = 4 so member alignments have 4 rows and can be PI.
    [Test]
    [CancelAfter(20000)]
    public void Sanity_SelectMarkers_RanksInformativeCoreClustersByDescendingPis()
    {
        // Four genomes g1..g4. Two single-copy core clusters:
        //   cHi: members AAAAA, AAACA, AACCG, ACCTG → PIS = 2 (the §7.1 alignment).
        //   cLo: members ACAA, ACCA, ACGA, ACTA → only col 3 = (A,C,G,T) four singletons
        //        ⇒ NOT informative... build a 1-PIS cluster instead:
        //        members GGTT, GGTT, GGAA, GGAA → col 3 = (T,T,A,A) two states ×2 ⇒ PIS 1;
        //        col 4 = (T,T,A,A) ⇒ another PIS ⇒ that's 2. Use GGTA,GGTA,GGAA,GGAA:
        //        col 3 (T,T,A,A) PI=1, col 4 (A,A,A,A) mono ⇒ total PIS = 1.
        var genomes = Genomes(
            ("g1", new[] { Gene("hi1", "AAAAA"), Gene("lo1", "GGTA") }),
            ("g2", new[] { Gene("hi2", "AAACA"), Gene("lo2", "GGTA") }),
            ("g3", new[] { Gene("hi3", "AACCG"), Gene("lo3", "GGAA") }),
            ("g4", new[] { Gene("hi4", "ACCTG"), Gene("lo4", "GGAA") }));

        var cHi = Cluster("cHi", new[] { "hi1", "hi2", "hi3", "hi4" },
            new[] { "g1", "g2", "g3", "g4" });               // single-copy core, PIS = 2
        var cLo = Cluster("cLo", new[] { "lo1", "lo2", "lo3", "lo4" },
            new[] { "g1", "g2", "g3", "g4" });               // single-copy core, PIS = 1
        // A paralog: two genes from the same genome (g1) ⇒ GeneIds.Count != GenomeCount
        // distinctly, and only 3 genomes ⇒ not core. Must be dropped (INV-04).
        var cPar = Cluster("cPar", new[] { "hi1", "lo1", "hi2", "hi3" },
            new[] { "g1", "g2", "g3" });

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(
            genomes, new[] { cHi, cLo, cPar }, totalGenomes: 4, maxMarkers: 100).ToList();

        markers.Select(m => m.ClusterId).Should().Equal(new[] { "cHi", "cLo" },
            "INV-04/05/06: only the two single-copy core informative clusters are kept, "
            + "ordered by descending PIS (cHi PIS=2 before cLo PIS=1); the paralog/below-core "
            + "cluster is dropped.");

        AssertWellFormedMarkers(markers, GeneIndex(genomes), totalGenomes: 4, maxMarkers: 100);
    }

    #endregion

    #region BE — Boundary: empty core (no candidate clusters / all filtered out)

    // §3.3 / §6.1: an EMPTY core (no candidate clusters) ⇒ an EMPTY marker sequence,
    // never a throw. This is the headline "empty core" boundary for the row.
    [Test]
    [CancelAfter(20000)]
    public void EmptyCore_NoCandidateClusters_ReturnsEmpty_NoThrow()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("a", "AAAA") }),
            ("g2", new[] { Gene("b", "AACA") }));

        var act = () => PanGenomeAnalyzer.SelectPhylogeneticMarkers(
            genomes, Array.Empty<PanGenomeAnalyzer.GeneCluster>(), totalGenomes: 2).ToList();
        act.Should().NotThrow("§3.3: empty core ⇒ empty marker sequence, never a throw.");

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(
            genomes, Array.Empty<PanGenomeAnalyzer.GeneCluster>(), totalGenomes: 2).ToList();
        markers.Should().BeEmpty("§6.1: no candidate clusters ⇒ no markers.");
    }

    // Core present in NAME but EMPTY in EFFECT: every candidate is filtered out — some
    // are below-core, some are paralogs, some are single-copy core but fully conserved
    // (0 PIS). All must be dropped ⇒ empty result (INV-04/INV-05).
    [Test]
    [CancelAfter(20000)]
    public void EmptyCore_AllCandidatesFilteredOut_ReturnsEmpty()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("con1", "AAAA"), Gene("bel1", "ACGT"), Gene("par1", "AAAA"), Gene("par1b", "CCCC") }),
            ("g2", new[] { Gene("con2", "AAAA"), Gene("bel2", "ACGT"), Gene("par2", "AAAA") }),
            ("g3", new[] { Gene("con3", "AAAA"), Gene("bel3", "ACGT"), Gene("par3", "AAAA") }));

        // Fully conserved single-copy core (all "AAAA") ⇒ 0 PIS ⇒ dropped (INV-05).
        var conserved = Cluster("conserved", new[] { "con1", "con2", "con3" },
            new[] { "g1", "g2", "g3" });
        // Below-core: present in only 2 of 3 genomes ⇒ GenomeCount != totalGenomes (INV-04).
        var belowCore = Cluster("below", new[] { "bel1", "bel2" }, new[] { "g1", "g2" });
        // Paralog: g1 contributes 2 genes ⇒ GeneIds.Count (4) != totalGenomes (3) (INV-04).
        var paralog = Cluster("paralog", new[] { "par1", "par1b", "par2", "par3" },
            new[] { "g1", "g2", "g3" });

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(
            genomes, new[] { conserved, belowCore, paralog }, totalGenomes: 3).ToList();

        markers.Should().BeEmpty(
            "INV-04/INV-05: a fully conserved core gene (0 PIS), a below-core cluster, and a "
            + "paralog are all dropped ⇒ an empty (in effect) core yields no markers.");
    }

    #endregion

    #region BE — Boundary: more markers than core (maxMarkers ≫ qualifying clusters)

    // §3 INV-06: maxMarkers far larger than the number of qualifying clusters (including
    // int.MaxValue) ⇒ ALL qualifying markers are returned; the Take cap never over-runs,
    // over-allocates, or throws. This is the headline "more markers than core" boundary.
    [Test]
    [CancelAfter(20000)]
    public void MoreMarkersThanCore_MaxMarkersHuge_ReturnsAllQualifying_NoOverrun()
    {
        // One single-copy core informative cluster over 4 genomes (the §7.1 alignment).
        var genomes = Genomes(
            ("g1", new[] { Gene("m1", "AAAAA") }),
            ("g2", new[] { Gene("m2", "AAACA") }),
            ("g3", new[] { Gene("m3", "AACCG") }),
            ("g4", new[] { Gene("m4", "ACCTG") }));
        var core = Cluster("c", new[] { "m1", "m2", "m3", "m4" },
            new[] { "g1", "g2", "g3", "g4" });

        foreach (int maxMarkers in new[] { 5, 1000, int.MaxValue })
        {
            var act = () => PanGenomeAnalyzer.SelectPhylogeneticMarkers(
                genomes, new[] { core }, totalGenomes: 4, maxMarkers: maxMarkers).ToList();
            act.Should().NotThrow(
                $"maxMarkers = {maxMarkers} ≫ 1 qualifying cluster ⇒ Take returns all, never over-runs.");

            var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(
                genomes, new[] { core }, totalGenomes: 4, maxMarkers: maxMarkers).ToList();
            markers.Should().HaveCount(1,
                $"maxMarkers = {maxMarkers} ⇒ all 1 qualifying marker returned (INV-06, ≤ cap).");
            AssertWellFormedMarkers(markers, GeneIndex(genomes), totalGenomes: 4, maxMarkers: maxMarkers);
        }
    }

    // maxMarkers ≤ 0 is invalid input ⇒ empty sequence (§3.3), never a throw / negative Take.
    [Test]
    [CancelAfter(20000)]
    public void MoreMarkersThanCore_MaxMarkersZeroAndNegative_ReturnsEmpty()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("m1", "AAAAA") }),
            ("g2", new[] { Gene("m2", "AAACA") }),
            ("g3", new[] { Gene("m3", "AACCG") }),
            ("g4", new[] { Gene("m4", "ACCTG") }));
        var core = Cluster("c", new[] { "m1", "m2", "m3", "m4" },
            new[] { "g1", "g2", "g3", "g4" });

        foreach (int maxMarkers in new[] { 0, -1, int.MinValue })
        {
            var act = () => PanGenomeAnalyzer.SelectPhylogeneticMarkers(
                genomes, new[] { core }, totalGenomes: 4, maxMarkers: maxMarkers).ToList();
            act.Should().NotThrow($"maxMarkers = {maxMarkers} ≤ 0 ⇒ empty result, never a throw (§3.3).");

            PanGenomeAnalyzer.SelectPhylogeneticMarkers(
                genomes, new[] { core }, totalGenomes: 4, maxMarkers: maxMarkers)
                .Should().BeEmpty($"§3.3: maxMarkers = {maxMarkers} ≤ 0 ⇒ empty marker sequence.");
        }
    }

    // maxMarkers = 1 with several qualifying clusters ⇒ exactly the TOP (highest-PIS)
    // marker is returned (INV-06 cap interacts with the descending order).
    [Test]
    [CancelAfter(20000)]
    public void MoreMarkersThanCore_CapOfOne_ReturnsHighestPisMarker()
    {
        var genomes = Genomes(
            ("g1", new[] { Gene("hi1", "AAAAA"), Gene("lo1", "GGTA") }),
            ("g2", new[] { Gene("hi2", "AAACA"), Gene("lo2", "GGTA") }),
            ("g3", new[] { Gene("hi3", "AACCG"), Gene("lo3", "GGAA") }),
            ("g4", new[] { Gene("hi4", "ACCTG"), Gene("lo4", "GGAA") }));
        var cHi = Cluster("cHi", new[] { "hi1", "hi2", "hi3", "hi4" },
            new[] { "g1", "g2", "g3", "g4" }); // PIS = 2
        var cLo = Cluster("cLo", new[] { "lo1", "lo2", "lo3", "lo4" },
            new[] { "g1", "g2", "g3", "g4" }); // PIS = 1

        var markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(
            genomes, new[] { cLo, cHi }, totalGenomes: 4, maxMarkers: 1).ToList();

        markers.Should().ContainSingle().Which.ClusterId.Should().Be("cHi",
            "INV-06: a cap of 1 keeps the single highest-PIS marker (cHi PIS=2 > cLo PIS=1).");
    }

    #endregion

    #region BE — Boundary: degenerate PIS-counter inputs (empty / sub-spec alignments)

    // §3.3 / §6.1 "< 2 aligned rows": 0 or 1 rows ⇒ PIS 0 (no state can occur in ≥ 2 rows).
    [Test]
    [CancelAfter(20000)]
    public void PisCounter_FewerThanTwoRows_ReturnsZero()
    {
        PanGenomeAnalyzer.CountParsimonyInformativeSites(Array.Empty<string>())
            .Should().Be(0, "§6.1: 0 rows ⇒ no countable alignment ⇒ PIS 0.");
        PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "ACGT" })
            .Should().Be(0, "§6.1: 1 row ⇒ no state can occur in ≥ 2 rows ⇒ PIS 0.");
    }

    // §3.3: null input, empty rows, and unequal-length rows ⇒ PIS 0, never a null-deref
    // or an IndexOutOfRange on a short row.
    [Test]
    [CancelAfter(20000)]
    public void PisCounter_NullEmptyAndUnequalLength_ReturnsZero_NoThrow()
    {
        var actNull = () => PanGenomeAnalyzer.CountParsimonyInformativeSites(null!);
        actNull.Should().NotThrow("§3.3: null alignment ⇒ PIS 0, never a throw.");
        PanGenomeAnalyzer.CountParsimonyInformativeSites(null!).Should().Be(0);

        // Empty rows (length 0) ⇒ no columns ⇒ PIS 0.
        PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "", "" })
            .Should().Be(0, "§3.3: empty rows ⇒ no columns ⇒ PIS 0.");

        // Unequal-length rows ⇒ no common alignment ⇒ PIS 0 (Assumption 1), never an
        // IndexOutOfRange reading the longer row past the short row.
        var actUneq = () => PanGenomeAnalyzer.CountParsimonyInformativeSites(
            new[] { "ACGTACGT", "ACG", "ACGTACGT" });
        actUneq.Should().NotThrow("§3.3: unequal-length rows ⇒ PIS 0, never IndexOutOfRange.");
        PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "ACGTACGT", "ACG", "ACGTACGT" })
            .Should().Be(0, "§6.1: unequal lengths ⇒ no common alignment ⇒ PIS 0.");

        // A null ROW among otherwise valid rows ⇒ PIS 0, never a null-deref.
        var actNullRow = () => PanGenomeAnalyzer.CountParsimonyInformativeSites(
            new[] { "ACGT", null!, "ACGT" });
        actNullRow.Should().NotThrow("§3.3: a null row ⇒ PIS 0, never a null-deref.");
        PanGenomeAnalyzer.CountParsimonyInformativeSites(new[] { "ACGT", null!, "ACGT" })
            .Should().Be(0);
    }

    // INV-07: PIS is invariant to row order and to a bijective state relabeling. A fuzz
    // relation: shuffling rows and applying a state bijection (A↔T, C↔G) must not change
    // the count — a column-content property, not a row-position one.
    [Test]
    [CancelAfter(20000)]
    public void PisCounter_InvariantToRowOrderAndStateRelabeling()
    {
        var aln = new[] { "AAAAA", "AAACA", "AACCG", "ACCTG" };
        int baseline = PanGenomeAnalyzer.CountParsimonyInformativeSites(aln); // 2

        // Reversed row order.
        PanGenomeAnalyzer.CountParsimonyInformativeSites(aln.Reverse().ToList())
            .Should().Be(baseline, "INV-07: PIS is invariant to row order.");

        // Bijective state relabeling A↔T, C↔G (a permutation of the alphabet).
        string Relabel(string s) => new string(s.Select(c => c switch
        {
            'A' => 'T', 'T' => 'A', 'C' => 'G', 'G' => 'C', _ => c
        }).ToArray());
        PanGenomeAnalyzer.CountParsimonyInformativeSites(aln.Select(Relabel).ToList())
            .Should().Be(baseline, "INV-07: PIS is invariant to bijective state relabeling.");
    }

    #endregion

    #region BE — Boundary: degenerate selector inputs (null / totalGenomes ≤ 0)

    // §3.3: null genomes, null coreClusters, or totalGenomes ≤ 0 ⇒ empty marker sequence,
    // never a throw.
    [Test]
    [CancelAfter(20000)]
    public void Selector_NullAndNonPositiveTotalGenomes_ReturnEmpty_NoThrow()
    {
        var genomes = Genomes(("g1", new[] { Gene("a", "AAAA") }));
        var core = new[] { Cluster("c", new[] { "a" }, new[] { "g1" }) };

        var actNullGenomes = () => PanGenomeAnalyzer.SelectPhylogeneticMarkers(
            null!, core, totalGenomes: 1).ToList();
        actNullGenomes.Should().NotThrow("§3.3: null genomes ⇒ empty, never a throw.");
        PanGenomeAnalyzer.SelectPhylogeneticMarkers(null!, core, totalGenomes: 1)
            .Should().BeEmpty();

        var actNullCore = () => PanGenomeAnalyzer.SelectPhylogeneticMarkers(
            genomes, null!, totalGenomes: 1).ToList();
        actNullCore.Should().NotThrow("§3.3: null coreClusters ⇒ empty, never a throw.");
        PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, null!, totalGenomes: 1)
            .Should().BeEmpty();

        foreach (int tg in new[] { 0, -1, int.MinValue })
        {
            PanGenomeAnalyzer.SelectPhylogeneticMarkers(genomes, core, totalGenomes: tg)
                .Should().BeEmpty($"§3.3: totalGenomes = {tg} ≤ 0 ⇒ empty marker sequence.");
        }
    }

    #endregion

    #region Randomized sweep — never throw / leak a non-core or 0-PIS marker / mis-order

    // Broad randomized fuzz over genome counts, cluster mixes (single-copy core, paralog,
    // below-core, conserved), random sequences, and degenerate maxMarkers including the
    // boundaries (empty core, maxMarkers ≫ core, maxMarkers ≤ 0). The result must ALWAYS
    // be well-formed: only single-copy core ≥1-PIS clusters, descending PIS, ≤ maxMarkers —
    // never a throw, a hang, a leaked non-core/0-PIS marker, or a mis-ordered result.
    [Test]
    [CancelAfter(60000)]
    public void Fuzz_RandomCoreSets_AlwaysWellFormedMarkers()
    {
        var rng = new Random(Seed);
        char[] alphabet = { 'A', 'C', 'G', 'T' };

        for (int iter = 0; iter < 400; iter++)
        {
            int totalGenomes = rng.Next(1, 7);
            int alnLen = rng.Next(0, 8); // include 0-length alignments

            // Build a gene pool: every genome contributes some genes.
            var genomes = new Dictionary<string, IReadOnlyList<(string GeneId, string Sequence)>>();
            string Seq() => new string(Enumerable.Range(0, alnLen).Select(_ => alphabet[rng.Next(4)]).ToArray());
            for (int g = 0; g < totalGenomes; g++)
            {
                var genes = new List<(string GeneId, string Sequence)>();
                int geneCount = rng.Next(0, 4);
                for (int k = 0; k < geneCount; k++)
                    genes.Add(Gene($"g{g}_k{k}_i{iter}", Seq()));
                genomes[$"g{g}"] = genes;
            }
            var geneIndex = GeneIndex(genomes);
            var allGeneIds = geneIndex.Keys.ToList();

            // Build a random mix of candidate clusters (may be empty ⇒ "empty core").
            int clusterCount = rng.Next(0, 6);
            var clusters = new List<PanGenomeAnalyzer.GeneCluster>();
            for (int c = 0; c < clusterCount && allGeneIds.Count > 0; c++)
            {
                int genomeMembers = rng.Next(1, totalGenomes + 1);
                int geneMembers = rng.Next(1, Math.Max(2, allGeneIds.Count + 1));
                var geneIds = Enumerable.Range(0, geneMembers)
                    .Select(_ => allGeneIds[rng.Next(allGeneIds.Count)])
                    .Distinct()
                    .ToArray();
                var genomeIds = Enumerable.Range(0, genomeMembers)
                    .Select(j => $"g{j}")
                    .ToArray();
                clusters.Add(Cluster($"c{iter}_{c}", geneIds, genomeIds));
            }

            int maxMarkers = rng.Next(4) switch
            {
                0 => rng.Next(1, 3),     // small cap
                1 => 1_000_000,          // huge cap (more markers than core)
                2 => int.MaxValue,       // MaxInt cap
                _ => rng.Next(-2, 1),    // ≤ 0 (invalid)
            };

            List<PanGenomeAnalyzer.GeneCluster> markers = null!;
            var act = () => markers = PanGenomeAnalyzer.SelectPhylogeneticMarkers(
                genomes, clusters, totalGenomes, maxMarkers).ToList();
            act.Should().NotThrow(
                $"iter {iter}: genomes={totalGenomes}, clusters={clusterCount}, "
                + $"alnLen={alnLen}, maxMarkers={maxMarkers} must never throw.");

            if (maxMarkers <= 0)
            {
                markers.Should().BeEmpty($"iter {iter}: maxMarkers ≤ 0 ⇒ empty (§3.3).");
                continue;
            }

            AssertWellFormedMarkers(markers, geneIndex, totalGenomes, maxMarkers);
        }
    }

    #endregion

    #endregion
}
