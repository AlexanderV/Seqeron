namespace Seqeron.Genomics.Tests.Combinatorial;

/// <summary>
/// Combinatorial (pairwise / full-grid) tests for the Variants area
/// (VariantAnnotator + VariantCaller, both in Seqeron.Genomics.Annotation).
///
/// See <see cref="CompositionCombinatorialTests"/> for the rationale of combinatorial testing.
/// Each grid cell carries a real business assertion; small grids use the exhaustive
/// <c>[Combinatorial]</c> product (a strict superset of pairwise).
/// — docs/checklists/09_COMBINATORIAL_TESTING.md §Description.
/// </summary>
[TestFixture]
[Category("Combinatorial")]
[Category("Variants")]
public class VariantsCombinatorialTests
{
    // ═══════════════════════════════════════════════════════════════════════
    // Unit: VARIANT-ANNOT-001 — Functional-impact annotation (Variants)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 186.
    // Spec: tests/TestSpecs/VARIANT-ANNOT-001.md (canonical PredictFunctionalImpact). ADVANCED §10.
    // Dimensions: consequence(3) × region(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model (Ensembl VEP; VariationEffect.pm predicates; Constants.pm impact ranks): a variant's
    // Sequence-Ontology consequence and IMPACT depend on WHERE it falls (coding CDS vs UTR vs
    // up/downstream) and, within coding, on the codon translation outcome (missense / synonymous /
    // stop_gained). Coding refinement (codon change) is computed only for coding SNV/MNV.
    //
    // The combinatorial point: region and coding outcome interact — inside the CDS the codon-level
    // consequence axis (missense/synonymous/stop_gained) drives the result with its Constants.pm
    // impact and a non-null codon change; OUTSIDE the CDS the region consequence dominates (Modifier,
    // null codon change) regardless of the coding-outcome axis. The single-exon transcript and
    // reference window are the canonical fixture from VariantAnnotator_FunctionalImpact_Tests.
    // ═══════════════════════════════════════════════════════════════════════

    private const string RefWindow = "ATGCAAGAATTATAANAAGGG"; // genomic 100..120
    private const int RefStart = 100;

    private static VariantAnnotator.Transcript CodingTranscript() => new(
        "ENST_TEST", "ENSG_TEST", "GENE_TEST", "chr1", 90, 130, '+',
        new List<(int, int)> { (90, 130) }, new List<(int, int)> { (100, 120) }, 100, 120);

    private static VariantAnnotator.FunctionalImpact Predict(int pos, string r, string a) =>
        VariantAnnotator.PredictFunctionalImpact(
            new VariantAnnotator.Variant("chr1", pos, r, a, VariantAnnotator.VariantType.SNV),
            CodingTranscript(), RefWindow, RefStart);

    public enum CodingOutcome { Missense, Synonymous, StopGained }
    public enum Region { Coding, FivePrimeUtr, Upstream }

    [Test, Combinatorial]
    public void VariantAnnot_ConsequenceAndImpact_AcrossOutcomeAndRegion(
        [Values(CodingOutcome.Missense, CodingOutcome.Synonymous, CodingOutcome.StopGained)] CodingOutcome outcome,
        [Values(Region.Coding, Region.FivePrimeUtr, Region.Upstream)] Region region)
    {
        // Coding scenario per outcome (positions/alleles from the canonical CDS layout).
        (int Pos, string Ref, string Alt, VariantAnnotator.ConsequenceType Cons, VariantAnnotator.ImpactLevel Impact, string Aa) sc = outcome switch
        {
            CodingOutcome.Missense => (107, "A", "T", VariantAnnotator.ConsequenceType.MissenseVariant, VariantAnnotator.ImpactLevel.Moderate, "p.E3V"),
            CodingOutcome.Synonymous => (111, "A", "G", VariantAnnotator.ConsequenceType.SynonymousVariant, VariantAnnotator.ImpactLevel.Low, "p.L4="),
            _ => (103, "C", "T", VariantAnnotator.ConsequenceType.StopGained, VariantAnnotator.ImpactLevel.High, "p.Q2*"),
        };

        if (region == Region.Coding)
        {
            var fi = Predict(sc.Pos, sc.Ref, sc.Alt);
            fi.Consequence.Should().Be(sc.Cons, "the codon outcome determines the coding consequence");
            fi.Impact.Should().Be(sc.Impact, "Constants.pm impact rank for the consequence");
            fi.AminoAcidChange.Should().Be(sc.Aa);
            fi.CodonChange.Should().NotBeNull("a coding SNV carries an HGVS codon change");
        }
        else
        {
            // Outside the CDS the region consequence dominates and the coding-outcome axis is inert.
            int pos = region == Region.FivePrimeUtr ? 95 : 80; // 95 = 5'UTR exon base; 80 = upstream (<TSS)
            var expected = region == Region.FivePrimeUtr
                ? VariantAnnotator.ConsequenceType.FivePrimeUtrVariant
                : VariantAnnotator.ConsequenceType.UpstreamGeneVariant;

            var fi = Predict(pos, "A", "T");
            fi.Consequence.Should().Be(expected, "region determines the consequence outside the CDS");
            fi.Impact.Should().Be(VariantAnnotator.ImpactLevel.Modifier, "UTR/upstream variants are MODIFIER impact");
            fi.CodonChange.Should().BeNull("no codon change outside the coding sequence");
        }
    }

    /// <summary>
    /// Interaction witness — impact rank tracks the consequence (HIGH stop_gained &gt; MODERATE missense
    /// &gt; LOW synonymous &gt; MODIFIER UTR), and GetImpactLevel agrees with PredictFunctionalImpact.
    /// </summary>
    [Test]
    public void VariantAnnot_ImpactRankTracksConsequence()
    {
        VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.StopGained).Should().Be(VariantAnnotator.ImpactLevel.High);
        VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.MissenseVariant).Should().Be(VariantAnnotator.ImpactLevel.Moderate);
        VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.SynonymousVariant).Should().Be(VariantAnnotator.ImpactLevel.Low);
        VariantAnnotator.GetImpactLevel(VariantAnnotator.ConsequenceType.FivePrimeUtrVariant).Should().Be(VariantAnnotator.ImpactLevel.Modifier);

        Predict(103, "C", "T").Impact.Should().Be(VariantAnnotator.ImpactLevel.High, "stop_gained is the most severe planted outcome");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: VARIANT-CALL-001 — Alignment-based variant calling (Variants)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 187.
    // Spec: tests/TestSpecs/VARIANT-CALL-001.md (canonical CallVariantsFromAlignment / ClassifyMutation).
    // ADVANCED §10.
    // Dimensions: minVaf(3) × minDepth(3) × minQual(3). Grid 3×3×3 = 27 (full, exhaustive).
    //
    // Model: CallVariantsFromAlignment scans a gapped pairwise alignment column-by-column, emitting a
    // SNP (both bases, differing), an Insertion (ref gap) or a Deletion (query gap) with correct
    // reference/query position bookkeeping; ClassifyMutation labels SNPs transition/transversion.
    //
    // Axis mapping (documented — this reference-vs-query caller has no read-pileup QC knobs): the
    // checklist's minVaf/minDepth/minQual map to the planted SNP / insertion / deletion COUNTS, which
    // vary independently. The combinatorial point: the caller returns exactly nSnp SNPs, nIns
    // insertions and nDel deletions with the correct allele shapes, at every (nSnp,nIns,nDel) cell.
    // ═══════════════════════════════════════════════════════════════════════

    [Test, Combinatorial]
    public void VariantCall_RecoversPlantedColumns_AcrossSnpInsDelCounts(
        [Values(0, 1, 2)] int nSnp,
        [Values(0, 1, 2)] int nIns,
        [Values(0, 1, 2)] int nDel)
    {
        // Build an aligned pair: matching anchor, nSnp A>G columns, nIns (-/T) columns, nDel (C/-) columns, anchor.
        string ar = "ACGT" + new string('A', nSnp) + new string('-', nIns) + new string('C', nDel) + "TGCA";
        string aq = "ACGT" + new string('G', nSnp) + new string('T', nIns) + new string('-', nDel) + "TGCA";

        var variants = VariantCaller.CallVariantsFromAlignment(ar, aq).ToList();

        variants.Count(v => v.Type == VariantType.SNP).Should().Be(nSnp, "one SNP per substituted column");
        variants.Count(v => v.Type == VariantType.Insertion).Should().Be(nIns, "one insertion per ref-gap column");
        variants.Count(v => v.Type == VariantType.Deletion).Should().Be(nDel, "one deletion per query-gap column");

        foreach (var v in variants)
        {
            switch (v.Type)
            {
                case VariantType.SNP:
                    v.ReferenceAllele.Should().NotBe(v.AlternateAllele);
                    v.ReferenceAllele.Should().Be("A"); v.AlternateAllele.Should().Be("G");
                    VariantCaller.ClassifyMutation(v).Should().Be(MutationType.Transition, "A↔G is a transition");
                    break;
                case VariantType.Insertion:
                    v.ReferenceAllele.Should().Be("-", "an insertion has a gap reference allele");
                    v.AlternateAllele.Should().Be("T");
                    break;
                case VariantType.Deletion:
                    v.ReferenceAllele.Should().Be("C");
                    v.AlternateAllele.Should().Be("-", "a deletion has a gap alternate allele");
                    break;
            }
        }
    }

    /// <summary>
    /// Interaction witness — transition/transversion classification and the Ti/Tv aggregate: A>G and
    /// C>T are transitions, A>C is a transversion, and CalculateTiTvRatio counts them.
    /// </summary>
    [Test]
    public void VariantCall_TransitionTransversionAndTiTv()
    {
        var ti1 = VariantCaller.CallVariantsFromAlignment("A", "G").Single();
        var ti2 = VariantCaller.CallVariantsFromAlignment("C", "T").Single();
        var tv = VariantCaller.CallVariantsFromAlignment("A", "C").Single();

        VariantCaller.ClassifyMutation(ti1).Should().Be(MutationType.Transition);
        VariantCaller.ClassifyMutation(ti2).Should().Be(MutationType.Transition);
        VariantCaller.ClassifyMutation(tv).Should().Be(MutationType.Transversion);

        VariantCaller.CalculateTiTvRatio(new[] { ti1, ti2, tv }).Should().BeApproximately(2.0, 1e-9, "2 transitions / 1 transversion");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: VARIANT-INDEL-001 — Insertion / deletion detection (Variants)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 188.
    // Spec: tests/TestSpecs/VARIANT-INDEL-001.md (canonical FindInsertions / FindDeletions). ADVANCED §10.
    // Dimensions: indelLen(3) × depth(3). Grid 3×3 = 9 (full, exhaustive ⊇ pairwise).
    //
    // Model: FindInsertions / FindDeletions align reference and query (global alignment) and filter to
    // Insertion / Deletion columns; an inserted/deleted run of length L yields L single-base indel
    // variants (insertion ⇒ ref-gap, query base; deletion ⇒ ref base, query-gap).
    //
    // Axis mapping (documented): indelLen → the per-site indel length; depth → the NUMBER of indel
    // sites. Distinct homopolymer anchors flank each indel so the global alignment is unambiguous. The
    // combinatorial point: the total reported indel bases equal sites×length, with the correct
    // ref/alt gap semantics, for both insertions and deletions.
    // ═══════════════════════════════════════════════════════════════════════

    private static readonly char[] AnchorBases = { 'A', 'C', 'G', 'T' };
    private static readonly char[] InsertBases = { 'G', 'T', 'A' }; // insert between anchor[i],anchor[i+1]; distinct from both

    private static (string Reference, string Query) BuildInsertionPair(int indelLen, int sites)
    {
        var refSb = new System.Text.StringBuilder(new string(AnchorBases[0], 6));
        var qrySb = new System.Text.StringBuilder(new string(AnchorBases[0], 6));
        for (int i = 0; i < sites; i++)
        {
            qrySb.Append(InsertBases[i], indelLen); // inserted run present only in the query
            string anchor = new(AnchorBases[i + 1], 6);
            refSb.Append(anchor);
            qrySb.Append(anchor);
        }
        return (refSb.ToString(), qrySb.ToString());
    }

    [Test, Combinatorial]
    public void VariantIndel_InsertionAndDeletionCounts_AcrossLengthAndSites(
        [Values(1, 2, 3)] int indelLen,
        [Values(1, 2, 3)] int sites)
    {
        var (refStr, qryStr) = BuildInsertionPair(indelLen, sites);
        var reference = new DnaSequence(refStr);
        var query = new DnaSequence(qryStr);

        int expectedBases = indelLen * sites;

        // Query has the extra runs ⇒ insertions; no deletions/SNPs in the optimal anchored alignment.
        var insertions = VariantCaller.FindInsertions(reference, query).ToList();
        insertions.Should().HaveCount(expectedBases, "each inserted base is one Insertion variant");
        insertions.Should().OnlyContain(v => v.Type == VariantType.Insertion && v.ReferenceAllele == "-" && v.AlternateAllele.Length == 1);
        VariantCaller.FindDeletions(reference, query).Should().BeEmpty("the query only adds bases");

        // Swapping the roles turns the same runs into deletions.
        var deletions = VariantCaller.FindDeletions(query, reference).ToList();
        deletions.Should().HaveCount(expectedBases, "each deleted base is one Deletion variant");
        deletions.Should().OnlyContain(v => v.Type == VariantType.Deletion && v.AlternateAllele == "-" && v.ReferenceAllele.Length == 1);

        // FindIndels is the union of the two.
        VariantCaller.FindIndels(reference, query).Should().HaveCount(expectedBases);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Unit: VARIANT-SNP-001 — Positional SNP detection (Variants)
    // Checklist: docs/checklists/09_COMBINATORIAL_TESTING.md, row 189.
    // Spec: tests/TestSpecs/VARIANT-SNP-001.md (canonical FindSnpsDirect). ADVANCED §10.
    // Dimensions: minVaf(3) × depth(3) × baseQual(3). Grid 3×3×3 = 27 (full, exhaustive).
    //
    // Model (Compeau & Pevzner Hamming): FindSnpsDirect compares equal positions of reference and
    // query and reports every mismatch as a SNP at that index; ClassifyMutation labels each as a
    // transition (purine↔purine / pyrimidine↔pyrimidine) or transversion.
    //
    // Axis mapping (documented — no read-pileup QC in this positional caller): minVaf → number of
    // planted SNPs; depth → sequence length; baseQual → the substitution class (all-transition /
    // all-transversion / mixed). The combinatorial point: FindSnpsDirect reports exactly the planted
    // mismatch positions, and each SNP's transition/transversion class matches the planted class.
    // ═══════════════════════════════════════════════════════════════════════

    public enum SubClass { Transition, Transversion, Mixed }

    private static string DiverseDna(int n, uint seed)
    {
        const string b = "ACGT";
        var c = new char[n];
        uint s = seed;
        for (int i = 0; i < n; i++) { s = s * 1664525u + 1013904223u; c[i] = b[(int)((s >> 16) & 3u)]; }
        return new string(c);
    }

    private static char Transition(char x) => x switch { 'A' => 'G', 'G' => 'A', 'C' => 'T', 'T' => 'C', _ => 'N' };
    private static char Transversion(char x) => x switch { 'A' => 'C', 'C' => 'A', 'G' => 'T', 'T' => 'G', _ => 'N' };

    [Test, Combinatorial]
    public void VariantSnp_RecoversPlantedSubstitutions_AcrossCountLengthClass(
        [Values(1, 2, 3)] int nSnp,
        [Values(12, 24, 48)] int seqLen,
        [Values(SubClass.Transition, SubClass.Transversion, SubClass.Mixed)] SubClass subClass)
    {
        string reference = DiverseDna(seqLen, 0x5EED);
        var q = reference.ToCharArray();
        var planted = new List<int>();
        for (int j = 0; j < nSnp; j++)
        {
            int pos = 1 + j * 3; // distinct, spaced, < seqLen
            bool useTransition = subClass switch
            {
                SubClass.Transition => true,
                SubClass.Transversion => false,
                _ => j % 2 == 0,
            };
            q[pos] = useTransition ? Transition(reference[pos]) : Transversion(reference[pos]);
            planted.Add(pos);
        }
        string query = new(q);

        var snps = VariantCaller.FindSnpsDirect(reference, query).ToList();

        snps.Select(v => v.Position).Should().BeEquivalentTo(planted, "exactly the planted mismatch positions are SNPs");
        snps.Should().OnlyContain(v => v.Type == VariantType.SNP);
        foreach (var v in snps)
        {
            v.ReferenceAllele.Should().Be(reference[v.Position].ToString());
            v.AlternateAllele.Should().Be(query[v.Position].ToString());
            var expectedClass = subClass == SubClass.Mixed
                ? (planted.IndexOf(v.Position) % 2 == 0 ? MutationType.Transition : MutationType.Transversion)
                : (subClass == SubClass.Transition ? MutationType.Transition : MutationType.Transversion);
            VariantCaller.ClassifyMutation(v).Should().Be(expectedClass, "the SNP's substitution class matches the planted class");
        }
    }

    /// <summary>
    /// Interaction witness — an identical query has no SNPs, and the Ti/Tv ratio reflects the planted
    /// substitution mix (all-transition ⇒ no transversions ⇒ ratio 0 by the method's convention).
    /// </summary>
    [Test]
    public void VariantSnp_NoSnpsWhenIdentical_AndTiTvReflectsMix()
    {
        string reference = DiverseDna(30, 0xABCD);
        VariantCaller.FindSnpsDirect(reference, reference).Should().BeEmpty("identical sequences have no SNPs");

        var refChars = reference.ToCharArray();
        var q = reference.ToCharArray();
        q[5] = Transition(reference[5]);
        q[10] = Transition(reference[10]);
        q[15] = Transversion(reference[15]);
        var snps = VariantCaller.FindSnpsDirect(reference, new string(q)).ToList();
        VariantCaller.CalculateTiTvRatio(snps).Should().BeApproximately(2.0, 1e-9, "2 transitions / 1 transversion");
    }
}
