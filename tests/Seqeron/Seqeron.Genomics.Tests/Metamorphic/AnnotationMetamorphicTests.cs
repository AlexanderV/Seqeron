using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using FluentAssertions;
using Seqeron.Genomics;
using Seqeron.Genomics.Annotation;
using Seqeron.Genomics.Core;

namespace Seqeron.Genomics.Tests.Metamorphic;

/// <summary>
/// Metamorphic tests for the Annotation area.
///
/// Each test encodes a metamorphic relation (MR) — a property relating the outputs of
/// multiple runs under an input transformation, with no hardcoded oracle. The relations
/// are derived from the ALGORITHM DEFINITION (and its documentation), not from the
/// implementation's observed output.
///
/// ───────────────────────────────────────────────────────────────────────────
/// Unit: ANNOT-ORF-001 — open reading frame detection (Annotation).
/// Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 28.
///
/// API under test (GenomeAnnotator.FindOrfs):
///   FindOrfs(dna, minLength = 100, searchBothStrands = true, requireStartCodon = true)
///     • Forward strand: in each of the 3 frames, scan codons; on a start codon
///       (ATG/GTG/TTG) push a pending start; on a stop codon (TAA/TAG/TGA) emit every
///       pending start s as the ORF [s, t+3) with aaLength = (t − s)/3, kept iff
///       aaLength ≥ minLength, then clear the pending starts.
///     • Reverse strand (searchBothStrands): the SAME scan runs on the reverse complement,
///       then each ORF is remapped to forward coordinates via
///       Start = |dna| − End_rev,  End = |dna| − Start_rev.
///   An ORF's identity is (Start, End, Frame, IsReverseComplement, Sequence, ProteinSequence);
///   Start/End are original-sequence coordinates (0-based, End exclusive).
///
/// Relations (derived from the definition / ORF_Detection.md §2.2–§2.4, NOT from output):
///   • MON   — minLength is ONLY a filter (aaLength ≥ minLength). Lowering it keeps every
///             ORF that already passed and may admit more, so the ORF set is a SUPERSET and
///             the count is non-decreasing along a decreasing-minLength chain.
///   • SHIFT — Start/End are forward-strand coordinates for both strands. Prepending an
///             IN-FRAME, non-coding flank F with |F| ≡ 0 (mod 3) preserves every reading
///             frame and creates no start/stop, so every ORF's Start and End advance by
///             exactly |F| with Frame, strand, Sequence and ProteinSequence preserved.
///             (Poly-C is start/stop-free; its reverse complement poly-G is too, so neither
///             strand gains an ORF; |F| ≡ 0 mod 3 keeps the frame labels.)
///   • INV   — a forward-strand ORF lying ENTIRELY upstream of an insertion point depends
///             only on codons before that point, so inserting ANY region at/after its End
///             leaves it byte-for-byte unchanged (same coordinates, sequence, protein).
/// ───────────────────────────────────────────────────────────────────────────
/// </summary>
[TestFixture]
[Category("Metamorphic")]
public class AnnotationMetamorphicTests
{
    #region Helpers

    /// <summary>Base seed so random inputs are reproducible.</summary>
    private const int RngSeed = 20260619;

    /// <summary>
    /// Generates a random DNA string of the given length over {A,C,G,T}.
    /// A FRESH, locally-seeded <see cref="Random"/> is used per call (seed = base ⊕ length)
    /// so the result is a pure, reproducible function of the length and is thread-safe under
    /// NUnit's parallel fixture execution — a single shared <see cref="Random"/> is not
    /// thread-safe and would yield non-deterministic, degenerate sequences when called
    /// concurrently from parallel test methods.
    /// </summary>
    private static string RandomDna(int length)
    {
        var rng = new Random(RngSeed + length);
        const string bases = "ACGT";
        var chars = new char[length];
        for (int i = 0; i < length; i++)
            chars[i] = bases[rng.Next(bases.Length)];
        return new string(chars);
    }

    /// <summary>Run-order-independent identity of an ORF.</summary>
    private static (int Start, int End, int Frame, bool Rc, string Seq, string Prot) OrfId(GenomeAnnotator.OpenReadingFrame o)
        => (o.Start, o.End, o.Frame, o.IsReverseComplement, o.Sequence, o.ProteinSequence);

    private static HashSet<(int, int, int, bool, string, string)> OrfSet(
        string dna, int minLength, bool bothStrands = true)
        => GenomeAnnotator.FindOrfs(dna, minLength, bothStrands).Select(OrfId).ToHashSet();

    /// <summary>
    /// Bodies carrying explicit forward-frame-0 ORFs (start codon … in-frame codons … stop),
    /// plus a fixed-seed random body. Codons avoid in-frame stops so the ORFs are unambiguous.
    /// </summary>
    private static IEnumerable<string> OrfBodies()
    {
        // ATG AAA CCC AAA TAA (aaLen 4) · spacer · ATG CCC AAA CCC TGA (aaLen 4)
        yield return "ATGAAACCCAAATAA" + "GGG" + "ATGCCCAAACCCTGA";
        // A single longer ORF: ATG + 8×AAA + TAA (aaLen 9)
        yield return "ATG" + string.Concat(Enumerable.Repeat("AAA", 8)) + "TAA";
        // Two ORFs sharing structure with extra GTG/TTG starts (alternate start codons)
        yield return "GTGAAATTTAAACCCTAA" + "TTT" + "TTGCCCAAATAA";
        // Fixed-seed random body (relations must hold for arbitrary input too)
        yield return RandomDna(120);
    }

    /// <summary>A small minimum length (in aa) so the short test ORFs qualify.</summary>
    private const int MinAa = 3;

    #endregion

    #region MON — lowering minLength yields a superset of ORFs (count non-decreasing)

    [Test]
    [Description("MON: along a DECREASING minLength chain the ORF set grows monotonically — each higher-threshold set is a subset of every lower-threshold one, count non-decreasing.")]
    public void FindOrfs_LoweringMinLength_YieldsSuperset_CountNonDecreasing()
    {
        int[] decreasingMinLength = { 10, 6, 4, 3, 2, 1 };

        foreach (var body in OrfBodies())
        {
            HashSet<(int, int, int, bool, string, string)>? previous = null;
            int previousCount = -1;

            foreach (int minLen in decreasingMinLength)
            {
                var orfs = OrfSet(body, minLen);

                if (previous is not null)
                {
                    orfs.IsSupersetOf(previous).Should().BeTrue(
                        because: $"minLength is only the filter aaLength ≥ minLength, so lowering it to {minLen} keeps every ORF that already passed and may add more");
                    orfs.Count.Should().BeGreaterThanOrEqualTo(previousCount,
                        because: $"a lower minLength ({minLen}) admits-or-keeps each ORF, never removes one — the count is non-decreasing");
                }

                previous = orfs;
                previousCount = orfs.Count;
            }
        }
    }

    [Test]
    [Description("MON: every ORF found under a stricter minLength is also found under a looser one (membership preserved when the threshold drops).")]
    public void FindOrfs_StrictThreshold_SubsetOfLooseThreshold()
    {
        foreach (var body in OrfBodies())
        {
            var strict = OrfSet(body, minLength: 8);
            var loose = OrfSet(body, minLength: 2);

            loose.IsSupersetOf(strict).Should().BeTrue(
                because: "an ORF whose aaLength ≥ 8 also satisfies aaLength ≥ 2, and no other rule changed, so the strict ORF set is a subset of the loose one");
        }
    }

    #endregion

    #region SHIFT — prepending an in-frame non-coding flank shifts every ORF by |F|

    [Test]
    [Description("SHIFT: prepending a poly-C flank of length ≡ 0 (mod 3) advances every ORF's Start and End by exactly |F|, preserving Frame, strand, Sequence and ProteinSequence.")]
    public void FindOrfs_PrependInFrameNonCodingFlank_ShiftsAllOrfsByFlankLength()
    {
        foreach (var body in OrfBodies())
        {
            var baseOrfs = GenomeAnnotator.FindOrfs(body, MinAa, searchBothStrands: true).ToList();
            baseOrfs.Should().NotBeEmpty(because: "each body embeds at least one ORF at the small minimum length");

            foreach (int flankLen in new[] { 3, 6, 12 })
            {
                string flank = new string('C', flankLen); // poly-C: no start/stop; revComp poly-G: no start/stop
                var shifted = OrfSet(flank + body, MinAa);

                var expected = baseOrfs
                    .Select(o => (o.Start + flankLen, o.End + flankLen, o.Frame, o.IsReverseComplement, o.Sequence, o.ProteinSequence))
                    .ToHashSet();

                shifted.SetEquals(expected).Should().BeTrue(
                    because: $"an in-frame (|F| = {flankLen} ≡ 0 mod 3) non-coding prefix advances every ORF by exactly {flankLen} on both strands " +
                             "while preserving frame, strand, sequence and protein — and creates no new ORF (poly-C / poly-G carry no start or stop codon)");
            }
        }
    }

    [Test]
    [Description("SHIFT anchor: a known forward ORF starting at index 0 moves to exactly |F| when an in-frame poly-C flank is prepended.")]
    public void FindOrfs_PrependFlank_KnownForwardOrf_ShiftsExactly()
    {
        const string body = "ATG" + "AAACCCAAACCC" + "TAA"; // ATG at 0; stop TAA at 15; End 18
        var baseOrf = GenomeAnnotator.FindOrfs(body, MinAa, searchBothStrands: false)
            .Single(o => o.Frame == 1 && o.Start == 0);
        baseOrf.End.Should().Be(18, because: "the ORF spans ATG…TAA over 18 nucleotides");

        const int flankLen = 9;
        var shiftedOrf = GenomeAnnotator.FindOrfs(new string('C', flankLen) + body, MinAa, searchBothStrands: false)
            .Single(o => o.Frame == 1 && o.IsReverseComplement == false && o.Sequence == baseOrf.Sequence);

        shiftedOrf.Start.Should().Be(0 + flankLen, because: "the in-frame poly-C prefix shifts the ORF start by exactly the flank length");
        shiftedOrf.End.Should().Be(18 + flankLen, because: "the ORF end shifts by the same flank length");
    }

    #endregion

    #region INV — inserting a region downstream does not change ORFs entirely upstream of it

    [Test]
    [Description("INV: inserting ANY region at a codon boundary leaves every forward-strand ORF that lies entirely upstream of the insertion point byte-for-byte unchanged.")]
    public void FindOrfs_InsertDownstream_UpstreamForwardOrfsUnchanged()
    {
        // prefix ends right after a complete forward ORF (ATG…TAA); insertion goes at |prefix|.
        const string prefix = "ATGAAACCCAAATAA";           // forward ORF [0,15)
        const string suffix = "GGGATGCCCAAACCCTGA";          // a downstream ORF (will move)
        int insertPos = prefix.Length;

        // Forward-only search isolates the "upstream" direction on the coding strand.
        var baseUpstream = GenomeAnnotator.FindOrfs(prefix + suffix, MinAa, searchBothStrands: false)
            .Where(o => o.End <= insertPos).Select(OrfId).ToHashSet();
        baseUpstream.Should().NotBeEmpty(because: "the prefix contains a complete forward ORF upstream of the insertion point");

        foreach (var insert in new[] { "CCC", "CCCCCCCCC", RandomDna(9), RandomDna(12), "ATGTTTGGGTAA" })
        {
            var changedUpstream = GenomeAnnotator.FindOrfs(prefix + insert + suffix, MinAa, searchBothStrands: false)
                .Where(o => o.End <= insertPos).Select(OrfId).ToHashSet();

            changedUpstream.SetEquals(baseUpstream).Should().BeTrue(
                because: "an ORF ending at or before the insertion point reads only codons before it, which the downstream insertion never touches — " +
                         "so its coordinates, sequence and protein are preserved exactly regardless of what is inserted");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ANNOT-GENE-001 — gene prediction (Annotation).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 29.
    //
    // API under test (GenomeAnnotator.PredictGenes):
    //   PredictGenes(dna, minOrfLength = 100, prefix = "gene")
    //     = FindOrfs(dna, minOrfLength, searchBothStrands: true, requireStartCodon: true),
    //       ORDERED by Start, each ORF emitted as a CDS GeneAnnotation whose span is the
    //       ORF span (Start, End), Strand = '+' forward / '-' reverse, with the ORF's frame,
    //       protein length and translation copied into Attributes. GeneId is the 1-based rank
    //       in the Start-ordered list ("{prefix}_{n:D4}").
    //
    // Relations (derived from this ORF→gene mapping, NOT from output):
    //   • COMP (gene ⊃ ORF): PredictGenes is a 1:1 image of FindOrfs — each gene's span
    //           equals an ORF's span and carries that ORF's translation, so the gene set is
    //           in bijection with the ORF set (gene count = ORF count, equal span sets).
    //   • INV  (non-coding insertion doesn't affect upstream): genes are ORFs, so a forward
    //           gene entirely upstream of an insertion point keeps its span/strand/product
    //           regardless of what is inserted downstream (the order-dependent GeneId is
    //           excluded, since downstream/other-strand shifts can renumber later genes).
    //   • MON  (longer seq → ≥ genes): genes are ORFs; an IN-FRAME append (|X| ≡ 0 mod 3)
    //           preserves every forward ORF and keeps the reverse-complement frames aligned,
    //           so the gene count is non-decreasing (and strictly grows when the appended
    //           region contributes a new ORF).
    // ───────────────────────────────────────────────────────────────────────────

    private static (int Start, int End, char Strand) GeneSpan(GenomeAnnotator.GeneAnnotation g)
        => (g.Start, g.End, g.Strand);

    private static (int Start, int End, char Strand) OrfSpan(GenomeAnnotator.OpenReadingFrame o)
        => (o.Start, o.End, o.IsReverseComplement ? '-' : '+');

    #region COMP — every predicted gene is exactly an ORF (1:1 image of FindOrfs)

    [Test]
    [Description("COMP: PredictGenes is a 1:1 image of FindOrfs — gene count equals ORF count, gene spans equal ORF spans, and each gene carries its ORF's translation.")]
    public void PredictGenes_IsBijectiveImageOfFindOrfs()
    {
        foreach (var body in OrfBodies())
        {
            var genes = GenomeAnnotator.PredictGenes(body, MinAa).ToList();
            var orfs = GenomeAnnotator.FindOrfs(body, MinAa, searchBothStrands: true, requireStartCodon: true).ToList();

            genes.Count.Should().Be(orfs.Count, because: "PredictGenes emits exactly one gene per ORF");
            genes.Select(GeneSpan).ToHashSet().SetEquals(orfs.Select(OrfSpan).ToHashSet()).Should().BeTrue(
                because: "each gene's [Start,End) and strand are copied verbatim from an ORF, so the span sets coincide (gene region ⊇ the ORF)");

            genes.Should().OnlyContain(g => g.End > g.Start && g.Type == "CDS",
                because: "a gene built from an ORF spans the coding interval and is annotated as a CDS");

            // Each gene carries its ORF's translation — the gene genuinely contains the ORF.
            foreach (var gene in genes)
            {
                var orf = orfs.Single(o => OrfSpan(o) == GeneSpan(gene));
                gene.Attributes["translation"].Should().Be(orf.ProteinSequence,
                    because: "the gene's translation attribute is the ORF's protein sequence");
            }
        }
    }

    #endregion

    #region INV — inserting downstream does not change genes entirely upstream of it

    [Test]
    [Description("INV: inserting any region at a codon boundary leaves every forward gene that lies entirely upstream of the insertion point unchanged in span, strand, type and product.")]
    public void PredictGenes_InsertDownstream_UpstreamForwardGenesUnchanged()
    {
        const string prefix = "ATGAAACCCAAATAA";           // forward CDS [0,15)
        const string suffix = "GGGATGCCCAAACCCTGA";          // a downstream CDS (will move)
        int insertPos = prefix.Length;

        (int, int, char, string, string) GeneKey(GenomeAnnotator.GeneAnnotation g)
            => (g.Start, g.End, g.Strand, g.Type, g.Product);

        var baseUpstream = GenomeAnnotator.PredictGenes(prefix + suffix, MinAa)
            .Where(g => g.Strand == '+' && g.End <= insertPos).Select(GeneKey).ToHashSet();
        baseUpstream.Should().NotBeEmpty(because: "the prefix contains a complete forward CDS upstream of the insertion point");

        foreach (var insert in new[] { "CCC", "CCCCCCCCC", RandomDna(9), "ATGTTTGGGTAA" })
        {
            var changedUpstream = GenomeAnnotator.PredictGenes(prefix + insert + suffix, MinAa)
                .Where(g => g.Strand == '+' && g.End <= insertPos).Select(GeneKey).ToHashSet();

            changedUpstream.SetEquals(baseUpstream).Should().BeTrue(
                because: "a forward gene ending at or before the insertion point reads only upstream codons, so a downstream insertion leaves its span/strand/product unchanged (GeneId aside)");
        }
    }

    #endregion

    #region MON — extending the sequence in frame never reduces the gene count

    [Test]
    [Description("MON: appending an in-frame region (|X| ≡ 0 mod 3) never reduces the predicted-gene count; appending a region that carries an ORF strictly increases it.")]
    public void PredictGenes_InFrameAppend_GeneCountNonDecreasing()
    {
        bool sawStrictIncrease = false;

        foreach (var body in OrfBodies())
        {
            int baseCount = GenomeAnnotator.PredictGenes(body, MinAa).Count();

            // All appends have length ≡ 0 (mod 3) so the reverse-complement frames stay aligned.
            foreach (var ext in new[] { "CCCCCC", "CCCCCCCCCCCC", RandomDna(9), RandomDna(12) })
            {
                int extCount = GenomeAnnotator.PredictGenes(body + ext, MinAa).Count();
                extCount.Should().BeGreaterThanOrEqualTo(baseCount,
                    because: $"an in-frame append preserves every forward ORF and keeps the reverse-strand frames aligned, so the gene count cannot drop (|ext| = {ext.Length})");
            }

            // Appending an explicit in-frame ORF adds at least that gene.
            int withOrf = GenomeAnnotator.PredictGenes(body + "ATGAAACCCTAA", MinAa).Count();
            withOrf.Should().BeGreaterThanOrEqualTo(baseCount);
            if (withOrf > baseCount) sawStrictIncrease = true;
        }

        sawStrictIncrease.Should().BeTrue(
            because: "appending a fresh ATG…TAA ORF introduces a new gene for at least one body — the monotone relation is exercised, not vacuous");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ANNOT-PROM-001 — promoter motif finding (Annotation).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 30.
    //
    // API under test (GenomeAnnotator.FindPromoterMotifs):
    //   FindPromoterMotifs(dna) enumerates every FORWARD-strand exact match of the bacterial
    //   −35 box variants {TTGACA, TTGAC, TGACA, TTGA} and −10 box variants
    //   {TATAAT, TATAA, ATAAT, TATA}, emitting (position, type, motif, score) with a FIXED
    //   probability-weighted score per variant (Harley & Reynolds 1987). There is no
    //   threshold parameter — the score cutoff is applied by the caller.
    //
    // Relations (derived from the exact-match definition, NOT from output):
    //   • MON   — a "promoter" is a hit with score ≥ threshold. Each hit's score is fixed by
    //             its variant, so LOWERING the threshold keeps every passing hit and may add
    //             more — the passing set is a SUPERSET, count non-decreasing.
    //   • SHIFT — positions are exact-substring indices on the forward strand. Prepending a
    //             poly-C flank (no motif starts with C, so neither an internal nor a junction
    //             match can form) advances every hit's position by exactly |F|, preserving
    //             type, motif and score.
    //   • INV   — a motif hit is a LOCAL exact substring, so changing the sequence strictly
    //             downstream of a hit (after its last base) leaves that hit unchanged: every
    //             motif ending at or before the change point is preserved exactly.
    // ───────────────────────────────────────────────────────────────────────────

    private static (int Position, string Type, string Motif, double Score) PromoterId(
        (int position, string type, string sequence, double score) h)
        => (h.position, h.type, h.sequence, h.score);

    /// <summary>Bodies embedding canonical −35 (TTGACA) and −10 (TATAAT) boxes, separated by neutral poly-C.</summary>
    private static IEnumerable<string> PromoterBodies()
    {
        yield return "TTGACA" + "CCCCC" + "TATAAT" + "CCCCC";
        yield return "ACGT" + "TTGACA" + "ACGTACGT" + "TATAAT" + "ACGT" + "TTGA" + "ACGT";
    }

    #region MON — lowering the score threshold yields a superset of promoter hits

    [Test]
    [Description("MON: along a DECREASING score-threshold chain the set of passing promoter hits grows monotonically — each higher cutoff's hits are a subset of every lower cutoff's.")]
    public void FindPromoterMotifs_LoweringScoreThreshold_YieldsSuperset_CountNonDecreasing()
    {
        double[] decreasingThresholds = { 1.0, 0.85, 0.81, 0.70, 0.60, 0.0 };

        foreach (var body in PromoterBodies().Append("ACTTGACATATAATGGGTATACCCTTGACGTAA").Append(RandomDna(120)))
        {
            var hits = GenomeAnnotator.FindPromoterMotifs(body).ToList();

            HashSet<(int, string, string, double)>? previous = null;
            int previousCount = -1;

            foreach (double t in decreasingThresholds)
            {
                var passing = hits.Where(h => h.score >= t).Select(PromoterId).ToHashSet();

                if (previous is not null)
                {
                    passing.IsSupersetOf(previous).Should().BeTrue(
                        because: $"each hit's score is fixed, so lowering the cutoff to {t} keeps every passing hit and may add more — the passing set is a superset");
                    passing.Count.Should().BeGreaterThanOrEqualTo(previousCount,
                        because: $"a lower score threshold ({t}) admits-or-keeps each hit, never removes one — count is non-decreasing");
                }

                previous = passing;
                previousCount = passing.Count;
            }
        }
    }

    #endregion

    #region SHIFT — prepending a neutral poly-C flank advances every promoter hit by |F|

    [Test]
    [Description("SHIFT: prepending a poly-C flank (no promoter motif starts with C) advances every promoter hit's position by exactly |F|, preserving type, motif and score.")]
    public void FindPromoterMotifs_PrependNeutralFlank_ShiftsAllHitsByFlankLength()
    {
        foreach (var body in PromoterBodies())
        {
            var baseHits = GenomeAnnotator.FindPromoterMotifs(body).ToList();
            baseHits.Should().NotBeEmpty(because: "each body embeds canonical −35 and −10 boxes");

            foreach (int flankLen in new[] { 1, 5, 17 })
            {
                string flank = new string('C', flankLen);
                GenomeAnnotator.FindPromoterMotifs(flank).Should().BeEmpty(
                    because: "a poly-C flank contains none of the −35/−10 box variants");

                var shifted = GenomeAnnotator.FindPromoterMotifs(flank + body).Select(PromoterId).ToHashSet();
                var expected = baseHits.Select(h => (h.position + flankLen, h.type, h.sequence, h.score)).ToHashSet();

                shifted.SetEquals(expected).Should().BeTrue(
                    because: $"positions are forward-strand exact-match indices, so a neutral length-{flankLen} prefix advances every hit by exactly that amount " +
                             "while preserving type, motif and score (no motif starts with C, so no junction match appears)");
            }
        }
    }

    #endregion

    #region INV — a change downstream of a promoter hit leaves it unchanged

    [Test]
    [Description("INV: changing the sequence downstream of a promoter motif leaves every hit ending at or before the change point unchanged — promoter detection is a local exact match.")]
    public void FindPromoterMotifs_DownstreamChange_UpstreamHitsUnchanged()
    {
        const string prefix = "TTGACA" + "CCC" + "TATAAT" + "CCC"; // −35 and −10 boxes, fully within the prefix
        int changePoint = prefix.Length;

        var baseUpstream = GenomeAnnotator.FindPromoterMotifs(prefix + "ACGTACGTACGT")
            .Where(h => h.position + h.sequence.Length <= changePoint).Select(PromoterId).ToHashSet();
        baseUpstream.Should().NotBeEmpty(because: "the prefix contains promoter boxes upstream of the change point");

        // Several distinct downstream regions, including motif-rich and fixed-seed random ones.
        foreach (var suffix in new[] { "GGGGGGGG", "TTGACATATAAT", RandomDna(20), RandomDna(40), "" })
        {
            var changedUpstream = GenomeAnnotator.FindPromoterMotifs(prefix + suffix)
                .Where(h => h.position + h.sequence.Length <= changePoint).Select(PromoterId).ToHashSet();

            changedUpstream.SetEquals(baseUpstream).Should().BeTrue(
                because: "a motif ending at or before the change point reads only upstream bases, so any downstream change leaves it exactly in place");
        }
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ANNOT-GFF-001 — GFF3 serialization round-trip (Annotation).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 31.
    //
    // API under test (GenomeAnnotator.ToGff3 / ParseGff3):
    //   ToGff3(annotations, seqId) emits "##gff-version 3" then one tab-delimited line per
    //     annotation: seqId . Type (Start+1) End . Strand phase attrs, where phase = 0 for
    //     CDS else '.', attrs = "ID={GeneId};product={Product};{other attrs except translation}"
    //     with GFF3 column-9 percent-encoding.
    //   ParseGff3(lines) skips blank/'#' lines and <9-column lines, then reads each feature:
    //     Start = col4 (kept 1-based), End = col5, Strand = col7[0], phase from col8,
    //     attributes split on ';' into a dict (order-independent, %-decoded), FeatureId = ID.
    //
    // Relations (derived from the format definition, NOT from output):
    //   • INV (line count): writing N annotations always yields N+1 lines (header + N), and
    //         parsing recovers exactly N features — the record count survives the round-trip.
    //   • COMP (parse∘write): each parsed feature reproduces its annotation's Type, End,
    //         Strand, ID and product, with Start in GFF 1-based form (annotation.Start + 1)
    //         and CDS phase 0; percent-encoded attribute values decode back to the originals.
    //   • INV (attribute order): column 9 parses into a dictionary, so permuting the
    //         attribute order yields an identical feature (same FeatureId and attribute map).
    // ───────────────────────────────────────────────────────────────────────────

    private static GenomeAnnotator.GeneAnnotation Gene(
        string id, int start, int end, char strand, string type, string product,
        params (string Key, string Value)[] attrs)
        => new(id, start, end, strand, type, product,
               attrs.ToDictionary(a => a.Key, a => a.Value));

    /// <summary>Annotations exercising CDS/non-CDS, both strands, extra attributes and special characters.</summary>
    private static List<GenomeAnnotator.GeneAnnotation> GffAnnotations() => new()
    {
        Gene("gene_0001", 10, 60, '+', "CDS", "hypothetical protein",
            ("frame", "1"), ("protein_length", "16"), ("translation", "MKLV")),
        Gene("gene_0002", 100, 130, '-', "gene", "regulatory; element=ABC",
            ("frame", "2"), ("note", "has,comma and=equals")),
        Gene("gene_0003", 200, 260, '+', "CDS", "enzyme alpha beta"),
    };

    #region INV — the record count survives write → parse

    [Test]
    [Description("INV: ToGff3 emits one header line plus one line per annotation, and ParseGff3 recovers exactly that many features — the record count is preserved.")]
    public void ToGff3_ThenParse_PreservesRecordCount()
    {
        var annotations = GffAnnotations();
        var lines = GenomeAnnotator.ToGff3(annotations).ToList();

        lines[0].Should().Be("##gff-version 3", because: "ToGff3 always starts with the GFF3 version pragma");
        lines.Count.Should().Be(annotations.Count + 1, because: "the output is the header followed by exactly one line per annotation");

        var features = GenomeAnnotator.ParseGff3(lines).ToList();
        features.Count.Should().Be(annotations.Count,
            because: "ParseGff3 skips the '#' header and yields one feature per data line, so the record count round-trips");
    }

    #endregion

    #region COMP — parse(write(x)) reproduces each annotation's fields (GFF 1-based coordinates)

    [Test]
    [Description("COMP: parsing the GFF3 written for each annotation reproduces its Type, End, Strand, ID and product, with Start in 1-based GFF form and CDS phase 0; encoded attribute values decode back to the originals.")]
    public void ParseGff3_OfWrittenAnnotations_ReproducesFields()
    {
        var annotations = GffAnnotations();
        var features = GenomeAnnotator.ParseGff3(GenomeAnnotator.ToGff3(annotations)).ToList();

        features.Count.Should().Be(annotations.Count);
        for (int i = 0; i < annotations.Count; i++)
        {
            var ann = annotations[i];
            var feat = features[i];

            feat.FeatureId.Should().Be(ann.GeneId, because: "the ID attribute carries the gene id");
            feat.Type.Should().Be(ann.Type);
            feat.Start.Should().Be(ann.Start + 1, because: "GFF coordinates are 1-based: the writer emits Start+1 and the parser keeps it 1-based");
            feat.End.Should().Be(ann.End);
            feat.Strand.Should().Be(ann.Strand);
            feat.Phase.Should().Be(ann.Type == "CDS" ? 0 : (int?)null, because: "phase is 0 for CDS features and '.' (null) otherwise");

            feat.Attributes["ID"].Should().Be(ann.GeneId);
            feat.Attributes["product"].Should().Be(ann.Product,
                because: "percent-encoded special characters in the product decode back to the original value");

            foreach (var (key, value) in ann.Attributes)
            {
                if (key == "translation")
                    feat.Attributes.ContainsKey(key).Should().BeFalse(because: "ToGff3 intentionally omits the large translation attribute");
                else
                    feat.Attributes[key].Should().Be(value, because: $"attribute '{key}' round-trips through encode/decode");
            }
        }
    }

    #endregion

    #region INV — attribute order in column 9 is irrelevant

    [Test]
    [Description("INV: two GFF3 lines that differ only in the order of their column-9 attributes parse to the same feature (same id, type, coordinates and attribute map).")]
    public void ParseGff3_AttributeOrder_DoesNotAffectFeature()
    {
        const string lineA = "seq1\t.\tgene\t10\t20\t.\t+\t.\tID=g1;product=foo;frame=1;color=red";
        const string lineB = "seq1\t.\tgene\t10\t20\t.\t+\t.\tcolor=red;frame=1;product=foo;ID=g1";

        var a = GenomeAnnotator.ParseGff3(new[] { lineA }).Single();
        var b = GenomeAnnotator.ParseGff3(new[] { lineB }).Single();

        b.FeatureId.Should().Be(a.FeatureId);
        b.Type.Should().Be(a.Type);
        b.Start.Should().Be(a.Start);
        b.End.Should().Be(a.End);
        b.Strand.Should().Be(a.Strand);

        b.Attributes.Should().BeEquivalentTo(a.Attributes,
            because: "column 9 parses into a key/value dictionary, so permuting the attribute order yields the same attribute map");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ANNOT-CODING-001 — coding-potential scoring (Annotation).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 216.
    //
    // API under test (GenomeAnnotator.CalculateCodingPotential):
    //   CPAT hexamer log-likelihood (Wang et al. 2013): the mean over in-frame hexamers of
    //   ln(coding[k]/noncoding[k]). Positive ⇒ coding-like, negative ⇒ non-coding-like.
    //
    // Relations (derived from the log-likelihood-ratio definition, NOT from output):
    //   • INV  (deterministic): a pure, case-insensitive function of (sequence, tables).
    //   • MON  (real ORF ⇒ higher score): a sequence built from coding-enriched hexamers (a CDS
    //          signature) scores strictly above a neutral one (coding = non-coding), which in turn
    //          scores above a background sequence built from non-coding-enriched hexamers.
    // ───────────────────────────────────────────────────────────────────────────

    #region ANNOT-CODING-001 — Helpers

    // In-frame (step-3) hexamer tables. "ACGACG" is coding-enriched (4:1), "TTTTTT" is
    // non-coding-enriched (1:4), and "GGGGGG" is neutral (2:2 ⇒ ln 1 = 0).
    private static readonly IReadOnlyDictionary<string, double> CodingHexamers = new Dictionary<string, double>
    {
        ["ACGACG"] = 4, ["GGGGGG"] = 2, ["TTTTTT"] = 1,
    };
    private static readonly IReadOnlyDictionary<string, double> NoncodingHexamers = new Dictionary<string, double>
    {
        ["ACGACG"] = 1, ["GGGGGG"] = 2, ["TTTTTT"] = 4,
    };

    // Period-3 sequences whose every step-3 hexamer window equals the repeated 6-mer.
    private const string CodingLikeSeq = "ACGACGACGACGACGACGACGACGACGACG"; // → "ACGACG" hexamers
    private const string NeutralSeq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";    // → "GGGGGG" hexamers
    private const string NoncodingLikeSeq = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"; // → "TTTTTT" hexamers

    private static double CodingPotential(string seq) =>
        GenomeAnnotator.CalculateCodingPotential(seq, CodingHexamers, NoncodingHexamers);

    #endregion

    #region ANNOT-CODING-001 INV — coding potential is deterministic and case-insensitive

    [Test]
    [Description("INV: CalculateCodingPotential is a pure, case-insensitive function of (sequence, tables).")]
    public void CodingPotential_SameInput_SameScore()
    {
        double first = CodingPotential(CodingLikeSeq);

        CodingPotential(CodingLikeSeq).Should().Be(first, because: "scoring has no hidden state");
        CodingPotential(CodingLikeSeq.ToLowerInvariant()).Should().Be(first,
            because: "the sequence is upper-cased before scoring, so case does not matter");
    }

    #endregion

    #region ANNOT-CODING-001 MON — coding-enriched sequence scores above background

    [Test]
    [Description("MON: a coding-enriched sequence scores strictly above a neutral one, which scores above a non-coding-enriched background.")]
    public void CodingPotential_MoreCodingLike_HigherScore()
    {
        double coding = CodingPotential(CodingLikeSeq);
        double neutral = CodingPotential(NeutralSeq);
        double noncoding = CodingPotential(NoncodingLikeSeq);

        coding.Should().BeGreaterThan(neutral, because: "coding-enriched hexamers have a positive log-likelihood ratio");
        neutral.Should().BeGreaterThan(noncoding, because: "non-coding-enriched hexamers have a negative log-likelihood ratio");

        coding.Should().BePositive(because: "a CDS-signature sequence is scored as coding");
        noncoding.Should().BeNegative(because: "a background sequence is scored as non-coding");
        neutral.Should().BeApproximately(0.0, 1e-12, because: "equal coding and non-coding frequencies give ln 1 = 0");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ANNOT-CODONUSAGE-001 — relative synonymous codon usage over a CDS set (Annotation).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 217.
    //
    // API under test (GenomeAnnotator.GetCodonUsage(IEnumerable<string>)):
    //   Pools codon counts across the coding sequences and returns RSCU per sense codon:
    //   RSCU_j = nᵢ·x_j/Σx over the synonymous family of amino acid i (Sharp & Li 1986).
    //
    // Relations (derived from the RSCU normalisation, NOT from output):
    //   • INV (codon order independent): RSCU depends only on the POOLED codon counts, so permuting
    //         the codons within a sequence — or reordering the input sequences — leaves it unchanged.
    //   • P   (per-AA sum = 1, frequency form): for each present amino acid ΣRSCU = nᵢ (its
    //         degeneracy), so the normalised frequency RSCU/nᵢ sums to exactly 1 over the family.
    // ───────────────────────────────────────────────────────────────────────────

    #region ANNOT-CODONUSAGE-001 — Helpers

    // A coding sequence exercising the F/I/V/L families (degeneracy 2/3/4/6), each with a 2:1
    // preference on its first codon so the RSCU values genuinely depart from 1.
    private const string AnnotCodingSeq =
        "TTTTTTTTC" +                            // F: TTT×2, TTC×1
        "ATTATTATCATA" +                         // I: ATT×2, ATC×1, ATA×1
        "GTTGTTGTCGTAGTG" +                      // V: GTT×2, GTC×1, GTA×1, GTG×1
        "CTTCTTCTCCTACTGTTATTG";                 // L: CTT×2, CTC×1, CTA×1, CTG×1, TTA×1, TTG×1

    private static List<string> SplitCodons(string seq)
    {
        var codons = new List<string>();
        for (int i = 0; i + 3 <= seq.Length; i += 3)
            codons.Add(seq.Substring(i, 3));
        return codons;
    }

    // Synonymous families (DNA codons) for each sense amino acid of the Standard code.
    private static IEnumerable<(char Aa, List<string> Codons)> SenseFamilies() =>
        GeneticCode.Standard.CodonTable
            .Where(kv => kv.Value != '*')
            .GroupBy(kv => kv.Value, kv => kv.Key.Replace('U', 'T'))
            .Select(g => (g.Key, g.ToList()));

    #endregion

    #region ANNOT-CODONUSAGE-001 INV — RSCU is independent of codon and sequence order

    [Test]
    [Description("INV: RSCU depends only on the pooled codon counts, so permuting codons within a sequence and reordering the input sequences leave it unchanged.")]
    public void CodonUsage_OrderInvariant()
    {
        var original = GenomeAnnotator.GetCodonUsage(new[] { AnnotCodingSeq });

        // Permute the codons of the sequence.
        var codons = SplitCodons(AnnotCodingSeq);
        var rng = new Random(20260620);
        for (int i = codons.Count - 1; i > 0; i--)
        {
            int j = rng.Next(i + 1);
            (codons[i], codons[j]) = (codons[j], codons[i]);
        }
        string shuffled = string.Concat(codons);

        GenomeAnnotator.GetCodonUsage(new[] { shuffled }).Should().BeEquivalentTo(original,
            because: "RSCU is a function of the pooled codon multiset, not codon order");

        // Splitting the same codons across two sequences, supplied in the opposite order, pools identically.
        int half = (codons.Count / 2) * 3;
        var partA = shuffled.Substring(0, half);
        var partB = shuffled.Substring(half);
        GenomeAnnotator.GetCodonUsage(new[] { partB, partA }).Should().BeEquivalentTo(original,
            because: "counts are pooled across all coding sequences before RSCU, so input order is irrelevant");
    }

    #endregion

    #region ANNOT-CODONUSAGE-001 P — each present amino acid's RSCU sums to its degeneracy (freq form sums to 1)

    [Test]
    [Description("P: for each present amino acid ΣRSCU = degeneracy nᵢ, so the normalised frequency RSCU/nᵢ sums to exactly 1 over the family.")]
    public void CodonUsage_PerAminoAcid_FrequencyFormSumsToOne()
    {
        var rscu = GenomeAnnotator.GetCodonUsage(new[] { AnnotCodingSeq });

        var present = new HashSet<char>();
        foreach (var (aa, family) in SenseFamilies())
        {
            double sum = family.Sum(c => rscu.GetValueOrDefault(c, 0.0));
            if (sum <= 0) continue; // amino acid absent from the CDS set

            present.Add(aa);
            int degeneracy = family.Count;
            sum.Should().BeApproximately(degeneracy, 1e-9,
                because: $"amino acid '{aa}' has {degeneracy} synonymous codons, so its RSCU values sum to {degeneracy}");

            double frequencySum = family.Sum(c => rscu.GetValueOrDefault(c, 0.0) / degeneracy);
            frequencySum.Should().BeApproximately(1.0, 1e-9,
                because: $"the frequency form RSCU/nᵢ sums to 1 over amino acid '{aa}'");
        }

        present.Should().Contain(new[] { 'F', 'I', 'V', 'L' }, because: "the fixture encodes those four amino acids");
        rscu.Values.Should().Contain(v => v > 1.0 + 1e-9, because: "the over-used preferred codons have RSCU > 1 — a non-vacuous fixture");
    }

    #endregion

    // ───────────────────────────────────────────────────────────────────────────
    // Unit: ANNOT-REPEAT-001 — repetitive-element detection (Annotation).
    // Checklist: docs/checklists/02_METAMORPHIC_TESTING.md, row 218.
    //
    // API under test (GenomeAnnotator.FindRepetitiveElements):
    //   Reports tandem repeats (≥ minCopies adjacent primitive-unit copies) and inverted repeats,
    //   filtered to a minimum span of minRepeatLength. Each element is (start, end, type, sequence).
    //
    // Relations (derived from the span filter + position-translation, NOT from output):
    //   • MON  (lower minLen ⇒ superset): minRepeatLength only filters candidates by span, so lowering
    //          it admits a superset of the identical elements.
    //   • SHIFT (prepend flank shifts elements): translating the sequence by a flank shifts every
    //          interior element's coordinates by the flank length (its (type, sequence) unchanged).
    // ───────────────────────────────────────────────────────────────────────────

    #region ANNOT-REPEAT-001 — Helpers

    // Two A-runs (span 4 and span 10) separated by a single G. Over {A,G} only: a reverse complement
    // would need T/C, which never occur, so there are NO inverted repeats — the result is purely the
    // two tandem A-arrays, making the relations clean.
    private const string RepeatFixture = "AAAAGAAAAAAAAAA";

    private static List<(int start, int end, string type, string sequence)> Repeats(string seq, int minLen) =>
        GenomeAnnotator.FindRepetitiveElements(seq, minRepeatLength: minLen, minCopies: 2).ToList();

    #endregion

    #region ANNOT-REPEAT-001 MON — lowering minRepeatLength admits a superset

    [Test]
    [Description("MON: minRepeatLength only filters candidates by span, so lowering it admits a superset of the identical elements.")]
    public void Repeats_LowerMinLength_SupersetOfElements()
    {
        var lenient = Repeats(RepeatFixture, 4);  // both A-runs (span 4 and 10)
        var strict = Repeats(RepeatFixture, 8);   // only the span-10 run

        lenient.Should().HaveCount(2, because: "both tandem A-runs clear a span-4 floor");
        strict.Should().NotBeEmpty(because: "the span-10 run clears the span-8 floor");

        strict.Should().BeSubsetOf(lenient, because: "raising minRepeatLength only drops whole elements, never alters them");
        strict.Should().HaveCountLessThan(lenient.Count, because: "the span-4 run is filtered at minLen = 8 but present at minLen = 4");
    }

    #endregion

    #region ANNOT-REPEAT-001 SHIFT — prepending a flank translates the elements

    [Test]
    [Description("SHIFT: prepending a flank shifts every interior repeat's coordinates by the flank length while preserving its type and sequence.")]
    public void Repeats_PrependFlank_ShiftsElements()
    {
        var baseline = Repeats(RepeatFixture, 4);
        baseline.Should().HaveCount(2, because: "the fixture has two tandem A-runs");

        // Flanks over {A,G} ending in G: G prevents the leading A-run from merging leftwards, and the
        // {A,G} alphabet keeps the whole sequence free of inverted repeats.
        foreach (string flank in new[] { "GAG", "GGAGAAG" })
        {
            int offset = flank.Length;
            var shifted = Repeats(flank + RepeatFixture, 4);

            var expected = baseline.Select(e => (e.start + offset, e.end + offset, e.type, e.sequence));

            // Elements lying inside the original region (start ≥ offset) are exactly the baseline,
            // translated by the flank length; any element the flank itself spawns sits before offset.
            shifted.Where(e => e.start >= offset).Should().BeEquivalentTo(expected,
                because: $"prepending {offset} bases shifts every interior repeat's coordinates by {offset}");
        }
    }

    #endregion
}
