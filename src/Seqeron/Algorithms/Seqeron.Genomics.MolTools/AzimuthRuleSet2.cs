using System.Buffers.Binary;
using System.Runtime.InteropServices;

namespace Seqeron.Genomics.MolTools;

// ---------------------------------------------------------------------------------------------
// Doench 2016 "Rule Set 2" / Azimuth on-target efficacy model (gradient-boosted regression trees).
//
// Unlike Doench 2014 "Rule Set 1" (a published linear model), Rule Set 2 is a trained
// scikit-learn GradientBoostingRegressor that cannot be reproduced from published coefficients.
// The two shipped models are faithful, sklearn-free reconstructions of Microsoft Research's
// Azimuth pickles, serialized to a compact binary blob by scripts/azimuth/extract_azimuth_model.py:
//
//   azimuth_rs2_nopos.bin  V3_model_nopos  (sequence only; the standard "Azimuth score" / Rule Set 2
//                                           used by CRISPOR -- needs only the 30-mer)
//   azimuth_rs2_full.bin   V3_model_full   (adds gene-context features: amino-acid cut position and
//                                           percent-peptide)
//
// Source: Doench, Fusi, Sullender, Hegde, et al. "Optimized sgRNA design to maximize activity and
//   minimize off-target effects of CRISPR-Cas9." Nat Biotechnol 34, 184-191 (2016). PMID 26780180.
//   Trained models: https://github.com/MicrosoftResearch/Azimuth (BSD-3-Clause).
//
// The featurization (incl. the Biopython Tm_NN melting-temperature features, the CPython-2.7 dict
// column ordering used at training time, and the alphabetical NGGX block) and the binary format are
// documented in the extractor script. The scores are validated end-to-end against the upstream
// regression fixture azimuth/tests/1000guides.csv (columns `truth nopos` / `truth pos`): all 947 rows
// of both models agree to 5e-10 (see the test suite).
// ---------------------------------------------------------------------------------------------

/// <summary>
/// Engine for the Doench 2016 Rule Set 2 / Azimuth on-target efficacy score. Internal; callers use
/// <see cref="CrisprDesigner.CalculateOnTargetRuleSet2(string)"/> and its gene-context overload.
/// </summary>
internal static class AzimuthRuleSet2
{
    /// <summary>One node of a regression tree (24 bytes, matches the on-disk blob layout exactly).</summary>
    [StructLayout(LayoutKind.Sequential)]
    private readonly struct Node
    {
        public readonly double ThresholdOrValue; // internal node: split threshold; leaf: pre-scaled value
        public readonly int Left;                 // child taken when feature <= threshold; -1 => leaf
        public readonly int Right;
        public readonly int Feature;              // split feature index; -1 => leaf
        public readonly int Pad;
    }

    private sealed class Model
    {
        public required Node[] Nodes;
        public required int[] TreeStart;
        public required double InitScore;
        public required int FeatureCount;
        public required bool HasGenePosition;
    }

    private const uint Magic = 0x32535241; // 'A','R','S','2'

    // Declared feature counts of the two Azimuth models. The featurizer below writes fixed column
    // offsets that are only valid for these shapes (see the column-order comment), so a blob whose
    // header declares anything else is malformed and is rejected by the loader rather than scored.
    private const int NoPosFeatureCount = 627;
    private const int FullFeatureCount = 630;

    // Gene-position column offsets within the full model's 630-column vector.
    private const int AaCutColumn = 1;
    private const int PercentPeptideBelow50Column = 607;
    private const int PercentPeptideColumn = 629;

    /// <summary>
    /// Azimuth's "no gene position" sentinel for <c>percent_peptide</c>. The reference answers a
    /// request carrying it with the sequence-only model: <c>azimuth/model_comparison.py</c> L539
    /// (<c>if np.any(percent_peptide == -1) ... model_name = 'V3_model_nopos.pickle'</c>) and L563
    /// (<c>if np.all(percent_peptide != -1) and (...)</c>, else the gene-position frame is filled
    /// with -1 and therefore unused).
    /// </summary>
    internal const double NoGenePositionSentinel = -1.0;

    private static readonly Lazy<Model> NoPosModel =
        new(() => Load("azimuth_rs2_nopos.bin", expectGenePosition: false),
            LazyThreadSafetyMode.ExecutionAndPublication);
    private static readonly Lazy<Model> FullModel =
        new(() => Load("azimuth_rs2_full.bin", expectGenePosition: true),
            LazyThreadSafetyMode.ExecutionAndPublication);

    private static Model Load(string fileName, bool expectGenePosition)
    {
        var resourceName = $"Seqeron.Genomics.MolTools.Resources.{fileName}";
        var assembly = typeof(AzimuthRuleSet2).Assembly;
        using var stream = assembly.GetManifestResourceStream(resourceName)
            ?? throw new InvalidOperationException($"Embedded Azimuth model '{resourceName}' not found.");

        var bytes = new byte[stream.Length];
        stream.ReadExactly(bytes);
        var span = (ReadOnlySpan<byte>)bytes;

        uint magic = BinaryPrimitives.ReadUInt32LittleEndian(span);
        if (magic != Magic)
            throw new InvalidOperationException($"Azimuth model '{fileName}' has bad magic 0x{magic:X8}.");
        // ushort version @4, ushort flags @6
        ushort flags = BinaryPrimitives.ReadUInt16LittleEndian(span[6..]);
        int treeCount = BinaryPrimitives.ReadInt32LittleEndian(span[8..]);
        int nodeCount = BinaryPrimitives.ReadInt32LittleEndian(span[12..]);
        int featureCount = BinaryPrimitives.ReadInt32LittleEndian(span[16..]);
        double init = BinaryPrimitives.ReadDoubleLittleEndian(span[20..]);
        // double learningRate @28 (informational; leaf values are already pre-scaled by it)

        // A malformed blob must fail as an InvalidOperationException (like the magic check above),
        // not as an IndexOutOfRange/ArgumentOutOfRange from the featurizer or the node slice.
        bool hasGenePosition = (flags & 1) != 0;
        if (hasGenePosition != expectGenePosition)
            throw new InvalidOperationException(
                $"Azimuth model '{fileName}' declares gene-position flag {hasGenePosition}, expected {expectGenePosition}.");
        int expectedFeatureCount = expectGenePosition ? FullFeatureCount : NoPosFeatureCount;
        if (featureCount != expectedFeatureCount)
            throw new InvalidOperationException(
                $"Azimuth model '{fileName}' declares {featureCount} features, expected {expectedFeatureCount}.");
        if (treeCount <= 0 || nodeCount <= 0)
            throw new InvalidOperationException(
                $"Azimuth model '{fileName}' declares {treeCount} trees / {nodeCount} nodes.");

        int offset = 36;
        if (span.Length < offset + treeCount * 4 + nodeCount * 24)
            throw new InvalidOperationException(
                $"Azimuth model '{fileName}' is truncated: {span.Length} bytes for {treeCount} trees / {nodeCount} nodes.");

        var treeStart = new int[treeCount];
        for (int i = 0; i < treeCount; i++)
            treeStart[i] = BinaryPrimitives.ReadInt32LittleEndian(span[(offset + i * 4)..]);
        offset += treeCount * 4;

        // Zero-copy reinterpretation of the node region, copied once into a managed array.
        var nodes = MemoryMarshal.Cast<byte, Node>(span.Slice(offset, nodeCount * 24)).ToArray();

        return new Model
        {
            Nodes = nodes,
            TreeStart = treeStart,
            InitScore = init,
            FeatureCount = featureCount,
            HasGenePosition = hasGenePosition,
        };
    }

    private static double Predict(Model model, ReadOnlySpan<double> features)
    {
        var nodes = model.Nodes;
        double score = model.InitScore;
        foreach (int root in model.TreeStart)
        {
            int n = root;
            while (nodes[n].Feature != -1)
            {
                ref readonly var node = ref nodes[n];
                n = features[node.Feature] <= node.ThresholdOrValue ? node.Left : node.Right;
            }
            score += nodes[n].ThresholdOrValue; // pre-scaled leaf value
        }
        return score;
    }

    // ----------------------------------------------------------------------------------- public

    /// <summary>Sequence-only Rule Set 2 score for a 30-mer (validated A/C/G/T, N<c>GG</c> PAM).</summary>
    public static double Score(string context30Mer) => ScoreNoPos(Validate(context30Mer));

    /// <summary>
    /// Rule Set 2 score using the gene-context (full) model -- except for Azimuth's
    /// <see cref="NoGenePositionSentinel"/>, which the reference answers with the sequence-only model
    /// (<c>model_comparison.py</c> L539/L563); see <see cref="NoGenePositionSentinel"/>.
    /// </summary>
    public static double Score(string context30Mer, int aminoAcidCutPosition, double percentPeptide)
    {
        // The reference validates the sequence first (model_comparison.py L522-L524), then the
        // gene-position arguments (L531-L534).
        var seq = Validate(context30Mer);

        // The reference only asserts percent_peptide is a real number (L534); the documented domain
        // of the feature is a percentage, and -1 is its own "no gene position" sentinel. Anything
        // else (including NaN) is a caller error rather than a silently-extrapolated tree walk.
        if (percentPeptide != NoGenePositionSentinel && !(percentPeptide >= 0.0 && percentPeptide <= 100.0))
            throw new ArgumentOutOfRangeException(
                nameof(percentPeptide), percentPeptide,
                "Rule Set 2 percent-peptide must be in [0, 100], or -1 for Azimuth's "
                + "\"no gene position\" sentinel (which selects the sequence-only model).");

        if (percentPeptide == NoGenePositionSentinel)
            return ScoreNoPos(seq);

        // "Amino Acid Cut position" is a residue count of the cut within the protein. The reference asserts
        // only that it is a real number (model_comparison.py, `np.isreal(aa_cut)`) and has no range check;
        // its own 1000-guide regression fixture (azimuth/tests/1000guides.csv, column "AA cut") holds
        // integers in [0, 499] -- 0 occurs, so the bound is >= 0, not >= 1. A negative count is not a
        // position; the reference's -1 sentinel belongs to percent_peptide (above, where the aa_cut is
        // then ignored), so a real percent_peptide paired with a negative aa_cut is a caller error rather
        // than a silently-extrapolated tree walk. There is no upper bound in the reference (proteins vary).
        if (aminoAcidCutPosition < 0)
            throw new ArgumentOutOfRangeException(
                nameof(aminoAcidCutPosition), aminoAcidCutPosition,
                "Rule Set 2 amino-acid cut position must be >= 0 (a residue count; Azimuth's own fixture spans "
                + "0..499). The \"no gene position\" sentinel -1 belongs to percentPeptide, which then selects "
                + "the sequence-only model and ignores this argument.");

        var model = FullModel.Value;
        if (!model.HasGenePosition || model.FeatureCount != FullFeatureCount)
            throw new InvalidOperationException(
                $"Azimuth full model has {model.FeatureCount} features / gene-position flag "
                + $"{model.HasGenePosition}; the gene-position columns are undefined for that shape.");
        Span<double> f = new double[model.FeatureCount];
        FeaturizeNoPos(seq, f, fullModel: true);
        FillGenePosition(f, aminoAcidCutPosition, percentPeptide);
        return Predict(model, f);
    }

    private static double ScoreNoPos(string seq)
    {
        var model = NoPosModel.Value;
        Span<double> f = new double[model.FeatureCount];
        FeaturizeNoPos(seq, f);
        return Predict(model, f);
    }

    private static string Validate(string context30Mer)
    {
        // A null reference is a missing argument; an empty string is a wrong value, and the length
        // check below reports it the same way as every other malformed context (the reference also
        // separates these: `assert isinstance(seq[0], basestring)` vs `assert len(seq[0]) > 0`,
        // model_comparison.py L523-L524).
        ArgumentNullException.ThrowIfNull(context30Mer);
        var seq = SequenceExtensions.ToUpperAscii(context30Mer)!;
        if (seq.Length != 30)
            throw new ArgumentException(
                $"Rule Set 2 requires a 30-nt context (4 upstream + 20 protospacer + 3 PAM + 3 downstream); got {seq.Length}.",
                nameof(context30Mer));
        foreach (var c in seq)
            if (c is not ('A' or 'C' or 'G' or 'T'))
                throw new ArgumentException(
                    $"Rule Set 2 context must contain only A/C/G/T; found '{c}'.", nameof(context30Mer));
        if (seq[25] != 'G' || seq[26] != 'G')
            throw new ArgumentException(
                "Rule Set 2 expects an SpCas9 NGG PAM at offsets 25-26 of the 30-nt context.",
                nameof(context30Mer));
        return seq;
    }

    // ------------------------------------------------------------------- featurization
    //
    // Column order is the CPython-2.7 dict iteration order used by azimuth at training time (see the
    // extractor). nopos (627):
    //   gc_count(1) | pd_Order2(464) | pd_Order1(120) | gc_above_10(1) | pi_Order1(4) | pi_Order2(16)
    //   | Tm(4) | gc_below_10(1) | NGGX(16)
    // full (630): same nucleotide/Tm blocks, with the three gene-position scalars interleaved:
    //   gc_count(1) | AA_cut(1) | pd_Order2(464) | pd_Order1(120) | gc_above_10(1) | pi_Order1(4)
    //   | pi_Order2(16) | pct_pep<50%(1) | Tm(4) | gc_below_10(1) | NGGX(16) | pct_pep(1)

    private static int NucIndex(char c) => c switch { 'A' => 0, 'T' => 1, 'C' => 2, 'G' => 3, _ => -1 };

    // The NGGX block is the one block that is NOT in azimuth's A,T,C,G alphabet order but in ALPHABETICAL
    // A,C,G,T order. Upstream NGGX_interaction_feature (azimuth/features/featurization.py L130-L145)
    // starts from an empty `pandas.DataFrame()` and `pandas.concat`s one 16-row one-hot frame per guide
    // along axis=1; the index union of non-equal indexes is what sorts the labels (the audit traced it to
    // pandas 0.19.2 `_union_indexes` -> `lib.fast_unique_multiple_list`; that pandas is not installable
    // here, so the mechanism is the explanation and the measurement below is the evidence), whereas the
    // pd_*/pi_* blocks are built by `Series.apply`, whose indexes are all equal and so keep A,T,C,G.
    // Measured against upstream's own fixture (azimuth/tests/1000guides.csv, `truth nopos` /
    // `truth pos`) with the shipped pickles: A,C,G,T gives 947/947 rows within 5e-10 for both models;
    // A,T,C,G leaves 362/947 (nopos) and 310/947 (full) rows off by up to 0.038 / 0.051.
    private static int NggxIndex(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, _ => -1 };

    /// <summary>Writes all sequence-derived blocks. Returns nothing meaningful; gene-position slots are
    /// left for <see cref="FillGenePosition"/> in the full model.</summary>
    private static void FeaturizeNoPos(string seq, Span<double> f, bool fullModel = false)
    {
        // Offsets differ between models because the gene-position scalars shift the nucleotide blocks.
        // GC over the 20-nt protospacer = seq[4:24], counted with the canonical GC primitive.
        int gcCount = seq.AsSpan(4, 20).CountGcAndValidNucleotides().GcCount;

        int oPd2, oPd1, oGcAbove, oPi1, oPi2, oTm, oGcBelow, oNggx, oGcCount;
        if (!fullModel)
        {
            oGcCount = 0; oPd2 = 1; oPd1 = 465; oGcAbove = 585; oPi1 = 586;
            oPi2 = 590; oTm = 606; oGcBelow = 610; oNggx = 611;
        }
        else
        {
            // gc_count(0), AA_cut(1), pd2(2), pd1(466), gc_above(586), pi1(587), pi2(591),
            // pct<50%(607), Tm(608), gc_below(612), NGGX(613), pct(629)
            oGcCount = 0; oPd2 = 2; oPd1 = 466; oGcAbove = 586; oPi1 = 587;
            oPi2 = 591; oTm = 608; oGcBelow = 612; oNggx = 613;
        }

        f[oGcCount] = gcCount;
        f[oGcAbove] = gcCount > 10 ? 1.0 : 0.0;
        f[oGcBelow] = gcCount < 10 ? 1.0 : 0.0;

        // Position-dependent one-hot (order 1: 30 positions x 4; order 2: 29 positions x 16).
        for (int pos = 0; pos < 30; pos++)
            f[oPd1 + pos * 4 + NucIndex(seq[pos])] = 1.0;
        for (int pos = 0; pos < 29; pos++)
        {
            int di = NucIndex(seq[pos]) * 4 + NucIndex(seq[pos + 1]);
            f[oPd2 + pos * 16 + di] = 1.0;
        }

        // Position-independent counts.
        for (int pos = 0; pos < 30; pos++)
            f[oPi1 + NucIndex(seq[pos])] += 1.0;
        for (int pos = 0; pos < 29; pos++)
            f[oPi2 + NucIndex(seq[pos]) * 4 + NucIndex(seq[pos + 1])] += 1.0;

        // NGGX: one-hot of the dinucleotide formed by the PAM 'N' (pos 24) and the +1 base (pos 27),
        // in alphabetical A,C,G,T order (see NggxIndex).
        f[oNggx + NggxIndex(seq[24]) * 4 + NggxIndex(seq[27])] = 1.0;

        // Melting temperatures: whole 30-mer and three sub-segments of the protospacer.
        f[oTm + 0] = MeltingTemp(seq, 0, 30);
        f[oTm + 1] = MeltingTemp(seq, 19, 24);
        f[oTm + 2] = MeltingTemp(seq, 11, 19);
        f[oTm + 3] = MeltingTemp(seq, 6, 11);
    }

    /// <summary>Writes the three gene-position columns of the full model. The loader has already
    /// rejected any blob whose declared shape differs from <see cref="FullFeatureCount"/>, so these
    /// offsets are in range.</summary>
    private static void FillGenePosition(Span<double> f, int aminoAcidCutPosition, double percentPeptide)
    {
        f[AaCutColumn] = aminoAcidCutPosition;
        // "Percent Peptide <50%" is `feature_sets["Percent Peptide"] < 50`
        // (azimuth/features/featurization.py L48).
        f[PercentPeptideBelow50Column] = percentPeptide < 50 ? 1.0 : 0.0;
        f[PercentPeptideColumn] = percentPeptide;
    }

    // ------------------------------------------------------------------- Tm (Biopython Tm_NN)
    //
    // Nearest-neighbor melting temperature, DNA_NN3 (Allawi & SantaLucia 1997), salt-correction
    // method 5, dnac1=dnac2=25 nM, Na=50 mM -- the exact parameters behind the Biopython call azimuth
    // makes. Upstream calls the OBSOLETE MeltingTemp.Tm_staluc (azimuth/features/featurization.py
    // L381-L384, Tm_feature), which Biopython implements as
    //   Tm_staluc(s, dnac=50, saltc=50) -> Tm_NN(s, dnac1=dnac/2.0, dnac2=dnac/2.0, Na=saltc)
    // (Bio/SeqUtils/MeltingTemp.py L1087/L1127 at tag biopython-176, where the function was last
    // shipped; it is gone in 1.88), i.e. Tm_NN's own defaults nn_table=DNA_NN3, saltcorr=5,
    // dnac1=dnac2=25 nM, Na=50 mM -- bit-equal on execution (1.76 Tm_staluc vs 1.88 Tm_NN, verified
    // for the five segments locked in the tests plus the 30-mer AAAAAAAAAAAAAAAAAAATGCAGCGGGAG
    // = 55.39603328329298) and the defaults of the canonical
    // ThermoConstants.CalculateNearestNeighborTm. The previously
    // inlined DNA_NN3 table + salt/concentration arithmetic was verified BIT-IDENTICAL to that
    // canonical implementation over 821,840 comparisons (200,000 random 30-mers x the four azimuth
    // windows [0,30) [19,24) [11,19) [6,11), plus every 2..7-mer exhaustively; 0 bitwise
    // differences), and both reproduce real Biopython 1.88 exactly, e.g.
    //   Tm_NN("ACGTACGTACGTACGTACGTACGTAGGACG") = 62.55983015288416
    //   Tm_NN("TAGG") = -67.90079099878676, Tm_NN("GCGCGCGC") = 36.45901389113385
    // so the duplicate table was removed in favour of the canonical call.

    private static double MeltingTemp(string seq, int start, int end) =>
        ThermoConstants.CalculateNearestNeighborTm(seq.Substring(start, end - start));
}
