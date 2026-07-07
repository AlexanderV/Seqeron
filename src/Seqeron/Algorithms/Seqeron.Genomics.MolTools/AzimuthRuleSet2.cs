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
// The featurization (incl. the Biopython Tm_NN melting-temperature features and the CPython-2.7 dict
// column ordering used at training time) and the binary format are documented in the extractor script
// and validated end-to-end against that verified Python reference (see the test suite).
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

    private static readonly Lazy<Model> NoPosModel =
        new(() => Load("azimuth_rs2_nopos.bin"), LazyThreadSafetyMode.ExecutionAndPublication);
    private static readonly Lazy<Model> FullModel =
        new(() => Load("azimuth_rs2_full.bin"), LazyThreadSafetyMode.ExecutionAndPublication);

    private static Model Load(string fileName)
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

        int offset = 36;
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
            HasGenePosition = (flags & 1) != 0,
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
    public static double Score(string context30Mer)
    {
        var seq = Validate(context30Mer);
        var model = NoPosModel.Value;
        Span<double> f = new double[model.FeatureCount];
        FeaturizeNoPos(seq, f);
        return Predict(model, f);
    }

    /// <summary>Rule Set 2 score using the gene-context (full) model.</summary>
    public static double Score(string context30Mer, int aminoAcidCutPosition, double percentPeptide)
    {
        var seq = Validate(context30Mer);
        var model = FullModel.Value;
        Span<double> f = new double[model.FeatureCount];
        FeaturizeNoPos(seq, f, fullModel: true);
        FillGenePosition(f, aminoAcidCutPosition, percentPeptide);
        return Predict(model, f);
    }

    private static string Validate(string context30Mer)
    {
        if (string.IsNullOrEmpty(context30Mer))
            throw new ArgumentNullException(nameof(context30Mer));
        var seq = context30Mer.ToUpperInvariant();
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

    /// <summary>Writes all sequence-derived blocks. Returns nothing meaningful; gene-position slots are
    /// left for <see cref="FillGenePosition"/> in the full model.</summary>
    private static void FeaturizeNoPos(string seq, Span<double> f, bool fullModel = false)
    {
        // Offsets differ between models because the gene-position scalars shift the nucleotide blocks.
        int gcCount = 0;
        for (int i = 4; i < 24; i++)
            if (seq[i] is 'G' or 'C') gcCount++;

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

        // NGGX: one-hot of the dinucleotide formed by the PAM 'N' (pos 24) and the +1 base (pos 27).
        f[oNggx + NucIndex(seq[24]) * 4 + NucIndex(seq[27])] = 1.0;

        // Melting temperatures: whole 30-mer and three sub-segments of the protospacer.
        f[oTm + 0] = MeltingTemp(seq, 0, 30);
        f[oTm + 1] = MeltingTemp(seq, 19, 24);
        f[oTm + 2] = MeltingTemp(seq, 11, 19);
        f[oTm + 3] = MeltingTemp(seq, 6, 11);
    }

    private static void FillGenePosition(Span<double> f, int aminoAcidCutPosition, double percentPeptide)
    {
        f[1] = aminoAcidCutPosition;             // Amino Acid Cut position
        f[607] = percentPeptide < 50 ? 1.0 : 0.0; // Percent Peptide <50%
        f[629] = percentPeptide;                 // Percent Peptide
    }

    // ------------------------------------------------------------------- Tm (Biopython Tm_NN)
    //
    // Nearest-neighbor melting temperature, DNA_NN3 (Allawi & SantaLucia 1997), salt-correction
    // method 5, dnac1=dnac2=25 nM, Na=50 mM -- the exact parameters azimuth passes to Biopython's
    // MeltingTemp.Tm_NN. Verified bit-comparable to real Biopython for whole and short AT-rich segments.

    private static double MeltingTemp(string seq, int start, int end)
    {
        int len = end - start;
        double dh = 0.0, ds = 0.0; // DNA_NN3 'init' is (0, 0)

        char first = seq[start], last = seq[end - 1];
        int at = (first is 'A' or 'T' ? 1 : 0) + (last is 'A' or 'T' ? 1 : 0);
        int gc = (first is 'G' or 'C' ? 1 : 0) + (last is 'G' or 'C' ? 1 : 0);
        dh += 2.3 * at + 0.1 * gc;   // init_A/T (2.3, 4.1), init_G/C (0.1, -2.8)
        ds += 4.1 * at + -2.8 * gc;

        for (int i = start; i < end - 1; i++)
        {
            var (h, s) = NearestNeighbor(seq[i], seq[i + 1]);
            dh += h; ds += s;
        }

        ds += 0.368 * (len - 1) * Math.Log(50e-3);     // salt correction, method 5
        const double k = (25.0 - 25.0 / 2.0) * 1e-9;   // (dnac1 - dnac2/2) * 1e-9
        return (1000.0 * dh) / (ds + 1.987 * Math.Log(k)) - 273.15;
    }

    /// <summary>DNA_NN3 nearest-neighbor (deltaH, deltaS); pair and its reverse map to the same value.</summary>
    private static (double H, double S) NearestNeighbor(char a, char b)
    {
        // Key is the Watson-Crick pair "XY/X'Y'"; only one orientation is stored, so we also try
        // the reverse, exactly as Biopython does.
        return (a, b) switch
        {
            ('A', 'A') or ('T', 'T') => (-7.9, -22.2),
            ('A', 'T') => (-7.2, -20.4),
            ('T', 'A') => (-7.2, -21.3),
            ('C', 'A') or ('T', 'G') => (-8.5, -22.7),
            ('G', 'T') or ('A', 'C') => (-8.4, -22.4),
            ('C', 'T') or ('A', 'G') => (-7.8, -21.0),
            ('G', 'A') or ('T', 'C') => (-8.2, -22.2),
            ('C', 'G') => (-10.6, -27.2),
            ('G', 'C') => (-9.8, -24.4),
            ('G', 'G') or ('C', 'C') => (-8.0, -19.9),
            _ => throw new ArgumentException($"Tm: non-ACGT dinucleotide '{a}{b}'."),
        };
    }
}
