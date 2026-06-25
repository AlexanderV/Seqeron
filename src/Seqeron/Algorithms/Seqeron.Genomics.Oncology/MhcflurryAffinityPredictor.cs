using System;
using System.Buffers.Binary;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Reflection;

namespace Seqeron.Genomics.Oncology;

/// <summary>
/// A faithful C# port of MHCflurry's Class I pan-allele binding-<b>affinity</b> (BA) neural-network
/// predictor (the <c>Class1AffinityPredictor</c> / <c>Class1NeuralNetwork</c> forward pass).
/// <para>
/// Given a peptide and an HLA allele, the predictor reproduces MHCflurry's pipeline:
/// </para>
/// <list type="number">
/// <item><description>encode the peptide with the BLOSUM62 <c>left_pad_centered_right_pad</c> scheme
/// (3 copies — left-aligned, centred, right-aligned — at <c>max_length</c>=15, giving 45×21 = 945 values);</description></item>
/// <item><description>encode the allele's 37-residue MHC-pocket <i>pseudosequence</i> with BLOSUM62 (37×21 = 777 values);</description></item>
/// <item><description>concatenate (peptide-flat then allele-flat) and run the feed-forward dense network
/// (<c>tanh</c> hidden layers, optional MHCflurry "with-skip-connections" densenet wiring, <c>sigmoid</c> output);</description></item>
/// <item><description>map the [0,1] output to IC50 (nM) via <c>IC50 = 50000^(1 − output)</c> (<c>regression_target.to_ic50</c>);</description></item>
/// <item><description>aggregate an ensemble of networks by the geometric mean of the per-network IC50s
/// (the default <c>centrality_measure="mean"</c> applied in log space: <c>exp(mean(log(ic50)))</c>).</description></item>
/// </list>
/// <para>
/// <b>Source.</b> MHCflurry (Apache License 2.0), github.com/openvax/mhcflurry — modules
/// <c>amino_acid.py</c> (BLOSUM62, amino-acid order <c>ACDEFGHIKLMNPQRSTVWYX</c>),
/// <c>encodable_sequences.py</c> (<c>left_pad_centered_right_pad</c>), <c>class1_neural_network.py</c>
/// (<c>make_network</c>, <c>predict</c>), <c>regression_target.py</c> (<c>to_ic50</c>),
/// <c>ensemble_centrality.py</c> (<c>mean</c>), and the <c>models_class1_pan</c> (release 20200610) network
/// architecture/weights. O'Donnell, Rubinsteyn, Laserson, "MHCflurry 2.0", Cell Systems 11(1):42-48.e7 (2020),
/// doi:10.1016/j.cels.2020.06.010. The bundled pseudosequence table and the Apache NOTICE are under
/// <c>Resources/</c>.
/// </para>
/// <para>
/// <b>Bundling.</b> The 37-residue allele pseudosequence table is bundled as an embedded resource. The trained
/// weight ensemble is large (~80 MB across the 10-network <c>models_class1_pan</c> ensemble), so the weights are
/// <b>not</b> embedded; they are loaded from a caller-supplied MHCflurry weight pack stream (see
/// <see cref="LoadWeightPack(Stream)"/> / the pack format documented there). This keeps the existing
/// <see cref="OncologyAnalyzer.ClassifyBindingAffinity(double)"/> classification and defaults unchanged; the
/// neural predictor is strictly opt-in.
/// </para>
/// </summary>
public static class MhcflurryAffinityPredictor
{
    #region Constants (MHCflurry encoding + transform)

    /// <summary>
    /// MHCflurry amino-acid alphabet order. Source: <c>amino_acid.py</c> —
    /// <c>COMMON_AMINO_ACIDS</c> sorted, then <c>"X"</c> (Unknown) appended; the resulting
    /// <c>AMINO_ACIDS</c> order used to index <c>BLOSUM62_MATRIX</c> is <c>ACDEFGHIKLMNPQRSTVWYX</c>.
    /// </summary>
    public const string AminoAcidOrder = "ACDEFGHIKLMNPQRSTVWYX";

    /// <summary>Unknown / padding amino-acid character. Source: <c>EncodableSequences.unknown_character = "X"</c>.</summary>
    public const char UnknownAminoAcid = 'X';

    /// <summary>Vector-encoding width (20 amino acids + X). Source: <c>amino_acid.vector_encoding_length("BLOSUM62") == 21</c>.</summary>
    public const int EncodingWidth = 21;

    /// <summary>
    /// Maximum supported peptide length for the pan-allele peptide encoding. Source: <c>models_class1_pan</c>
    /// hyperparameter <c>peptide_encoding.max_length = 15</c>.
    /// </summary>
    public const int PeptideMaxLength = 15;

    /// <summary>
    /// Minimum supported peptide length for the <c>left_pad_centered_right_pad</c> encoding. Source:
    /// <c>encodable_sequences.py</c> — "We arbitrarily set a minimum length of 5" for this alignment method.
    /// </summary>
    public const int PeptideMinLength = 5;

    /// <summary>Length of an MHC pseudosequence in the bundled table (residues). Source: <c>allele_sequences.csv</c> (all rows length 37).</summary>
    public const int PseudosequenceLength = 37;

    /// <summary>
    /// Flattened peptide-encoding length: 3 × <see cref="PeptideMaxLength"/> × <see cref="EncodingWidth"/> = 945.
    /// Source: <c>left_pad_centered_right_pad</c> produces <c>3 × max_length</c> positions, each a 21-vector.
    /// </summary>
    public const int PeptideFlatLength = 3 * PeptideMaxLength * EncodingWidth;

    /// <summary>Flattened allele-encoding length: 37 × 21 = 777. Source: <c>AlleleEncoding.allele_representations</c> flattened.</summary>
    public const int AlleleFlatLength = PseudosequenceLength * EncodingWidth;

    /// <summary>Network input width = peptide-flat ‖ allele-flat = 945 + 777 = 1722. Source: <c>make_network</c> concatenate.</summary>
    public const int NetworkInputLength = PeptideFlatLength + AlleleFlatLength;

    /// <summary>
    /// IC50 ceiling used in the regression-target transform. Source: <c>regression_target.to_ic50</c>:
    /// <c>ic50 = max_ic50 ** (1 - x)</c> with <c>max_ic50 = 50000.0</c>.
    /// </summary>
    public const double MaxIc50Nm = 50000.0;

    private const string PseudosequenceResourceName =
        "Seqeron.Genomics.Oncology.Resources.mhcflurry.allele_sequences.csv";

    // Weight-pack binary format magic + version. See LoadWeightPack for the layout.
    private static readonly byte[] WeightPackMagic = "MHCF"u8.ToArray();
    private const int WeightPackVersion = 1;

    // Topology codes in the weight pack (mirrors MHCflurry hyperparameter "topology").
    private const int TopologyFeedforward = 0;          // make_network: densenet_layers = None
    private const int TopologyWithSkipConnections = 1;  // make_network: densenet topology

    #endregion

    #region BLOSUM62 (amino_acid.py, order ACDEFGHIKLMNPQRSTVWYX)

    // BLOSUM62_MATRIX from mhcflurry/amino_acid.py, re-indexed to AminoAcidOrder. Row i = encoding vector for
    // residue AminoAcidOrder[i]. The X row is all-zero except the X column (=1), and all other rows have 0 in
    // the X column (per the source matrix).
    private static readonly sbyte[][] Blosum62 =
    [
        /* A */ [ 4, 0,-2,-1,-2, 0,-2,-1,-1,-1,-1,-2,-1,-1,-1, 1, 0, 0,-3,-2, 0],
        /* C */ [ 0, 9,-3,-4,-2,-3,-3,-1,-3,-1,-1,-3,-3,-3,-3,-1,-1,-1,-2,-2, 0],
        /* D */ [-2,-3, 6, 2,-3,-1,-1,-3,-1,-4,-3, 1,-1, 0,-2, 0,-1,-3,-4,-3, 0],
        /* E */ [-1,-4, 2, 5,-3,-2, 0,-3, 1,-3,-2, 0,-1, 2, 0, 0,-1,-2,-3,-2, 0],
        /* F */ [-2,-2,-3,-3, 6,-3,-1, 0,-3, 0, 0,-3,-4,-3,-3,-2,-2,-1, 1, 3, 0],
        /* G */ [ 0,-3,-1,-2,-3, 6,-2,-4,-2,-4,-3, 0,-2,-2,-2, 0,-2,-3,-2,-3, 0],
        /* H */ [-2,-3,-1, 0,-1,-2, 8,-3,-1,-3,-2, 1,-2, 0, 0,-1,-2,-3,-2, 2, 0],
        /* I */ [-1,-1,-3,-3, 0,-4,-3, 4,-3, 2, 1,-3,-3,-3,-3,-2,-1, 3,-3,-1, 0],
        /* K */ [-1,-3,-1, 1,-3,-2,-1,-3, 5,-2,-1, 0,-1, 1, 2, 0,-1,-2,-3,-2, 0],
        /* L */ [-1,-1,-4,-3, 0,-4,-3, 2,-2, 4, 2,-3,-3,-2,-2,-2,-1, 1,-2,-1, 0],
        /* M */ [-1,-1,-3,-2, 0,-3,-2, 1,-1, 2, 5,-2,-2, 0,-1,-1,-1, 1,-1,-1, 0],
        /* N */ [-2,-3, 1, 0,-3, 0, 1,-3, 0,-3,-2, 6,-2, 0, 0, 1, 0,-3,-4,-2, 0],
        /* P */ [-1,-3,-1,-1,-4,-2,-2,-3,-1,-3,-2,-2, 7,-1,-2,-1,-1,-2,-4,-3, 0],
        /* Q */ [-1,-3, 0, 2,-3,-2, 0,-3, 1,-2, 0, 0,-1, 5, 1, 0,-1,-2,-2,-1, 0],
        /* R */ [-1,-3,-2, 0,-3,-2, 0,-3, 2,-2,-1, 0,-2, 1, 5,-1,-1,-3,-3,-2, 0],
        /* S */ [ 1,-1, 0, 0,-2, 0,-1,-2, 0,-2,-1, 1,-1, 0,-1, 4, 1,-2,-3,-2, 0],
        /* T */ [ 0,-1,-1,-1,-2,-2,-2,-1,-1,-1,-1, 0,-1,-1,-1, 1, 5, 0,-2,-2, 0],
        /* V */ [ 0,-1,-3,-2,-1,-3,-3, 3,-2, 1, 1,-3,-2,-2,-3,-2, 0, 4,-3,-1, 0],
        /* W */ [-3,-2,-4,-3, 1,-2,-2,-3,-3,-2,-1,-4,-4,-2,-3,-3,-2,-3,11, 2, 0],
        /* Y */ [-2,-2,-3,-2, 3,-3, 2,-1,-2,-1,-1,-2,-3,-1,-2,-2,-2,-1, 2, 7, 0],
        /* X */ [ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    ];

    private static readonly int[] AminoAcidIndex = BuildAminoAcidIndex();

    private static int[] BuildAminoAcidIndex()
    {
        // Map ASCII upper/lower-case letters to the row index in AminoAcidOrder; unmapped -> X index.
        int xIndex = AminoAcidOrder.IndexOf(UnknownAminoAcid);
        var map = new int[128];
        Array.Fill(map, xIndex);
        for (int i = 0; i < AminoAcidOrder.Length; i++)
        {
            char c = AminoAcidOrder[i];
            map[char.ToUpperInvariant(c)] = i;
            map[char.ToLowerInvariant(c)] = i;
        }
        return map;
    }

    private static int IndexOfResidue(char residue)
    {
        // Non-canonical / out-of-range residues are treated as X, mirroring
        // EncodableSequences with allow_unsupported_amino_acids and the X fallback for alleles.
        return residue < 128 ? AminoAcidIndex[residue] : AminoAcidIndex[UnknownAminoAcid];
    }

    #endregion

    #region Pseudosequence table (bundled, Apache-2.0)

    private static readonly Lazy<IReadOnlyDictionary<string, string>> PseudosequenceTable =
        new(LoadPseudosequenceTable);

    /// <summary>
    /// Loads the bundled MHCflurry allele → 37-residue pseudosequence table (release 20200610), parsing the
    /// embedded <c>allele_sequences.csv</c> (comment header lines beginning with <c>#</c> are skipped). Apache-2.0;
    /// see <c>Resources/MHCFLURRY_NOTICE.txt</c>.
    /// </summary>
    /// <returns>A read-only map from MHCflurry allele name (e.g. <c>HLA-A*02:01</c>) to its pseudosequence.</returns>
    public static IReadOnlyDictionary<string, string> GetAllelePseudosequences() => PseudosequenceTable.Value;

    /// <summary>
    /// Gets the 37-residue MHC pseudosequence for <paramref name="allele"/> from the bundled table.
    /// </summary>
    /// <param name="allele">MHCflurry allele name exactly as keyed in the table (e.g. <c>HLA-A*02:01</c>).</param>
    /// <returns>The 37-residue pseudosequence.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="allele"/> is null.</exception>
    /// <exception cref="KeyNotFoundException">The allele is not present in the bundled table.</exception>
    public static string GetPseudosequence(string allele)
    {
        ArgumentNullException.ThrowIfNull(allele);
        if (!PseudosequenceTable.Value.TryGetValue(allele, out string? seq))
        {
            throw new KeyNotFoundException(
                $"Allele '{allele}' is not in the bundled MHCflurry pseudosequence table. " +
                "Use the exact MHCflurry allele name (e.g. 'HLA-A*02:01').");
        }
        return seq;
    }

    private static IReadOnlyDictionary<string, string> LoadPseudosequenceTable()
    {
        var asm = typeof(MhcflurryAffinityPredictor).Assembly;
        using Stream stream = asm.GetManifestResourceStream(PseudosequenceResourceName)
            ?? throw new InvalidOperationException(
                $"Bundled MHCflurry pseudosequence resource '{PseudosequenceResourceName}' was not found.");
        using var reader = new StreamReader(stream);

        var table = new Dictionary<string, string>(StringComparer.Ordinal);
        bool headerSeen = false;
        string? line;
        while ((line = reader.ReadLine()) is not null)
        {
            if (line.Length == 0 || line[0] == '#')
            {
                continue;
            }
            int comma = line.IndexOf(',');
            if (comma <= 0)
            {
                continue;
            }
            string key = line[..comma];
            string value = line[(comma + 1)..].Trim();
            if (!headerSeen)
            {
                // The first non-comment line is the "allele,sequence" header.
                headerSeen = true;
                if (string.Equals(key, "allele", StringComparison.OrdinalIgnoreCase))
                {
                    continue;
                }
            }
            table[key] = value;
        }
        return table;
    }

    #endregion

    #region Encoders (encodable_sequences.py / allele_encoding.py)

    /// <summary>
    /// Encodes a peptide into the flattened BLOSUM62 <c>left_pad_centered_right_pad</c> representation used by
    /// the pan-allele network: three copies of the peptide at <see cref="PeptideMaxLength"/> positions —
    /// left-aligned, centred, right-aligned — each residue a 21-wide BLOSUM62 vector, in row-major
    /// (position, channel) order, giving <see cref="PeptideFlatLength"/> values. Padding positions use the X
    /// (all-but-X-zero) vector. Source: <c>encodable_sequences.py</c> (<c>left_pad_centered_right_pad</c>) +
    /// <c>amino_acid.fixed_vectors_encoding</c>.
    /// </summary>
    /// <param name="peptide">Peptide sequence (length in [<see cref="PeptideMinLength"/>, <see cref="PeptideMaxLength"/>]).</param>
    /// <returns>A <see cref="PeptideFlatLength"/>-length vector of BLOSUM62 values.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="peptide"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The peptide length is unsupported by this encoding.</exception>
    public static double[] EncodePeptide(string peptide)
    {
        ArgumentNullException.ThrowIfNull(peptide);
        int len = peptide.Length;
        if (len < PeptideMinLength || len > PeptideMaxLength)
        {
            throw new ArgumentOutOfRangeException(
                nameof(peptide), len,
                $"Peptide length must be in [{PeptideMinLength}, {PeptideMaxLength}] for the MHCflurry pan-allele encoding.");
        }

        // Index-encode the peptide first.
        Span<int> idx = stackalloc int[len];
        for (int i = 0; i < len; i++)
        {
            idx[i] = IndexOfResidue(peptide[i]);
        }
        int xIdx = AminoAcidIndex[UnknownAminoAcid];

        // Build a (3*max_length) integer layout filled with X, then place the three copies.
        int positions = 3 * PeptideMaxLength;
        Span<int> layout = stackalloc int[positions];
        for (int i = 0; i < positions; i++)
        {
            layout[i] = xIdx;
        }

        // Block 0 [0, max_length): left-aligned  -> result[:length]
        for (int i = 0; i < len; i++)
        {
            layout[i] = idx[i];
        }
        // Block 2 [2*max_length, 3*max_length): right-aligned -> result[-length:] (of the WHOLE 3*max_length array)
        for (int i = 0; i < len; i++)
        {
            layout[positions - len + i] = idx[i];
        }
        // Block 1 [max_length, 2*max_length): centred. center_left_padding = floor((max_length - length)/2);
        // offset = max_length + center_left_padding.
        int centerLeftPadding = (PeptideMaxLength - len) / 2; // floor for non-negative
        int centerOffset = PeptideMaxLength + centerLeftPadding;
        for (int i = 0; i < len; i++)
        {
            layout[centerOffset + i] = idx[i];
        }

        // Vector-encode the integer layout via BLOSUM62.
        var result = new double[PeptideFlatLength];
        int w = 0;
        for (int p = 0; p < positions; p++)
        {
            sbyte[] row = Blosum62[layout[p]];
            for (int c = 0; c < EncodingWidth; c++)
            {
                result[w++] = row[c];
            }
        }
        return result;
    }

    /// <summary>
    /// Encodes an MHC allele pseudosequence into the flattened BLOSUM62 representation (37×21 =
    /// <see cref="AlleleFlatLength"/> values, position-major). Source: <c>allele_encoding.py</c>
    /// (<c>allele_representations</c>) + <c>amino_acid.fixed_vectors_encoding</c>.
    /// </summary>
    /// <param name="pseudosequence">The allele's pseudosequence (typically <see cref="PseudosequenceLength"/> residues).</param>
    /// <returns>A vector of length <c>pseudosequence.Length × 21</c> of BLOSUM62 values.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="pseudosequence"/> is null.</exception>
    public static double[] EncodePseudosequence(string pseudosequence)
    {
        ArgumentNullException.ThrowIfNull(pseudosequence);
        var result = new double[pseudosequence.Length * EncodingWidth];
        int w = 0;
        foreach (char residue in pseudosequence)
        {
            sbyte[] row = Blosum62[IndexOfResidue(residue)];
            for (int c = 0; c < EncodingWidth; c++)
            {
                result[w++] = row[c];
            }
        }
        return result;
    }

    #endregion

    #region Network forward pass + weight pack

    /// <summary>
    /// A single MHCflurry feed-forward affinity network: an ordered list of dense layers operating on the
    /// 1722-wide peptide‖allele input. All hidden layers use <c>tanh</c>; the final (single-unit) layer uses
    /// <c>sigmoid</c>; the [0,1] output is converted to IC50(nM). Supports both MHCflurry topologies
    /// (<c>feedforward</c> and <c>with-skip-connections</c>).
    /// </summary>
    public sealed class Network
    {
        private readonly DenseLayer[] _layers;
        private readonly bool _withSkipConnections;

        internal Network(DenseLayer[] layers, bool withSkipConnections)
        {
            _layers = layers;
            _withSkipConnections = withSkipConnections;
        }

        /// <summary>Number of dense layers (hidden layers + the single output layer).</summary>
        public int LayerCount => _layers.Length;

        /// <summary>
        /// Runs the forward pass on a pre-concatenated <see cref="NetworkInputLength"/>-wide input and returns
        /// the raw network output in [0,1] (before the IC50 transform). The hidden layers apply <c>tanh</c>; the
        /// output layer applies <c>sigmoid</c>. For <c>with-skip-connections</c> the i-th hidden layer (i ≥ 1)
        /// receives <c>concat(input[i-1], input[i])</c> where the tracked inputs are the merged input followed by
        /// each hidden activation (MHCflurry <c>make_network</c>: <c>densenet_layers[-2:]</c>); the output layer
        /// has no skip. Source: <c>class1_neural_network.make_network</c>.
        /// </summary>
        /// <param name="input">The concatenated peptide‖allele input (length <see cref="NetworkInputLength"/>).</param>
        /// <returns>The raw sigmoid output in [0,1].</returns>
        public double ForwardRaw(ReadOnlySpan<double> input)
        {
            if (input.Length != NetworkInputLength)
            {
                throw new ArgumentException(
                    $"Network input length must be {NetworkInputLength}.", nameof(input));
            }

            int hiddenCount = _layers.Length - 1;

            if (!_withSkipConnections)
            {
                double[] current = input.ToArray();
                for (int li = 0; li < hiddenCount; li++)
                {
                    current = _layers[li].Apply(current, Activation.Tanh);
                }
                double[] outVec = _layers[hiddenCount].Apply(current, Activation.Sigmoid);
                return outVec[0];
            }

            // with-skip-connections: track the inputs that were fed to each dense layer.
            // appended[0] = merged input; appended[k] = activation of hidden layer k-1.
            var appended = new List<double[]>(hiddenCount + 1);
            double[] merged = input.ToArray();
            double[] cur = merged;
            for (int li = 0; li < hiddenCount; li++)
            {
                appended.Add(cur);
                double[] layerInput = appended.Count > 1
                    ? Concat(appended[^2], appended[^1])
                    : appended[0];
                cur = _layers[li].Apply(layerInput, Activation.Tanh);
            }
            double[] output = _layers[hiddenCount].Apply(cur, Activation.Sigmoid);
            return output[0];
        }

        /// <summary>
        /// Predicts the IC50 (nM) for the given peptide / allele-pseudosequence: encodes, runs the forward pass,
        /// then applies <c>IC50 = 50000^(1 − output)</c> (<see cref="ToIc50"/>).
        /// </summary>
        /// <param name="peptide">Peptide sequence.</param>
        /// <param name="pseudosequence">Allele pseudosequence.</param>
        /// <returns>The predicted IC50 in nM.</returns>
        public double PredictIc50(string peptide, string pseudosequence)
        {
            double[] pep = EncodePeptide(peptide);
            double[] allele = EncodePseudosequence(pseudosequence);
            double[] merged = Concat(pep, allele);
            return ToIc50(ForwardRaw(merged));
        }

        private static double[] Concat(double[] a, double[] b)
        {
            var result = new double[a.Length + b.Length];
            Array.Copy(a, 0, result, 0, a.Length);
            Array.Copy(b, 0, result, a.Length, b.Length);
            return result;
        }
    }

    internal enum Activation
    {
        Tanh,
        Sigmoid,
    }

    internal sealed class DenseLayer
    {
        // Kernel stored row-major as [inDim][outDim]; bias length outDim. Matches Keras Dense:
        // output[j] = sum_i input[i] * kernel[i, j] + bias[j].
        private readonly float[] _kernel;
        private readonly float[] _bias;
        private readonly int _inDim;
        private readonly int _outDim;

        internal DenseLayer(float[] kernel, float[] bias, int inDim, int outDim)
        {
            _kernel = kernel;
            _bias = bias;
            _inDim = inDim;
            _outDim = outDim;
        }

        internal int InDim => _inDim;

        internal double[] Apply(double[] input, Activation activation)
        {
            if (input.Length != _inDim)
            {
                throw new ArgumentException(
                    $"Dense layer expected input length {_inDim} but got {input.Length}.", nameof(input));
            }
            var output = new double[_outDim];
            for (int j = 0; j < _outDim; j++)
            {
                output[j] = _bias[j];
            }
            for (int i = 0; i < _inDim; i++)
            {
                double xi = input[i];
                if (xi == 0.0)
                {
                    continue;
                }
                int rowBase = i * _outDim;
                for (int j = 0; j < _outDim; j++)
                {
                    output[j] += xi * _kernel[rowBase + j];
                }
            }
            for (int j = 0; j < _outDim; j++)
            {
                output[j] = activation == Activation.Tanh
                    ? Math.Tanh(output[j])
                    : Sigmoid(output[j]);
            }
            return output;
        }

        private static double Sigmoid(double x) => 1.0 / (1.0 + Math.Exp(-x));
    }

    /// <summary>
    /// Converts a network output in [0,1] to an IC50 in nM via the MHCflurry regression-target transform
    /// <c>IC50 = 50000^(1 − x)</c>. Source: <c>regression_target.to_ic50</c>.
    /// </summary>
    /// <param name="output">The raw sigmoid output in [0,1].</param>
    /// <returns>The IC50 in nM.</returns>
    public static double ToIc50(double output) => Math.Pow(MaxIc50Nm, 1.0 - output);

    /// <summary>
    /// Loads an MHCflurry weight pack (an ensemble of pan-allele affinity networks) from a binary stream.
    /// <para>
    /// <b>Pack format</b> (little-endian): magic <c>"MHCF"</c> (4 bytes); int32 version (=1); int32 networkCount;
    /// then, per network: int32 topology (0 = feedforward, 1 = with-skip-connections); int32 denseLayerCount
    /// (hidden layers + output); then, per dense layer: int32 inDim, int32 outDim, then <c>inDim×outDim</c>
    /// float32 kernel values (row-major), then <c>outDim</c> float32 bias values. The first network's first
    /// hidden layer must have <c>inDim == </c><see cref="NetworkInputLength"/>. This format is produced from the
    /// MHCflurry <c>models_class1_pan</c> <c>weights_*.npz</c> dense kernels/biases (the empty embedding array is
    /// dropped because the allele pseudosequence is supplied via <see cref="EncodePseudosequence"/>).
    /// </para>
    /// </summary>
    /// <param name="stream">The weight-pack stream.</param>
    /// <returns>The ensemble of networks (in pack order).</returns>
    /// <exception cref="ArgumentNullException"><paramref name="stream"/> is null.</exception>
    /// <exception cref="InvalidDataException">The stream is not a valid weight pack.</exception>
    public static IReadOnlyList<Network> LoadWeightPack(Stream stream)
    {
        ArgumentNullException.ThrowIfNull(stream);
        using var reader = new BinaryReader(stream, System.Text.Encoding.ASCII, leaveOpen: true);

        Span<byte> magic = stackalloc byte[4];
        if (reader.Read(magic) != 4 || !magic.SequenceEqual(WeightPackMagic))
        {
            throw new InvalidDataException("Not an MHCflurry weight pack (bad magic).");
        }
        int version = ReadInt32(reader);
        if (version != WeightPackVersion)
        {
            throw new InvalidDataException($"Unsupported MHCflurry weight-pack version {version} (expected {WeightPackVersion}).");
        }
        int networkCount = ReadInt32(reader);
        if (networkCount <= 0)
        {
            throw new InvalidDataException("Weight pack contains no networks.");
        }

        var networks = new Network[networkCount];
        for (int n = 0; n < networkCount; n++)
        {
            int topology = ReadInt32(reader);
            if (topology != TopologyFeedforward && topology != TopologyWithSkipConnections)
            {
                throw new InvalidDataException($"Unknown topology code {topology} in weight pack.");
            }
            int layerCount = ReadInt32(reader);
            if (layerCount < 2)
            {
                throw new InvalidDataException("A network must have at least one hidden layer plus an output layer.");
            }
            var layers = new DenseLayer[layerCount];
            for (int li = 0; li < layerCount; li++)
            {
                int inDim = ReadInt32(reader);
                int outDim = ReadInt32(reader);
                if (inDim <= 0 || outDim <= 0)
                {
                    throw new InvalidDataException("Invalid dense-layer dimensions in weight pack.");
                }
                var kernel = new float[inDim * outDim];
                for (int k = 0; k < kernel.Length; k++)
                {
                    kernel[k] = ReadSingle(reader);
                }
                var bias = new float[outDim];
                for (int b = 0; b < outDim; b++)
                {
                    bias[b] = ReadSingle(reader);
                }
                layers[li] = new DenseLayer(kernel, bias, inDim, outDim);
            }
            if (layers[0].InDim != NetworkInputLength)
            {
                throw new InvalidDataException(
                    $"First dense layer input dim {layers[0].InDim} != expected {NetworkInputLength}.");
            }
            networks[n] = new Network(layers, topology == TopologyWithSkipConnections);
        }
        return networks;
    }

    private static int ReadInt32(BinaryReader reader)
    {
        Span<byte> buf = stackalloc byte[4];
        if (reader.Read(buf) != 4)
        {
            throw new InvalidDataException("Unexpected end of weight pack.");
        }
        return BinaryPrimitives.ReadInt32LittleEndian(buf);
    }

    private static float ReadSingle(BinaryReader reader)
    {
        Span<byte> buf = stackalloc byte[4];
        if (reader.Read(buf) != 4)
        {
            throw new InvalidDataException("Unexpected end of weight pack.");
        }
        return BinaryPrimitives.ReadSingleLittleEndian(buf);
    }

    #endregion

    #region Ensemble prediction + predict→classify chain

    /// <summary>
    /// Predicts the ensemble IC50 (nM) for a peptide / allele pair using an MHCflurry network ensemble: each
    /// network predicts an IC50, and the ensemble value is the <b>geometric mean</b> of the per-network IC50s
    /// (<c>exp(mean(log(ic50)))</c>), exactly as MHCflurry combines models with the default
    /// <c>centrality_measure="mean"</c> in log space. Source:
    /// <c>class1_affinity_predictor.predict_to_dataframe</c> + <c>ensemble_centrality.py</c>.
    /// </summary>
    /// <param name="networks">The network ensemble (from <see cref="LoadWeightPack(Stream)"/>).</param>
    /// <param name="peptide">Peptide sequence.</param>
    /// <param name="allele">MHCflurry allele name (looked up in the bundled pseudosequence table).</param>
    /// <returns>The ensemble IC50 in nM.</returns>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">The ensemble is empty.</exception>
    /// <exception cref="KeyNotFoundException">The allele is not in the bundled pseudosequence table.</exception>
    public static double PredictIc50(IReadOnlyList<Network> networks, string peptide, string allele)
    {
        ArgumentNullException.ThrowIfNull(networks);
        ArgumentNullException.ThrowIfNull(peptide);
        ArgumentNullException.ThrowIfNull(allele);
        string pseudosequence = GetPseudosequence(allele);
        return PredictIc50WithPseudosequence(networks, peptide, pseudosequence);
    }

    /// <summary>
    /// As <see cref="PredictIc50(IReadOnlyList{Network}, string, string)"/> but with a caller-supplied
    /// pseudosequence (bypassing the bundled allele table), for alleles not in the table or for testing.
    /// </summary>
    /// <param name="networks">The network ensemble.</param>
    /// <param name="peptide">Peptide sequence.</param>
    /// <param name="pseudosequence">Allele pseudosequence.</param>
    /// <returns>The ensemble IC50 in nM (geometric mean of the per-network IC50s).</returns>
    public static double PredictIc50WithPseudosequence(
        IReadOnlyList<Network> networks, string peptide, string pseudosequence)
    {
        ArgumentNullException.ThrowIfNull(networks);
        ArgumentNullException.ThrowIfNull(peptide);
        ArgumentNullException.ThrowIfNull(pseudosequence);
        if (networks.Count == 0)
        {
            throw new ArgumentException("The network ensemble is empty.", nameof(networks));
        }

        double[] pep = EncodePeptide(peptide);
        double[] allele = EncodePseudosequence(pseudosequence);
        var merged = new double[pep.Length + allele.Length];
        Array.Copy(pep, 0, merged, 0, pep.Length);
        Array.Copy(allele, 0, merged, pep.Length, allele.Length);

        double logSum = 0.0;
        for (int i = 0; i < networks.Count; i++)
        {
            double ic50 = ToIc50(networks[i].ForwardRaw(merged));
            logSum += Math.Log(ic50);
        }
        return Math.Exp(logSum / networks.Count);
    }

    /// <summary>
    /// End-to-end MHCflurry prediction → classification: predicts the ensemble IC50 with
    /// <see cref="PredictIc50(IReadOnlyList{Network}, string, string)"/> and classifies it into a
    /// <see cref="OncologyAnalyzer.BindingStrength"/> via the existing
    /// <see cref="OncologyAnalyzer.ClassifyBindingAffinity(double)"/> (strong &lt; 50 nM, weak &lt; 500 nM).
    /// The classification cutoffs and defaults are unchanged.
    /// </summary>
    /// <param name="networks">The MHCflurry network ensemble.</param>
    /// <param name="peptide">Peptide sequence.</param>
    /// <param name="allele">MHCflurry allele name.</param>
    /// <returns>The predicted IC50 (nM) and its binding-strength classification.</returns>
    public static (double Ic50Nm, OncologyAnalyzer.BindingStrength Strength) PredictAndClassify(
        IReadOnlyList<Network> networks, string peptide, string allele)
    {
        double ic50 = PredictIc50(networks, peptide, allele);
        return (ic50, OncologyAnalyzer.ClassifyBindingAffinity(ic50));
    }

    #endregion
}
