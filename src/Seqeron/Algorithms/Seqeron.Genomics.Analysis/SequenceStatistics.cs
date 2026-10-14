namespace Seqeron.Genomics.Analysis;

/// <summary>
/// Provides comprehensive statistical analysis of biological sequences.
/// Includes composition, thermodynamics, and physical properties.
/// </summary>
public static class SequenceStatistics
{
    #region Basic Composition

    /// <summary>
    /// DNA/RNA nucleotide composition statistics.
    /// </summary>
    public readonly record struct NucleotideComposition(
        int Length,
        int CountA,
        int CountT,
        int CountG,
        int CountC,
        int CountU,
        int CountN,
        int CountOther,
        double GcContent,
        double AtContent,
        double GcSkew,
        double AtSkew);

    /// <summary>
    /// Amino acid composition statistics.
    /// </summary>
    public readonly record struct AminoAcidComposition(
        int Length,
        IReadOnlyDictionary<char, int> Counts,
        double MolecularWeight,
        double IsoelectricPoint,
        double Hydrophobicity,
        double ChargedResidueRatio,
        double AromaticResidueRatio);

    /// <summary>
    /// Calculates nucleotide composition of a DNA/RNA sequence: per-symbol counts, GC/AT content
    /// and GC/AT skew.
    /// </summary>
    /// <remarks>
    /// <para>Counting is case-insensitive. A, T, G, C, U and N are counted individually; every other
    /// character (IUPAC ambiguity codes R/Y/S/W/K/M/B/D/H/V, gaps, digits, ...) is counted as
    /// <c>CountOther</c>, so the counts partition <c>Length</c>.</para>
    /// <para><c>GcContent</c> = (G+C)/(A+T+G+C+U), delegated to the canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char})"/> (Biopython
    /// <c>gc_fraction</c> "remove" mode restricted to the unambiguous alphabet; S/W are not counted —
    /// use <see cref="SequenceExtensions.CalculateGcFraction(string, GcAmbiguityMode)"/> for exact
    /// Biopython parity on ambiguity codes). <c>AtContent</c> = (A+T+U)/(A+T+G+C+U).</para>
    /// <para><c>GcSkew</c> = (G−C)/(G+C) and <c>AtSkew</c> = (A−T)/(A+T) (Lobry 1996; Biopython
    /// <c>GC_skew</c>), delegated to the canonical <see cref="GcSkewCalculator.CalculateGcSkew(string)"/>
    /// and <see cref="GcSkewCalculator.CalculateAtSkew(string)"/>; each is 0 when its denominator is 0.
    /// AT skew is the DNA definition: U is not paired with A (an all-RNA sequence therefore has
    /// AtSkew = +1 whenever it contains A).</para>
    /// <para>Null or empty input returns an all-zero composition.</para>
    /// </remarks>
    public static NucleotideComposition CalculateNucleotideComposition(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new NucleotideComposition(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);
        }

        int a = 0, t = 0, g = 0, c = 0, u = 0, n = 0, other = 0;

        foreach (char ch in sequence.ToUpperInvariant())
        {
            switch (ch)
            {
                case 'A': a++; break;
                case 'T': t++; break;
                case 'G': g++; break;
                case 'C': c++; break;
                case 'U': u++; break;
                case 'N': n++; break;
                default: other++; break;
            }
        }

        int total = a + t + g + c + u;
        int at = a + t + u;

        // Canonical implementations (no re-implementation): GC fraction lives in
        // Core/SequenceExtensions, GC/AT skew in Analysis/GcSkewCalculator.
        double gcContent = sequence.AsSpan().CalculateGcFraction();
        double atContent = total > 0 ? (double)at / total : 0;
        double gcSkew = GcSkewCalculator.CalculateGcSkew(sequence);
        double atSkew = GcSkewCalculator.CalculateAtSkew(sequence);

        return new NucleotideComposition(
            Length: sequence.Length,
            CountA: a,
            CountT: t,
            CountG: g,
            CountC: c,
            CountU: u,
            CountN: n,
            CountOther: other,
            GcContent: gcContent,
            AtContent: atContent,
            GcSkew: gcSkew,
            AtSkew: atSkew);
    }

    /// <summary>
    /// Calculates amino acid composition of a protein sequence.
    /// </summary>
    /// <remarks>
    /// <para><c>Counts</c> holds the case-insensitive count of every ASCII letter (the 20 standard residues
    /// plus any ambiguity/extended codes such as B, Z, X, U, O, J), folded ASCII-only as in
    /// <see cref="ProteinPhysicochemistry"/>; non-ASCII letters (e.g. U+017F 'ſ', U+0131 'ı') are not
    /// residue codes (EMBOSS pepstats skips them). <c>Length</c> is the number of
    /// letters; non-letter symbols (stop <c>*</c>, gap <c>-</c>, digits, whitespace) are not residues
    /// and are excluded from <c>Counts</c>, <c>Length</c> and the ratio denominators. This differs from
    /// Biopython <c>ProteinAnalysis</c>, whose denominator is <c>len(seq)</c> including such symbols
    /// (opt in with <see cref="CalculateAminoAcidComposition(string, AminoAcidCompositionDenominator)"/>,
    /// <see cref="AminoAcidCompositionDenominator.SequenceLength"/>); for letter-only input the two agree
    /// exactly. The residue denominator is the EMBOSS pepstats "Mole%" convention.</para>
    /// <para><c>AromaticResidueRatio</c> = (F+W+Y)/Length — the aromaticity of Lobry &amp; Gautier
    /// (1994) as implemented by Biopython <c>ProteinAnalysis.aromaticity()</c> (EMBOSS pepstats'
    /// "Aromatic" class additionally includes H). <c>ChargedResidueRatio</c> = (D+E+H+K+R)/Length —
    /// the EMBOSS pepstats "Charged" class (B+D+E+H+K+R+Z) restricted to unambiguous residues.</para>
    /// <para><c>MolecularWeight</c>, <c>IsoelectricPoint</c> and <c>Hydrophobicity</c> delegate to
    /// <see cref="CalculateMolecularWeight(string)"/>, <see cref="CalculateIsoelectricPoint"/> and
    /// <see cref="CalculateHydrophobicity"/>. Null/empty input returns Length 0, empty counts,
    /// ratios 0 and the neutral-pH pI default of those methods.</para>
    /// </remarks>
    public static AminoAcidComposition CalculateAminoAcidComposition(string sequence) =>
        CalculateAminoAcidComposition(sequence, AminoAcidCompositionDenominator.Residues);

    /// <summary>
    /// Denominator of the amino-acid composition ratios
    /// (<see cref="CalculateAminoAcidComposition(string, AminoAcidCompositionDenominator)"/>,
    /// <see cref="CalculateAminoAcidFractions"/>).
    /// </summary>
    public enum AminoAcidCompositionDenominator
    {
        /// <summary>
        /// Number of residues (ASCII letters) — EMBOSS 6.6.0 pepstats "Mole%" convention: stop <c>*</c>,
        /// gaps, digits and whitespace are not residues (pepstats "MKW*" → Residues = 3, Mole% 33.333).
        /// Default.
        /// </summary>
        Residues = 0,

        /// <summary>
        /// Full input length <c>sequence.Length</c>, every character included — Biopython 1.88
        /// <c>ProteinAnalysis</c> convention (<c>self.length = len(self.sequence)</c>, used by
        /// <c>amino_acids_percent</c> and <c>aromaticity()</c>): <c>ProteinAnalysis("MKW*")</c> →
        /// percent 25.0 each, aromaticity 0.25. For letter-only input both denominators agree.
        /// </summary>
        SequenceLength = 1,
    }

    /// <summary>
    /// Amino-acid composition with a selectable ratio denominator. With
    /// <see cref="AminoAcidCompositionDenominator.Residues"/> this is exactly
    /// <see cref="CalculateAminoAcidComposition(string)"/>; with
    /// <see cref="AminoAcidCompositionDenominator.SequenceLength"/> <c>ChargedResidueRatio</c> and
    /// <c>AromaticResidueRatio</c> divide by <c>sequence.Length</c> (Biopython 1.88
    /// <c>ProteinAnalysis.aromaticity()</c>, e.g. "mkwfy**" → 3/7). <c>Length</c>, <c>Counts</c>,
    /// molecular weight, pI and GRAVY are unaffected by the denominator.
    /// </summary>
    /// <param name="sequence">Protein sequence (case-insensitive).</param>
    /// <param name="denominator">Ratio denominator (default convention: residues).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="denominator"/> is not a defined value.</exception>
    public static AminoAcidComposition CalculateAminoAcidComposition(
        string sequence,
        AminoAcidCompositionDenominator denominator)
    {
        ValidateCompositionDenominator(denominator);
        if (string.IsNullOrEmpty(sequence))
        {
            return new AminoAcidComposition(0, new Dictionary<char, int>(), 0, NeutralPhDefault, 0, 0, 0);
        }

        var counts = CountAsciiResidues(sequence);
        int residues = counts.Values.Sum();
        int ratioDenominator = denominator == AminoAcidCompositionDenominator.SequenceLength ? sequence.Length : residues;
        double mw = CalculateMolecularWeight(sequence);
        double pi = CalculateIsoelectricPoint(sequence);
        double hydro = CalculateHydrophobicity(sequence);

        // Charged residues: D, E (negative), K, R, H (positive)
        int charged = counts.GetValueOrDefault('D') + counts.GetValueOrDefault('E') +
                     counts.GetValueOrDefault('K') + counts.GetValueOrDefault('R') +
                     counts.GetValueOrDefault('H');
        double chargedRatio = ratioDenominator > 0 ? (double)charged / ratioDenominator : 0;

        // Aromatic residues: F, Y, W
        int aromatic = counts.GetValueOrDefault('F') + counts.GetValueOrDefault('Y') +
                      counts.GetValueOrDefault('W');
        double aromaticRatio = ratioDenominator > 0 ? (double)aromatic / ratioDenominator : 0;

        return new AminoAcidComposition(
            Length: residues,
            Counts: counts,
            MolecularWeight: mw,
            IsoelectricPoint: pi,
            Hydrophobicity: hydro,
            ChargedResidueRatio: chargedRatio,
            AromaticResidueRatio: aromaticRatio);
    }

    /// <summary>
    /// Per-residue composition fractions count(X) / denominator for every residue code (ASCII letter,
    /// case-insensitive, ASCII-only fold) that occurs in <paramref name="sequence"/>.
    /// </summary>
    /// <remarks>
    /// <see cref="AminoAcidCompositionDenominator.Residues"/> (default) reproduces EMBOSS 6.6.0 pepstats
    /// "Mole%"/100 (denominator = residues; "MKW*" → 1/3 each); <see cref="AminoAcidCompositionDenominator.SequenceLength"/>
    /// reproduces Biopython 1.88 <c>ProteinAnalysis.amino_acids_percent</c>/100 (denominator = len(seq);
    /// "MKW*" → 0.25 each, "MKW-A" → 0.2 each). Absent residues are not keys (fraction 0). Unlike
    /// Biopython, which reports only the 20 standard residues, ambiguity/extended codes (B, X, Z, U, O, J)
    /// are keys too — as the per-letter rows of pepstats. Null/empty input returns an empty map.
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="denominator"/> is not a defined value.</exception>
    public static IReadOnlyDictionary<char, double> CalculateAminoAcidFractions(
        string sequence,
        AminoAcidCompositionDenominator denominator = AminoAcidCompositionDenominator.Residues)
    {
        ValidateCompositionDenominator(denominator);
        var fractions = new Dictionary<char, double>();
        if (string.IsNullOrEmpty(sequence))
            return fractions;

        var counts = CountAsciiResidues(sequence);
        int total = denominator == AminoAcidCompositionDenominator.SequenceLength
            ? sequence.Length
            : counts.Values.Sum();
        foreach (var (residue, count) in counts)
            fractions[residue] = (double)count / total;

        return fractions;
    }

    // Residue codes are ASCII letters, folded ASCII-only (as ProteinPhysicochemistry):
    // non-ASCII letters such as U+017F 'ſ' are not residues (EMBOSS pepstats skips them).
    private static Dictionary<char, int> CountAsciiResidues(string sequence)
    {
        var counts = new Dictionary<char, int>();
        foreach (char raw in sequence)
        {
            if (char.IsAsciiLetter(raw))
            {
                char ch = ToUpperAscii(raw);
                counts[ch] = counts.GetValueOrDefault(ch) + 1;
            }
        }

        return counts;
    }

    private static void ValidateCompositionDenominator(AminoAcidCompositionDenominator denominator)
    {
        if (!Enum.IsDefined(denominator))
            throw new ArgumentOutOfRangeException(nameof(denominator), denominator, "Undefined composition denominator.");
    }

    #endregion

    #region Molecular Weight

    // Average isotopic mass of one water molecule (Da): Biopython SeqUtils water = 18.0153
    // (Expasy FindMod H2O = 18.01524). Single source: ProteinPhysicochemistry.AverageWaterMass.
    private const double AverageWaterMass = ProteinPhysicochemistry.AverageWaterMass;

    // Average molecular masses of DNA mononucleotides (5'-monophosphate, Da).
    // Source: Biopython Bio/Data/IUPACData.py `unambiguous_dna_weights`.
    private static readonly Dictionary<char, double> DnaNucleotideWeights = new()
    {
        { 'A', 331.2218 }, { 'C', 307.1971 }, { 'G', 347.2212 }, { 'T', 322.2085 }
    };

    // Average molecular masses of RNA mononucleotides (5'-monophosphate, Da).
    // Source: Biopython Bio/Data/IUPACData.py `unambiguous_rna_weights`.
    private static readonly Dictionary<char, double> RnaNucleotideWeights = new()
    {
        { 'A', 347.2212 }, { 'C', 323.1965 }, { 'G', 363.2206 }, { 'U', 324.1813 }
    };

    // Monoisotopic mass of one water molecule (Da): Biopython 1.88 Bio.SeqUtils.molecular_weight
    // (monoisotopic=True) water = 18.010565.
    private const double MonoisotopicWaterMass = 18.010565;

    // Monoisotopic masses of the free amino acids (20 standard + U Sec + O Pyl, Da).
    // Source: Biopython 1.88 Bio/Data/IUPACData.py `monoisotopic_protein_weights`.
    // (Proposed home: Infrastructure/ProteinPhysicochemistry next to AverageAminoAcidMasses —
    // B03 cross-batch request R26 to B02.)
    private static readonly Dictionary<char, double> MonoisotopicAminoAcidWeights = new()
    {
        { 'A', 89.047678 },  { 'C', 121.019749 }, { 'D', 133.037508 }, { 'E', 147.053158 },
        { 'F', 165.078979 }, { 'G', 75.032028 },  { 'H', 155.069477 }, { 'I', 131.094629 },
        { 'K', 146.105528 }, { 'L', 131.094629 }, { 'M', 149.051049 }, { 'N', 132.053492 },
        { 'O', 255.158292 }, { 'P', 115.063329 }, { 'Q', 146.069142 }, { 'R', 174.111676 },
        { 'S', 105.042593 }, { 'T', 119.058243 }, { 'U', 168.964203 }, { 'V', 117.078979 },
        { 'W', 204.089878 }, { 'Y', 181.073893 }
    };

    // Monoisotopic masses of DNA mononucleotides (5'-monophosphate, Da).
    // Source: Biopython 1.88 Bio/Data/IUPACData.py `monoisotopic_unambiguous_dna_weights`.
    private static readonly Dictionary<char, double> MonoisotopicDnaNucleotideWeights = new()
    {
        { 'A', 331.06817 }, { 'C', 307.056936 }, { 'G', 347.063084 }, { 'T', 322.056602 }
    };

    // Monoisotopic masses of RNA mononucleotides (5'-monophosphate, Da).
    // Source: Biopython 1.88 Bio/Data/IUPACData.py `monoisotopic_unambiguous_rna_weights`.
    private static readonly Dictionary<char, double> MonoisotopicRnaNucleotideWeights = new()
    {
        { 'A', 347.063084 }, { 'C', 323.051851 }, { 'G', 363.057999 }, { 'U', 324.035867 }
    };

    /// <summary>
    /// Options of the molecular-weight calculations
    /// (<see cref="CalculateMolecularWeight(string, MolecularWeightOptions)"/>,
    /// <see cref="CalculateNucleotideMolecularWeight(string, bool, bool, MolecularWeightOptions)"/>),
    /// mirroring the keyword arguments of Biopython 1.88 <c>Bio.SeqUtils.molecular_weight</c>.
    /// The defaults reproduce the existing overloads (average masses, linear, lenient).
    /// </summary>
    public sealed record MolecularWeightOptions
    {
        /// <summary>Default options: average masses, linear molecule, unknown symbols skipped.</summary>
        public static MolecularWeightOptions Default { get; } = new();

        /// <summary>
        /// True to use monoisotopic monomer masses and water = 18.010565 Da (Biopython
        /// <c>monoisotopic=True</c>); false (default) for average masses and water = 18.0153 Da.
        /// </summary>
        public bool Monoisotopic { get; init; }

        /// <summary>
        /// True for a circular molecule: one further water is lost per strand for the ring-closing
        /// bond (Biopython <c>circular=True</c>; for a protein, a head-to-tail cyclic peptide).
        /// </summary>
        public bool Circular { get; init; }

        /// <summary>
        /// True to reject any symbol that has no mass in the selected table (ambiguity codes
        /// B/Z/J/X, '*', gaps, N, U in DNA, T in RNA, …) with <see cref="ArgumentException"/>, as
        /// Biopython raises <c>ValueError</c>; whitespace is ignored (Biopython strips it).
        /// False (default) skips such symbols (no mass, no bond).
        /// </summary>
        public bool Strict { get; init; }
    }

    // Shared accumulation of the MW of one strand: sum of monomer masses − bonds × water, where a
    // linear strand of n monomers has n − 1 bonds and a circular one n bonds (Biopython
    // molecular_weight). Case folding is ASCII-only. In strict mode a non-whitespace symbol without
    // a mass throws; otherwise it is skipped. Returns (sum of masses, monomer count).
    private static (double Sum, int Monomers) SumMonomerMasses(
        string sequence, IReadOnlyDictionary<char, double> table, bool strict, string paramName, string kind)
    {
        double sum = 0;
        int monomers = 0;
        foreach (char raw in sequence)
        {
            char ch = ToUpperAscii(raw);
            if (table.TryGetValue(ch, out double mass))
            {
                sum += mass;
                monomers++;
            }
            else if (strict && !char.IsWhiteSpace(raw))
            {
                throw new ArgumentException(
                    $"'{raw}' is not a valid unambiguous letter for {kind} (strict mode).", paramName);
            }
        }
        return (sum, monomers);
    }

    /// <summary>
    /// Calculates the average-isotopic molecular weight of a protein sequence (Da).
    /// </summary>
    /// <remarks>
    /// Implements the Expasy Compute pI/Mw definition: the sum of the average isotopic
    /// masses of the amino acids plus the average isotopic mass of one water molecule.
    /// Equivalently (Biopython): sum(free amino-acid masses) − (n − 1) × water, removing
    /// one water per peptide bond. Recognized alphabet: the 20 standard amino acids plus
    /// selenocysteine (U) and pyrrolysine (O), as in Biopython <c>protein_weights</c>.
    /// Unknown/ambiguous symbols (B, Z, X, J, '*', gaps, …) are skipped (contribute no mass and
    /// no bond) — Biopython instead raises <c>ValueError</c>. Linear chain, average masses only.
    /// Delegates to the canonical <see cref="ProteinPhysicochemistry.MolecularWeight"/>; case folding
    /// is ASCII-only, so a non-ASCII letter such as 'ſ' is skipped like any unknown symbol (EMBOSS
    /// pepstats; Python's Unicode <c>str.upper()</c> in Biopython would read it as S).
    /// </remarks>
    /// <param name="proteinSequence">Protein sequence (case-insensitive, one-letter codes).</param>
    /// <returns>Molecular weight in daltons; 0 for null/empty input.</returns>
    public static double CalculateMolecularWeight(string proteinSequence) =>
        ProteinPhysicochemistry.MolecularWeight(proteinSequence);

    /// <summary>
    /// Calculates the molecular weight of a protein (Da) with Biopython 1.88
    /// <c>molecular_weight(seq, "protein", circular, monoisotopic)</c> options and an optional strict mode.
    /// </summary>
    /// <remarks>
    /// Σ free amino-acid masses − bonds × water, bonds = n − 1 (linear) or n (circular, head-to-tail
    /// cyclic peptide: one more water lost). Masses: average (<see cref="ProteinPhysicochemistry.AverageAminoAcidMasses"/>,
    /// water 18.0153) or monoisotopic (Biopython <c>IUPACData.monoisotopic_protein_weights</c>,
    /// water 18.010565). With <see cref="MolecularWeightOptions.Default"/> this is exactly
    /// <see cref="CalculateMolecularWeight(string)"/>. Alphabet: 20 standard residues + U + O.
    /// </remarks>
    /// <param name="proteinSequence">Protein sequence (case-insensitive, one-letter codes).</param>
    /// <param name="options">Mass set, topology and strictness; null = <see cref="MolecularWeightOptions.Default"/>.</param>
    /// <returns>Molecular weight in daltons; 0 for null/empty input or when no residue has a mass.</returns>
    /// <exception cref="ArgumentException">Strict mode and a symbol without a mass (e.g. B, Z, X, '*').</exception>
    public static double CalculateMolecularWeight(string proteinSequence, MolecularWeightOptions? options)
    {
        options ??= MolecularWeightOptions.Default;
        if (string.IsNullOrEmpty(proteinSequence))
            return 0;
        if (!options.Monoisotopic && !options.Circular && !options.Strict)
            return ProteinPhysicochemistry.MolecularWeight(proteinSequence);

        IReadOnlyDictionary<char, double> table = options.Monoisotopic
            ? MonoisotopicAminoAcidWeights
            : ProteinPhysicochemistry.AverageAminoAcidMasses;
        double water = options.Monoisotopic ? MonoisotopicWaterMass : AverageWaterMass;

        var (sum, residues) = SumMonomerMasses(proteinSequence, table, options.Strict, nameof(proteinSequence), "protein");
        if (residues == 0)
            return 0;

        int bonds = options.Circular ? residues : residues - 1;
        return sum - bonds * water;
    }

    /// <summary>
    /// Calculates the average-isotopic molecular weight of a DNA or RNA sequence (Da).
    /// </summary>
    /// <remarks>
    /// Uses average monophosphate (5'-phosphate) mononucleotide masses and removes one
    /// water per phosphodiester bond: sum(monophosphate masses) − (n − 1) × water
    /// (Biopython Bio.SeqUtils.molecular_weight). Single-stranded, linear molecule; for the
    /// double-stranded and/or circular molecule use
    /// <see cref="CalculateNucleotideMolecularWeight(string, bool, bool, bool)"/>.
    /// Unknown symbols are skipped (Biopython raises <c>ValueError</c> instead).
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence (case-insensitive).</param>
    /// <param name="isDna">True for DNA (uses A/C/G/T table); false for RNA (A/C/G/U table).</param>
    /// <returns>Molecular weight in daltons; 0 for null/empty input.</returns>
    public static double CalculateNucleotideMolecularWeight(string sequence, bool isDna = true)
        => CalculateNucleotideMolecularWeight(sequence, isDna, doubleStranded: false, circular: false);

    /// <summary>
    /// Calculates the average-isotopic molecular weight of a DNA or RNA molecule (Da),
    /// optionally double-stranded and/or circular.
    /// </summary>
    /// <remarks>
    /// Realises Biopython <c>Bio.SeqUtils.molecular_weight(seq, seq_type, double_stranded, circular)</c>:
    /// each strand weighs sum(monophosphate masses) − (n − 1) × water; a circular strand loses one
    /// further water (the ring-closing phosphodiester bond); when <paramref name="doubleStranded"/> is
    /// true the complementary strand (canonical <see cref="Seqeron.Genomics.Core.SequenceExtensions.GetComplementBase"/> /
    /// <see cref="Seqeron.Genomics.Core.SequenceExtensions.GetRnaComplementBase"/>) is added with the
    /// same rule. Unknown symbols are skipped on both strands (Biopython raises <c>ValueError</c>).
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence of one strand (case-insensitive).</param>
    /// <param name="isDna">True for DNA (A/C/G/T table); false for RNA (A/C/G/U table).</param>
    /// <param name="doubleStranded">True to add the Watson–Crick complementary strand.</param>
    /// <param name="circular">True for a circular molecule (one extra water lost per strand).</param>
    /// <returns>Molecular weight in daltons; 0 when no recognized nucleotide is present.</returns>
    public static double CalculateNucleotideMolecularWeight(
        string sequence, bool isDna, bool doubleStranded, bool circular = false)
        => CalculateNucleotideMolecularWeight(sequence, isDna, doubleStranded,
            circular ? CircularAverageOptions : MolecularWeightOptions.Default);

    private static readonly MolecularWeightOptions CircularAverageOptions = new() { Circular = true };

    /// <summary>
    /// Calculates the molecular weight of a DNA or RNA molecule (Da) with Biopython 1.88
    /// <c>molecular_weight(seq, seq_type, double_stranded, circular, monoisotopic)</c> options and an
    /// optional strict mode.
    /// </summary>
    /// <remarks>
    /// Each strand weighs Σ monophosphate masses − bonds × water (bonds = n − 1 linear, n circular);
    /// when <paramref name="doubleStranded"/> is true the complementary strand (canonical
    /// <see cref="Seqeron.Genomics.Core.SequenceExtensions.GetComplementBase"/> /
    /// <see cref="Seqeron.Genomics.Core.SequenceExtensions.GetRnaComplementBase"/>) is added by the same
    /// rule. Masses: average (Biopython <c>unambiguous_dna/rna_weights</c>, water 18.0153) or
    /// monoisotopic (<c>monoisotopic_unambiguous_dna/rna_weights</c>, water 18.010565).
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence of one strand (case-insensitive).</param>
    /// <param name="isDna">True for DNA (A/C/G/T table); false for RNA (A/C/G/U table).</param>
    /// <param name="doubleStranded">True to add the Watson–Crick complementary strand.</param>
    /// <param name="options">Mass set, topology and strictness; null = <see cref="MolecularWeightOptions.Default"/>.</param>
    /// <returns>Molecular weight in daltons; 0 when no recognized nucleotide is present.</returns>
    /// <exception cref="ArgumentException">Strict mode and a symbol without a mass (e.g. N; U in DNA; T in RNA).</exception>
    public static double CalculateNucleotideMolecularWeight(
        string sequence, bool isDna, bool doubleStranded, MolecularWeightOptions? options)
    {
        options ??= MolecularWeightOptions.Default;
        if (string.IsNullOrEmpty(sequence))
            return 0;

        Dictionary<char, double> table;
        if (options.Monoisotopic)
            table = isDna ? MonoisotopicDnaNucleotideWeights : MonoisotopicRnaNucleotideWeights;
        else
            table = isDna ? DnaNucleotideWeights : RnaNucleotideWeights;
        double water = options.Monoisotopic ? MonoisotopicWaterMass : AverageWaterMass;

        var (weight, monomers) = SumMonomerMasses(sequence, table, options.Strict, nameof(sequence), isDna ? "DNA" : "RNA");
        if (monomers == 0)
            return 0;

        // One water is lost per phosphodiester bond: n monomers form (n − 1) bonds in a linear
        // strand and n bonds in a circular one.
        int bondsPerStrand = options.Circular ? monomers : monomers - 1;
        double result = weight - bondsPerStrand * water;

        if (doubleStranded)
        {
            double complementWeight = 0;
            foreach (char raw in sequence)
            {
                char ch = ToUpperAscii(raw);
                if (!table.ContainsKey(ch))
                    continue;
                char complement = isDna
                    ? Seqeron.Genomics.Core.SequenceExtensions.GetComplementBase(ch)
                    : Seqeron.Genomics.Core.SequenceExtensions.GetRnaComplementBase(ch);
                complementWeight += table[complement];
            }
            result += complementWeight - bondsPerStrand * water;
        }

        return result;
    }

    #endregion

    #region Isoelectric Point

    /// <summary>
    /// pK set used by <see cref="CalculateIsoelectricPoint(string, PkaScale)"/> and
    /// <see cref="CalculateNetCharge(string, double, PkaScale)"/>.
    /// </summary>
    public enum PkaScale
    {
        /// <summary>
        /// EMBOSS <c>iep</c> default <c>Epk.dat</c> (EMBOSS 6.6.0): N-terminus 7.5, C-terminus 3.6,
        /// C 8.5, D 3.9, E 4.1, H 6.5, K 10.8, R 12.5, Y 10.1; ambiguity codes B/Z are split into
        /// D/N and E/Q by Dayhoff frequencies exactly as <c>embIepCompC</c> does.
        /// </summary>
        Emboss,

        /// <summary>
        /// Bjellqvist et al. 1993/1994 (ExPASy Compute pI/Mw; Biopython
        /// <c>Bio.SeqUtils.IsoelectricPoint</c>): side chains C 9.0, D 4.05, E 4.45, H 5.98,
        /// K 10.0, R 12.0, Y 10.0; N-terminus 7.5 unless the N-terminal residue is
        /// A 7.59 / M 7.0 / S 6.93 / P 8.36 / T 6.82 / V 7.44 / E 7.7; C-terminus 3.55 unless the
        /// C-terminal residue is D 4.55 / E 4.75.
        /// </summary>
        Bjellqvist
    }

    // The pK tables (EMBOSS 6.6.0 Epk.dat; Bjellqvist 1993/1994 incl. terminal-residue pKs), the
    // EMBOSS B/Z Dayhoff split, the [0, 14] bisection to 1e-9 pH and the 2-dp rounding live in the
    // canonical ProteinPhysicochemistry (review 2026-09, B02 request; FIN-B03 A1).

    // pI returned for empty/null input: pI is undefined for a zero-length protein (a real
    // protein always has both termini); neutral 7.0 is used as a documented input-guard sentinel.
    private const double NeutralPhDefault = 7.0;

    /// <summary>
    /// Calculates the theoretical isoelectric point (pI) of a protein on the EMBOSS <c>iep</c>
    /// pK scale. Equivalent to <see cref="CalculateIsoelectricPoint(string, PkaScale)"/> with
    /// <see cref="PkaScale.Emboss"/>.
    /// </summary>
    /// <param name="proteinSequence">Single-letter amino-acid sequence (case-insensitive).</param>
    /// <returns>The isoelectric point in [0, 14], rounded to two decimal places; 7.0 for null/empty.</returns>
    public static double CalculateIsoelectricPoint(string proteinSequence) =>
        CalculateIsoelectricPoint(proteinSequence, PkaScale.Emboss);

    /// <summary>
    /// Calculates the theoretical isoelectric point (pI) of a protein: the pH at which the
    /// Henderson–Hasselbalch net charge (<see cref="CalculateNetCharge"/>) is zero.
    /// The root is located by bisection over [0, 14] to 1e-9 pH and rounded to two decimals.
    /// </summary>
    /// <param name="proteinSequence">Single-letter amino-acid sequence (case-insensitive). Residues
    /// without an ionizable side chain in the chosen scale are ignored. Null or empty returns the
    /// neutral sentinel 7.0.</param>
    /// <param name="scale">pK set: <see cref="PkaScale.Emboss"/> (EMBOSS iep, default) or
    /// <see cref="PkaScale.Bjellqvist"/> (ExPASy Compute pI/Mw, Biopython).</param>
    /// <returns>The isoelectric point in [0, 14], rounded to two decimal places.</returns>
    /// <remarks>
    /// <para>EMBOSS: reproduces EMBOSS 6.6.0 <c>iep</c> (nucleus/embiep.c, data/Epk.dat) with default
    /// options (both termini charged, no disulphides, no modified lysines). EMBOSS searches pH [1, 14]
    /// and reports "none" when the charge does not change sign there; this method searches [0, 14]
    /// (the EMBOSS window, "none" and the <c>iep</c> qualifiers are available through
    /// <see cref="TryCalculateIsoelectricPoint"/>).</para>
    /// <para>Bjellqvist: reproduces Biopython <c>IsoelectricPoint.pi()</c> / ExPASy Compute pI; Biopython
    /// clamps its bisection to [4.05, 12], so for extremely acidic/basic peptides whose true root lies
    /// outside that window Biopython returns the window edge while this method returns the true root
    /// (<see cref="IsoelectricPointWindow.Biopython"/> reproduces Biopython's window).</para>
    /// </remarks>
    public static double CalculateIsoelectricPoint(string proteinSequence, PkaScale scale) =>
        ProteinPhysicochemistry.IsoelectricPoint(proteinSequence, ToProteinPkaScale(scale));

    /// <summary>
    /// Net charge of a protein at a given pH (Henderson–Hasselbalch): basic groups (N-terminus,
    /// K, R, H) contribute +1/(1+10^(pH−pKa)), acidic groups (C-terminus, D, E, C, Y) contribute
    /// −1/(1+10^(pKa−pH)); each terminus is counted once.
    /// Equivalent to EMBOSS <c>embIepGetCharge</c> and Biopython <c>IsoelectricPoint.charge_at_pH</c>.
    /// </summary>
    /// <param name="proteinSequence">Single-letter amino-acid sequence (case-insensitive, ASCII fold).</param>
    /// <param name="pH">pH at which to evaluate the charge.</param>
    /// <param name="scale">pK set (default EMBOSS).</param>
    /// <returns>Net charge (elementary charges); 0 for null/empty input.</returns>
    public static double CalculateNetCharge(string proteinSequence, double pH, PkaScale scale = PkaScale.Emboss) =>
        ProteinPhysicochemistry.NetCharge(proteinSequence, pH, ToProteinPkaScale(scale));

    /// <summary>
    /// pH search convention used by <see cref="TryCalculateIsoelectricPoint"/>.
    /// </summary>
    public enum IsoelectricPointWindow
    {
        /// <summary>
        /// Library convention (= <see cref="CalculateIsoelectricPoint(string, PkaScale)"/>): bisection
        /// over [0, 14] to 1e-9 pH, result rounded to 2 decimals; a root is always reported (the window
        /// edge when the charge does not change sign); null/empty → 7.0.
        /// </summary>
        Default,

        /// <summary>
        /// EMBOSS 6.6.0 <c>iep</c> (<c>embIepPhConverge</c>): bisection over [1, 14] while the bracket
        /// is wider than 1e-4, returning its lower bound (or the midpoint where the charge is exactly
        /// 0), unrounded (<c>iep</c> prints it with 4 decimals). When the net charge has the same strict
        /// sign at pH 1 and pH 14 there is no pI ("Isoelectric Point = None"): the Try method returns
        /// false. Null/empty → false.
        /// </summary>
        Emboss,

        /// <summary>
        /// Biopython 1.88 <c>IsoelectricPoint.pi(pH=7.775, min_=4.05, max_=12)</c>: bisection that
        /// starts at pH 7.775 inside [4.05, 12] and stops once the bracket is ≤ 1e-4, returning the
        /// last midpoint, unrounded. There is no range extension: a root outside the window is reported
        /// at the window edge (4.0500284… for very acidic, 11.9999677… for very basic peptides).
        /// Null/empty → false (Biopython raises <c>IndexError</c>).
        /// </summary>
        Biopython
    }

    /// <summary>
    /// Options for <see cref="TryCalculateIsoelectricPoint"/> and
    /// <see cref="CalculateNetCharge(string, double, IsoelectricPointOptions)"/>, mirroring the
    /// EMBOSS 6.6.0 <c>iep</c> qualifiers <c>-[no]termini</c>, <c>-lysinemodified</c> and
    /// <c>-disulphides</c> and the EMBOSS / Biopython pI search windows. The defaults reproduce
    /// <see cref="CalculateIsoelectricPoint(string)"/> and <see cref="CalculateNetCharge(string, double, PkaScale)"/>.
    /// </summary>
    /// <remarks>
    /// The charge model is always the canonical <see cref="ProteinPhysicochemistry.NetCharge"/>;
    /// options only adjust the residue counts handed to it (EMBOSS <c>embIepCompC</c>) or drop the
    /// termini. Custom pK tables (<c>-pkdata</c>) and N-/C-terminus multiplicities
    /// (<c>-amino</c>/<c>-carboxyl</c>) need a parameterised charge model in
    /// <c>ProteinPhysicochemistry</c> and are not offered here (review 2026-09 request R27 to B02).
    /// </remarks>
    public sealed record IsoelectricPointOptions
    {
        /// <summary>Default options: EMBOSS pK scale, library window, both termini charged, no modifications.</summary>
        public static IsoelectricPointOptions Default { get; } = new();

        /// <summary>pK set (default <see cref="PkaScale.Emboss"/>).</summary>
        public PkaScale Scale { get; init; } = PkaScale.Emboss;

        /// <summary>pH search convention (default <see cref="IsoelectricPointWindow.Default"/>).</summary>
        public IsoelectricPointWindow Window { get; init; } = IsoelectricPointWindow.Default;

        /// <summary>
        /// True (default) to include the charged N- and C-terminus; false = EMBOSS <c>-notermini</c>
        /// (side chains only).
        /// </summary>
        public bool Termini { get; init; } = true;

        /// <summary>
        /// Number of chemically modified (uncharged) lysines, EMBOSS <c>-lysinemodified N</c>: N lysines
        /// are removed from the charged-group count; if the sequence has fewer than N lysines, all are
        /// removed (<c>embIepCompC</c>). Must be ≥ 0.
        /// </summary>
        public int LysinesModified { get; init; }

        /// <summary>
        /// Number of disulphide bonds, EMBOSS <c>-disulphides N</c>: 2N cysteines are removed from the
        /// charged-group count; if the sequence has fewer than 2N cysteines, all are removed
        /// (<c>embIepCompC</c>). Must be ≥ 0.
        /// </summary>
        public int Disulphides { get; init; }
    }

    // Residue used to neutralise a modified K / bonded C in the string handed to
    // ProteinPhysicochemistry: G has no ionizable side chain in either scale and, like K and C, has
    // no residue-specific Bjellqvist terminal pK (N-term A/M/S/P/T/V/E, C-term D/E), so substituting
    // it changes the side-chain counts only and leaves both terminal pKs untouched.
    private const char NeutralResidue = 'G';

    // EMBOSS embIepPhConverge window and bracket width; Biopython IsoelectricPoint.pi defaults.
    private const double EmbossPiLowPh = 1.0;
    private const double EmbossPiHighPh = 14.0;
    private const double EmbossPiBracket = 0.0001;
    private const double BiopythonPiStartPh = 7.775;
    private const double BiopythonPiLowPh = 4.05;
    private const double BiopythonPiHighPh = 12.0;
    private const double BiopythonPiBracket = 0.0001;

    // Library window (= ProteinPhysicochemistry.IsoelectricPoint): [0, 14], 1e-9 pH, 2 dp.
    private const double DefaultPiLowPh = 0.0;
    private const double DefaultPiHighPh = 14.0;
    private const double DefaultPiPrecision = 1e-9;
    private const int DefaultPiDecimalPlaces = 2;

    /// <summary>
    /// Net charge at a given pH with EMBOSS <c>iep</c> options (termini on/off, modified lysines,
    /// disulphide-bonded cysteines). With <see cref="IsoelectricPointOptions.Default"/> (or any
    /// options that keep the termini and modify nothing) this is exactly
    /// <see cref="CalculateNetCharge(string, double, PkaScale)"/>. <see cref="IsoelectricPointOptions.Window"/> is ignored.
    /// </summary>
    /// <param name="proteinSequence">Single-letter amino-acid sequence (case-insensitive, ASCII fold).</param>
    /// <param name="pH">pH at which to evaluate the charge.</param>
    /// <param name="options">Options; null = <see cref="IsoelectricPointOptions.Default"/>.</param>
    /// <returns>Net charge (elementary charges); 0 for null/empty input.</returns>
    /// <exception cref="ArgumentOutOfRangeException">Negative <c>LysinesModified</c>/<c>Disulphides</c> or an unknown scale.</exception>
    public static double CalculateNetCharge(string proteinSequence, double pH, IsoelectricPointOptions? options)
    {
        options ??= IsoelectricPointOptions.Default;
        ValidateIsoelectricPointOptions(options);
        if (string.IsNullOrEmpty(proteinSequence))
            return 0.0;

        return BuildOptionsCharge(proteinSequence, options)(pH);
    }

    /// <summary>
    /// Isoelectric point with EMBOSS <c>iep</c> / Biopython conventions selected by
    /// <paramref name="options"/>: pK scale, search window (<see cref="IsoelectricPointWindow"/>),
    /// <c>-notermini</c>, <c>-lysinemodified N</c> and <c>-disulphides N</c>.
    /// </summary>
    /// <param name="proteinSequence">Single-letter amino-acid sequence (case-insensitive, ASCII fold).</param>
    /// <param name="options">Options; null = <see cref="IsoelectricPointOptions.Default"/>.</param>
    /// <param name="isoelectricPoint">The pI, or <see cref="double.NaN"/> when there is none.</param>
    /// <returns>
    /// False when the selected convention reports no pI: <see cref="IsoelectricPointWindow.Emboss"/>
    /// when the net charge has the same strict sign at pH 1 and 14 (<c>iep</c> "None"), or a
    /// null/empty sequence under the EMBOSS / Biopython windows. Always true for the default window.
    /// </returns>
    /// <remarks>
    /// <para>Default window + default options returns exactly <see cref="CalculateIsoelectricPoint(string, PkaScale)"/>.</para>
    /// <para>EMBOSS window: <c>embIepPhConverge</c> over the <c>embIepCompC</c>-adjusted counts; the
    /// result is unrounded (compare with <c>iep</c>'s 4-dp output). Biopython window:
    /// <c>IsoelectricPoint.pi()</c> with its defaults, unrounded (Bjellqvist scale reproduces
    /// <c>ProteinAnalysis.isoelectric_point()</c>).</para>
    /// <para>The modifications and <c>-notermini</c> are EMBOSS features; combined with the Bjellqvist
    /// scale they apply the same count adjustment (terminal-residue pKs are unaffected).</para>
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException">Negative <c>LysinesModified</c>/<c>Disulphides</c>, unknown scale or window.</exception>
    public static bool TryCalculateIsoelectricPoint(string proteinSequence, IsoelectricPointOptions? options, out double isoelectricPoint)
    {
        options ??= IsoelectricPointOptions.Default;
        ValidateIsoelectricPointOptions(options);
        isoelectricPoint = double.NaN;

        if (string.IsNullOrEmpty(proteinSequence))
        {
            if (options.Window != IsoelectricPointWindow.Default)
                return false;
            isoelectricPoint = NeutralPhDefault;
            return true;
        }

        bool unmodified = options.Termini && options.LysinesModified == 0 && options.Disulphides == 0;
        if (options.Window == IsoelectricPointWindow.Default && unmodified)
        {
            isoelectricPoint = CalculateIsoelectricPoint(proteinSequence, options.Scale);
            return true;
        }

        Func<double, double> charge = BuildOptionsCharge(proteinSequence, options);
        switch (options.Window)
        {
            case IsoelectricPointWindow.Default:
                isoelectricPoint = BisectDefaultWindow(charge);
                return true;
            case IsoelectricPointWindow.Emboss:
                return TryBisectEmbossWindow(charge, out isoelectricPoint);
            default: // Biopython (validated above)
                isoelectricPoint = BisectBiopythonWindow(charge);
                return true;
        }
    }

    private static void ValidateIsoelectricPointOptions(IsoelectricPointOptions options)
    {
        if (options.LysinesModified < 0)
            throw new ArgumentOutOfRangeException(nameof(options), options.LysinesModified, "LysinesModified must be non-negative.");
        if (options.Disulphides < 0)
            throw new ArgumentOutOfRangeException(nameof(options), options.Disulphides, "Disulphides must be non-negative.");
        if (!Enum.IsDefined(options.Window))
            throw new ArgumentOutOfRangeException(nameof(options), options.Window, "Unknown pI window.");
        _ = ToProteinPkaScale(options.Scale);
    }

    // Charge function composed over the canonical ProteinPhysicochemistry.NetCharge: modified
    // lysines / bonded cysteines are neutralised in the string (embIepCompC count adjustment); with
    // Termini = false the side-chain charge is NetCharge("G"+s+"G") − NetCharge("GG"), whose default
    // termini cancel (G is not ionizable and not a terminal-specific residue in either scale).
    private static Func<double, double> BuildOptionsCharge(string proteinSequence, IsoelectricPointOptions options)
    {
        var scale = ToProteinPkaScale(options.Scale);
        string adjusted = NeutraliseResidues(proteinSequence, 'K', options.LysinesModified);
        adjusted = NeutraliseResidues(adjusted, 'C', 2L * options.Disulphides);

        if (options.Termini)
            return pH => ProteinPhysicochemistry.NetCharge(adjusted, pH, scale);

        string flanked = NeutralResidue + adjusted + NeutralResidue;
        string terminiOnly = new(NeutralResidue, 2);
        return pH => ProteinPhysicochemistry.NetCharge(flanked, pH, scale)
                     - ProteinPhysicochemistry.NetCharge(terminiOnly, pH, scale);
    }

    // embIepCompC: remove `count` residues of `residue` from the charged-group count, or all of them
    // when fewer are present. Case folding is ASCII-only (as ProteinPhysicochemistry).
    private static string NeutraliseResidues(string sequence, char residue, long count)
    {
        if (count == 0)
            return sequence;

        char lower = (char)(residue | 0x20);
        var chars = sequence.ToCharArray();
        long remaining = count;
        for (int i = 0; i < chars.Length && remaining > 0; i++)
        {
            if (chars[i] == residue || chars[i] == lower)
            {
                chars[i] = NeutralResidue;
                remaining--;
            }
        }

        return new string(chars);
    }

    // Same bracket as ProteinPhysicochemistry.IsoelectricPoint (net charge is non-increasing in pH).
    private static double BisectDefaultWindow(Func<double, double> charge)
    {
        double low = DefaultPiLowPh, high = DefaultPiHighPh;
        while (high - low > DefaultPiPrecision)
        {
            double mid = (low + high) / 2.0;
            if (charge(mid) > 0)
                low = mid;
            else
                high = mid;
        }

        return Math.Round((low + high) / 2.0, DefaultPiDecimalPlaces);
    }

    // EMBOSS 6.6.0 nucleus/embiep.c embIepPhConverge (and embIepIepC's "!*iep → none").
    private static bool TryBisectEmbossWindow(Func<double, double> charge, out double isoelectricPoint)
    {
        isoelectricPoint = double.NaN;
        double tph = EmbossPiLowPh, bph = EmbossPiHighPh;
        double top = charge(tph), bot = charge(bph);
        if ((top > 0.0 && bot > 0.0) || (top < 0.0 && bot < 0.0))
            return false;

        while (bph - tph > EmbossPiBracket)
        {
            double mid = ((bph - tph) / 2.0) + tph;
            double c = charge(mid);
            if (c > 0.0)
                tph = mid;
            else if (c < 0.0)
                bph = mid;
            else
            {
                isoelectricPoint = mid;
                return true;
            }
        }

        isoelectricPoint = tph;
        return true;
    }

    // Biopython 1.88 Bio/SeqUtils/IsoelectricPoint.py IsoelectricPoint.pi (recursion unrolled).
    private static double BisectBiopythonWindow(Func<double, double> charge)
    {
        double pH = BiopythonPiStartPh, min = BiopythonPiLowPh, max = BiopythonPiHighPh;
        while (max - min > BiopythonPiBracket)
        {
            if (charge(pH) > 0.0)
                min = pH;
            else
                max = pH;
            pH = (min + max) / 2;
        }

        return pH;
    }

    // Maps the (backward-compatible) public SequenceStatistics.PkaScale onto the canonical scale.
    private static ProteinPkaScale ToProteinPkaScale(PkaScale scale) => scale switch
    {
        PkaScale.Emboss => ProteinPkaScale.Emboss,
        PkaScale.Bjellqvist => ProteinPkaScale.Bjellqvist,
        _ => throw new ArgumentOutOfRangeException(nameof(scale), scale, "Unknown pK scale."),
    };

    #endregion

    #region Hydrophobicity

    // Kyte-Doolittle (1982) hydropathy index of the 20 standard residues (= Biopython
    // ProtParamData.kd); single source: ProteinPhysicochemistry.KyteDoolittle.
    private static readonly IReadOnlyDictionary<char, double> HydrophobicityScale = ProteinPhysicochemistry.KyteDoolittle;

    // Default sliding-window size; Kyte & Doolittle (1982) report window 9 gives the
    // best results for surface regions of globular proteins (per GCAT/Davidson summary).
    private const int DefaultHydropathyWindow = 9;

    /// <summary>
    /// Calculates the grand average of hydropathy (GRAVY) index: the sum of Kyte-Doolittle
    /// hydropathy values of all recognized residues divided by the residue count.
    /// Positive values indicate hydrophobic, negative hydrophilic. Non-standard residues
    /// are skipped (not counted). Returns 0 for null/empty input or no recognized residues.
    /// Delegates to the canonical <see cref="ProteinPhysicochemistry.Gravy"/> (ASCII-only case fold).
    /// </summary>
    public static double CalculateHydrophobicity(string proteinSequence) =>
        ProteinPhysicochemistry.Gravy(proteinSequence);

    /// <summary>
    /// Calculates the sliding-window Kyte-Doolittle hydropathy profile (Kyte &amp; Doolittle 1982;
    /// ExPASy ProtScale / Biopython <c>ProteinAnalysis.protein_scale(kd, window, edge)</c>).
    /// </summary>
    /// <remarks>
    /// <para>Yields exactly N − W + 1 values; value <c>k</c> (0-based) belongs to the window
    /// <c>[k, k + W − 1]</c> and, for odd W, is the score of its central residue <c>k + (W − 1)/2</c>
    /// (ProtScale convention). Yields nothing when W exceeds the sequence length or the input is
    /// null/empty.</para>
    /// <para>Weighting (ProtScale "linear" weight-variation model, as implemented by Biopython
    /// <c>_weight_list</c>): the central residue has weight 1, the two window ends have weight
    /// <paramref name="edgeWeight"/>, and weights vary linearly in between
    /// (<c>w_j = edge + j·2(1 − edge)/(W − 1)</c> for the j-th position from either end); each value is
    /// Σ w·kd / Σ w. With the default <paramref name="edgeWeight"/> = 1 this is the unweighted window
    /// mean of the original Kyte-Doolittle method. A weighted window needs a central residue, so
    /// <paramref name="edgeWeight"/> &lt; 1 requires an odd window (ProtScale accepts odd windows only).
    /// For an even window with edge 1 the plain mean over the W residues is returned (Biopython instead
    /// counts residue W/2 twice and divides by W + 1 — an artefact of its odd-window loop; reproduce it with
    /// <see cref="HydropathyWindowConvention.Biopython"/> in the options overload).</para>
    /// <para>Non-standard residues (B, Z, X, gaps, stop) have no scale value and contribute 0 to the
    /// weighted sum while keeping their weight in the divisor — identical to Biopython for such a
    /// residue at the window centre (Biopython also drops the symmetric partner of an off-centre
    /// unknown residue; this library does not — use <see cref="HydropathyWindowConvention.Biopython"/> for that, or
    /// <see cref="HydrophobicityProfileOptions.Strict"/> to reject unknown residues).</para>
    /// <para>Other scales, strict mode and ProtScale 0..1 scale normalisation:
    /// <see cref="CalculateHydrophobicityProfile(string, int, double, HydrophobicityProfileOptions)"/>.</para>
    /// </remarks>
    /// <param name="proteinSequence">One-letter amino-acid sequence (case-insensitive, ASCII-only fold:
    /// non-ASCII characters such as 'ſ' are non-residues).</param>
    /// <param name="windowSize">Window length W (≥ 1; default 9).</param>
    /// <param name="edgeWeight">Relative weight of the window edges, in [0, 1] (default 1 = unweighted).</param>
    /// <exception cref="ArgumentOutOfRangeException">W &lt; 1, or <paramref name="edgeWeight"/> outside [0, 1].</exception>
    /// <exception cref="ArgumentException"><paramref name="edgeWeight"/> &lt; 1 with an even W.</exception>
    public static IEnumerable<double> CalculateHydrophobicityProfile(
        string proteinSequence,
        int windowSize = DefaultHydropathyWindow,
        double edgeWeight = 1.0)
    {
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1.");
        if (double.IsNaN(edgeWeight) || edgeWeight < 0.0 || edgeWeight > 1.0)
            throw new ArgumentOutOfRangeException(nameof(edgeWeight), edgeWeight, "Edge weight must be in [0, 1].");
        if (edgeWeight < 1.0 && windowSize % 2 == 0)
            throw new ArgumentException("An edge-weighted window must have an odd size (a central residue).", nameof(windowSize));

        if (string.IsNullOrEmpty(proteinSequence) || windowSize > proteinSequence.Length)
            return Array.Empty<double>();

        return HydrophobicityProfileIterator(proteinSequence, windowSize, edgeWeight);
    }

    /// <summary>
    /// Window convention of <see cref="CalculateHydrophobicityProfile(string, int, double, HydrophobicityProfileOptions)"/>.
    /// </summary>
    public enum HydropathyWindowConvention
    {
        /// <summary>Default: every window position j has its own weight w_j, the value is Σ w·v / Σ w over the W
        /// residues; an unknown residue contributes 0 and keeps its weight; even windows need edge = 1.</summary>
        Standard = 0,

        /// <summary>Line-for-line Biopython 1.88 <c>ProteinAnalysis.protein_scale(param_dict, window, edge)</c>:
        /// positions j and W−1−j (j &lt; ⌊W/2⌋) are scored as a pair with the half-list weight
        /// <c>edge + j·2(1 − edge)/(W − 1)</c>, then residue ⌊W/2⌋ is added with weight 1, and the sum is divided by
        /// 2·Σ half-weights + 1. Hence for an even W residue W/2 is counted twice (once in its pair, once as the
        /// "middle") and the divisor is W + 1 at edge 1; even windows with edge &lt; 1 are accepted; a pair with an
        /// unknown residue is dropped as a whole (its known partner too), an unknown middle is dropped; the
        /// divisor never changes. W = 1 is rejected (Biopython raises ZeroDivisionError).</summary>
        Biopython = 1,
    }

    /// <summary>
    /// Options for <see cref="CalculateHydrophobicityProfile(string, int, double, HydrophobicityProfileOptions)"/>.
    /// <see cref="Default"/> reproduces <see cref="CalculateHydrophobicityProfile(string, int, double)"/> exactly.
    /// </summary>
    public sealed record HydrophobicityProfileOptions
    {
        /// <summary>Kyte-Doolittle scale, standard window convention, lenient, unnormalised.</summary>
        public static HydrophobicityProfileOptions Default { get; } = new();

        /// <summary>Window convention (default <see cref="HydropathyWindowConvention.Standard"/>).</summary>
        public HydropathyWindowConvention Convention { get; init; } = HydropathyWindowConvention.Standard;

        /// <summary>Residue scale keyed by upper-case one-letter code (e.g. Biopython <c>ProtParamData.hw</c>,
        /// <c>em</c>, <c>ja</c>, <c>Flex</c>); null (default) = Kyte-Doolittle
        /// (<see cref="ProteinPhysicochemistry.KyteDoolittle"/>). Sequence letters are upper-cased (ASCII) before lookup.</summary>
        public IReadOnlyDictionary<char, double>? Scale { get; init; }

        /// <summary>True: throw <see cref="ArgumentException"/> when the sequence contains a symbol absent from the
        /// scale (the Biopython <c>gravy()</c> KeyError behaviour). False (default): such symbols are handled by the
        /// window convention (no score).</summary>
        public bool Strict { get; init; }

        /// <summary>True: normalise the scale to 0..1 before profiling, v' = (v − min)/(max − min) over the scale's
        /// entries (ExPASy ProtScale "normalize the scale from 0 to 1", Gasteiger et al. 2005). Default false.</summary>
        public bool NormalizeScale { get; init; }
    }

    /// <summary>
    /// Sliding-window amino-acid-scale profile with an explicit scale, window convention, strictness and
    /// ProtScale 0..1 scale normalisation (see <see cref="HydrophobicityProfileOptions"/>). With
    /// <see cref="HydrophobicityProfileOptions.Default"/> it equals
    /// <see cref="CalculateHydrophobicityProfile(string, int, double)"/>.
    /// </summary>
    /// <param name="proteinSequence">One-letter amino-acid sequence (ASCII case fold).</param>
    /// <param name="windowSize">Window length W (≥ 1; ≥ 2 in the Biopython convention).</param>
    /// <param name="edgeWeight">Relative weight of the window edges, in [0, 1].</param>
    /// <param name="options">Scale / convention / strictness / normalisation; null = <see cref="HydrophobicityProfileOptions.Default"/>.</param>
    /// <exception cref="ArgumentOutOfRangeException">Bad W or edge weight; W = 1 in the Biopython convention.</exception>
    /// <exception cref="ArgumentException">Standard convention with edge &lt; 1 and even W; an empty or constant scale
    /// to normalise; strict mode and a symbol absent from the scale.</exception>
    public static IEnumerable<double> CalculateHydrophobicityProfile(
        string proteinSequence,
        int windowSize,
        double edgeWeight,
        HydrophobicityProfileOptions? options)
    {
        options ??= HydrophobicityProfileOptions.Default;
        bool biopython = options.Convention == HydropathyWindowConvention.Biopython;
        if (!Enum.IsDefined(options.Convention))
            throw new ArgumentOutOfRangeException(nameof(options), options.Convention, "Unknown window convention.");
        if (windowSize < (biopython ? 2 : 1))
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize,
                biopython ? "The Biopython convention needs a window of at least 2 (W = 1 divides by W − 1 = 0)." : "Window size must be at least 1.");
        if (double.IsNaN(edgeWeight) || edgeWeight < 0.0 || edgeWeight > 1.0)
            throw new ArgumentOutOfRangeException(nameof(edgeWeight), edgeWeight, "Edge weight must be in [0, 1].");
        if (!biopython && edgeWeight < 1.0 && windowSize % 2 == 0)
            throw new ArgumentException("An edge-weighted window must have an odd size (a central residue).", nameof(windowSize));

        IReadOnlyDictionary<char, double> scale = options.Scale ?? HydrophobicityScale;
        if (options.NormalizeScale)
            scale = NormalizeScaleZeroToOne(scale, nameof(options));

        if (string.IsNullOrEmpty(proteinSequence))
            return Array.Empty<double>();

        if (options.Strict)
        {
            foreach (char c in proteinSequence)
            {
                if (!scale.ContainsKey(ToUpperAscii(c)))
                    throw new ArgumentException($"Symbol '{c}' has no value in the scale (strict mode).", nameof(proteinSequence));
            }
        }

        if (windowSize > proteinSequence.Length)
            return Array.Empty<double>();

        return biopython
            ? BiopythonScaleProfileIterator(proteinSequence, windowSize, edgeWeight, scale)
            : HydrophobicityProfileIterator(proteinSequence, windowSize, edgeWeight, scale);
    }

    /// <summary>
    /// Sliding-window profile of a caller-supplied amino-acid scale (e.g. Biopython <c>ProtParamData.hw</c>
    /// Hopp-Woods, <c>em</c> Eisenberg, <c>ja</c> Janin, <c>Flex</c> Vihinen flexibility) with the standard
    /// window convention; shorthand for the options overload with <see cref="HydrophobicityProfileOptions.Scale"/>.
    /// </summary>
    public static IEnumerable<double> CalculateHydrophobicityProfile(
        string proteinSequence,
        IReadOnlyDictionary<char, double> scale,
        int windowSize = DefaultHydropathyWindow,
        double edgeWeight = 1.0)
    {
        ArgumentNullException.ThrowIfNull(scale);
        return CalculateHydrophobicityProfile(proteinSequence, windowSize, edgeWeight,
            new HydrophobicityProfileOptions { Scale = scale });
    }

    // ExPASy ProtScale "normalize the scale from 0 to 1": min–max over the scale's own entries.
    private static IReadOnlyDictionary<char, double> NormalizeScaleZeroToOne(IReadOnlyDictionary<char, double> scale, string paramName)
    {
        if (scale.Count == 0)
            throw new ArgumentException("Cannot normalise an empty scale.", paramName);
        double min = scale.Values.Min(), max = scale.Values.Max();
        if (!(max > min))
            throw new ArgumentException("Cannot normalise a constant scale (max = min).", paramName);
        var normalized = new Dictionary<char, double>(scale.Count);
        foreach (var (key, value) in scale)
            normalized[key] = (value - min) / (max - min);
        return normalized;
    }

    private static IEnumerable<double> HydrophobicityProfileIterator(
        string sequence, int windowSize, double edgeWeight, IReadOnlyDictionary<char, double>? scale = null)
    {
        scale ??= HydrophobicityScale;
        // Position weights (all 1 for the unweighted Kyte-Doolittle mean).
        var weights = new double[windowSize];
        double weightSum = 0;
        int half = windowSize / 2;
        for (int j = 0; j < windowSize; j++)
        {
            int fromEdge = Math.Min(j, windowSize - 1 - j);
            weights[j] = edgeWeight >= 1.0 || fromEdge >= half
                ? 1.0
                : edgeWeight + fromEdge * 2.0 * (1.0 - edgeWeight) / (windowSize - 1);
            weightSum += weights[j];
        }

        for (int i = 0; i <= sequence.Length - windowSize; i++)
        {
            double sum = 0;
            for (int j = 0; j < windowSize; j++)
            {
                if (scale.TryGetValue(ToUpperAscii(sequence[i + j]), out double value))
                    sum += weights[j] * value;
            }
            yield return sum / weightSum;
        }
    }

    // Line-for-line port of Biopython 1.88 ProteinAnalysis.protein_scale / _weight_list (W ≥ 2).
    private static IEnumerable<double> BiopythonScaleProfileIterator(
        string sequence, int windowSize, double edgeWeight, IReadOnlyDictionary<char, double> scale)
    {
        int half = windowSize / 2;
        double unit = 2.0 * (1.0 - edgeWeight) / (windowSize - 1);
        var weights = new double[half];
        double halfSum = 0;
        for (int j = 0; j < half; j++)
        {
            weights[j] = edgeWeight + unit * j;
            halfSum += weights[j];
        }
        double sumOfWeights = halfSum * 2 + 1;

        for (int i = 0; i <= sequence.Length - windowSize; i++)
        {
            double score = 0.0;
            for (int j = 0; j < half; j++)
            {
                if (scale.TryGetValue(ToUpperAscii(sequence[i + j]), out double front)
                    && scale.TryGetValue(ToUpperAscii(sequence[i + windowSize - j - 1]), out double back))
                    score += weights[j] * front + weights[j] * back;
            }
            if (scale.TryGetValue(ToUpperAscii(sequence[i + half]), out double middle))
                score += middle;
            yield return score / sumOfWeights;
        }
    }

    // ASCII-only case folding, the convention of the canonical ProteinPhysicochemistry (B02
    // TRANS-PROT-001): char.ToUpperInvariant folds U+017F 'ſ' → 'S' and U+0131 'ı' → 'I', which
    // EMBOSS treats as non-residues (pepstats "MSKſ" → 3 residues, MW 364.46 = MW("MSK")).
    private static char ToUpperAscii(char c) => char.IsAsciiLetterLower(c) ? (char)(c - 32) : c;

    #endregion

    #region DNA Thermodynamics

    // The nearest-neighbour model is the canonical Infrastructure core
    // ThermoConstants.CalculateNearestNeighborDuplex (a line-by-line port of Biopython 1.88
    // Bio.SeqUtils.MeltingTemp.Tm_NN; 2026-10 FIN-B03 F26). SequenceStatistics only maps its
    // mol/L inputs onto Tm_NN's nM / mM arguments, keeps its "< 2 bases → zeros" sentinel and
    // derives ΔG°₃₇ = ΔH° − 310.15·ΔS°/1000.

    // Reference temperature for ΔG° calculation: 37 °C = 310.15 K.
    private const double ReferenceTemperatureKelvin = 310.15;

    /// <summary>
    /// DNA thermodynamic properties.
    /// </summary>
    public readonly record struct ThermodynamicProperties(
        double DeltaH,
        double DeltaS,
        double DeltaG,
        double MeltingTemperature);

    /// <summary>
    /// Reaction conditions beyond [Na⁺] / C_T for
    /// <see cref="CalculateThermodynamics(string, double, double, bool, ThermodynamicsOptions)"/>:
    /// the Biopython <c>MeltingTemp.salt_correction</c> method (0–7) with its K⁺ / Tris / Mg²⁺ / dNTP
    /// inputs, and whether the result is rounded. The defaults (method 5, no other ions, rounded) give
    /// exactly <see cref="CalculateThermodynamics(string, double, double, bool)"/>.
    /// </summary>
    public sealed record ThermodynamicsOptions
    {
        /// <summary>Default conditions: SantaLucia (1998) method-5 entropy correction, Na⁺ only, rounded output.</summary>
        public static ThermodynamicsOptions Default { get; } = new();

        /// <summary>Biopython <c>saltcorr</c> method (default 5, SantaLucia 1998 entropy). 1–4 add to Tm,
        /// 5 adds to ΔS°, 6 (Owczarzy 2004) and 7 (Owczarzy 2008, Mg²⁺/dNTPs) add to 1/Tm; 0 = none.</summary>
        public NnSaltCorrection SaltCorrection { get; init; } = NnSaltCorrection.SantaLucia1998Entropy;

        /// <summary>[K⁺] in mol/L (default 0); ≥ 0.</summary>
        public double PotassiumConcentration { get; init; }

        /// <summary>[Tris] in mol/L (default 0); ≥ 0. Counts as [Tris]/2 monovalent ions (Biopython).</summary>
        public double TrisConcentration { get; init; }

        /// <summary>[Mg²⁺] in mol/L (default 0); ≥ 0.</summary>
        public double MagnesiumConcentration { get; init; }

        /// <summary>[dNTPs] in mol/L (default 0); ≥ 0.</summary>
        public double DntpConcentration { get; init; }

        /// <summary>Round ΔH°/ΔS°/ΔG° to 2 decimals and Tm to 1 decimal (default true); false returns
        /// the unrounded values (Biopython <c>Tm_NN</c> returns an unrounded Tm).</summary>
        public bool RoundResults { get; init; } = true;
    }

    /// <summary>
    /// Calculates thermodynamic properties (ΔH°, ΔS°, ΔG°₃₇ and Tm) of a DNA duplex formed by
    /// two non-self-complementary strands in equal amount, using the nearest-neighbor model of
    /// Allawi &amp; SantaLucia (1997) (Biopython <c>DNA_NN3</c>) with the SantaLucia (1998)
    /// "method 5" Na+ entropy correction. Reproduces Biopython
    /// <c>MeltingTemp.Tm_NN(seq, nn_table=DNA_NN3, Na=1000·[Na+], dnac1=dnac2=C_T/2, saltcorr=5)</c>.
    /// Equivalent to <see cref="CalculateThermodynamics(string, double, double, bool)"/> with
    /// <c>selfComplementary: false</c>.
    /// </summary>
    /// <param name="dnaSequence">DNA sequence (5'→3'). Normalised like Biopython <c>Tm_NN</c>
    /// (<c>_check</c>): case-insensitive, whitespace removed, RNA U read as T, and every character
    /// other than A/C/G/T and inosine I removed before the model is applied. Inosine is scored with
    /// Biopython's <c>DNA_IMM1</c> parameters (SantaLucia 2005); a neighbour pair without parameters
    /// (e.g. "II") throws, as Biopython's default strict mode does. Fewer than 2 remaining bases → all zero.</param>
    /// <param name="naConcentration">Na+ concentration in mol/L (default 0.05 = 50 mM); must be &gt; 0.</param>
    /// <param name="primerConcentration">
    /// Total strand concentration C_T in mol/L (default 2.5e-7 = 250 nM; must be &gt; 0); the Tm
    /// equation divides this by F = 4 for two non-self-complementary strands in equal amount.
    /// </param>
    /// <returns>
    /// ΔH° (kcal/mol), ΔS° (cal/(mol·K), salt-corrected), ΔG°₃₇ (kcal/mol) rounded to 2 decimals
    /// and Tm (°C) rounded to 1 decimal. For null input or fewer than 2 A/C/G/T/I (U) bases all four
    /// fields are 0.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">A concentration is not a positive finite number.</exception>
    /// <exception cref="ArgumentException">A neighbour pair has no nearest-neighbour parameter (only possible
    /// with inosine, e.g. "II"; Biopython <c>Tm_NN</c> raises <c>ValueError</c>).</exception>
    public static ThermodynamicProperties CalculateThermodynamics(
        string dnaSequence,
        double naConcentration = 0.05, // 50 mM
        double primerConcentration = 0.00000025) // 250 nM
        => CalculateThermodynamics(dnaSequence, naConcentration, primerConcentration, selfComplementary: false);

    /// <summary>
    /// Calculates ΔH°, ΔS°, ΔG°₃₇ and Tm of a DNA duplex with the Allawi &amp; SantaLucia (1997)
    /// nearest-neighbor model (Biopython <c>DNA_NN3</c>), SantaLucia (1998) method-5 salt
    /// correction, and an explicit choice of the duplex stoichiometry.
    /// </summary>
    /// <param name="dnaSequence">DNA sequence (5'→3'); normalised as in
    /// <see cref="CalculateThermodynamics(string, double, double)"/>.</param>
    /// <param name="naConcentration">Na+ concentration in mol/L; must be &gt; 0.</param>
    /// <param name="primerConcentration">Strand concentration C_T in mol/L; must be &gt; 0.
    /// Non-self-complementary: total of two equimolar strands (Tm uses C_T/4).
    /// Self-complementary: concentration of the single strand (Tm uses C_T/1).</param>
    /// <param name="selfComplementary">
    /// <c>true</c> for a self-complementary (homo)duplex: adds the symmetry correction
    /// ΔS° −1.4 cal/(mol·K) and uses x = 1 in the Tm equation (SantaLucia 1998; Biopython
    /// <c>Tm_NN(selfcomp=True)</c>, which takes <c>k = dnac1</c>). Like Biopython, the flag is the
    /// caller's statement of the experiment and is not inferred from the sequence.
    /// </param>
    /// <exception cref="ArgumentOutOfRangeException">A concentration is not a positive finite number.</exception>
    /// <exception cref="ArgumentException">A neighbour pair has no nearest-neighbour parameter (inosine only).</exception>
    public static ThermodynamicProperties CalculateThermodynamics(
        string dnaSequence,
        double naConcentration,
        double primerConcentration,
        bool selfComplementary)
    {
        // Biopython salt_correction raises for a zero ion concentration and math.log for a
        // negative one; C_T enters ln(C_T/x). Reject non-physical inputs instead of returning
        // NaN / ±∞ / −273.15 °C.
        if (!(naConcentration > 0) || double.IsInfinity(naConcentration))
            throw new ArgumentOutOfRangeException(nameof(naConcentration), naConcentration,
                "Na+ concentration must be a positive finite value in mol/L.");
        return CalculateThermodynamics(dnaSequence, naConcentration, primerConcentration, selfComplementary,
            ThermodynamicsOptions.Default);
    }

    /// <summary>
    /// Calculates ΔH°, ΔS°, ΔG°₃₇ and Tm of a DNA duplex (Allawi &amp; SantaLucia 1997, Biopython
    /// <c>DNA_NN3</c>) under any Biopython <c>MeltingTemp.salt_correction</c> method (0–7) with
    /// Na⁺ / K⁺ / Tris / Mg²⁺ / dNTPs, optionally unrounded. Reproduces Biopython 1.88
    /// <c>Tm_NN(seq, Na=1000·[Na⁺], K=, Tris=, Mg=, dNTPs= (mM), dnac1=dnac2=C_T/2 (nM; self-complementary:
    /// dnac1=C_T, dnac2=0), selfcomp=, saltcorr=method)</c> through the canonical
    /// <see cref="ThermoConstants.CalculateNearestNeighborDuplex"/>.
    /// </summary>
    /// <param name="dnaSequence">DNA sequence (5'→3'); normalised as in
    /// <see cref="CalculateThermodynamics(string, double, double)"/>.</param>
    /// <param name="naConcentration">Na⁺ concentration in mol/L; ≥ 0 and finite (a zero total ion
    /// concentration is rejected by methods 1–6, as in Biopython).</param>
    /// <param name="primerConcentration">Strand concentration C_T in mol/L; must be &gt; 0 (see the 4-argument overload).</param>
    /// <param name="selfComplementary">Self-complementary (homo)duplex (see the 4-argument overload).</param>
    /// <param name="options">Salt-correction method, other ions and rounding; <c>null</c> = <see cref="ThermodynamicsOptions.Default"/>.</param>
    /// <returns>
    /// ΔH° (kcal/mol, 1 M NaCl); ΔS° (cal/(mol·K)) — salt-corrected for method 5, the 1 M NaCl value for every
    /// other method (methods 1–4 / 6–7 correct Tm or 1/Tm, not ΔS°); ΔG°₃₇ = ΔH° − 310.15·ΔS°/1000 from those two;
    /// Tm (°C) with the chosen correction. Rounded to 2 / 1 decimals unless
    /// <see cref="ThermodynamicsOptions.RoundResults"/> is false. All zero for fewer than 2 A/C/G/T/I (U) bases.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">C_T not positive finite, an ion concentration negative or
    /// not finite, or a method outside 0–7.</exception>
    /// <exception cref="ArgumentException">A neighbour pair without parameters (inosine) or a salt input Biopython
    /// rejects (zero total ion concentration for methods 1–6, ln of a non-positive concentration).</exception>
    public static ThermodynamicProperties CalculateThermodynamics(
        string dnaSequence,
        double naConcentration,
        double primerConcentration,
        bool selfComplementary,
        ThermodynamicsOptions? options)
    {
        options ??= ThermodynamicsOptions.Default;
        if (!(primerConcentration > 0) || double.IsInfinity(primerConcentration))
            throw new ArgumentOutOfRangeException(nameof(primerConcentration), primerConcentration,
                "Strand concentration must be a positive finite value in mol/L.");
        RequireNonNegativeFinite(naConcentration, nameof(naConcentration));
        RequireNonNegativeFinite(options.PotassiumConcentration, nameof(options.PotassiumConcentration));
        RequireNonNegativeFinite(options.TrisConcentration, nameof(options.TrisConcentration));
        RequireNonNegativeFinite(options.MagnesiumConcentration, nameof(options.MagnesiumConcentration));
        RequireNonNegativeFinite(options.DntpConcentration, nameof(options.DntpConcentration));
        if (options.SaltCorrection is < NnSaltCorrection.None or > NnSaltCorrection.Owczarzy2008)
            throw new ArgumentOutOfRangeException(nameof(options), options.SaltCorrection,
                "Allowed salt-correction methods are 0–7.");

        // Sentinel kept from the original API: fewer than two scorable bases → zeros (Biopython has no
        // duplex either: a 1-nt Tm_NN has no neighbour pair).
        if (ThermoConstants.NormalizeForNearestNeighbor(dnaSequence).Length < 2)
            return new ThermodynamicProperties(0, 0, 0, 0);

        // mol/L → Tm_NN units: strands in nM, ions in mM. Equimolar strands dnac1 = dnac2 = C_T/2
        // give k = dnac1 − dnac2/2 = C_T/4; self-complementary k = dnac1 = C_T.
        double ctNanomolar = primerConcentration * 1e9;
        var r = ThermoConstants.CalculateNearestNeighborDuplex(
            dnaSequence,
            parameterSet: NnParameterSet.AllawiSantaLucia1997,
            dnac1: selfComplementary ? ctNanomolar : ctNanomolar / 2,
            dnac2: selfComplementary ? 0 : ctNanomolar / 2,
            selfComplementary: selfComplementary,
            sodium: naConcentration * 1000,
            potassium: options.PotassiumConcentration * 1000,
            tris: options.TrisConcentration * 1000,
            magnesium: options.MagnesiumConcentration * 1000,
            dntps: options.DntpConcentration * 1000,
            saltCorrection: options.SaltCorrection);

        // ΔG°₃₇ = ΔH° − T·ΔS° (ΔS° converted from cal to kcal).
        double dG = r.DeltaH - (ReferenceTemperatureKelvin * r.SaltCorrectedDeltaS / 1000.0);
        return options.RoundResults
            ? new ThermodynamicProperties(
                DeltaH: Math.Round(r.DeltaH, 2),
                DeltaS: Math.Round(r.SaltCorrectedDeltaS, 2),
                DeltaG: Math.Round(dG, 2),
                MeltingTemperature: Math.Round(r.MeltingTemperature, 1))
            : new ThermodynamicProperties(r.DeltaH, r.SaltCorrectedDeltaS, dG, r.MeltingTemperature);

        static void RequireNonNegativeFinite(double value, string name)
        {
            if (!(value >= 0) || double.IsInfinity(value))
                throw new ArgumentOutOfRangeException(name, value,
                    "Ion concentration must be a non-negative finite value in mol/L.");
        }
    }

    /// <summary>
    /// Calculates the "basic" melting temperature of an oligonucleotide (OligoCalc, Kibbe 2007,
    /// NAR 35:W43): Wallace rule Tm = 2(A+T) + 4(G+C) (Thein &amp; Wallace 1986) for fewer than 14
    /// bases, otherwise Tm = 64.9 + 41·(G+C − 16.4)/N (customarily attributed to Marmur &amp; Doty 1962),
    /// both at OligoCalc's fixed standard conditions (50 nM primer, 50 mM Na+, pH 7.0).
    /// </summary>
    /// <remarks>
    /// Only A, C, G, T and U are counted (case-insensitive); RNA uracil is read as T, as Biopython
    /// <c>MeltingTemp._check</c> back-transcribes RNA for <c>Tm_Wallace</c>/<c>Tm_GC</c> (and as
    /// <see cref="CalculateThermodynamics(string, double, double)"/> does). Every other character
    /// (N, IUPAC codes, gaps, whitespace) is ignored. As in OligoCalc — whose GC formula divides by
    /// N = wA+xT+yG+zC — the length that selects the formula is the number of counted bases, so for
    /// <paramref name="useWallaceRule"/> = true delegates to the canonical
    /// <see cref="ThermoConstants.CalculateBasicTm"/> (also behind <c>PrimerDesigner.CalculateMeltingTemperature</c>);
    /// false uses <see cref="ThermoConstants.CountBasicTmBases"/> + <see cref="ThermoConstants.CalculateMarmurDotyTm"/>.
    /// With <paramref name="useWallaceRule"/> = false the GC formula is applied at any length; below
    /// 14 bases this is outside its published domain and can be negative. Returns 0 for null/empty
    /// input or when no A/C/G/T base is present.
    /// </remarks>
    public static double CalculateMeltingTemperature(string dnaSequence, bool useWallaceRule = true)
    {
        // Canonical basic Tm (Infrastructure, B07): same A/C/G/T/U count with U read as T
        // (Biopython MeltingTemp._check; 2026-09 B03 F20) and the same 14-base Wallace switch
        // (2026-10 FIN-B03 F25).
        if (useWallaceRule)
            return ThermoConstants.CalculateBasicTm(dnaSequence);

        // GC formula (Marmur-Doty / OligoCalc basic Tm) at any length.
        var (at, gc) = ThermoConstants.CountBasicTmBases(dnaSequence);
        return ThermoConstants.CalculateMarmurDotyTm(gc, at + gc);
    }

    /// <summary>
    /// Basic melting temperature with optional IUPAC-ambiguity weighting as Biopython 1.88 non-strict
    /// <c>Bio.SeqUtils.MeltingTemp.Tm_Wallace(seq, strict=False)</c> /
    /// <c>Tm_GC(seq, userset=(64.9, 0.41, 672.4, 0), saltcorr=0, strict=False)</c>.
    /// </summary>
    /// <remarks>
    /// <paramref name="weightAmbiguousBases"/> = false is exactly
    /// <see cref="CalculateMeltingTemperature(string, bool)"/>. When true:
    /// <list type="bullet">
    /// <item>The formula is selected as in the default: Wallace when <paramref name="useWallaceRule"/> is
    /// true and fewer than 14 A/C/G/T/U bases are present (OligoCalc's N = wA+xT+yG+zC), the GC formula
    /// otherwise — so weighting never changes which formula applies.</item>
    /// <item>Wallace = Biopython <c>Tm_Wallace(strict=False)</c>: 2·(A+T+W) + 4·(C+G+S) + 3·(K+M+N+R+Y)
    /// + 10/3·(B+V) + 8/3·(D+H) (U read as T; other characters ignored).</item>
    /// <item>GC formula = Biopython <c>Tm_GC(strict=False, mismatch=True)</c> with the OligoCalc userset:
    /// the sequence is reduced to the <c>Tm_GC</c> alphabet A B C D G H I K M N R S T V W X Y (U read as T),
    /// N is its length, %GC = 100·<c>gc_fraction(seq, "weighted")</c> − 50·X/N (canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(string, SequenceExtensions.GcAmbiguityMode)"/>,
    /// Weighted mode), Tm = 64.9 + 0.41·%GC − 672.4/N. X and I therefore count in N but add no GC.</item>
    /// </list>
    /// Case folding is ASCII-only (Python <c>str.upper</c> also maps U+0131 'ı' to I). Returns 0 for
    /// null/empty input and when the selected formula has no base to score.
    /// </remarks>
    /// <param name="dnaSequence">Oligonucleotide (DNA or RNA), case-insensitive.</param>
    /// <param name="useWallaceRule">True: Wallace below 14 A/C/G/T/U bases, GC formula otherwise; false: GC formula.</param>
    /// <param name="weightAmbiguousBases">True to weight IUPAC ambiguity codes as Biopython non-strict.</param>
    /// <returns>Tm in °C.</returns>
    public static double CalculateMeltingTemperature(string dnaSequence, bool useWallaceRule, bool weightAmbiguousBases)
    {
        if (!weightAmbiguousBases)
            return CalculateMeltingTemperature(dnaSequence, useWallaceRule);
        if (string.IsNullOrEmpty(dnaSequence))
            return 0;

        var (at, gc) = ThermoConstants.CountBasicTmBases(dnaSequence);
        return useWallaceRule && at + gc < ThermoConstants.WallaceMaxLength
            ? AmbiguityWeightedWallaceTm(dnaSequence)
            : AmbiguityWeightedGcTm(dnaSequence);
    }

    // Biopython 1.88 MeltingTemp.Tm_Wallace(seq, check=True, strict=False), term for term.
    private static double AmbiguityWeightedWallaceTm(string sequence)
    {
        int weak = 0, strong = 0, half = 0, twoThirds = 0, oneThird = 0;
        foreach (char raw in sequence)
        {
            switch (ToUpperAscii(raw))
            {
                case 'A': case 'T': case 'U': case 'W': weak++; break;
                case 'C': case 'G': case 'S': strong++; break;
                case 'K': case 'M': case 'N': case 'R': case 'Y': half++; break;
                case 'B': case 'V': twoThirds++; break;
                case 'D': case 'H': oneThird++; break;
            }
        }
        double meltingTemp = ThermoConstants.WallaceAtContribution * weak + ThermoConstants.WallaceGcContribution * strong;
        double ambiguous = 3 * half + 10 / 3.0 * twoThirds + 8 / 3.0 * oneThird;
        return meltingTemp + ambiguous;
    }

    // Biopython 1.88 MeltingTemp.Tm_GC(seq, userset=(64.9, 0.41, 672.4, 0), saltcorr=0, strict=False):
    // _check keeps the Tm_GC alphabet (RNA back-transcribed), %GC = gc_fraction(weighted)·100 − X·50/N.
    private static double AmbiguityWeightedGcTm(string sequence)
    {
        var kept = new System.Text.StringBuilder(sequence.Length);
        int countX = 0;
        foreach (char raw in sequence)
        {
            char c = ToUpperAscii(raw);
            if (c == 'U') c = 'T';
            if ("ABCDGHIKMNRSTVWXY".IndexOf(c) < 0)
                continue;
            if (c == 'X') countX++;
            kept.Append(c);
        }
        int length = kept.Length;
        if (length == 0)
            return 0;

        double percentGc = kept.ToString().CalculateGcFraction(SequenceExtensions.GcAmbiguityMode.Weighted) * 100;
        percentGc -= countX * 50.0 / length;
        return ThermoConstants.MarmurDotyBase
            + ThermoConstants.MarmurDotyGcCoefficient / 100 * percentGc
            - ThermoConstants.MarmurDotyGcCoefficient * ThermoConstants.MarmurDotyGcOffset / length;
    }

    #endregion

    #region Sequence Patterns

    /// <summary>
    /// Calculates normalized dinucleotide frequencies f_XY = count(XY) / (number of counted
    /// overlapping dinucleotides) over the alphabet {A,T,G,C,U}. Pairs containing any other symbol
    /// (N, IUPAC ambiguity, gap) are excluded from both the count and the denominator, so the
    /// denominator equals N−1 only when every base is in the alphabet. This is seqinr
    /// <c>count(seq, 2, freq = TRUE)</c> (denominator = sum of in-alphabet word counts) and, for
    /// pure-alphabet input, EMBOSS <c>compseq -word 2</c> (Obs Frequency over N−1 words).
    /// Frequency normalization follows the Karlin genomic-signature convention
    /// (Karlin S., "Pervasive properties of the genomic signature", PMC126251).
    /// Case-insensitive; only observed dinucleotides are keys (absent key ⇒ count 0).
    /// Counting delegates to the canonical overlapping k-mer counter
    /// <see cref="KmerAnalyzer.CountKmers(string, int)"/> with k = 2.
    /// </summary>
    public static IReadOnlyDictionary<string, double> CalculateDinucleotideFrequencies(string sequence)
    {
        var freq = new Dictionary<string, double>();

        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return freq;

        var counts = CountAlphabetDinucleotides(sequence, DinucleotideAlphabet);
        int total = counts.Values.Sum();

        foreach (var (dinuc, count) in counts)
            freq[dinuc] = (double)count / total;

        return freq;
    }

    // Single-strand dinucleotide alphabet (DNA + RNA) and the double-stranded DNA alphabet used for
    // the strand-symmetrized ρ* (Karlin; seqinr default alphabet "acgt").
    private const string DinucleotideAlphabet = "ACGTU";
    private const string DoubleStrandedDinucleotideAlphabet = "ACGT";

    /// <summary>
    /// Overlapping dinucleotide counts restricted to pairs whose two symbols are in
    /// <paramref name="alphabet"/> (upper case). Delegates to the canonical k-mer counter.
    /// </summary>
    private static Dictionary<string, int> CountAlphabetDinucleotides(string sequence, string alphabet)
    {
        var counts = new Dictionary<string, int>();
        foreach (var (kmer, count) in KmerAnalyzer.CountKmers(sequence, 2))
        {
            if (alphabet.Contains(kmer[0]) && alphabet.Contains(kmer[1]))
                counts[kmer] = count;
        }

        return counts;
    }

    /// <summary>
    /// Calculates the dinucleotide relative abundance (odds ratio) ρ_XY = f_XY / (f_X · f_Y),
    /// where f_XY is the dinucleotide frequency and f_X, f_Y are single-base frequencies.
    /// ρ = 1 indicates no bias (observed equals the product of base frequencies); ρ &gt; 1
    /// over-representation and ρ &lt; 1 under-representation
    /// (Karlin S., PMC126251; Karlin &amp; Burge 1995, Trends Genet 11(7):283-290).
    /// Single-strand form; equivalent to <see cref="CalculateDinucleotideRatios(string, bool)"/>
    /// with <c>strandSymmetric: false</c>.
    /// </summary>
    public static IReadOnlyDictionary<string, double> CalculateDinucleotideRatios(string sequence)
        => CalculateDinucleotideRatios(sequence, strandSymmetric: false);

    /// <summary>
    /// Calculates dinucleotide relative abundances (odds ratios), either single-strand
    /// ρ_XY = f_XY / (f_X · f_Y) or Karlin's strand-symmetrized genomic-signature form ρ*_XY.
    /// </summary>
    /// <remarks>
    /// <para><b>Single strand</b> (<paramref name="strandSymmetric"/> = false): alphabet {A,C,G,T,U};
    /// f_XY = count(XY) / (counted in-alphabet dinucleotides) and f_X = count(X) / (count of
    /// in-alphabet bases); other symbols (N, ambiguity codes, gaps) are excluded from counts and
    /// denominators. Identical to seqinr <c>rho(seq)</c> (wordcount / (Σwordcount · f_X f_Y)) and,
    /// for pure-ACGT input, to EMBOSS <c>compseq -word 2 -calcfreq</c> Obs/Exp.</para>
    /// <para><b>Strand-symmetrized</b> (true): Karlin's ρ*_XY for double-stranded DNA, computed from
    /// the sequence together with its inverted complement: f*_A = f*_T = (f_A + f_T)/2,
    /// f*_C = f*_G = (f_C + f_G)/2, f*_XY = (f_XY + f_{X̄'Ȳ'})/2 where X̄'Ȳ' is the reverse complement of
    /// XY (e.g. f*_GT = (f_GT + f_AC)/2), and ρ*_XY = f*_XY / (f*_X f*_Y) (Karlin &amp; Mrázek 1997;
    /// Karlin 1998, PMC126251). No junction dinucleotide is created (equals seqinr <c>rho</c> of the
    /// sequence, a separator and its reverse complement). Alphabet {A,C,G,T}; U is not a
    /// double-stranded DNA base and is excluded like N.</para>
    /// <para>Only dinucleotides that occur (on either strand, for ρ*) are keys. A dinucleotide that
    /// occurs always has both constituent bases present, so every returned ratio is finite and &gt; 0;
    /// an absent key means ρ = 0 when both bases occur (seqinr/compseq report 0) and undefined
    /// (0/0) otherwise. This odds ratio normalizes the dinucleotide count by the counted
    /// dinucleotides (≈ N−1); the Gardiner-Garden &amp; Frommer CpG O/E used for CpG-island
    /// calling normalizes by the window length N — see
    /// <c>EpigeneticsAnalyzer.CalculateCpGObservedExpected</c>.</para>
    /// </remarks>
    /// <param name="sequence">Nucleotide sequence (case-insensitive).</param>
    /// <param name="strandSymmetric">true for Karlin's strand-symmetrized ρ*; false for single-strand ρ.</param>
    /// <returns>Map dinucleotide → odds ratio; empty for null/empty/length &lt; 2 or no counted dinucleotide.</returns>
    public static IReadOnlyDictionary<string, double> CalculateDinucleotideRatios(string sequence, bool strandSymmetric)
    {
        var ratios = new Dictionary<string, double>();

        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return ratios;

        var comp = CalculateNucleotideComposition(sequence);

        if (!strandSymmetric)
        {
            int total = comp.CountA + comp.CountT + comp.CountG + comp.CountC + comp.CountU;
            if (total == 0) return ratios;

            var singleFreq = new Dictionary<char, double>
            {
                { 'A', (double)comp.CountA / total },
                { 'T', (double)comp.CountT / total },
                { 'G', (double)comp.CountG / total },
                { 'C', (double)comp.CountC / total },
                { 'U', (double)comp.CountU / total }
            };

            foreach (var (dinuc, observed) in CalculateDinucleotideFrequencies(sequence))
            {
                double expected = singleFreq[dinuc[0]] * singleFreq[dinuc[1]];
                ratios[dinuc] = observed / expected;
            }

            return ratios;
        }

        // Karlin ρ*: counts of the sequence plus its inverted complement.
        var counts = CountAlphabetDinucleotides(sequence, DoubleStrandedDinucleotideAlphabet);
        int dinucTotal = counts.Values.Sum();
        int baseTotal = comp.CountA + comp.CountT + comp.CountG + comp.CountC;
        if (dinucTotal == 0 || baseTotal == 0) return ratios;

        double fAT = (double)(comp.CountA + comp.CountT) / (2.0 * baseTotal);
        double fCG = (double)(comp.CountC + comp.CountG) / (2.0 * baseTotal);
        double SymmetricBaseFrequency(char b) => b is 'A' or 'T' ? fAT : fCG;

        var symmetricCounts = new Dictionary<string, int>();
        foreach (var (dinuc, count) in counts)
        {
            string reverseComplement = new(new[]
            {
                Seqeron.Genomics.Core.SequenceExtensions.GetComplementBase(dinuc[1]),
                Seqeron.Genomics.Core.SequenceExtensions.GetComplementBase(dinuc[0])
            });
            symmetricCounts[dinuc] = symmetricCounts.GetValueOrDefault(dinuc) + count;
            symmetricCounts[reverseComplement] = symmetricCounts.GetValueOrDefault(reverseComplement) + count;
        }

        foreach (var (dinuc, count) in symmetricCounts)
        {
            double observed = count / (2.0 * dinucTotal);
            ratios[dinuc] = observed / (SymmetricBaseFrequency(dinuc[0]) * SymmetricBaseFrequency(dinuc[1]));
        }

        return ratios;
    }

    /// <summary>
    /// Karlin's representation class of a dinucleotide relative abundance (odds ratio) ρ / ρ*.
    /// </summary>
    public enum DinucleotideRepresentation
    {
        /// <summary>0.78 &lt; ρ &lt; 1.23 — within the normal range.</summary>
        Normal = 0,

        /// <summary>ρ ≤ 0.78 — under-represented (Karlin &amp; Burge 1995).</summary>
        UnderRepresented = 1,

        /// <summary>ρ ≥ 1.23 — over-represented (Karlin &amp; Burge 1995).</summary>
        OverRepresented = 2,
    }

    /// <summary>Karlin under-representation threshold: ρ ≤ 0.78 (Karlin &amp; Burge 1995).</summary>
    public const double KarlinUnderRepresentedThreshold = 0.78;

    /// <summary>Karlin over-representation threshold: ρ ≥ 1.23 (Karlin &amp; Burge 1995).</summary>
    public const double KarlinOverRepresentedThreshold = 1.23;

    /// <summary>
    /// Classifies a dinucleotide relative abundance with Karlin's criterion: under-represented when
    /// ρ ≤ 0.78, over-represented when ρ ≥ 1.23, normal otherwise (both bounds inclusive).
    /// </summary>
    /// <remarks>
    /// Criterion of Karlin &amp; Burge (1995), Trends Genet 11(7):283–290, as quoted verbatim in
    /// MBE 19(6):964: "the XY dinucleotide was considered to be underrepresented if ρ*XY ≤ 0.78 and
    /// overrepresented if ρ*XY ≥ 1.23" (also Karlin, Mrázek &amp; Campbell 1997; Karlin 1998). Karlin
    /// states it for the strand-symmetrized ρ*; it applies equally to single-strand ρ.
    /// </remarks>
    /// <param name="ratio">Odds ratio ρ (finite or +∞, ≥ 0).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ratio"/> is NaN or negative.</exception>
    public static DinucleotideRepresentation ClassifyDinucleotideRatio(double ratio)
    {
        if (double.IsNaN(ratio) || ratio < 0)
            throw new ArgumentOutOfRangeException(nameof(ratio), ratio, "Relative abundance must be a non-negative number.");

        if (ratio <= KarlinUnderRepresentedThreshold) return DinucleotideRepresentation.UnderRepresented;
        if (ratio >= KarlinOverRepresentedThreshold) return DinucleotideRepresentation.OverRepresented;
        return DinucleotideRepresentation.Normal;
    }

    /// <summary>
    /// Karlin representation class of every dinucleotide returned by
    /// <see cref="CalculateDinucleotideRatios(string, bool)"/> (same keys; see
    /// <see cref="ClassifyDinucleotideRatio(double)"/>). A dinucleotide absent from the ratio map but
    /// whose two bases occur has ρ = 0 and is under-represented; it is not a key here, as in the ratio map.
    /// </summary>
    /// <param name="sequence">Nucleotide sequence (case-insensitive).</param>
    /// <param name="strandSymmetric">true for Karlin's strand-symmetrized ρ* (default); false for single-strand ρ.</param>
    public static IReadOnlyDictionary<string, DinucleotideRepresentation> ClassifyDinucleotideRatios(
        string sequence,
        bool strandSymmetric = true)
    {
        var classes = new Dictionary<string, DinucleotideRepresentation>();
        foreach (var (dinucleotide, ratio) in CalculateDinucleotideRatios(sequence, strandSymmetric))
            classes[dinucleotide] = ClassifyDinucleotideRatio(ratio);

        return classes;
    }

    /// <summary>
    /// Calculates codon usage frequencies by reading consecutive, non-overlapping triplets from the
    /// given reading frame: frequency = count(codon) / total counted codons. Input is case-insensitive
    /// and may be DNA (T) or RNA (U): U is read as T and codons are reported in DNA spelling (EMBOSS
    /// cusp, CodonW <c>ident_codon</c>; same contract as the canonical
    /// <c>CodonUsageAnalyzer.CountCodons(string)</c> in MolTools). Triplets containing any other symbol
    /// (N, IUPAC ambiguity codes, gaps, …) are excluded from both count and total without shifting the
    /// frame, and trailing 1-2 leftover bases are ignored (EMBOSS <c>ajCodSetTripletsS</c>). This is the
    /// count/total fraction used by the Kazusa Codon Usage Database (CUTG) and the EMBOSS cusp
    /// "/1000" column divided by 1000; it is distinct from the per-amino-acid "Fraction" column of cusp.
    /// Input shorter than 3 bases, or with no valid codon (total = 0), yields an empty table.
    /// Sources: Nakamura, Gojobori, Ikemura (2000), Nucleic Acids Res 28(1):292,
    /// DOI 10.1093/nar/28.1.292; Kazusa CUTG README, https://www.kazusa.or.jp/codon/readme_codon.html;
    /// EMBOSS 6.6.0 cusp / ajcod.c.
    /// </summary>
    /// <param name="dnaSequence">DNA or RNA coding sequence; case-insensitive, U read as T, triplets with other symbols excluded.</param>
    /// <param name="readingFrame">0-based offset of the first codon (0, 1, or 2 in practice; larger values are
    /// plain offsets, as EMBOSS compseq <c>-frame</c>). Must be non-negative.</param>
    /// <returns>Map of codon (DNA spelling) to its frequency (count / total counted codons); empty if no valid codon.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="readingFrame"/> is negative.</exception>
    public static IReadOnlyDictionary<string, double> CalculateCodonFrequencies(
        string dnaSequence,
        int readingFrame = 0)
    {
        // Canonical in-frame splitter (Core Translator.CountCodons, review 2026-09 R16): U read as T,
        // ASCII case fold, triplets with any non-ACGT symbol skipped without shifting the frame,
        // trailing partial triplet dropped, negative frame throws ArgumentOutOfRangeException.
        var counts = Translator.CountCodons(dnaSequence, readingFrame);
        int total = counts.Values.Sum();

        // total == 0 (no valid codon) leaves counts empty, so this yields an empty table and
        // avoids any division by zero.
        var freq = new Dictionary<string, double>(counts.Count);
        foreach (var (codon, count) in counts)
        {
            freq[codon] = (double)count / total;
        }

        return freq;
    }

    /// <summary>
    /// Normalisation of <see cref="CalculateCodonFrequencies(string, int, CodonFrequencyNormalization, GeneticCode?)"/>,
    /// matching the columns of EMBOSS 6.6.0 <c>cusp</c>.
    /// </summary>
    public enum CodonFrequencyNormalization
    {
        /// <summary>count / total counted codons (Kazusa CUTG; cusp "Frequency"/1000). Default.</summary>
        Total = 0,

        /// <summary>count / total counted codons × 1000 — cusp "Frequency" (per thousand) column.</summary>
        PerThousand = 1,

        /// <summary>
        /// count / Σ counts of the codons that encode the same amino acid (stop codons form one
        /// '*' class) under the genetic code — cusp "Fraction" column.
        /// </summary>
        PerAminoAcid = 2,
    }

    /// <summary>
    /// Codon usage with a selectable normalisation — the "Frequency" (per thousand) and per-amino-acid
    /// "Fraction" columns of EMBOSS 6.6.0 <c>cusp</c>.
    /// </summary>
    /// <remarks>
    /// <para>Codons are counted exactly as by <see cref="CalculateCodonFrequencies(string, int)"/> (canonical
    /// <c>Translator.CountCodons</c>: U read as T, case-insensitive, triplets with non-ACGT symbols skipped
    /// without shifting the frame, trailing partial triplet dropped). <see cref="CodonFrequencyNormalization.Total"/>
    /// equals that method; <see cref="CodonFrequencyNormalization.PerThousand"/> multiplies by 1000.</para>
    /// <para><see cref="CodonFrequencyNormalization.PerAminoAcid"/> divides each count by the total count of
    /// its synonymous family, the amino acid being <see cref="GeneticCode.Translate(string)"/> of the codon
    /// (stop codons share the '*' family). cusp always uses the standard code (NCBI table 1); other tables
    /// re-group the families (e.g. table 2: AGA/AGG join '*', TGA joins W). Example (cusp-executed):
    /// <c>…AAAAAGAAA…</c> → AAA 0.667, AAG 0.333; TAA/TAG/TGA 0.333 each.</para>
    /// <para>Only observed codons are keys (absent codon ⇒ 0, as cusp's zero rows). Difference from cusp:
    /// cusp maps a triplet with an ambiguous base to the lowest constituent (NNN → AAA) and counts it;
    /// here it is skipped (B02 F9 / <c>ajCodSetTripletsS</c> contract).</para>
    /// </remarks>
    /// <param name="dnaSequence">DNA or RNA coding sequence.</param>
    /// <param name="readingFrame">0-based offset of the first codon (non-negative).</param>
    /// <param name="normalization">Output normalisation.</param>
    /// <param name="geneticCode">Genetic code for <see cref="CodonFrequencyNormalization.PerAminoAcid"/>
    /// (default <see cref="GeneticCode.Standard"/>, as cusp); ignored by the other normalisations.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="readingFrame"/> is negative or
    /// <paramref name="normalization"/> is undefined.</exception>
    public static IReadOnlyDictionary<string, double> CalculateCodonFrequencies(
        string dnaSequence,
        int readingFrame,
        CodonFrequencyNormalization normalization,
        GeneticCode? geneticCode = null)
    {
        if (!Enum.IsDefined(normalization))
            throw new ArgumentOutOfRangeException(nameof(normalization), normalization, "Undefined codon normalisation.");

        if (normalization == CodonFrequencyNormalization.Total)
            return CalculateCodonFrequencies(dnaSequence, readingFrame);

        var counts = Translator.CountCodons(dnaSequence, readingFrame);
        var result = new Dictionary<string, double>(counts.Count);

        if (normalization == CodonFrequencyNormalization.PerThousand)
        {
            int total = counts.Values.Sum();
            foreach (var (codon, count) in counts)
                result[codon] = count * PerThousandScale / total;
            return result;
        }

        var code = geneticCode ?? GeneticCode.Standard;
        var familyTotals = new Dictionary<char, int>();
        var familyOf = new Dictionary<string, char>(counts.Count);
        foreach (var (codon, count) in counts)
        {
            char aminoAcid = code.Translate(codon);
            familyOf[codon] = aminoAcid;
            familyTotals[aminoAcid] = familyTotals.GetValueOrDefault(aminoAcid) + count;
        }

        foreach (var (codon, count) in counts)
            result[codon] = (double)count / familyTotals[familyOf[codon]];

        return result;
    }

    // cusp "Frequency" column scale: codons per thousand counted codons.
    private const double PerThousandScale = 1000.0;

    #endregion

    #region Entropy and Complexity

    /// <summary>
    /// Calculates the Shannon entropy H = −Σ pᵢ·log₂ pᵢ (bits per symbol; Shannon 1948) of the
    /// per-letter composition of a sequence.
    /// </summary>
    /// <remarks>
    /// Alphabet: every letter (<see cref="char.IsLetter(char)"/>) after upper-casing is its own symbol, so
    /// the method serves DNA, RNA and protein alike; non-letters (gaps, digits, '*', whitespace) are
    /// excluded from numerator and denominator. N and IUPAC codes are counted as distinct symbols and
    /// T and U are not merged — use <c>SequenceComplexity.CalculateShannonEntropy</c> for the
    /// nucleotide-only {A, C, G, T/U} alphabet (max 2 bits).
    /// The kernel is the canonical <see cref="StatisticsHelper.ShannonIndex"/> (natural log) converted to
    /// bits by dividing by ln 2 — the same computation as <c>scipy.stats.entropy(counts, base=2)</c>
    /// and scikit-bio <c>shannon(counts, base=2)</c>. Null, empty or letter-free input returns 0.
    /// </remarks>
    public static double CalculateShannonEntropy(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        var counts = new Dictionary<char, int>();
        foreach (char ch in sequence.ToUpperInvariant())
        {
            if (char.IsLetter(ch))
                counts[ch] = counts.GetValueOrDefault(ch) + 1;
        }

        if (counts.Count == 0)
            return 0;

        return StatisticsHelper.ShannonIndex(counts.Values.ToArray()) / Ln2;
    }

    // Base conversion ln → log₂ (bits): H₂ = H_e / ln 2 (scipy.stats.entropy divides by log(base)).
    private static readonly double Ln2 = Math.Log(2.0);

    /// <summary>
    /// Calculates linguistic complexity LC = Σ_{k=1..m} V_k / Σ_{k=1..m} min(a^k, N − k + 1)
    /// (V_k = distinct k-words, m = min(maxK, N), a = alphabet size) — Orlov &amp; Potapov (2004) summation form, equal to
    /// Troyanskaya et al. (2002) / Rosalind LING when m ≥ N. Delegates to the canonical
    /// <see cref="SequenceComplexity.CalculateLinguisticComplexity(string, int)"/>.
    /// </summary>
    /// <remarks>
    /// <paramref name="maxK"/> is the canonical maximum word length m (same meaning: word lengths
    /// 1..min(maxK, N)); only the default differs (6 here, 10 in <c>SequenceComplexity</c>).
    /// Before the 2026-09 review (B03 F21) this method returned the unsourced arithmetic mean of the
    /// per-k usages U_k (e.g. ATTTGGATT, m = 6: 293/336 instead of 29/34). Input is upper-cased; every character
    /// is a word symbol and a = |{A, C, G, T (U for RNA)} ∪ observed symbols| (inferred as in <c>SequenceComplexity</c>,
    /// B04 F20), so a ≥ 4 and LC ≤ 1 for any input (e.g. ACGTN: a = 5 → 15/15 = 1.0, formerly 15/14).
    /// Null/empty input or <paramref name="maxK"/> &lt; 1 returns 0.
    /// </remarks>
    public static double CalculateLinguisticComplexity(string sequence, int maxK = 6) =>
        SequenceComplexity.CalculateLinguisticComplexity(sequence, maxK);

    #endregion

    #region Protein Secondary Structure

    // Chou-Fasman (1978) conformational parameters and β-turn bend frequencies.
    // Pa/Pb/Pt are the published integer parameters (propensity × 100); f(i)..f(i+3) are the
    // positional bend frequencies of residues at the four positions of a β-turn tetrapeptide.
    // Values from Chou PY, Fasman GD (1978) "Empirical predictions of protein conformation"
    // Annu Rev Biochem 47:251-276 / Adv Enzymol 47:45-148, as reproduced in the
    // prowl.rockefeller.edu Chou-Fasman table (copied by residue NAME in
    // ravihansa3000/ChouFasman and by one-letter code in hassan11196/Chou-Fasman; both agree).
    // See docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md.
    private readonly record struct ChouFasmanParameters(
        int Pa, int Pb, int Pt, double F0, double F1, double F2, double F3);

    private static readonly Dictionary<char, ChouFasmanParameters> ChouFasmanTable = new()
    {
        ['A'] = new(142, 83, 66, 0.060, 0.076, 0.035, 0.058),
        ['R'] = new(98, 93, 95, 0.070, 0.106, 0.099, 0.085),
        ['N'] = new(67, 89, 156, 0.161, 0.083, 0.191, 0.091),
        ['D'] = new(101, 54, 146, 0.147, 0.110, 0.179, 0.081),
        ['C'] = new(70, 119, 119, 0.149, 0.050, 0.117, 0.128),
        ['E'] = new(151, 37, 74, 0.056, 0.060, 0.077, 0.064),
        ['Q'] = new(111, 110, 98, 0.074, 0.098, 0.037, 0.098),
        ['G'] = new(57, 75, 156, 0.102, 0.085, 0.190, 0.152),
        ['H'] = new(100, 87, 95, 0.140, 0.047, 0.093, 0.054),
        ['I'] = new(108, 160, 47, 0.043, 0.034, 0.013, 0.056),
        ['L'] = new(121, 130, 59, 0.061, 0.025, 0.036, 0.070),
        ['K'] = new(114, 74, 101, 0.055, 0.115, 0.072, 0.095),
        ['M'] = new(145, 105, 60, 0.068, 0.082, 0.014, 0.055),
        ['F'] = new(113, 138, 60, 0.059, 0.041, 0.065, 0.065),
        ['P'] = new(57, 55, 152, 0.102, 0.301, 0.034, 0.068),
        ['S'] = new(77, 75, 143, 0.120, 0.139, 0.125, 0.106),
        ['T'] = new(83, 119, 96, 0.086, 0.108, 0.065, 0.079),
        ['W'] = new(108, 137, 96, 0.077, 0.013, 0.064, 0.167),
        ['Y'] = new(69, 147, 114, 0.082, 0.065, 0.114, 0.125),
        ['V'] = new(106, 170, 50, 0.062, 0.048, 0.028, 0.053),
    };

    // Integer parameters are propensity × 100; 142 / 100.0 is bit-identical to the literal 1.42.
    private const double ChouFasmanScale = 100.0;

    private static readonly Dictionary<char, (double Helix, double Sheet, double Turn)> SecondaryStructurePropensity =
        ChouFasmanTable.ToDictionary(
            kv => kv.Key,
            kv => (kv.Value.Pa / ChouFasmanScale, kv.Value.Pb / ChouFasmanScale, kv.Value.Pt / ChouFasmanScale));

    // Default sliding-window length. Chou & Fasman (1978) scan a hexapeptide window for
    // helix nucleation (4 of 6) and a pentapeptide window for sheet nucleation (3 of 5);
    // this profile uses a single configurable window whose default is the helix window + 1
    // residue of context. The window length is a caller parameter, not a Chou-Fasman constant.
    private const int DefaultWindowSize = 7;

    /// <summary>
    /// Computes per-window mean Chou-Fasman conformational propensities (helix Pa,
    /// sheet Pb, turn Pt) along a protein sequence. Each emitted tuple is the average
    /// propensity over a sliding window of <paramref name="windowSize"/> residues, stepping
    /// one residue at a time. A window mean &gt; 1.0 indicates a residue stretch that favours
    /// the corresponding conformation. Unknown residues (e.g. X, B, Z, gaps) are skipped and
    /// excluded from the per-window average.
    /// </summary>
    /// <param name="proteinSequence">Amino-acid sequence in one-letter code; case-insensitive.</param>
    /// <param name="windowSize">Sliding-window length (residues). Must be positive and
    /// not exceed the sequence length.</param>
    /// <returns>One (Helix, Sheet, Turn) mean-propensity tuple per window position, in
    /// order from the N-terminus. Empty when the sequence is null/empty, the window is
    /// larger than the sequence, or the window is non-positive.</returns>
    public static IEnumerable<(double Helix, double Sheet, double Turn)> PredictSecondaryStructure(
        string proteinSequence,
        int windowSize = DefaultWindowSize)
    {
        if (string.IsNullOrEmpty(proteinSequence) || windowSize < 1 || windowSize > proteinSequence.Length)
            yield break;

        string upper = proteinSequence.ToUpperInvariant();

        for (int i = 0; i <= upper.Length - windowSize; i++)
        {
            double helixSum = 0, sheetSum = 0, turnSum = 0;
            int count = 0;

            for (int j = 0; j < windowSize; j++)
            {
                if (SecondaryStructurePropensity.TryGetValue(upper[i + j], out var prop))
                {
                    helixSum += prop.Helix;
                    sheetSum += prop.Sheet;
                    turnSum += prop.Turn;
                    count++;
                }
            }

            if (count > 0)
            {
                yield return (helixSum / count, sheetSum / count, turnSum / count);
            }
        }
    }

    // Chou & Fasman (1978) assignment rules, as stated by Chen, Gu & Huang (2006) BMC
    // Bioinformatics 7(Suppl 4):S14 "Methods" rules 1-3, and the β-turn rule (p(t) > 7.5e-5,
    // <Pt> > 1.00, <Pa> < <Pt> > <Pb>) of the 1978 method.
    private const int HelixNucleationWindow = 6;    // 4 of 6 helix formers
    private const int HelixNucleationFormers = 4;
    private const int SheetNucleationWindow = 5;    // 3 of 5 sheet formers
    private const int SheetNucleationFormers = 3;
    private const int FormerThreshold = 100;        // former: P > 1.00
    private const int ExtensionTetrapeptide = 4;    // extend until tetrapeptide <P> < 1.00
    private const int ExtensionMinimumSum = 400;    // 4 × 1.00 (integer units)
    private const int HelixAcceptThreshold = 103;   // segment <Pa> > 1.03
    private const int SheetAcceptThreshold = 105;   // segment <Pb> > 1.05
    private const int TurnLength = 4;
    private const int TurnMinimumPtSum = 400;       // tetrapeptide <Pt> > 1.00
    private const double TurnBendProbabilityThreshold = 7.5e-5; // p(t) = f(i)f(i+1)f(i+2)f(i+3)

    /// <summary>
    /// Nucleation rule set of <see cref="PredictSecondaryStructureChouFasman(string, ChouFasmanRuleSet)"/>.
    /// </summary>
    public enum ChouFasmanRuleSet
    {
        /// <summary>Default: a residue is a former when its parameter exceeds 1.00 (Pα &gt; 1.00 /
        /// Pβ &gt; 1.00); nucleus = ≥ 4 formers in 6 (helix) / ≥ 3 in 5 (sheet); breakers are not counted
        /// (Chen, Gu &amp; Huang 2006 statement of the rules). Same output as the one-argument overload.</summary>
        Standard = 0,

        /// <summary>Chou &amp; Fasman (1978) graded classes for the 1978 parameters (Prevelige &amp; Fasman 1989
        /// Table II): helix Hα = E M A L, hα = K F Q W I V, Iα = D H, iα = R T S C, bα = Y N, Bα = P G;
        /// sheet Hβ = V I Y, hβ = F W L C T Q M, iβ = R N H A, bβ = S G K, Bβ = P D E (no Iβ residue).
        /// Nucleus: Hα/hα (Hβ/hβ) count 1, Iα (Iβ) ½; helix ≥ 4 of 6, sheet ≥ 3 of 5, and at most one
        /// breaker (b or B; ≥ ⅓ breakers is unfavourable). Boundary rules: Pro may occur only in the first three
        /// helix residues (it is never added at the C-terminal end or inside the helix) and counts as Iα in the
        /// N-terminal turn of a helix nucleus; Arg counts as Iα in the C-terminal three residues of a helix
        /// nucleus; Trp counts as bβ in the C-terminal two residues of a sheet nucleus. Extension (until a
        /// tetrapeptide mean &lt; 1.00), acceptance (⟨Pα⟩ &gt; 1.03 / ⟨Pβ⟩ &gt; 1.05 and greater than the
        /// competitor), overlap and β-turn steps are those of <see cref="Standard"/>.</summary>
        Graded1978 = 1,
    }

    /// <summary>Per-residue state code: α-helix.</summary>
    public const char ChouFasmanHelix = 'H';
    /// <summary>Per-residue state code: β-strand (sheet).</summary>
    public const char ChouFasmanSheet = 'E';
    /// <summary>Per-residue state code: β-turn.</summary>
    public const char ChouFasmanTurn = 'T';
    /// <summary>Per-residue state code: coil (no assignment).</summary>
    public const char ChouFasmanCoil = 'C';

    /// <summary>
    /// Assigns a discrete secondary-structure state to every residue with the Chou &amp; Fasman
    /// (1978) prediction rules, using the 1978 conformational parameters (Pα, Pβ, Pt) and
    /// β-turn bend frequencies f(i)..f(i+3):
    /// <list type="number">
    /// <item><description>Helix nucleation: any 6-residue window with ≥ 4 helix formers (Pα &gt; 1.00).
    /// Sheet nucleation: any 5-residue window with ≥ 3 sheet formers (Pβ &gt; 1.00).</description></item>
    /// <item><description>Extension: each nucleus is extended residue by residue in both directions
    /// while the tetrapeptide formed by the new residue and the three adjacent segment residues has
    /// mean propensity ≥ 1.00 (extension stops when it drops below 1.00).</description></item>
    /// <item><description>Acceptance: an extended helix is kept if ⟨Pα⟩ &gt; 1.03 and ⟨Pα⟩ &gt; ⟨Pβ⟩
    /// over the segment; an extended strand if ⟨Pβ⟩ &gt; 1.05 and ⟨Pβ⟩ &gt; ⟨Pα⟩.</description></item>
    /// <item><description>Overlap: each maximal run of residues covered by both a helix and a strand
    /// is assigned helix if ⟨Pα⟩ &gt; ⟨Pβ⟩ over the run, otherwise strand.</description></item>
    /// <item><description>β-turn: a tetrapeptide i..i+3 is a turn if
    /// p(t) = f(i)·f(i+1)·f(i+2)·f(i+3) &gt; 7.5×10⁻⁵, ⟨Pt⟩ &gt; 1.00 and ⟨Pα⟩ &lt; ⟨Pt⟩ &gt; ⟨Pβ⟩;
    /// all four residues are marked turn and take precedence over helix/strand.</description></item>
    /// </list>
    /// Residues other than the 20 standard amino acids have no parameters: they are never
    /// formers, windows/tetrapeptides containing them are ineligible, and they stay coil.
    /// Unlike <see cref="PredictSecondaryStructure(string, int)"/> (a windowed mean-propensity
    /// profile), this is the discrete Chou-Fasman assignment. Q3 accuracy is ~50-60%.
    /// </summary>
    /// <param name="proteinSequence">Amino-acid sequence in one-letter code; case-insensitive.</param>
    /// <returns>A string of the same length as the input with one state code per residue:
    /// <see cref="ChouFasmanHelix"/> ('H'), <see cref="ChouFasmanSheet"/> ('E'),
    /// <see cref="ChouFasmanTurn"/> ('T') or <see cref="ChouFasmanCoil"/> ('C').
    /// Empty for null/empty input.</returns>
    public static string PredictSecondaryStructureChouFasman(string proteinSequence) =>
        PredictSecondaryStructureChouFasman(proteinSequence, ChouFasmanRuleSet.Standard);

    /// <summary>
    /// Discrete Chou &amp; Fasman assignment with a selectable nucleation rule set.
    /// <see cref="ChouFasmanRuleSet.Standard"/> is identical to
    /// <see cref="PredictSecondaryStructureChouFasman(string)"/>;
    /// <see cref="ChouFasmanRuleSet.Graded1978"/> applies the graded former/breaker classes
    /// (Hα/hα/Iα/iα/bα/Bα, Hβ/hβ/Iβ/iβ/bβ/Bβ) of Chou &amp; Fasman (1978), fractional former counting,
    /// the breaker limit, the proline rule and the helix/sheet boundary reassignments — see
    /// <see cref="ChouFasmanRuleSet.Graded1978"/>. Extension, acceptance, overlap and β-turn steps are shared.
    /// </summary>
    /// <param name="proteinSequence">Amino-acid sequence in one-letter code; case-insensitive.</param>
    /// <param name="ruleSet">Nucleation rule set (default behaviour = <see cref="ChouFasmanRuleSet.Standard"/>).</param>
    /// <returns>One state code (H/E/T/C) per residue; empty for null/empty input.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="ruleSet"/> is not a defined value.</exception>
    public static string PredictSecondaryStructureChouFasman(string proteinSequence, ChouFasmanRuleSet ruleSet)
    {
        if (ruleSet is not (ChouFasmanRuleSet.Standard or ChouFasmanRuleSet.Graded1978))
            throw new ArgumentOutOfRangeException(nameof(ruleSet), ruleSet, "Unknown Chou-Fasman rule set.");

        if (string.IsNullOrEmpty(proteinSequence))
            return string.Empty;

        bool graded = ruleSet == ChouFasmanRuleSet.Graded1978;

        int n = proteinSequence.Length;
        var parameters = new ChouFasmanParameters?[n];
        for (int i = 0; i < n; i++)
        {
            if (ChouFasmanTable.TryGetValue(char.ToUpperInvariant(proteinSequence[i]), out var p))
                parameters[i] = p;
        }

        bool[] helix, sheet;
        if (graded)
        {
            string upper = proteinSequence.ToUpperInvariant();
            helix = FindChouFasmanSegments(parameters, static p => p.Pa, static p => p.Pb,
                HelixNucleationWindow, HelixNucleationFormers, HelixAcceptThreshold,
                (start, window) => IsGradedHelixNucleus(upper, start, window), upper);
            sheet = FindChouFasmanSegments(parameters, static p => p.Pb, static p => p.Pa,
                SheetNucleationWindow, SheetNucleationFormers, SheetAcceptThreshold,
                (start, window) => IsGradedSheetNucleus(upper, start, window), null);
        }
        else
        {
            helix = FindChouFasmanSegments(parameters, static p => p.Pa, static p => p.Pb,
                HelixNucleationWindow, HelixNucleationFormers, HelixAcceptThreshold, null, null);
            sheet = FindChouFasmanSegments(parameters, static p => p.Pb, static p => p.Pa,
                SheetNucleationWindow, SheetNucleationFormers, SheetAcceptThreshold, null, null);
        }

        var states = new char[n];
        Array.Fill(states, ChouFasmanCoil);

        for (int i = 0; i < n;)
        {
            if (helix[i] && sheet[i])
            {
                int j = i;
                while (j + 1 < n && helix[j + 1] && sheet[j + 1])
                    j++;

                long sumPa = 0, sumPb = 0;
                for (int k = i; k <= j; k++)
                {
                    sumPa += parameters[k]!.Value.Pa;
                    sumPb += parameters[k]!.Value.Pb;
                }

                char winner = sumPa > sumPb ? ChouFasmanHelix : ChouFasmanSheet;
                for (int k = i; k <= j; k++)
                    states[k] = winner;
                i = j + 1;
            }
            else
            {
                if (helix[i])
                    states[i] = ChouFasmanHelix;
                else if (sheet[i])
                    states[i] = ChouFasmanSheet;
                i++;
            }
        }

        for (int i = 0; i + TurnLength <= n; i++)
        {
            if (IsChouFasmanTurn(parameters, i))
            {
                for (int k = i; k < i + TurnLength; k++)
                    states[k] = ChouFasmanTurn;
            }
        }

        return new string(states);
    }

    private static bool[] FindChouFasmanSegments(
        ChouFasmanParameters?[] parameters,
        Func<ChouFasmanParameters, int> own,
        Func<ChouFasmanParameters, int> competitor,
        int window,
        int minimumFormers,
        int acceptThreshold,
        Func<int, int, bool>? gradedNucleus,
        string? prolineRuleSequence)
    {
        int n = parameters.Length;
        var covered = new bool[n];

        for (int start = 0; start + window <= n; start++)
        {
            int formers = 0;
            bool allKnown = true;
            for (int k = start; k < start + window; k++)
            {
                if (parameters[k] is not { } p) { allKnown = false; break; }
                if (own(p) > FormerThreshold) formers++;
            }
            if (!allKnown)
                continue;
            if (gradedNucleus is null ? formers < minimumFormers : !gradedNucleus(start, window))
                continue;

            int first = start, last = start + window - 1;

            // C-terminal extension: tetrapeptide = last three segment residues + the new residue.
            // Proline rule (graded helix only): Pro is never added at the C-terminal end.
            while (last + 1 < n && TetrapeptideSum(parameters, last - 2, own) >= ExtensionMinimumSum
                   && !(prolineRuleSequence is not null && prolineRuleSequence[last + 1] == 'P'))
                last++;
            // N-terminal extension: tetrapeptide = the new residue + first three segment residues.
            // Proline rule (graded helix only): every Pro must stay within the first three helix residues.
            while (first - 1 >= 0 && TetrapeptideSum(parameters, first - 1, own) >= ExtensionMinimumSum
                   && (prolineRuleSequence is null || ProlinesWithinNTerminalTurn(prolineRuleSequence, first - 1, last)))
                first--;

            long sumOwn = 0, sumCompetitor = 0;
            for (int k = first; k <= last; k++)
            {
                sumOwn += own(parameters[k]!.Value);
                sumCompetitor += competitor(parameters[k]!.Value);
            }

            int length = last - first + 1;
            if (sumOwn > (long)acceptThreshold * length && sumOwn > sumCompetitor)
            {
                for (int k = first; k <= last; k++)
                    covered[k] = true;
            }
        }

        return covered;
    }

    // Sum of a parameter over the tetrapeptide starting at 'from'; int.MinValue if any residue
    // lacks parameters (so it can never satisfy an extension threshold).
    private static int TetrapeptideSum(
        ChouFasmanParameters?[] parameters, int from, Func<ChouFasmanParameters, int> selector)
    {
        int sum = 0;
        for (int k = from; k < from + ExtensionTetrapeptide; k++)
        {
            if (parameters[k] is not { } p)
                return int.MinValue;
            sum += selector(p);
        }
        return sum;
    }

    // ---- Graded 1978 rule set (ChouFasmanRuleSet.Graded1978) ----
    // Class assignments for the 1978 (29-protein) parameters: Chou & Fasman (1978) as reproduced in
    // Prevelige & Fasman (1989) Table II (brackets Hα/hα/Iα/iα/bα/Bα and Hβ/hβ/iβ/bβ/Bβ; no residue is Iβ).
    private enum ChouFasmanClass { StrongFormer, Former, WeakFormer, Indifferent, Breaker, StrongBreaker }

    private const string HelixStrongFormers = "EMAL";     // Hα
    private const string HelixFormers = "KFQWIV";         // hα
    private const string HelixWeakFormers = "DH";         // Iα
    private const string HelixBreakers = "YN";            // bα
    private const string HelixStrongBreakers = "PG";      // Bα   (iα = R T S C)
    private const string SheetStrongFormers = "VIY";      // Hβ
    private const string SheetFormers = "FWLCTQM";        // hβ
    private const string SheetBreakers = "SGK";           // bβ
    private const string SheetStrongBreakers = "PDE";     // Bβ   (iβ = R N H A; Iβ empty)

    // Fractional counting in half-former units: Hα/hα (Hβ/hβ) = 2, Iα (Iβ) = 1 ("Iα counts as ½ hα").
    private const int HalfUnitsPerFormer = 2;
    private const int MaximumNucleusBreakers = 1;         // ≥ 1/3 breakers (bα/Bα, bβ/Bβ) is unfavourable
    // Boundary zone of a nucleation window: residues "near the N-terminal end" are window offsets
    // 0..2 (the first helical turn, where Pro is allowed); "near the C-terminal end" are the last three
    // offsets of the helix window and the last two of the β window.
    private const int HelixBoundaryZone = 3;
    private const int SheetBoundaryZone = 2;

    private static ChouFasmanClass HelixClass(char residue) => residue switch
    {
        _ when HelixStrongFormers.Contains(residue) => ChouFasmanClass.StrongFormer,
        _ when HelixFormers.Contains(residue) => ChouFasmanClass.Former,
        _ when HelixWeakFormers.Contains(residue) => ChouFasmanClass.WeakFormer,
        _ when HelixBreakers.Contains(residue) => ChouFasmanClass.Breaker,
        _ when HelixStrongBreakers.Contains(residue) => ChouFasmanClass.StrongBreaker,
        _ => ChouFasmanClass.Indifferent,
    };

    private static ChouFasmanClass SheetClass(char residue) => residue switch
    {
        _ when SheetStrongFormers.Contains(residue) => ChouFasmanClass.StrongFormer,
        _ when SheetFormers.Contains(residue) => ChouFasmanClass.Former,
        _ when SheetBreakers.Contains(residue) => ChouFasmanClass.Breaker,
        _ when SheetStrongBreakers.Contains(residue) => ChouFasmanClass.StrongBreaker,
        _ => ChouFasmanClass.Indifferent,
    };

    private static bool IsGradedNucleus(IEnumerable<ChouFasmanClass> classes, int requiredFormers)
    {
        int halfUnits = 0, breakers = 0;
        foreach (var c in classes)
        {
            halfUnits += c switch
            {
                ChouFasmanClass.StrongFormer or ChouFasmanClass.Former => HalfUnitsPerFormer,
                ChouFasmanClass.WeakFormer => 1,
                _ => 0,
            };
            if (c is ChouFasmanClass.Breaker or ChouFasmanClass.StrongBreaker)
                breakers++;
        }
        return halfUnits >= requiredFormers * HalfUnitsPerFormer && breakers <= MaximumNucleusBreakers;
    }

    // Helix nucleus: Pro may only occupy the N-terminal turn (offsets 0..2) and is then assigned Iα;
    // Arg near the C-terminal end (last three offsets) is assigned Iα (Table I footnote b).
    private static bool IsGradedHelixNucleus(string upper, int start, int window)
    {
        var classes = new ChouFasmanClass[window];
        for (int k = 0; k < window; k++)
        {
            char r = upper[start + k];
            bool nTerminal = k < HelixBoundaryZone;
            bool cTerminal = k >= window - HelixBoundaryZone;
            if (r == 'P' && !nTerminal)
                return false;
            classes[k] = (r == 'P' && nTerminal) || (r == 'R' && cTerminal)
                ? ChouFasmanClass.WeakFormer
                : HelixClass(r);
        }
        return IsGradedNucleus(classes, HelixNucleationFormers);
    }

    // Sheet nucleus: Trp near the C-terminal end (last two offsets) is assigned bβ (Table I footnote c).
    private static bool IsGradedSheetNucleus(string upper, int start, int window)
    {
        var classes = new ChouFasmanClass[window];
        for (int k = 0; k < window; k++)
        {
            char r = upper[start + k];
            classes[k] = r == 'W' && k >= window - SheetBoundaryZone ? ChouFasmanClass.Breaker : SheetClass(r);
        }
        return IsGradedNucleus(classes, SheetNucleationFormers);
    }

    // True if every Pro in [first, last] lies within the first three residues of the segment.
    private static bool ProlinesWithinNTerminalTurn(string upper, int first, int last)
    {
        for (int k = first + HelixBoundaryZone; k <= last; k++)
        {
            if (upper[k] == 'P')
                return false;
        }
        return true;
    }

    private static bool IsChouFasmanTurn(ChouFasmanParameters?[] parameters, int i)
    {
        if (parameters[i] is not { } a || parameters[i + 1] is not { } b ||
            parameters[i + 2] is not { } c || parameters[i + 3] is not { } d)
            return false;

        double bendProbability = a.F0 * b.F1 * c.F2 * d.F3;
        int sumPt = a.Pt + b.Pt + c.Pt + d.Pt;
        int sumPa = a.Pa + b.Pa + c.Pa + d.Pa;
        int sumPb = a.Pb + b.Pb + c.Pb + d.Pb;

        return bendProbability > TurnBendProbabilityThreshold
            && sumPt > TurnMinimumPtSum
            && sumPt > sumPa
            && sumPt > sumPb;
    }

    #endregion

    #region Sequence Windows

    // GC content is expressed as a percentage per the GC-content definition
    // GC% = (G + C) / (A + T + G + C) × 100  [Wikipedia GC-content; Biopython gc_fraction ×100].
    private const double PercentScale = 100.0;
    // Default sliding-window width for GC profiling (bp); a window of 100 is a common
    // local-GC resolution and matches the sibling profile methods in this class.
    private const int DefaultGcProfileWindow = 100;

    /// <summary>
    /// Calculates the GC-content profile of a sequence: the GC content of each sliding
    /// window of length <paramref name="windowSize"/> advanced by <paramref name="stepSize"/>,
    /// expressed as a percentage GC% = (G + C) / (A + T + G + C) × 100.
    /// </summary>
    /// <param name="sequence">Input nucleotide sequence (DNA or RNA; case-insensitive).</param>
    /// <param name="windowSize">Window width W in bases (≥ 1; default 100). Must be ≤ sequence length for any window to be produced.</param>
    /// <param name="stepSize">Window advance in bases (≥ 1; default 1).</param>
    /// <param name="fraction">true → report a fraction in [0, 1] (Biopython <c>gc_fraction</c>); false (default) → percentage.</param>
    /// <returns>
    /// One GC% value per window position, in order, for offsets 0, stepSize, 2·stepSize, …
    /// up to (length − windowSize); empty when the sequence is null/empty or
    /// <paramref name="windowSize"/> exceeds its length. The denominator counts only the
    /// standard bases A/T/U/G/C in the window (ambiguous symbols such as N are excluded,
    /// matching Biopython gc_fraction's default <c>ambiguous="remove"</c>); a window with no
    /// standard base yields 0.
    /// </returns>
    /// <remarks>
    /// GC-content definition: (G + C) / (A + T + G + C) × 100 — Wikipedia, GC-content
    /// (citing primary literature); Biopython <c>Bio.SeqUtils.gc_fraction</c> returns the
    /// same quantity as a fraction in [0, 1] (×100 here). U is treated as a non-GC base
    /// equivalent to T. Per-window value = canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char})"/> (S/W excluded like
    /// other ambiguity codes). Only complete windows are reported (as EMBOSS <c>isochore</c> and
    /// <see cref="GcSkewCalculator.CalculateWindowedGcSkew(string,int,int)"/>; Biopython
    /// <c>GC_skew</c> instead appends a trailing partial window — opt in via
    /// <see cref="CalculateGcContentProfile(string,int,int,bool,bool)"/>). Values carry no positions: window
    /// k starts at 0-based offset k·stepSize.
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly).</exception>
    public static IEnumerable<double> CalculateGcContentProfile(
        string sequence,
        int windowSize = DefaultGcProfileWindow,
        int stepSize = 1,
        bool fraction = false)
    {
        ValidateWindowAndStep(windowSize, stepSize);

        if (string.IsNullOrEmpty(sequence) || windowSize > sequence.Length)
            return Array.Empty<double>();

        // Canonical sliding-GC driver (B01): same complete-window loop and per-window
        // SequenceExtensions.CalculateGcFraction, ×100 unless fraction (2026-10 FIN-B03 F24).
        return GcSkewCalculator.CalculateWindowedGcContent(sequence, windowSize, stepSize, fraction)
            .Select(p => p.GcContent);
    }

    /// <summary>
    /// GC-content profile with Biopython <c>gc_fraction(window, ambiguous=…)</c> IUPAC-ambiguity handling:
    /// each complete window is scored by the canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char},SequenceExtensions.GcAmbiguityMode)"/>
    /// (<c>Remove</c>: S counts as GC, S/W in the denominator, other codes excluded; <c>Ignore</c>:
    /// denominator = window length; <c>Weighted</c>: ambiguity codes add their mean GC, e.g. N = 0.5).
    /// Window/step semantics are those of <see cref="CalculateGcContentProfile(string,int,int,bool)"/>.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1.</exception>
    public static IEnumerable<double> CalculateGcContentProfile(
        string sequence,
        int windowSize,
        int stepSize,
        bool fraction,
        SequenceExtensions.GcAmbiguityMode ambiguityMode)
    {
        ValidateWindowAndStep(windowSize, stepSize);

        if (string.IsNullOrEmpty(sequence) || windowSize > sequence.Length)
            return Array.Empty<double>();

        return GcSkewCalculator.CalculateWindowedGcContent(sequence, windowSize, stepSize, fraction, ambiguityMode)
            .Select(p => p.GcContent);
    }

    /// <summary>
    /// GC-content profile with an opt-in trailing partial window, Biopython <c>GC_skew</c> windowing:
    /// window starts are <c>range(0, len(seq), stepSize)</c> and each window is the slice
    /// <c>seq[i:i+windowSize]</c>, truncated at the sequence end.
    /// </summary>
    /// <param name="sequence">Input nucleotide sequence (DNA or RNA; case-insensitive).</param>
    /// <param name="windowSize">Window width W in bases (≥ 1).</param>
    /// <param name="stepSize">Window advance in bases (≥ 1).</param>
    /// <param name="fraction">true → fraction in [0, 1]; false → percentage.</param>
    /// <param name="includePartialWindow">false → identical to
    /// <see cref="CalculateGcContentProfile(string,int,int,bool)"/> (complete windows only).
    /// true → after the complete windows, every further start &lt; length is emitted with the
    /// truncated window (so a sequence shorter than W yields its partial windows instead of nothing).</param>
    /// <remarks>
    /// Biopython 1.88 has no windowed GC-content function; its only windowed nucleotide statistic,
    /// <c>Bio.SeqUtils.GC_skew(seq, window)</c>, iterates <c>for i in range(0, len(seq), window): s = seq[i:i+window]</c>
    /// and so emits the trailing partial window. This option applies that convention (generalised to
    /// step ≠ W) to <c>gc_fraction</c>: oracle <c>[gc_fraction(s[i:i+w]) for i in range(0, len(s), step)]</c>.
    /// Complete windows come from the canonical <see cref="GcSkewCalculator.CalculateWindowedGcContent(string,int,int,bool)"/>;
    /// each partial window is scored by the canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char})"/> (0 for a window without A/C/G/T/U).
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly).</exception>
    public static IEnumerable<double> CalculateGcContentProfile(
        string sequence,
        int windowSize,
        int stepSize,
        bool fraction,
        bool includePartialWindow)
    {
        if (!includePartialWindow)
            return CalculateGcContentProfile(sequence, windowSize, stepSize, fraction);

        ValidateWindowAndStep(windowSize, stepSize);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<double>();

        double scale = fraction ? 1.0 : PercentScale;
        return CalculateGcContentProfile(sequence, windowSize, stepSize, fraction)
            .Concat(PartialWindowStarts(sequence.Length, windowSize, stepSize)
                .Select(i => sequence.AsSpan(i).CalculateGcFraction() * scale));
    }

    /// <summary>
    /// <see cref="CalculateGcContentProfile(string,int,int,bool,bool)"/> with Biopython
    /// <c>gc_fraction(window, ambiguous=…)</c> handling (see
    /// <see cref="CalculateGcContentProfile(string,int,int,bool,SequenceExtensions.GcAmbiguityMode)"/>);
    /// partial windows are scored by the canonical
    /// <see cref="SequenceExtensions.CalculateGcFraction(ReadOnlySpan{char},SequenceExtensions.GcAmbiguityMode)"/>.
    /// </summary>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1.</exception>
    public static IEnumerable<double> CalculateGcContentProfile(
        string sequence,
        int windowSize,
        int stepSize,
        bool fraction,
        SequenceExtensions.GcAmbiguityMode ambiguityMode,
        bool includePartialWindow)
    {
        if (!includePartialWindow)
            return CalculateGcContentProfile(sequence, windowSize, stepSize, fraction, ambiguityMode);

        ValidateWindowAndStep(windowSize, stepSize);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<double>();

        double scale = fraction ? 1.0 : PercentScale;
        return CalculateGcContentProfile(sequence, windowSize, stepSize, fraction, ambiguityMode)
            .Concat(PartialWindowStarts(sequence.Length, windowSize, stepSize)
                .Select(i => sequence.AsSpan(i).CalculateGcFraction(ambiguityMode) * scale));
    }

    // Starts k·step (k ≥ 0) with start < length but start + windowSize > length — the truncated
    // windows that Biopython's range(0, len, step) slicing adds after the complete ones. long
    // arithmetic avoids int overflow for huge steps/windows.
    private static IEnumerable<int> PartialWindowStarts(int length, int windowSize, int stepSize)
    {
        long first = length >= windowSize
            ? ((long)(length - windowSize) / stepSize + 1) * stepSize
            : 0;
        for (long i = first; i < length; i += stepSize)
            yield return (int)i;
    }

    // Shared eager validation for the sliding-window profiles (as GcSkewCalculator.CalculateWindowedGcSkew):
    // a window below 1 has no content and a step below 1 never advances (infinite loop).
    private static void ValidateWindowAndStep(int windowSize, int stepSize)
    {
        if (windowSize < 1)
            throw new ArgumentOutOfRangeException(nameof(windowSize), windowSize, "Window size must be at least 1.");
        if (stepSize < 1)
            throw new ArgumentOutOfRangeException(nameof(stepSize), stepSize, "Step size must be at least 1.");
    }

    /// <summary>
    /// Calculates the Shannon entropy profile of a sequence: the per-symbol Shannon
    /// entropy H = -Σ pᵢ log₂ pᵢ (in bits) of each sliding window of length
    /// <paramref name="windowSize"/>, advanced by <paramref name="stepSize"/>.
    /// Per-window entropy is delegated to <see cref="CalculateShannonEntropy"/>.
    /// </summary>
    /// <param name="sequence">Input sequence; symbol frequencies are taken over its letters (case-folded).</param>
    /// <param name="windowSize">Window width W in symbols (≥ 1; default 50). Must be ≤ sequence length for any window to be produced.</param>
    /// <param name="stepSize">Window advance in symbols (≥ 1; default 1).</param>
    /// <returns>
    /// One entropy value (bits) per window position, in order, for offsets
    /// 0, stepSize, 2·stepSize, … up to (length − windowSize); empty when the
    /// sequence is null/empty or <paramref name="windowSize"/> exceeds its length.
    /// </returns>
    /// <remarks>
    /// Shannon C. E. (1948), A Mathematical Theory of Communication, Bell Syst. Tech. J.
    /// 27(3):379–423. Base-2 logarithm yields bits; maximum is log₂k for k distinct symbols
    /// (2 bits for the 4-letter DNA alphabet). Windows are taken over raw characters; only
    /// complete windows are reported; window k starts at 0-based offset k·stepSize.
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly).</exception>
    public static IEnumerable<double> CalculateEntropyProfile(
        string sequence,
        int windowSize = 50,
        int stepSize = 1)
    {
        ValidateWindowAndStep(windowSize, stepSize);

        if (string.IsNullOrEmpty(sequence) || windowSize > sequence.Length)
            return Array.Empty<double>();

        return EntropyProfileIterator(sequence, windowSize, stepSize);
    }

    private static IEnumerable<double> EntropyProfileIterator(string sequence, int windowSize, int stepSize)
    {
        for (int i = 0; i <= sequence.Length - windowSize; i += stepSize)
            yield return CalculateShannonEntropy(sequence.Substring(i, windowSize));
    }

    /// <summary>
    /// Entropy profile with optional Biopython-style trailing partial windows (letter-symbol entropy as
    /// <see cref="CalculateEntropyProfile(string, int, int)"/>).
    /// </summary>
    /// <remarks>
    /// With <paramref name="includePartialWindow"/> = true the window starts are
    /// <c>range(0, len(seq), stepSize)</c> as in Biopython 1.88 <c>GC_skew</c> slicing: after the complete
    /// windows, every start with start + W &gt; length yields the entropy of the truncated window
    /// <c>seq[start:]</c> — the same convention as
    /// <see cref="CalculateGcContentProfile(string, int, int, bool, bool)"/> (B03 F24). false equals the
    /// three-argument overload. Oracle: <c>scipy.stats.entropy(letter counts of seq[i:i+W], base=2)</c>.
    /// A window with no letter yields 0 (scipy: nan).
    /// </remarks>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/> or
    /// <paramref name="stepSize"/> is less than 1 (validated eagerly).</exception>
    public static IEnumerable<double> CalculateEntropyProfile(
        string sequence,
        int windowSize,
        int stepSize,
        bool includePartialWindow)
    {
        if (!includePartialWindow)
            return CalculateEntropyProfile(sequence, windowSize, stepSize);

        ValidateWindowAndStep(windowSize, stepSize);
        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<double>();

        return CalculateEntropyProfile(sequence, windowSize, stepSize)
            .Concat(PartialWindowStarts(sequence.Length, windowSize, stepSize)
                .Select(i => CalculateShannonEntropy(sequence.Substring(i))));
    }

    /// <summary>
    /// Higher-order (k-mer, block) entropy profile: the Shannon entropy in bits of the overlapping k-mer
    /// distribution of each window, delegated per window to the canonical
    /// <see cref="SequenceComplexity.CalculateKmerEntropy(string, int)"/> (no new kernel).
    /// </summary>
    /// <remarks>
    /// <para>Order-k block entropy H_k (Herzel, Ebeling &amp; Schmitt 1994; Schmitt &amp; Herzel 1997) per
    /// window of length W (W − k + 1 overlapping k-mers), plug-in estimate, maximum
    /// min(2k, log₂(W − k + 1)) bits for DNA. Kernel conventions apply: the window is upper-cased and every
    /// length-k substring is a symbol (N, IUPAC codes, gaps and non-letters are symbols, unlike the
    /// letter-only <see cref="CalculateEntropyProfile(string, int, int)"/>; for k = 1 the two coincide on
    /// letter-only input). A window shorter than k yields 0.</para>
    /// <para>Windows: offsets 0, step, 2·step, … ≤ length − W; with <paramref name="includePartialWindow"/>
    /// also the Biopython-style trailing truncated windows <c>seq[start:]</c> (as
    /// <see cref="CalculateEntropyProfile(string, int, int, bool)"/>). Oracle:
    /// <c>scipy.stats.entropy(Counter(w[j:j+k] for j in range(len(w)-k+1)).values(), base=2)</c>.</para>
    /// </remarks>
    /// <param name="sequence">Input sequence.</param>
    /// <param name="windowSize">Window width W (≥ 1).</param>
    /// <param name="stepSize">Window advance (≥ 1).</param>
    /// <param name="k">K-mer length (≥ 1).</param>
    /// <param name="includePartialWindow">Also emit trailing partial windows (default false).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="windowSize"/>, <paramref name="stepSize"/>
    /// or <paramref name="k"/> is less than 1 (validated eagerly).</exception>
    public static IEnumerable<double> CalculateEntropyProfile(
        string sequence,
        int windowSize,
        int stepSize,
        int k,
        bool includePartialWindow = false)
    {
        ValidateWindowAndStep(windowSize, stepSize);
        ArgumentOutOfRangeException.ThrowIfLessThan(k, 1);

        if (string.IsNullOrEmpty(sequence))
            return Array.Empty<double>();

        var starts = CompleteWindowStarts(sequence.Length, windowSize, stepSize);
        if (includePartialWindow)
            starts = starts.Concat(PartialWindowStarts(sequence.Length, windowSize, stepSize));

        return starts.Select(i => SequenceComplexity.CalculateKmerEntropy(
            sequence.Substring(i, Math.Min(windowSize, sequence.Length - i)), k));
    }

    // Offsets 0, step, 2·step, … ≤ length − windowSize of the complete windows (none when W > length).
    private static IEnumerable<int> CompleteWindowStarts(int length, int windowSize, int stepSize)
    {
        for (long i = 0; i <= (long)length - windowSize; i += stepSize)
            yield return (int)i;
    }

    #endregion

    #region Summary Statistics

    /// <summary>
    /// Comprehensive sequence statistics.
    /// </summary>
    public readonly record struct SequenceSummary(
        int Length,
        double GcContent,
        double Entropy,
        double Complexity,
        double MeltingTemperature,
        IReadOnlyDictionary<char, int> Composition)
    {
        /// <summary>
        /// Average-isotopic single-stranded linear molecular weight (Da) from
        /// <see cref="CalculateNucleotideMolecularWeight(string, bool)"/> (Biopython <c>molecular_weight</c>);
        /// RNA table when the sequence has U and no T, DNA table otherwise; unknown symbols skipped.
        /// Non-positional (2026-10 FIN-B03 F33) so the positional record shape is unchanged.
        /// </summary>
        public double MolecularWeight { get; init; }
    }

    /// <summary>
    /// Generates comprehensive summary statistics for a DNA/RNA sequence.
    /// </summary>
    /// <remarks>
    /// Pure aggregation — every field is the return value of its component method on the same input:
    /// <c>Length</c>/<c>GcContent</c>/<c>Composition</c> from <see cref="CalculateNucleotideComposition"/>
    /// (GC over A+C+G+T+U, Biopython <c>gc_fraction</c> "remove" mode without S/W),
    /// <c>Entropy</c> from <see cref="CalculateShannonEntropy"/> (bits, every letter a symbol),
    /// <c>Complexity</c> from <see cref="CalculateLinguisticComplexity"/> (maxK = 6; Σ V_k / Σ V_max,k, canonical <c>SequenceComplexity</c>) and
    /// <c>MeltingTemperature</c> from <see cref="CalculateMeltingTemperature"/> with the Wallace rule
    /// enabled (U read as T). <c>Composition</c> holds A, T, G, C, U, N only; other symbols are counted
    /// in <c>Length</c> but not listed. Null is treated as empty (all-zero summary).
    /// </remarks>
    public static SequenceSummary SummarizeNucleotideSequence(string? sequence)
        => SummarizeNucleotideSequence(sequence, weightAmbiguousBasesInTm: false);

    /// <summary>
    /// Generates the summary of <see cref="SummarizeNucleotideSequence(string?)"/>, optionally with the
    /// melting temperature weighted for IUPAC ambiguity codes as Biopython non-strict
    /// <c>Tm_Wallace</c>/<c>Tm_GC</c> (<see cref="CalculateMeltingTemperature(string, bool, bool)"/>).
    /// </summary>
    /// <param name="sequence">DNA/RNA sequence; null is treated as empty.</param>
    /// <param name="weightAmbiguousBasesInTm">False (default behaviour): Tm counts A/C/G/T/U only.</param>
    public static SequenceSummary SummarizeNucleotideSequence(string? sequence, bool weightAmbiguousBasesInTm)
    {
        // Treat null and empty identically (each per-metric method guards IsNullOrEmpty);
        // normalizing null to empty avoids dereferencing a null length when selecting the
        // Tm formula branch and keeps every per-metric call non-null.
        string seq = sequence ?? string.Empty;
        var comp = CalculateNucleotideComposition(seq);
        double entropy = CalculateShannonEntropy(seq);
        double complexity = CalculateLinguisticComplexity(seq);
        // Wallace rule applies to short oligos (< WallaceMaxLength = 14 A/C/G/T/U bases); the GC/Marmur-Doty
        // formula applies otherwise (SEQ-TM-001).
        // CalculateMeltingTemperature(useWallaceRule: true) itself selects the formula from the
        // number of A/C/G/T/U bases (OligoCalc), so N/gaps cannot push a short oligo into the GC formula.
        double tm = CalculateMeltingTemperature(seq, useWallaceRule: true, weightAmbiguousBases: weightAmbiguousBasesInTm);
        // Canonical nucleotide MW (Biopython molecular_weight, average masses, ss linear); RNA table for
        // U-only spellings so an RNA summary is not reduced to its non-U residues.
        bool isRna = comp.CountU > 0 && comp.CountT == 0;
        double molecularWeight = CalculateNucleotideMolecularWeight(seq, isDna: !isRna);

        var composition = new Dictionary<char, int>
        {
            { 'A', comp.CountA },
            { 'T', comp.CountT },
            { 'G', comp.CountG },
            { 'C', comp.CountC },
            { 'U', comp.CountU },
            { 'N', comp.CountN }
        };

        return new SequenceSummary(
            Length: comp.Length,
            GcContent: comp.GcContent,
            Entropy: entropy,
            Complexity: complexity,
            MeltingTemperature: tm,
            Composition: composition)
        {
            MolecularWeight = molecularWeight
        };
    }

    #endregion
}
