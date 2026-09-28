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
    /// <para><c>Counts</c> holds the case-insensitive count of every letter (the 20 standard residues
    /// plus any ambiguity/extended codes such as B, Z, X, U, O, J). <c>Length</c> is the number of
    /// letters; non-letter symbols (stop <c>*</c>, gap <c>-</c>, digits, whitespace) are not residues
    /// and are excluded from <c>Counts</c>, <c>Length</c> and the ratio denominators. This differs from
    /// Biopython <c>ProteinAnalysis</c>, whose denominator is <c>len(seq)</c> including such symbols;
    /// for letter-only input the two agree exactly.</para>
    /// <para><c>AromaticResidueRatio</c> = (F+W+Y)/Length — the aromaticity of Lobry &amp; Gautier
    /// (1994) as implemented by Biopython <c>ProteinAnalysis.aromaticity()</c> (EMBOSS pepstats'
    /// "Aromatic" class additionally includes H). <c>ChargedResidueRatio</c> = (D+E+H+K+R)/Length —
    /// the EMBOSS pepstats "Charged" class (B+D+E+H+K+R+Z) restricted to unambiguous residues.</para>
    /// <para><c>MolecularWeight</c>, <c>IsoelectricPoint</c> and <c>Hydrophobicity</c> delegate to
    /// <see cref="CalculateMolecularWeight"/>, <see cref="CalculateIsoelectricPoint"/> and
    /// <see cref="CalculateHydrophobicity"/>. Null/empty input returns Length 0, empty counts,
    /// ratios 0 and the neutral-pH pI default of those methods.</para>
    /// </remarks>
    public static AminoAcidComposition CalculateAminoAcidComposition(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
        {
            return new AminoAcidComposition(0, new Dictionary<char, int>(), 0, NeutralPhDefault, 0, 0, 0);
        }

        var counts = new Dictionary<char, int>();
        foreach (char ch in sequence.ToUpperInvariant())
        {
            if (char.IsLetter(ch))
            {
                counts[ch] = counts.GetValueOrDefault(ch) + 1;
            }
        }

        int length = counts.Values.Sum();
        double mw = CalculateMolecularWeight(sequence);
        double pi = CalculateIsoelectricPoint(sequence);
        double hydro = CalculateHydrophobicity(sequence);

        // Charged residues: D, E (negative), K, R, H (positive)
        int charged = counts.GetValueOrDefault('D') + counts.GetValueOrDefault('E') +
                     counts.GetValueOrDefault('K') + counts.GetValueOrDefault('R') +
                     counts.GetValueOrDefault('H');
        double chargedRatio = length > 0 ? (double)charged / length : 0;

        // Aromatic residues: F, Y, W
        int aromatic = counts.GetValueOrDefault('F') + counts.GetValueOrDefault('Y') +
                      counts.GetValueOrDefault('W');
        double aromaticRatio = length > 0 ? (double)aromatic / length : 0;

        return new AminoAcidComposition(
            Length: length,
            Counts: counts,
            MolecularWeight: mw,
            IsoelectricPoint: pi,
            Hydrophobicity: hydro,
            ChargedResidueRatio: chargedRatio,
            AromaticResidueRatio: aromaticRatio);
    }

    #endregion

    #region Molecular Weight

    // Average isotopic mass of one water molecule (Da).
    // Expasy FindMod "Other mass values": H2O = 18.01524; Biopython SeqUtils uses water = 18.0153.
    // Source: https://web.expasy.org/findmod/findmod_masses.html ;
    //         Biopython Bio/SeqUtils/__init__.py (molecular_weight).
    private const double AverageWaterMass = 18.0153;

    // Average molecular masses of the 20 standard free amino acids plus the two
    // genetically encoded non-standard ones, selenocysteine (U) and pyrrolysine (O) (Da).
    // Source: Biopython Bio/Data/IUPACData.py `protein_weights` (master, 22 entries incl. O/U),
    // "Mass data taken from PubChem"; consistent with Expasy FindMod average residue masses +
    // AverageWaterMass. Expasy ProtParam / Compute pI/Mw likewise accept U and O.
    private static readonly Dictionary<char, double> AminoAcidWeights = new()
    {
        { 'A', 89.0932 },  { 'C', 121.1582 }, { 'D', 133.1027 }, { 'E', 147.1293 },
        { 'F', 165.1891 }, { 'G', 75.0666 },  { 'H', 155.1546 }, { 'I', 131.1729 },
        { 'K', 146.1876 }, { 'L', 131.1729 }, { 'M', 149.2113 }, { 'N', 132.1179 },
        { 'O', 255.3134 }, { 'P', 115.1305 }, { 'Q', 146.1445 }, { 'R', 174.201 },
        { 'S', 105.0926 }, { 'T', 119.1192 }, { 'U', 168.0532 }, { 'V', 117.1463 },
        { 'W', 204.2252 }, { 'Y', 181.1885 }
    };

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
    /// </remarks>
    /// <param name="proteinSequence">Protein sequence (case-insensitive, one-letter codes).</param>
    /// <returns>Molecular weight in daltons; 0 for null/empty input.</returns>
    public static double CalculateMolecularWeight(string proteinSequence)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            return 0;

        double weight = 0;
        int residues = 0;

        foreach (char aa in proteinSequence.ToUpperInvariant())
        {
            if (AminoAcidWeights.TryGetValue(aa, out double aaWeight))
            {
                weight += aaWeight;
                residues++;
            }
        }

        if (residues == 0)
            return 0;

        // One water is lost per peptide bond; n residues form (n − 1) bonds.
        return weight - (residues - 1) * AverageWaterMass;
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
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        Dictionary<char, double> table = isDna ? DnaNucleotideWeights : RnaNucleotideWeights;

        double weight = 0;
        double complementWeight = 0;
        int monomers = 0;

        foreach (char ch in sequence.ToUpperInvariant())
        {
            if (table.TryGetValue(ch, out double ntWeight))
            {
                weight += ntWeight;
                monomers++;
                if (doubleStranded)
                {
                    char complement = isDna
                        ? Seqeron.Genomics.Core.SequenceExtensions.GetComplementBase(ch)
                        : Seqeron.Genomics.Core.SequenceExtensions.GetRnaComplementBase(ch);
                    complementWeight += table[complement];
                }
            }
        }

        if (monomers == 0)
            return 0;

        // One water is lost per phosphodiester bond: n monomers form (n − 1) bonds in a linear
        // strand and n bonds in a circular one.
        int bondsPerStrand = circular ? monomers : monomers - 1;
        double result = weight - bondsPerStrand * AverageWaterMass;

        if (doubleStranded)
            result += complementWeight - bondsPerStrand * AverageWaterMass;

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

    // EMBOSS 6.6.0 emboss/data/Epk.dat (the file shipped with and read by `iep`; embiep.c
    // embIepPkReadFile also defaults amino = 7.50 when the file has no "Amino" line).
    // NOTE: the Epk.dat listing printed on the EMBOSS iep web page (Amino 8.6) is stale — the
    // page's own worked outputs (LACI_ECOLI pI 6.8385, IFNA2_HUMAN pI 5.7240) are only reproduced
    // with Amino 7.5 (verified against the EMBOSS 6.6.0 `iep` binary).
    private const double EmbossNTerminusPka = 7.5;
    private const double EmbossCTerminusPka = 3.6;

    // Ionizable side chains: pKa and sign (+1 basic, -1 acidic).
    private static readonly Dictionary<char, (double pKa, int charge)> EmbossSideChains = new()
    {
        { 'C', (8.5, -1) }, { 'D', (3.9, -1) }, { 'E', (4.1, -1) }, { 'Y', (10.1, -1) },
        { 'H', (6.5, 1) },  { 'K', (10.8, 1) }, { 'R', (12.5, 1) }
    };

    // Bjellqvist et al. 1993 (Electrophoresis 14:1023) / 1994 (Electrophoresis 15:529) pK set,
    // as implemented by Biopython Bio/SeqUtils/IsoelectricPoint.py (positive_pKs, negative_pKs,
    // pKnterminal, pKcterminal) and ExPASy Compute pI/Mw.
    private static readonly Dictionary<char, (double pKa, int charge)> BjellqvistSideChains = new()
    {
        { 'C', (9.0, -1) }, { 'D', (4.05, -1) }, { 'E', (4.45, -1) }, { 'Y', (10.0, -1) },
        { 'H', (5.98, 1) }, { 'K', (10.0, 1) },  { 'R', (12.0, 1) }
    };

    private const double BjellqvistNTerminusPka = 7.5;
    private const double BjellqvistCTerminusPka = 3.55;

    private static readonly Dictionary<char, double> BjellqvistNTerminalResiduePka = new()
    {
        { 'A', 7.59 }, { 'M', 7.0 }, { 'S', 6.93 }, { 'P', 8.36 }, { 'T', 6.82 }, { 'V', 7.44 }, { 'E', 7.7 }
    };

    private static readonly Dictionary<char, double> BjellqvistCTerminalResiduePka = new()
    {
        { 'D', 4.55 }, { 'E', 4.75 }
    };

    // Bisection window and convergence. The root is located to 1e-9 pH and then rounded, so the
    // returned value is the correctly rounded pI (a 0.01-wide final bracket could round wrongly).
    private const double MinPh = 0.0;
    private const double MaxPh = 14.0;
    private const double PiBisectionPrecision = 1e-9;

    // pI returned for empty/null input: pI is undefined for a zero-length protein (a real
    // protein always has both termini); neutral 7.0 is used as a documented input-guard sentinel.
    private const double NeutralPhDefault = 7.0;

    private const int PiDecimalPlaces = 2;

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
    /// and reports "none" when the charge does not change sign there; this method searches [0, 14].</para>
    /// <para>Bjellqvist: reproduces Biopython <c>IsoelectricPoint.pi()</c> / ExPASy Compute pI; Biopython
    /// clamps its bisection to [4.05, 12], so for extremely acidic/basic peptides whose true root lies
    /// outside that window Biopython returns the window edge while this method returns the true root.</para>
    /// </remarks>
    public static double CalculateIsoelectricPoint(string proteinSequence, PkaScale scale)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            return NeutralPhDefault;

        var model = BuildChargeModel(proteinSequence, scale);

        double pHLow = MinPh;
        double pHHigh = MaxPh;

        while (pHHigh - pHLow > PiBisectionPrecision)
        {
            double pH = (pHLow + pHHigh) / 2.0;

            // Net charge is monotonically non-increasing in pH: positive ⇒ pI is higher.
            if (model.NetCharge(pH) > 0)
                pHLow = pH;
            else
                pHHigh = pH;
        }

        return Math.Round((pHLow + pHHigh) / 2.0, PiDecimalPlaces);
    }

    /// <summary>
    /// Net charge of a protein at a given pH (Henderson–Hasselbalch): basic groups (N-terminus,
    /// K, R, H) contribute +1/(1+10^(pH−pKa)), acidic groups (C-terminus, D, E, C, Y) contribute
    /// −1/(1+10^(pKa−pH)); each terminus is counted once.
    /// Equivalent to EMBOSS <c>embIepGetCharge</c> and Biopython <c>IsoelectricPoint.charge_at_pH</c>.
    /// </summary>
    /// <param name="proteinSequence">Single-letter amino-acid sequence (case-insensitive).</param>
    /// <param name="pH">pH at which to evaluate the charge.</param>
    /// <param name="scale">pK set (default EMBOSS).</param>
    /// <returns>Net charge (elementary charges); 0 for null/empty input.</returns>
    public static double CalculateNetCharge(string proteinSequence, double pH, PkaScale scale = PkaScale.Emboss)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            return 0.0;

        return BuildChargeModel(proteinSequence, scale).NetCharge(pH);
    }

    private readonly record struct ChargeModel(
        double NTerminusPka,
        double CTerminusPka,
        Dictionary<char, (double pKa, int charge)> SideChains,
        Dictionary<char, int> Counts)
    {
        public double NetCharge(double pH)
        {
            double charge = 1.0 / (1.0 + Math.Pow(10, pH - NTerminusPka));
            charge -= 1.0 / (1.0 + Math.Pow(10, CTerminusPka - pH));

            foreach (var (aa, count) in Counts)
            {
                var (pKa, sign) = SideChains[aa];
                if (sign > 0)
                    charge += count / (1.0 + Math.Pow(10, pH - pKa));   // basic group
                else
                    charge -= count / (1.0 + Math.Pow(10, pKa - pH));   // acidic group
            }

            return charge;
        }
    }

    private static ChargeModel BuildChargeModel(string proteinSequence, PkaScale scale)
    {
        string upper = proteinSequence.ToUpperInvariant();
        var sideChains = scale == PkaScale.Bjellqvist ? BjellqvistSideChains : EmbossSideChains;

        var counts = new Dictionary<char, int>();
        int countB = 0, countZ = 0;
        foreach (char aa in upper)
        {
            if (sideChains.ContainsKey(aa))
                counts[aa] = counts.GetValueOrDefault(aa) + 1;
            else if (aa == 'B')
                countB++;
            else if (aa == 'Z')
                countZ++;
        }

        if (scale == PkaScale.Bjellqvist)
        {
            // Terminal-residue-specific pKs (Bjellqvist 1994; Biopython _update_pKs_tables).
            double nPka = BjellqvistNTerminalResiduePka.GetValueOrDefault(upper[0], BjellqvistNTerminusPka);
            double cPka = BjellqvistCTerminalResiduePka.GetValueOrDefault(upper[^1], BjellqvistCTerminusPka);
            return new ChargeModel(nPka, cPka, sideChains, counts);
        }

        // EMBOSS embIepCompC: B = D or N, Z = E or Q, split by Dayhoff frequencies
        // (D 5.5 / N 4.3; E 6.0 / Q 3.9), rounding half up via (int)(0.5 + x).
        if (countB > 0)
        {
            int asp = (int)(0.5 + countB * 5.5 / 9.8);
            if (asp > 0) counts['D'] = counts.GetValueOrDefault('D') + asp;
        }
        if (countZ > 0)
        {
            int glu = (int)(0.5 + countZ * 6.0 / 9.9);
            if (glu > 0) counts['E'] = counts.GetValueOrDefault('E') + glu;
        }

        return new ChargeModel(EmbossNTerminusPka, EmbossCTerminusPka, sideChains, counts);
    }

    #endregion

    #region Hydrophobicity

    // Kyte-Doolittle hydropathy scale (kcal-derived index, the 20 standard residues).
    // Values per Kyte J, Doolittle RF (1982) J Mol Biol 157:105-132; identical to
    // Biopython Bio.SeqUtils.ProtParamData.kd.
    private static readonly Dictionary<char, double> HydrophobicityScale = new()
    {
        { 'A', 1.8 },  { 'R', -4.5 }, { 'N', -3.5 }, { 'D', -3.5 },
        { 'C', 2.5 },  { 'E', -3.5 }, { 'Q', -3.5 }, { 'G', -0.4 },
        { 'H', -3.2 }, { 'I', 4.5 },  { 'L', 3.8 },  { 'K', -3.9 },
        { 'M', 1.9 },  { 'F', 2.8 },  { 'P', -1.6 }, { 'S', -0.8 },
        { 'T', -0.7 }, { 'W', -0.9 }, { 'Y', -1.3 }, { 'V', 4.2 }
    };

    // Default sliding-window size; Kyte & Doolittle (1982) report window 9 gives the
    // best results for surface regions of globular proteins (per GCAT/Davidson summary).
    private const int DefaultHydropathyWindow = 9;

    /// <summary>
    /// Calculates the grand average of hydropathy (GRAVY) index: the sum of Kyte-Doolittle
    /// hydropathy values of all recognized residues divided by the residue count.
    /// Positive values indicate hydrophobic, negative hydrophilic. Non-standard residues
    /// are skipped (not counted). Returns 0 for null/empty input or no recognized residues.
    /// </summary>
    public static double CalculateHydrophobicity(string proteinSequence)
    {
        if (string.IsNullOrEmpty(proteinSequence))
            return 0;

        double sum = 0;
        int count = 0;

        foreach (char aa in proteinSequence.ToUpperInvariant())
        {
            if (HydrophobicityScale.TryGetValue(aa, out double value))
            {
                sum += value;
                count++;
            }
        }

        return count > 0 ? sum / count : 0;
    }

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
    /// counts residue W/2 twice and divides by W + 1 — an artefact of its odd-window loop).</para>
    /// <para>Non-standard residues (B, Z, X, gaps, stop) have no scale value and contribute 0 to the
    /// weighted sum while keeping their weight in the divisor — identical to Biopython for such a
    /// residue at the window centre (Biopython also drops the symmetric partner of an off-centre
    /// unknown residue; this library does not).</para>
    /// </remarks>
    /// <param name="proteinSequence">One-letter amino-acid sequence (case-insensitive).</param>
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

        return HydrophobicityProfileIterator(proteinSequence.ToUpperInvariant(), windowSize, edgeWeight);
    }

    private static IEnumerable<double> HydrophobicityProfileIterator(string upper, int windowSize, double edgeWeight)
    {
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

        for (int i = 0; i <= upper.Length - windowSize; i++)
        {
            double sum = 0;
            for (int j = 0; j < windowSize; j++)
            {
                if (HydrophobicityScale.TryGetValue(upper[i + j], out double value))
                    sum += weights[j] * value;
            }
            yield return sum / weightSum;
        }
    }

    #endregion

    #region DNA Thermodynamics

    // Unified nearest-neighbor thermodynamic parameters for Watson-Crick base pairs
    // in 1 M NaCl, ΔH in kcal/mol and ΔS in cal/(mol·K).
    // Source: Allawi HT, SantaLucia J (1997), Biochemistry 36(34):10581-10594,
    //         Table 1 (also tabulated as DNA_NN3 in Biopython Bio.SeqUtils.MeltingTemp).
    // Each entry and its Watson-Crick complement share the same values
    // (e.g. AA/TT, CA/TG), so all 16 dinucleotides are listed explicitly.
    private static readonly Dictionary<string, (double dH, double dS)> NearestNeighborParams = new()
    {
        { "AA", (-7.9, -22.2) }, { "TT", (-7.9, -22.2) },
        { "AT", (-7.2, -20.4) },
        { "TA", (-7.2, -21.3) },
        { "CA", (-8.5, -22.7) }, { "TG", (-8.5, -22.7) },
        { "GT", (-8.4, -22.4) }, { "AC", (-8.4, -22.4) },
        { "CT", (-7.8, -21.0) }, { "AG", (-7.8, -21.0) },
        { "GA", (-8.2, -22.2) }, { "TC", (-8.2, -22.2) },
        { "CG", (-10.6, -27.2) },
        { "GC", (-9.8, -24.4) },
        { "GG", (-8.0, -19.9) }, { "CC", (-8.0, -19.9) }
    };

    // Helix-initiation parameters applied once at EACH duplex terminus,
    // selected by whether that terminal base pair is G·C or A·T.
    // Source: Allawi & SantaLucia (1997), Table 1; "init. w/ term. G·C" and
    // "init. w/ term. A·T" (Biopython DNA_NN3 keys init_G/C, init_A/T).
    private const double InitTerminalGcDeltaH = 0.1;   // kcal/mol
    private const double InitTerminalGcDeltaS = -2.8;  // cal/(mol·K)
    private const double InitTerminalAtDeltaH = 2.3;   // kcal/mol
    private const double InitTerminalAtDeltaS = 4.1;   // cal/(mol·K)

    // Salt (Na+) entropy correction per SantaLucia (1998) "method 5":
    // ΔS(salt) = 0.368 * (N-1) * ln[Na+], [Na+] in mol/L.
    // Source: SantaLucia J (1998), PNAS 95(4):1460-1465; Biopython salt_correction method 5.
    private const double SaltEntropyCoefficient = 0.368;

    // Gas constant R in cal/(mol·K) used in the Tm equation.
    // Source: Allawi & SantaLucia (1997); Biopython Tm_NN uses R = 1.987.
    private const double GasConstantCalPerMolK = 1.987;

    // Reference temperature for ΔG° calculation: 37 °C = 310.15 K.
    private const double ReferenceTemperatureKelvin = 310.15;

    // Kelvin-to-Celsius offset.
    private const double KelvinToCelsiusOffset = 273.15;

    // Total-strand-concentration divisor F (x) for the Tm equation Tm = ΔH°/(ΔS° + R·ln(C_T/x)).
    // x = 4 for two non-self-complementary strands in equal amount, x = 1 for a
    // self-complementary (homo)duplex. Source: SantaLucia (1998) PNAS 95:1460 Eq. 3;
    // Biopython Tm_NN: k = dnac1 − dnac2/2 (non-self-complementary), k = dnac1 (selfcomp).
    private const double NonSelfComplementaryFactor = 4.0;
    private const double SelfComplementaryFactor = 1.0;

    // Symmetry correction for a self-complementary duplex (ΔH° = 0, ΔS° = −1.4 cal/(mol·K)).
    // Source: SantaLucia (1998) Table 2 "symmetry correction"; Biopython DNA_NN3 key "sym" (0, −1.4).
    private const double SymmetryCorrectionDeltaS = -1.4;

    /// <summary>
    /// DNA thermodynamic properties.
    /// </summary>
    public readonly record struct ThermodynamicProperties(
        double DeltaH,
        double DeltaS,
        double DeltaG,
        double MeltingTemperature);

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
    /// other than A/C/G/T removed before the model is applied. Fewer than 2 remaining bases → all zero.</param>
    /// <param name="naConcentration">Na+ concentration in mol/L (default 0.05 = 50 mM); must be &gt; 0.</param>
    /// <param name="primerConcentration">
    /// Total strand concentration C_T in mol/L (default 2.5e-7 = 250 nM; must be &gt; 0); the Tm
    /// equation divides this by F = 4 for two non-self-complementary strands in equal amount.
    /// </param>
    /// <returns>
    /// ΔH° (kcal/mol), ΔS° (cal/(mol·K), salt-corrected), ΔG°₃₇ (kcal/mol) rounded to 2 decimals
    /// and Tm (°C) rounded to 1 decimal. For null input or fewer than 2 A/C/G/T(U) bases all four
    /// fields are 0.
    /// </returns>
    /// <exception cref="ArgumentOutOfRangeException">A concentration is not a positive finite number.</exception>
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
        if (!(primerConcentration > 0) || double.IsInfinity(primerConcentration))
            throw new ArgumentOutOfRangeException(nameof(primerConcentration), primerConcentration,
                "Strand concentration must be a positive finite value in mol/L.");

        string seq = NormalizeForNearestNeighbor(dnaSequence);
        if (seq.Length < 2)
            return new ThermodynamicProperties(0, 0, 0, 0);

        // Calculate ΔH and ΔS using the nearest-neighbor method.
        double dH = 0;
        double dS = 0;

        // Helix-initiation parameters are applied at BOTH duplex termini
        // (first and last base pair) per Allawi & SantaLucia (1997), Table 1.
        AddTerminalInitiation(seq[0], ref dH, ref dS);
        AddTerminalInitiation(seq[^1], ref dH, ref dS);

        // Sum nearest-neighbor contributions over each overlapping dinucleotide
        // (the normalised sequence is pure A/C/G/T, so every step is in the table).
        for (int i = 0; i < seq.Length - 1; i++)
        {
            var param = NearestNeighborParams[seq.Substring(i, 2)];
            dH += param.dH;
            dS += param.dS;
        }

        double strandFactor = NonSelfComplementaryFactor;
        if (selfComplementary)
        {
            dS += SymmetryCorrectionDeltaS;
            strandFactor = SelfComplementaryFactor;
        }

        // Salt correction for ΔS (SantaLucia 1998, method 5); N = number of bases.
        dS += SaltEntropyCoefficient * (seq.Length - 1) * Math.Log(naConcentration);

        // ΔG° at 37 °C: ΔG° = ΔH° - T·ΔS° (ΔS° converted from cal to kcal).
        double dG = dH - (ReferenceTemperatureKelvin * dS / 1000.0);

        // Tm = ΔH° / (ΔS° + R · ln(C_T / x)) - 273.15, with ΔH° converted to cal.
        double tm = (dH * 1000) /
                    (dS + GasConstantCalPerMolK * Math.Log(primerConcentration / strandFactor))
                    - KelvinToCelsiusOffset;

        return new ThermodynamicProperties(
            DeltaH: Math.Round(dH, 2),
            DeltaS: Math.Round(dS, 2),
            DeltaG: Math.Round(dG, 2),
            MeltingTemperature: Math.Round(tm, 1));
    }

    // Biopython MeltingTemp._check(seq, "Tm_NN"): upper-case, drop whitespace, back-transcribe
    // (U → T), keep only bases the NN table can score. DNA_NN3 has no inosine (I) parameters
    // (Biopython raises on I), so only A/C/G/T are kept.
    private static string NormalizeForNearestNeighbor(string? sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return string.Empty;

        var sb = new System.Text.StringBuilder(sequence.Length);
        foreach (char ch in sequence)
        {
            char b = char.ToUpperInvariant(ch);
            if (b == 'U') b = 'T';
            if (b is 'A' or 'C' or 'G' or 'T')
                sb.Append(b);
        }
        return sb.ToString();
    }

    // Adds the helix-initiation contribution for one terminal base (A/C/G/T).
    private static void AddTerminalInitiation(char terminalBase, ref double dH, ref double dS)
    {
        if (terminalBase is 'G' or 'C')
        {
            dH += InitTerminalGcDeltaH;
            dS += InitTerminalGcDeltaS;
        }
        else
        {
            dH += InitTerminalAtDeltaH;
            dS += InitTerminalAtDeltaS;
        }
    }

    /// <summary>
    /// Calculates the "basic" melting temperature of an oligonucleotide (OligoCalc, Kibbe 2007,
    /// NAR 35:W43): Wallace rule Tm = 2(A+T) + 4(G+C) (Thein &amp; Wallace 1986) for fewer than 14
    /// bases, otherwise Tm = 64.9 + 41·(G+C − 16.4)/N (customarily attributed to Marmur &amp; Doty 1962),
    /// both at OligoCalc's fixed standard conditions (50 nM primer, 50 mM Na+, pH 7.0).
    /// </summary>
    /// <remarks>
    /// Only A, C, G, T are counted (case-insensitive); every other character (N, IUPAC codes, U,
    /// gaps, whitespace) is ignored. As in OligoCalc — whose GC formula divides by
    /// N = wA+xT+yG+zC — the length that selects the formula is the number of counted bases, so the
    /// result with <paramref name="useWallaceRule"/> = true equals the canonical
    /// <c>PrimerDesigner.CalculateMeltingTemperature</c> (MolTools; not callable from this assembly).
    /// With <paramref name="useWallaceRule"/> = false the GC formula is applied at any length; below
    /// 14 bases this is outside its published domain and can be negative. Returns 0 for null/empty
    /// input or when no A/C/G/T base is present.
    /// </remarks>
    public static double CalculateMeltingTemperature(string dnaSequence, bool useWallaceRule = true)
    {
        if (string.IsNullOrEmpty(dnaSequence))
            return 0;

        var comp = CalculateNucleotideComposition(dnaSequence);
        int at = comp.CountA + comp.CountT;
        int gc = comp.CountG + comp.CountC;
        int validLength = at + gc;
        if (validLength == 0)
            return 0;

        if (useWallaceRule && validLength < ThermoConstants.WallaceMaxLength)
        {
            // Wallace rule for short oligos: Tm = 2(A+T) + 4(G+C)
            return ThermoConstants.CalculateWallaceTm(at, gc);
        }

        // GC formula (Marmur-Doty / OligoCalc basic Tm)
        return ThermoConstants.CalculateMarmurDotyTm(gc, validLength);
    }

    #endregion

    #region Sequence Patterns

    /// <summary>
    /// Calculates normalized dinucleotide frequencies f_XY = count(XY) / (number of dinucleotide
    /// positions, i.e. N-1) over the alphabet {A,T,G,C,U}. Non-alphabet dinucleotides are excluded.
    /// Frequency normalization follows the Karlin genomic-signature convention
    /// (Karlin S., "Pervasive properties of the genomic signature", PMC126251).
    /// </summary>
    public static IReadOnlyDictionary<string, double> CalculateDinucleotideFrequencies(string sequence)
    {
        var counts = new Dictionary<string, int>();
        var freq = new Dictionary<string, double>();

        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return freq;

        string upper = sequence.ToUpperInvariant();

        // Count all dinucleotides
        int total = 0;
        for (int i = 0; i < upper.Length - 1; i++)
        {
            string dinuc = upper.Substring(i, 2);
            if (dinuc.All(c => "ATGCU".Contains(c)))
            {
                counts[dinuc] = counts.GetValueOrDefault(dinuc) + 1;
                total++;
            }
        }

        // Convert to frequencies
        foreach (var (dinuc, count) in counts)
        {
            freq[dinuc] = (double)count / total;
        }

        return freq;
    }

    /// <summary>
    /// Calculates the dinucleotide relative abundance (odds ratio) ρ_XY = f_XY / (f_X · f_Y),
    /// where f_XY is the dinucleotide frequency and f_X, f_Y are single-base frequencies.
    /// ρ = 1 indicates no bias (observed equals the product of base frequencies); ρ &gt; 1
    /// over-representation and ρ &lt; 1 under-representation
    /// (Karlin S., PMC126251; Karlin &amp; Burge 1995, Trends Genet 11(7):283-290).
    /// When a constituent base is absent the expected frequency is 0 and the ratio is reported as 0.
    /// </summary>
    public static IReadOnlyDictionary<string, double> CalculateDinucleotideRatios(string sequence)
    {
        var ratios = new Dictionary<string, double>();

        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return ratios;

        var comp = CalculateNucleotideComposition(sequence);
        var dinucFreq = CalculateDinucleotideFrequencies(sequence);

        int total = comp.CountA + comp.CountT + comp.CountG + comp.CountC + comp.CountU;
        if (total == 0) return ratios;

        // Single nucleotide frequencies
        var singleFreq = new Dictionary<char, double>
        {
            { 'A', (double)comp.CountA / total },
            { 'T', (double)comp.CountT / total },
            { 'G', (double)comp.CountG / total },
            { 'C', (double)comp.CountC / total },
            { 'U', (double)comp.CountU / total }
        };

        // Calculate observed/expected ratios
        foreach (var (dinuc, observed) in dinucFreq)
        {
            double expected = singleFreq.GetValueOrDefault(dinuc[0]) *
                             singleFreq.GetValueOrDefault(dinuc[1]);
            ratios[dinuc] = expected > 0 ? observed / expected : 0;
        }

        return ratios;
    }

    /// <summary>
    /// Calculates codon usage frequencies by reading consecutive, non-overlapping triplets from the
    /// given reading frame: frequency = count(codon) / total counted codons. Triplets containing any
    /// non-ACGT base are excluded and trailing 1-2 leftover bases are ignored. This is the count/total
    /// fraction used by the Kazusa Codon Usage Database (CUTG); it equals the CUTG per-thousand
    /// frequency divided by 1000, and is distinct from the per-amino-acid "fraction" reported by
    /// EMBOSS cusp. Input shorter than 3 bases, or with no valid codon (total = 0), yields an empty
    /// table. Sources: Nakamura, Gojobori, Ikemura (2000), Nucleic Acids Res 28(1):292,
    /// DOI 10.1093/nar/28.1.292; Kazusa CUTG README, https://www.kazusa.or.jp/codon/readme_codon.html.
    /// </summary>
    /// <param name="dnaSequence">DNA coding sequence; case-insensitive, non-ACGT bases excluded.</param>
    /// <param name="readingFrame">0-based offset of the first codon (0, 1, or 2 in practice).</param>
    /// <returns>Map of codon to its frequency (count / total counted codons); empty if no valid codon.</returns>
    public static IReadOnlyDictionary<string, double> CalculateCodonFrequencies(
        string dnaSequence,
        int readingFrame = 0)
    {
        // Codon length is fixed by the genetic code (non-overlapping triplets), per Kazusa CUTG.
        const int CodonLength = 3;

        var counts = new Dictionary<string, int>();
        var freq = new Dictionary<string, double>();

        if (string.IsNullOrEmpty(dnaSequence) || dnaSequence.Length < CodonLength)
            return freq;

        string upper = dnaSequence.ToUpperInvariant();
        int total = 0;

        for (int i = readingFrame; i <= upper.Length - CodonLength; i += CodonLength)
        {
            string codon = upper.Substring(i, CodonLength);
            // Kazusa CUTG: codons containing an ambiguous (non-ACGT) base are excluded from the count.
            if (codon.All(c => "ATGC".Contains(c)))
            {
                counts[codon] = counts.GetValueOrDefault(codon) + 1;
                total++;
            }
        }

        // total == 0 (no valid codon) leaves counts empty, so this loop yields an empty table
        // and avoids any division by zero.
        foreach (var (codon, count) in counts)
        {
            freq[codon] = (double)count / total;
        }

        return freq;
    }

    #endregion

    #region Entropy and Complexity

    /// <summary>
    /// Calculates Shannon entropy of a sequence.
    /// </summary>
    public static double CalculateShannonEntropy(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        var counts = new Dictionary<char, int>();
        int total = 0;

        foreach (char ch in sequence.ToUpperInvariant())
        {
            if (char.IsLetter(ch))
            {
                counts[ch] = counts.GetValueOrDefault(ch) + 1;
                total++;
            }
        }

        if (total == 0) return 0;

        double entropy = 0;
        foreach (int count in counts.Values)
        {
            double freq = (double)count / total;
            if (freq > 0)
            {
                entropy -= freq * Math.Log2(freq);
            }
        }

        return entropy;
    }

    /// <summary>
    /// Calculates linguistic complexity of a sequence.
    /// </summary>
    public static double CalculateLinguisticComplexity(string sequence, int maxK = 6)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        string upper = sequence.ToUpperInvariant();
        int n = upper.Length;
        double totalRatio = 0;
        int kCount = 0;

        for (int k = 1; k <= Math.Min(maxK, n); k++)
        {
            var observedKmers = new HashSet<string>();
            for (int i = 0; i <= n - k; i++)
            {
                observedKmers.Add(upper.Substring(i, k));
            }

            // Maximum possible k-mers
            int maxPossible = Math.Min((int)Math.Pow(4, k), n - k + 1);
            if (maxPossible > 0)
            {
                totalRatio += (double)observedKmers.Count / maxPossible;
                kCount++;
            }
        }

        return kCount > 0 ? totalRatio / kCount : 0;
    }

    #endregion

    #region Protein Secondary Structure

    // Chou-Fasman conformational parameters (Pa = helix, Pb = sheet, Pt = turn),
    // expressed as propensities (the published integer parameters / 100).
    // Values verbatim from Chou PY, Fasman GD (1978) "Empirical predictions of protein
    // conformation" Annu Rev Biochem 47:251-276, as reproduced in the cited academic
    // tables and reference implementation. See docs/Evidence/SEQ-SECSTRUCT-001-Evidence.md.
    private static readonly Dictionary<char, (double Helix, double Sheet, double Turn)> SecondaryStructurePropensity = new()
    {
        { 'A', (1.42, 0.83, 0.66) }, { 'R', (0.98, 0.93, 0.95) },
        { 'N', (0.67, 0.89, 1.56) }, { 'D', (1.01, 0.54, 1.46) },
        { 'C', (0.70, 1.19, 1.19) }, { 'E', (1.51, 0.37, 0.74) },
        { 'Q', (1.11, 1.10, 0.98) }, { 'G', (0.57, 0.75, 1.56) },
        { 'H', (1.00, 0.87, 0.95) }, { 'I', (1.08, 1.60, 0.47) },
        { 'L', (1.21, 1.30, 0.59) }, { 'K', (1.14, 0.74, 1.01) },
        { 'M', (1.45, 1.05, 0.60) }, { 'F', (1.13, 1.38, 0.60) },
        { 'P', (0.57, 0.55, 1.52) }, { 'S', (0.77, 0.75, 1.43) },
        { 'T', (0.83, 1.19, 0.96) }, { 'W', (1.08, 1.37, 0.96) },
        { 'Y', (0.69, 1.47, 1.14) }, { 'V', (1.06, 1.70, 0.50) }
    };

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
    /// <param name="windowSize">Window width W in bases (default 100). Must be ≤ sequence length for any window to be produced.</param>
    /// <param name="stepSize">Window advance in bases (default 1).</param>
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
    /// equivalent to T.
    /// </remarks>
    public static IEnumerable<double> CalculateGcContentProfile(
        string sequence,
        int windowSize = DefaultGcProfileWindow,
        int stepSize = 1,
        bool fraction = false)
    {
        if (string.IsNullOrEmpty(sequence) || windowSize > sequence.Length)
            yield break;

        string upper = sequence.ToUpperInvariant();

        // Opt-in Biopython convention: when fraction == true, emit GC in [0,1] (matching
        // Bio.SeqUtils.gc_fraction) instead of the default percentage [0,100]. The default
        // (false) is unchanged.
        double scale = fraction ? 1.0 : PercentScale;

        for (int i = 0; i <= upper.Length - windowSize; i += stepSize)
        {
            int gc = 0;
            int total = 0;

            for (int j = 0; j < windowSize; j++)
            {
                char ch = upper[i + j];
                if (ch == 'G' || ch == 'C')
                {
                    gc++;
                    total++;
                }
                else if (ch == 'A' || ch == 'T' || ch == 'U')
                {
                    total++;
                }
            }

            yield return total > 0 ? (double)gc / total * scale : 0;
        }
    }

    /// <summary>
    /// Calculates the Shannon entropy profile of a sequence: the per-symbol Shannon
    /// entropy H = -Σ pᵢ log₂ pᵢ (in bits) of each sliding window of length
    /// <paramref name="windowSize"/>, advanced by <paramref name="stepSize"/>.
    /// Per-window entropy is delegated to <see cref="CalculateShannonEntropy"/>.
    /// </summary>
    /// <param name="sequence">Input sequence; symbol frequencies are taken over its letters (case-folded).</param>
    /// <param name="windowSize">Window width W in symbols (default 50). Must be ≤ sequence length for any window to be produced.</param>
    /// <param name="stepSize">Window advance in symbols (default 1).</param>
    /// <returns>
    /// One entropy value (bits) per window position, in order, for offsets
    /// 0, stepSize, 2·stepSize, … up to (length − windowSize); empty when the
    /// sequence is null/empty or <paramref name="windowSize"/> exceeds its length.
    /// </returns>
    /// <remarks>
    /// Shannon C. E. (1948), A Mathematical Theory of Communication, Bell Syst. Tech. J.
    /// 27(3):379–423. Base-2 logarithm yields bits; maximum is log₂k for k distinct symbols
    /// (2 bits for the 4-letter DNA alphabet).
    /// </remarks>
    public static IEnumerable<double> CalculateEntropyProfile(
        string sequence,
        int windowSize = 50,
        int stepSize = 1)
    {
        if (string.IsNullOrEmpty(sequence) || windowSize > sequence.Length)
            yield break;

        for (int i = 0; i <= sequence.Length - windowSize; i += stepSize)
        {
            string window = sequence.Substring(i, windowSize);
            yield return CalculateShannonEntropy(window);
        }
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
        IReadOnlyDictionary<char, int> Composition);

    /// <summary>
    /// Generates comprehensive summary statistics for a DNA/RNA sequence.
    /// </summary>
    public static SequenceSummary SummarizeNucleotideSequence(string? sequence)
    {
        // Treat null and empty identically (each per-metric method guards IsNullOrEmpty);
        // normalizing null to empty avoids dereferencing a null length when selecting the
        // Tm formula branch and keeps every per-metric call non-null.
        string seq = sequence ?? string.Empty;
        var comp = CalculateNucleotideComposition(seq);
        double entropy = CalculateShannonEntropy(seq);
        double complexity = CalculateLinguisticComplexity(seq);
        // Wallace rule applies to short oligos (< WallaceMaxLength = 14 A/C/G/T bases); the GC/Marmur-Doty
        // formula applies otherwise (SEQ-TM-001).
        // CalculateMeltingTemperature(useWallaceRule: true) itself selects the formula from the
        // number of A/C/G/T bases (OligoCalc), so N/gaps cannot push a short oligo into the GC formula.
        double tm = CalculateMeltingTemperature(seq, useWallaceRule: true);

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
            Composition: composition);
    }

    #endregion
}
