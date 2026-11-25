// SPDX-License-Identifier: GPL-2.0-or-later
// Port of Primer3 (thal.c ntthal, dpal.c, masker.c, libprimer3; primer3-org/primer3 and primer3-py 2.3.1; GNU General Public License v2 or later).
// Part of Seqeron.Genomics.OligoThermodynamics, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.OligoThermodynamics;

/// <summary>Thermodynamic (PrimerDesigner) members of the Primer3 ntthal family; see the licence notes of this assembly.</summary>
public static partial class ThermodynamicPrimerDesigner
{
    private static bool UsesStructureTerms(Primer3PenaltyWeights w, PrimerStructureScreen screen) => screen switch
    {
        PrimerStructureScreen.Primer3Alignment => w.SelfAny != 0 || w.SelfEnd != 0,
        PrimerStructureScreen.Primer3Thermodynamic => w.SelfAnyTh != 0 || w.SelfEndTh != 0 || w.HairpinTh != 0,
        _ => false,
    };

    private static StructureValues ComputeStructureValues(string seq, PrimerParameters param) => param.StructureScreen switch
    {
        PrimerStructureScreen.Primer3Alignment => new StructureValues(null,
            CalculatePrimerSelfAnyComplementarity(seq), CalculatePrimerSelfEndComplementarity(seq)),
        PrimerStructureScreen.Primer3Thermodynamic => new StructureValues(CalculatePrimer3OligoStructure(seq,
            param.EffectiveMonovalentMillimolar, param.EffectiveDivalentMillimolar, param.EffectiveDntpMillimolar,
            param.EffectiveDnaConcentrationNanomolar), null, null),
        _ => default,
    };

    /// <summary>
    /// Thermodynamic structure model of <see cref="FindMostStableHairpin(string, StructureModel, double?, double?, double?, double?, int?)"/>,
    /// <see cref="CalculateHairpinMeltingTemperature(string, StructureModel, double?, double?, double?, double?, int?)"/> and
    /// <see cref="FindMostStableDimer(string, string, StructureModel, double?, double?, double?, double?, double?, int?)"/>.
    /// </summary>
    public enum StructureModel
    {
        /// <summary>The default single-helix models: the hairpin single-stem core (SantaLucia &amp; Hicks 2004
        /// Table 1 stem stacks + Table 4 loop initiation, 1 M NaCl, no terminal-mismatch / special-loop /
        /// bulge / internal-loop terms) and the gapless contiguous Watson–Crick dimer scorer.</summary>
        SingleHelix = 0,

        /// <summary>The complete Primer3 <c>ntthal</c> (thal.c) dynamic program — SantaLucia &amp; Hicks (2004)
        /// stacks with terminal mismatches (<c>tstack</c>/<c>tstack2</c>), dangling ends, single mismatches,
        /// bulges, internal loops, size-keyed hairpin loops with the triloop/tetraloop bonus tables and the
        /// Mg²⁺/dNTP salt model — identical to primer3-py 2.3.1 <c>calc_hairpin</c> /
        /// <c>calc_homodimer</c> / <c>calc_heterodimer</c> (delegates to
        /// <see cref="CalculateHairpinStructureNtthal(string, double, double, double, double, int)"/> /
        /// <see cref="CalculateDimerStructureNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int)"/>).</summary>
        Ntthal = 1,
    }

    /// <summary>
    /// <see cref="FindMostStableHairpin(string, int, double)"/> with a selectable <paramref name="model"/>.
    /// <see cref="StructureModel.SingleHelix"/> is exactly <c>FindMostStableHairpin(sequence)</c> (the SantaLucia &amp;
    /// Hicks 2004 single-stem core at 1 M NaCl; every condition argument must then be <c>null</c>).
    /// <see cref="StructureModel.Ntthal"/> folds the oligo with the full thal.c hairpin DP (primer3-py
    /// <c>calc_hairpin</c>) and maps the optimal structure onto <see cref="HairpinResult"/>:
    /// <c>StemStart</c> = the 5′-most paired base, <c>StemEnd</c> = its partner, <c>StemLength</c> = base pairs of the
    /// stem closed by that pair (counted across bulges / internal loops), <c>LoopSize</c> = unpaired bases enclosed by
    /// its innermost pair, and ΔH°/ΔS°/ΔG = the ntthal values of the whole structure (kcal/mol, cal/(K·mol) incl. the
    /// ntthal salt term, kcal/mol at <paramref name="temperatureCelsius"/>) — exactly primer3-py <c>dh</c>/1000,
    /// <c>ds</c>, <c>dg</c>/1000. When ntthal places more than one hairpin on the oligo, the span / stem / loop fields
    /// describe the 5′-most one.
    /// </summary>
    /// <param name="sequence">DNA oligo (5′→3′), ACGT only (case-insensitive).</param>
    /// <param name="model">Structure model.</param>
    /// <param name="sodiumMolar">Ntthal only: mv, mol/L (null = 0.05, the calc_hairpin default).</param>
    /// <param name="divalentMolar">Ntthal only: dv, mol/L (null = 0.0015).</param>
    /// <param name="dntpMolar">Ntthal only: dNTP, mol/L (null = 0.0006).</param>
    /// <param name="temperatureCelsius">Ntthal only: temp_c (null = 37).</param>
    /// <param name="maxLoop">Ntthal only: max_loop 0–30 (null = 30).</param>
    /// <returns>The hairpin, or <c>null</c> for invalid input / no hairpin (ntthal <c>no_structure</c>).</returns>
    /// <exception cref="ArgumentException"><see cref="StructureModel.SingleHelix"/> with a non-null condition argument,
    /// or <see cref="StructureModel.Ntthal"/> with an oligo longer than 60 nt (thal.c <c>THAL_MAX_ALIGN</c>).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="model"/> undefined or max_loop outside 0–30.</exception>
    public static HairpinResult? FindMostStableHairpin(
        string sequence,
        StructureModel model,
        double? sodiumMolar = null,
        double? divalentMolar = null,
        double? dntpMolar = null,
        double? temperatureCelsius = null,
        int? maxLoop = null)
    {
        if (!UseNtthal(model, sodiumMolar, divalentMolar, dntpMolar, null, temperatureCelsius, maxLoop))
            return PrimerDesigner.FindMostStableHairpin(sequence);

        var r = RunHairpinNtthal(sequence, sodiumMolar ?? 0.05, divalentMolar ?? 0.0015, dntpMolar ?? 0.0006,
            temperatureCelsius ?? NtthalDefaultTemperatureCelsius, maxLoop ?? NtthalDefaultMaxLoop,
            withStructure: false, NtthalMaxAlignLength);
        if (r is null)
            return null;
        var v = r.Value;
        int[] bp = v.PairPartners!;
        int first = Array.FindIndex(bp, p => p > 0); // 0-based 5'-most paired base
        int a = first, b = bp[first] - 1, stem = 1;
        while (true)
        {
            // Next pair nested inside (a, b): the first paired base after a (a hairpin stem has no branches).
            int k = a + 1;
            while (k < b && bp[k] == 0) k++;
            if (k >= b) break;
            a = k;
            b = bp[k] - 1;
            stem++;
        }
        return new HairpinResult(first, bp[first] - 1, stem, b - a - 1,
            v.DeltaH / 1000.0, v.DeltaS, v.DeltaG37 / 1000.0);
    }

    /// <summary>
    /// <see cref="CalculateHairpinMeltingTemperature(string, int, double)"/> with a selectable <paramref name="model"/>:
    /// <see cref="StructureModel.SingleHelix"/> is exactly <c>CalculateHairpinMeltingTemperature(sequence)</c> (every
    /// condition argument must be <c>null</c>); <see cref="StructureModel.Ntthal"/> returns the ntthal hairpin Tm
    /// (primer3-py <c>calc_hairpin(...).tm</c>, <see cref="CalculateHairpinThermodynamicsNtthal(string, double, double, double, double, int)"/>).
    /// Parameters as in <see cref="FindMostStableHairpin(string, StructureModel, double?, double?, double?, double?, int?)"/>.
    /// </summary>
    /// <returns>The hairpin Tm in °C, or <c>double.NaN</c> for invalid input / no hairpin.</returns>
    /// <exception cref="ArgumentException">As <see cref="FindMostStableHairpin(string, StructureModel, double?, double?, double?, double?, int?)"/>.</exception>
    /// <exception cref="ArgumentOutOfRangeException">As <see cref="FindMostStableHairpin(string, StructureModel, double?, double?, double?, double?, int?)"/>.</exception>
    public static double CalculateHairpinMeltingTemperature(
        string sequence,
        StructureModel model,
        double? sodiumMolar = null,
        double? divalentMolar = null,
        double? dntpMolar = null,
        double? temperatureCelsius = null,
        int? maxLoop = null)
    {
        if (!UseNtthal(model, sodiumMolar, divalentMolar, dntpMolar, null, temperatureCelsius, maxLoop))
            return PrimerDesigner.CalculateHairpinMeltingTemperature(sequence);

        var r = RunHairpinNtthal(sequence, sodiumMolar ?? 0.05, divalentMolar ?? 0.0015, dntpMolar ?? 0.0006,
            temperatureCelsius ?? NtthalDefaultTemperatureCelsius, maxLoop ?? NtthalDefaultMaxLoop,
            withStructure: false, NtthalMaxAlignLength);
        return r is null ? double.NaN : r.Value.TmCelsius;
    }

    // Validates a StructureModel and its ntthal-only condition arguments: true = Ntthal; false = SingleHelix, which
    // takes no ntthal condition (its own conditions, if any, are passed separately by the caller).
    private static bool UseNtthal(StructureModel model, double? mv, double? dv, double? dntp, double? dnaConc,
        double? tempC, int? maxLoop)
    {
        switch (model)
        {
            case StructureModel.Ntthal:
                return true;
            case StructureModel.SingleHelix:
                if (mv is not null || dv is not null || dntp is not null || dnaConc is not null || tempC is not null || maxLoop is not null)
                    throw new ArgumentException(
                        "The SingleHelix model has fixed conditions (hairpin: 1 M NaCl, 37 °C; dimer: monovalent salt and " +
                        "strand concentration only); ntthal condition arguments require StructureModel.Ntthal.");
                return false;
            default:
                throw new ArgumentOutOfRangeException(nameof(model), model, "Unknown StructureModel.");
        }
    }

    /// <summary>
    /// <see cref="FindMostStableDimer(string, string, double, double)"/> with a selectable <paramref name="model"/>.
    /// <see cref="StructureModel.SingleHelix"/> is exactly <c>FindMostStableDimer(strand1, strand2, sodiumMolar,
    /// strandConcentrationMolar)</c> (the gapless contiguous Watson–Crick scorer; null = its defaults 50 mM / 50 nM;
    /// <paramref name="divalentMolar"/>, <paramref name="dntpMolar"/>, <paramref name="temperatureCelsius"/> and
    /// <paramref name="maxLoop"/> must be <c>null</c>). <see cref="StructureModel.Ntthal"/> runs the full thal.c dimer
    /// DP in mode ANY (primer3-py <c>calc_heterodimer(strand1, strand2, …)</c>; a self-dimer —
    /// <c>calc_homodimer</c> — by passing the same sequence twice) and maps its optimal structure onto
    /// <see cref="DimerResult"/>: <c>Strand1Start</c> / <c>Strand2Start</c> = the 0-based 5′-most paired base on each
    /// strand, <c>BasePairs</c> = paired bases of the structure (ntthal reports single-pair structures too), and
    /// ΔH°/ΔS°/ΔG = the ntthal values (kcal/mol; cal/(K·mol) incl. the ntthal salt term; kcal/mol at
    /// <paramref name="temperatureCelsius"/>) — exactly primer3-py <c>dh</c>/1000, <c>ds</c>, <c>dg</c>/1000. The
    /// structure may contain mismatches, bulges, internal loops and dangling ends; its Tm is
    /// <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int)"/>.
    /// </summary>
    /// <param name="strand1">First DNA oligo (5′→3′), ACGT only (case-insensitive).</param>
    /// <param name="strand2">Second DNA oligo (5′→3′); the same string for a self-dimer.</param>
    /// <param name="model">Structure model.</param>
    /// <param name="sodiumMolar">mv, mol/L (null = 0.05 for both models).</param>
    /// <param name="divalentMolar">Ntthal only: dv, mol/L (null = 0.0015, the primer3-py default).</param>
    /// <param name="dntpMolar">Ntthal only: dNTP, mol/L (null = 0.0006).</param>
    /// <param name="strandConcentrationMolar">Oligo concentration, mol/L (null = 50 nM for both models).</param>
    /// <param name="temperatureCelsius">Ntthal only: temp_c (null = 37).</param>
    /// <param name="maxLoop">Ntthal only: max_loop 0–30 (null = 30).</param>
    /// <returns>The dimer, or <c>null</c> for invalid input / no duplex (ntthal <c>no_structure</c>).</returns>
    /// <exception cref="ArgumentException"><see cref="StructureModel.SingleHelix"/> with a non-null ntthal-only argument,
    /// or <see cref="StructureModel.Ntthal"/> with both strands longer than 60 nt / either longer than 10 000 nt
    /// (thal.c <c>THAL_MAX_ALIGN</c> / <c>THAL_MAX_SEQ</c>).</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="model"/> undefined or max_loop outside 0–30.</exception>
    public static DimerResult? FindMostStableDimer(
        string strand1,
        string strand2,
        StructureModel model,
        double? sodiumMolar = null,
        double? divalentMolar = null,
        double? dntpMolar = null,
        double? strandConcentrationMolar = null,
        double? temperatureCelsius = null,
        int? maxLoop = null)
    {
        double mv = sodiumMolar ?? ThermoConstants.DefaultNaConcentration;
        double dnaConc = strandConcentrationMolar ?? DefaultDimerStrandConcentrationMolar;
        if (!UseNtthal(model, null, divalentMolar, dntpMolar, null, temperatureCelsius, maxLoop))
            return PrimerDesigner.FindMostStableDimer(strand1, strand2, mv, dnaConc);

        var r = RunDimerNtthal(strand1, strand2, NtthalAlignmentMode.Any, mv, divalentMolar ?? 0.0015,
            dntpMolar ?? 0.0006, dnaConc, temperatureCelsius ?? NtthalDefaultTemperatureCelsius,
            maxLoop ?? NtthalDefaultMaxLoop, withStructure: false, NtthalMaxAlignLength);
        if (r is null)
            return null;
        var v = r.Value;
        // thal.c pairs strand 1 (1-based i) with the REVERSED strand 2 (1-based j); the traceback starts at the
        // 3'-most strand-1 pair (Strand1End, Strand2End), whose partner is the 5'-most paired base of strand 2.
        return new DimerResult(v.Strand1Start - 1, strand2.Length - v.Strand2End, v.BasePairs,
            v.DeltaH / 1000.0, v.DeltaS, v.DeltaG37 / 1000.0);
    }

    /// <summary>
    /// Computes the intermolecular <b>dimer</b> melting temperature (°C) between two DNA oligos
    /// (a self-dimer when both arguments are the same sequence; a hetero/cross-dimer otherwise),
    /// as the bimolecular Tm of the most stable duplex found by the Primer3 / <c>ntthal</c>
    /// thermodynamic alignment (SantaLucia &amp; Hicks 2004 unified NN):
    /// Tm = ΔH°·1000/(ΔS° + R·ln(C_T/x)) − 273.15, x = 1 if both oligos are reverse-complement
    /// palindromes else x = 4 (C_T = the total strand concentration, default 50 nM). <b>Opt-in</b>: existing Tm methods are unchanged.
    /// </summary>
    /// <param name="strand1">First DNA oligo (5'→3'); ACGT only (case-insensitive).</param>
    /// <param name="strand2">Second DNA oligo (5'→3'); ACGT only (case-insensitive).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The dimer Tm in °C (ntthal also reports weak / single-pair structures, whose Tm may
    /// be negative), or <c>double.NaN</c> if ntthal finds no structure / either sequence is null,
    /// empty or non-ACGT.</returns>
    /// <exception cref="ArgumentException">Both strands are longer than 60 nt (thal.c THAL_MAX_ALIGN).</exception>
    public static double CalculateDimerMeltingTemperature(
        string strand1,
        string strand2,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar)
    {
        var dimer = CalculateDimerThermodynamicsNtthal(strand1, strand2, sodiumMolar, strandConcentrationMolar);
        return dimer is null ? double.NaN : dimer.Value.TmCelsius;
    }

    /// <summary>
    /// Full <c>ntthal</c> dimer thermodynamics (ΔH°, ΔS°, ΔG°37 and bimolecular Tm) for the most
    /// stable intermolecular DNA duplex between two oligos, computed by the complete Primer3
    /// <c>ntthal</c> dynamic program (mode ANY): matched nearest-neighbour stacks, single internal
    /// mismatches, internal loops, single- and multi-base bulges, and terminal overhangs /
    /// dangling ends (the <c>tstack2</c> terminal table + 5′/3′ dangling-end tables + interior /
    /// bulge loop-length parameters). Unlike <see cref="FindMostStableDimer"/> (which scores only
    /// the best contiguous Watson–Crick run), this reproduces primer3-py's
    /// <c>calc_homodimer</c>/<c>calc_heterodimer</c> for dimers whose optimum is <b>non-contiguous</b>.
    /// This overload uses monovalent salt only (primer3-py with <c>dv_conc=0, dntp_conc=0</c>; note the
    /// primer3-py defaults are 1.5 mM Mg²⁺ / 0.6 mM dNTP — use the overload taking
    /// <see cref="NtthalAlignmentMode"/> for those).
    /// <b>Opt-in</b>: all other Tm methods and defaults are unchanged.
    /// </summary>
    /// <exception cref="ArgumentException">Both strands are longer than 60 nt (thal.c THAL_MAX_ALIGN).</exception>
    /// <param name="strand1">First DNA oligo (5′→3′); ≥ 1 ACGT base.</param>
    /// <param name="strand2">Second DNA oligo (5′→3′); the same string for a self-dimer.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The most stable dimer's thermodynamics, or <c>null</c> if either strand is
    /// null/empty/contains a non-ACGT character, or no duplex can be formed (ntthal
    /// <c>no_structure</c>).</returns>
    public static DimerThermodynamics? CalculateDimerThermodynamicsNtthal(
        string strand1,
        string strand2,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar)
    {
        if (string.IsNullOrEmpty(strand1) || string.IsNullOrEmpty(strand2))
            return null;
        string s1 = strand1.ToUpperInvariant();
        string s2 = strand2.ToUpperInvariant();
        foreach (char c in s1)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;
        foreach (char c in s2)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;

        var r = NtthalDimer.Run(s1, s2, sodiumMolar, strandConcentrationMolar);
        if (r is null)
            return null;

        var v = r.Value;
        // Convert ntthal native cal/mol → the library's kcal/mol convention for ΔH/ΔG.
        return new DimerThermodynamics(
            DeltaH: v.DeltaH / 1000.0,
            DeltaS: v.DeltaS,
            DeltaG37: v.DeltaG37 / 1000.0,
            TmCelsius: v.TmCelsius,
            BasePairs: v.BasePairs);
    }

    /// <summary>
    /// Computes the full <b>ntthal</b> intramolecular-hairpin thermodynamics (ΔH°, ΔS°, ΔG°37, Tm)
    /// of a DNA oligo, reproducing primer3-py's <c>calc_hairpin</c>. This runs the complete
    /// Primer3 <c>ntthal</c> monomer dynamic program (a single stem with internal
    /// mismatches/loops, bulges, terminal mismatch / dangling-end terminal contributions and the
    /// size-keyed hairpin-loop initiation) and — unlike <see cref="FindMostStableHairpin"/>'s
    /// SantaLucia &amp; Hicks (2004) Table 4 model — <b>automatically applies the bundled
    /// sequence-specific special triloop / tetraloop stability bonuses</b> (the primer3
    /// <c>triloop.dh/.ds</c> + <c>tetraloop.dh/.ds</c> tables, keyed on the full loop string
    /// including the closing base pair). No caller-supplied loop bonus is required.
    /// <b>Opt-in</b>: <see cref="FindMostStableHairpin"/>, the duplex/dimer Tm methods and all
    /// defaults are unchanged.
    /// </summary>
    /// <param name="sequence">The DNA oligo (5′→3'); ≥ 1 ACGT base.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM, the
    /// primer3 <c>calc_hairpin</c> default of <c>mv=50</c>). This overload has no divalent cations or
    /// dNTPs: it equals <c>calc_hairpin(seq, mv_conc, dv_conc=0, dntp_conc=0)</c> at 37 °C, max loop 30.</param>
    /// <returns>The most stable hairpin's thermodynamics, or <c>null</c> if the sequence is
    /// null/empty/contains a non-ACGT character, or no hairpin can form (ntthal
    /// <c>no_structure</c>, e.g. a homopolymer).</returns>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt (thal.c <c>THAL_MAX_ALIGN</c>).</exception>
    public static HairpinThermodynamics? CalculateHairpinThermodynamicsNtthal(
        string sequence,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration)
    {
        if (string.IsNullOrEmpty(sequence))
            return null;
        string seq = sequence.ToUpperInvariant();
        foreach (char c in seq)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;

        var r = NtthalHairpin.Run(seq, sodiumMolar);
        if (r is null)
            return null;

        var v = r.Value;
        // Convert ntthal native cal/mol → the library's kcal/mol convention for ΔH/ΔG.
        return new HairpinThermodynamics(
            DeltaH: v.DeltaH / 1000.0,
            DeltaS: v.DeltaS,
            DeltaG37: v.DeltaG37 / 1000.0,
            TmCelsius: v.TmCelsius,
            BasePairs: v.BasePairs);
    }

    /// <summary>
    /// ntthal dimer alignment type (Primer3 <c>thal_alignment_type</c>, ntthal <c>-a</c>).
    /// </summary>
    public enum NtthalAlignmentMode
    {
        /// <summary>THAL_ANY: the most stable duplex anywhere (primer3-py <c>calc_heterodimer</c>).</summary>
        Any,
        /// <summary>THAL_END1: the duplex must contain the 3′-terminal base of strand 1
        /// (primer3-py <c>calc_end_stability(strand1, strand2)</c>).</summary>
        End1,
        /// <summary>THAL_END2: the duplex must contain the 3′-terminal base of strand 2
        /// (= END1 with the strands swapped).</summary>
        End2,
    }

    /// <summary>
    /// Full <c>ntthal</c> dimer thermodynamics with an explicit alignment type and the complete
    /// ntthal salt model (<c>saltCorrectS</c>: 0.368·ln((mv + 120·√max(0, dv − dntp))/1000), mM),
    /// reproducing primer3-py <c>calc_heterodimer</c> (mode <see cref="NtthalAlignmentMode.Any"/>)
    /// and <c>calc_end_stability</c> (mode <see cref="NtthalAlignmentMode.End1"/>) at any
    /// mv/dv/dntp/dna_conc (temperature 37 °C, max loop 30 — the primer3-py defaults).
    /// </summary>
    /// <param name="strand1">First DNA oligo (5′→3′), ACGT only.</param>
    /// <param name="strand2">Second DNA oligo (5′→3′), ACGT only.</param>
    /// <param name="mode">ntthal alignment type.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration, mol/L.</param>
    /// <param name="divalentMolar">Mg²⁺ concentration, mol/L.</param>
    /// <param name="dntpMolar">dNTP concentration, mol/L.</param>
    /// <param name="strandConcentrationMolar">Oligo concentration, mol/L (ntthal dna_conc).</param>
    /// <returns>The thermodynamics, or <c>null</c> for invalid input or when no duplex forms.</returns>
    /// <exception cref="ArgumentException">Both strands are longer than 60 nt, or either is longer
    /// than 10 000 nt (thal.c <c>THAL_MAX_ALIGN</c> / <c>THAL_MAX_SEQ</c>; primer3-py raises).</exception>
    public static DimerThermodynamics? CalculateDimerThermodynamicsNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double strandConcentrationMolar) =>
        CalculateDimerThermodynamicsNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, NtthalDefaultTemperatureCelsius, NtthalDefaultMaxLoop);

    /// <summary>primer3-py / ntthal default analysis temperature (°C) at which ΔG is reported.</summary>
    public const double NtthalDefaultTemperatureCelsius = 37.0;

    /// <summary>primer3-py / ntthal default (and maximum) internal-loop / bulge size.</summary>
    public const int NtthalDefaultMaxLoop = 30;

    /// <summary>
    /// thal.h <c>THAL_MAX_ALIGN</c> = 60: Primer3 / primer3-py refuse a dimer whose two strands are both longer, and a
    /// hairpin of a longer oligo ("At least one sequence must be equal to or shorter than 60bp …"). It is a
    /// compile-time constant of thal.c (<c>#ifndef THAL_MAX_ALIGN</c>, used only by <c>thal_check_errors</c>; the DP
    /// tables are allocated from the actual lengths) chosen by the Primer3 authors as "the maximum reasonable length
    /// for nearest neighbor models … only two states of melting" (thal.h). The ntthal overloads taking
    /// <c>maxAlignLength</c> raise it (opt-in) and then reproduce ntthal built with <c>-DTHAL_MAX_ALIGN=…</c>.
    /// </summary>
    public const int NtthalMaxAlignLength = 60;

    /// <summary>thal.h <c>THAL_MAX_SEQ</c> = 10 000: maximum length of either strand (and the largest
    /// <c>maxAlignLength</c> accepted).</summary>
    public const int NtthalMaxSequenceLength = 10000;

    /// <summary>
    /// Full <c>ntthal</c> dimer thermodynamics with every primer3-py <c>calc_heterodimer</c> /
    /// <c>calc_end_stability</c> argument: alignment type, mv/dv/dntp/dna_conc, <c>temp_c</c> (the
    /// temperature at which ΔG is evaluated: ΔG = ΔH − (temp_c + 273.15)·ΔS, thal.c <c>calcDimer</c>;
    /// the DP ranking and Tm do not depend on it) and <c>max_loop</c> (largest internal loop / bulge
    /// considered, 0–30). Bit-faithful port of primer3-py 2.3.1 <c>thal.c</c>.
    /// </summary>
    /// <param name="strand1">First DNA oligo (5′→3′), ACGT only.</param>
    /// <param name="strand2">Second DNA oligo (5′→3′), ACGT only.</param>
    /// <param name="mode">ntthal alignment type.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration, mol/L.</param>
    /// <param name="divalentMolar">Mg²⁺ concentration, mol/L.</param>
    /// <param name="dntpMolar">dNTP concentration, mol/L.</param>
    /// <param name="strandConcentrationMolar">Oligo concentration, mol/L (ntthal dna_conc).</param>
    /// <param name="temperatureCelsius">primer3-py <c>temp_c</c>; the returned
    /// <see cref="DimerThermodynamics.DeltaG37"/> is ΔG at this temperature.</param>
    /// <param name="maxLoop">primer3-py <c>max_loop</c>, 0–30.</param>
    /// <returns>The thermodynamics, or <c>null</c> for invalid input or when no duplex forms.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30.</exception>
    /// <exception cref="ArgumentException">Both strands are longer than 60 nt, or either is longer
    /// than 10 000 nt (thal.c <c>THAL_MAX_ALIGN</c> / <c>THAL_MAX_SEQ</c>).</exception>
    public static DimerThermodynamics? CalculateDimerThermodynamicsNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double strandConcentrationMolar,
        double temperatureCelsius,
        int maxLoop) =>
        CalculateDimerStructureNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: false,
            NtthalMaxAlignLength)?.Thermodynamics;

    /// <summary>
    /// As <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int)"/>
    /// with thal.h <c>THAL_MAX_ALIGN</c> raised to <paramref name="maxAlignLength"/> (opt-in; see
    /// <see cref="NtthalMaxAlignLength"/>): the unchanged ntthal recursions run when at least one strand is at most
    /// <paramref name="maxAlignLength"/> nt — identical to ntthal / thal.c compiled with
    /// <c>-DTHAL_MAX_ALIGN=maxAlignLength</c>. With 60 it is the 9-argument overload. Cost is O(len1·len2·maxLoop²).
    /// </summary>
    /// <param name="maxAlignLength">THAL_MAX_ALIGN, <see cref="NtthalMaxAlignLength"/> (60) …
    /// <see cref="NtthalMaxSequenceLength"/> (10 000) — smaller values would refuse what Primer3 accepts.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30 or
    /// <paramref name="maxAlignLength"/> outside 60–10 000.</exception>
    /// <exception cref="ArgumentException">Both strands are longer than <paramref name="maxAlignLength"/>, or either
    /// is longer than 10 000 nt.</exception>
    public static DimerThermodynamics? CalculateDimerThermodynamicsNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double strandConcentrationMolar,
        double temperatureCelsius,
        int maxLoop,
        int maxAlignLength) =>
        CalculateDimerStructureNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: false,
            CheckedMaxAlignLength(maxAlignLength))?.Thermodynamics;

    // Public THAL_MAX_ALIGN override range: 60 (Primer3) … 10 000 (THAL_MAX_SEQ).
    private static int CheckedMaxAlignLength(int maxAlignLength)
    {
        if (maxAlignLength < NtthalMaxAlignLength || maxAlignLength > NtthalMaxSequenceLength)
            throw new ArgumentOutOfRangeException(nameof(maxAlignLength), maxAlignLength,
                "THAL_MAX_ALIGN override must be in 60..10000 (Primer3 default .. THAL_MAX_SEQ).");
        return maxAlignLength;
    }

    /// <summary>
    /// As <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int)"/>,
    /// plus the optimal duplex drawn exactly as thal.c <c>drawDimer</c> / primer3-py
    /// <c>ThermoResult.ascii_structure_lines</c> (<c>output_structure=True</c>): four lines
    /// "SEQ	…" (unpaired strand-1 bases), "SEQ	…" (paired strand-1 bases), "STR	…" (paired
    /// strand-2 bases), "STR	…" (unpaired strand-2 bases); strand 2 runs 3′→5′ and '-' pads the
    /// shorter side of a loop. In mode <see cref="NtthalAlignmentMode.End2"/> the strands are drawn
    /// swapped, as ntthal does.
    /// </summary>
    /// <returns>The thermodynamics and structure lines, or <c>null</c> for invalid input or when no
    /// duplex forms.</returns>
    public static NtthalDimerStructure? CalculateDimerStructureNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode = NtthalAlignmentMode.Any,
        double sodiumMolar = 0.05,
        double divalentMolar = 0.0015,
        double dntpMolar = 0.0006,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar,
        double temperatureCelsius = NtthalDefaultTemperatureCelsius,
        int maxLoop = NtthalDefaultMaxLoop) =>
        CalculateDimerStructureNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: true, NtthalMaxAlignLength);

    /// <summary>
    /// As <see cref="CalculateDimerStructureNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int)"/>
    /// with THAL_MAX_ALIGN raised to <paramref name="maxAlignLength"/> (opt-in, 60–10 000; see
    /// <see cref="CalculateDimerThermodynamicsNtthal(string, string, NtthalAlignmentMode, double, double, double, double, double, int, int)"/>).
    /// </summary>
    public static NtthalDimerStructure? CalculateDimerStructureNtthal(
        string strand1,
        string strand2,
        NtthalAlignmentMode mode,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double strandConcentrationMolar,
        double temperatureCelsius,
        int maxLoop,
        int maxAlignLength) =>
        CalculateDimerStructureNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: true,
            CheckedMaxAlignLength(maxAlignLength));

    private static NtthalDimerStructure? CalculateDimerStructureNtthal(
        string strand1, string strand2, NtthalAlignmentMode mode, double sodiumMolar, double divalentMolar,
        double dntpMolar, double strandConcentrationMolar, double temperatureCelsius, int maxLoop, bool withStructure,
        int maxAlign)
    {
        var r = RunDimerNtthal(strand1, strand2, mode, sodiumMolar, divalentMolar, dntpMolar,
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure, maxAlign);
        if (r is null)
            return null;
        var v = r.Value;
        return new NtthalDimerStructure(
            new DimerThermodynamics(v.DeltaH / 1000.0, v.DeltaS, v.DeltaG37 / 1000.0, v.TmCelsius, v.BasePairs),
            v.AsciiStructure ?? Array.Empty<string>());
    }

    // The one entry into the ntthal dimer engine (ACGT check, upper-casing, mode mapping, °C → K).
    private static NtthalDimer.Result? RunDimerNtthal(
        string strand1, string strand2, NtthalAlignmentMode mode, double sodiumMolar, double divalentMolar,
        double dntpMolar, double strandConcentrationMolar, double temperatureCelsius, int maxLoop, bool withStructure,
        int maxAlign)
    {
        if (!IsAcgtOnly(strand1) || !IsAcgtOnly(strand2))
            return null;
        var type = mode switch
        {
            NtthalAlignmentMode.Any => NtthalDimer.AlignmentType.Any,
            NtthalAlignmentMode.End1 => NtthalDimer.AlignmentType.End1,
            NtthalAlignmentMode.End2 => NtthalDimer.AlignmentType.End2,
            _ => throw new ArgumentOutOfRangeException(nameof(mode)),
        };
        return NtthalDimer.Run(strand1.ToUpperInvariant(), strand2.ToUpperInvariant(),
            sodiumMolar, strandConcentrationMolar, type, divalentMolar, dntpMolar,
            temperatureCelsius + KelvinOffset, maxLoop, withStructure, maxAlign);
    }

    /// <summary>
    /// ntthal dimer thermodynamics plus the thal.c <c>drawDimer</c> ASCII duplex
    /// (primer3-py <c>ThermoResult.ascii_structure_lines</c>).
    /// </summary>
    /// <param name="Thermodynamics">ΔH/ΔS/ΔG/Tm of the optimal duplex.</param>
    /// <param name="AsciiStructureLines">The four "SEQ	"/"SEQ	"/"STR	"/"STR	" lines.</param>
    public sealed record NtthalDimerStructure(DimerThermodynamics Thermodynamics, IReadOnlyList<string> AsciiStructureLines);

    /// <summary>
    /// Full <c>ntthal</c> hairpin thermodynamics with the complete ntthal salt model (divalent
    /// cations and dNTPs enter through <c>saltCorrectS</c>), reproducing primer3-py
    /// <c>calc_hairpin</c> at any mv/dv/dntp (temperature 37 °C, max loop 30 — the primer3-py defaults).
    /// </summary>
    /// <param name="sequence">DNA oligo (5′→3′), ACGT only, at most 60 nt.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration, mol/L.</param>
    /// <param name="divalentMolar">Mg²⁺ concentration, mol/L.</param>
    /// <param name="dntpMolar">dNTP concentration, mol/L.</param>
    /// <returns>The thermodynamics, or <c>null</c> for invalid input or when no hairpin forms.</returns>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt (thal.c
    /// <c>THAL_MAX_ALIGN</c>; primer3-py raises).</exception>
    public static HairpinThermodynamics? CalculateHairpinThermodynamicsNtthal(
        string sequence,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar) =>
        CalculateHairpinThermodynamicsNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            NtthalDefaultTemperatureCelsius, NtthalDefaultMaxLoop);

    /// <summary>
    /// Full <c>ntthal</c> hairpin thermodynamics with every primer3-py <c>calc_hairpin</c> argument:
    /// mv/dv/dntp, <c>temp_c</c> (the analysis temperature: it sets the reported ΔG = ΔH − (temp_c +
    /// 273.15)·ΔS and, as in thal.c <c>calc_terminal_bp</c>, the exterior-loop acceptance test
    /// ΔH − T·ΔS &lt; 0, so it can change the selected structure) and <c>max_loop</c> (largest internal
    /// loop / bulge considered, 0–30). Bit-faithful port of primer3-py 2.3.1 <c>thal.c</c> (type 4).
    /// </summary>
    /// <param name="sequence">DNA oligo (5′→3′), ACGT only (case-insensitive), at most 60 nt.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration, mol/L.</param>
    /// <param name="divalentMolar">Mg²⁺ concentration, mol/L.</param>
    /// <param name="dntpMolar">dNTP concentration, mol/L.</param>
    /// <param name="temperatureCelsius">primer3-py <c>temp_c</c>; the returned
    /// <see cref="HairpinThermodynamics.DeltaG37"/> is ΔG at this temperature.</param>
    /// <param name="maxLoop">primer3-py <c>max_loop</c>, 0–30.</param>
    /// <returns>The thermodynamics, or <c>null</c> for invalid input or when no hairpin forms.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30.</exception>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt.</exception>
    public static HairpinThermodynamics? CalculateHairpinThermodynamicsNtthal(
        string sequence,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double temperatureCelsius,
        int maxLoop) =>
        CalculateHairpinStructureNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure: false, NtthalMaxAlignLength)?.Thermodynamics;

    /// <summary>
    /// As <see cref="CalculateHairpinThermodynamicsNtthal(string, double, double, double, double, int)"/> with thal.h
    /// <c>THAL_MAX_ALIGN</c> raised to <paramref name="maxAlignLength"/> (opt-in; see <see cref="NtthalMaxAlignLength"/>):
    /// oligos up to <paramref name="maxAlignLength"/> nt are folded by the unchanged ntthal hairpin recursions —
    /// identical to ntthal / thal.c compiled with <c>-DTHAL_MAX_ALIGN=maxAlignLength</c> (verified on 61–120-mers
    /// against thal.c from primer3-py 2.3.1 built with a larger THAL_MAX_ALIGN). With 60 it is the 6-argument
    /// overload. Cost is O(n²·maxLoop²). Note that the Primer3 authors chose 60 as the limit of the two-state
    /// nearest-neighbour model (thal.h); longer oligos are folded with the same single-structure model.
    /// </summary>
    /// <param name="maxAlignLength">THAL_MAX_ALIGN, <see cref="NtthalMaxAlignLength"/> (60) …
    /// <see cref="NtthalMaxSequenceLength"/> (10 000).</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30 or
    /// <paramref name="maxAlignLength"/> outside 60–10 000.</exception>
    /// <exception cref="ArgumentException">The sequence is longer than <paramref name="maxAlignLength"/>.</exception>
    public static HairpinThermodynamics? CalculateHairpinThermodynamicsNtthal(
        string sequence,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double temperatureCelsius,
        int maxLoop,
        int maxAlignLength) =>
        CalculateHairpinStructureNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure: false, CheckedMaxAlignLength(maxAlignLength))?.Thermodynamics;

    /// <summary>
    /// As <see cref="CalculateHairpinThermodynamicsNtthal(string, double, double, double, double, int)"/>,
    /// plus the optimal hairpin drawn exactly as thal.c <c>drawHairpin</c> / primer3-py
    /// <c>ThermoResult.ascii_structure_lines</c> (<c>output_structure=True</c>): two lines,
    /// "SEQ\t" followed by one character per base ('-' unpaired; for each base pair the 5′ partner
    /// is drawn '/' and the 3′ partner '\') and "STR\t" followed by the (upper-case) oligo.
    /// Defaults are the primer3-py <c>calc_hairpin</c> defaults (mv 50 mM, dv 1.5 mM, dNTP 0.6 mM,
    /// 37 °C, max loop 30).
    /// </summary>
    /// <returns>The thermodynamics and structure lines, or <c>null</c> for invalid input or when no
    /// hairpin forms.</returns>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxLoop"/> outside 0–30.</exception>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt.</exception>
    public static NtthalHairpinStructure? CalculateHairpinStructureNtthal(
        string sequence,
        double sodiumMolar = 0.05,
        double divalentMolar = 0.0015,
        double dntpMolar = 0.0006,
        double temperatureCelsius = NtthalDefaultTemperatureCelsius,
        int maxLoop = NtthalDefaultMaxLoop) =>
        CalculateHairpinStructureNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure: true, NtthalMaxAlignLength);

    /// <summary>
    /// As <see cref="CalculateHairpinStructureNtthal(string, double, double, double, double, int)"/> with THAL_MAX_ALIGN
    /// raised to <paramref name="maxAlignLength"/> (opt-in, 60–10 000; see
    /// <see cref="CalculateHairpinThermodynamicsNtthal(string, double, double, double, double, int, int)"/>).
    /// </summary>
    public static NtthalHairpinStructure? CalculateHairpinStructureNtthal(
        string sequence,
        double sodiumMolar,
        double divalentMolar,
        double dntpMolar,
        double temperatureCelsius,
        int maxLoop,
        int maxAlignLength) =>
        CalculateHairpinStructureNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure: true, CheckedMaxAlignLength(maxAlignLength));

    private static NtthalHairpinStructure? CalculateHairpinStructureNtthal(
        string sequence, double sodiumMolar, double divalentMolar, double dntpMolar,
        double temperatureCelsius, int maxLoop, bool withStructure, int maxAlign)
    {
        var r = RunHairpinNtthal(sequence, sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius, maxLoop, withStructure, maxAlign);
        if (r is null)
            return null;
        var v = r.Value;
        return new NtthalHairpinStructure(
            new HairpinThermodynamics(v.DeltaH / 1000.0, v.DeltaS, v.DeltaG37 / 1000.0, v.TmCelsius, v.BasePairs),
            v.AsciiStructure ?? Array.Empty<string>());
    }

    // The one entry into the ntthal hairpin engine (ACGT check, upper-casing, °C → K).
    private static NtthalHairpin.Result? RunHairpinNtthal(
        string sequence, double sodiumMolar, double divalentMolar, double dntpMolar,
        double temperatureCelsius, int maxLoop, bool withStructure, int maxAlign) =>
        IsAcgtOnly(sequence)
            ? NtthalHairpin.Run(sequence.ToUpperInvariant(), sodiumMolar, divalentMolar, dntpMolar,
                temperatureCelsius + KelvinOffset, maxLoop, withStructure, maxAlign)
            : null;

    /// <summary>
    /// ntthal hairpin thermodynamics plus the thal.c <c>drawHairpin</c> ASCII structure
    /// (primer3-py <c>ThermoResult.ascii_structure_lines</c>).
    /// </summary>
    /// <param name="Thermodynamics">ΔH/ΔS/ΔG/Tm of the optimal hairpin.</param>
    /// <param name="AsciiStructureLines">The "SEQ\t…" and "STR\t…" lines.</param>
    public sealed record NtthalHairpinStructure(HairpinThermodynamics Thermodynamics, IReadOnlyList<string> AsciiStructureLines);

    /// <summary>
    /// Computes a primer's Primer3 thermodynamic secondary-structure Tm values exactly as Primer3
    /// (default <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1</c>) does in <c>oligo_compl_thermod</c> /
    /// <c>oligo_hairpin</c>: self_any = ntthal ANY(primer, primer), self_end = ntthal END1(primer,
    /// primer), hairpin = ntthal HAIRPIN(primer), each the Tm (°C) of the most stable structure, 0
    /// when none forms (ntthal <c>no_structure</c>) or the Tm is negative (<c>align_thermod</c>). Conditions default to Primer3's primer
    /// conditions (50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo). Cross-checked against
    /// primer3-py 2.3.1 <c>calc_homodimer</c>/<c>calc_end_stability</c>/<c>calc_hairpin</c> and the
    /// <c>PRIMER_LEFT_0_SELF_ANY_TH</c>/<c>_SELF_END_TH</c>/<c>_HAIRPIN_TH</c> values of
    /// <c>design_primers</c>.
    /// </summary>
    /// <param name="primer">Primer (5′→3′), case-insensitive, ACGT only.</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM.</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM.</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM.</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM.</param>
    /// <returns>The three Tm values, or <c>null</c> when the primer is null/empty or contains a
    /// non-ACGT character.</returns>
    public static Primer3OligoStructure? CalculatePrimer3OligoStructure(
        string primer,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar,
        double dnaConcentrationNanomolar = Primer3DnaConcentrationNanomolar) =>
        CalculatePrimer3OligoStructureCore(primer, monovalentMillimolar, divalentMillimolar, dntpMillimolar,
            dnaConcentrationNanomolar, NtthalMaxAlignLength);

    private static Primer3OligoStructure? CalculatePrimer3OligoStructureCore(
        string primer, double monovalentMillimolar, double divalentMillimolar, double dntpMillimolar,
        double dnaConcentrationNanomolar, int maxAlign)
    {
        if (!IsAcgtOnly(primer))
            return null;
        string p = primer.ToUpperInvariant();
        double mv = monovalentMillimolar / 1000.0, dv = divalentMillimolar / 1000.0, dntp = dntpMillimolar / 1000.0;
        double conc = dnaConcentrationNanomolar * 1e-9;
        double any = TmOrZero(NtthalDimer.Run(p, p, mv, conc, NtthalDimer.AlignmentType.Any, dv, dntp, maxAlign: maxAlign));
        double end = TmOrZero(NtthalDimer.Run(p, p, mv, conc, NtthalDimer.AlignmentType.End1, dv, dntp, maxAlign: maxAlign));
        var h = NtthalHairpin.Run(p, mv, dv, dntp, maxAlign: maxAlign);
        return new Primer3OligoStructure(any, end, h is null ? 0.0 : Math.Max(0.0, h.Value.TmCelsius));
    }

    // characterize_pair (thermodynamic mode): s1 = left, s2_rev = right, s2 = revcomp(right),
    // s1_rev = revcomp(left); compl_any = ANY(s1, s2_rev), compl_end = max(END1/END2(s1, s2_rev),
    // END1/END2(s2, s1_rev)). Stops as soon as a value exceeds stopAbove (the remaining values are
    // then irrelevant to a pass/fail decision).
    private static (double ComplAny, double ComplEnd) Primer3PairTms(
        string l, string r, double mv, double dv, double dntp, double conc, double stopAbove)
    {
        double any = TmOrZero(NtthalDimer.Run(l, r, mv, conc, NtthalDimer.AlignmentType.Any, dv, dntp));
        if (any > stopAbove)
            return (any, 0.0);
        string rcL = DnaSequence.GetReverseComplementString(l), rcR = DnaSequence.GetReverseComplementString(r);
        double end = 0.0;
        foreach (var (a, b, t) in new[]
                 {
                     (l, r, NtthalDimer.AlignmentType.End1), (l, r, NtthalDimer.AlignmentType.End2),
                     (rcR, rcL, NtthalDimer.AlignmentType.End1), (rcR, rcL, NtthalDimer.AlignmentType.End2),
                 })
        {
            end = Math.Max(end, TmOrZero(NtthalDimer.Run(a, b, mv, conc, t, dv, dntp)));
            if (end > stopAbove)
                break;
        }
        return (any, end);
    }

    // libprimer3.cc align_thermod: Tm of the structure, 0 when none forms or when Tm < 0 °C.
    private static double TmOrZero(NtthalDimer.Result? r) => r is null ? 0.0 : Math.Max(0.0, r.Value.TmCelsius);

    // Primer3 pair search (make_detection_primer_lists + choose_pair_or_triple + characterize_pair + obj_fn).
    private sealed class PrimerPairSearch
    {
        private readonly int _targetStart, _targetEnd, _incStart, _incEnd;
        private readonly PrimerParameters _param;
        private readonly PrimerPairOptions _opt;
        private readonly Primer3PairWeights _w;
        private readonly List<(PrimerCandidate C, double Tm)> _fwd = new(), _rev = new();
        private readonly bool?[] _fwdOk, _revOk;
        // Primer3 primer_rec.overlaps: excluded by PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE after a pair was selected.
        private readonly bool[] _fwdUsed, _revUsed;
        private readonly Dictionary<string, bool> _structureBySequence = new(StringComparer.Ordinal);
        private readonly Dictionary<(string F, string R), (bool Fails, double? Any, double? End)> _dimer = new();
        private readonly int[] _gcPrefix;
        // long_seq_tm salt for the product Tm: Primer3 passes the primer conditions (p_args) — "skewed" by its own comment.
        private readonly double _productMonovalentEq;
        private readonly ThermodynamicProbeDesigner.Primer3ProbeSettings? _intlSettings;
        private readonly List<ProbeDesigner.Primer3Probe>? _intl;
        private readonly ProbeDesigner.Primer3Probe?[]? _intlChecked;
        private readonly bool?[]? _intlOk;
        private bool _sawCharacterized, _sawTm, _sawDimer, _sawProductTm, _sawInternal, _sawLibrary;
        // PRIMER_MISPRIMING_LIBRARY (null when absent / empty) and each candidate's repeat_sim scores.
        private readonly PrimerMisprimingLibrary? _library;
        private readonly Dictionary<PrimerCandidate, LibraryMispriming> _libraryByCandidate = new(ReferenceEqualityComparer.Instance);
        // Template mispriming (null when Primer3 never computes it: no limit / weight set in the active mode).
        private readonly TemplateContext? _template;
        private readonly bool _needPairTemplate;
        private readonly Dictionary<PrimerCandidate, TemplateMisprimingScore> _templateByCandidate = new(ReferenceEqualityComparer.Instance);
        private bool _sawTemplate;
        // Non-default PRIMER_INSIDE/OUTSIDE_PENALTY: primers are scored against the target (compute_position_penalty)
        // instead of being kept off it; null under Primer3's defaults.
        private readonly TargetPosition? _position;
        // SEQUENCE_QUALITY (null when absent).
        private readonly QualityContext? _quality;
        // PRIMER_MASK_TEMPLATE: the masked copies of the included region (null without masking).
        private readonly MaskContext? _mask;
        // PRIMER_LOWERCASE_MASKING (forced on by PRIMER_MASK_TEMPLATE): the case-preserving template (trimmed_orig_seq);
        // null when lower-case masking is off or the template has no case information.
        private readonly string? _lowercase;

        // seq: the upper-case template; caseTemplate: the same template as given (null = upper case only).
        public PrimerPairSearch(string seq, string? caseTemplate, int targetStart, int targetEnd,
            PrimerParameters param, PrimerPairOptions opt)
        {
            ArgumentNullException.ThrowIfNull(opt);
            int n = seq.Length;
            if (targetStart < 0 || targetEnd >= n || targetStart >= targetEnd)
                throw new ArgumentException("Invalid target region.");
            _targetStart = targetStart;
            _targetEnd = targetEnd;
            _param = param;
            _opt = opt;
            _w = opt.Weights ?? throw new ArgumentException("Pair weights cannot be null.", nameof(opt));
            _library = param.ActiveLibrary;
            ValidateOptions(n);
            _quality = QualityContext.Create(opt.SequenceQuality, param.EffectiveQualityRangeMax);
            // _pr_need_pair_template_mispriming[_thermod] / _pr_need_template_mispriming[_thermod].
            bool thermoTemplate = param.ThermodynamicTemplateAlignment;
            _needPairTemplate = thermoTemplate
                ? opt.MaxTemplateMisprimingTh >= 0 || _w.TemplateMisprimingTh > 0
                : opt.MaxTemplateMispriming >= 0 || _w.TemplateMispriming > 0;
            if (_needPairTemplate || param.ActiveMaxTemplateMispriming >= 0 || param.ActiveTemplateMisprimingWeight > 0)
                _template = new TemplateContext(seq, thermoTemplate, param.EffectiveMonovalentMillimolar,
                    param.EffectiveDivalentMillimolar, param.EffectiveDntpMillimolar, param.EffectiveDnaConcentrationNanomolar);
            _productMonovalentEq = Primer3MonovalentEquivalent(param.EffectiveMonovalentMillimolar,
                param.EffectiveDivalentMillimolar, param.EffectiveDntpMillimolar);
            (_incStart, _incEnd) = opt.IncludedRegion is { } inc ? (inc.Start, inc.Start + inc.Length) : (0, n);
            // primer3-py sets lowercase_masking = mask_template when masking (thermoanalysis.pyx), so PRIMER_MASK_TEMPLATE implies
            // PRIMER_LOWERCASE_MASKING.
            if ((opt.LowercaseMasking || opt.MaskTemplate) && caseTemplate is not null)
                _lowercase = caseTemplate;
            // Primer3 masks the included region (trimmed_orig_seq, original case: input lower case stays masked) once per
            // design, before the primer lists are built.
            if (opt.MaskTemplate)
                _mask = new MaskContext(_lowercase ?? seq, _incStart, _incEnd, opt);

            _gcPrefix = new int[n + 1];
            for (int i = 0; i < n; i++)
                _gcPrefix[i + 1] = _gcPrefix[i] + (seq[i] is 'G' or 'C' ? 1 : 0);

            int minProduct = int.MaxValue, maxProduct = 0;
            foreach (var r in opt.ProductSizeRanges)
            {
                minProduct = Math.Min(minProduct, r.Min);
                maxProduct = Math.Max(maxProduct, r.Max);
            }

            // Default position penalties: a primer may not overlap the target (oligo_overlaps_interval). Otherwise
            // (compute_position_penalty) a left primer's 3' end may lie anywhere up to the target's last base and a
            // right primer's from the target's first base on (make_detection_primer_lists then searches the whole
            // included region); the pair must still span the target (pair_spans_target).
            if (!opt.DefaultPositionPenalties)
                _position = new TargetPosition(targetStart, targetEnd, opt.InsidePenalty, opt.OutsidePenalty);
            int leftLimit = _position is null ? targetStart : targetEnd;   // left primer end (exclusive) ≤ leftLimit
            int rightLimit = _position is null ? targetEnd : targetStart;  // right primer start ≥ rightLimit

            // Forward candidates: 3' end before the target; Primer3 drops a left primer that starts past
            // n − min product (pick_primer_range); starts that cannot reach any product range are skipped.
            int fStartLo = Math.Max(_incStart, rightLimit + param.MinLength - maxProduct);
            int fStartHi = Math.Min(leftLimit - param.MinLength, _incEnd - minProduct);
            for (int start = fStartLo; start <= fStartHi; start++)
            {
                for (int len = param.MinLength; len <= param.MaxLength && start + len <= leftLimit; len++)
                {
                    var (candidate, tm, lib, tmp) = EvaluatePrimerCore(seq.Substring(start, len), start, true, param, evaluateStructure: false, template: _template, target: _position, quality: _quality, mask: _mask, lowercase: _lowercase);
                    if (candidate.IsValid)
                    {
                        _fwd.Add((candidate, tm));
                        if (lib is not null) _libraryByCandidate[candidate] = lib;
                        if (tmp is { } t) _templateByCandidate[candidate] = t;
                    }
                }
            }

            // Reverse candidates (evaluated as the reverse complement): 5' end on the top strand = end − 1.
            int rEndLo = Math.Max(rightLimit + param.MinLength, _incStart + minProduct);
            int rEndHi = Math.Min(_incEnd, leftLimit - param.MinLength + maxProduct);
            for (int end = rEndLo; end <= rEndHi; end++)
            {
                for (int len = param.MinLength; len <= param.MaxLength && end - len >= rightLimit; len++)
                {
                    int start = end - len;
                    var revComp = DnaSequence.GetReverseComplementString(seq.Substring(start, len));
                    var (candidate, tm, lib, tmp) = EvaluatePrimerCore(revComp, start, false, param, evaluateStructure: false, template: _template, target: _position, quality: _quality, mask: _mask, lowercase: _lowercase);
                    if (candidate.IsValid)
                    {
                        _rev.Add((candidate, tm));
                        if (lib is not null) _libraryByCandidate[candidate] = lib;
                        if (tmp is { } t) _templateByCandidate[candidate] = t;
                    }
                }
            }

            // Primer3 examines primers in increasing penalty order (sort_primer_array).
            _fwd.Sort((a, b) => CompareLeft(a.C, b.C));
            _rev.Sort((a, b) => CompareRight(a.C, b.C));
            _fwdOk = new bool?[_fwd.Count];
            _revOk = new bool?[_rev.Count];
            _fwdUsed = new bool[_fwd.Count];
            _revUsed = new bool[_rev.Count];

            if (opt.PickInternalOligo)
            {
                // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT is global in Primer3: the internal oligo follows the primer screen.
                _intlSettings = (opt.InternalOligo ?? new ThermodynamicProbeDesigner.Primer3ProbeSettings()) with
                {
                    ThermodynamicOligoAlignment = param.StructureScreen != PrimerStructureScreen.Primer3Alignment,
                    // PRIMER_LIB_AMBIGUITY_CODES_CONSENSUS and PRIMER_ANNEALING_TEMP are global too.
                    LibraryAmbiguityCodesConsensus = param.EffectiveLibraryAmbiguityCodesConsensus,
                    AnnealingTemperature = param.EffectiveAnnealingTemperature,
                    // PRIMER_QUALITY_RANGE_MIN/MAX are global as well.
                    QualityRangeMin = param.EffectiveQualityRangeMin,
                    QualityRangeMax = param.EffectiveQualityRangeMax,
                };
                // make_internal_oligo_list over the included region; choose_internal_oligo takes the
                // lowest-penalty oligo (first in enumeration order among equals) → stable sort.
                var list = ThermodynamicProbeDesigner.EnumeratePrimer3InternalOligos(seq, _incStart, _incEnd, _intlSettings, screenStructure: false, _quality?.Values, _lowercase);
                _intl = list
                    .Select((p, i) => (p, i))
                    .OrderBy(t => t.p.Penalty).ThenBy(t => t.i)
                    .Select(t => t.p)
                    .ToList();
                _intlChecked = new ProbeDesigner.Primer3Probe?[_intl.Count];
                _intlOk = new bool?[_intl.Count];
            }
        }

        private void ValidateOptions(int n)
        {
            var o = _opt;
            if (o.ProductSizeRanges is null || o.ProductSizeRanges.Count == 0)
                throw new ArgumentException("At least one product size range is required (PRIMER_PRODUCT_SIZE_RANGE).");
            int minProduct = int.MaxValue;
            foreach (var r in o.ProductSizeRanges)
            {
                if (r.Min < 1 || r.Max < r.Min)
                    throw new ArgumentException($"Invalid product size range {r.Min}-{r.Max}: need 1 ≤ min ≤ max.");
                minProduct = Math.Min(minProduct, r.Min);
            }
            _param.ValidateConditions("parameters");
            if (_param.MinLength < 1 || _param.MaxLength < _param.MinLength)
                throw new ArgumentException("Primer sizes must satisfy 1 ≤ MinLength ≤ MaxLength.");
            if (_param.MaxLength > minProduct)
                throw new ArgumentException("PRIMER_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE (Primer3 _pr_data_control).");
            if (!(o.MaxTmDifference >= 0))
                throw new ArgumentException("MaxTmDifference must be ≥ 0 °C.");
            if (o.NumReturn < 1)
                throw new ArgumentException("PRIMER_NUM_RETURN < 1 (Primer3 _pr_data_control).");
            if (o.MinLeftThreePrimeDistance < -1 || o.MinRightThreePrimeDistance < -1)
                throw new ArgumentException("Minimum 3' distance must be >= -1 (min_*_three_prime_distance) (Primer3 _pr_data_control).");
            // _pr_data_control sequence-quality checks (they precede the objective-function checks; the PRIMER_INTERNAL_*
            // values count whether or not an internal oligo is picked).
            var io = o.InternalOligo ?? new ThermodynamicProbeDesigner.Primer3ProbeSettings();
            ValidatePrimer3Quality(o.SequenceQuality, n, _param.MinQuality, io.MinQuality, _param.EffectiveQualityRangeMin,
                _param.EffectiveQualityRangeMax, _param.PenaltyWeights?.SequenceQuality ?? 0.0, io.WeightSequenceQuality, nameof(o.SequenceQuality));
            if ((_w.ProductTmLt != 0 || _w.ProductTmGt != 0) && o.ProductOptTm is null)
                throw new ArgumentException("Product temperature is part of objective function while optimum temperature is not defined (Primer3 _pr_data_control).");
            if ((_w.ProductSizeLt != 0 || _w.ProductSizeGt != 0) && o.ProductOptSize is null)
                throw new ArgumentException("Product size is part of objective function while optimum size is not defined (Primer3 _pr_data_control).");
            // _pr_data_control: GC weights without PRIMER_[INTERNAL_]OPT_GC_PERCENT (the internal-oligo check applies whether or
            // not an internal oligo is picked).
            _param.ValidateGcOptimum("parameters");
            ThermodynamicProbeDesigner.ValidatePrimer3ProbeGcOptimum(io, nameof(o.InternalOligo));
            // primer3-py: PRIMER_MASK_TEMPLATE needs the k-mer lists (PRIMER_MASK_KMERLIST_PATH).
            if (o.MaskTemplate && o.MaskKmerLists is null)
                throw new ArgumentException("masking template chosen, but path to PRIMER_MASK_KMERLIST_PATH not specified (primer3-py).");
            ValidateMaskSettings(o.MaskFailureRate, o.MaskFivePrimeDirection, o.MaskThreePrimeDirection, nameof(o.MaskFailureRate));
            if (double.IsNaN(_param.PenaltyWeights?.MaskFailureRate ?? 0.0))
                throw new ArgumentException("PRIMER_WT_MASK_FAILURE_RATE must not be NaN.");
            // _pr_data_control: PRIMER_PAIR_WT_IO_PENALTY without PRIMER_PICK_INTERNAL_OLIGO.
            if (_w.InternalOligoPenalty != 0 && !o.PickInternalOligo)
                throw new ArgumentException("Internal oligo quality is part of objective function while internal oligo choice is not required (Primer3 _pr_data_control).");
            if (!(o.MaxComplAny >= 0 && o.MaxComplAny <= short.MaxValue && o.MaxComplEnd >= 0 && o.MaxComplEnd <= short.MaxValue
                  && _param.EffectiveMaxSelfAny >= 0 && _param.EffectiveMaxSelfAny <= short.MaxValue
                  && _param.EffectiveMaxSelfEnd >= 0 && _param.EffectiveMaxSelfEnd <= short.MaxValue))
                throw new ArgumentException("Illegal value for primer complementarity restrictions (Primer3 _pr_data_control: 0 ≤ limit ≤ 32767).");
            if (double.IsNaN(o.MaxLibraryMispriming)
                || (o.MaxLibraryMispriming > short.MaxValue && _param.StructureScreen == PrimerStructureScreen.Primer3Alignment))
                throw new ArgumentException("Value too large at tag PRIMER_PAIR_MAX_LIBRARY_MISPRIMING (Primer3 _pr_data_control).");
            if (_w.LibraryMispriming != 0 && _library is null)
                throw new ArgumentException("Mispriming score is part of objective function, but mispriming library is not defined (Primer3 _pr_data_control).");
            if (double.IsNaN(o.MaxTemplateMispriming) || double.IsNaN(o.MaxTemplateMisprimingTh)
                || (o.MaxTemplateMispriming > short.MaxValue && !_param.ThermodynamicTemplateAlignment))
                throw new ArgumentException("Value too large at tag PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING (Primer3 _pr_data_control).");
            if (_w.PrimerPenalty < 0 || _w.InternalOligoPenalty < 0 || _w.DiffTm < 0 || _w.ComplAnyTh < 0 || _w.ComplEndTh < 0
                || _w.ComplAny < 0 || _w.ComplEnd < 0 || _w.LibraryMispriming < 0
                || !(_w.TemplateMispriming >= 0) || !(_w.TemplateMisprimingTh >= 0)
                || _w.ProductTmLt < 0 || _w.ProductTmGt < 0 || _w.ProductSizeLt < 0 || _w.ProductSizeGt < 0)
                throw new ArgumentException("Pair weights must be ≥ 0.");
            if (o.IncludedRegion is { } inc)
            {
                if (inc.Start < 0 || inc.Length < 1 || inc.Start + inc.Length > n)
                    throw new ArgumentException("Included region outside the template.");
                if (_targetStart < inc.Start || _targetEnd > inc.Start + inc.Length)
                    throw new ArgumentException("TARGET outside of INCLUDED_REGION (Primer3 _check_and_adjust_1_interval).");
            }
            if (o.PickInternalOligo)
            {
                var s = (o.InternalOligo ?? new ThermodynamicProbeDesigner.Primer3ProbeSettings()) with
                {
                    ThermodynamicOligoAlignment = _param.StructureScreen != PrimerStructureScreen.Primer3Alignment,
                    AnnealingTemperature = _param.EffectiveAnnealingTemperature,
                };
                ThermodynamicProbeDesigner.ValidatePrimer3ProbeSettings(s, nameof(o.InternalOligo));
                if (s.MaxSize > minProduct)
                    throw new ArgumentException("PRIMER_INTERNAL_MAX_SIZE > min PRIMER_PRODUCT_SIZE_RANGE (Primer3 _pr_data_control).");
            }
        }

        // Primer3 primer_rec_comp (sort_primer_array): penalty ascending, then Primer3 "start"
        // descending, then shorter first. A left primer's start is its 5' end (Position); a right
        // primer's start is its 5' end on the top strand (Position + Length − 1).
        private static int CompareLeft(PrimerCandidate a, PrimerCandidate b)
        {
            int c = a.Penalty.CompareTo(b.Penalty);
            if (c != 0) return c;
            c = b.Position.CompareTo(a.Position);
            return c != 0 ? c : a.Length.CompareTo(b.Length);
        }

        private static int CompareRight(PrimerCandidate a, PrimerCandidate b)
        {
            int c = a.Penalty.CompareTo(b.Penalty);
            if (c != 0) return c;
            c = (b.Position + b.Length).CompareTo(a.Position + a.Length);
            return c != 0 ? c : a.Length.CompareTo(b.Length);
        }

        // Primer3 compare_primer_pair (libprimer3.cc): quality (±1e-6), then left start descending,
        // right start (its 5' end) ascending, left length ascending, right length ascending.
        private static int ComparePair(
            double q1, PrimerCandidate l1, PrimerCandidate r1,
            double q2, PrimerCandidate l2, PrimerCandidate r2)
        {
            if (q1 + PairQualityEpsilon < q2) return -1;
            if (q1 > q2 + PairQualityEpsilon) return 1;
            int c = l2.Position.CompareTo(l1.Position);
            if (c != 0) return c;
            c = (r1.Position + r1.Length).CompareTo(r2.Position + r2.Length);
            if (c != 0) return c;
            c = l1.Length.CompareTo(l2.Length);
            return c != 0 ? c : r1.Length.CompareTo(r2.Length);
        }

        // Pair-level complementarity screen used by DesignPrimers (Primer3 characterize_pair): whether the
        // pair fails, and (thermodynamic screen) PRIMER_PAIR_COMPL_ANY_TH / _COMPL_END_TH.
        private static (bool Fails, double? Any, double? End) PairScreen(
            string forward, string reverse, PrimerParameters param, PrimerPairOptions opt)
        {
            if (param.StructureScreen == PrimerStructureScreen.Heuristic)
                return (HasPrimerDimer(forward, reverse), null, null);
            if (param.StructureScreen == PrimerStructureScreen.Primer3Alignment)
                return AlignmentPairScreen(forward, reverse, param, opt);
            if (!IsAcgtOnly(forward) || !IsAcgtOnly(reverse))
                return (false, null, null);
            // Same values as CalculatePrimer3PairComplementarity(...), stopping at the first alignment over the
            // limit (characterize_pair also fails the pair on compl_any first); complete when the pair passes.
            double max = param.EffectiveMaxStructureTm;
            // characterize_pair uses the primer conditions (thal_arg_to_use = create_thal_arg_holder(p_args)).
            var (any, end) = Primer3PairTms(forward.ToUpperInvariant(), reverse.ToUpperInvariant(),
                param.EffectiveMonovalentMillimolar / 1000.0, param.EffectiveDivalentMillimolar / 1000.0,
                param.EffectiveDntpMillimolar / 1000.0, param.EffectiveDnaConcentrationNanomolar * 1e-9, max);
            return (any > max || end > max, any, end);
        }

        // characterize_pair, PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0 (s1 = left, s2 = revcomp(right),
        // s1_rev = revcomp(left), s2_rev = right): compl_any = align(s1, s2, LOCAL) > PRIMER_PAIR_MAX_COMPL_ANY fails;
        // compl_end = align(s1, s2, GLOBAL_END) > PRIMER_PAIR_MAX_COMPL_END fails; then, when
        // align(s2_rev, s1_rev, GLOBAL_END) is larger, it fails above PRIMER_MAX_SELF_END (Primer3 uses the
        // per-primer limit there) and becomes compl_end.
        private static (bool Fails, double? Any, double? End) AlignmentPairScreen(
            string forward, string reverse, PrimerParameters param, PrimerPairOptions opt)
        {
            string l = forward.ToUpperInvariant(), r = reverse.ToUpperInvariant();
            string rcL = ReverseComplementPrimer3(l), rcR = ReverseComplementPrimer3(r);
            double any = DpalLocalScore(l, rcR);
            if (any > opt.MaxComplAny)
                return (true, any, null);
            double end = DpalGlobalEndScore(l, rcR);
            if (end > opt.MaxComplEnd)
                return (true, any, end);
            double end2 = DpalGlobalEndScore(r, rcL);
            if (end2 > end)
            {
                if (end2 > param.EffectiveMaxSelfEnd)
                    return (true, any, end2);
                end = end2;
            }
            return (false, any, end);
        }

        private sealed record PairEval(
            int Fi, int Ri, double Penalty, int ProductSize, double ProductTm,
            double? ComplAny, double? ComplEnd, ProbeDesigner.Primer3Probe? Internal,
            int? LibraryScore = null, string? LibraryName = null, double? TemplateScore = null);

        public List<PrimerPairResult> Run(int numReturn)
        {
            var results = new List<PrimerPairResult>();
            if (_fwd.Count == 0 || _rev.Count == 0)
                return results;
            // PRIMER_PICK_INTERNAL_OLIGO with no acceptable internal oligo: make_internal_oligo_list fails and Primer3 returns
            // before examining any pair.
            if (_intl is { Count: 0 })
            {
                _sawCharacterized = _sawInternal = true;
                return results;
            }

            var ranges = _opt.ProductSizeRanges;
            var cache = new Dictionary<(int Ri, int Fi), PairEval?>();
            double wq = _w.PrimerPenalty;
            int rangeIndex = 0;
            while (true)
            {
                PairEval? best = null;
                for (int i = 0; i < _rev.Count; i++)
                {
                    // Only primers that are (still) legal; their expensive checks run in characterize_pair.
                    if (_revOk[i] == false)
                        continue;
                    // No pair with this or any later reverse primer can beat the best pair.
                    if (wq * (_rev[i].C.Penalty + _fwd[0].C.Penalty) > BestQuality(best))
                        break;
                    if (_revUsed[i])
                        continue; // 3′ end too close to a right primer of a selected pair

                    for (int j = 0; j < _fwd.Count; j++)
                    {
                        if (_revOk[i] == false)
                            break;
                        if (_fwdOk[j] == false)
                            continue;
                        if (wq * (_fwd[j].C.Penalty + _rev[i].C.Penalty) > BestQuality(best))
                            break;
                        if (_fwdUsed[j])
                            continue;

                        int product = ProductSize(j, i);
                        var range = ranges[rangeIndex];
                        if (product < range.Min || product > range.Max)
                            continue;

                        if (!cache.TryGetValue((i, j), out var e))
                        {
                            e = Characterize(j, i, product);
                            cache[(i, j)] = e;
                        }
                        if (e is null)
                            continue; // illegal, or already selected

                        if (best is null || ComparePair(e.Penalty, _fwd[j].C, _rev[i].C,
                                best.Penalty, _fwd[best.Fi].C, _rev[best.Ri].C) < 0)
                            best = e;
                        if (best.Penalty == 0)
                            break; // there cannot be a better pair
                    }
                    if (best is { Penalty: 0 })
                        break;
                }

                if (best is null)
                {
                    // No pair in this product-size range: try the next one, if any.
                    if (++rangeIndex >= ranges.Count)
                        break;
                    continue;
                }

                results.Add(ToResult(best));
                cache[(best.Ri, best.Fi)] = null; // mark as selected
                MarkUsedPrimers(best);
                if (results.Count == numReturn)
                    break;
            }
            return results;
        }

        private static double BestQuality(PairEval? best) => best?.Penalty ?? double.MaxValue;

        // Primer3 left/right_oligo_in_pair_overlaps_used_oligo after a pair is selected: min distance −1 → nothing;
        // 0 → the identical primer (same start and length); d > 0 → every primer whose 3′ end is < d bases from the
        // selected primer's 3′ end (left 3′ end = Position + Length − 1; right 3′ end on the top strand = Position).
        private void MarkUsedPrimers(PairEval best)
        {
            MarkUsed(_fwd, _fwdUsed, _fwd[best.Fi].C, _opt.MinLeftThreePrimeDistance, c => c.Position + c.Length - 1);
            MarkUsed(_rev, _revUsed, _rev[best.Ri].C, _opt.MinRightThreePrimeDistance, c => c.Position);

            static void MarkUsed(List<(PrimerCandidate C, double Tm)> list, bool[] used, PrimerCandidate chosen,
                int minDistance, Func<PrimerCandidate, int> threePrime)
            {
                if (minDistance == -1)
                    return;
                int chosenEnd = threePrime(chosen);
                for (int k = 0; k < list.Count; k++)
                {
                    var c = list[k].C;
                    if (minDistance == 0
                        ? c.Position == chosen.Position && c.Length == chosen.Length
                        : Math.Abs(threePrime(c) - chosenEnd) < minDistance)
                        used[k] = true;
                }
            }
        }

        private int ProductSize(int fi, int ri) =>
            _rev[ri].C.Position + _rev[ri].C.Length - _fwd[fi].C.Position;

        // characterize_pair (+ choose_internal_oligo and obj_fn for a legal pair).
        private PairEval? Characterize(int fi, int ri, int product)
        {
            _sawCharacterized = true;
            var f = _fwd[fi];
            var r = _rev[ri];

            // pair_spans_target: the left primer's 3' end must precede the right primer's (always true when neither
            // may overlap the target).
            if (_position is not null && f.C.Position + f.C.Length - 1 >= r.C.Position)
                return null;

            int fStart = f.C.Position;
            int gc = _gcPrefix[fStart + product] - _gcPrefix[fStart];
            double productTm = LongSeqTm(gc, product, _productMonovalentEq);
            if ((_opt.ProductMinTm is { } minTm && productTm < minTm)
                || (_opt.ProductMaxTm is { } maxTm && productTm > maxTm))
            {
                _sawProductTm = true;
                return null;
            }

            double diffTm = Math.Abs(f.Tm - r.Tm);
            if (diffTm > _opt.MaxTmDifference)
            {
                _sawTm = true;
                return null;
            }

            if (!StructureOk(_fwd, _fwdOk, fi) || !StructureOk(_rev, _revOk, ri))
                return null;

            var (fails, any, end) = PairComplementarity(f.C.Sequence, r.C.Sequence);
            if (fails)
            {
                _sawDimer = true;
                return null;
            }

            // pair_repeat_sim > PRIMER_PAIR_MAX_LIBRARY_MISPRIMING fails the pair.
            int? libScore = null;
            string? libName = null;
            if (_library is not null)
            {
                (int score, libName) = PairLibraryMispriming(LibraryOf(f.C), LibraryOf(r.C), _library);
                if (score > _opt.MaxLibraryMispriming)
                {
                    _sawLibrary = true;
                    return null;
                }
                libScore = score;
            }

            // Pair template mispriming (end of characterize_pair): max(left.T + right.T_r, left.T_r + right.T); fails
            // above PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING when that is ≥ 0, or (thermodynamic mode) above
            // PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH when that is non-zero.
            double? templateScore = null;
            if (_needPairTemplate)
            {
                if (!_templateByCandidate.TryGetValue(f.C, out var lt) || !_templateByCandidate.TryGetValue(r.C, out var rt))
                    throw new InvalidOperationException(
                        "Primer3 PR_ASSERT(template_mispriming != ALIGN_SCORE_UNDEF): with a mispriming library scored at pick time (PRIMER_WT_LIBRARY_MISPRIMING ≠ 0) Primer3 never scores template mispriming for the pair unless PRIMER_WT_TEMPLATE_MISPRIMING[_TH] ≠ 0.");
                double v = lt.SameStrand + rt.OtherStrand;
                if (lt.OtherStrand + rt.SameStrand > v)
                    v = lt.OtherStrand + rt.SameStrand;
                bool templateFails = _param.ThermodynamicTemplateAlignment
                    ? _opt.MaxTemplateMisprimingTh != 0 && v > _opt.MaxTemplateMisprimingTh
                    : _opt.MaxTemplateMispriming >= 0 && v > _opt.MaxTemplateMispriming;
                if (templateFails)
                {
                    _sawTemplate = true;
                    return null;
                }
                templateScore = v;
            }

            ProbeDesigner.Primer3Probe? intl = null;
            if (_intl is not null)
            {
                intl = ChooseInternalOligo(f.C, r.C);
                if (intl is null)
                {
                    _sawInternal = true;
                    return null;
                }
            }

            double penalty = ObjectiveFunction(f, r, diffTm, any, end, productTm, product, intl, libScore ?? 0, templateScore ?? 0);
            return new PairEval(fi, ri, penalty, product, productTm, any, end, intl, libScore, libName, templateScore);
        }

        // The candidate's repeat_sim scores (computed at pick time when weighted, otherwise now).
        private LibraryMispriming LibraryOf(PrimerCandidate c)
        {
            if (!_libraryByCandidate.TryGetValue(c, out var lib))
            {
                lib = ComputeLibraryMispriming(c.Sequence, c.IsForward, _library!,
                    _param.EffectiveLibraryAmbiguityCodesConsensus, _param.EffectiveMaxLibraryMispriming);
                _libraryByCandidate[c] = lib;
            }
            return lib;
        }

        // obj_fn (libprimer3.cc): thermodynamic mode (ComplAnyTh/ComplEndTh terms), alignment mode (linear
        // ComplAny/ComplEnd terms) or, under the heuristic screen, without complementarity terms.
        private double ObjectiveFunction(
            (PrimerCandidate C, double Tm) f, (PrimerCandidate C, double Tm) r, double diffTm,
            double? complAny, double? complEnd, double productTm, int product, ProbeDesigner.Primer3Probe? intl,
            int libraryScore, double templateScore)
        {
            double sum = 0.0;
            double lowerTm = r.Tm;
            if (f.Tm < r.Tm) lowerTm = f.Tm;

            if (_w.PrimerPenalty != 0)
                sum += _w.PrimerPenalty * (f.C.Penalty + r.C.Penalty);
            if (_w.InternalOligoPenalty != 0 && intl is { } io)
                sum += _w.InternalOligoPenalty * io.Penalty;
            if (_w.DiffTm != 0)
                sum += _w.DiffTm * diffTm;
            if (_param.StructureScreen == PrimerStructureScreen.Primer3Alignment)
            {
                // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0: linear PRIMER_PAIR_WT_COMPL_ANY / _COMPL_END terms.
                if (_w.ComplAny != 0)
                    sum += _w.ComplAny * complAny!.Value;
                if (_w.ComplEnd != 0)
                    sum += _w.ComplEnd * complEnd!.Value;
            }
            else
            {
                if (complAny is { } a)
                    sum += ThermodynamicStructurePenalty(_w.ComplAnyTh, lowerTm, a);
                if (complEnd is { } e)
                    sum += ThermodynamicStructurePenalty(_w.ComplEndTh, lowerTm, e);
            }
            if (_w.ProductTmLt != 0 && productTm < _opt.ProductOptTm!.Value)
                sum += _w.ProductTmLt * (_opt.ProductOptTm.Value - productTm);
            if (_w.ProductTmGt != 0 && productTm > _opt.ProductOptTm!.Value)
                sum += _w.ProductTmGt * (productTm - _opt.ProductOptTm.Value);
            if (_w.ProductSizeLt != 0 && product < _opt.ProductOptSize!.Value)
                sum += _w.ProductSizeLt * (_opt.ProductOptSize.Value - product);
            if (_w.ProductSizeGt != 0 && product > _opt.ProductOptSize!.Value)
                sum += _w.ProductSizeGt * (product - _opt.ProductOptSize.Value);
            if (_w.LibraryMispriming != 0)
                sum += _w.LibraryMispriming * libraryScore;
            if (!_param.ThermodynamicTemplateAlignment)
            {
                if (_w.TemplateMispriming != 0)
                    sum += _w.TemplateMispriming * templateScore;
            }
            else
                sum += ThermodynamicStructurePenalty(_w.TemplateMisprimingTh, lowerTm, templateScore);
            // obj_fn ends with PR_ASSERT(sum >= 0.0): Primer3 aborts on a negative pair penalty, which only negative
            // per-primer terms produce (e.g. the default PRIMER_INSIDE_PENALTY −1 with a changed PRIMER_OUTSIDE_PENALTY).
            if (sum < 0.0)
                throw new InvalidOperationException(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Negative primer-pair penalty {sum} (Primer3 obj_fn PR_ASSERT(sum >= 0.0) aborts): negative per-primer penalties — e.g. PRIMER_INSIDE_PENALTY < 0 (default −1) with a non-default PRIMER_OUTSIDE_PENALTY — are not usable."));
            return sum;
        }

        // choose_internal_oligo: the lowest-penalty internal oligo lying strictly between the primers
        // (start after the left primer's 3' end, end before the right primer's 5'-most top-strand base)
        // whose postponed self-any / self-end / hairpin checks pass; a failing oligo stays rejected.
        private ProbeDesigner.Primer3Probe? ChooseInternalOligo(PrimerCandidate left, PrimerCandidate right)
        {
            int leftEnd = left.Position + left.Length - 1;
            for (int k = 0; k < _intl!.Count; k++)
            {
                var h = _intl[k];
                if (h.Start <= leftEnd || h.Start + h.Length - 1 >= right.Position)
                    continue;
                if (_intlOk![k] is null)
                {
                    _intlChecked![k] = ThermodynamicProbeDesigner.PassesPrimer3ProbeStructure(h, _intlSettings!);
                    _intlOk[k] = _intlChecked[k] is not null;
                }
                if (_intlOk[k] == true)
                    return _intlChecked![k];
            }
            return null;
        }

        private bool StructureOk(List<(PrimerCandidate C, double Tm)> list, bool?[] cache, int i)
        {
            if (cache[i] is { } known)
                return known;
            string seq = list[i].C.Sequence;
            if (!_structureBySequence.TryGetValue(seq, out bool ok))
            {
                var issues = new List<string>();
                AddStructureIssues(seq, _param, issues);
                ok = issues.Count == 0;
                _structureBySequence[seq] = ok;
            }
            // characterize_pair then checks the primer against the mispriming library (oligo_repeat_library_mispriming).
            if (ok && _library is not null && LibraryOf(list[i].C).Exceeds)
            {
                ok = false;
                _sawLibrary = true;
            }
            // ... and, when its library scores were not computed at pick time (repeat_sim.score == NULL: no library, or
            // PRIMER_WT_LIBRARY_MISPRIMING = 0), against the template (oligo_template_mispriming).
            if (ok && _template is not null
                && (_library is null || EffectivePenaltyWeights(_param).LibraryMispriming == 0)
                && TemplateMisprimingExceeds(TemplateOf(list[i].C), _param))
            {
                ok = false;
                _sawTemplate = true;
            }
            cache[i] = ok;
            return ok;
        }

        private (bool Fails, double? Any, double? End) PairComplementarity(string f, string r)
        {
            if (!_dimer.TryGetValue((f, r), out var v))
            {
                v = PairScreen(f, r, _param, _opt);
                _dimer[(f, r)] = v;
            }
            return v;
        }

        // The candidate's template mispriming scores (computed at pick time when weighted, otherwise now).
        private TemplateMisprimingScore TemplateOf(PrimerCandidate c)
        {
            if (!_templateByCandidate.TryGetValue(c, out var t))
            {
                t = _template!.Score(c.Sequence, c.Position, c.IsForward);
                _templateByCandidate[c] = t;
            }
            return t;
        }

        private PrimerPairResult ToResult(PairEval e)
        {
            var forward = WithTemplate(Reevaluate(_fwd[e.Fi].C), _fwd[e.Fi].C);
            var reverse = WithTemplate(Reevaluate(_rev[e.Ri].C), _rev[e.Ri].C);
            bool alignment = _param.StructureScreen == PrimerStructureScreen.Primer3Alignment;
            return new PrimerPairResult(
                Forward: forward,
                Reverse: reverse,
                IsValid: true,
                Message: "Valid primer pair found.",
                ProductSize: e.ProductSize,
                PairPenalty: e.Penalty,
                ProductTm: e.ProductTm,
                ComplAnyTh: alignment ? null : e.ComplAny,
                ComplEndTh: alignment ? null : e.ComplEnd,
                InternalOligo: e.Internal)
            {
                ComplAny = alignment ? e.ComplAny : null,
                ComplEnd = alignment ? e.ComplEnd : null,
                LibraryMispriming = e.LibraryScore,
                LibraryMisprimingName = e.LibraryName,
                TemplateMispriming = e.TemplateScore,
            };
        }

        private PrimerCandidate WithTemplate(PrimerCandidate full, PrimerCandidate picked) =>
            _templateByCandidate.TryGetValue(picked, out var t) ? full with { TemplateMispriming = t.Max } : full;

        // Full evaluation (including the structure values) of a chosen primer.
        private PrimerCandidate Reevaluate(PrimerCandidate c) =>
            EvaluatePrimerCore(c.Sequence, c.Position, c.IsForward, _param, template: _template, target: _position, quality: _quality, mask: _mask, lowercase: _lowercase).Candidate;

        // Result when no pair qualifies: the individually lowest-penalty primers that pass their own
        // (structure) constraints, with the violated pair constraint.
        public PrimerPairResult Failure()
        {
            const string noPrimers = "Could not find valid primers for the target region.";
            string noProduct = $"No primer pair with a product size in {string.Join(", ", _opt.ProductSizeRanges.Select(x => $"{x.Min}-{x.Max}"))} bp (PRIMER_PRODUCT_SIZE_RANGE).";
            if (_fwd.Count == 0 || _rev.Count == 0)
            {
                // A search region shorter than the smallest product cannot hold any pair (Primer3
                // _pr_data_control per-sequence error).
                int minProduct = _opt.ProductSizeRanges.Min(x => x.Min);
                return new PrimerPairResult(null, null, false,
                    _incEnd - _incStart < minProduct ? "SEQUENCE_INCLUDED_REGION length < min PRIMER_PRODUCT_SIZE_RANGE" : noPrimers, 0);
            }
            int f0i = FirstStructurallyValid(_fwd, _fwdOk);
            int r0i = FirstStructurallyValid(_rev, _revOk);
            if (f0i < 0 || r0i < 0)
                return new PrimerPairResult(null, null, false, noPrimers, 0);

            var f0 = Reevaluate(_fwd[f0i].C);
            var r0 = Reevaluate(_rev[r0i].C);
            double maxDiff = _opt.MaxTmDifference;
            string reason;
            if (!_sawCharacterized)
                reason = noProduct;
            else if (_sawProductTm || _sawInternal || _sawLibrary || _sawTemplate)
            {
                var parts = new List<string>();
                if (_sawLibrary) parts.Add("the mispriming-library limits");
                if (_sawTemplate) parts.Add("the template-mispriming limits");
                if (_sawProductTm) parts.Add("product Tm limits");
                if (_sawTm) parts.Add($"the {maxDiff:0.##}°C Tm-difference limit");
                if (_sawDimer) parts.Add("primer-dimer avoidance");
                if (_sawInternal) parts.Add("an acceptable internal oligo between the primers");
                reason = $"No primer pair satisfies {string.Join(", ", parts)}.";
            }
            else if (_sawTm && !_sawDimer)
                reason = $"No primer pair within the {maxDiff:0.##}°C Tm-difference limit (best primers: Tm {f0.MeltingTemperature:F1}/{r0.MeltingTemperature:F1}°C).";
            else if (_sawDimer && !_sawTm)
                reason = "Every candidate primer pair forms a primer-dimer.";
            else
                reason = $"No primer pair satisfies both the {maxDiff:0.##}°C Tm-difference limit and primer-dimer avoidance.";
            return new PrimerPairResult(
                Forward: f0,
                Reverse: r0,
                IsValid: false,
                Message: reason,
                ProductSize: r0.Position + r0.Length - f0.Position);

            int FirstStructurallyValid(List<(PrimerCandidate C, double Tm)> list, bool?[] cache)
            {
                for (int i = 0; i < list.Count; i++)
                    if (StructureOk(list, cache, i))
                        return i;
                return -1;
            }
        }
    }

    // Evaluates a candidate and also returns its unrounded Tm (Primer3 compares unrounded Tm values).
    // With evaluateStructure = false the secondary-structure screen is skipped (DesignPrimers runs it
    // lazily, like Primer3's characterize_pair, and re-evaluates the chosen primers in full).
    // template (pair search only): Primer3 scores template mispriming at pick time when the active template weight is
    // non-zero (calc_and_check_oligo_features), for an oligo that passed every earlier check.
    private static (PrimerCandidate Candidate, double Tm, LibraryMispriming? Library, TemplateMisprimingScore? Template) EvaluatePrimerCore(
        string sequence,
        int position,
        bool isForward,
        PrimerParameters param,
        bool evaluateStructure = true,
        TemplateContext? template = null,
        TargetPosition? target = null,
        QualityContext? quality = null,
        MaskContext? mask = null,
        string? lowercase = null)
    {
        var seq = sequence.ToUpperInvariant();

        double gcContent = CalculateGcContent(seq);
        // seqtm (MAX_NN_TM_LENGTH 36) with the fraction bound at PRIMER_ANNEALING_TEMP.
        double annealing = param.EffectiveAnnealingTemperature;
        var (tmRaw, boundRaw) = Primer3SeqTm(seq, param.EffectiveDnaConcentrationNanomolar,
            Primer3MonovalentEquivalent(param.EffectiveMonovalentMillimolar, param.EffectiveDivalentMillimolar,
                param.EffectiveDntpMillimolar), annealing);
        bool tmComputable = !double.IsNaN(tmRaw);
        double tm = tmComputable ? tmRaw : 0.0;
        int homopolymer = FindLongestHomopolymer(seq);
        int dinucRepeat = FindLongestDinucleotideRepeat(seq);
        double stability3Prime = Calculate3PrimeStability(seq);

        var issues = new List<string>();

        // Validate against parameters
        if (seq.Length < param.MinLength || seq.Length > param.MaxLength)
            issues.Add($"Length {seq.Length} outside range [{param.MinLength}-{param.MaxLength}]");

        if (gcContent < param.MinGcContent || gcContent > param.MaxGcContent)
            issues.Add($"GC content {gcContent:F1}% outside range [{param.MinGcContent}-{param.MaxGcContent}]%");

        if (!tmComputable)
            issues.Add("Tm not computable: sequence contains a non-ACGT base");
        else if (tm < param.MinTm || tm > param.MaxTm)
            issues.Add($"Tm {tm:F1}°C outside range [{param.MinTm}-{param.MaxTm}]°C");

        // Primer3 fraction bound (calc_and_check_oligo_features, after the Tm checks): only when PRIMER_ANNEALING_TEMP > 0.
        if (annealing > 0.0)
        {
            if (boundRaw < param.EffectiveMinBound)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Fraction bound {BoundText(boundRaw)} below {param.EffectiveMinBound:0.##}% (Primer3 PRIMER_MIN_BOUND)"));
            if (boundRaw > param.EffectiveMaxBound)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Fraction bound {BoundText(boundRaw)} above {param.EffectiveMaxBound:0.##}% (Primer3 PRIMER_MAX_BOUND)"));
        }

        // PRIMER_MASK_TEMPLATE (is_lowercase_masked, among the first checks of calc_and_check_oligo_features): the primer's
        // 3' base must not be masked on its strand's copy of the template.
        // With masking the masked copies (built from the case-preserving template) carry the input's lower case too.
        if (mask is not null && mask.ThreePrimeEndMasked(position, seq.Length, isForward))
            issues.Add("3' end overlaps masked sequence (Primer3 PRIMER_MASK_TEMPLATE)");
        else if (mask is null && lowercase is not null
                 && IsLowercaseMaskedBase(lowercase[isForward ? position + seq.Length - 1 : position]))
            issues.Add("3' end overlaps lower-case masked sequence (Primer3 PRIMER_LOWERCASE_MASKING)");

        // Primer3 position penalty relative to the target (non-default PRIMER_INSIDE/OUTSIDE_PENALTY, pair search only):
        // a 3' end past the target is an infinite position penalty (OP_OVERLAPS_TARGET).
        double? positionPenalty = null;
        if (target is { } tp)
        {
            positionPenalty = CalculatePositionPenaltyPrimer3(position, seq.Length, isForward, tp.Start, tp.End,
                tp.InsidePenalty, tp.OutsidePenalty);
            if (positionPenalty is null)
                issues.Add("3' end beyond the target (Primer3 infinite position penalty, overlaps target)");
        }

        if (homopolymer > param.MaxHomopolymer)
            issues.Add($"Homopolymer run of {homopolymer} exceeds max {param.MaxHomopolymer}");

        if (dinucRepeat > param.MaxDinucleotideRepeats)
            issues.Add($"Dinucleotide repeat of {dinucRepeat} exceeds max {param.MaxDinucleotideRepeats}");

        bool hasHairpin = false;
        StructureValues structure = default;
        if (evaluateStructure)
            (hasHairpin, structure) = AddStructureIssues(seq, param, issues);

        // Primer3 3′-end checks of left/right primers (calc_and_check_oligo_features): PRIMER_GC_CLAMP — the
        // gc_clamp 3′-most bases must all be G/C; PRIMER_MAX_END_GC — checked only when < 5, at most max_end_gc
        // G/C among the five 3′-most bases (a shorter primer counts its own bases); PRIMER_MAX_END_STABILITY —
        // end_stability = end_oligodg(seq, 5) = −ΔG of the 3′ pentamer must not exceed the limit.
        int gcClamp = param.GcClamp;
        for (int i = 0; i < gcClamp; i++)
        {
            int k = seq.Length - 1 - i;
            if (k < 0 || seq[k] is not ('G' or 'C'))
            {
                issues.Add($"No GC clamp: the {gcClamp} 3'-most bases must be G/C (PRIMER_GC_CLAMP)");
                break;
            }
        }

        int maxEndGc = param.EffectiveMaxEndGc;
        if (maxEndGc < Primer3MaxEndGc)
        {
            int endGc = 0;
            for (int k = Math.Max(0, seq.Length - 5); k < seq.Length; k++)
                if (seq[k] is 'G' or 'C') endGc++;
            if (endGc > maxEndGc)
                issues.Add($"{endGc} G/C in the last 5 bases exceeds max {maxEndGc} (PRIMER_MAX_END_GC)");
        }

        // Primer3 sequence_quality_is_ok (after the end-GC check; only with SEQUENCE_QUALITY): the minimum base quality,
        // then the minimum over the five 3'-most bases.
        int? minQuality = null;
        if (quality is not null)
        {
            var (qMin, qEnd) = CalculateSequenceQualityPrimer3(quality.Values, position, seq.Length, isForward, quality.RangeMax);
            minQuality = qMin;
            if (qMin < param.MinQuality)
                issues.Add($"Minimum base quality {qMin} below {param.MinQuality} (Primer3 PRIMER_MIN_QUALITY)");
            else if (qEnd < param.MinEndQuality)
                issues.Add($"Minimum 3'-end base quality {qEnd} below {param.MinEndQuality} (Primer3 PRIMER_MIN_END_QUALITY)");
        }

        double endStability = -stability3Prime;
        if (endStability > param.EffectiveMaxEndStability)
            issues.Add($"3' end stability {endStability:F2} kcal/mol exceeds max {param.EffectiveMaxEndStability} (PRIMER_MAX_END_STABILITY)");

        // Library rule kept for source compatibility (not Primer3; superseded by GcClamp / MaxEndGc).
        if (param.Avoid3PrimeGC && seq.Length >= 2)
        {
            string last2 = seq.Substring(seq.Length - 2);
            int gcCount = last2.Count(c => c == 'G' || c == 'C');
            if (gcCount == 0)
                issues.Add("No GC clamp at 3' end");
        }

        // Primer3 mispriming library (oligo_repeat_library_mispriming): at pick time only when weighted
        // (PRIMER_WT_LIBRARY_MISPRIMING ≠ 0), otherwise in characterize_pair — the pair search does that lazily.
        var weights = EffectivePenaltyWeights(param);
        LibraryMispriming? library = null;
        if (param.ActiveLibrary is { } lib && (evaluateStructure || weights.LibraryMispriming != 0))
        {
            library = ComputeLibraryMispriming(seq, isForward, lib, param.EffectiveLibraryAmbiguityCodesConsensus,
                param.EffectiveMaxLibraryMispriming);
            if (library.Exceeds)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Library mispriming score {library.Scores.Max():0.00} exceeds {param.EffectiveMaxLibraryMispriming:0.00} (Primer3 PRIMER_MAX_LIBRARY_MISPRIMING)"));
        }

        TemplateMisprimingScore? templateScore = null;
        if (template is not null && issues.Count == 0 && param.ActiveTemplateMisprimingWeight != 0)
        {
            templateScore = template.Score(seq, position, isForward);
            if (TemplateMisprimingExceeds(templateScore.Value, param))
                issues.Add(TemplateMisprimingIssue(templateScore.Value, param));
        }
        // Without a template (EvaluatePrimer) the template mispriming score is undefined: no template term
        // (callers combine CalculateTemplateMispriming with CalculatePrimer3Penalty themselves).
        if (templateScore is null)
            weights = weights with { TemplateMispriming = 0.0, TemplateMisprimingTh = 0.0 };

        bool isValid = issues.Count == 0;

        // Masker failure rate (end of calc_and_check_oligo_features, 16-nt window): only with PRIMER_MASK_TEMPLATE.
        double? failureRate = mask is null ? null : CalculateMaskFailureRatePrimer3(seq, mask.Lists);

        // Informational heuristic score and the Primer3 ranking penalty.
        double score = CalculatePrimerScore(seq, gcContent, tm, homopolymer, param);
        if (!evaluateStructure && UsesStructureTerms(weights, param.StructureScreen))
            structure = ComputeStructureValues(seq, param); // Primer3 computes them at pick time when weighted
        double penalty = CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(tm, seq.Length, gcContent,
                SelfAny: (weights.ThermodynamicOligoAlignment ? structure.Thermo?.SelfAnyTh : structure.SelfAny) ?? 0.0,
                SelfEnd: (weights.ThermodynamicOligoAlignment ? structure.Thermo?.SelfEndTh : structure.SelfEnd) ?? 0.0,
                HairpinTh: structure.Thermo?.HairpinTh ?? 0.0,
                // p_obj_fn end_stability term (left/right primers only): h->end_stability =
                // end_oligodg(seq, 5, santalucia) = −ΔG of the 3′ pentamer (positive magnitude).
                EndStability: double.IsNaN(stability3Prime) ? 0.0 : -stability3Prime)
            {
                LibraryMispriming = library?.MaxScore ?? 0.0,
                TemplateMispriming = templateScore?.Max ?? 0.0,
                Bound = annealing > 0.0 ? boundRaw : null,
                PositionPenalty = positionPenalty ?? 0.0,
                SequenceQuality = minQuality,
                QualityRangeMax = param.EffectiveQualityRangeMax,
                MaskFailureRate = failureRate ?? 0.0,
            },
            weights,
            new Primer3Optima(param.OptimalTm, param.OptimalLength, param.OptimalGcPercent) { OptBound = param.EffectiveOptBound });

        var candidate = new PrimerCandidate(
            Sequence: seq,
            Position: position,
            IsForward: isForward,
            Length: seq.Length,
            GcContent: Math.Round(gcContent, 1),
            MeltingTemperature: Math.Round(tm, 1),
            HomopolymerLength: homopolymer,
            HasHairpin: hasHairpin,
            Stability3Prime: Math.Round(stability3Prime, 1),
            IsValid: isValid,
            Issues: issues.AsReadOnly(),
            Score: Math.Round(score, 2),
            Penalty: penalty,
            SelfAnyTh: structure.Thermo?.SelfAnyTh,
            SelfEndTh: structure.Thermo?.SelfEndTh,
            HairpinTh: structure.Thermo?.HairpinTh
        )
        {
            SelfAny = structure.SelfAny,
            SelfEnd = structure.SelfEnd,
            LibraryMispriming = library?.MaxScore,
            LibraryMisprimingName = library?.Name,
            TemplateMispriming = templateScore?.Max,
            Bound = annealing > 0.0 && boundRaw != Primer3OligoTmError ? boundRaw : null,
            PositionPenalty = positionPenalty,
            MinSequenceQuality = minQuality,
            MaskFailureRate = failureRate,
        };
        return (candidate, tm, library, templateScore);

        static string BoundText(double b) => b == Primer3OligoTmError
            ? "undefined (Primer3 OLIGOTM_ERROR: > 36 bases or non-ACGT)"
            : string.Create(System.Globalization.CultureInfo.InvariantCulture, $"{b:0.##}%");
    }

    private static PrimerPairResult DesignPrimersCore(string seq, string? caseTemplate, int targetStart, int targetEnd,
        PrimerParameters? parameters, PrimerPairOptions? pairOptions)
    {
        var search = new PrimerPairSearch(seq, caseTemplate, targetStart, targetEnd,
            parameters ?? DefaultParameters, pairOptions ?? PrimerPairOptions.Default);
        var pairs = search.Run(1);
        return pairs.Count > 0 ? pairs[0] : search.Failure();
    }

    /// <summary>
    /// Returns up to <see cref="PrimerPairOptions.NumReturn"/> primer pairs (Primer3 PRIMER_NUM_RETURN,
    /// default 5), best first, exactly as Primer3's <c>choose_pair_or_triple</c> picks them (see
    /// <see cref="DesignPrimers"/>): after a pair is selected it is removed and the search is repeated;
    /// by default primers may be reused in later pairs (Primer3 PRIMER_MIN_LEFT/RIGHT_THREE_PRIME_DISTANCE = −1),
    /// otherwise a selected pair excludes every left / right primer whose 3′ end lies closer than
    /// <see cref="PrimerPairOptions.MinLeftThreePrimeDistance"/> / <see cref="PrimerPairOptions.MinRightThreePrimeDistance"/>
    /// to its own (0: only the identical primer) from all later pairs; product-size ranges are tried in order. Each result is valid and carries
    /// <see cref="PrimerPairResult.PairPenalty"/> (PRIMER_PAIR_k_PENALTY),
    /// <see cref="PrimerPairResult.ProductTm"/> (PRIMER_PAIR_k_PRODUCT_TM), the pair complementarity Tm
    /// values and, with <see cref="PrimerPairOptions.PickInternalOligo"/>, the internal oligo
    /// (PRIMER_INTERNAL_k_*). An empty list means no pair satisfies the constraints
    /// (<see cref="DesignPrimers"/> then reports why).
    /// </summary>
    /// <param name="template">The DNA template sequence.</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Per-primer design parameters (default <see cref="DefaultParameters"/>).</param>
    /// <param name="pairOptions">Pair options (default <see cref="PrimerPairOptions.Default"/>).</param>
    /// <returns>The selected pairs, rank 0 first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentException">Invalid target or options.</exception>
    public static IReadOnlyList<PrimerPairResult> DesignPrimerPairs(
        DnaSequence template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null,
        PrimerPairOptions? pairOptions = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        var opts = pairOptions ?? PrimerPairOptions.Default;
        return new PrimerPairSearch(template.Sequence, null, targetStart, targetEnd, parameters ?? DefaultParameters, opts)
            .Run(opts.NumReturn);
    }

    /// <summary>
    /// <see cref="DesignPrimerPairs(DnaSequence, int, int, PrimerParameters?, PrimerPairOptions?)"/> on a case-preserving
    /// template string: with <see cref="PrimerPairOptions.LowercaseMasking"/> (PRIMER_LOWERCASE_MASKING) or
    /// <see cref="PrimerPairOptions.MaskTemplate"/> primers and internal oligos whose 3′-terminal template base is lower case
    /// are rejected (Primer3 <c>is_lowercase_masked</c>); otherwise case is ignored.
    /// </summary>
    /// <param name="template">Template (A/C/G/T, any case; lower case marks masked bases).</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Per-primer design parameters (default <see cref="DefaultParameters"/>).</param>
    /// <param name="pairOptions">Pair options (default <see cref="PrimerPairOptions.Default"/>).</param>
    /// <returns>The selected pairs, rank 0 first.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentException">A character other than A/C/G/T, an invalid target or options.</exception>
    public static IReadOnlyList<PrimerPairResult> DesignPrimerPairs(
        string template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null,
        PrimerPairOptions? pairOptions = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        var opts = pairOptions ?? PrimerPairOptions.Default;
        return new PrimerPairSearch(new DnaSequence(template).Sequence, template, targetStart, targetEnd,
            parameters ?? DefaultParameters, opts).Run(opts.NumReturn);
    }

    // Primer3 compare_primer_pair epsilon on pair_quality.
    private const double PairQualityEpsilon = 1e-6;

    /// <summary>
    /// Evaluates a single primer candidate against the per-primer constraints of
    /// <paramref name="parameters"/> (length, GC%, Tm, homopolymer, dinucleotide repeat, secondary
    /// structure, Primer3's 3′-end checks PRIMER_GC_CLAMP / PRIMER_MAX_END_GC / PRIMER_MAX_END_STABILITY
    /// (<see cref="PrimerParameters.GcClamp"/>, <see cref="PrimerParameters.MaxEndGc"/>,
    /// <see cref="PrimerParameters.MaxEndStability"/>), the deprecated library rule <c>Avoid3PrimeGC</c>). The secondary-structure screen is
    /// <see cref="PrimerParameters.StructureScreen"/>: by default Primer3's thermodynamic limits
    /// (ntthal self-dimer, 3′ self-dimer and hairpin Tm ≤ 47 °C, reported in
    /// <see cref="PrimerCandidate.SelfAnyTh"/>/<see cref="PrimerCandidate.SelfEndTh"/>/<see cref="PrimerCandidate.HairpinTh"/>),
    /// Primer3's alignment-mode limits (<see cref="PrimerStructureScreen.Primer3Alignment"/>: dpal self_any ≤ 8,
    /// self_end ≤ 3, reported in <see cref="PrimerCandidate.SelfAny"/>/<see cref="PrimerCandidate.SelfEnd"/>),
    /// or the sequence-only <see cref="HasHairpinPotential"/>. The Tm is Primer3's primer Tm
    /// (<see cref="CalculateMeltingTemperaturePrimer3"/>: SantaLucia 1998 nearest-neighbour,
    /// SantaLucia salt correction) at the reaction conditions of <paramref name="parameters"/>
    /// (<see cref="PrimerParameters.MonovalentMillimolar"/>, <see cref="PrimerParameters.DivalentMillimolar"/>,
    /// <see cref="PrimerParameters.DntpMillimolar"/>, <see cref="PrimerParameters.DnaConcentrationNanomolar"/>; Primer3's
    /// defaults 50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo), the scale on which the Primer3-sourced Tm
    /// window 57–63 °C (opt 60) is defined; the ntthal structure values use the same conditions (Primer3
    /// <c>create_thal_arg_holder(p_args)</c>). A sequence
    /// containing a non-ACGT base has no computable Tm (Primer3 PRIMER_MAX_NS_ACCEPTED = 0): it is
    /// reported with Tm 0 and an issue. <see cref="PrimerCandidate.Penalty"/> is the Primer3
    /// per-primer penalty (<see cref="CalculatePrimer3Penalty"/>, <see cref="PrimerParameters.PenaltyWeights"/>, optima
    /// OptimalTm / OptimalLength / <see cref="PrimerParameters.OptimalGcPercent"/>) used by <see cref="DesignPrimers"/> for ranking;
    /// <see cref="PrimerCandidate.Score"/> is an informational additive quality score (0–100,
    /// higher is better) that does not drive selection.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="sequence"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException">Illegal reaction conditions (Primer3 <c>_pr_data_control</c>:
    /// salt or DNA concentration ≤ 0, negative divalent / dNTP concentration, NaN / ∞).</exception>
    /// <exception cref="ArgumentException">PRIMER_MIN_QUALITY / a quality weight without a template, or a GC weight without
    /// <see cref="PrimerParameters.OptimalGcPercent"/> (Primer3 <c>_pr_data_control</c>).</exception>
    public static PrimerCandidate EvaluatePrimer(
        string sequence,
        int position,
        bool isForward,
        PrimerParameters? parameters = null)
    {
        ArgumentNullException.ThrowIfNull(sequence);
        var param = parameters ?? DefaultParameters;
        param.ValidateConditions(nameof(parameters));
        // No template, hence no SEQUENCE_QUALITY: Primer3 rejects PRIMER_MIN_QUALITY / PRIMER_WT_SEQ_QUAL without it.
        ValidatePrimer3Quality(null, 0, param.MinQuality, 0, param.EffectiveQualityRangeMin, param.EffectiveQualityRangeMax,
            param.PenaltyWeights?.SequenceQuality ?? 0.0, 0.0, nameof(parameters));
        param.ValidateGcOptimum(nameof(parameters));
        return EvaluatePrimerCore(sequence, position, isForward, param).Candidate;
    }

    // The single target of a pair search when the position penalties are not Primer3's defaults
    // (compute_position_penalty): [Start, End) 0-based, PRIMER_INSIDE_PENALTY / PRIMER_OUTSIDE_PENALTY.
    private readonly record struct TargetPosition(int Start, int End, double InsidePenalty, double OutsidePenalty);

    // The per-primer weights with the structure mode of the screen (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT).
    private static Primer3PenaltyWeights EffectivePenaltyWeights(PrimerParameters param) =>
        (param.PenaltyWeights ?? DefaultPrimer3Weights) with
        {
            ThermodynamicOligoAlignment = param.StructureScreen != PrimerStructureScreen.Primer3Alignment,
            ThermodynamicTemplateAlignment = param.ThermodynamicTemplateAlignment,
        };

    // primer_mispriming_to_template[_thermod]: max_template_mispriming[_th] ≥ 0 and the max score strictly greater.
    private static bool TemplateMisprimingExceeds(TemplateMisprimingScore score, PrimerParameters param)
    {
        double max = param.ActiveMaxTemplateMispriming;
        return max >= 0 && score.Max > max;
    }

    private static string TemplateMisprimingIssue(TemplateMisprimingScore score, PrimerParameters param) =>
        string.Create(System.Globalization.CultureInfo.InvariantCulture,
            $"Template mispriming score {score.Max:0.00} exceeds {param.ActiveMaxTemplateMispriming:0.00} (Primer3 PRIMER_MAX_TEMPLATE_MISPRIMING{(param.ThermodynamicTemplateAlignment ? "_TH" : "")})");

    // Per-primer secondary-structure screen. Primer3Thermodynamic: Primer3's default
    // (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1) ntthal limits PRIMER_MAX_SELF_ANY_TH / _SELF_END_TH /
    // _HAIRPIN_TH = 47 °C (a non-ACGT primer has no ntthal structure; it is already invalid by Tm).
    // Primer3Alignment: PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0 oligo_compl — dpal self_any > PRIMER_MAX_SELF_ANY,
    // then self_end > PRIMER_MAX_SELF_END (no hairpin value in this mode).
    // Heuristic: the sequence-only library stem-loop screen HasHairpinPotential (default stem 4 = unsourced library
    // threshold, loop 3 = thal.c min_hrpn_loop); not Primer3's hairpin screen.
    private static (bool HasHairpin, StructureValues Structure) AddStructureIssues(
        string seq, PrimerParameters param, List<string> issues)
    {
        if (param.StructureScreen == PrimerStructureScreen.Heuristic)
        {
            bool hp = HasHairpinPotential(seq);
            if (hp)
                issues.Add("Potential hairpin structure detected");
            return (hp, default);
        }

        var values = ComputeStructureValues(seq, param);
        if (param.StructureScreen == PrimerStructureScreen.Primer3Alignment)
        {
            double maxAny = param.EffectiveMaxSelfAny, maxEnd = param.EffectiveMaxSelfEnd;
            if (values.SelfAny > maxAny)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"Self-complementarity {values.SelfAny:0.00} exceeds {maxAny:0.00} (Primer3 PRIMER_MAX_SELF_ANY)"));
            if (values.SelfEnd > maxEnd)
                issues.Add(string.Create(System.Globalization.CultureInfo.InvariantCulture,
                    $"3' self-complementarity {values.SelfEnd:0.00} exceeds {maxEnd:0.00} (Primer3 PRIMER_MAX_SELF_END)"));
            return (false, values);
        }

        if (values.Thermo is not { } v)
            return (false, values);
        double max = param.EffectiveMaxStructureTm;
        bool hasHairpin = v.HairpinTh > max;
        if (hasHairpin)
            issues.Add($"Hairpin melting temperature {v.HairpinTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_HAIRPIN_TH)");
        if (v.SelfAnyTh > max)
            issues.Add($"Self-dimer melting temperature {v.SelfAnyTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_SELF_ANY_TH)");
        if (v.SelfEndTh > max)
            issues.Add($"3' self-dimer melting temperature {v.SelfEndTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_SELF_END_TH)");
        return (hasHairpin, values);
    }

    /// <summary>
    /// Computes the intermolecular <b>self-dimer</b> melting temperature (°C) of a DNA oligo: the
    /// bimolecular Tm of the most stable duplex it forms with a second copy of itself, via the
    /// Primer3 / <c>ntthal</c> thermodynamic alignment (SantaLucia &amp; Hicks 2004 unified NN).
    /// <b>Opt-in</b>: the perfect-match duplex / hairpin / default Tm methods are unchanged.
    /// </summary>
    /// <param name="sequence">DNA oligo (5'→3'); ACGT only (case-insensitive).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The self-dimer Tm in °C (ntthal also reports weak / single-pair structures, whose Tm
    /// may be negative), or <c>double.NaN</c> if ntthal finds no structure / the sequence is
    /// null, empty or non-ACGT.</returns>
    /// <exception cref="ArgumentException">The sequence is longer than 60 nt (thal.c THAL_MAX_ALIGN).</exception>
    public static double CalculateSelfDimerMeltingTemperature(
        string sequence,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar) =>
        CalculateDimerMeltingTemperature(sequence, sequence, sodiumMolar, strandConcentrationMolar);

    /// <summary>
    /// As <see cref="CalculatePrimer3OligoStructure(string, double, double, double, double)"/> with thal.h
    /// <c>THAL_MAX_ALIGN</c> raised to <paramref name="maxAlignLength"/> (opt-in; see <see cref="NtthalMaxAlignLength"/>):
    /// self_any / self_end / hairpin of oligos up to <paramref name="maxAlignLength"/> nt by the unchanged ntthal
    /// recursions (= thal.c compiled with <c>-DTHAL_MAX_ALIGN=maxAlignLength</c>). Primer3 itself never sees such
    /// oligos (PRIMER_MAX_SIZE ≤ 35 / internal ≤ 36 nt and THAL_MAX_ALIGN = 60).
    /// </summary>
    /// <param name="maxAlignLength">THAL_MAX_ALIGN, 60 … 10 000.</param>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="maxAlignLength"/> outside 60–10 000.</exception>
    /// <exception cref="ArgumentException">The oligo is longer than <paramref name="maxAlignLength"/>.</exception>
    public static Primer3OligoStructure? CalculatePrimer3OligoStructure(
        string primer,
        double monovalentMillimolar,
        double divalentMillimolar,
        double dntpMillimolar,
        double dnaConcentrationNanomolar,
        int maxAlignLength) =>
        CalculatePrimer3OligoStructureCore(primer, monovalentMillimolar, divalentMillimolar, dntpMillimolar,
            dnaConcentrationNanomolar, CheckedMaxAlignLength(maxAlignLength));

    /// <summary>
    /// Computes Primer3's thermodynamic pair complementarity between a left (forward) and a right
    /// (reverse) primer, both given 5′→3′ as synthesised, exactly as <c>characterize_pair</c> does in
    /// the default thermodynamic mode: compl_any = ntthal ANY(left, right); compl_end = the maximum of
    /// ntthal END1(left, right), END2(left, right), END1(rc(right), rc(left)) and END2(rc(right), rc(left))
    /// (Primer3 evaluates the last two on the reverse complements, <c>align_thermod(s2, s1_rev, …)</c>).
    /// Tm in °C, 0 when no structure or Tm &lt; 0 °C. Conditions default to Primer3's primer conditions.
    /// </summary>
    /// <param name="leftPrimer">Forward primer (5′→3′), ACGT only.</param>
    /// <param name="rightPrimer">Reverse primer (5′→3′), ACGT only.</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM.</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM.</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM.</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM.</param>
    /// <returns>The pair values, or <c>null</c> when either primer is null/empty or contains a
    /// non-ACGT character.</returns>
    public static Primer3PairComplementarity? CalculatePrimer3PairComplementarity(
        string leftPrimer,
        string rightPrimer,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar,
        double dnaConcentrationNanomolar = Primer3DnaConcentrationNanomolar)
    {
        if (!IsAcgtOnly(leftPrimer) || !IsAcgtOnly(rightPrimer))
            return null;
        string l = leftPrimer.ToUpperInvariant(), r = rightPrimer.ToUpperInvariant();
        double mv = monovalentMillimolar / 1000.0, dv = divalentMillimolar / 1000.0, dntp = dntpMillimolar / 1000.0;
        double conc = dnaConcentrationNanomolar * 1e-9;
        var (any, end) = Primer3PairTms(l, r, mv, dv, dntp, conc, double.PositiveInfinity);
        return new Primer3PairComplementarity(any, end);
    }

    private static double CalculatePrimerScore(string seq, double gc, double tm, int homopolymer, PrimerParameters param)
    {
        double score = 100;

        // Penalize for deviation from optimal length
        score -= Math.Abs(seq.Length - param.OptimalLength) * 2;

        // Penalize for deviation from optimal Tm
        score -= Math.Abs(tm - param.OptimalTm) * 2;

        // Penalize for deviation from 50% GC
        score -= Math.Abs(gc - 50) * 0.5;

        // Penalize for homopolymers
        score -= homopolymer * 5;

        // Bonus for GC clamp at 3' end
        if (seq.Length >= 2)
        {
            char last = seq[^1];
            if (last == 'G' || last == 'C')
                score += 5;
        }

        return Math.Max(0, score);
    }

    /// <summary>
    /// Designs the best forward/reverse primer pair for a target region with Primer3's pair-selection
    /// algorithm (Untergasser et al. 2012; <c>libprimer3.cc</c> <c>make_detection_primer_lists</c>,
    /// <c>choose_pair_or_triple</c>, <c>characterize_pair</c>, <c>obj_fn</c>, <c>compare_primer_pair</c>)
    /// — the first pair (rank 0) of <see cref="DesignPrimerPairs"/>:
    /// <list type="number">
    /// <item><b>Search region</b> (Primer3 <c>SEQUENCE_TARGET</c> / <c>SEQUENCE_INCLUDED_REGION</c> /
    /// <c>PRIMER_PRODUCT_SIZE_RANGE</c>): forward candidates end at or before <paramref name="targetStart"/>,
    /// reverse candidates start at or after <paramref name="targetEnd"/> (primers never overlap the target),
    /// both inside <see cref="PrimerPairOptions.IncludedRegion"/> (default: the whole template); the
    /// product-size ranges (default 100–300 bp, Primer3's default) bound the product. Every candidate is
    /// evaluated by <see cref="EvaluatePrimer"/> (per-primer limits of <paramref name="parameters"/>) and
    /// carries its Primer3 per-primer penalty (<see cref="PrimerCandidate.Penalty"/>).</item>
    /// <item><b>Pair search</b>: candidates sorted as Primer3's <c>sort_primer_array</c>; pairs are examined
    /// with Primer3's pruning and must satisfy, in <c>characterize_pair</c> order, the product size range
    /// (ranges are tried in order: the next range is used only when no pair fits the current one), the
    /// product Tm limits (<see cref="PrimerPairOptions.ProductMinTm"/>/<see cref="PrimerPairOptions.ProductMaxTm"/>,
    /// Primer3 <c>long_seq_tm</c>), |Tm_f − Tm_r| ≤ <see cref="PrimerPairOptions.MaxTmDifference"/>
    /// (PRIMER_PAIR_MAX_DIFF_TM), the per-primer secondary-structure screen and the pair complementarity
    /// screen — with the default <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> Primer3's ntthal
    /// limits (self-any/self-end/hairpin and pair compl-any/compl-end Tm ≤ 47 °C), evaluated lazily as
    /// <c>characterize_pair</c> does; with <see cref="PrimerStructureScreen.Heuristic"/>
    /// <see cref="HasHairpinPotential"/> and <see cref="HasPrimerDimer"/> — and, with
    /// <see cref="PrimerPairOptions.PickInternalOligo"/>, a hybridization oligo between the primers
    /// (<c>choose_internal_oligo</c>).</item>
    /// <item><b>Objective</b>: the pair penalty is Primer3's <c>obj_fn</c> with
    /// <see cref="PrimerPairOptions.Weights"/> (default: PRIMER_PAIR_WT_PR_PENALTY = 1, all other pair weights
    /// 0, i.e. the sum of the two primer penalties); ties within 1e-6 are broken as
    /// <c>compare_primer_pair</c> (left primer further 3′, then right primer further 5′, then shorter left,
    /// then shorter right).</item>
    /// </list>
    /// Coordinates: the target is the half-open interval [<paramref name="targetStart"/>,
    /// <paramref name="targetEnd"/>) (0-based). Reverse-primer <see cref="PrimerCandidate.Position"/> is
    /// the leftmost template coordinate of its binding site; the product size is
    /// <c>reverse.Position + reverse.Length − forward.Position</c> (Primer3 PRIMER_PAIR_PRODUCT_SIZE).
    /// When candidates exist but no pair satisfies the pair constraints, the individually lowest-penalty
    /// forward and reverse candidates are returned with <c>IsValid = false</c> and a message naming the
    /// violated constraint. Library defaults that differ from Primer3 (documented): the per-primer limits of
    /// <see cref="DefaultParameters"/> and PRIMER_PAIR_MAX_DIFF_TM = 5 °C (Primer3: 100); pass
    /// <see cref="Primer3DefaultParameters"/> and <see cref="PrimerPairOptions.Primer3Defaults"/> for
    /// Primer3's defaults. Reaction conditions: the primer Tm, the ntthal structure / pair complementarity values and
    /// the product Tm use the primer conditions of <paramref name="parameters"/> (PRIMER_SALT_MONOVALENT,
    /// PRIMER_SALT_DIVALENT, PRIMER_DNTP_CONC, PRIMER_DNA_CONC; Primer3 defaults 50 mM / 1.5 mM / 0.6 mM / 50 nM), the
    /// internal oligo those of <see cref="PrimerPairOptions.InternalOligo"/> (PRIMER_INTERNAL_*; defaults 50 mM / 0 / 0 /
    /// 50 nM). Verified against primer3-py 2.3.1 <c>design_primers</c>, including random reaction conditions.
    /// </summary>
    /// <param name="template">The DNA template sequence.</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Per-primer design parameters (default <see cref="DefaultParameters"/>).</param>
    /// <param name="pairOptions">Pair options (default <see cref="PrimerPairOptions.Default"/>).</param>
    /// <returns>Primer pair result.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentException">Invalid target or options (see <see cref="PrimerPairOptions"/>).</exception>
    public static PrimerPairResult DesignPrimers(
        DnaSequence template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null,
        PrimerPairOptions? pairOptions = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        return DesignPrimersCore(template.Sequence, null, targetStart, targetEnd, parameters, pairOptions);
    }

    /// <summary>
    /// <see cref="DesignPrimers(DnaSequence, int, int, PrimerParameters?, PrimerPairOptions?)"/> on a case-preserving
    /// template string (Primer3 SEQUENCE_TEMPLATE as given): with <see cref="PrimerPairOptions.LowercaseMasking"/>
    /// (PRIMER_LOWERCASE_MASKING) or <see cref="PrimerPairOptions.MaskTemplate"/> a primer or internal oligo whose 3′-terminal
    /// template base is lower case is rejected (Primer3 <c>is_lowercase_masked</c>); otherwise case is ignored and the result
    /// equals the <see cref="DnaSequence"/> overload. All other computations use the upper-cased template.
    /// </summary>
    /// <param name="template">Template (A/C/G/T, any case; lower case marks masked bases).</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Per-primer design parameters (default <see cref="DefaultParameters"/>).</param>
    /// <param name="pairOptions">Pair options (default <see cref="PrimerPairOptions.Default"/>).</param>
    /// <returns>Primer pair result.</returns>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentException">A character other than A/C/G/T, an invalid target or options.</exception>
    public static PrimerPairResult DesignPrimers(
        string template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null,
        PrimerPairOptions? pairOptions = null)
    {
        ArgumentNullException.ThrowIfNull(template);
        return DesignPrimersCore(new DnaSequence(template).Sequence, template, targetStart, targetEnd, parameters, pairOptions);
    }

    // Primer3 secondary-structure values of one primer: ntthal Tm values (thermodynamic screen) or dpal
    // self_any / self_end scores (alignment screen).
    private readonly record struct StructureValues(Primer3OligoStructure? Thermo, double? SelfAny, double? SelfEnd);

    private static IEnumerable<PrimerCandidate> GeneratePrimerCandidatesIterator(
        DnaSequence template, int regionStart, int regionEnd, bool forward, PrimerParameters param)
    {
        for (int start = regionStart; start + param.MinLength <= regionEnd; start++)
        {
            for (int len = param.MinLength; len <= param.MaxLength && start + len <= regionEnd; len++)
            {
                var seq = template.Sequence.Substring(start, len);
                if (!forward)
                    seq = DnaSequence.GetReverseComplementString(seq);

                yield return EvaluatePrimer(seq, start, forward, param);
            }
        }
    }

    /// <summary>
    /// Generates all possible primers for a region.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="template"/> is null.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="regionStart"/> &lt; 0, or
    /// <paramref name="regionStart"/> / <paramref name="regionEnd"/> beyond the end of the template
    /// (<paramref name="regionEnd"/> is exclusive; an empty region yields no candidates).</exception>
    public static IEnumerable<PrimerCandidate> GeneratePrimerCandidates(
        DnaSequence template,
        int regionStart,
        int regionEnd,
        bool forward = true,
        PrimerParameters? parameters = null)
    {
        // Validated eagerly (audit round 7, A7-3/A7-4): out-of-template regions crashed in Substring on enumeration.
        ArgumentNullException.ThrowIfNull(template);
        ArgumentOutOfRangeException.ThrowIfNegative(regionStart);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(regionStart, template.Length);
        ArgumentOutOfRangeException.ThrowIfGreaterThan(regionEnd, template.Length);
        return GeneratePrimerCandidatesIterator(template, regionStart, regionEnd, forward, parameters ?? DefaultParameters);
    }
}
