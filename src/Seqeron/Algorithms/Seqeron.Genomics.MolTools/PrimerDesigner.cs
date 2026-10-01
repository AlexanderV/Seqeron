using System.Text;

namespace Seqeron.Genomics.MolTools;

/// <summary>
/// Designs PCR primers for DNA sequences with various quality criteria.
/// </summary>
public static class PrimerDesigner
{
    /// <summary>
    /// Default primer design parameters.
    /// </summary>
    public static readonly PrimerParameters DefaultParameters = new(
        MinLength: 18,
        MaxLength: 25,
        OptimalLength: 20,
        MinGcContent: 40,
        MaxGcContent: 60,
        MinTm: 57,
        MaxTm: 63,
        OptimalTm: 60,
        MaxHomopolymer: 4,
        MaxDinucleotideRepeats: 4,
        Avoid3PrimeGC: false,
        Check3PrimeStability: true
    );

    /// <summary>
    /// Designs a forward/reverse primer pair flanking a target region, following Primer3's
    /// pair-selection semantics (Untergasser et al. 2012; <c>libprimer3.cc</c>
    /// <c>choose_pair_or_triple</c> / <c>compare_primer_pair</c>):
    /// <list type="number">
    /// <item>every forward candidate that ends at or before <paramref name="targetStart"/> (within
    /// <see cref="PrimerSearchFlank"/> bp) and every reverse candidate that starts at or after
    /// <paramref name="targetEnd"/> (within <see cref="PrimerSearchFlank"/> bp) is evaluated by
    /// <see cref="EvaluatePrimer"/>; only candidates passing every per-primer constraint are kept
    /// (Primer3 <c>SEQUENCE_TARGET</c>: primers never overlap the target);</item>
    /// <item>each kept candidate carries its Primer3 per-primer penalty (<see cref="PrimerCandidate.Penalty"/>,
    /// Primer3 <c>p_obj_fn</c> with default weights = |Tm − OptimalTm| + |length − OptimalLength|);</item>
    /// <item>the pair returned is the one with the <b>lowest pair penalty</b> (sum of the two primer
    /// penalties, Primer3 default <c>PRIMER_PAIR_WT_PR_PENALTY = 1</c>, all other pair weights 0) among
    /// all pairs satisfying the pair constraints |Tm_f − Tm_r| ≤ <see cref="MaxPairTmDifference"/> °C and
    /// no primer-dimer, and whose primers pass the per-primer secondary-structure screen. With the default
    /// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> these are Primer3's default ntthal limits
    /// (self-any/self-end/hairpin Tm and pair compl-any/compl-end Tm ≤ <see cref="Primer3MaxStructureTm"/>
    /// = 47 °C; <see cref="CalculatePrimer3OligoStructure"/>, <see cref="CalculatePrimer3PairComplementarity"/>),
    /// evaluated lazily in the pair loop as Primer3's <c>characterize_pair</c> does; with
    /// <see cref="PrimerStructureScreen.Heuristic"/> they are <see cref="HasHairpinPotential"/> and
    /// <see cref="HasPrimerDimer"/>. Ties (within 1e-6) are broken exactly as
    /// Primer3's <c>compare_primer_pair</c>: left primer further 3' (right), then right primer further
    /// 5' (left), then shorter left, then shorter right.</item>
    /// </list>
    /// Coordinates: the target is the half-open interval [<paramref name="targetStart"/>,
    /// <paramref name="targetEnd"/>) (0-based). Reverse-primer <see cref="PrimerCandidate.Position"/> is
    /// the leftmost template coordinate of its binding site; the product size is
    /// <c>reverse.Position + reverse.Length − forward.Position</c> (Primer3 PRIMER_PAIR_PRODUCT_SIZE).
    /// When candidates exist on both sides but no pair satisfies the pair constraints, the
    /// individually lowest-penalty forward and reverse candidates are returned with
    /// <c>IsValid = false</c> and a message naming the violated constraint.
    /// Deviations from Primer3 defaults (documented): pair ΔTm ≤ 5 °C (Primer3 PRIMER_PAIR_MAX_DIFF_TM
    /// = 100) and no PRIMER_PRODUCT_SIZE_RANGE (the ±200 bp search flanks bound the product instead);
    /// the per-primer limits are those of <paramref name="parameters"/>.
    /// </summary>
    /// <param name="template">The DNA template sequence.</param>
    /// <param name="targetStart">0-based inclusive start of the target region.</param>
    /// <param name="targetEnd">0-based exclusive end of the target region; must be &lt; template length.</param>
    /// <param name="parameters">Primer design parameters (optional).</param>
    /// <returns>Primer pair result.</returns>
    public static PrimerPairResult DesignPrimers(
        DnaSequence template,
        int targetStart,
        int targetEnd,
        PrimerParameters? parameters = null)
    {
        var param = parameters ?? DefaultParameters;

        if (targetStart < 0 || targetEnd >= template.Length || targetStart >= targetEnd)
            throw new ArgumentException("Invalid target region.");

        // Forward candidates (upstream of target, may not overlap it).
        var forwardCandidates = new List<(PrimerCandidate C, double Tm)>();
        int forwardSearchStart = Math.Max(0, targetStart - PrimerSearchFlank);
        for (int start = forwardSearchStart; start < targetStart; start++)
        {
            for (int len = param.MinLength; len <= param.MaxLength && start + len <= targetStart; len++)
            {
                var (candidate, tm) = EvaluatePrimerCore(template.Sequence.Substring(start, len), start, true, param, evaluateStructure: false);
                if (candidate.IsValid)
                    forwardCandidates.Add((candidate, tm));
            }
        }

        // Reverse candidates (downstream of target, evaluated as the reverse complement).
        var reverseCandidates = new List<(PrimerCandidate C, double Tm)>();
        int reverseSearchEnd = Math.Min(template.Length, targetEnd + PrimerSearchFlank);
        for (int end = targetEnd + param.MinLength; end <= reverseSearchEnd; end++)
        {
            for (int len = param.MinLength; len <= param.MaxLength && end - len >= targetEnd; len++)
            {
                int start = end - len;
                var revComp = DnaSequence.GetReverseComplementString(template.Sequence.Substring(start, len));
                var (candidate, tm) = EvaluatePrimerCore(revComp, start, false, param, evaluateStructure: false);
                if (candidate.IsValid)
                    reverseCandidates.Add((candidate, tm));
            }
        }

        if (forwardCandidates.Count == 0 || reverseCandidates.Count == 0)
        {
            return new PrimerPairResult(
                null, null, false,
                "Could not find valid primers for the target region.",
                0
            );
        }

        // Primer3 examines primers in increasing penalty order (sort_primer_array).
        forwardCandidates.Sort((a, b) => CompareLeft(a.C, b.C));
        reverseCandidates.Sort((a, b) => CompareRight(a.C, b.C));

        // Secondary-structure screen of individual primers, run lazily and cached exactly where
        // Primer3 runs it (characterize_pair: the "expensive" per-primer checks are postponed until a
        // primer takes part in a pair that passed the cheaper pair checks).
        var forwardStructureOk = new bool?[forwardCandidates.Count];
        var reverseStructureOk = new bool?[reverseCandidates.Count];
        // Results depend only on the sequences, so they are also cached by sequence (templates with
        // repeats yield many candidates with the same sequence at different positions).
        var structureBySequence = new Dictionary<string, bool>(StringComparer.Ordinal);
        var dimerBySequencePair = new Dictionary<(string F, string R), bool>();
        bool StructureOk(List<(PrimerCandidate C, double Tm)> list, bool?[] cache, int i)
        {
            if (cache[i] is { } known)
                return known;
            string seq = list[i].C.Sequence;
            if (!structureBySequence.TryGetValue(seq, out bool ok))
            {
                var issues = new List<string>();
                AddStructureIssues(seq, param, issues);
                ok = issues.Count == 0;
                structureBySequence[seq] = ok;
            }
            cache[i] = ok;
            return ok;
        }
        bool FormsDimer(string f, string r)
        {
            if (!dimerBySequencePair.TryGetValue((f, r), out bool dimer))
            {
                dimer = PairFormsDimer(f, r, param);
                dimerBySequencePair[(f, r)] = dimer;
            }
            return dimer;
        }

        int bestFi = -1, bestRi = -1;
        double bestQuality = double.PositiveInfinity;
        bool sawTmFailure = false, sawDimerFailure = false;

        for (int ri = 0; ri < reverseCandidates.Count; ri++)
        {
            var r = reverseCandidates[ri];
            for (int fi = 0; fi < forwardCandidates.Count; fi++)
            {
                var f = forwardCandidates[fi];
                double quality = f.C.Penalty + r.C.Penalty;
                // choose_pair_or_triple: no later forward primer can improve on the best pair.
                if (quality > bestQuality)
                    break;

                if (Math.Abs(f.Tm - r.Tm) > MaxPairTmDifference)
                {
                    sawTmFailure = true;
                    continue;
                }
                if (!StructureOk(forwardCandidates, forwardStructureOk, fi))
                    continue;
                if (!StructureOk(reverseCandidates, reverseStructureOk, ri))
                    break; // this reverse primer fails on its own; no pair with it can be valid
                if (FormsDimer(f.C.Sequence, r.C.Sequence))
                {
                    sawDimerFailure = true;
                    continue;
                }

                if (bestFi < 0 || ComparePair(quality, f.C, r.C, bestQuality,
                        forwardCandidates[bestFi].C, reverseCandidates[bestRi].C) < 0)
                {
                    bestFi = fi;
                    bestRi = ri;
                    bestQuality = quality;
                }
            }
        }

        if (bestFi < 0)
        {
            // The individually lowest-penalty primers that pass their own (structure) constraints.
            int f0i = FirstStructurallyValid(forwardCandidates, forwardStructureOk);
            int r0i = FirstStructurallyValid(reverseCandidates, reverseStructureOk);
            if (f0i < 0 || r0i < 0)
            {
                return new PrimerPairResult(
                    null, null, false,
                    "Could not find valid primers for the target region.",
                    0
                );
            }
            var f0 = Reevaluate(forwardCandidates[f0i].C);
            var r0 = Reevaluate(reverseCandidates[r0i].C);
            string reason;
            if (sawTmFailure && !sawDimerFailure)
                reason = $"No primer pair within the {MaxPairTmDifference:F0}°C Tm-difference limit (best primers: Tm {f0.MeltingTemperature:F1}/{r0.MeltingTemperature:F1}°C).";
            else if (sawDimerFailure && !sawTmFailure)
                reason = "Every candidate primer pair forms a primer-dimer.";
            else
                reason = $"No primer pair satisfies both the {MaxPairTmDifference:F0}°C Tm-difference limit and primer-dimer avoidance.";
            return new PrimerPairResult(
                Forward: f0,
                Reverse: r0,
                IsValid: false,
                Message: reason,
                ProductSize: r0.Position + r0.Length - f0.Position);
        }

        var forward = Reevaluate(forwardCandidates[bestFi].C);
        var reverse = Reevaluate(reverseCandidates[bestRi].C);
        return new PrimerPairResult(
            Forward: forward,
            Reverse: reverse,
            IsValid: true,
            Message: "Valid primer pair found.",
            ProductSize: reverse.Position + reverse.Length - forward.Position
        );

        int FirstStructurallyValid(List<(PrimerCandidate C, double Tm)> list, bool?[] cache)
        {
            for (int i = 0; i < list.Count; i++)
                if (StructureOk(list, cache, i))
                    return i;
            return -1;
        }

        // Full evaluation (including the structure values) of a chosen primer.
        PrimerCandidate Reevaluate(PrimerCandidate c) =>
            EvaluatePrimerCore(c.Sequence, c.Position, c.IsForward, param).Candidate;
    }

    /// <summary>Flank (bp) searched on each side of the target for primer candidates.</summary>
    public const int PrimerSearchFlank = 200;

    /// <summary>
    /// Maximum |Tm_forward − Tm_reverse| (°C) accepted for a primer pair (Addgene "within 5 °C";
    /// stricter than Primer3's PRIMER_PAIR_MAX_DIFF_TM default of 100).
    /// </summary>
    public const double MaxPairTmDifference = 5.0;

    // Primer3 compare_primer_pair epsilon on pair_quality.
    private const double PairQualityEpsilon = 1e-6;

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

    /// <summary>
    /// Evaluates a single primer candidate against the per-primer constraints of
    /// <paramref name="parameters"/> (length, GC%, Tm, homopolymer, dinucleotide repeat, secondary
    /// structure, 3'-end stability, optional GC clamp). The secondary-structure screen is
    /// <see cref="PrimerParameters.StructureScreen"/>: by default Primer3's thermodynamic limits
    /// (ntthal self-dimer, 3′ self-dimer and hairpin Tm ≤ 47 °C, reported in
    /// <see cref="PrimerCandidate.SelfAnyTh"/>/<see cref="PrimerCandidate.SelfEndTh"/>/<see cref="PrimerCandidate.HairpinTh"/>),
    /// or the sequence-only <see cref="HasHairpinPotential"/>. The Tm is Primer3's default primer Tm
    /// (<see cref="CalculateMeltingTemperaturePrimer3"/>: SantaLucia 1998 nearest-neighbour,
    /// SantaLucia salt correction, 50 mM monovalent, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM oligo), the
    /// scale on which the Primer3-sourced Tm window 57–63 °C (opt 60) is defined. A sequence
    /// containing a non-ACGT base has no computable Tm (Primer3 PRIMER_MAX_NS_ACCEPTED = 0): it is
    /// reported with Tm 0 and an issue. <see cref="PrimerCandidate.Penalty"/> is the Primer3
    /// per-primer penalty (<see cref="CalculatePrimer3Penalty"/>, default weights, optima
    /// OptimalTm / OptimalLength) used by <see cref="DesignPrimers"/> for ranking;
    /// <see cref="PrimerCandidate.Score"/> is an informational additive quality score (0–100,
    /// higher is better) that does not drive selection.
    /// </summary>
    public static PrimerCandidate EvaluatePrimer(
        string sequence,
        int position,
        bool isForward,
        PrimerParameters? parameters = null) =>
        EvaluatePrimerCore(sequence, position, isForward, parameters ?? DefaultParameters).Candidate;

    // Evaluates a candidate and also returns its unrounded Tm (Primer3 compares unrounded Tm values).
    // With evaluateStructure = false the secondary-structure screen is skipped (DesignPrimers runs it
    // lazily, like Primer3's characterize_pair, and re-evaluates the chosen primers in full).
    private static (PrimerCandidate Candidate, double Tm) EvaluatePrimerCore(
        string sequence,
        int position,
        bool isForward,
        PrimerParameters param,
        bool evaluateStructure = true)
    {
        var seq = sequence.ToUpperInvariant();

        double gcContent = CalculateGcContent(seq);
        double tmRaw = CalculateMeltingTemperaturePrimer3(seq);
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

        if (homopolymer > param.MaxHomopolymer)
            issues.Add($"Homopolymer run of {homopolymer} exceeds max {param.MaxHomopolymer}");

        if (dinucRepeat > param.MaxDinucleotideRepeats)
            issues.Add($"Dinucleotide repeat of {dinucRepeat} exceeds max {param.MaxDinucleotideRepeats}");

        bool hasHairpin = false;
        Primer3OligoStructure? structure = null;
        if (evaluateStructure)
            (hasHairpin, structure) = AddStructureIssues(seq, param, issues);

        if (param.Check3PrimeStability && stability3Prime < -9)
            issues.Add($"3' end too stable (ΔG = {stability3Prime:F1} kcal/mol)");

        // Check 3' end for GC clamp
        if (param.Avoid3PrimeGC && seq.Length >= 2)
        {
            string last2 = seq.Substring(seq.Length - 2);
            int gcCount = last2.Count(c => c == 'G' || c == 'C');
            if (gcCount == 0)
                issues.Add("No GC clamp at 3' end");
        }

        bool isValid = issues.Count == 0;

        // Informational heuristic score and the Primer3 ranking penalty.
        double score = CalculatePrimerScore(seq, gcContent, tm, homopolymer, param);
        double penalty = CalculatePrimer3Penalty(
            new Primer3PenaltyInputs(tm, seq.Length, gcContent),
            DefaultPrimer3Weights,
            new Primer3Optima(param.OptimalTm, param.OptimalLength, DefaultPrimer3Optima.OptGcPercent));

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
            SelfAnyTh: structure?.SelfAnyTh,
            SelfEndTh: structure?.SelfEndTh,
            HairpinTh: structure?.HairpinTh
        );
        return (candidate, tm);
    }

    // Per-primer secondary-structure screen. Primer3Thermodynamic: Primer3's default
    // (PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1) ntthal limits PRIMER_MAX_SELF_ANY_TH / _SELF_END_TH /
    // _HAIRPIN_TH = 47 °C (a non-ACGT primer has no ntthal structure; it is already invalid by Tm).
    // Heuristic: the sequence-only stem-loop screen HasHairpinPotential (default stem 4, loop 3).
    private static (bool HasHairpin, Primer3OligoStructure? Structure) AddStructureIssues(
        string seq, PrimerParameters param, List<string> issues)
    {
        if (param.StructureScreen == PrimerStructureScreen.Heuristic)
        {
            bool hp = HasHairpinPotential(seq);
            if (hp)
                issues.Add("Potential hairpin structure detected");
            return (hp, null);
        }

        var st = CalculatePrimer3OligoStructure(seq);
        if (st is null)
            return (false, null);
        var v = st.Value;
        double max = param.EffectiveMaxStructureTm;
        bool hasHairpin = v.HairpinTh > max;
        if (hasHairpin)
            issues.Add($"Hairpin melting temperature {v.HairpinTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_HAIRPIN_TH)");
        if (v.SelfAnyTh > max)
            issues.Add($"Self-dimer melting temperature {v.SelfAnyTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_SELF_ANY_TH)");
        if (v.SelfEndTh > max)
            issues.Add($"3' self-dimer melting temperature {v.SelfEndTh:F1}°C exceeds {max:0.##}°C (Primer3 PRIMER_MAX_SELF_END_TH)");
        return (hasHairpin, st);
    }

    // Pair-level complementarity screen used by DesignPrimers (Primer3 characterize_pair).
    private static bool PairFormsDimer(string forward, string reverse, PrimerParameters param)
    {
        if (param.StructureScreen == PrimerStructureScreen.Heuristic)
            return HasPrimerDimer(forward, reverse);
        if (!IsAcgtOnly(forward) || !IsAcgtOnly(reverse))
            return false;
        // Same values as CalculatePrimer3PairComplementarity(...).Exceeds(max), stopping at the first
        // alignment over the limit (characterize_pair also fails the pair on compl_any first).
        var (any, end) = Primer3PairTms(forward.ToUpperInvariant(), reverse.ToUpperInvariant(),
            0.050, 0.0015, 0.0006, 50e-9, param.EffectiveMaxStructureTm);
        return any > param.EffectiveMaxStructureTm || end > param.EffectiveMaxStructureTm;
    }

    /// <summary>
    /// Calculates the "basic" melting temperature for DNA primers (OligoCalc basic Tm,
    /// Kibbe 2007, NAR 35:W43): Wallace rule Tm = 2(A+T) + 4(G+C) (Thein &amp; Wallace 1986) for
    /// &lt; 14 counted bases, and Tm = 64.9 + 41·(G+C − 16.4)/N (customarily attributed to
    /// Marmur &amp; Doty 1962) for ≥ 14 counted bases. Both formulas assume fixed standard conditions
    /// (50 nM primer, 50 mM Na+, pH 7.0); use <see cref="CalculateMeltingTemperatureWithSalt"/>
    /// for another [Na+], or <see cref="CalculateMeltingTemperatureNN"/> for a nearest-neighbor Tm.
    /// A, C, G, T and U are counted (case-insensitive; U read as T, as Biopython <c>MeltingTemp._check</c>
    /// back-transcribes RNA for <c>Tm_Wallace</c>/<c>Tm_GC</c>); all other characters are ignored.
    /// Delegates to the canonical <see cref="ThermoConstants.CalculateBasicTm"/>.
    /// </summary>
    public static double CalculateMeltingTemperature(string primer) => ThermoConstants.CalculateBasicTm(primer);

    /// <summary>
    /// Calculates the OligoCalc "salt adjusted" melting temperature (Kibbe 2007, NAR 35:W43),
    /// the [Na+]-aware counterpart of <see cref="CalculateMeltingTemperature(string)"/>:
    /// <list type="bullet">
    /// <item>&lt; 14 valid bases: Tm = 2(A+T) + 4(G+C) − 16.6·log10(0.050) + 16.6·log10([Na+]) —
    /// the Wallace rule (defined at 50 mM Na+) shifted by the Schildkraut–Lifson relative correction;</item>
    /// <item>≥ 14 valid bases: Tm = 100.5 + 41·(G+C)/N − 820/N + 16.6·log10([Na+]).</item>
    /// </list>
    /// [Na+] enters in mol/L (converted from the mM argument). The basic formulas already assume
    /// 50 mM Na+, so the salt term is never added on top of them. Result rounded to one decimal.
    /// A/C/G/T/U are counted as in <see cref="CalculateMeltingTemperature(string)"/> (U read as T);
    /// returns 0 for null/empty input or no counted bases.
    /// </summary>
    /// <param name="primer">Primer sequence.</param>
    /// <param name="naConcentration">Na+ concentration in mM (default: 50). Must be &gt; 0.</param>
    /// <returns>Salt-adjusted melting temperature in °C.</returns>
    /// <exception cref="ArgumentOutOfRangeException">When <paramref name="naConcentration"/> is not a positive finite number.</exception>
    public static double CalculateMeltingTemperatureWithSalt(string primer, double naConcentration = 50)
    {
        if (!(naConcentration > 0) || double.IsInfinity(naConcentration))
            throw new ArgumentOutOfRangeException(nameof(naConcentration), naConcentration,
                "Na+ concentration must be a positive finite value in mM.");

        if (string.IsNullOrEmpty(primer))
            return 0;

        var (at, gc) = ThermoConstants.CountBasicTmBases(primer);
        if (at + gc == 0)
            return 0;

        return Math.Round(
            ThermoConstants.CalculateOligoCalcSaltAdjustedTm(at, gc, naConcentration / 1000.0), 1);
    }

    // ---- Primer3 default primer Tm (PRIMER-DESIGN-001) -------------------------------------
    // Reproduces Primer3's seqtm()/oligotm() with PRIMER_TM_FORMULA = 1 (SantaLucia 1998) and
    // PRIMER_SALT_CORRECTIONS = 1 (SantaLucia 1998), the Primer3 ≥ 2.0 defaults, including the
    // divalent→monovalent equivalence of von Ahsen et al. (2001) and the long_seq_tm fallback
    // for oligos longer than MAX_NN_TM_LENGTH (36). Source: primer3 src/oligotm.c (oligotm,
    // seqtm, long_seq_tm, divalent_to_monovalent, symmetry; SantaLucia_1998_dH/dS tables);
    // SantaLucia (1998) PNAS 95:1460 Table 2 (NN + terminal-initiation terms); von Ahsen et al.
    // (2001) Clin Chem 47:1956. Cross-checked against primer3-py 2.3.1 primer3.calc_tm.

    /// <summary>Primer3 MAX_NN_TM_LENGTH: longer oligos use the long_seq_tm GC/length formula.</summary>
    public const int Primer3MaxNnTmLength = 36;

    // SantaLucia (1998) Table 2 unified NN parameters exactly as tabulated in Primer3 oligotm.c
    // (SantaLucia_1998_dH in units of −100 cal/mol, SantaLucia_1998_dS in units of −0.1 cal/(K·mol);
    // index order A, C, G, T). Primer3 sums these as integers, so primers with the same
    // nearest-neighbour multiset get bit-identical Tm values (this matters for tie-breaking in
    // DesignPrimers). Differs from the SantaLucia & Hicks (2004) set (ThermoConstants NnParameterSet.SantaLuciaHicks2004) in
    // AA/TT (−7.9/−22.2 vs −7.6/−21.3).
    private static readonly int[,] Primer3SantaLucia1998Dh =
    {
        { 79, 84, 78, 72 },
        { 85, 80, 106, 78 },
        { 82, 98, 80, 84 },
        { 72, 82, 85, 79 }
    };
    private static readonly int[,] Primer3SantaLucia1998Ds =
    {
        { 222, 224, 210, 204 },
        { 227, 199, 272, 210 },
        { 222, 244, 199, 224 },
        { 213, 222, 227, 222 }
    };

    // SantaLucia (1998) Table 2 initiation per duplex end, in the same integer units
    // (oligotm.c "Terminal penalty"): terminal A·T ΔH +2.3 / ΔS +4.1 → (−23, −41);
    // terminal G·C ΔH +0.1 / ΔS −2.8 → (−1, +28); symmetry ΔS −1.4 → +14.
    private const int Primer3TerminalAtDh = -23, Primer3TerminalAtDs = -41;
    private const int Primer3TerminalGcDh = -1, Primer3TerminalGcDs = 28;
    private const int Primer3SymmetryDs = 14;

    // oligotm.c constants: R = 1.987 cal/(K·mol); divalent_to_monovalent factor 120 (mM units).
    private const double Primer3GasConstant = 1.987;
    private const double Primer3DivalentFactor = 120.0;

    /// <summary>
    /// Primer3's default primer melting temperature (°C) — Primer3 <c>seqtm</c> with
    /// PRIMER_TM_FORMULA = SantaLucia 1998 and PRIMER_SALT_CORRECTIONS = SantaLucia 1998:
    /// <para>
    /// [Mon]_eq = [Mon] + 120·√([Mg²⁺] − [dNTP]) (mM; von Ahsen 2001, 0 when [Mg²⁺] ≤ [dNTP]);
    /// ΔS = ΔS°(1 M) + 0.368·(N − 1)·ln([Mon]_eq / 1000);
    /// Tm = ΔH / (ΔS + R·ln(C/4)) − 273.15, with C/1 for a self-complementary oligo,
    /// R = 1.987, C = DNA concentration in mol/L;
    /// for N &gt; <see cref="Primer3MaxNnTmLength"/>: Tm = 81.5 + 16.6·log10([Mon]_eq/1000) + 41·GC/N − 600/N.
    /// </para>
    /// Defaults are Primer3's PRIMER_DNA_CONC = 50 nM, PRIMER_SALT_MONOVALENT = 50 mM,
    /// PRIMER_SALT_DIVALENT = 1.5 mM, PRIMER_DNTP_CONC = 0.6 mM. Matches primer3-py
    /// <c>calc_tm</c> (e.g. AGCTAGCTAGCTAGCTAGCT → 58.101 °C).
    /// </summary>
    /// <param name="primer">Primer sequence (case-insensitive).</param>
    /// <param name="dnaConcentrationNanomolar">Oligo concentration, nM (&gt; 0).</param>
    /// <param name="monovalentMillimolar">Monovalent cation concentration, mM (≥ 0).</param>
    /// <param name="divalentMillimolar">Mg²⁺ concentration, mM (≥ 0).</param>
    /// <param name="dntpMillimolar">dNTP concentration, mM (≥ 0).</param>
    /// <returns>Tm in °C, or <c>double.NaN</c> when the sequence is null/shorter than 2 bases or contains a
    /// non-ACGT character.</returns>
    /// <exception cref="ArgumentOutOfRangeException">For a non-positive DNA concentration, negative ion
    /// concentrations, or a zero total monovalent-equivalent concentration.</exception>
    public static double CalculateMeltingTemperaturePrimer3(
        string primer,
        double dnaConcentrationNanomolar = 50.0,
        double monovalentMillimolar = 50.0,
        double divalentMillimolar = 1.5,
        double dntpMillimolar = 0.6)
    {
        if (!(dnaConcentrationNanomolar > 0) || double.IsInfinity(dnaConcentrationNanomolar))
            throw new ArgumentOutOfRangeException(nameof(dnaConcentrationNanomolar), dnaConcentrationNanomolar,
                "DNA concentration must be a positive finite value in nM.");
        if (!(monovalentMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), monovalentMillimolar,
                "Monovalent cation concentration must be ≥ 0 mM.");
        if (!(divalentMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(divalentMillimolar), divalentMillimolar,
                "Divalent cation concentration must be ≥ 0 mM.");
        if (!(dntpMillimolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(dntpMillimolar), dntpMillimolar,
                "dNTP concentration must be ≥ 0 mM.");

        // divalent_to_monovalent (oligotm.c): no divalent ⇒ dNTP ignored; Mg ≤ dNTP ⇒ no contribution.
        double freeDivalent = divalentMillimolar == 0 ? 0 : Math.Max(0, divalentMillimolar - dntpMillimolar);
        double monovalentEq = monovalentMillimolar + Primer3DivalentFactor * Math.Sqrt(freeDivalent);
        if (!(monovalentEq > 0))
            throw new ArgumentOutOfRangeException(nameof(monovalentMillimolar), monovalentMillimolar,
                "Total monovalent-equivalent cation concentration must be > 0 mM.");

        if (string.IsNullOrEmpty(primer) || primer.Length < 2)
            return double.NaN;

        string seq = primer.ToUpperInvariant();
        int n = seq.Length;
        int gc = 0;
        foreach (char c in seq)
        {
            if (c is 'G' or 'C') gc++;
            else if (c is not ('A' or 'T')) return double.NaN;
        }

        if (n > Primer3MaxNnTmLength)
            return ThermoConstants.CalculateSaltAdjustedTm((double)gc / n, n, monovalentEq / 1000.0);

        // oligotm(): integer accumulation, then ΔH = dh·(−100) cal/mol, ΔS = ds·(−0.1) cal/(K·mol).
        bool symmetric = IsSelfComplementary(seq);
        int dh = 0, ds = symmetric ? Primer3SymmetryDs : 0;
        foreach (char end in new[] { seq[0], seq[^1] })
        {
            if (end is 'A' or 'T') { dh += Primer3TerminalAtDh; ds += Primer3TerminalAtDs; }
            else { dh += Primer3TerminalGcDh; ds += Primer3TerminalGcDs; }
        }
        for (int i = 0; i < n - 1; i++)
        {
            int x = BaseIndex(seq[i]), y = BaseIndex(seq[i + 1]);
            dh += Primer3SantaLucia1998Dh[x, y];
            ds += Primer3SantaLucia1998Ds[x, y];
        }

        double deltaH = dh * -100.0;
        double deltaS = ds * -0.1;
        deltaS += SantaLuciaEntropySaltCoefficient * (n - 1) * Math.Log(monovalentEq / 1000.0);
        // Equation A (self-complementary, C_T/1) or Equation B (C_T/4) of oligotm.c.
        double strandDivisor = symmetric ? 1000000000.0 : 4000000000.0;
        return deltaH / (deltaS + Primer3GasConstant * Math.Log(dnaConcentrationNanomolar / strandDivisor)) - KelvinOffset;

        static int BaseIndex(char c) => c switch { 'A' => 0, 'C' => 1, 'G' => 2, _ => 3 };
    }

    /// <summary>
    /// Calculates GC content as a percentage.
    /// </summary>
    public static double CalculateGcContent(string sequence) =>
        string.IsNullOrEmpty(sequence) ? 0 : sequence.CalculateGcContentFast();

    /// <summary>
    /// Finds the longest homopolymer run (consecutive identical nucleotides, case-insensitive) — the
    /// quantity Primer3 limits with <c>PRIMER_MAX_POLY_X</c> (a primer fails when the run exceeds
    /// it; <c>libprimer3.cc</c> <c>_pr_violates_poly_x</c>). As in Primer3, N is a wildcard that takes
    /// the worst case: a forward scan assigns each N to the preceding non-N base, a reverse scan to
    /// the following one, and the longer run is reported (e.g. ANA → 3, GNGNG → 5, ANGNG → 4).
    /// </summary>
    /// <returns>0 for null/empty input, otherwise ≥ 1.</returns>
    public static int FindLongestHomopolymer(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        string seq = sequence.ToUpperInvariant();
        int len = seq.Length;

        // Forward scan (N counts as the last seen non-N base; leading Ns as the first non-N base).
        char lastNonN = seq[0];
        if (lastNonN == 'N')
        {
            for (int i = 1; i < len; i++)
                if (seq[i] != 'N') { lastNonN = seq[i]; break; }
        }
        int run = 1, maxRun = 1;
        bool hasN = false;
        for (int i = 1; i < len; i++)
        {
            if (seq[i] == 'N') { hasN = true; run++; }
            else if (seq[i] == lastNonN) run++;
            else { run = 1; lastNonN = seq[i]; }
            if (run > maxRun) maxRun = run;
        }
        if (!hasN)
            return maxRun;

        // Reverse scan (N counts as the next non-N base).
        lastNonN = seq[len - 1];
        if (lastNonN == 'N')
        {
            for (int i = len - 2; i >= 0; i--)
                if (seq[i] != 'N') { lastNonN = seq[i]; break; }
        }
        run = 1;
        for (int i = len - 2; i >= 0; i--)
        {
            if (seq[i] == 'N' || seq[i] == lastNonN) run++;
            else { run = 1; lastNonN = seq[i]; }
            if (run > maxRun) maxRun = run;
        }
        return maxRun;
    }

    /// <summary>
    /// Finds the longest dinucleotide repeat (e.g., ATATAT).
    /// </summary>
    public static int FindLongestDinucleotideRepeat(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 4)
            return 0;

        var seq = sequence.ToUpperInvariant();
        int maxRepeats = 0;

        for (int i = 0; i < seq.Length - 3; i++)
        {
            string dinuc = seq.Substring(i, 2);
            int repeats = 1;
            int j = i + 2;

            while (j + 1 < seq.Length && seq.Substring(j, 2) == dinuc)
            {
                repeats++;
                j += 2;
            }

            maxRepeats = Math.Max(maxRepeats, repeats);
        }

        return maxRepeats;
    }

    /// <summary>
    /// Sequence-only stem-loop screen: <c>true</c> when the sequence contains two non-overlapping
    /// segments of <paramref name="minStemLength"/> bases that are exact Watson–Crick reverse
    /// complements of each other (an antiparallel stem) separated by at least
    /// <paramref name="minLoopLength"/> unpaired bases (hairpin loops shorter than 3 nt are sterically
    /// excluded; SantaLucia &amp; Hicks 2004). Case-insensitive; no G·T wobble, mismatches or energies.
    /// This is a structural screen, not a thermodynamic model: the Primer3 hairpin Tm
    /// (PRIMER_HAIRPIN_TH) is <see cref="CalculatePrimer3OligoStructure"/> /
    /// <see cref="CalculateHairpinThermodynamicsNtthal(string, double)"/>, which
    /// <see cref="EvaluatePrimer"/> uses by default. Uses an O(n²) scan below 100 nt and a suffix tree
    /// at ≥ 100 nt (identical results).
    /// </summary>
    /// <param name="sequence">DNA sequence to check.</param>
    /// <param name="minStemLength">Minimum stem length (default 4).</param>
    /// <param name="minLoopLength">Minimum loop length (default 3).</param>
    /// <returns>True if hairpin potential detected.</returns>
    public static bool HasHairpinPotential(string sequence, int minStemLength = 4, int minLoopLength = 3)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < minStemLength * 2 + minLoopLength)
            return false;

        var seq = sequence.ToUpperInvariant();

        // For short sequences (typical primers), use simple O(n²) approach
        // Break-even point is ~100bp based on suffix tree construction overhead
        if (seq.Length < 100)
        {
            return HasHairpinPotentialSimple(seq, minStemLength, minLoopLength);
        }

        // For longer sequences, use suffix tree for O(n) lookup
        return HasHairpinPotentialWithSuffixTree(seq, minStemLength, minLoopLength);
    }

    /// <summary>
    /// Simple O(n²) hairpin detection for short sequences.
    /// </summary>
    private static bool HasHairpinPotentialSimple(string seq, int minStemLength, int minLoopLength)
    {
        // Check for self-complementary regions
        for (int i = 0; i <= seq.Length - minStemLength; i++)
        {
            string fragment = seq.Substring(i, minStemLength);
            // Look for complementary sequence at least minLoopLength positions away
            for (int j = i + minStemLength + minLoopLength; j <= seq.Length - minStemLength; j++)
            {
                string target = seq.Substring(j, minStemLength);
                if (AreComplementary(fragment, Reverse(target)))
                    return true;
            }
        }

        return false;
    }

    /// <summary>
    /// Suffix tree-based O(n) hairpin detection for long sequences.
    /// 
    /// Algorithm: A hairpin forms when a substring S at position i is complementary
    /// to a substring at position j (in reverse). This is equivalent to:
    /// - seq[i..i+k] being present in revComp at some position p
    /// - The positions must satisfy: j = n - p - k, and j >= i + k + minLoopLength
    /// 
    /// We build a suffix tree on seq and search for all substrings of revComp,
    /// checking if any match satisfies the loop constraint.
    /// </summary>
    private static bool HasHairpinPotentialWithSuffixTree(string seq, int minStemLength, int minLoopLength)
    {
        var revComp = DnaSequence.GetReverseComplementString(seq);
        var tree = global::SuffixTree.SuffixTree.Build(seq);

        // For each position in revComp, find matches in seq via suffix tree
        // and check if they form valid hairpin (sufficient loop distance)
        int n = seq.Length;

        // Slide through revComp looking for stems
        for (int p = 0; p <= n - minStemLength; p++)
        {
            var pattern = revComp.AsSpan(p, minStemLength);
            var matches = tree.FindAllOccurrences(pattern);

            foreach (int i in matches)
            {
                // Position in revComp p corresponds to position (n - p - minStemLength) in seq
                // when we reverse complement back
                int j = n - p - minStemLength;

                // Check if positions form valid hairpin: j >= i + minStemLength + minLoopLength
                // Also check i and j don't overlap with the stem itself
                if (j >= i + minStemLength + minLoopLength && j + minStemLength <= n)
                {
                    return true;
                }
                // Also check the reverse case: i is the 3' stem, j is the 5' stem
                if (i >= j + minStemLength + minLoopLength && i + minStemLength <= n)
                {
                    return true;
                }
            }
        }

        return false;
    }

    /// <summary>
    /// Checks whether two primers can form a 3′-end primer-dimer using Primer3's alignment-based
    /// (non-thermodynamic, <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>) pair 3′-complementarity:
    /// returns <c>true</c> when <see cref="CalculatePrimerDimerEndComplementarity"/> ≥
    /// <paramref name="minComplementarity"/>. For ACGT primers the score is an integer, so the
    /// default threshold 4 is exactly Primer3's default <c>PRIMER_PAIR_MAX_COMPL_END = 3.00</c>
    /// (a pair fails when <c>compl_end &gt; 3</c>). The thermodynamic counterpart (Primer3's default
    /// mode) is <see cref="CalculatePrimer3PairComplementarity"/>.
    /// </summary>
    /// <param name="primer1">First primer (5′→3′).</param>
    /// <param name="primer2">Second primer (5′→3′).</param>
    /// <param name="minComplementarity">Minimum 3′-anchored complementarity score that flags a dimer (default 4).</param>
    /// <returns><c>false</c> for null/empty primers.</returns>
    public static bool HasPrimerDimer(string primer1, string primer2, int minComplementarity = 4)
    {
        if (string.IsNullOrEmpty(primer1) || string.IsNullOrEmpty(primer2))
            return false;
        return CalculatePrimerDimerEndComplementarity(primer1, primer2) >= minComplementarity;
    }

    /// <summary>
    /// Primer3 alignment-mode pair 3′-complementarity <c>compl_end</c> (Primer3 <c>libprimer3.cc</c>
    /// <c>characterize_pair</c>, <c>PRIMER_PAIR_COMPL_END</c> with
    /// <c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=0</c>; Rozen &amp; Skaletsky 2000; Untergasser et al. 2012):
    /// the maximum of <c>align(p1, revcomp(p2))</c> and <c>align(p2, revcomp(p1))</c>, where
    /// <c>align</c> is the <c>dpal</c> end-anchored (<c>DPAL_GLOBAL_END</c>) alignment — it must end
    /// at the 3′-terminal base of the first sequence — scored +1 per complementary base pair, −1 per
    /// mismatch, −0.25 against N (any non-ACGT character is scored as N), −2 per single-base gap
    /// (max gap 1), and floored at 0. Two primers whose 3′-terminal k bases are reverse complements
    /// score k; a homopolymer against itself (e.g. A₈/A₈) scores 0 because it cannot pair.
    /// </summary>
    /// <param name="primer1">First primer (5′→3′), case-insensitive.</param>
    /// <param name="primer2">Second primer (5′→3′), case-insensitive.</param>
    /// <returns>The Primer3 <c>compl_end</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerDimerEndComplementarity(string primer1, string primer2)
    {
        if (string.IsNullOrEmpty(primer1) || string.IsNullOrEmpty(primer2))
            return 0;
        string p1 = primer1.ToUpperInvariant();
        string p2 = primer2.ToUpperInvariant();
        // Primer3 compares s1 with s2 taken from the same (top) strand, i.e. the right primer's
        // template-strand copy = revcomp(right primer); then also s2_rev against s1_rev.
        double a = DpalGlobalEndScore(p1, ReverseComplementPrimer3(p2));
        double b = DpalGlobalEndScore(p2, ReverseComplementPrimer3(p1));
        return Math.Max(a, b);
    }

    /// <summary>
    /// Primer3 alignment-mode self 3′-complementarity <c>self_end</c>
    /// (<c>oligo_compl</c>: <c>align(oligo, revcomp(oligo), DPAL_GLOBAL_END)</c>, Primer3
    /// <c>PRIMER_LEFT/RIGHT_SELF_END</c>; default limit <c>PRIMER_MAX_SELF_END = 3.00</c>).
    /// </summary>
    /// <param name="primer">Primer (5′→3′), case-insensitive.</param>
    /// <returns>The Primer3 <c>self_end</c> score (≥ 0); 0 for null/empty input.</returns>
    public static double CalculatePrimerSelfEndComplementarity(string primer)
    {
        if (string.IsNullOrEmpty(primer))
            return 0;
        string p = primer.ToUpperInvariant();
        return DpalGlobalEndScore(p, ReverseComplementPrimer3(p));
    }

    // Primer3 p3_reverse_complement: ACGT complemented, every other character becomes N.
    private static string ReverseComplementPrimer3(string seq)
    {
        var chars = new char[seq.Length];
        for (int i = 0; i < seq.Length; i++)
        {
            chars[seq.Length - 1 - i] = seq[i] switch
            {
                'A' => 'T',
                'T' => 'A',
                'G' => 'C',
                'C' => 'G',
                _ => 'N',
            };
        }
        return new string(chars);
    }

    // dpal default primer-picking scoring (dpal.c set_dpal_args): identity matrix ×100.
    private const int DpalMatch = 100, DpalMismatch = -100, DpalN = -25, DpalGap = -200;

    private static int DpalSsm(char x, char y)
    {
        bool xn = x is not ('A' or 'C' or 'G' or 'T');
        bool yn = y is not ('A' or 'C' or 'G' or 'T');
        if (xn || yn) return DpalN;
        return x == y ? DpalMatch : DpalMismatch;
    }

    /// <summary>
    /// Primer3 <c>align(X, Y, DPAL_GLOBAL_END)</c> score / 100 floored at 0 (libprimer3.cc
    /// <c>align</c>), computed by a line-by-line port of dpal.c
    /// <c>_dpal_long_nopath_maxgap1_global_end</c> (the routine Primer3 runs in DPM_FAST mode with
    /// max_gap = 1); inputs too short for that routine (|X| ≤ 3 or |Y| = 1, where the C code reads
    /// past the sequence end) use dpal.c's <c>_dpal_generic</c> recurrence for GLOBAL_END.
    /// </summary>
    private static double DpalGlobalEndScore(string x, string y)
    {
        int xlen = x.Length, ylen = y.Length;
        // The fast routine reads Y[ylen] (C's NUL terminator) when |X| ≤ 3 and Y[1] when |Y| = 1.
        int smax = xlen < 4 || ylen < 2 ? DpalGlobalEndGeneric(x, y) : DpalGlobalEndFast(x, y);
        return smax < 0 ? 0.0 : smax / 100.0;
    }

    private static int DpalGlobalEndFast(string x, string y)
    {
        int xlen = x.Length, ylen = y.Length;
        const int gap = DpalGap;
        var s0 = new int[xlen];
        var s1 = new int[xlen];
        var s2 = new int[xlen];
        int score, a;

        int smax = DpalSsm(x[xlen - 1], y[0]);
        for (int j = 0; j < xlen; j++) s0[j] = DpalSsm(x[j], y[0]);

        s1[0] = DpalSsm(x[0], y[1]);
        for (int j = 1; j < xlen; j++)
        {
            score = s0[j - 1];
            if (j > 1 && (a = s0[j - 2] + gap) > score) score = a;
            score += DpalSsm(x[j], y[1]);
            if (score > smax && j == xlen - 1) smax = score;
            s1[j] = score;
        }

        int k = ylen - xlen / 2 + 1;
        if (k < 1) k = 1;

        // Rectangular part.
        for (int j = 2; j < k + 1; j++)
        {
            s2[0] = DpalSsm(x[0], y[j]);
            score = s1[0];
            if ((a = s0[0] + gap) > score) score = a;
            score += DpalSsm(x[1], y[j]);
            s2[1] = score;
            for (int i = 2; i < xlen - 1; i++)
            {
                score = s1[i - 2];
                if ((a = s0[i - 1]) > score) score = a;
                score += gap;
                if ((a = s1[i - 1]) > score) score = a;
                score += DpalSsm(x[i], y[j]);
                s2[i] = score;
            }
            score = s1[xlen - 3];
            if ((a = s0[xlen - 2]) > score) score = a;
            score += gap;
            if ((a = s1[xlen - 2]) > score) score = a;
            score += DpalSsm(x[xlen - 1], y[j]);
            s2[xlen - 1] = score;
            if (score > smax) smax = score;
            (s0, s1, s2) = (s1, s2, s0);
        }

        // Triangular part (cells left of the band are not recomputed, exactly as dpal.c).
        int t = 2;
        for (int j = k + 1; j < ylen; j++)
        {
            for (int i = t; i < xlen - 1; i++)
            {
                score = s1[i - 2];
                if ((a = s0[i - 1]) > score) score = a;
                score += gap;
                if ((a = s1[i - 1]) > score) score = a;
                score += DpalSsm(x[i], y[j]);
                s2[i] = score;
            }
            t += 2;
            score = s1[xlen - 3];
            if ((a = s0[xlen - 2]) > score) score = a;
            score += gap;
            if ((a = s1[xlen - 2]) > score) score = a;
            score += DpalSsm(x[xlen - 1], y[j]);
            s2[xlen - 1] = score;
            if (score > smax) smax = score;
            (s0, s1, s2) = (s1, s2, s0);
        }
        return smax;
    }

    // dpal.c _dpal_generic with flag DPAL_GLOBAL_END and max_gap = 1 (score only).
    private static int DpalGlobalEndGeneric(string x, string y)
    {
        int xlen = x.Length, ylen = y.Length;
        var sm = new int[xlen, ylen];
        for (int i = 0; i < xlen; i++) sm[i, 0] = DpalSsm(x[i], y[0]);
        int smax = sm[xlen - 1, 0];
        for (int j = 0; j < ylen; j++) sm[0, j] = DpalSsm(x[0], y[j]);
        for (int i = 1; i < xlen; i++)
        {
            for (int j = 1; j < ylen; j++)
            {
                long a = sm[i - 1, j - 1];
                long b = i > 1 ? (long)sm[i - 2, j - 1] + DpalGap : long.MinValue;
                long c = j > 1 ? (long)sm[i - 1, j - 2] + DpalGap : long.MinValue;
                long best;
                if (a >= b && a >= c) best = a;
                else if (b > a && b >= c) best = b;
                else best = c;
                int score = (int)(best + DpalSsm(x[i], y[j]));
                if (score >= smax && i == xlen - 1) smax = score;
                sm[i, j] = score;
            }
        }
        return smax;
    }

    /// <summary>
    /// 3′-end stability: the duplex ΔG°37 (kcal/mol, 1 M NaCl) of the last five bases (the whole
    /// primer when it is shorter than five), computed exactly as Primer3's
    /// <c>end_oligodg(seq, 5, santalucia)</c> (<c>oligotm.c</c> <c>oligodg</c>; the value Primer3 reports
    /// as <c>PRIMER_{LEFT,RIGHT}_n_END_STABILITY</c> and limits with <c>PRIMER_MAX_END_STABILITY</c>):
    /// SantaLucia (1998) Table 1 unified NN ΔG°37 values summed over the steps, plus initiation
    /// +1.96, +0.05 per terminal A·T base pair and +0.43 for a self-complementary (even-length)
    /// sequence. For a 5-mer this is identical to SantaLucia (1998)'s "initiation with terminal G·C
    /// +0.98 / terminal A·T +1.03" form (0.98+0.98 = 1.96, 1.03−0.98 = 0.05). N is accepted with
    /// Primer3's N-row values (NA 0.58, NC 1.30, NG 1.28, NT 0.88, NN 0.58 negated, etc.) and no A·T
    /// penalty.
    /// <para><b>Sign convention:</b> this method returns the physical ΔG (negative = stable); Primer3's
    /// END_STABILITY is the same quantity with the opposite sign (GCGCG → −6.86 here, 6.86 in
    /// Primer3; TATAT → −0.86 / 0.86).</para>
    /// </summary>
    /// <param name="sequence">Primer sequence (5′→3′), case-insensitive.</param>
    /// <returns>ΔG°37 in kcal/mol; 0 for null/empty input; <c>double.NaN</c> when the 3′ window
    /// contains a character other than A, C, G, T, N (Primer3 <c>OLIGOTM_ERROR</c>).</returns>
    public static double Calculate3PrimeStability(string sequence)
    {
        if (string.IsNullOrEmpty(sequence))
            return 0;

        string window = (sequence.Length > 5 ? sequence[^5..] : sequence).ToUpperInvariant();
        var idx = new int[window.Length];
        for (int i = 0; i < window.Length; i++)
        {
            int b = window[i] switch { 'A' => 0, 'C' => 1, 'G' => 2, 'T' => 3, 'N' => 4, _ => -1 };
            if (b < 0) return double.NaN;
            idx[i] = b;
        }

        // oligodg (santalucia): dg = −1960 [−430 if symmetric] [−50 per terminal A/T] + Σ table, in cal/mol
        // of −ΔG; the method returns ΔG = −dg/1000.
        int dg = -1960;
        if (IsPrimer3Symmetric(window)) dg += -430;
        if (window[0] is 'A' or 'T') dg += -50;
        for (int i = 0; i + 1 < idx.Length; i++)
            dg += Primer3SantaLucia1998Dg[idx[i], idx[i + 1]];
        if (window[^1] is 'A' or 'T') dg += -50;

        return -dg / 1000.0;
    }

    // oligotm.c SantaLucia_1998_dG (−ΔG°37, cal/mol), rows/cols A, C, G, T, N.
    private static readonly int[,] Primer3SantaLucia1998Dg =
    {
        { 1000, 1440, 1280,  880,  880 },
        { 1450, 1840, 2170, 1280, 1450 },
        { 1300, 2240, 1840, 1440, 1300 },
        {  580, 1300, 1450, 1000,  580 },
        {  580, 1300, 1280,  880,  580 },
    };

    // oligotm.c symmetry(): even length and every A/T and C/G position Watson-Crick paired with its mirror.
    private static bool IsPrimer3Symmetric(string seq)
    {
        int n = seq.Length;
        if (n % 2 == 1) return false;
        for (int i = 0; i < n / 2; i++)
        {
            char s = seq[i], e = seq[n - 1 - i];
            if ((s == 'A' && e != 'T') || (s == 'T' && e != 'A') || (e == 'A' && s != 'T') || (e == 'T' && s != 'A'))
                return false;
            if ((s == 'C' && e != 'G') || (s == 'G' && e != 'C') || (e == 'C' && s != 'G') || (e == 'G' && s != 'C'))
                return false;
        }
        return true;
    }

    // ---- Nearest-neighbour salt-corrected Tm (PRIMER-NNTM-001, opt-in) ----------
    // SantaLucia & Hicks (2004) Watson-Crick NN ΔH°/ΔS° (1 M NaCl; Biopython DNA_NN4 — NOT the
    // SantaLucia 1998 / Allawi & SantaLucia 1997 set: AA/TT is −7.6/−21.3 here, −7.9/−22.2 in 1998,
    // and 2004 uses one duplex-initiation term plus a terminal A·T penalty instead of per-end
    // initiation) with the bimolecular Tm equation and published salt corrections. This is an OPT-IN
    // design Tm: the default CalculateMeltingTemperature (Wallace / Marmur-Doty) is unchanged.
    // All arithmetic is delegated to the canonical NN core ThermoConstants.CalculateNearestNeighborDuplex
    // (a line-by-line port of Biopython Tm_NN); this class only maps its API onto it:
    //   parameter set = SantaLuciaHicks2004 (DNA_NN4); R = 1.9872 (SantaLucia & Hicks 2004 Eq. 3;
    //   Biopython uses 1.987); C_T total → dnac1 = dnac2 = C_T/2 (k = C_T/4), self-complementary →
    //   dnac1 = C_T (k = C_T); self-complementarity detected from the sequence; salt mode → Biopython
    //   salt_correction method 0 / 5 / 6 / 7 with [Na⁺] (and, for method 7 only, [Mg²⁺], [dNTPs]).
    // Sources: SantaLucia J, Hicks D (2004) Annu Rev Biophys Biomol Struct 33:415, Table 1 + Eq. 3/5;
    //   Owczarzy R et al. (2004) Biochemistry 43:3537 (monovalent); Owczarzy R et al. (2008)
    //   Biochemistry 47:5336 (Mg²⁺/dNTP); Biopython 1.88 Bio.SeqUtils.MeltingTemp (reference).

    /// <summary>R of SantaLucia &amp; Hicks (2004) Eq. 3, used by the NN Tm methods of this class.</summary>
    private const double GasConstant = ThermoConstants.NnGasConstantSantaLuciaHicks2004;

    // Strand-concentration divisor x in Tm = ΔH°/(ΔS° + R·ln(C_T/x)) (SantaLucia & Hicks 2004 Eq. 3):
    // x = 4 for non-self-complementary, x = 1 for self-complementary duplexes (used by the dimer Tm).
    private const double NonSelfComplementaryFactor = 4.0;
    private const double SelfComplementaryFactor = 1.0;

    // Kelvin-to-Celsius offset.
    private const double KelvinOffset = 273.15;

    // Default total strand concentration C_T = 0.5 µM (a common PCR primer working
    // concentration). Exposed as a parameter; the caller may override.
    private const double DefaultStrandConcentrationMolar = 0.5e-6;

    // Default total strand concentration for the intermolecular dimer Tm. Primer3 / ntthal
    // uses dna_conc = 50 nM (thal.c default a->dna_conc = 50, lines 829/844); the opt-in dimer
    // methods adopt the same convention so they reproduce the ntthal reference out of the box.
    private const double DefaultDimerStrandConcentrationMolar = 50e-9;

    // The NN parameter set of this class's NN Tm, hairpin and dimer helpers (SantaLucia & Hicks 2004).
    private const NnParameterSet DesignNnParameterSet = NnParameterSet.SantaLuciaHicks2004;

    // SantaLucia & Hicks (2004) Table 1 initiation-type terms of DesignNnParameterSet
    // (duplex initiation +0.2/−5.7, terminal A·T +2.2/+6.9, symmetry 0/−1.4).
    private static readonly NnInitiationTerms DesignNnInitiation =
        ThermoConstants.GetNearestNeighborInitiation(DesignNnParameterSet);

    // Watson–Crick stack of the top-strand dinucleotide (upper-case ACGT) from DesignNnParameterSet.
    private static bool TryGetDesignStack(string dinucleotide, out (double DeltaH, double DeltaS) p) =>
        ThermoConstants.TryGetNearestNeighborStack(DesignNnParameterSet, dinucleotide, out p);

    /// <summary>Salt-correction mode for <see cref="CalculateMeltingTemperatureNN"/>.</summary>
    public enum SaltCorrectionMode
    {
        /// <summary>No correction — Tm at the SantaLucia 1 M NaCl reference state (Biopython saltcorr 0).</summary>
        None,

        /// <summary>
        /// SantaLucia (1998) / SantaLucia &amp; Hicks (2004) Eq. 5 entropy correction (Biopython saltcorr 5):
        /// ΔS°[Na] = ΔS°[1 M] + 0.368·(N − 1)·ln[Na⁺], N = oligo length (N − 1 = half the duplex phosphates);
        /// applied to ΔS° before the Tm equation.
        /// </summary>
        SantaLuciaEntropy,

        /// <summary>
        /// Owczarzy et al. (2004) monovalent 1/Tm correction (Biochemistry 43:3537; Biopython saltcorr 6):
        /// 1/Tm = 1/Tm(1 M) + (4.29·f(GC) − 3.95)·10⁻⁵·ln[Na⁺] + 9.40·10⁻⁶·ln²[Na⁺]. [Mg²⁺] and dNTPs are ignored.
        /// </summary>
        Owczarzy2004Monovalent,

        /// <summary>
        /// Owczarzy et al. (2008) divalent Mg²⁺ (and dNTP-adjusted) correction (Biochemistry 47:5336; Biopython
        /// saltcorr 7); reduces to the 2004 monovalent form when √[Mg²⁺]/[Na⁺] &lt; 0.22.
        /// </summary>
        Owczarzy2008Divalent
    }

    // SantaLucia (1998) / SantaLucia & Hicks (2004) Eq. 5 entropy salt-correction coefficient (0.368),
    // used by the Primer3 seqtm port and the dimer Tm.
    private const double SantaLuciaEntropySaltCoefficient = 0.368;

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M NaCl) of a DNA oligonucleotide and its
    /// perfect complement with the SantaLucia &amp; Hicks (2004) nearest-neighbour parameters (Biopython
    /// <c>DNA_NN4</c>): duplex initiation, a terminal A·T penalty per A·T-closed end, the stacks and (for a
    /// self-complementary sequence) the symmetry correction. Equals Biopython
    /// <c>Tm_NN(seq, nn_table=DNA_NN4, selfcomp=…)</c>'s ΔH°/ΔS°; computed by
    /// <see cref="ThermoConstants.CalculateNearestNeighborThermodynamics"/>.
    /// </summary>
    /// <param name="sequence">DNA sequence (one strand, 5'→3'; case-insensitive).</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> if the sequence is empty,
    /// shorter than 2 bases, or contains a non-ACGT character.</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamics(string sequence)
    {
        if (string.IsNullOrEmpty(sequence) || sequence.Length < 2)
            return null;

        string seq = sequence.ToUpperInvariant();
        if (!IsAcgtOnly(seq))
            return null; // non-ACGT base present

        bool selfComp = IsSelfComplementary(seq);
        var (dH, dS) = ThermoConstants.CalculateNearestNeighborThermodynamics(
            seq, parameterSet: DesignNnParameterSet, selfComplementary: selfComp, check: false);
        return (dH, dS, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) of a primer/oligonucleotide using the
    /// SantaLucia &amp; Hicks (2004) nearest-neighbour thermodynamics (Biopython <c>DNA_NN4</c>) and the
    /// bimolecular Tm equation, with an optional published salt correction. <b>Opt-in</b>: the default
    /// <see cref="CalculateMeltingTemperature(string)"/> (Wallace / Marmur-Doty) is unchanged.
    /// <para>
    /// Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15, with R = 1.9872 cal/(K·mol),
    /// x = 4 for a non-self-complementary duplex and x = 1 for a self-complementary one
    /// (SantaLucia &amp; Hicks 2004, Eq. 3); salt corrections per <paramref name="saltMode"/>.
    /// Identical to Biopython <c>Tm_NN(seq, nn_table=DNA_NN4, dnac1=dnac2=C_T/2 (self-complementary:
    /// dnac1=C_T), selfcomp, Na, Mg, dNTPs, saltcorr=0/5/6/7)</c> with Biopython's R = 1.987 replaced by
    /// 1.9872 (≈ +0.005 °C); computed by the canonical <see cref="ThermoConstants.CalculateNearestNeighborDuplex"/>.
    /// </para>
    /// </summary>
    /// <param name="primer">DNA primer sequence (5'→3'). Must be ≥ 2 ACGT bases.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L
    /// (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation ([Na⁺]+[K⁺]+[Tris]/2) concentration in
    /// mol/L (default 0.05 M = 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the
    /// <see cref="SaltCorrectionMode.Owczarzy2008Divalent"/> mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0; sequesters Mg²⁺
    /// in the divalent mode).</param>
    /// <param name="saltMode">Which salt correction to apply (default
    /// <see cref="SaltCorrectionMode.Owczarzy2004Monovalent"/>).</param>
    /// <returns>The nearest-neighbour Tm in °C, or <c>double.NaN</c> if the sequence is
    /// empty, shorter than 2 bases, or contains a non-ACGT character.</returns>
    public static double CalculateMeltingTemperatureNN(
        string primer,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        ValidateNnConditions(strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar);

        var thermo = CalculateNearestNeighborThermodynamics(primer);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        return NnTm(dH, dS, primer.ToUpperInvariant(), selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);
    }

    // Parameter-domain guards: the Tm equation takes R·ln(C_T/x) and the salt corrections take
    // ln([Na⁺]) / ln([Mg²⁺]); a zero or negative concentration would make ln undefined. Reject these
    // out-of-domain inputs explicitly (Biopython raises ValueError for them too).
    private static void ValidateNnConditions(
        double strandConcentrationMolar, double sodiumMolar, double magnesiumMolar, double dntpMolar)
    {
        if (!(strandConcentrationMolar > 0))
            throw new ArgumentOutOfRangeException(nameof(strandConcentrationMolar),
                strandConcentrationMolar, "Strand concentration C_T must be > 0 mol/L.");
        if (!(sodiumMolar > 0))
            throw new ArgumentOutOfRangeException(nameof(sodiumMolar),
                sodiumMolar, "Monovalent cation concentration [Na⁺] must be > 0 mol/L.");
        if (!(magnesiumMolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(magnesiumMolar),
                magnesiumMolar, "[Mg²⁺] must be ≥ 0 mol/L.");
        if (!(dntpMolar >= 0))
            throw new ArgumentOutOfRangeException(nameof(dntpMolar),
                dntpMolar, "Total dNTP concentration must be ≥ 0 mol/L.");
    }

    // Maps this class's NN Tm API onto the canonical Biopython Tm_NN core (see the region comment):
    // C_T → dnac1 = dnac2 = C_T/2 (k = C_T/4) or dnac1 = C_T for a self-complementary duplex (k = C_T);
    // mol/L → mM; salt mode → method 0/5/6/7 ([Mg²⁺]/[dNTPs] only for method 7); R = 1.9872.
    private static double NnTm(
        double deltaH, double deltaS, string saltSequence, bool selfComp,
        double strandConcentrationMolar, double sodiumMolar, double magnesiumMolar, double dntpMolar,
        SaltCorrectionMode saltMode)
    {
        double dnac = selfComp ? strandConcentrationMolar * 1e9 : strandConcentrationMolar * 1e9 / 2.0;
        var method = saltMode switch
        {
            SaltCorrectionMode.None => NnSaltCorrection.None,
            SaltCorrectionMode.SantaLuciaEntropy => NnSaltCorrection.SantaLucia1998Entropy,
            SaltCorrectionMode.Owczarzy2004Monovalent => NnSaltCorrection.Owczarzy2004,
            SaltCorrectionMode.Owczarzy2008Divalent => NnSaltCorrection.Owczarzy2008,
            _ => throw new ArgumentOutOfRangeException(nameof(saltMode), saltMode, "Unknown salt-correction mode.")
        };
        bool divalent = method == NnSaltCorrection.Owczarzy2008;
        try
        {
            return ThermoConstants.CalculateNearestNeighborTmFromThermodynamics(
                deltaH, deltaS, saltSequence, dnac1: dnac, dnac2: dnac, selfComplementary: selfComp,
                sodium: sodiumMolar * 1000.0,
                magnesium: divalent ? magnesiumMolar * 1000.0 : 0,
                dntps: divalent ? dntpMolar * 1000.0 : 0,
                saltCorrection: method, gasConstant: GasConstant).MeltingTemperature;
        }
        catch (ArgumentException)
        {
            return double.NaN; // degenerate duplex the reference rejects (Tm_NN ZeroDivisionError)
        }
    }

    // ---- NN internal-mismatch + dangling-end Tm (PRIMER-NNTM-001, opt-in extension) -------
    // Extends the perfect-match NN model to a probe–target duplex with internal mismatches,
    // terminal mismatches and/or a single unpaired dangling base at either end, exactly as
    // Biopython Tm_NN(seq, c_seq, shift, nn_table=DNA_NN4, tmm_table=DNA_TMM1, imm_table=DNA_IMM1,
    // de_table=DNA_DE1) does (the core ThermoConstants.CalculateNearestNeighborDuplex is a port of it).
    //
    // Convention of this API: the top strand is 5'→3'; the bottom strand is supplied 3'→5' (the
    // complement of the top read in the SAME left-to-right order, NOT the reverse complement), so
    // column i pairs top[i] with bottom[i]. A '.' as the first or last character of a strand marks
    // the missing partner of a single dangling base on the other strand. It maps onto Tm_NN as
    // seq = top without '.', c_seq = bottom without '.', shift = (leading '.' of top) − (leading '.'
    // of bottom). Tm_NN then scores, in this order: dangling ends (Bommarito et al. 2000), terminal
    // mismatches (SantaLucia & Peyret 2001, Biopython DNA_TMM1), the initiation terms — the terminal
    // A·T penalty is taken from the first and last base of the TOP strand without '.' (Tm_NN's
    // `ends = seq[0] + seq[-1]`) — and every remaining stack from the internal-mismatch table
    // (Allawi & SantaLucia 1997/1998, Peyret et al. 1999) or the Watson-Crick table.

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol), 1 M NaCl) of a probe–target DNA duplex
    /// that may contain internal single mismatches (Allawi &amp; SantaLucia 1997/1998; Peyret et al. 1999),
    /// terminal mismatches (SantaLucia &amp; Peyret 2001) and/or a single dangling end at either end
    /// (Bommarito et al. 2000), with the SantaLucia &amp; Hicks (2004) Watson–Crick parameters. Identical to
    /// Biopython <c>Tm_NN(seq, c_seq, shift, nn_table=DNA_NN4)</c>'s ΔH°/ΔS° (see the mapping in the
    /// remarks of <see cref="CalculateMeltingTemperatureNNMismatch"/>).
    /// </summary>
    /// <param name="topStrand">Top strand 5'→3'. May start/end with a single '.' marking a
    /// dangling base of the bottom strand.</param>
    /// <param name="bottomStrand">Bottom strand written 3'→5' (the complement of the top
    /// read left-to-right, NOT the reverse complement), so base i pairs with top base i.
    /// May start/end with a single '.' marking a dangling base of the top strand.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) or <c>null</c> if the strands are null,
    /// unequal length, shorter than two columns, contain a character other than A/C/G/T (case-insensitive)
    /// or a '.' that is not a single terminal marker facing a base, or contain a stack with no parameter
    /// (e.g. adjacent mismatches other than the tandem G·T motifs).</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsMismatch(
        string topStrand, string bottomStrand)
    {
        var mapped = MapMismatchDuplex(topStrand, bottomStrand);
        if (mapped is null)
            return null;

        var (seq, cSeq, shift, selfComp) = mapped.Value;
        try
        {
            var (dH, dS) = ThermoConstants.CalculateNearestNeighborThermodynamics(
                seq, cSeq, shift, DesignNnParameterSet, selfComp, check: false, strict: true);
            return (dH, dS, selfComp);
        }
        catch (ArgumentException)
        {
            return null; // a neighbour pair with no thermodynamic parameter
        }
    }

    // Validates the column-aligned (top, bottom) pair and maps it to Tm_NN's (seq, c_seq, shift).
    private static (string Seq, string CSeq, int Shift, bool SelfComp)? MapMismatchDuplex(string top, string bottom)
    {
        if (top is null || bottom is null)
            return null;
        string t = top.ToUpperInvariant();
        string b = bottom.ToUpperInvariant();
        int n = t.Length;
        if (n != b.Length || n < 2)
            return null;

        for (int i = 0; i < n; i++)
        {
            char x = t[i], y = b[i];
            bool xDot = x == '.', yDot = y == '.';
            if ((xDot || yDot) && i != 0 && i != n - 1)
                return null; // a dangling marker must be terminal
            if (xDot && yDot)
                return null; // a column with no base at all
            if ((!xDot && x is not ('A' or 'C' or 'G' or 'T')) || (!yDot && y is not ('A' or 'C' or 'G' or 'T')))
                return null;
        }

        string seq = t.Replace(".", string.Empty, StringComparison.Ordinal);
        string cSeq = b.Replace(".", string.Empty, StringComparison.Ordinal);
        if (seq.Length < 1 || cSeq.Length < 1)
            return null;
        int shift = (t[0] == '.' ? 1 : 0) - (b[0] == '.' ? 1 : 0);

        // Symmetry term only for a fully paired, self-complementary duplex (Tm_NN's selfcomp flag).
        bool hasDangling = t.Contains('.') || b.Contains('.');
        bool selfComp = !hasDangling && IsSelfComplementary(t)
                        && string.Equals(b, Complement(t), StringComparison.Ordinal);
        return (seq, cSeq, shift, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) for a probe–target DNA duplex that may contain
    /// internal and terminal mismatches and/or a single dangling end at either end, with the
    /// SantaLucia &amp; Hicks (2004) Watson–Crick parameters plus the internal-mismatch, terminal-mismatch and
    /// Bommarito (2000) dangling-end terms (<see cref="CalculateNearestNeighborThermodynamicsMismatch"/>), the
    /// same bimolecular Tm equation and salt corrections as <see cref="CalculateMeltingTemperatureNN"/>.
    /// <b>Opt-in extension</b>: a fully paired duplex through this path equals
    /// <see cref="CalculateMeltingTemperatureNN"/>.
    /// </summary>
    /// <remarks>
    /// Equals Biopython <c>Tm_NN(seq, c_seq=…, shift=…, nn_table=DNA_NN4, dnac1=dnac2=C_T/2, Na, Mg, dNTPs,
    /// saltcorr=0/5/6/7)</c> with R = 1.9872, where seq/c_seq are the strands without '.', shift = (leading '.'
    /// of the top) − (leading '.' of the bottom). The salt correction uses the top strand without '.'
    /// (its length and GC fraction), as Tm_NN does.
    /// </remarks>
    /// <param name="topStrand">Top strand 5'→3' (may carry a leading/trailing '.' dangling-end marker).</param>
    /// <param name="bottomStrand">Bottom strand 3'→5', aligned base-for-base under the top
    /// (complement direction, NOT reverse complement; may carry a '.' dangling-end marker).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction to apply (default Owczarzy2004Monovalent).</param>
    /// <returns>The NN Tm in °C, or <c>double.NaN</c> if the duplex is not computable
    /// (see <see cref="CalculateNearestNeighborThermodynamicsMismatch"/>).</returns>
    /// <exception cref="ArgumentOutOfRangeException">A non-positive C_T or [Na⁺], or a negative [Mg²⁺]/[dNTPs].</exception>
    public static double CalculateMeltingTemperatureNNMismatch(
        string topStrand,
        string bottomStrand,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        ValidateNnConditions(strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar);

        var thermo = CalculateNearestNeighborThermodynamicsMismatch(topStrand, bottomStrand);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        string topPaired = topStrand.ToUpperInvariant().Replace(".", string.Empty, StringComparison.Ordinal);
        return NnTm(dH, dS, topPaired, selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);
    }

    // ---- LNA (locked nucleic acid)-adjusted NN Tm (PROBE-DESIGN-001, opt-in extension) ----
    // Extends the perfect-match DNA NN model with the McTigue, Peterson & Kahn (2004) LNA-DNA
    // nearest-neighbour increments so the NN Tm can be computed for a DNA oligo carrying one or
    // more INTERNAL LNA substitutions on one strand. The LNA value is an ADDITIVE increment
    // (ΔΔH°, ΔΔS°) added to the underlying DNA NN stack for each step containing the LNA base,
    // exactly as the MELTING 5 reference implementation realises it (McTigue04LockedAcid.java:
    // DNA NN sum, then `enthalpy += lockedAcidValue.getEnthalpy()`). The perfect-match
    // CalculateMeltingTemperatureNN above is UNCHANGED; this is opt-in.
    //
    // Convention: the increment for a step is keyed by the DNA dinucleotide step (e.g. "TT")
    // and which of the two bases of that step is locked (0 = the 5' base, 1 = the 3' base),
    // matching the paper's MX_L / X_L N notation (XML keys e.g. "TTL/AA" = step TT, 3' base
    // locked; "TLG/AC" = step TG, 5' base locked). Terminal LNA positions are NOT parameterised
    // by McTigue (2004) and are rejected (return not-computable), per MELTING isApplicable.
    //
    // Source (retrieved & cross-checked this session, 2026-06-24):
    //   McTigue PM, Peterson RJ, Kahn JD (2004) Biochemistry 43:5388-5405, DOI 10.1021/bi035976d
    //   — ΔΔH°/ΔΔS° for all 32 LNA+DNA:DNA nearest neighbours.
    //   Parameter values transcribed VERBATIM from the MELTING 5 data file
    //   "McTigue2004lockedmn.xml" (Dumousseau et al. 2012, BMC Bioinformatics 13:101; mirrored in
    //   aravind-j/rmelting). MELTING stores ΔΔH°/ΔΔS° in cal/mol and cal/(mol·K); the kcal/mol
    //   values below are XML_value / 1000. Worked example reproduced: CCATT(L)GCTACC at C=1e-4,
    //   Na=1 → ΔH°=-80.014 kcal/mol, ΔS°=-216.6 cal/(mol·K), Tm=63.528 °C (MELTING mct04: 63.614).

    /// <summary>Which base of a nearest-neighbour dinucleotide step is the LNA monomer.</summary>
    private enum LnaStepPosition
    {
        /// <summary>The 5' (first) base of the step is locked — McTigue X_L N (key e.g. "TLG/AC").</summary>
        FivePrime = 0,

        /// <summary>The 3' (second) base of the step is locked — McTigue MX_L (key e.g. "TTL/AA").</summary>
        ThreePrime = 1
    }

    /// <summary>
    /// McTigue, Peterson &amp; Kahn (2004) LNA-DNA nearest-neighbour increments
    /// (ΔΔH° in kcal/mol, ΔΔS° in cal/(K·mol)), added to the base DNA NN stack for the step
    /// containing the LNA base. Key = (DNA dinucleotide step 5'→3', which base is locked).
    /// All 32 nearest neighbours are present (16 with the 5' base locked, 16 with the 3' base
    /// locked). Values transcribed verbatim from MELTING 5 <c>McTigue2004lockedmn.xml</c>
    /// (cal/mol ÷ 1000 = kcal/mol); the XML <c>sequence</c> key is shown in the comment.
    /// Source: McTigue et al. (2004) Biochemistry 43:5388-5405 (DOI 10.1021/bi035976d).
    /// </summary>
    private static readonly Dictionary<(string Step, LnaStepPosition Locked), (double DeltaH, double DeltaS)> McTigueLnaIncrements = new()
    {
        // 5'-base locked (X_L N): XML key "XLY/comp".
        [("AA", LnaStepPosition.FivePrime)] = (0.707, 2.5),    // ALA/TT
        [("AT", LnaStepPosition.FivePrime)] = (2.282, 7.5),    // ALT/TA
        [("AG", LnaStepPosition.FivePrime)] = (0.264, 2.6),    // ALG/TC
        [("AC", LnaStepPosition.FivePrime)] = (1.131, 4.1),    // ALC/TG
        [("TA", LnaStepPosition.FivePrime)] = (-0.046, 1.6),   // TLA/AT
        [("TT", LnaStepPosition.FivePrime)] = (1.528, 5.3),    // TLT/AA
        [("TG", LnaStepPosition.FivePrime)] = (-1.540, -3.0),  // TLG/AC
        [("TC", LnaStepPosition.FivePrime)] = (1.893, 6.7),    // TLC/AG
        [("GA", LnaStepPosition.FivePrime)] = (3.162, 10.5),   // GLA/CT
        [("GT", LnaStepPosition.FivePrime)] = (-0.212, 0.1),   // GLT/CA
        [("GG", LnaStepPosition.FivePrime)] = (-2.844, -6.7),  // GLG/CC
        [("GC", LnaStepPosition.FivePrime)] = (-0.360, -0.3),  // GLC/CG
        [("CA", LnaStepPosition.FivePrime)] = (1.049, 4.3),    // CLA/GT
        [("CT", LnaStepPosition.FivePrime)] = (0.708, 4.2),    // CLT/GA
        [("CG", LnaStepPosition.FivePrime)] = (0.785, 3.7),    // CLG/GC
        [("CC", LnaStepPosition.FivePrime)] = (2.096, 8.0),    // CLC/GG

        // 3'-base locked (MX_L): XML key "MXL/comp".
        [("AA", LnaStepPosition.ThreePrime)] = (0.992, 4.1),   // AAL/TT
        [("AT", LnaStepPosition.ThreePrime)] = (1.816, 6.9),   // ATL/TA
        [("AG", LnaStepPosition.ThreePrime)] = (-1.200, -1.8), // AGL/TC
        [("AC", LnaStepPosition.ThreePrime)] = (2.890, 10.6),  // ACL/TG
        [("TA", LnaStepPosition.ThreePrime)] = (1.591, 5.3),   // TAL/AT
        [("TT", LnaStepPosition.ThreePrime)] = (2.326, 8.1),   // TTL/AA
        [("TG", LnaStepPosition.ThreePrime)] = (2.165, 7.2),   // TGL/AC
        [("TC", LnaStepPosition.ThreePrime)] = (0.609, 3.2),   // TCL/AG
        [("GA", LnaStepPosition.ThreePrime)] = (0.444, 2.9),   // GAL/CT
        [("GT", LnaStepPosition.ThreePrime)] = (-0.635, -0.3), // GTL/CA
        [("GG", LnaStepPosition.ThreePrime)] = (-0.943, -0.9), // GGL/CC
        [("GC", LnaStepPosition.ThreePrime)] = (-0.925, -1.1), // GCL/CG
        [("CA", LnaStepPosition.ThreePrime)] = (1.358, 4.4),   // CAL/GT
        [("CT", LnaStepPosition.ThreePrime)] = (-1.671, -4.1), // CTL/GA
        [("CG", LnaStepPosition.ThreePrime)] = (-0.276, -0.7), // CGL/GC
        [("CC", LnaStepPosition.ThreePrime)] = (2.063, 7.6)    // CCL/GG
    };

    /// <summary>
    /// Computes the duplex ΔH° (kcal/mol) and ΔS° (cal/(K·mol)) of a DNA oligonucleotide that
    /// carries one or more <b>internal</b> LNA (locked nucleic acid) substitutions, by adding the
    /// McTigue, Peterson &amp; Kahn (2004) LNA-DNA nearest-neighbour increments to the SantaLucia
    /// &amp; Hicks (2004) DNA NN stack (initiation + terminal-A·T + symmetry are computed on the underlying
    /// DNA sequence, unchanged). <b>Opt-in</b>: the perfect-match
    /// <see cref="CalculateNearestNeighborThermodynamics"/> is unchanged, and an empty
    /// <paramref name="lnaPositions"/> reproduces it exactly.
    /// </summary>
    /// <param name="sequence">DNA sequence (one strand, 5'→3'); the LNA monomers are at the given
    /// positions of this sequence. Must be ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based positions of the LNA monomers within
    /// <paramref name="sequence"/>. Order and duplicates are tolerated. A <b>terminal</b> position
    /// (0 or length−1) is not parameterised by McTigue (2004) and makes the result not computable.</param>
    /// <returns>(ΔH°, ΔS°, IsSelfComplementary) of the LNA-substituted duplex, or <c>null</c> if the
    /// sequence is empty/&lt; 2 bases/contains a non-ACGT base, or any LNA position is out of range
    /// or terminal.</returns>
    public static (double DeltaH, double DeltaS, bool IsSelfComplementary)? CalculateNearestNeighborThermodynamicsLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions)
    {
        ArgumentNullException.ThrowIfNull(lnaPositions);

        var dna = CalculateNearestNeighborThermodynamics(sequence);
        if (dna is null)
            return null;

        string seq = sequence.ToUpperInvariant();
        var locked = new HashSet<int>();
        foreach (int pos in lnaPositions)
        {
            // McTigue (2004) parameters are for internal LNA only — reject terminal/out-of-range.
            if (pos <= 0 || pos >= seq.Length - 1)
                return null;
            locked.Add(pos);
        }

        var (dH, dS, selfComp) = dna.Value;

        // Add the McTigue increment to each NN step (i, i+1) that contains an LNA base.
        for (int i = 0; i < seq.Length - 1; i++)
        {
            string step = seq.Substring(i, 2);
            if (locked.Contains(i)
                && McTigueLnaIncrements.TryGetValue((step, LnaStepPosition.FivePrime), out var inc5))
            {
                dH += inc5.DeltaH; dS += inc5.DeltaS;
            }
            if (locked.Contains(i + 1)
                && McTigueLnaIncrements.TryGetValue((step, LnaStepPosition.ThreePrime), out var inc3))
            {
                dH += inc3.DeltaH; dS += inc3.DeltaS;
            }
        }

        return (dH, dS, selfComp);
    }

    /// <summary>
    /// Computes the design melting temperature (°C) of a DNA oligonucleotide carrying one or more
    /// <b>internal</b> LNA substitutions, using the McTigue (2004) LNA-DNA nearest-neighbour
    /// increments on top of the SantaLucia &amp; Hicks (2004) DNA NN model, with the same bimolecular Tm
    /// equation and optional salt corrections as <see cref="CalculateMeltingTemperatureNN"/>.
    /// <b>Opt-in</b>: the perfect-match <see cref="CalculateMeltingTemperatureNN"/> is unchanged,
    /// and an empty <paramref name="lnaPositions"/> equals it exactly.
    /// </summary>
    /// <param name="sequence">DNA sequence (5'→3'). Must be ≥ 2 ACGT bases.</param>
    /// <param name="lnaPositions">Zero-based positions of the internal LNA monomers (see
    /// <see cref="CalculateNearestNeighborThermodynamicsLna"/>).</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L (default 0.5 µM).</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM).</param>
    /// <param name="magnesiumMolar">[Mg²⁺] in mol/L (default 0; only used by the divalent mode).</param>
    /// <param name="dntpMolar">Total dNTP concentration in mol/L (default 0).</param>
    /// <param name="saltMode">Salt correction to apply (default Owczarzy2004Monovalent).</param>
    /// <returns>The LNA-adjusted NN Tm in °C, or <c>double.NaN</c> if the duplex is not computable
    /// (empty/&lt; 2 bases/non-ACGT, or an out-of-range/terminal LNA position).</returns>
    public static double CalculateMeltingTemperatureNNLna(
        string sequence,
        IReadOnlyCollection<int> lnaPositions,
        double strandConcentrationMolar = DefaultStrandConcentrationMolar,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double magnesiumMolar = 0.0,
        double dntpMolar = 0.0,
        SaltCorrectionMode saltMode = SaltCorrectionMode.Owczarzy2004Monovalent)
    {
        var thermo = CalculateNearestNeighborThermodynamicsLna(sequence, lnaPositions);
        if (thermo is null)
            return double.NaN;

        var (dH, dS, selfComp) = thermo.Value;
        return NnTm(dH, dS, sequence.ToUpperInvariant(), selfComp,
            strandConcentrationMolar, sodiumMolar, magnesiumMolar, dntpMolar, saltMode);
    }

    // ---- DNA hairpin folding + secondary-structure (hairpin) Tm (PRIMER-TM-001, opt-in) ----
    // Finds the most stable intramolecular hairpin (a single stem closing one hairpin loop)
    // of a DNA oligo and computes its ΔH°/ΔS°/ΔG°37 and unimolecular melting temperature.
    // The perfect-match CalculateMeltingTemperatureNN and the default Wallace/Marmur-Doty Tm
    // are UNCHANGED; this is a new opt-in capability.
    //
    // Model (SantaLucia & Hicks 2004, Annu Rev Biophys 33:415, "Hairpin Loops", Eqs 8–11):
    //   ΔG°37(hairpin) = Σ stem NN stacks (Table 1)  +  ΔG°37(hairpin loop of N) (Table 4)
    //   "To compute the stability of a complete hairpin + stem, one simply adds the
    //    salt-corrected base pair NN contributions (Table 1; Equation 3) to the loop energy
    //    from Equations 8–10." (paper, p.428).
    //   The stem contributes only its nearest-neighbour STACKS (Table 1 propagation terms) —
    //   the bimolecular duplex-initiation term (+0.2/−5.7) is a TWO-strand nucleation cost and
    //   does NOT apply to a unimolecular hairpin (the loop-initiation term is the nucleation
    //   cost instead). Loop ΔH° = 0 for every loop size, and the loop ΔS° increment is
    //   ΔS° = −ΔG°37 × 1000 / 310.15 (Table 4 note: "ΔS° = ΔG°37 × 1000/310.15"; the loop is
    //   destabilising so ΔG°37 > 0 → ΔS° < 0).
    //   Two-state hairpin Tm (Eq 11) is UNIMOLECULAR/concentration-independent:
    //       Tm = ΔH° × 1000 / ΔS° − 273.15        (NO R·ln(C_T/x) strand-concentration term).
    //
    // NOT bundled (honest residual): the supplementary triloop/tetraloop bonus tables (length-3
    // and length-4 special loops) and the terminal-mismatch increment (the first mismatch stack
    // closing loops of length ≥4) are separate Annual-Reviews supplementary tables not embedded
    // here. They are exposed as an OPT-IN caller-supplied additive ΔG°37/ΔH° adjustment
    // (default 0) so a caller who has those tables can supply the increment; without it the
    // result is the stem-stack + loop-initiation core, which is exact and fully sourced.

    /// <summary>Hairpin loop ΔG°37 increment (kcal/mol, 1 M NaCl) by loop size (number of
    /// unpaired loop nucleotides). SantaLucia &amp; Hicks (2004) Table 4 "Hairpin loops" column;
    /// sizes 3–30 are tabulated. ΔH° = 0 for all sizes; ΔS° = −ΔG°37·1000/310.15.</summary>
    private static readonly Dictionary<int, double> HairpinLoopInitiationDeltaG = new()
    {
        [3] = 3.5, [4] = 3.5, [5] = 3.3, [6] = 4.0, [7] = 4.2, [8] = 4.3, [9] = 4.5,
        [10] = 4.6, [12] = 5.0, [14] = 5.1, [16] = 5.3, [18] = 5.5, [20] = 5.7,
        [25] = 6.1, [30] = 6.3
    };

    // SantaLucia & Hicks (2004): minimum sterically allowed hairpin loop size is 3 nt
    // ("Hairpin loops with lengths shorter than 3 are sterically prohibited.").
    private const int MinHairpinLoopSize = 3;

    // Reference temperature for the ΔG°37 ↔ ΔS° conversion (310.15 K = 37 °C).
    // SantaLucia & Hicks (2004) Table 4 note: ΔS° = ΔG°37 × 1000/310.15.
    private const double ReferenceTemperatureKelvin = 310.15;

    // Jacobson-Stockmayer entropic extrapolation coefficient for loop sizes beyond the
    // tabulated lengths. SantaLucia & Hicks (2004) Eq. 7: ΔG°37(loop-n) =
    // ΔG°37(loop-x) + 2.44·R·310.15·ln(n/x); the 2.44 coefficient is from recent DNA
    // kinetics measurements (ref 22), preferred over the older 1.75.
    private const double JacobsonStockmayerCoefficient = 2.44;

    /// <summary>
    /// Hairpin loop ΔG°37 (kcal/mol) for a loop of <paramref name="loopSize"/> unpaired
    /// nucleotides. Tabulated sizes return Table 4 directly; non-tabulated sizes are filled by
    /// the Jacobson-Stockmayer extrapolation from the largest tabulated size ≤ n
    /// (SantaLucia &amp; Hicks 2004 Eq. 7). Loop sizes &lt; 3 are sterically prohibited.
    /// </summary>
    private static double HairpinLoopDeltaG(int loopSize)
    {
        if (HairpinLoopInitiationDeltaG.TryGetValue(loopSize, out double dg))
            return dg;

        // Jacobson-Stockmayer from the largest tabulated x ≤ loopSize.
        int x = 0;
        foreach (int size in HairpinLoopInitiationDeltaG.Keys)
            if (size <= loopSize && size > x) x = size;

        return HairpinLoopInitiationDeltaG[x]
               + JacobsonStockmayerCoefficient * GasConstant * ReferenceTemperatureKelvin
                 * 1e-3 * Math.Log((double)loopSize / x);
    }

    private static bool IsWatsonCrickPair(char a, char b) =>
        (a == 'A' && b == 'T') || (a == 'T' && b == 'A') ||
        (a == 'G' && b == 'C') || (a == 'C' && b == 'G');

    /// <summary>
    /// The most stable intramolecular DNA hairpin found in an oligo: the closing stem span,
    /// stem length (base pairs), loop size, and the hairpin ΔH° (kcal/mol), ΔS° (cal/(K·mol)),
    /// and ΔG°37 (kcal/mol).
    /// </summary>
    /// <param name="StemStart">5'-most index (0-based) of the stem on the input strand.</param>
    /// <param name="StemEnd">3'-most index (0-based) of the stem on the input strand.</param>
    /// <param name="StemLength">Number of base pairs in the stem.</param>
    /// <param name="LoopSize">Number of unpaired loop nucleotides closed by the stem.</param>
    /// <param name="DeltaH">Hairpin ΔH° in kcal/mol.</param>
    /// <param name="DeltaS">Hairpin ΔS° in cal/(K·mol).</param>
    /// <param name="DeltaG37">Hairpin ΔG°37 in kcal/mol (negative = stable).</param>
    public readonly record struct HairpinResult(
        int StemStart, int StemEnd, int StemLength, int LoopSize,
        double DeltaH, double DeltaS, double DeltaG37);

    /// <summary>
    /// Finds the most stable (minimum ΔG°37) intramolecular DNA <b>hairpin</b> — a single
    /// Watson-Crick stem closing one hairpin loop — in <paramref name="sequence"/>, using the
    /// SantaLucia &amp; Hicks (2004) Table 1 nearest-neighbour stem stacks and their
    /// Table 4 hairpin-loop initiation increments. <b>Opt-in</b>: the duplex Tm methods
    /// are unchanged. Returns <c>null</c> when the sequence is empty, contains a non-ACGT
    /// character, or admits no hairpin at all (no stem of ≥ 2 bp can close a loop of ≥ 3 nt,
    /// e.g. a homopolymer such as poly-A).
    /// <para>
    /// Model: ΔG°37 = Σ stem NN stacks (Table 1) + ΔG°37(loop of N) (Table 4); the bimolecular
    /// duplex-initiation term is intentionally excluded for this unimolecular structure. Loop
    /// ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15. The supplementary triloop/tetraloop and
    /// terminal-mismatch increments are not bundled (see <paramref name="loopBonusDeltaG37"/>).
    /// </para>
    /// </summary>
    /// <param name="sequence">DNA oligo (5'→3').</param>
    /// <param name="minStemLength">Minimum stem length in base pairs (default 2 → at least one
    /// NN stack). Must be ≥ 2.</param>
    /// <param name="loopBonusDeltaG37">Optional caller-supplied additive ΔG°37 increment
    /// (kcal/mol) for the terminal-mismatch / special triloop-tetraloop bonus that is NOT
    /// bundled (default 0). Added to the loop free energy; its ΔS° contribution follows the
    /// same −ΔG·1000/310.15 rule (ΔH° contribution 0), consistent with the Table 4 loop model.</param>
    /// <returns>The most stable hairpin, or <c>null</c> if none exists / invalid input.</returns>
    public static HairpinResult? FindMostStableHairpin(
        string sequence,
        int minStemLength = 2,
        double loopBonusDeltaG37 = 0.0)
    {
        if (string.IsNullOrEmpty(sequence) || minStemLength < 2)
            return null;

        string seq = sequence.ToUpperInvariant();
        int n = seq.Length;
        foreach (char c in seq)
            if (c is not ('A' or 'C' or 'G' or 'T'))
                return null; // non-ACGT base present

        HairpinResult? best = null;

        // For every candidate outermost closing pair (i, j), extend the stem inward as far as
        // Watson-Crick pairing allows, then close the remaining inner bases as a hairpin loop.
        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (!IsWatsonCrickPair(seq[i], seq[j]))
                    continue;

                double dH = 0.0, dS = 0.0;
                // Extend the stem from the outermost pair (i, j) inward.
                int a = i, b = j;
                int stemPairs = 0;
                while (a < b && IsWatsonCrickPair(seq[a], seq[b]))
                {
                    if (stemPairs > 0)
                    {
                        // NN stack between pair (a-1, b+1) and (a, b): key = 5'-strand dinucleotide seq[a-1..a].
                        string step = seq.Substring(a - 1, 2);
                        if (!TryGetDesignStack(step, out var p))
                            break;
                        dH += p.DeltaH;
                        dS += p.DeltaS;
                    }
                    stemPairs++;
                    a++;
                    b--;
                }

                if (stemPairs < minStemLength)
                    continue;

                // Innermost pair is (a-1, b+1); loop is the bases strictly between them.
                int innerLeft = a - 1;
                int innerRight = b + 1;
                int loopSize = innerRight - innerLeft - 1;
                if (loopSize < MinHairpinLoopSize)
                    continue;

                double loopDg = HairpinLoopDeltaG(loopSize) + loopBonusDeltaG37;
                // Loop ΔH° = 0; loop ΔS° = −ΔG°37·1000/310.15 (destabilising loop).
                double loopDs = -loopDg * 1000.0 / ReferenceTemperatureKelvin;

                double totalDh = dH;                       // loop ΔH° contribution is 0
                double totalDs = dS + loopDs;
                double dG37 = totalDh - ReferenceTemperatureKelvin * totalDs / 1000.0;

                if (best is null || dG37 < best.Value.DeltaG37)
                    best = new HairpinResult(i, j, stemPairs, loopSize, totalDh, totalDs, dG37);
            }
        }

        return best;
    }

    /// <summary>
    /// Computes the secondary-structure (hairpin) melting temperature (°C) of a DNA oligo:
    /// finds its most stable intramolecular hairpin (<see cref="FindMostStableHairpin"/>) and
    /// returns the <b>unimolecular</b> two-state Tm = ΔH°·1000/ΔS° − 273.15
    /// (SantaLucia &amp; Hicks 2004, Eq. 11). A hairpin is intramolecular, so the Tm is
    /// concentration-independent: there is <b>no</b> R·ln(C_T/x) strand-concentration term.
    /// <b>Opt-in</b>: the duplex Tm methods (<see cref="CalculateMeltingTemperatureNN"/>) and the
    /// default <see cref="CalculateMeltingTemperature(string)"/> are unchanged.
    /// </summary>
    /// <param name="sequence">DNA oligo (5'→3').</param>
    /// <param name="minStemLength">Minimum stem length in base pairs (default 2).</param>
    /// <param name="loopBonusDeltaG37">Optional caller-supplied terminal-mismatch / special-loop
    /// ΔG°37 increment (default 0; not bundled — see <see cref="FindMostStableHairpin"/>).</param>
    /// <returns>The hairpin Tm in °C, or <c>double.NaN</c> if no hairpin exists / invalid input.</returns>
    public static double CalculateHairpinMeltingTemperature(
        string sequence,
        int minStemLength = 2,
        double loopBonusDeltaG37 = 0.0)
    {
        var hairpin = FindMostStableHairpin(sequence, minStemLength, loopBonusDeltaG37);
        if (hairpin is null)
            return double.NaN;

        var h = hairpin.Value;
        // Unimolecular: NO concentration term (Eq. 11).
        return (h.DeltaH * 1000.0) / h.DeltaS - KelvinOffset;
    }

    // ---- Self-dimer / hetero-dimer (intermolecular) Tm via thermodynamic alignment ----
    // PRIMER-TM-001, opt-in. Finds the most stable INTERMOLECULAR antiparallel duplex
    // between two oligonucleotides (self-dimer = an oligo against a second copy of itself;
    // hetero/cross-dimer = two different oligos) and returns its NN ΔH°/ΔS° and the
    // bimolecular Tm. Reuses the SantaLucia & Hicks (2004) Table 1 NN stacking table
    // (ThermoConstants NnParameterSet.SantaLuciaHicks2004), the terminal-A·T penalty, the duplex-initiation term and the
    // 0.368 entropy salt coefficient already used by CalculateMeltingTemperatureNN.
    // The duplex Tm / hairpin Tm / default Tm methods and their defaults are UNCHANGED.
    //
    // Model (Primer3 / ntthal — SantaLucia & Hicks 2004 unified NN):
    //   For each gapless antiparallel offset of strand2 (read 3'→5') under strand1 (5'→3'),
    //   each maximal contiguous run of Watson-Crick pairs (≥ 1 NN stack) is a candidate
    //   duplex with
    //     ΔH° = ΔH°_init + Σ stacks ΔH° + ΔH°_AT-penalty(per A·T-closed end),
    //     ΔS° = ΔS°_init + Σ stacks ΔS° + ΔS°_AT-penalty(per A·T-closed end)
    //           + 0.368·N_stacks·ln[Na⁺]   (salt correction baked into ΔS°, ntthal saltCorrectS),
    //   and bimolecular
    //     Tm = ΔH°·1000 / (ΔS° + R·ln(C_T / x)) − 273.15,
    //   with x = 1 when BOTH oligos are reverse-complement palindromes (ntthal symmetry_thermo),
    //   else x = 4. ntthal keeps the candidate with the highest Tm — so does this method.
    //
    // Sources (retrieved & extracted this session, 2026-06-25):
    //   SantaLucia J, Hicks D (2004) Annu Rev Biophys 33:415-440 — unified NN parameters
    //     (Table 1, ThermoConstants NnParameterSet.SantaLuciaHicks2004) + the bimolecular Tm Eq. 3 + Eq. 5 entropy
    //     salt correction (0.368 coefficient).
    //   Untergasser A et al. (2012) Nucleic Acids Res 40:e115 (Primer3 2.0) — the ntthal
    //     thermodynamic-alignment engine for oligo dimers.
    //   Primer3 `thal.c` (primer3-py vendored libprimer3, retrieved
    //     https://raw.githubusercontent.com/libnano/primer3-py/master/primer3/src/libprimer3/thal.c):
    //     dplx_init_H=200 cal, dplx_init_S=−5.7 (lines 588-589); AT_H=2200, AT_S=6.9
    //     (lines 128-129); saltCorrectS = 0.368·ln((mv+120·√max(0,dv−dntp))/1000) per stack
    //     (lines 623-624, 1042); RC = R·ln(dna_conc/1e9) when both strands symmetric else
    //     R·ln(dna_conc/4e9) (lines 590-593); symmetry_thermo = reverse-complement palindrome
    //     (line 2771). dna_conc is in nM, so /1e9 → mol/L with x=1 and /4e9 → x=4.
    //   Cross-checked against primer3-py 2.3.0 calc_homodimer / calc_heterodimer
    //     (mv=50, dv=0, dntp=0, dna_conc=50 nM): this method reproduces ntthal's ΔH°, ΔS°
    //     and Tm to machine precision for every case whose optimal structure is a contiguous
    //     Watson-Crick duplex (e.g. GCGCGCGC, ACGTACGTACGT, ATCGATCGATCG/CGATCGATCGAT,
    //     CGATCGATCG self-dimer, GCATGC, GGGGCCCC). ntthal's extra terminal-stack /
    //     overhang-extension terms for some sequences are NOT modelled here (documented limit).

    /// <summary>
    /// The most stable intermolecular DNA duplex (self- or hetero-dimer) found between two
    /// oligonucleotides: the aligned spans on each strand, the number of base pairs, and the
    /// dimer ΔH° (kcal/mol), ΔS° (cal/(K·mol)) and ΔG°37 (kcal/mol).
    /// </summary>
    /// <param name="Strand1Start">5'-most aligned index (0-based) on strand 1.</param>
    /// <param name="Strand2Start">5'-most aligned index (0-based) on strand 2.</param>
    /// <param name="BasePairs">Number of contiguous Watson-Crick base pairs in the duplex.</param>
    /// <param name="DeltaH">Dimer ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Dimer ΔS° in cal/(K·mol), including the 0.368 salt correction.</param>
    /// <param name="DeltaG37">Dimer ΔG°37 = ΔH° − 310.15·ΔS°/1000 in kcal/mol (negative = stable).</param>
    public readonly record struct DimerResult(
        int Strand1Start, int Strand2Start, int BasePairs,
        double DeltaH, double DeltaS, double DeltaG37);

    /// <summary>
    /// Finds the most stable (highest-Tm) intermolecular DNA duplex between two oligonucleotides
    /// using the Primer3 / <c>ntthal</c> thermodynamic alignment over the SantaLucia &amp; Hicks
    /// (2004) unified nearest-neighbour model. A <b>self-dimer</b> is obtained by passing the same
    /// sequence as both strands; a <b>hetero/cross-dimer</b> by passing two different sequences.
    /// <b>Opt-in</b>: the duplex (<see cref="CalculateMeltingTemperatureNN"/>) and hairpin Tm
    /// methods, and the default <see cref="CalculateMeltingTemperature(string)"/>, are unchanged.
    /// <para>
    /// The two strands are aligned antiparallel (strand 2 read 3'→5' under strand 1 5'→3') over
    /// every gapless offset; each maximal contiguous Watson-Crick run of ≥ 2 bp is scored as a
    /// duplex with ΔH° = init + Σ stacks + terminal-A·T penalty per A·T-closed end, and ΔS° the
    /// same plus the 0.368·N<sub>stacks</sub>·ln[Na⁺] salt correction. The candidate with the
    /// highest bimolecular Tm is returned (the ΔS° in the result already includes the salt term
    /// for the supplied <paramref name="sodiumMolar"/>).
    /// </para>
    /// </summary>
    /// <param name="strand1">First DNA oligo (5'→3'). Must contain ≥ 2 ACGT bases.</param>
    /// <param name="strand2">Second DNA oligo (5'→3'); the same string as
    /// <paramref name="strand1"/> for a self-dimer. Must contain ≥ 2 ACGT bases.</param>
    /// <param name="sodiumMolar">Monovalent cation concentration in mol/L (default 50 mM); only
    /// the entropy salt correction depends on it.</param>
    /// <param name="strandConcentrationMolar">Total strand concentration C_T in mol/L for the
    /// bimolecular Tm used to rank candidates (default 50 nM, the Primer3/ntthal convention).</param>
    /// <returns>The most stable dimer, or <c>null</c> if either strand is null/&lt; 2 bases/contains
    /// a non-ACGT character, or no duplex of ≥ 2 contiguous base pairs exists between them.</returns>
    public static DimerResult? FindMostStableDimer(
        string strand1,
        string strand2,
        double sodiumMolar = ThermoConstants.DefaultNaConcentration,
        double strandConcentrationMolar = DefaultDimerStrandConcentrationMolar)
    {
        if (string.IsNullOrEmpty(strand1) || string.IsNullOrEmpty(strand2))
            return null;

        string s1 = strand1.ToUpperInvariant();
        string s2 = strand2.ToUpperInvariant();
        if (s1.Length < 2 || s2.Length < 2)
            return null;
        foreach (char c in s1)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;
        foreach (char c in s2)
            if (c is not ('A' or 'C' or 'G' or 'T')) return null;

        // x = 1 only when both strands are reverse-complement palindromes (ntthal symmetry_thermo).
        bool symmetric = IsSelfComplementary(s1) && IsSelfComplementary(s2);
        double x = symmetric ? SelfComplementaryFactor : NonSelfComplementaryFactor;
        // Strand-concentration term R·ln(C_T / x) for the bimolecular Tm (constant over candidates).
        double rcTerm = GasConstant * Math.Log(strandConcentrationMolar / x);
        double saltPerStack = SantaLuciaEntropySaltCoefficient * Math.Log(sodiumMolar);

        // Strand 2 read 3'→5' so its index i pairs base-for-base under strand 1 read 5'→3'.
        string s2Rev = Reverse(s2);
        int n = s1.Length, m = s2.Length;

        DimerResult? best = null;
        double bestTm = double.NegativeInfinity;

        // Slide strand 2 across strand 1 over every gapless antiparallel offset.
        for (int offset = -(m - 1); offset < n; offset++)
        {
            // Within this offset, walk the overlap and split it into maximal contiguous
            // Watson-Crick runs; each run of ≥ 2 bp is a candidate duplex.
            int runStart = -1; // strand-1 index where the current WC run started
            for (int i = 0; i <= n; i++)
            {
                int j = i - offset; // index into s2Rev paired with s1[i]
                bool paired = i < n && j >= 0 && j < m && IsWatsonCrickPair(s1[i], s2Rev[j]);

                if (paired && runStart < 0)
                    runStart = i;

                if (!paired && runStart >= 0)
                {
                    int runEnd = i - 1; // inclusive strand-1 index of the run end
                    EvaluateRun(s1, s2Rev, offset, runStart, runEnd, saltPerStack, rcTerm,
                                ref best, ref bestTm);
                    runStart = -1;
                }
            }
        }

        return best;

        // Scores one contiguous Watson-Crick run [runStart..runEnd] on strand 1 and keeps it
        // if its bimolecular Tm exceeds the best found so far.
        static void EvaluateRun(
            string s1, string s2Rev, int offset, int runStart, int runEnd,
            double saltPerStack, double rcTerm, ref DimerResult? best, ref double bestTm)
        {
            int basePairs = runEnd - runStart + 1;
            if (basePairs < MinDimerBasePairs)
                return;

            double dH = DesignNnInitiation.Initiation.DeltaH;
            double dS = DesignNnInitiation.Initiation.DeltaS;
            for (int k = runStart; k < runEnd; k++)
            {
                string step = s1.Substring(k, 2);
                // A contiguous WC run pairs every column; the stack is the perfect-match NN
                // keyed by the strand-1 dinucleotide (its complement is the strand-2 stack).
                if (!TryGetDesignStack(step, out var p))
                    return;
                dH += p.DeltaH;
                dS += p.DeltaS;
            }

            // Terminal A·T penalty per duplex end closing with an A·T pair.
            if (s1[runStart] is 'A' or 'T') { dH += DesignNnInitiation.TerminalAT.DeltaH; dS += DesignNnInitiation.TerminalAT.DeltaS; }
            if (s1[runEnd] is 'A' or 'T') { dH += DesignNnInitiation.TerminalAT.DeltaH; dS += DesignNnInitiation.TerminalAT.DeltaS; }

            int stacks = basePairs - 1;
            dS += stacks * saltPerStack; // 0.368·N_stacks·ln[Na⁺] (ntthal saltCorrectS)

            double tmKelvin = (dH * 1000.0) / (dS + rcTerm);
            double tmCelsius = tmKelvin - KelvinOffset;

            if (tmCelsius > bestTm)
            {
                bestTm = tmCelsius;
                int strand2Start5 = s2Rev.Length - 1 - (runEnd - offset); // → 5'→3' index on strand 2
                double dG37 = dH - ReferenceTemperatureKelvin * dS / 1000.0;
                best = new DimerResult(runStart, strand2Start5, basePairs, dH, dS, dG37);
            }
        }
    }

    // Minimum base pairs for a dimer duplex (at least one NN stack). ntthal requires a paired
    // region; a single base pair has no stacking energy and is not a duplex.
    private const int MinDimerBasePairs = 2;

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
    /// Full <c>ntthal</c> dimer thermodynamics of the most stable intermolecular duplex.
    /// </summary>
    /// <param name="DeltaH">Dimer ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Dimer ΔS° in cal/(K·mol), including the N·saltCorrection term.</param>
    /// <param name="DeltaG37">Dimer ΔG° = ΔH° − T·ΔS°/1000 in kcal/mol (negative = stable) at the
    /// analysis temperature T (310.15 K = 37 °C unless an overload with <c>temperatureCelsius</c> is used;
    /// ntthal / primer3-py <c>temp_c</c>). ntthal may report a positive ΔG for its optimal structure.</param>
    /// <param name="TmCelsius">Bimolecular melting temperature in °C.</param>
    /// <param name="BasePairs">Number of paired bases in the optimal structure.</param>
    public readonly record struct DimerThermodynamics(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius, int BasePairs);

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
    /// Full <c>ntthal</c> intramolecular-hairpin thermodynamics of the most stable hairpin
    /// (reproduces primer3-py <c>calc_hairpin</c>; special tri/tetraloop bonuses applied).
    /// </summary>
    /// <param name="DeltaH">Hairpin ΔH° in kcal/mol (salt-independent).</param>
    /// <param name="DeltaS">Hairpin ΔS° in cal/(K·mol), including the (N/2−1)·saltCorrection term.</param>
    /// <param name="DeltaG37">Hairpin ΔG = ΔH° − T·ΔS°/1000 in kcal/mol (negative = stable) at T = 310.15 K
    /// (37 °C), or at <c>temperatureCelsius</c> + 273.15 for the overloads taking it (primer3-py <c>temp_c</c>).</param>
    /// <param name="TmCelsius">Unimolecular melting temperature in °C (no strand-concentration term).</param>
    /// <param name="BasePairs">ntthal N/2: half the number of paired positions among bases 1..len−1 of the
    /// optimal structure (the count thal.c uses in the (N/2 − 1)·saltCorrection term).</param>
    public readonly record struct HairpinThermodynamics(
        double DeltaH, double DeltaS, double DeltaG37, double TmCelsius, int BasePairs);

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
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: false)?.Thermodynamics;

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
            strandConcentrationMolar, temperatureCelsius, maxLoop, withStructure: true);

    private static NtthalDimerStructure? CalculateDimerStructureNtthal(
        string strand1, string strand2, NtthalAlignmentMode mode, double sodiumMolar, double divalentMolar,
        double dntpMolar, double strandConcentrationMolar, double temperatureCelsius, int maxLoop, bool withStructure)
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
        var r = NtthalDimer.Run(strand1.ToUpperInvariant(), strand2.ToUpperInvariant(),
            sodiumMolar, strandConcentrationMolar, type, divalentMolar, dntpMolar,
            temperatureCelsius + KelvinOffset, maxLoop, withStructure);
        if (r is null)
            return null;
        var v = r.Value;
        return new NtthalDimerStructure(
            new DimerThermodynamics(v.DeltaH / 1000.0, v.DeltaS, v.DeltaG37 / 1000.0, v.TmCelsius, v.BasePairs),
            v.AsciiStructure ?? Array.Empty<string>());
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
            temperatureCelsius, maxLoop, withStructure: false)?.Thermodynamics;

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
            temperatureCelsius, maxLoop, withStructure: true);

    private static NtthalHairpinStructure? CalculateHairpinStructureNtthal(
        string sequence, double sodiumMolar, double divalentMolar, double dntpMolar,
        double temperatureCelsius, int maxLoop, bool withStructure)
    {
        if (!IsAcgtOnly(sequence))
            return null;
        var r = NtthalHairpin.Run(sequence.ToUpperInvariant(), sodiumMolar, divalentMolar, dntpMolar,
            temperatureCelsius + KelvinOffset, maxLoop, withStructure);
        if (r is null)
            return null;
        var v = r.Value;
        return new NtthalHairpinStructure(
            new HairpinThermodynamics(v.DeltaH / 1000.0, v.DeltaS, v.DeltaG37 / 1000.0, v.TmCelsius, v.BasePairs),
            v.AsciiStructure ?? Array.Empty<string>());
    }

    /// <summary>
    /// ntthal hairpin thermodynamics plus the thal.c <c>drawHairpin</c> ASCII structure
    /// (primer3-py <c>ThermoResult.ascii_structure_lines</c>).
    /// </summary>
    /// <param name="Thermodynamics">ΔH/ΔS/ΔG/Tm of the optimal hairpin.</param>
    /// <param name="AsciiStructureLines">The "SEQ\t…" and "STR\t…" lines.</param>
    public sealed record NtthalHairpinStructure(HairpinThermodynamics Thermodynamics, IReadOnlyList<string> AsciiStructureLines);

    private static bool IsAcgtOnly(string? s)
    {
        if (string.IsNullOrEmpty(s)) return false;
        foreach (char c in s)
            if (c is not ('A' or 'C' or 'G' or 'T' or 'a' or 'c' or 'g' or 't')) return false;
        return true;
    }

    /// <summary>
    /// Primer3's default (thermodynamic) limit, in °C, on every secondary-structure Tm:
    /// <c>PRIMER_MAX_SELF_ANY_TH</c>, <c>PRIMER_MAX_SELF_END_TH</c>, <c>PRIMER_MAX_HAIRPIN_TH</c>,
    /// <c>PRIMER_PAIR_MAX_COMPL_ANY_TH</c>, <c>PRIMER_PAIR_MAX_COMPL_END_TH</c> = 47.0
    /// (<c>libprimer3.cc</c> <c>pr_set_default_global_args_2</c>). A value is a violation when it is
    /// strictly greater than the limit.
    /// </summary>
    public const double Primer3MaxStructureTm = 47.0;

    /// <summary>
    /// Primer3 thermodynamic secondary-structure values of one primer (Tm in °C; 0 when ntthal finds
    /// no structure or the Tm is below 0 °C, as <c>align_thermod</c> reports): <c>PRIMER_*_SELF_ANY_TH</c>, <c>PRIMER_*_SELF_END_TH</c>, <c>PRIMER_*_HAIRPIN_TH</c>.
    /// </summary>
    /// <param name="SelfAnyTh">Self-dimer Tm, ntthal ANY of (primer, primer).</param>
    /// <param name="SelfEndTh">3′-anchored self-dimer Tm, ntthal END1 of (primer, primer).</param>
    /// <param name="HairpinTh">Hairpin Tm, ntthal HAIRPIN of the primer.</param>
    public readonly record struct Primer3OligoStructure(double SelfAnyTh, double SelfEndTh, double HairpinTh)
    {
        /// <summary>True when any value exceeds <paramref name="maxTm"/> (Primer3 rejects the primer).</summary>
        public bool Exceeds(double maxTm = Primer3MaxStructureTm) =>
            SelfAnyTh > maxTm || SelfEndTh > maxTm || HairpinTh > maxTm;
    }

    /// <summary>
    /// Primer3 thermodynamic pair complementarity (Tm in °C; 0 when no structure or Tm &lt; 0 °C):
    /// <c>PRIMER_PAIR_COMPL_ANY_TH</c> and <c>PRIMER_PAIR_COMPL_END_TH</c>.
    /// </summary>
    /// <param name="ComplAnyTh">Hetero-dimer Tm, ntthal ANY of (left, right).</param>
    /// <param name="ComplEndTh">Max of ntthal END1/END2 of (left, right) and of (rc(right), rc(left)).</param>
    public readonly record struct Primer3PairComplementarity(double ComplAnyTh, double ComplEndTh)
    {
        /// <summary>True when either value exceeds <paramref name="maxTm"/> (Primer3 rejects the pair).</summary>
        public bool Exceeds(double maxTm = Primer3MaxStructureTm) => ComplAnyTh > maxTm || ComplEndTh > maxTm;
    }

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
        double monovalentMillimolar = 50.0,
        double divalentMillimolar = 1.5,
        double dntpMillimolar = 0.6,
        double dnaConcentrationNanomolar = 50.0)
    {
        if (!IsAcgtOnly(primer))
            return null;
        string p = primer.ToUpperInvariant();
        double mv = monovalentMillimolar / 1000.0, dv = divalentMillimolar / 1000.0, dntp = dntpMillimolar / 1000.0;
        double conc = dnaConcentrationNanomolar * 1e-9;
        double any = TmOrZero(NtthalDimer.Run(p, p, mv, conc, NtthalDimer.AlignmentType.Any, dv, dntp));
        double end = TmOrZero(NtthalDimer.Run(p, p, mv, conc, NtthalDimer.AlignmentType.End1, dv, dntp));
        var h = NtthalHairpin.Run(p, mv, dv, dntp);
        return new Primer3OligoStructure(any, end, h is null ? 0.0 : Math.Max(0.0, h.Value.TmCelsius));
    }

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
        double monovalentMillimolar = 50.0,
        double divalentMillimolar = 1.5,
        double dntpMillimolar = 0.6,
        double dnaConcentrationNanomolar = 50.0)
    {
        if (!IsAcgtOnly(leftPrimer) || !IsAcgtOnly(rightPrimer))
            return null;
        string l = leftPrimer.ToUpperInvariant(), r = rightPrimer.ToUpperInvariant();
        double mv = monovalentMillimolar / 1000.0, dv = divalentMillimolar / 1000.0, dntp = dntpMillimolar / 1000.0;
        double conc = dnaConcentrationNanomolar * 1e-9;
        var (any, end) = Primer3PairTms(l, r, mv, dv, dntp, conc, double.PositiveInfinity);
        return new Primer3PairComplementarity(any, end);
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

    // Watson-Crick complement (same left-to-right order) via the canonical Core per-base complement
    // (SequenceExtensions.TryGetComplement); identical to the former local ACGT switch on the
    // validated ACGT duplexes it is applied to.
    private static string Complement(string seq) =>
        string.Create(seq.Length, seq, static (dest, src) => src.AsSpan().TryGetComplement(dest));

    /// <summary>True if the sequence equals its own reverse complement (self-complementary).</summary>
    private static bool IsSelfComplementary(string seq)
    {
        int n = seq.Length;
        if (n % 2 != 0) return false; // odd-length cannot be self-complementary
        for (int i = 0; i < n; i++)
        {
            char a = seq[i];
            char b = seq[n - 1 - i];
            bool pair = (a == 'A' && b == 'T') || (a == 'T' && b == 'A')
                     || (a == 'G' && b == 'C') || (a == 'C' && b == 'G');
            if (!pair) return false;
        }
        return true;
    }

    /// <summary>
    /// Generates all possible primers for a region.
    /// </summary>
    public static IEnumerable<PrimerCandidate> GeneratePrimerCandidates(
        DnaSequence template,
        int regionStart,
        int regionEnd,
        bool forward = true,
        PrimerParameters? parameters = null)
    {
        var param = parameters ?? DefaultParameters;

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

    // ---- Primer3 weighted penalty objective (PRIMER-TM-001) -------------------
    // Reproduces the per-primer objective function `p_obj_fn` (left/right primer
    // branch) from Primer3's reference source `libprimer3.cc`, with the documented
    // default weights and optima. Lower penalty = better primer, exactly as Primer3.
    // Source: Primer3 source libprimer3.cc p_obj_fn / pr_set_default_global_args_2;
    //         Primer3 manual §19 "HOW PRIMER3 CALCULATES THE PENALTY VALUE";
    //         Untergasser et al. (2012) NAR 40(15):e115; Koressaar & Remm (2007).

    /// <summary>
    /// Primer3 default per-primer objective weights, taken verbatim from
    /// <c>pr_set_default_global_args_2</c> in Primer3's <c>libprimer3.cc</c>:
    /// length/Tm weights = 1; GC, self-complementarity and N weights = 0.
    /// Source: Primer3 source (branch main), function pr_set_default_global_args_2.
    /// </summary>
    public static readonly Primer3PenaltyWeights DefaultPrimer3Weights = new(
        TmGt: 1.0,           // PRIMER_WT_TM_GT  (libprimer3.cc: weights.temp_gt = 1)
        TmLt: 1.0,           // PRIMER_WT_TM_LT  (weights.temp_lt = 1)
        SizeGt: 1.0,         // PRIMER_WT_SIZE_GT (weights.length_gt = 1)
        SizeLt: 1.0,         // PRIMER_WT_SIZE_LT (weights.length_lt = 1)
        GcGt: 0.0,           // PRIMER_WT_GC_PERCENT_GT (weights.gc_content_gt = 0)
        GcLt: 0.0,           // PRIMER_WT_GC_PERCENT_LT (weights.gc_content_lt = 0)
        SelfAny: 0.0,        // PRIMER_WT_SELF_ANY (weights.compl_any = 0)
        SelfEnd: 0.0,        // PRIMER_WT_SELF_END (weights.compl_end = 0)
        NumNs: 0.0,          // PRIMER_WT_NUM_NS  (weights.num_ns = 0)
        SelfAnyTh: 0.0,      // PRIMER_WT_SELF_ANY_TH (weights.compl_any_th = 0)
        SelfEndTh: 0.0,      // PRIMER_WT_SELF_END_TH (weights.compl_end_th = 0)
        HairpinTh: 0.0,      // PRIMER_WT_HAIRPIN_TH  (weights.hairpin_th = 0)
        EndStability: 0.0,   // PRIMER_WT_END_STABILITY (weights.end_stability = 0)
        ThermodynamicOligoAlignment: true // PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT = 1 (pr_set_default_global_args_2)
    );

    /// <summary>
    /// Primer3's fixed <c>weights.temp_cutoff</c> (= 5 °C, <c>pr_set_default_global_args_1</c>; not a
    /// user-settable tag). In thermodynamic mode a secondary-structure Tm within this many degrees of
    /// the primer Tm is penalised linearly, otherwise by a reciprocal term.
    /// </summary>
    public const double Primer3TempCutoff = 5.0;

    /// <summary>
    /// Primer3 default per-primer optima. OPT_TM = 60 °C and OPT_SIZE = 20 bases are
    /// from <c>libprimer3.cc</c> (opt_tm = 60.0, opt_size = 20); OPT_GC_PERCENT = 50.0 %
    /// is the published default in the Primer3 manual (PRIMER_OPT_GC_PERCENT, default 50.0).
    /// </summary>
    public static readonly Primer3Optima DefaultPrimer3Optima = new(
        OptTm: 60.0,         // PRIMER_OPT_TM (libprimer3.cc: opt_tm = 60.0) °C
        OptSize: 20,         // PRIMER_OPT_SIZE (libprimer3.cc: opt_size = 20) bases
        OptGcPercent: 50.0   // PRIMER_OPT_GC_PERCENT (manual default 50.0) %
    );

    /// <summary>
    /// Computes the Primer3 per-primer penalty (objective function value) for a single
    /// primer, faithfully reproducing the left/right-primer branch of Primer3's
    /// <c>p_obj_fn</c>. The penalty is the weighted sum of one-sided deviations of Tm,
    /// length and GC% from their optima, plus the weighted secondary-structure terms
    /// (alignment-score mode: <c>compl_any</c>/<c>compl_end</c>; thermodynamic mode, Primer3's
    /// default: <c>compl_any_th</c>/<c>compl_end_th</c>/<c>hairpin_th</c> with the fixed 5 °C
    /// <c>temp_cutoff</c>), the number-of-Ns term and the 3'-end-stability term. Each term is
    /// added only when its weight is non-zero (and, for Tm/GC/size, the deviation has the
    /// matching sign), so the result is always ≥ 0; <b>lower is better</b>, exactly as Primer3
    /// sorts candidates. Cross-checked against primer3-py 2.3.1 <c>design_primers</c>
    /// PRIMER_LEFT/RIGHT_n_PENALTY in both alignment modes.
    /// <para>Not modelled (all zero under Primer3 defaults): the annealing-temperature
    /// <c>bound</c> term (only when PRIMER_ANNEALING_TEMP &gt; 0), <c>failure_rate</c>,
    /// <c>repeat_sim</c> (needs a mispriming library), <c>pos_penalty</c> (needs
    /// PRIMER_INSIDE/OUTSIDE_PENALTY), <c>seq_quality</c> (needs base qualities) and
    /// <c>template_mispriming</c>; callers needing them add weight·value themselves.</para>
    /// </summary>
    /// <param name="inputs">Measured primer properties (Tm in °C, length in bases, GC in
    /// percent 0–100, self/3' local-alignment scores, count of N bases).</param>
    /// <param name="weights">Objective weights; defaults to <see cref="DefaultPrimer3Weights"/>.</param>
    /// <param name="optima">Parameter optima; defaults to <see cref="DefaultPrimer3Optima"/>.</param>
    /// <returns>The Primer3 objective-function value (penalty); 0 means every term is at its optimum.</returns>
    public static double CalculatePrimer3Penalty(
        Primer3PenaltyInputs inputs,
        Primer3PenaltyWeights? weights = null,
        Primer3Optima? optima = null)
    {
        var w = weights ?? DefaultPrimer3Weights;
        var o = optima ?? DefaultPrimer3Optima;

        double sum = 0.0;

        // Tm term: one-sided, separate _gt / _lt weights (p_obj_fn temp_gt / temp_lt).
        if (w.TmGt != 0 && inputs.Tm > o.OptTm)
            sum += w.TmGt * (inputs.Tm - o.OptTm);
        if (w.TmLt != 0 && inputs.Tm < o.OptTm)
            sum += w.TmLt * (o.OptTm - inputs.Tm);

        // GC% term (gc_content is a percentage 0–100 in libprimer3.cc, line 3856).
        if (w.GcGt != 0 && inputs.GcPercent > o.OptGcPercent)
            sum += w.GcGt * (inputs.GcPercent - o.OptGcPercent);
        if (w.GcLt != 0 && inputs.GcPercent < o.OptGcPercent)
            sum += w.GcLt * (o.OptGcPercent - inputs.GcPercent);

        // Length/size term (p_obj_fn length_lt / length_gt).
        if (w.SizeLt != 0 && inputs.Length < o.OptSize)
            sum += w.SizeLt * (o.OptSize - inputs.Length);
        if (w.SizeGt != 0 && inputs.Length > o.OptSize)
            sum += w.SizeGt * (inputs.Length - o.OptSize);

        // Secondary-structure terms: p_obj_fn switches on thermodynamic_oligo_alignment.
        if (!w.ThermodynamicOligoAlignment)
        {
            // Mode 0: local-alignment scores, linear weights (compl_any, compl_end).
            if (w.SelfAny != 0)
                sum += w.SelfAny * inputs.SelfAny;
            if (w.SelfEnd != 0)
                sum += w.SelfEnd * inputs.SelfEnd;
        }
        else
        {
            // Mode 1 (Primer3 default): SelfAny / SelfEnd / HairpinTh are structure Tm values (°C)
            // (compl_any_th, compl_end_th, hairpin_th).
            sum += ThermodynamicStructurePenalty(w.SelfAnyTh, inputs.Tm, inputs.SelfAny);
            sum += ThermodynamicStructurePenalty(w.SelfEndTh, inputs.Tm, inputs.SelfEnd);
            sum += ThermodynamicStructurePenalty(w.HairpinTh, inputs.Tm, inputs.HairpinTh);
        }

        // Number-of-Ns term (num_ns).
        if (w.NumNs != 0)
            sum += w.NumNs * inputs.NumNs;

        // 3'-end stability term (end_stability): weight · ΔG magnitude (kcal/mol, as Primer3 reports it).
        if (w.EndStability != 0)
            sum += w.EndStability * inputs.EndStability;

        return sum;
    }

    // p_obj_fn thermodynamic secondary-structure term (libprimer3.cc):
    //   if (Tm − temp_cutoff) ≤ s : w · (s − (Tm − temp_cutoff − 1))
    //   else                       : w · 1 / (Tm − temp_cutoff + 1 − s)
    private static double ThermodynamicStructurePenalty(double weight, double primerTm, double structureTm)
    {
        if (weight == 0)
            return 0;
        double threshold = primerTm - Primer3TempCutoff;
        return threshold <= structureTm
            ? weight * (structureTm - (threshold - 1.0))
            : weight * (1.0 / (threshold + 1.0 - structureTm));
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

    private static string Reverse(string s)
    {
        var chars = s.ToCharArray();
        Array.Reverse(chars);
        return new string(chars);
    }

    private static bool AreComplementary(string s1, string s2)
    {
        if (s1.Length != s2.Length) return false;
        for (int i = 0; i < s1.Length; i++)
        {
            if (!IsComplementary(s1[i], s2[i]))
                return false;
        }
        return true;
    }

    private static bool IsComplementary(char c1, char c2) =>
        (c1 == 'A' && c2 == 'T') || (c1 == 'T' && c2 == 'A') ||
        (c1 == 'G' && c2 == 'C') || (c1 == 'C' && c2 == 'G');
}

/// <summary>
/// Parameters for primer design.
/// </summary>
public readonly record struct PrimerParameters(
    int MinLength,
    int MaxLength,
    int OptimalLength,
    double MinGcContent,
    double MaxGcContent,
    double MinTm,
    double MaxTm,
    double OptimalTm,
    int MaxHomopolymer,
    int MaxDinucleotideRepeats,
    bool Avoid3PrimeGC,
    bool Check3PrimeStability,
    PrimerStructureScreen StructureScreen = PrimerStructureScreen.Primer3Thermodynamic,
    double MaxStructureTm = PrimerDesigner.Primer3MaxStructureTm)
{
    /// <summary>
    /// Limit (°C) for every Primer3 thermodynamic structure value under
    /// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> — Primer3's PRIMER_MAX_SELF_ANY_TH,
    /// PRIMER_MAX_SELF_END_TH, PRIMER_MAX_HAIRPIN_TH, PRIMER_PAIR_MAX_COMPL_ANY_TH and
    /// PRIMER_PAIR_MAX_COMPL_END_TH, all 47 °C by default. A value strictly greater fails.
    /// The value <c>0</c> (e.g. from <c>default(PrimerParameters)</c>) also means 47 °C.
    /// </summary>
    public double EffectiveMaxStructureTm => MaxStructureTm > 0 ? MaxStructureTm : PrimerDesigner.Primer3MaxStructureTm;
}

/// <summary>
/// Secondary-structure screen applied by <see cref="PrimerDesigner.EvaluatePrimer"/> and
/// <see cref="PrimerDesigner.DesignPrimers"/>.
/// </summary>
public enum PrimerStructureScreen
{
    /// <summary>
    /// Primer3's default thermodynamic screen (<c>PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1</c>): a
    /// primer is rejected when its ntthal self-dimer, 3′ self-dimer or hairpin Tm exceeds 47 °C
    /// (<see cref="PrimerDesigner.CalculatePrimer3OligoStructure"/>), a pair when its ntthal
    /// hetero-dimer or 3′ hetero-dimer Tm exceeds 47 °C
    /// (<see cref="PrimerDesigner.CalculatePrimer3PairComplementarity"/>).
    /// </summary>
    Primer3Thermodynamic = 0,

    /// <summary>
    /// Sequence-only screen: <see cref="PrimerDesigner.HasHairpinPotential"/> (a ≥ 4-bp
    /// Watson–Crick stem closing a ≥ 3-nt loop) per primer and <see cref="PrimerDesigner.HasPrimerDimer"/>
    /// (Primer3 alignment-mode pair 3′ complementarity ≥ 4) per pair.
    /// </summary>
    Heuristic = 1,
}

/// <summary>
/// A primer candidate with quality metrics. <see cref="MeltingTemperature"/> is the Primer3-default
/// Tm rounded to 0.1 °C; <see cref="Score"/> is an informational heuristic (higher is better);
/// <see cref="Penalty"/> is the unrounded Primer3 per-primer penalty (lower is better) that
/// <see cref="PrimerDesigner.DesignPrimers"/> ranks by. Under the default
/// <see cref="PrimerStructureScreen.Primer3Thermodynamic"/> screen <see cref="SelfAnyTh"/>,
/// <see cref="SelfEndTh"/> and <see cref="HairpinTh"/> carry Primer3's PRIMER_*_SELF_ANY_TH /
/// _SELF_END_TH / _HAIRPIN_TH (°C) and <see cref="HasHairpin"/> means HairpinTh &gt; 47 °C; under the
/// heuristic screen they are <c>null</c> and <see cref="HasHairpin"/> is
/// <see cref="PrimerDesigner.HasHairpinPotential"/>.
/// </summary>
public sealed record PrimerCandidate(
    string Sequence,
    int Position,
    bool IsForward,
    int Length,
    double GcContent,
    double MeltingTemperature,
    int HomopolymerLength,
    bool HasHairpin,
    double Stability3Prime,
    bool IsValid,
    IReadOnlyList<string> Issues,
    double Score,
    double Penalty = 0.0,
    double? SelfAnyTh = null,
    double? SelfEndTh = null,
    double? HairpinTh = null);

/// <summary>
/// Measured properties of a single primer used as input to the Primer3 penalty
/// objective (<see cref="PrimerDesigner.CalculatePrimer3Penalty"/>). Units mirror
/// Primer3's <c>p_obj_fn</c>: Tm in °C, length in bases, GC in percent (0–100),
/// self/3'-complementarity as local-alignment scores, and the count of N bases.
/// </summary>
/// <param name="Tm">Primer melting temperature in °C.</param>
/// <param name="Length">Primer length in bases.</param>
/// <param name="GcPercent">GC content as a percentage in [0, 100].</param>
/// <param name="SelfAny">Self-complementarity: the local-alignment score PRIMER_SELF_ANY when
/// <see cref="Primer3PenaltyWeights.ThermodynamicOligoAlignment"/> is false, or the self-dimer Tm in °C
/// (PRIMER_SELF_ANY_TH) when it is true — Primer3 stores both in the same field.</param>
/// <param name="SelfEnd">3'-self-complementarity: PRIMER_SELF_END score (alignment mode) or
/// PRIMER_SELF_END_TH Tm in °C (thermodynamic mode).</param>
/// <param name="NumNs">Number of ambiguous N bases in the primer.</param>
/// <param name="HairpinTh">Hairpin Tm in °C (PRIMER_HAIRPIN_TH); used only in thermodynamic mode.</param>
/// <param name="EndStability">3'-end stability ΔG magnitude in kcal/mol (PRIMER_END_STABILITY).</param>
public readonly record struct Primer3PenaltyInputs(
    double Tm,
    int Length,
    double GcPercent,
    double SelfAny = 0.0,
    double SelfEnd = 0.0,
    int NumNs = 0,
    double HairpinTh = 0.0,
    double EndStability = 0.0);

/// <summary>
/// Weights for the Primer3 per-primer penalty objective (the <c>PRIMER_WT_*</c>
/// parameters). Tm, GC and size each have separate "greater-than" (_gt) and
/// "less-than" (_lt) weights, applied one-sidedly relative to the optimum.
/// Defaults are <see cref="PrimerDesigner.DefaultPrimer3Weights"/>.
/// </summary>
/// <param name="TmGt">PRIMER_WT_TM_GT.</param>
/// <param name="TmLt">PRIMER_WT_TM_LT.</param>
/// <param name="SizeGt">PRIMER_WT_SIZE_GT.</param>
/// <param name="SizeLt">PRIMER_WT_SIZE_LT.</param>
/// <param name="GcGt">PRIMER_WT_GC_PERCENT_GT.</param>
/// <param name="GcLt">PRIMER_WT_GC_PERCENT_LT.</param>
/// <param name="SelfAny">PRIMER_WT_SELF_ANY (alignment mode only).</param>
/// <param name="SelfEnd">PRIMER_WT_SELF_END (alignment mode only).</param>
/// <param name="NumNs">PRIMER_WT_NUM_NS.</param>
/// <param name="SelfAnyTh">PRIMER_WT_SELF_ANY_TH (thermodynamic mode only).</param>
/// <param name="SelfEndTh">PRIMER_WT_SELF_END_TH (thermodynamic mode only).</param>
/// <param name="HairpinTh">PRIMER_WT_HAIRPIN_TH (thermodynamic mode only).</param>
/// <param name="EndStability">PRIMER_WT_END_STABILITY.</param>
/// <param name="ThermodynamicOligoAlignment">PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT: false = alignment-score
/// secondary-structure terms (SelfAny/SelfEnd weights), true = thermodynamic terms (…Th weights). Defaults to
/// false for an explicitly constructed weight set; <see cref="PrimerDesigner.DefaultPrimer3Weights"/> uses
/// Primer3's default (true).</param>
public readonly record struct Primer3PenaltyWeights(
    double TmGt,
    double TmLt,
    double SizeGt,
    double SizeLt,
    double GcGt,
    double GcLt,
    double SelfAny,
    double SelfEnd,
    double NumNs,
    double SelfAnyTh = 0.0,
    double SelfEndTh = 0.0,
    double HairpinTh = 0.0,
    double EndStability = 0.0,
    bool ThermodynamicOligoAlignment = false);

/// <summary>
/// Optimal parameter values for the Primer3 penalty objective
/// (<c>PRIMER_OPT_TM</c>, <c>PRIMER_OPT_SIZE</c>, <c>PRIMER_OPT_GC_PERCENT</c>).
/// Defaults are <see cref="PrimerDesigner.DefaultPrimer3Optima"/>.
/// </summary>
public readonly record struct Primer3Optima(
    double OptTm,
    int OptSize,
    double OptGcPercent);

/// <summary>
/// Result of primer pair design.
/// </summary>
public sealed record PrimerPairResult(
    PrimerCandidate? Forward,
    PrimerCandidate? Reverse,
    bool IsValid,
    string Message,
    int ProductSize);
