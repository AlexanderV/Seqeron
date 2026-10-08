namespace Seqeron.Genomics.MolTools;

// Primer3 template mispriming (PRIMER_MAX_TEMPLATE_MISPRIMING[_TH], PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING[_TH],
// PRIMER_WT_TEMPLATE_MISPRIMING[_TH], PRIMER_PAIR_WT_TEMPLATE_MISPRIMING[_TH], PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT):
// libprimer3.c oligo_template_mispriming, primer_mispriming_to_template (dpal DPAL_LOCAL_END, the F41 port),
// primer_mispriming_to_template_thermod (ntthal THAL_END1, NtthalDimer), _pr_need_[pair_]template_mispriming[_thermod],
// characterize_pair (pair score), p_obj_fn / obj_fn terms. Bit-exact to primer3-py 2.3.1 design_primers.
public static partial class PrimerDesigner
{
    /// <summary>
    /// Primer3 default PRIMER_MAX_TEMPLATE_MISPRIMING, PRIMER_MAX_TEMPLATE_MISPRIMING_TH,
    /// PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING and PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH (<c>pr_set_default_global_args_1</c>:
    /// <c>PR_UNDEFINED_ALIGN_OPT</c> = −100): a negative per-primer or alignment-mode pair limit is not checked.
    /// </summary>
    public const double Primer3UndefinedTemplateMispriming = -100.0;

    /// <summary>
    /// Primer3 template mispriming of one primer (<c>libprimer3.c</c> <c>oligo_template_mispriming</c>): the primer —
    /// 5′→3′ as synthesised, a left primer equal to the top strand at <paramref name="position"/>, a right primer the
    /// reverse complement of the top strand at <paramref name="position"/> — is aligned with the template strand it
    /// lies on, once against the part 5′ of its own site and once against the part 3′ of it
    /// (<see cref="TemplateMisprimingScore.SameStrand"/> = the larger: Primer3 <c>template_mispriming</c>), and with the
    /// whole opposite strand (<see cref="TemplateMisprimingScore.OtherStrand"/>, <c>template_mispriming_r</c>).
    /// <para>Alignment mode (PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = 0, Primer3's default,
    /// <c>primer_mispriming_to_template</c>): dpal local alignment anchored at the primer's 3′ end
    /// (<c>DPAL_LOCAL_END</c>, +1 match, −1 mismatch, −2 per single-base gap, max gap 1; a segment shorter than 3 nt
    /// scores its length), the same port as the mispriming library (<see cref="CalculateLibraryMispriming"/>).</para>
    /// <para>Thermodynamic mode (<paramref name="thermodynamic"/> = true, <c>primer_mispriming_to_template_thermod</c>):
    /// the ntthal THAL_END1 Tm (°C; Primer3 <c>use_end_for_th_template_mispriming</c> = 1) of the primer against each
    /// segment read on the complementary strand, at the primer reaction conditions (Primer3
    /// <c>create_thal_arg_holder(p_args)</c>); 0 for no structure, a negative Tm or an empty segment. Primer3 swaps the
    /// strands of the two calls relative to alignment mode (the 5′/3′ split is taken on the complementary strand), which
    /// is reproduced. The template must then be ≤ 10000 nt (ntthal THAL_MAX_SEQ).</para>
    /// </summary>
    /// <param name="template">Template (top strand, A/C/G/T).</param>
    /// <param name="position">0-based leftmost top-strand base of the primer's site (<see cref="PrimerCandidate.Position"/>).</param>
    /// <param name="length">Primer length.</param>
    /// <param name="isForward">True for a left primer, false for a right primer.</param>
    /// <param name="thermodynamic">PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT (default false = 0, Primer3's default).</param>
    /// <param name="monovalentMillimolar">PRIMER_SALT_MONOVALENT, mM (thermodynamic mode).</param>
    /// <param name="divalentMillimolar">PRIMER_SALT_DIVALENT, mM (thermodynamic mode).</param>
    /// <param name="dntpMillimolar">PRIMER_DNTP_CONC, mM (thermodynamic mode).</param>
    /// <param name="dnaConcentrationNanomolar">PRIMER_DNA_CONC, nM (thermodynamic mode).</param>
    /// <returns>The two Primer3 scores; <see cref="TemplateMisprimingScore.Max"/> is PRIMER_LEFT/RIGHT_n_TEMPLATE_MISPRIMING.</returns>
    /// <exception cref="ArgumentNullException">Null template.</exception>
    /// <exception cref="ArgumentOutOfRangeException">The site lies outside the template or <paramref name="length"/> &lt; 1.</exception>
    /// <exception cref="ArgumentException">Thermodynamic mode with a template longer than 10000 nt.</exception>
    public static TemplateMisprimingScore CalculateTemplateMispriming(
        DnaSequence template,
        int position,
        int length,
        bool isForward,
        bool thermodynamic = false,
        double monovalentMillimolar = Primer3MonovalentMillimolar,
        double divalentMillimolar = Primer3DivalentMillimolar,
        double dntpMillimolar = Primer3DntpMillimolar,
        double dnaConcentrationNanomolar = Primer3DnaConcentrationNanomolar)
    {
        ArgumentNullException.ThrowIfNull(template);
        string seq = template.Sequence.ToUpperInvariant();
        if (length < 1 || position < 0 || position > seq.Length - length)
            throw new ArgumentOutOfRangeException(nameof(position), "The primer site must lie inside the template.");
        ValidatePrimer3Conditions(monovalentMillimolar, divalentMillimolar, dntpMillimolar, dnaConcentrationNanomolar,
            nameof(monovalentMillimolar));
        var ctx = new TemplateContext(seq, thermodynamic, monovalentMillimolar, divalentMillimolar, dntpMillimolar,
            dnaConcentrationNanomolar);
        string site = seq.Substring(position, length);
        string primer = isForward ? site : DnaSequence.GetReverseComplementString(site);
        return ctx.Score(primer, position, isForward);
    }

    // The template strands (sa->upcased_seq / upcased_seq_r) and the alignment mode / thal conditions.
    internal sealed class TemplateContext
    {
        private readonly string _seq, _rc;
        private readonly bool _thermodynamic;
        private readonly double _mv, _dv, _dntp, _conc;

        public TemplateContext(string upperSeq, bool thermodynamic, double mvMm, double dvMm, double dntpMm, double dnaNm)
        {
            _seq = upperSeq;
            _rc = DnaSequence.GetReverseComplementString(upperSeq);
            _thermodynamic = thermodynamic;
            _mv = mvMm / 1000.0;
            _dv = dvMm / 1000.0;
            _dntp = dntpMm / 1000.0;
            _conc = dnaNm * 1e-9;
            if (thermodynamic && upperSeq.Length > NtthalDimer.ThalMaxSeq)
                throw new ArgumentException(
                    "Target sequence length > maximum allowed (10000) in thermodynamic alignment (PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = 1).");
        }

        // primer_mispriming_to_template / _thermod. first/last = the site on the top strand; for the strand on which
        // the oligo is matched against the template, the coordinates are flipped (n − 1 − last, n − 1 − first).
        public TemplateMisprimingScore Score(string oligo, int position, bool isForward)
        {
            int n = _seq.Length;
            int first = position, last = position + oligo.Length - 1;
            // Alignment mode: left → target = seq (unflipped), target_r = rc; right → target = rc (flipped), target_r = seq.
            // Thermodynamic mode: left → target = rc (flipped), target_r = seq; right → target = seq (unflipped), target_r = rc.
            bool targetIsRc = _thermodynamic ? isForward : !isForward;
            string target = targetIsRc ? _rc : _seq;
            string targetR = targetIsRc ? _seq : _rc;
            if (targetIsRc)
                (first, last) = (n - 1 - last, n - 1 - first);

            // 1. 5′ of the oligo: target[0, first); 2. 3′ of the oligo: target[last + 1, n); 3. the max; 4. other strand.
            double five = Align(oligo, target.Substring(0, first));
            double same = Align(oligo, target.Substring(last + 1));
            if (five > same) same = five;
            double other = Align(oligo, targetR);
            return new TemplateMisprimingScore(same, other);
        }

        // libprimer3.c align(…, local_end) / align_thermod(…, thal_args.end1).
        private double Align(string oligo, string target)
        {
            if (!_thermodynamic)
                return AlignLibrary(oligo, target, DpalPrimerMatrix, localEnd: true);
            if (target.Length == 0)
                return 0.0; // thal: "Empty second sequence" → temp 0
            return TmOrZero(NtthalDimer.Run(oligo, target, _mv, _conc, NtthalDimer.AlignmentType.End1, _dv, _dntp));
        }
    }
}

/// <summary>
/// Primer3 template mispriming scores of one primer (<c>oligo_template_mispriming</c>): dpal scores in alignment mode,
/// ntthal Tm (°C) in thermodynamic mode (PRIMER_THERMODYNAMIC_TEMPLATE_ALIGNMENT = 1).
/// </summary>
/// <param name="SameStrand">Best ectopic site on the strand the primer is matched against, outside its own site
/// (Primer3 <c>template_mispriming</c>).</param>
/// <param name="OtherStrand">Best site on the opposite strand (Primer3 <c>template_mispriming_r</c>).</param>
public readonly record struct TemplateMisprimingScore(double SameStrand, double OtherStrand)
{
    /// <summary>
    /// PRIMER_LEFT/RIGHT_n_TEMPLATE_MISPRIMING (Primer3 <c>oligo_max_template_mispriming</c>): the larger of the two.
    /// </summary>
    public double Max => SameStrand > OtherStrand ? SameStrand : OtherStrand;
}
