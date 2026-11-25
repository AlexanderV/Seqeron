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
