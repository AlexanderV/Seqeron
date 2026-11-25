// SPDX-License-Identifier: GPL-2.0-or-later
// Port of Primer3 (thal.c ntthal, dpal.c, masker.c, libprimer3; primer3-org/primer3 and primer3-py 2.3.1; GNU General Public License v2 or later).
// Part of Seqeron.Genomics.OligoThermodynamics, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Text;
using Seqeron.Genomics.MolTools;

namespace Seqeron.Genomics.OligoThermodynamics;

public static partial class ThermodynamicProbeDesigner
{
    // _pr_data_control: PRIMER_INTERNAL_WT_GC_PERCENT_GT/_LT need PRIMER_INTERNAL_OPT_GC_PERCENT.
    internal static void ValidatePrimer3ProbeGcOptimum(Primer3ProbeSettings s, string paramName)
    {
        if ((s.WeightGcPercentGt != 0 || s.WeightGcPercentLt != 0) && s.OptGcPercent is null)
            throw new ArgumentException(PrimerDesigner.ProbeGcOptimumUndefinedMessage, paramName);
    }

    // Primer3 settings check shared by the hybridization-probe picker and PrimerDesigner's
    // PRIMER_PICK_INTERNAL_OLIGO path.
    internal static void ValidatePrimer3ProbeSettings(Primer3ProbeSettings s, string paramName)
    {
        if (s.MinSize < 1 || s.MaxSize < s.MinSize || s.MaxSize > Primer3MaxOligoLength)
            throw new ArgumentOutOfRangeException(paramName,
                $"Probe sizes must satisfy 1 ≤ MinSize ≤ MaxSize ≤ {Primer3MaxOligoLength}.");
        if (!(s.MaxSelfAny >= 0 && s.MaxSelfAny <= short.MaxValue && s.MaxSelfEnd >= 0 && s.MaxSelfEnd <= short.MaxValue))
            throw new ArgumentOutOfRangeException(paramName,
                "Illegal value for internal oligo complementarity restrictions (Primer3: 0 ≤ PRIMER_INTERNAL_MAX_SELF_ANY/_END ≤ 32767).");
        PrimerDesigner.ValidatePrimer3Conditions(
            s.MonovalentMillimolar, s.DivalentMillimolar, s.DntpMillimolar, s.DnaConcentrationNanomolar, paramName);
        if (s.OptGcPercent is { } opt && !double.IsFinite(opt))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_INTERNAL_OPT_GC_PERCENT must be finite.");
        ValidatePrimer3ProbeGcOptimum(s, paramName);
        // _pr_data_control: PRIMER_INTERNAL_MAX_LIBRARY_MISHYB > SHRT_MAX (alignment mode); a mishyb weight without a library.
        if (double.IsNaN(s.MaxLibraryMishyb) || Math.Abs(s.MaxLibraryMishyb) >= int.MaxValue
            || (s.MaxLibraryMishyb > short.MaxValue && !s.ThermodynamicOligoAlignment))
            throw new ArgumentOutOfRangeException(paramName, "Value too large at tag PRIMER_INTERNAL_MAX_LIBRARY_MISHYB.");
        if (s.WeightLibraryMishyb != 0 && s.ActiveMishybLibrary is null)
            throw new ArgumentException(
                "Internal oligo mispriming score is part of objective function while mishyb library is not defined (Primer3 _pr_data_control).",
                paramName);
        // _pr_data_control: PRIMER_INTERNAL_OPT_BOUND within [MIN, MAX]; PRIMER_ANNEALING_TEMP ≤ 100.
        if (double.IsNaN(s.WeightBoundGt) || double.IsNaN(s.WeightBoundLt))
            throw new ArgumentOutOfRangeException(paramName, "PRIMER_INTERNAL_WT_BOUND_GT/_LT must not be NaN.");
        ThermodynamicPrimerDesigner.ValidatePrimer3Bound(s.AnnealingTemperature, s.MinBound, s.MaxBound, s.OptBound,
            internalOligo: true, paramName);
    }
}
