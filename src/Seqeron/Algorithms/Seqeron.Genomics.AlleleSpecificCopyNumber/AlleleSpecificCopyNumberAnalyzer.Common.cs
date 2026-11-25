// SPDX-License-Identifier: GPL-3.0-only
// Port of the R packages ASCAT 3.2.0 and Battenberg 3.0.0 (VanLoo-lab / Wedge-lab), both GPL-3; the R Mersenne-Twister generator (RMersenneTwister.cs) follows R (GPL-2 or GPL-3).
// Part of Seqeron.Genomics.AlleleSpecificCopyNumber, licensed under the GNU General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.AlleleSpecificCopyNumber;

/// <summary>Allele-specific copy-number analysis: shared constants and label helpers.</summary>
public static partial class AlleleSpecificCopyNumberAnalyzer
{
    /// <summary>Normal (germline) total copy number of an autosomal locus (diploid).</summary>
    private const double NormalDiploidCopyNumber = 2.0;

    /// <summary>R <c>mad</c> consistency constant 1.4826 = 1/Φ⁻¹(3/4).</summary>
    private const double MadConsistencyConstant = 1.4826;

    /// <summary>
    /// Removes a leading case-insensitive "chr" from a chromosome label ("chr7" → "7", "ChrX" → "X"); a label that is
    /// exactly "chr" (or shorter) is returned unchanged.
    /// </summary>
    private static ReadOnlySpan<char> StripChrPrefix(ReadOnlySpan<char> name) =>
        name.Length > 3 && name.StartsWith("chr", StringComparison.OrdinalIgnoreCase) ? name[3..] : name;
}
