// SPDX-License-Identifier: AGPL-3.0-only
// Port of RSAT oligo-analysis v1.169 / dyad-analysis v1.78 / RSA.disco.lib (rsa-tools/rsat-code, GNU Affero General Public License v3).
// Part of Seqeron.Genomics.OverrepresentedWords, licensed under the GNU Affero General Public License, version 3 (only) (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.OverrepresentedWords;

/// <summary>
/// Over-represented word statistics of RSAT <c>oligo-analysis</c> / <c>dyad-analysis</c>, and the exact PWM p-value under an
/// <see cref="OligoBackgroundModel"/> (the model's explicit word distribution is handed to the MIT p-value engine of
/// <see cref="MotifFinder"/>).
/// </summary>
public static partial class OverrepresentedWordAnalyzer
{
    /// <summary>
    /// Exact p-value of a PWM score under the Markov background <paramref name="background"/> of an RSAT model:
    /// P(S ≥ <paramref name="score"/>) for a random word whose probability is the RSAT <c>segment_proba</c> of the model.
    /// Supported models: <see cref="OligoBackgroundModel.Equiprobable"/>, <see cref="OligoBackgroundModel.Bernoulli"/> and
    /// <see cref="OligoBackgroundModel.MarkovFromOligoFrequencies"/> (order ≤ 10); the input-estimated and lexicon models throw.
    /// </summary>
    /// <exception cref="ArgumentNullException">An argument is null.</exception>
    /// <exception cref="ArgumentException">Non-finite matrix cell, unsupported model kind, or invalid options.</exception>
    /// <exception cref="ArgumentOutOfRangeException"><paramref name="score"/> is NaN or infinite.</exception>
    public static PwmPValueResult PwmMarkovScorePValue(
        PositionWeightMatrix pwm, double score, OligoBackgroundModel background, PwmPValueOptions? options = null)
    {
        if (!double.IsFinite(score))
            throw new ArgumentOutOfRangeException(nameof(score), score, "Score must be finite.");
        ArgumentNullException.ThrowIfNull(pwm);
        ArgumentNullException.ThrowIfNull(background);
        return MotifFinder.PwmMarkovScorePValue(pwm, score, background.ToPValueChain(nameof(background)), options);
    }

    /// <summary>
    /// Exact score threshold of a p-value under the Markov background <paramref name="background"/> of an RSAT model (the
    /// <see cref="MotifFinder.PwmScoreThresholdForPValue(PositionWeightMatrix, double, IReadOnlyList{double}?)"/> contract with the
    /// word distribution of <see cref="PwmMarkovScorePValue"/>).
    /// </summary>
    public static PwmPValueResult PwmMarkovScoreThresholdForPValue(
        PositionWeightMatrix pwm, double pValue, OligoBackgroundModel background, PwmPValueOptions? options = null)
    {
        if (!(pValue >= 0 && pValue <= 1))
            throw new ArgumentOutOfRangeException(nameof(pValue), pValue, "P-value must be in [0, 1].");
        ArgumentNullException.ThrowIfNull(pwm);
        ArgumentNullException.ThrowIfNull(background);
        return MotifFinder.PwmMarkovScoreThresholdForPValue(pwm, pValue, background.ToPValueChain(nameof(background)), options);
    }
}
