// SPDX-License-Identifier: GPL-2.0-or-later
// Port of EMBOSS 6.6.0 (application sources; the GNU General Public License, version 2 or later).
// Part of Seqeron.Genomics.Emboss.PairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using System.Globalization;
using System.Text;
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.Emboss.PairwiseAlignment;

using Seqeron.Genomics.Alignment;

/// <summary>The arguments of <see cref="EmbossEst2genomePrograms.Run"/>: the inputs, the options of EMBOSS <c>est2genome</c> 6.6.0 and the sequence qualifiers.</summary>
/// <param name="EstFasta">The spliced EST sequences (<c>-estsequence</c>, a <c>seqall</c>): FASTA text.</param>
/// <param name="GenomeFasta">The unspliced genomic sequence (<c>-genomesequence</c>): FASTA text (the first record is used).</param>
/// <param name="Match"><c>-match</c>: score for matching two bases (default 1).</param>
/// <param name="Mismatch"><c>-mismatch</c>: cost for mismatching two bases (default 1; the absolute value is used).</param>
/// <param name="GapPenalty"><c>-gappenalty</c>: cost for deleting a single base in either sequence, excluding introns (default 2; absolute value).</param>
/// <param name="IntronPenalty"><c>-intronpenalty</c>: cost for an intron, independent of length (default 40; absolute value).</param>
/// <param name="SplicePenalty"><c>-splicepenalty</c>: cost for an intron between donor and acceptor sites (default 20; absolute value).</param>
/// <param name="MinScore"><c>-minscore</c>: exclude alignments with a score below this (default 30).</param>
/// <param name="Reverse"><c>-reverse</c>: print the EST positions of the reversed EST on the original strand (default false).</param>
/// <param name="UseSplice"><c>-usesplice</c>: use donor and acceptor splice sites (default true).</param>
/// <param name="Mode"><c>-mode</c>: <c>both</c> (default), <c>forward</c> or <c>reverse</c> (case-insensitive, unique abbreviations as the ACD menu accepts).</param>
/// <param name="Best"><c>-best</c>: print only the best alignment (default true; <c>-nobest</c> with a single strand mode makes the program crash).</param>
/// <param name="Space"><c>-space</c>: space threshold in megabytes for the linear-space recursion (default 10.0).</param>
/// <param name="Shuffle"><c>-shuffle</c>: number of shuffled comparisons (the program's shuffle does not change the sequence, see <see cref="EstEngine"/>; any value above 0 raises <c>-minscore</c> to the best score + 1).</param>
/// <param name="Seed"><c>-seed</c>: random number seed (no effect on the output).</param>
/// <param name="Align"><c>-align</c>: show the alignment (default false).</param>
/// <param name="Width"><c>-width</c>: alignment width (default 50; 0 or less never ends in the program).</param>
/// <param name="SequenceQualifiers">The sequence qualifiers as typed after the files (<c>-sbegin1 5 -send2 90 -sreverse1</c>; input 1 = the ESTs, 2 = the genomic sequence).</param>
public sealed record EmbossEst2genomeRequest(string EstFasta, string GenomeFasta, int Match = 1, int Mismatch = 1, int GapPenalty = 2, int IntronPenalty = 40, int SplicePenalty = 20, int MinScore = 30,
    bool Reverse = false, bool UseSplice = true, string Mode = "both", bool Best = true, float Space = 10.0f, int Shuffle = 0, int Seed = 20825, bool Align = false, int Width = 50, string? SequenceQualifiers = null);

/// <summary>Output of <see cref="EmbossEst2genomePrograms.Run"/>.</summary>
/// <param name="Text">The <c>-outfile</c> (the Note / Exon / Intron / Span / Segment lines and, with <c>-align</c>, the alignment).</param>
/// <param name="Messages">What the program writes to standard error (the space warning).</param>
public sealed record EmbossEst2genomeOutput(string Text, string Messages);

/// <summary>The EMBOSS 6.6.0 program <c>est2genome</c> (<c>emboss/est2genome.c</c> over <c>nucleus/embest.c</c>; B09 A2-28d): the spliced alignment of ESTs to a genomic sequence on the forward and the reverse strand with the donor / acceptor splice site model.</summary>
public static class EmbossEst2genomePrograms
{
    /// <summary>Runs <c>est2genome</c>.</summary>
    /// <exception cref="ArgumentException">An input has no sequence, the mode is not one of both / forward / reverse, or <see cref="EmbossEst2genomeRequest.Width"/> is not positive with <see cref="EmbossEst2genomeRequest.Align"/> (the program never ends).</exception>
    /// <exception cref="InvalidOperationException">The program crashes: <c>-nobest</c> with <c>-mode forward</c> / <c>reverse</c> (a null pointer is dereferenced), or an inconsistent table of remembered intron starts (<c>ajFatal</c>).</exception>
    public static EmbossEst2genomeOutput Run(EmbossEst2genomeRequest r)
    {
        ArgumentNullException.ThrowIfNull(r);
        ArgumentNullException.ThrowIfNull(r.EstFasta);
        ArgumentNullException.ThrowIfNull(r.GenomeFasta);
        string mode = ResolveMode(r.Mode);
        if (r.Align && r.Width <= 0)
            throw new ArgumentException("est2genome -align with a width of 0 or less never ends (the alignment loop adds the width).", nameof(r));
        var qualifiers = EmbossSequenceQualifiers.ParseInputs("est2genome", r.SequenceQualifiers);
        var ests = EmbossFastaReader.Read(r.EstFasta);
        var genomes = EmbossFastaReader.Read(r.GenomeFasta);
        if (ests.Count == 0 || genomes.Count == 0)
            throw new ArgumentException("Unable to read sequence.", nameof(r));
        var set = ests.Select(e => new EmbossSequence(e.Name, e.Text, e.Description, qualifiers[0])).ToList();
        var (_, estViews) = EmbossSetInput.ReadAll(set, EmbossSequenceType.Dna, EmbossInputKind.Stream);
        var genomeView = EmbossInputSequence.Read(genomes[0].Text, qualifiers[1], EmbossSequenceType.Dna);

        int mismatch = Math.Abs(r.Mismatch), gap = Math.Abs(r.GapPenalty), intron = Math.Abs(r.IntronPenalty), splice = Math.Abs(r.SplicePenalty);
        float megabytes = r.Space;
        var messages = new StringBuilder();
        string genome = genomeView.Residues;
        if (megabytes < genome.Length * 1.5e-6)
        {
            float adjusted = (float)1.5e-6 * genome.Length;
            messages.Append("Warning: increasing space from ").Append(CFormat.Fixed(megabytes, 3)).Append(" to ").Append(CFormat.Fixed(adjusted, 3)).Append(" Mb\n");
            megabytes = adjusted;
        }

        var engine = new EstEngine(r.Match, mismatch, gap);
        var spliceSites = r.UseSplice ? EstEngine.FindSpliceSites(genome, true) : null;
        var reversedSpliceSites = mode == "both" && r.UseSplice ? EstEngine.FindSpliceSites(genome, false) : null;
        string genomeName = genomes[0].Name;
        int genomeOffset = genomeView.Offset;
        int minScore = r.MinScore;
        int maxScore = 0;
        var o = new StringBuilder();

        for (int n = 0; n < estViews.Count; n++)
        {
            string est = estViews[n].Residues, estName = ests[n].Name, estDesc = ests[n].Description;
            int estOffset = estViews[n].Offset;
            if (r.Shuffle > 0)
            {
                // embEstShuffleSeq leaves the sequence unchanged: every shuffled comparison is the local alignment of the EST itself
                int score = engine.AlignNonRecursive(est, genome, gap, intron, splice, spliceSites, false, false, EstEngine.Diagonal).Score;
                for (int k = 0; k < r.Shuffle; k++)
                    if (score > maxScore)
                        maxScore = score;
                minScore = maxScore + 1;
            }

            EstAlign? fge = null, rge = null, bge = null;
            string reversedEst = string.Empty;
            int reversedOffset = 0;
            if (mode != "reverse")
                fge = engine.AlignLinearSpace(est, genome, r.Match, mismatch, gap, intron, splice, spliceSites, megabytes);
            if (mode != "forward")
            {
                reversedEst = EmbossInputSequence.ReverseComplement(est);
                reversedOffset = estViews[n].Offend;                    // ajSeqReverseForce swaps Offset and Offend
                rge = engine.AlignLinearSpace(reversedEst, genome, r.Match, mismatch, gap, intron, splice, spliceSites, megabytes);
            }

            void Output(string estSeq, int estOff, EstAlign? ge, bool reverse)
            {
                if (ge == null)
                    throw new InvalidOperationException("est2genome: segmentation fault (the alignment of the other strand does not exist: -nobest with a single strand mode).");
                if (ge.Score >= minScore)
                {
                    engine.OutBlastStyle(o, genome, genomeName, genomeOffset, estSeq, estName, estDesc, estOff, ge, gap, intron, splice, true, reverse);
                    o.Append('\n');
                    engine.OutBlastStyle(o, genome, genomeName, genomeOffset, estSeq, estName, estDesc, estOff, ge, gap, intron, splice, false, reverse);
                    if (r.Align)
                    {
                        o.Append("\n\n").Append(genomeName).Append(" vs ").Append(estName).Append(":\n");
                        engine.PrintAlign(o, genome, genomeName, estSeq, estName, ge, r.Width);
                    }
                }
            }

            if (mode == "both")
            {
                if (fge!.Score > rge!.Score)
                {
                    bge = engine.AlignLinearSpace(est, genome, r.Match, mismatch, gap, intron, splice, reversedSpliceSites, megabytes);
                    if (bge.Score > fge.Score)
                    {
                        o.Append("Note Best alignment is between forward est and forward genome, but splice sites imply REVERSED GENE\n");
                        Output(est, estOffset, bge, false);
                        if (!r.Best)
                            Output(est, estOffset, fge, false);
                    }
                    else
                    {
                        o.Append("Note Best alignment is between forward est and forward genome, and splice sites imply forward gene\n");
                        Output(est, estOffset, fge, false);
                        if (!r.Best)
                            Output(est, estOffset, bge, false);
                    }
                }
                else
                {
                    bge = engine.AlignLinearSpace(reversedEst, genome, r.Match, mismatch, gap, intron, splice, reversedSpliceSites, megabytes);
                    if (bge.Score > rge.Score)
                    {
                        o.Append("Note Best alignment is between reversed est and forward genome, but splice sites imply REVERSED GENE\n");
                        Output(reversedEst, reversedOffset, bge, r.Reverse);
                        if (!r.Best)
                            Output(reversedEst, reversedOffset, rge, r.Reverse);
                    }
                    else
                    {
                        o.Append("Note Best alignment is between reversed est and forward genome, and splice sites imply forward gene\n");
                        Output(reversedEst, reversedOffset, rge, r.Reverse);
                        if (!r.Best)
                            Output(reversedEst, reversedOffset, bge, r.Reverse);
                    }
                }
            }
            else if (mode == "forward")
            {
                o.Append("Note requested forward est and forward genome\n");
                Output(est, estOffset, fge, false);
                if (!r.Best)
                    Output(est, estOffset, bge, false);
            }
            else
            {
                o.Append("Note requested reversed est and forward genome\n");
                Output(reversedEst, reversedOffset, rge, r.Reverse);
                if (!r.Best)
                    Output(reversedEst, reversedOffset, bge, r.Reverse);
            }
        }

        return new EmbossEst2genomeOutput(o.ToString(), messages.ToString());
    }

    private static string ResolveMode(string? mode)
    {
        string m = (mode ?? string.Empty).Trim().ToLowerInvariant();
        string[] menu = ["both", "forward", "reverse"];
        var hits = m.Length == 0 ? [] : menu.Where(x => x.StartsWith(m, StringComparison.Ordinal)).ToArray();
        if (hits.Length != 1)
            throw new ArgumentException("est2genome: '" + mode + "' is not a valid menu option (both, forward, reverse).", nameof(mode));
        return hits[0];
    }
}
