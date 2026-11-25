// SPDX-License-Identifier: GPL-2.0-or-later
// Tests of the MCP tools of the port of EMBOSS 6.6.0 (the GNU General Public License, version 2 or later).
// Part of Seqeron.Mcp.EmbossPairwiseAlignment, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using NUnit.Framework;
using Seqeron.Mcp.Alignment.Tools;
using Seqeron.Mcp.EmbossPairwiseAlignment.Tools;

namespace Seqeron.Mcp.EmbossPairwiseAlignment.Tests;

/// <summary>
/// R14-7 (B09 F300): the <c>significance</c> option of <c>local_align_affine</c> / <c>local_align</c> / <c>water_align</c>: bit score, E-value, P-value and the Karlin-Altschul parameters of the score. Truth: the installed BLAST+
/// 2.12.0 on the same pairs (protein: <c>blastp -comp_based_stats 0 -seg no -subject</c>, 150 x 160, score 256: bit score 103.219, E 1.24008e-33; nucleotide: <c>blastn -task blastn -dust no</c>, 430 x 880, score 145:
/// bit score 132.031, E 1.21e-34), the same numbers the library test locks over 2 283 HSPs.
/// </summary>
[TestFixture]
public class LocalAlignSignificanceTests
{

    [Test]
    public void WaterAlign_Protein_IntegerGapCosts_UseExistenceMinusExtension()
    {
        // EMBOSS gapopen 12 / gapextend 1 = BLAST existence 11 extension 1 (a gap of length k costs 12 + (k - 1)): the same row of the BLOSUM62 table
        var r = EmbossPairwiseAlignmentTools.WaterAlign(ProteinQuery, ProteinSubject, matrix: "EBLOSUM62", gapOpen: 12, gapExtend: 1, significance: true);
        Assert.That(r.Significance, Is.Not.Null);
        Assert.That((r.Significance!.Lambda, r.Significance.K), Is.EqualTo((0.267, 0.041)));
        Assert.That(r.Significance.Method, Is.EqualTo("SpougeFiniteSizeCorrection"));
    }

    private const string ProteinQuery = "MDLNLSMGSSGCKANPAVQNPWARCGYHEIRNTNTKREWNYKNTARICGENSIPVEWIAHQKGPGDFYYRFFAAHHGGLMHVHGHPLANQGFKDMLYWAYLCQDHMTYNFMKVDLMLGDFLSGCDYVPCIYNKRQFCCSMHFWFQEGQNF";

    private const string ProteinSubject = "CQLFRYGTRSMSKLSPFEPVGSMGDSKTVTNDNWCLNVKSKLMGWASVKMKRLTNNKNQNGRNMTFTGHNYKNTAWICGENSIWVEGISHFGGPGDFHNPFFAAHHGGLMHVHQHPLYNQFFKDELYWAYWPQQISLSPPGYYKLSKQAMLSLFSAEYRI";
}

/// <summary>R14-7 (B09 F300): <c>local_align_zscore</c>, the shuffle Z-score of a local score (EMBOSS prss / FASTA idea; seeded, so reproducible).</summary>
[TestFixture]
public class LocalAlignZScoreTests
{
    private const string Core = "ACGTTGCAAGGCTTAACCGGTTAACGTACGGTACCATGCA";

    [Test]
    public void ZScore_OfASharedMotif_IsLarge_AndReproducible()
    {
        var a = AlignmentTools.LocalAlignZScore("TTTT" + Core + "GGGG", "CCCCAAAA" + Core + "TTTTCCCC", shuffles: 200, seed: 7);
        var b = AlignmentTools.LocalAlignZScore("TTTT" + Core + "GGGG", "CCCCAAAA" + Core + "TTTTCCCC", shuffles: 200, seed: 7);
        Assert.That(a, Is.EqualTo(b));
        Assert.That(a.Score, Is.EqualTo(AlignmentTools.LocalAlign("TTTT" + Core + "GGGG", "CCCCAAAA" + Core + "TTTTCCCC").Score), "the score of the pair is the local_align score");
        Assert.That(a.ZScore, Is.GreaterThan(5));
        Assert.That(a.StandardDeviation, Is.GreaterThan(0));
        Assert.That((a.Shuffles, a.Seed), Is.EqualTo((200, 7)));
        Assert.That(AlignmentTools.LocalAlignZScore("TTTT" + Core + "GGGG", "CCCCAAAA" + Core + "TTTTCCCC", shuffles: 200, seed: 8).Mean, Is.Not.EqualTo(a.Mean), "another seed, other shuffles");
    }

    [Test]
    public void ZScore_Validation()
    {
        Assert.Throws<ArgumentException>(() => AlignmentTools.LocalAlignZScore(null!, "ACGT"));
        Assert.Throws<ArgumentOutOfRangeException>(() => AlignmentTools.LocalAlignZScore("ACGT", "ACGT", shuffles: 1));
        using var cts = new CancellationTokenSource();
        cts.Cancel();
        Assert.Throws<OperationCanceledException>(() => AlignmentTools.LocalAlignZScore("ACGTACGT", "ACGTACGT", cancellationToken: cts.Token));
    }
}
