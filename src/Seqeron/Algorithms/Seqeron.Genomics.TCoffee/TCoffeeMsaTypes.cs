// SPDX-License-Identifier: GPL-2.0-or-later
// Derived from T-Coffee 13.41 (https://github.com/cbcrg/tcoffee), Copyright (C) Centre National de la Recherche Scientifique (CNRS), Cedric Notredame; GPL-2.0-or-later. ProbCons 1.12 (public domain) is hosted here too.
// Part of Seqeron.Genomics.TCoffee, licensed under the GNU General Public License, version 2 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
namespace Seqeron.Genomics.TCoffee;

using Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Infrastructure;

/// <summary>Sequence redundancy weights of the consistency-based aligner (T-Coffee <c>-seq_weight</c>).</summary>
public enum MsaSequenceWeights
{
    /// <summary>No weights (every sequence counts 1).</summary>
    None = 0,

    /// <summary>
    /// <c>-seq_weight t_coffee</c> (<c>compute_t_coffee_weight</c>): w = 2 / (1 + sum over the other sequences of (10 x identity)^3), rescaled to a mean of 1; an extended library weight of
    /// the pair (Si.p, Sj.q) is multiplied by min(w_i, w_j) (<c>re_weight_constraint_list</c>, applied after the extension). The identities are those of the program's default k-tuple
    /// distance matrix; the residue type is decided by the program's rule (more than 98 % nucleotide letters = DNA).
    /// </summary>
    TCoffee = 1,
}

/// <summary>Progressive stage of the consistency-based (T-Coffee) multiple aligner.</summary>
public enum MsaConsistencyEngine
{
    /// <summary>The paper's objective: extended library weight (direct + sum of triplet minima) summed over every cross pair of two columns, zero gap penalty, full DP (default).</summary>
    PaperExtension = 0,

    /// <summary>
    /// The T-Coffee 13.41 program's progressive stage (<c>SequenceAligner.MultipleAlignTCoffee</c>: library relaxation, sparse linked-pair-wise DP).
    /// Identical to the program for the same primary library and guide tree; the library is this repository's (global + local alignments, weight 10 x percent identity) and the tree the configured guide tree.
    /// </summary>
    TCoffeeLinked = 1,
}

/// <summary>Options of the consistency-based aligner (<c>MultipleAlignConsistency</c>): the options of the progressive family plus the T-Coffee library settings.</summary>
public sealed record ConsistencyMsaOptions : ProgressiveMsaOptions
{
    /// <summary>
    /// Progressive stage of <c>MultipleAlignConsistency</c> (default <see cref="MsaConsistencyEngine.PaperExtension"/>, the published formula).
    /// </summary>
    public MsaConsistencyEngine ConsistencyEngine { get; init; } = MsaConsistencyEngine.PaperExtension;

    /// <summary>
    /// <c>MultipleAlignConsistency</c> only: the number of non-intersecting local alignments per pair that enter the primary library (default 1 = the single best Smith-Waterman alignment of the
    /// published formula). With N &gt; 1 the Lalign library of the paper is built as T-Coffee's <c>lalign_id_pair</c> does: the Waterman-Eggert top-N alignments
    /// (<c>SequenceAligner.LocalAlignTopN</c>, the program's <c>-lalign_n_top</c> default is 10), each weighted by its own percent identity and
    /// signal-added with the global library.
    /// </summary>
    public int PrimaryLocalAlignments { get; init; } = 1;

    /// <summary>
    /// <c>MultipleAlignConsistency</c> only: redundancy weights of the sequences (T-Coffee <c>-seq_weight</c>; default <see cref="MsaSequenceWeights.None"/>).
    /// </summary>
    public MsaSequenceWeights SequenceWeights { get; init; } = MsaSequenceWeights.None;
}
