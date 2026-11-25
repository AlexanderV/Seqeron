// SPDX-License-Identifier: LGPL-3.0-or-later
// Derived from ClustalW 2.1 (http://www.clustal.org/), Copyright (C) Des Higgins, Julie Thompson, Toby Gibson; LGPL-3.0-or-later.
// Part of Seqeron.Genomics.ClustalW, licensed under the GNU Lesser General Public License, version 3 or (at your option) any later version (see LICENSE in the project directory).
// This file is NOT covered by the MIT licence of the other Seqeron assemblies.
using Seqeron.Genomics.Infrastructure;

namespace Seqeron.Genomics.ClustalW;

using Seqeron.Genomics.Alignment;
using static Seqeron.Genomics.Alignment.SequenceAligner;

public static partial class ClustalWAligner
{
    // ClustalW 2.1 command-level entry points (audit A2-14b, B09 F104): the alignment commands with NAMED sequences, so that guide-tree files can be written (-NEWTREE, -NEWTREE1/2) and read
    // (-USETREE, -USETREE1/2) and the rows can be listed in the program's -OUTORDER=ALIGNED order. The alignment itself is the port of F56 / F62 / F103 (ClustalWMsaState).

    /// <summary>
    /// <c>clustalw -ALIGN</c> / <c>-NEWTREE</c> / <c>-USETREE</c> for DNA sequences: the multiple alignment (<c>Clustal::align</c>; with <see cref="ClustalWRunOptions.UseTree"/> the guide tree comes from the text
    /// instead of the pairwise stage, <c>Clustal::doAlignUseOldTree</c>), the guide-tree file the program writes, and the order of the rows.
    /// </summary>
    /// <param name="sequences">The named sequences (at least two; the program writes nothing for one).</param>
    /// <param name="options">ClustalW DNA parameters (default: ClustalW's DNA defaults).</param>
    /// <param name="run">Guide-tree and row-order parameters.</param>
    /// <exception cref="ClustalWTreeFileException">The user tree is refused by the program.</exception>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static ClustalWRunResult RunClustalW(IReadOnlyList<ClustalWSequence> sequences, ClustalWOptions? options = null, ClustalWRunOptions? run = null)
    {
        options ??= new ClustalWOptions();
        run ??= new ClustalWRunOptions();
        ValidateClustalWOptions(options);
        var names = ValidateClustalWSequences(sequences, minimum: 2, nameof(sequences));
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var raw = sequences.Select(s => s.Residues).ToArray();
        var codes = raw.Select(r => EncodeClustalWDna(r, "sequences")).ToArray();
        RequireNonEmpty(codes);
        var state = new ClustalWMsaState(codes, CwConfig.FromDna(options)) { Names = names, UserTree = run.UseTree, GuideTreesOnly = run.GuideTreeOnly && run.UseTree == null };
        return FinishMain(state, names, run, () => ClustalWPairwiseDistanceMatrix(raw, options), protein: false);
    }

    /// <summary>The PROTEIN counterpart of <see cref="RunClustalW(IReadOnlyList{ClustalWSequence}, ClustalWOptions?, ClustalWRunOptions?)"/> (<c>clustalw -TYPE=PROTEIN</c>).</summary>
    public static ClustalWRunResult RunClustalW(IReadOnlyList<ClustalWSequence> sequences, ClustalWProteinOptions options, ClustalWRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(options);
        run ??= new ClustalWRunOptions();
        ValidateClustalWProteinMultipleOptions(options);
        var names = ValidateClustalWSequences(sequences, minimum: 2, nameof(sequences));
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var raw = sequences.Select(s => s.Residues).ToArray();
        var codes = raw.Select(r => EncodeClustalWProtein(r, "sequences")).ToArray();
        RequireNonEmpty(codes);
        var state = new ClustalWMsaState(codes, ProteinConfig(options)) { Names = names, UserTree = run.UseTree, GuideTreesOnly = run.GuideTreeOnly && run.UseTree == null };
        return FinishMain(state, names, run, () => ClustalWProteinDistanceMatrix(raw.Select(r => r.ToUpperInvariant()).ToArray(), options), protein: true);
    }

    /// <summary>
    /// <c>clustalw -PROFILE1=a -PROFILE2=b -PROFILE</c> for DNA with the guide-tree files of the profiles (<see cref="ClustalWRunOptions.UseTree1"/> / <see cref="ClustalWRunOptions.UseTree2"/>,
    /// <see cref="ClustalWRunResult.GuideTree1"/> / <see cref="ClustalWRunResult.GuideTree2"/>); the rows are profile 1 then profile 2.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static ClustalWRunResult RunClustalWProfile(IReadOnlyList<ClustalWSequence> profile1, IReadOnlyList<ClustalWSequence> profile2, ClustalWOptions? options = null, ClustalWRunOptions? run = null)
    {
        options ??= new ClustalWOptions();
        run ??= new ClustalWRunOptions();
        ValidateClustalWOptions(options);
        ValidateClustalWProfileOptions(options.Iteration);
        var (names, all) = ValidateClustalWProfiles(profile1, profile2);
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var codes = all.Select(s => EncodeClustalWDna(s.Residues, "profiles", allowGaps: true)).ToArray();
        var state = new ClustalWMsaState(codes, CwConfig.FromDna(options)) { Names = names, UserTree1 = run.UseTree1, UserTree2 = run.UseTree2, GuideTreesOnly = run.GuideTreeOnly, StructureOptions = run.Structure, Structure1Input = run.ProfileStructure1, Structure2Input = run.ProfileStructure2 };
        return FinishProfile(state, names, run, profile1.Count, protein: false);
    }

    /// <summary>The PROTEIN counterpart of <see cref="RunClustalWProfile(IReadOnlyList{ClustalWSequence}, IReadOnlyList{ClustalWSequence}, ClustalWOptions?, ClustalWRunOptions?)"/>.</summary>
    public static ClustalWRunResult RunClustalWProfile(IReadOnlyList<ClustalWSequence> profile1, IReadOnlyList<ClustalWSequence> profile2, ClustalWProteinOptions options, ClustalWRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(options);
        run ??= new ClustalWRunOptions();
        ValidateClustalWProteinMultipleOptions(options);
        ValidateClustalWProfileOptions(options.Iteration);
        var (names, all) = ValidateClustalWProfiles(profile1, profile2);
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var codes = all.Select(s => EncodeClustalWProtein(s.Residues, "profiles", allowGaps: true)).ToArray();
        var state = new ClustalWMsaState(codes, ProteinConfig(options)) { Names = names, UserTree1 = run.UseTree1, UserTree2 = run.UseTree2, GuideTreesOnly = run.GuideTreeOnly, StructureOptions = run.Structure, Structure1Input = run.ProfileStructure1, Structure2Input = run.ProfileStructure2 };
        return FinishProfile(state, names, run, profile1.Count, protein: true);
    }

    /// <summary>
    /// <c>clustalw -PROFILE1=a -PROFILE2=b -SEQUENCES</c> for DNA with the guide tree of all sequences (<see cref="ClustalWRunOptions.UseTree"/>: <c>getWeightsFromGuideTree</c>, no pairwise stage;
    /// <see cref="ClustalWRunResult.GuideTree"/> otherwise); the rows are the profile, then the sequences.
    /// </summary>
    [System.Runtime.CompilerServices.OverloadResolutionPriority(1)]
    public static ClustalWRunResult RunClustalWSequences(IReadOnlyList<ClustalWSequence> profile, IReadOnlyList<ClustalWSequence> sequences, ClustalWOptions? options = null, ClustalWRunOptions? run = null)
    {
        options ??= new ClustalWOptions();
        run ??= new ClustalWRunOptions();
        ValidateClustalWOptions(options);
        ValidateClustalWProfileOptions(options.Iteration);
        var (names, all) = ValidateClustalWProfiles(profile, sequences);
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var codes = all.Select(s => EncodeClustalWDna(s.Residues, "profiles", allowGaps: true)).ToArray();
        var state = new ClustalWMsaState(codes, CwConfig.FromDna(options)) { Names = names, UserTree = run.UseTree, GuideTreesOnly = run.GuideTreeOnly, StructureOptions = run.Structure, Structure1Input = run.ProfileStructure1, Structure2Input = run.ProfileStructure2 };
        return FinishSequences(state, names, run, profile.Count, ClustalWSequencesPairwise(options, codes), protein: false);
    }

    /// <summary>The PROTEIN counterpart of <see cref="RunClustalWSequences(IReadOnlyList{ClustalWSequence}, IReadOnlyList{ClustalWSequence}, ClustalWOptions?, ClustalWRunOptions?)"/>.</summary>
    public static ClustalWRunResult RunClustalWSequences(IReadOnlyList<ClustalWSequence> profile, IReadOnlyList<ClustalWSequence> sequences, ClustalWProteinOptions options, ClustalWRunOptions? run = null)
    {
        ArgumentNullException.ThrowIfNull(options);
        run ??= new ClustalWRunOptions();
        ValidateClustalWProteinMultipleOptions(options);
        ValidateClustalWProfileOptions(options.Iteration);
        var (names, all) = ValidateClustalWProfiles(profile, sequences);
        using var cancelScope = AlignmentCancellation.Enter(options.CancellationToken);
        var codes = all.Select(s => EncodeClustalWProtein(s.Residues, "profiles", allowGaps: true)).ToArray();
        var state = new ClustalWMsaState(codes, ProteinConfig(options)) { Names = names, UserTree = run.UseTree, GuideTreesOnly = run.GuideTreeOnly, StructureOptions = run.Structure, Structure1Input = run.ProfileStructure1, Structure2Input = run.ProfileStructure2 };
        return FinishSequences(state, names, run, profile.Count, ClustalWSequencesPairwise(options, codes), protein: true);
    }

    private static CwConfig ProteinConfig(ClustalWProteinOptions options)
    {
        string hyd = new string(options.HydrophilicResidues.Take(ClustalWMaxHydrophilicResidues).TakeWhile(char.IsAsciiLetter).Select(char.ToUpperInvariant).ToArray());
        return CwConfig.FromProtein(options) with { HydrophilicResidues = hyd.Length > 0 ? hyd : "GPSNDQEKR" };
    }

    private static void RequireNonEmpty(int[][] codes)
    {
        for (int i = 0; i < codes.Length; i++)
        {
            if (codes[i].Length == 0)
                throw new ArgumentException(EmptyMemberMessage("ClustalW", i), nameof(codes));
        }
    }

    // FileReader / Alignment: names are the first word of the header (<= MAXNAMES), all different (exact comparison), no empty sequences
    private static List<string> ValidateClustalWSequences(IReadOnlyList<ClustalWSequence> sequences, int minimum, string paramName)
    {
        ArgumentNullException.ThrowIfNull(sequences);
        var names = new List<string>(sequences.Count);
        foreach (var s in sequences)
        {
            ArgumentNullException.ThrowIfNull(s);
            ArgumentNullException.ThrowIfNull(s.Name);
            ArgumentNullException.ThrowIfNull(s.Residues);
            if (s.Name.Length == 0 || s.Name.Length > CwTreeMaxNames || s.Name.Any(char.IsWhiteSpace))
                throw new ArgumentException($"Sequence name '{s.Name}': a ClustalW name is the first word of the FASTA header (1 to {CwTreeMaxNames} characters, no blanks).", paramName);
            names.Add(s.Name);
        }

        if (names.Distinct(StringComparer.Ordinal).Count() != names.Count)
            throw new ArgumentException("Multiple sequences found with same name (ClustalW: not all sequence names are different).", paramName);
        if (sequences.Count < minimum)
            throw new ArgumentException("ClustalW needs at least two sequences for a multiple alignment (the program prints 'Only 1 sequence, cannot do multiple alignment' and writes nothing).", paramName);
        return names;
    }

    private static (List<string> Names, List<ClustalWSequence> All) ValidateClustalWProfiles(IReadOnlyList<ClustalWSequence> profile1, IReadOnlyList<ClustalWSequence> profile2)
    {
        ArgumentNullException.ThrowIfNull(profile1);
        ArgumentNullException.ThrowIfNull(profile2);
        if (profile1.Count == 0 || profile2.Count == 0)
            throw new ArgumentException("A profile alignment needs at least one sequence in each profile (ClustalW: 'No sequences in memory').", profile1.Count == 0 ? nameof(profile1) : nameof(profile2));
        var all = profile1.Concat(profile2).ToList();
        var names = ValidateClustalWSequences(all, minimum: 2, "profiles");
        if (all.Any(s => s.Residues.Replace("-", "").Length == 0))
            throw new ArgumentException("Empty sequences found.", nameof(profile1));
        return (names, all);
    }

    private static ClustalWRunResult FinishMain(ClustalWMsaState state, List<string> names, ClustalWRunOptions run, Func<double[,]> distances, bool protein)
    {
        int n = names.Count;
        // Clustal::doAlignUseOldTree runs no pairwise stage
        var dist = run.UseTree != null ? new double[n, n] : distances();
        string[] rows = state.Run(dist);
        if (state.StoppedAfterTrees)
            return new ClustalWRunResult(null, state.GuideTreeText, null, null);

        return new ClustalWRunResult(BuildAlignment(state, names, rows, run, protein), state.GuideTreeText, null, null);
    }

    private static ClustalWRunResult FinishProfile(ClustalWMsaState state, List<string> names, ClustalWRunOptions run, int profile1Count, bool protein)
    {
        string[] rows = state.RunProfileAlignment(profile1Count);
        if (state.StoppedAfterTrees)
            return new ClustalWRunResult(null, null, state.GuideTree1Text, state.GuideTree2Text);

        return new ClustalWRunResult(BuildAlignment(state, names, rows, run, protein, profile1Count), null, state.GuideTree1Text, state.GuideTree2Text);
    }

    private static ClustalWRunResult FinishSequences(ClustalWMsaState state, List<string> names, ClustalWRunOptions run, int profile1Count, Func<int[], int[], double> pairwise, bool protein)
    {
        string[] rows = state.RunSequencesToProfile(profile1Count, pairwise);
        if (state.StoppedAfterTrees)
            return new ClustalWRunResult(null, state.GuideTreeText, null, null);

        return new ClustalWRunResult(BuildAlignment(state, names, rows, run, protein, profile1Count), state.GuideTreeText, null, null);
    }

    private static ClustalWAlignment BuildAlignment(ClustalWMsaState state, List<string> names, string[] rows, ClustalWRunOptions run, bool protein, int profile1Count = 0)
    {
        int n = names.Count;
        IReadOnlyList<int> order = run.OutputOrder == ClustalWOutputOrder.Input ? Enumerable.Range(0, n).ToList() : state.AlignedOrder.ToList();
        return new ClustalWAlignment(names, rows, order, protein, state.FinalSequenceWeights.Take(n).ToList())
        {
            Structure1 = state.FinalStructure1,
            Structure2 = state.FinalStructure2,
            Profile1Count = profile1Count,
        };
    }
}
