// SPDX-License-Identifier: AGPL-3.0-only
// This file is part of Seqeron.Mcp.OverrepresentedWords, licensed under the terms stated in LICENSE (this directory).
using System.ComponentModel;
using System.Text;
using ModelContextProtocol.Server;
using Seqeron.Mcp.Analysis.Tools;
using static Seqeron.Mcp.Analysis.Tools.AnalysisTools;

namespace Seqeron.Mcp.OverrepresentedWords.Tools;

public partial class OverrepresentedWordsTools
{
    [McpServerTool(Name = "oligo_analysis", Title = "Motifs — RSAT oligo-analysis Significance", ReadOnly = true)]
    [Description("RSAT oligo-analysis over-representation of the k-mers of one DNA sequence (or a set with extraSequences): occurrences, expected frequency under a background model (equiprobable, input Bernoulli, given Bernoulli, input Markov order m, Markov table, or lexicon), observed/expected ratio and binomial significance occ_P / occ_E / occ_sig; single strand or reverse-complement pairs (-2str); overlapping or non-overlapping (-noov) counts. RSAT options: z-scores with the expected variance and overlap coefficient (-return zscore), pseudo-frequency on expected frequencies (-pseudo), degenerate words with one N or one IUPAC code (-oneN / -onedeg), calibration tables with negative-binomial / Poisson P-values (-calibN / -calib1); protein or free-text input with sequenceType (-seqtype prot|other).")]
    public static OligoAnalysisResultDto OligoAnalysis(
        [Description("DNA sequence.")] string sequence,
        [Description("Oligonucleotide length k (>= 1, default 6).")] int k = 6,
        [Description("Occurrence threshold (RSAT -lth occ, default 2).")] int minCount = 2,
        [Description("Background model: 'input' (default, Bernoulli from input), 'equiprobable', 'bernoulli' (needs residueFrequencies), 'markov' (order markovOrder from input), 'markov_table' (needs oligoFrequencies) or 'lexicon' (RSAT -lexicon, k >= 2).")] string background = "input",
        [Description("Markov order for background 'markov' (>= 0, <= k-2 when > 0; default 1).")] int markovOrder = 1,
        [Description("Residue probabilities A,C,G,T for background 'bernoulli'.")] double[]? residueFrequencies = null,
        [Description("(m+1)-mer frequency table for background 'markov_table' (RSAT -bgfile oligos).")] Dictionary<string, double>? oligoFrequencies = null,
        [Description("Pseudo-frequency in [0,1] for background 'markov_table' (default 0.01).")] double pseudoFrequency = 0.01,
        [Description("The 'markov_table' frequencies are strand-insensitive pair frequencies (default false).")] bool strandInsensitive = false,
        [Description("Strands: 'single' (-1str, default) or 'both' (-2str, reverse-complement pairs).")] string strands = "single",
        [Description("Count overlapping occurrences (-ovlp, default true); false = -noov.")] bool countOverlapping = true,
        [Description("Additional DNA sequences analysed together with 'sequence' (RSAT multi-sequence input; positions carry sequenceIndices, 0 = 'sequence').")] string[]? extraSequences = null,
        [Description("Also report RSAT z-scores: expectedVariance, overlapCoefficient, zScore (-return zscore). Default false.")] bool zscore = false,
        [Description("RSAT -pseudo: pseudo-frequency psi in [0,1] on the expected frequencies, exp_freq = (1 - psi)·exp_freq + psi / possibleOligos per strand (default 0 = off).")] double expectedFrequencyPseudo = 0,
        [Description("Degenerate words: 'none' (default), 'oneN' (-oneN, one N per word) or 'onedeg' (-onedeg, one IUPAC code R Y W S M K H B V D N per word).")] string degenerate = "none",
        [Description("RSAT calibration table text (-calibN / -calib1; RSAT fit-distribution layout: word [id] mean sd variance ...). When given, P-values use a negative binomial (mean < variance) or Poisson fitted to the table.")] string? calibrationTable = null,
        [Description("Calibration mode: 'set' (-calibN, default) or 'sequence' (-calib1, mean and variance multiplied by the number of sequences).")] string calibrationMode = "set",
        [Description("RSAT -seqtype: 'dna' (default), 'protein' (20-amino-acid alphabet; windows with other residues such as X or * are discarded) or 'other' (any text; alphabet = the observed residues). Protein / other: single strand, no degenerate words or calibration, background input / equiprobable / markov / lexicon; case and white space are folded.")] string sequenceType = "dna")
    {
        var seqType = (sequenceType ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "dna" => global::Seqeron.Genomics.OverrepresentedWords.OligoSequenceType.Dna,
            "protein" or "prot" => global::Seqeron.Genomics.OverrepresentedWords.OligoSequenceType.Protein,
            "other" => global::Seqeron.Genomics.OverrepresentedWords.OligoSequenceType.Other,
            _ => throw new ArgumentException("sequenceType must be 'dna', 'protein' or 'other'.", nameof(sequenceType)),
        };
        if (seqType != global::Seqeron.Genomics.OverrepresentedWords.OligoSequenceType.Dna)
            return OligoAnalysisResidueAlphabet(sequence, k, minCount, background, markovOrder, strands, countOverlapping,
                extraSequences, zscore, expectedFrequencyPseudo, degenerate, calibrationTable, seqType);

        var dna = RequireDna(sequence, nameof(sequence));
        if (k < 1)
            throw new ArgumentOutOfRangeException(nameof(k), "k must be >= 1.");
        var bg = ToOligoBackground(background, markovOrder, residueFrequencies, oligoFrequencies, pseudoFrequency, strandInsensitive);
        var mode = ToOligoStrandMode(strands);
        var degeneracy = ToOligoDegeneracy(degenerate);
        if (!(expectedFrequencyPseudo >= 0 && expectedFrequencyPseudo <= 1))
            throw new ArgumentOutOfRangeException(nameof(expectedFrequencyPseudo), "expectedFrequencyPseudo must be in [0, 1].");

        bool extended = (extraSequences is { Length: > 0 }) || zscore || expectedFrequencyPseudo > 0
                        || degeneracy != global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy.None || calibrationTable is not null
                        || ReferenceEquals(bg, global::Seqeron.Genomics.OverrepresentedWords.OligoBackgroundModel.Lexicon);
        if (!extended)
        {
            var r = global::Seqeron.Genomics.OverrepresentedWords.OverrepresentedWordAnalyzer.DiscoverMotifs(dna, k, minCount, bg, mode, countOverlapping);
            var items = r.Motifs
                .Select(m => new OligoMotifItem(
                    m.Sequence, m.ReverseComplement, m.Count, m.Positions.ToArray(),
                    m.ExpectedFrequency, m.ExpectedOccurrences, FiniteOrNull(m.Ratio),
                    m.OccurrenceProbability, m.OccurrenceEValue, m.OccurrenceSignificance))
                .ToArray();
            return new OligoAnalysisResultDto(items, r.OligoLength, StrandName(r.Strands), r.CountOverlapping,
                r.TotalOccurrences, r.TestedPatterns, FiniteOrNull(r.PossibleOligos));
        }

        var dnaList = new List<DnaSequence> { dna };
        foreach (var extra in extraSequences ?? Array.Empty<string>())
            dnaList.Add(RequireDna(extra, nameof(extraSequences)));
        global::Seqeron.Genomics.OverrepresentedWords.OligoCalibration? calibration = null;
        if (calibrationTable is not null)
        {
            var calMode = (calibrationMode ?? string.Empty).Trim().ToLowerInvariant() switch
            {
                "set" => global::Seqeron.Genomics.OverrepresentedWords.OligoCalibrationMode.PerSet,
                "sequence" => global::Seqeron.Genomics.OverrepresentedWords.OligoCalibrationMode.PerSequence,
                _ => throw new ArgumentException("calibrationMode must be 'set' or 'sequence'.", nameof(calibrationMode)),
            };
            try
            {
                calibration = global::Seqeron.Genomics.OverrepresentedWords.OligoCalibration.Parse(calibrationTable, calMode);
            }
            catch (FormatException e)
            {
                throw new ArgumentException(e.Message, nameof(calibrationTable), e);
            }
        }

        var report = global::Seqeron.Genomics.OverrepresentedWords.OverrepresentedWordAnalyzer.AnalyzeOligos(dnaList, k,
            new global::Seqeron.Genomics.OverrepresentedWords.OligoAnalysisOptions
            {
                MinCount = minCount,
                Background = bg,
                Strands = mode,
                CountOverlapping = countOverlapping,
                PseudoFrequency = expectedFrequencyPseudo,
                Degeneracy = degeneracy,
                Calibration = calibration,
            });
        // Patterns without a valid expected frequency (RSAT "NA", not tested) are omitted.
        var extendedItems = ToOligoMotifItems(report, zscore, includeVariance: zscore || calibration is not null);
        return new OligoAnalysisResultDto(extendedItems, report.OligoLength, StrandName(report.Strands), report.CountOverlapping,
            report.TotalOccurrences, report.TestedPatterns, FiniteOrNull(report.PossibleOligos))
        {
            SequenceCount = report.SequenceCount,
            Degenerate = DegeneracyName(degeneracy) ?? "none",
        };
    }
    [McpServerTool(Name = "dyad_analysis", Title = "Motifs — RSAT dyad-analysis (Spaced Dyads)", ReadOnly = true)]
    [Description("RSAT dyad-analysis: over-represented spaced dyads M1 n{s} M2 (two monads of length m separated by s unspecified residues, s in [minSpacing, maxSpacing]) in DNA sequences — any dyad, direct repeats (dr), inverted repeats (ir) or both (rep); expected frequency from the input monad frequencies (or a monad / dyad frequency table), expected occurrences, z-score, and binomial significance occ_P / occ_E / occ_sig over the T_s possible positions. RSAT defaults: m = 3, spacing 0-20, both strands, non-overlapping counts.")]
    public static DyadAnalysisResultDto DyadAnalysis(
        [Description("DNA sequences.")] string[] sequences,
        [Description("Monad length m (RSAT -l, >= 1, default 3).")] int monadLength = 3,
        [Description("Minimal spacing (>= 0, default 0).")] int minSpacing = 0,
        [Description("Maximal spacing (>= minSpacing, default 20).")] int maxSpacing = 20,
        [Description("Dyad type: 'any' (default), 'dr' (direct repeats), 'ir' (inverted repeats) or 'rep' (both).")] string dyadType = "any",
        [Description("Strands: 'both' (-2str, default; reverse-complement pairs) or 'single' (-1str).")] string strands = "both",
        [Description("Count overlapping occurrences (-ovlp); default false = RSAT -noov.")] bool countOverlapping = false,
        [Description("Occurrence threshold (RSAT -lth occ, default 1); <= 0 tests every combination of observed monads (RSAT without a threshold).")] int minCount = 1,
        [Description("Background: 'monads' (default, input monad frequencies), 'monad_table' (needs monadFrequencies, RSAT -mncf) or 'dyad_table' (needs dyadFrequencies, RSAT -expfreq).")] string background = "monads",
        [Description("Monad -> expected frequency, for background 'monad_table'.")] Dictionary<string, double>? monadFrequencies = null,
        [Description("Dyad (RSAT notation, e.g. 'acgn{4}tta') -> expected single-strand frequency, for background 'dyad_table'.")] Dictionary<string, double>? dyadFrequencies = null)
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required", nameof(sequences));
        if (monadLength < 1)
            throw new ArgumentOutOfRangeException(nameof(monadLength), "monadLength must be >= 1.");
        if (minSpacing < 0 || maxSpacing < minSpacing)
            throw new ArgumentOutOfRangeException(nameof(maxSpacing), "Spacings must satisfy 0 <= minSpacing <= maxSpacing.");
        var dnaList = sequences.Select(s => RequireDna(s, nameof(sequences))).ToList();
        var type = (dyadType ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "any" => global::Seqeron.Genomics.OverrepresentedWords.DyadType.Any,
            "dr" => global::Seqeron.Genomics.OverrepresentedWords.DyadType.DirectRepeat,
            "ir" => global::Seqeron.Genomics.OverrepresentedWords.DyadType.InvertedRepeat,
            "rep" => global::Seqeron.Genomics.OverrepresentedWords.DyadType.Repeat,
            _ => throw new ArgumentException("dyadType must be 'any', 'dr', 'ir' or 'rep'.", nameof(dyadType)),
        };
        var bg = (background ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "monads" => global::Seqeron.Genomics.OverrepresentedWords.DyadBackgroundModel.Monads,
            "monad_table" => monadFrequencies is { Count: > 0 }
                ? global::Seqeron.Genomics.OverrepresentedWords.DyadBackgroundModel.MonadFrequencies(monadFrequencies)
                : throw new ArgumentException("Background 'monad_table' requires monadFrequencies.", nameof(monadFrequencies)),
            "dyad_table" => dyadFrequencies is { Count: > 0 }
                ? global::Seqeron.Genomics.OverrepresentedWords.DyadBackgroundModel.DyadFrequencies(dyadFrequencies)
                : throw new ArgumentException("Background 'dyad_table' requires dyadFrequencies.", nameof(dyadFrequencies)),
            _ => throw new ArgumentException("Background must be 'monads', 'monad_table' or 'dyad_table'.", nameof(background)),
        };

        var r = global::Seqeron.Genomics.OverrepresentedWords.OverrepresentedWordAnalyzer.AnalyzeDyads(dnaList, new global::Seqeron.Genomics.OverrepresentedWords.DyadAnalysisOptions
        {
            MonadLength = monadLength,
            MinSpacing = minSpacing,
            MaxSpacing = maxSpacing,
            Type = type,
            Strands = ToOligoStrandMode(strands),
            CountOverlapping = countOverlapping,
            MinCount = minCount,
            Background = bg,
        });
        var items = r.Dyads.Select(d => new DyadItem(
                d.Pattern, d.FirstMonad, d.Spacing, d.SecondMonad, d.ReverseComplement, d.Occurrences, d.Overlaps,
                d.Positions.Select(p => p.SequenceIndex).ToArray(), d.Positions.Select(p => p.Position).ToArray(),
                d.ObservedFrequency, FiniteOrNull(d.ExpectedFrequency), FiniteOrNull(d.ExpectedOccurrences),
                FiniteOrNull(d.ExpectedVariance), d.OverlapCoefficient, FiniteOrNull(d.ZScore), d.Ratio,
                FiniteOrNull(d.OccurrenceProbability), FiniteOrNull(d.OccurrenceEValue), FiniteOrNull(d.OccurrenceSignificance),
                d.IsDirectRepeat, d.IsReversePalindrome, d.ExpectedFromMonads))
            .ToArray();
        var spacings = r.Spacings.Select(s => new DyadSpacingItem(s.Spacing, s.PossiblePositions, s.Occurrences, s.Overlaps)).ToArray();
        return new DyadAnalysisResultDto(items, r.MonadLength, r.MinSpacing, r.MaxSpacing, dyadType!.Trim().ToLowerInvariant(),
            StrandName(r.Strands), r.CountOverlapping, r.SequenceCount, r.MonadOccurrences, spacings, r.TestedPatterns,
            FiniteOrNull(r.PossibleDyads));
    }
    [McpServerTool(Name = "shared_motifs_significance", Title = "Motifs — Shared Motifs Significance (RSAT)", ReadOnly = true)]
    [Description("RSAT oligo-analysis matching-sequence statistics: k-mers (or reverse-complement pairs) present in at least minSequences of the input DNA sequences, with expected matching sequences exp_ms and binomial significance ms_P / ms_E / ms_sig under a background model; optionally degenerate words with one N or one IUPAC code (-oneN / -onedeg).")]
    public static SharedMotifSignificanceResult SharedMotifsSignificance(
        [Description("DNA sequences.")] string[] sequences,
        [Description("Word length k (>= 1, default 6).")] int k = 6,
        [Description("Matching-sequence quorum (RSAT -lth mseq, >= 1, default 2).")] int minSequences = 2,
        [Description("Background model: 'input' (default), 'equiprobable', 'bernoulli' (needs residueFrequencies), 'markov' (order markovOrder) or 'markov_table' (needs oligoFrequencies).")] string background = "input",
        [Description("Markov order for background 'markov' (default 1).")] int markovOrder = 1,
        [Description("Residue probabilities A,C,G,T for background 'bernoulli'.")] double[]? residueFrequencies = null,
        [Description("(m+1)-mer frequency table for background 'markov_table'.")] Dictionary<string, double>? oligoFrequencies = null,
        [Description("Pseudo-frequency in [0,1] for background 'markov_table' (default 0.01).")] double pseudoFrequency = 0.01,
        [Description("The 'markov_table' frequencies are strand-insensitive pair frequencies (default false).")] bool strandInsensitive = false,
        [Description("Strands: 'single' (default) or 'both' (reverse-complement pairs).")] string strands = "single",
        [Description("Degenerate words: 'none' (default), 'oneN' (-oneN, one N per word) or 'onedeg' (-onedeg, one IUPAC code R Y W S M K H B V D N per word); a sequence matches a degenerate word when one of its k-mers matches it.")] string degenerate = "none")
    {
        if (sequences is null || sequences.Length == 0)
            throw new ArgumentException("At least one sequence is required", nameof(sequences));
        if (k < 1)
            throw new ArgumentOutOfRangeException(nameof(k), "k must be >= 1.");
        if (minSequences < 1)
            throw new ArgumentOutOfRangeException(nameof(minSequences), "minSequences must be >= 1.");
        var degeneracy = ToOligoDegeneracy(degenerate);

        var dnaList = sequences.Select(s => RequireDna(s, nameof(sequences))).ToList();
        var bg = ToOligoBackground(background, markovOrder, residueFrequencies, oligoFrequencies, pseudoFrequency, strandInsensitive);
        var mode = ToOligoStrandMode(strands);

        var r = degeneracy == global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy.None
            ? global::Seqeron.Genomics.OverrepresentedWords.OverrepresentedWordAnalyzer.FindSharedMotifs(dnaList, k, minSequences, bg, mode)
            : global::Seqeron.Genomics.OverrepresentedWords.OverrepresentedWordAnalyzer.FindSharedMotifs(dnaList, k, minSequences, bg, mode, 0.0, degeneracy);
        var items = r.Motifs
            .Select(m => new SharedMotifSignificanceItem(
                m.Sequence, m.ReverseComplement, m.SequenceIndices.ToArray(), m.Prevalence,
                m.ExpectedFrequency, m.ExpectedMatchingSequences, m.MatchingSequenceProbability,
                FiniteOrNull(m.MatchingSequenceEValue), FiniteOrNull(m.MatchingSequenceSignificance)))
            .ToArray();
        return new SharedMotifSignificanceResult(items, r.OligoLength, StrandName(r.Strands), r.SequenceCount,
            r.PossiblePositions, FiniteOrNull(r.PossibleOligos))
        {
            Degenerate = DegeneracyName(degeneracy),
        };
    }

    // oligo_analysis with sequenceType 'protein' / 'other' (RSAT -seqtype prot|other): MotifFinder.AnalyzeOligoStrings.
    private static OligoAnalysisResultDto OligoAnalysisResidueAlphabet(
        string sequence, int k, int minCount, string background, int markovOrder, string strands, bool countOverlapping,
        string[]? extraSequences, bool zscore, double expectedFrequencyPseudo, string degenerate, string? calibrationTable,
        global::Seqeron.Genomics.OverrepresentedWords.OligoSequenceType seqType)
    {
        if (string.IsNullOrWhiteSpace(sequence))
            throw new ArgumentException("Sequence cannot be null or empty.", nameof(sequence));
        if (k < 1)
            throw new ArgumentOutOfRangeException(nameof(k), "k must be >= 1.");
        if (!(expectedFrequencyPseudo >= 0 && expectedFrequencyPseudo <= 1))
            throw new ArgumentOutOfRangeException(nameof(expectedFrequencyPseudo), "expectedFrequencyPseudo must be in [0, 1].");
        string type = seqType == global::Seqeron.Genomics.OverrepresentedWords.OligoSequenceType.Protein ? "protein" : "other";
        if (calibrationTable is not null)
            throw new ArgumentException($"Calibration tables are DNA word tables; not available for sequenceType '{type}'.", nameof(calibrationTable));
        if (!string.Equals((degenerate ?? "none").Trim(), "none", StringComparison.OrdinalIgnoreCase))
            throw new ArgumentException($"Degenerate words use IUPAC nucleotide codes; not available for sequenceType '{type}'.", nameof(degenerate));
        var bgName = (background ?? string.Empty).Trim().ToLowerInvariant();
        if (bgName is not ("input" or "equiprobable" or "markov" or "lexicon"))
            throw new ArgumentException($"sequenceType '{type}' supports background 'input', 'equiprobable', 'markov' or 'lexicon'.", nameof(background));
        var bg = ToOligoBackground(bgName, markovOrder, null, null, 0.01, false);
        var mode = ToOligoStrandMode(strands);
        if (mode == global::Seqeron.Genomics.OverrepresentedWords.OligoStrandMode.Both)
            throw new ArgumentException($"sequenceType '{type}' counts a single strand (RSAT -seqtype {type}).", nameof(strands));

        var list = new List<string> { sequence };
        foreach (var extra in extraSequences ?? Array.Empty<string>())
        {
            if (extra is null)
                throw new ArgumentException("extraSequences cannot contain null.", nameof(extraSequences));
            list.Add(extra);
        }

        var report = global::Seqeron.Genomics.OverrepresentedWords.OverrepresentedWordAnalyzer.AnalyzeOligoStrings(list, k,
            new global::Seqeron.Genomics.OverrepresentedWords.OligoAnalysisOptions
            {
                MinCount = minCount,
                Background = bg,
                CountOverlapping = countOverlapping,
                PseudoFrequency = expectedFrequencyPseudo,
                SequenceType = seqType,
            });
        // Single strand only (checked above), so every ReverseComplement is null.
        var items = ToOligoMotifItems(report, zscore, includeVariance: zscore);
        return new OligoAnalysisResultDto(items, report.OligoLength, StrandName(report.Strands), report.CountOverlapping,
            report.TotalOccurrences, report.TestedPatterns, FiniteOrNull(report.PossibleOligos))
        {
            SequenceCount = report.SequenceCount,
            Degenerate = "none",
            SequenceType = type,
            AlphabetSize = report.AlphabetSize,
        };
    }
    private static global::Seqeron.Genomics.OverrepresentedWords.OligoBackgroundModel ToOligoBackground(
        string background, int markovOrder, double[]? residueFrequencies,
        Dictionary<string, double>? oligoFrequencies, double pseudoFrequency, bool strandInsensitive)
    {
        switch ((background ?? string.Empty).Trim().ToLowerInvariant())
        {
            case "input":
                return global::Seqeron.Genomics.OverrepresentedWords.OligoBackgroundModel.BernoulliFromInput;
            case "equiprobable":
                return global::Seqeron.Genomics.OverrepresentedWords.OligoBackgroundModel.Equiprobable;
            case "bernoulli":
                if (residueFrequencies is null || residueFrequencies.Length != 4)
                    throw new ArgumentException("Background 'bernoulli' requires 4 residueFrequencies (A,C,G,T).", nameof(residueFrequencies));
                return global::Seqeron.Genomics.OverrepresentedWords.OligoBackgroundModel.Bernoulli(residueFrequencies);
            case "markov":
                if (markovOrder < 0)
                    throw new ArgumentOutOfRangeException(nameof(markovOrder), "markovOrder must be >= 0.");
                return global::Seqeron.Genomics.OverrepresentedWords.OligoBackgroundModel.MarkovFromInput(markovOrder);
            case "lexicon":
                return global::Seqeron.Genomics.OverrepresentedWords.OligoBackgroundModel.Lexicon;
            case "markov_table":
                if (oligoFrequencies is null || oligoFrequencies.Count == 0)
                    throw new ArgumentException("Background 'markov_table' requires oligoFrequencies.", nameof(oligoFrequencies));
                return global::Seqeron.Genomics.OverrepresentedWords.OligoBackgroundModel.MarkovFromOligoFrequencies(
                    oligoFrequencies, pseudoFrequency, strandInsensitive);
            default:
                throw new ArgumentException(
                    "Background must be 'input', 'equiprobable', 'bernoulli', 'markov', 'markov_table' or 'lexicon'.", nameof(background));
        }
    }
    private static global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy ToOligoDegeneracy(string degenerate)
        => (degenerate ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "none" => global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy.None,
            "onen" => global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy.OneN,
            "onedeg" => global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy.OneDegenerate,
            _ => throw new ArgumentException("degenerate must be 'none', 'oneN' or 'onedeg'.", nameof(degenerate)),
        };
    private static string? DegeneracyName(global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy degeneracy) => degeneracy switch
    {
        global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy.OneN => "oneN",
        global::Seqeron.Genomics.OverrepresentedWords.OligoDegeneracy.OneDegenerate => "onedeg",
        _ => null,
    };
    // AnalyzeOligos / AnalyzeOligoStrings rows → DTO; rows without a valid expected frequency (RSAT "NA", not tested) are omitted.
    private static OligoMotifItem[] ToOligoMotifItems(
        global::Seqeron.Genomics.OverrepresentedWords.OligoAnalysisReport report, bool zscore, bool includeVariance)
        => report.Patterns
            .Where(p => p.FittedDistribution != global::Seqeron.Genomics.OverrepresentedWords.OligoFittedDistribution.None)
            .Select(p => new OligoMotifItem(
                p.Pattern, p.ReverseComplement, p.Occurrences, p.Positions.Select(o => o.Position).ToArray(),
                p.ExpectedFrequency, p.ExpectedOccurrences, FiniteOrNull(p.Ratio),
                p.OccurrenceProbability, p.OccurrenceEValue, p.OccurrenceSignificance)
            {
                SequenceIndices = p.Positions.Select(o => o.SequenceIndex).ToArray(),
                Overlaps = p.Overlaps,
                ObservedFrequency = p.ObservedFrequency,
                ExpectedVariance = includeVariance ? FiniteOrNull(p.ExpectedVariance) : null,
                OverlapCoefficient = zscore ? FiniteOrNull(p.OverlapCoefficient) : null,
                ZScore = zscore ? FiniteOrNull(p.ZScore) : null,
                FittedDistribution = p.FittedDistribution.ToString(),
                LexiconSegmentation = p.LexiconSegmentation is { } seg ? $"{seg.Prefix}|{seg.Suffix}" : null,
            })
            .ToArray();
    private static global::Seqeron.Genomics.OverrepresentedWords.OligoStrandMode ToOligoStrandMode(string strands)
        => (strands ?? string.Empty).Trim().ToLowerInvariant() switch
        {
            "single" => global::Seqeron.Genomics.OverrepresentedWords.OligoStrandMode.Single,
            "both" => global::Seqeron.Genomics.OverrepresentedWords.OligoStrandMode.Both,
            _ => throw new ArgumentException("Strands must be 'single' or 'both'.", nameof(strands)),
        };
    private static string StrandName(global::Seqeron.Genomics.OverrepresentedWords.OligoStrandMode mode)
        => mode == global::Seqeron.Genomics.OverrepresentedWords.OligoStrandMode.Both ? "both" : "single";

}
