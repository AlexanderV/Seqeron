using System.Collections.Concurrent;
using System.Globalization;

namespace Seqeron.Genomics.Infrastructure;

/// <summary>
/// Named built-in matrices (audit A2-15): every EMBOSS 6.6.0 residue-pair comparison matrix (<c>EBLOSUM30</c> ... <c>EBLOSUM90</c>, <c>EBLOSUM62-12</c>, <c>EBLOSUMN</c>, <c>EPAM10</c> ... <c>EPAM500</c>,
/// <c>EDNAFULL</c>, <c>EDNAMAT</c>, <c>EDNASIMPLE</c>, <c>ENUC.4.4</c>) read from the EMBOSS data directory by <c>scripts/alignment/generate_emboss_matrices.py</c>, with the NCBI aliases
/// (<c>BLOSUM45</c>, <c>BLOSUM50</c>, <c>BLOSUM62</c>, <c>BLOSUM80</c>, <c>PAM30</c>, <c>PAM70</c>, <c>PAM250</c>, ... = the same cells: 0 differing cells against Biopython 1.88 <c>substitution_matrices</c>).
/// </summary>
public sealed partial class SubstitutionMatrix
{
    private static readonly ConcurrentDictionary<string, SubstitutionMatrix> Registry = new(StringComparer.Ordinal);
    private static readonly Lazy<Dictionary<string, int>> EmbossIndex = new(() =>
    {
        var d = new Dictionary<string, int>(StringComparer.Ordinal);
        for (int i = 0; i < EmbossData.Length; i++)
            d[EmbossData[i].Name] = i;
        return d;
    });

    /// <summary>The EMBOSS file names of the built-in matrices (e.g. <c>EBLOSUM62</c>, <c>EPAM120</c>, <c>EDNAFULL</c>), in the order protein BLOSUM, protein PAM, nucleotide.</summary>
    public static IReadOnlyList<string> BuiltInNames { get; } = EmbossData.Select(d => d.Name).ToArray();

    /// <summary>
    /// Looks up a built-in matrix by name (case-insensitive): an EMBOSS data file name (<c>EBLOSUM62</c>, <c>EPAM250</c>, <c>EDNAFULL</c>, ...) or its NCBI name without the leading <c>E</c>
    /// (<c>BLOSUM62</c>, <c>PAM250</c>, <c>BLOSUM45</c>, <c>PAM30</c>, <c>DNAFULL</c>, <c>NUC.4.4</c> / <c>NUC44</c> = EDNAFULL as before). <c>EBLOSUM62</c> / <c>EPAM250</c> / <c>EDNAFULL</c> return the
    /// instances <see cref="Blosum62"/> / <see cref="Pam250"/> / <see cref="EdnaFull"/>. Returns false for an unknown name.
    /// </summary>
    public static bool TryGet(string? name, out SubstitutionMatrix matrix)
    {
        matrix = null!;
        if (string.IsNullOrWhiteSpace(name))
            return false;
        string key = name.Trim().ToUpperInvariant();
        if (key is "NUC.4.4" or "NUC44" or "DNAFULL")
            key = "EDNAFULL";
        else if (!EmbossIndex.Value.ContainsKey(key) && EmbossIndex.Value.ContainsKey("E" + key))
            key = "E" + key;
        if (!EmbossIndex.Value.TryGetValue(key, out int idx))
            return false;
        matrix = key switch
        {
            "EBLOSUM62" => Blosum62,
            "EPAM250" => Pam250,
            "EDNAFULL" => EdnaFull,
            _ => Registry.GetOrAdd(key, _ => FromEmbossRow(EmbossData[idx])),
        };
        return true;
    }

    /// <summary>Looks up a built-in matrix (<see cref="TryGet"/>).</summary>
    /// <exception cref="ArgumentException">Unknown name; the message lists the names.</exception>
    public static SubstitutionMatrix Get(string name) =>
        TryGet(name, out var m)
            ? m
            : throw new ArgumentException($"Unknown matrix '{name}'. Built-in: {string.Join(", ", BuiltInNames)} (NCBI names without the leading E also accepted; ENUC.4.2 is rectangular and not supported).", nameof(name));

    /// <summary>True for a nucleotide matrix: every symbol is an IUPAC nucleotide code (A C G T U and the ambiguity codes R Y S W K M B D H V N, X).</summary>
    public bool IsNucleotide => Alphabet.All(c => "ACGTURYSWKMBDHVNX".Contains(c));

    private static SubstitutionMatrix FromEmbossRow((string Name, string Description, string Alphabet, string Scores) row)
    {
        int n = row.Alphabet.Length;
        var scores = new int[n, n];
        int k = 0;
        foreach (var token in row.Scores.Split(' '))
        {
            scores[k / n, k % n] = int.Parse(token, CultureInfo.InvariantCulture);
            k++;
        }

        if (k != n * n)
            throw new InvalidOperationException($"Matrix {row.Name}: {k} scores for a {n} x {n} table.");
        return new SubstitutionMatrix(row.Name, row.Alphabet, scores);
    }
}
