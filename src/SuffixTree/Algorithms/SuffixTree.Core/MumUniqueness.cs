namespace SuffixTree;

/// <summary>
/// Uniqueness requirement for maximal unique matches (MUMs), as in MUMmer 3
/// (Kurtz et al. 2004, Genome Biol 5:R12; <c>mummer -mum</c> / <c>-mumreference</c>).
/// </summary>
public enum MumUniqueness
{
    /// <summary>
    /// The matched string occurs exactly once in the reference (tree text) and exactly once in the
    /// query — <c>mummer -mum</c> (MUMmer 3 <c>findmumcandidates</c> + <c>mumuniqueinquery</c>).
    /// </summary>
    Both = 0,

    /// <summary>
    /// The matched string occurs exactly once in the reference; it may repeat in the query —
    /// <c>mummer -mumreference</c> / <c>-mumcand</c> (MUMmer's default; "MUM-candidates",
    /// called maximal almost-unique matches, MAMs, in MUMmer 4 / essaMEM).
    /// </summary>
    Reference = 1,
}
