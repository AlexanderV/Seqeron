using System.ComponentModel;
using System.Globalization;
using ModelContextProtocol.Server;

namespace Seqeron.Mcp.Alignment.Tools;

// B09 A2-28c / A2-29b: EMBOSS supermatcher and needleall, and the ksw2 extension functions; delegate to EmbossPairwisePrograms.Supermatcher / NeedleAll and Ksw2Aligner.
public partial class AlignmentTools
{

    [McpServerTool(Name = "ksw2_extend", Title = "Pairwise — ksw2 Extension Alignment (z-drop, end bonus, band)", ReadOnly = true)]
    [Description("The ksw2 / minimap2 extension-alignment functions: global or extension alignment of a query against a target with affine (q + k*e) or dual-affine (min(q + k*e, q2 + k*e2)) gap cost, a diagonal band, the Z-drop rule, the end bonus, the ambiguous-base (N) penalty and CIGAR output. variant: 'extz' (ksw_extz, scalar affine), 'extz2' (ksw_extz2_sse, minimap2 affine), 'extd' (ksw_extd, scalar dual affine), 'extd2' (ksw_extd2_sse, minimap2 -O/-E dual affine; the default), 'exts2' (ksw_exts2_sse, the spliced extension of minimap2 -x splice: introns of at least (q2 - q) / e - 1 bases are N operations of the CIGAR, charged q2 plus a donor / acceptor score: spliceStrand forward = GT-AG, reverse = CT-AC, nonCanonical = minimap2 -C, junctionBonus with a junctions string (digits: bit 0 forward donor, 1 forward acceptor, 2 reverse acceptor, 3 reverse donor; the program's --junc-bed), spliceModel 'miniprot' (-J 1) or 'old' (-J 0, with spliceFlank), reverseCigar), 'll' (ksw_ll_i16, the 16-bit striped local alignment minimap2 uses to test a Z-drop for an inversion: returns the best local score in max and its end coordinates in maxQ / maxT, gap cost gapOpen + k*gapExtend). Bases are encoded A,C,G,T = 0..3 and anything else = N (4); match / mismatch / ambiguous give the matrix as minimap2 builds it (-A, -B as a penalty, -N penalty of the wildcard). Penalties are POSITIVE (minimap2 -O q,q2 -E e,e2; a gap of length k costs q + k*e; the second extension penalty e2 may be 0). w = band half-width (-1: none; the SIMD variants need w >= 1), zdrop = off-diagonal drop-off to stop an extension (-1: none; minimap2 -z), endBonus (the extension end bonus; SIMD variants, extensionOnly). extensionOnly = the path ends at the best cell (or at the query end when mqe + endBonus > max); otherwise it is global (corner to corner). Returned: max, zdropped, maxQ/maxT, mqe (best score reaching the query end) with mqeT, mte with mteQ, score (corner score; -1073741824 = not reached), reachEnd and the CIGAR (M, I = query only, D = target only, N = intron). VERIFICATION: the SIMD variants extz2, extd2, exts2 and ll equal the C functions of minimap2 v2.26 compiled with SSE4.1 on 500 random inputs each (every field of ksw_extz_t and the CIGAR: band, Z-drop, end bonus, extension only, left / right gaps, reversed CIGAR, approximate maximum / drop, score only); extz / extd (scalar) equal an independent dynamic program on 330 cases per variant (max, mqe, mte, score; every CIGAR re-scored along its path).")]
    public static Ksw2ExtendResult Ksw2Extend(
        [Description("Query sequence (ACGT, other letters = N).")] string query,
        [Description("Target sequence.")] string target,
        [Description("'extd2' (default), 'extz2', 'extd', 'extz', 'exts2' (spliced extension) or 'll' (striped local alignment).")] string variant = "extd2",
        [Description("Score of a match (minimap2 -A).")] int match = 2,
        [Description("Penalty of a mismatch, positive (minimap2 -B).")] int mismatch = 4,
        [Description("Penalty of the ambiguous base N, positive (minimap2 -N / sc_ambi); 0 gives -e.")] int ambiguous = 1,
        [Description("Gap open penalty q (minimap2 -O first value), >= 0.")] int gapOpen = 4,
        [Description("Gap extension penalty e (minimap2 -E first value), >= 1.")] int gapExtend = 2,
        [Description("Second gap open penalty q2 (dual variants; minimap2 -O second value).")] int gapOpen2 = 24,
        [Description("Second gap extension penalty e2 (dual variants; minimap2 -E second value).")] int gapExtend2 = 1,
        [Description("Band half-width (-1 = none).")] int band = -1,
        [Description("Z-drop (-1 = none; minimap2 -z).")] int zdrop = -1,
        [Description("End bonus (SIMD variants, extensionOnly).")] int endBonus = 0,
        [Description("Extension alignment: the path ends at the best cell (or the query end with the end bonus).")] bool extensionOnly = false,
        [Description("Right-align gaps (KSW_EZ_RIGHT).")] bool rightAlignGaps = false,
        [Description("Score only, no CIGAR.")] bool scoreOnly = false,
        [Description("exts2: 'none' (no splice signal), 'forward' (GT-AG, minimap2 -u f) or 'reverse' (CT-AC, the reverse complement).")] string spliceStrand = "none",
        [Description("exts2: the penalty of a non-canonical splice site (minimap2 -C, noncan).")] int nonCanonical = 0,
        [Description("exts2: the bonus of an annotated junction (minimap2 --junc-bonus).")] int junctionBonus = 0,
        [Description("exts2: 'miniprot' (minimap2 -J 1, the default) or 'old' (-J 0).")] string spliceModel = "miniprot",
        [Description("exts2, old model: the flanking-sequence bonus (--splice-flank).")] bool spliceFlank = false,
        [Description("exts2: the CIGAR is reported last operation first (KSW_EZ_REV_CIGAR; the left extension of minimap2).")] bool reverseCigar = false,
        [Description("exts2: junction flags of the target positions, one digit 0-9 per target base (bit 0 forward donor, bit 1 forward acceptor, bit 2 reverse acceptor, bit 3 reverse donor), or null.")] string? junctions = null,
        CancellationToken cancellationToken = default)
    {
        using var cancelScope = global::Seqeron.Genomics.Alignment.SequenceAligner.BeginCancellationScope(cancellationToken);
        RequireSequence(query, nameof(query));
        RequireSequence(target, nameof(target));
        if (gapOpen < 0 || gapExtend < 1 || gapOpen > 127 || gapExtend > 127)
            throw new ArgumentOutOfRangeException(nameof(gapOpen), "Gap open must be 0..127 and gap extension 1..127 (they are bytes in ksw2).");
        var flags = global::Seqeron.Genomics.Alignment.Ksw2Option.None;
        if (extensionOnly)
            flags |= global::Seqeron.Genomics.Alignment.Ksw2Option.ExtensionOnly;
        if (rightAlignGaps)
            flags |= global::Seqeron.Genomics.Alignment.Ksw2Option.Right;
        if (scoreOnly)
            flags |= global::Seqeron.Genomics.Alignment.Ksw2Option.ScoreOnly;
        byte[] q = global::Seqeron.Genomics.Alignment.Ksw2Aligner.EncodeNucleotides(query);
        byte[] t = global::Seqeron.Genomics.Alignment.Ksw2Aligner.EncodeNucleotides(target);
        var mat = global::Seqeron.Genomics.Alignment.Ksw2Aligner.SimpleMatrix(5, match, -mismatch, ambiguous);
        string v = (variant ?? "extd2").Trim().ToLowerInvariant();
        if (v == "ll")
        {
            var ll = global::Seqeron.Genomics.Alignment.Ksw2Aligner.LocalScore(q, t, 5, mat, gapOpen, gapExtend);
            return new Ksw2ExtendResult(ll.Score, false, ll.QueryEnd, ll.TargetEnd, global::Seqeron.Genomics.Alignment.Ksw2Aligner.NegInf, -1, global::Seqeron.Genomics.Alignment.Ksw2Aligner.NegInf, -1, global::Seqeron.Genomics.Alignment.Ksw2Aligner.NegInf, false, "*");
        }

        if (v == "exts2")
        {
            var sflags = flags;
            string strand = (spliceStrand ?? "none").Trim().ToLowerInvariant();
            if (strand == "forward")
                sflags |= global::Seqeron.Genomics.Alignment.Ksw2Option.SpliceForward;
            else if (strand == "reverse")
                sflags |= global::Seqeron.Genomics.Alignment.Ksw2Option.SpliceReverse;
            else if (strand != "none")
                throw new ArgumentException($"spliceStrand must be none, forward or reverse, got '{spliceStrand}'.", nameof(spliceStrand));
            string model = (spliceModel ?? "miniprot").Trim().ToLowerInvariant();
            if (model == "miniprot")
                sflags |= global::Seqeron.Genomics.Alignment.Ksw2Option.SpliceComplex;
            else if (model != "old")
                throw new ArgumentException($"spliceModel must be miniprot or old, got '{spliceModel}'.", nameof(spliceModel));
            if (spliceFlank)
                sflags |= global::Seqeron.Genomics.Alignment.Ksw2Option.SpliceFlank;
            if (reverseCigar)
                sflags |= global::Seqeron.Genomics.Alignment.Ksw2Option.ReverseCigar;
            byte[]? junc = null;
            if (!string.IsNullOrEmpty(junctions))
            {
                if (junctions.Length != target.Length || junctions.Any(c => c < '0' || c > '9'))
                    throw new ArgumentException("junctions must have one digit 0-9 per target base.", nameof(junctions));
                junc = [.. junctions.Select(c => (byte)(c - '0'))];
            }

            var sr = global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendS2(q, t, 5, mat, gapOpen, gapExtend, gapOpen2, nonCanonical, zdrop, junctionBonus, sflags, junc);
            return new Ksw2ExtendResult(sr.Max, sr.ZDropped, sr.MaxQ, sr.MaxT, sr.Mqe, sr.MqeT, sr.Mte, sr.MteQ, sr.Score, sr.ReachEnd, sr.CigarString);
        }

        var r = v switch
        {
            "extz" => global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendZ(q, t, 5, mat, gapOpen, gapExtend, band, zdrop, flags),
            "extz2" => global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendZ2(q, t, 5, mat, gapOpen, gapExtend, band, zdrop, endBonus, flags),
            "extd" => global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendD(q, t, 5, mat, gapOpen, gapExtend, gapOpen2, gapExtend2, band, zdrop, flags),
            "extd2" => global::Seqeron.Genomics.Alignment.Ksw2Aligner.ExtendD2(q, t, 5, mat, gapOpen, gapExtend, gapOpen2, gapExtend2, band, zdrop, endBonus, flags),
            _ => throw new ArgumentException($"variant must be extz, extz2, extd, extd2, exts2 or ll, got '{variant}'.", nameof(variant)),
        };
        return new Ksw2ExtendResult(r.Max, r.ZDropped, r.MaxQ, r.MaxT, r.Mqe, r.MqeT, r.Mte, r.MteQ, r.Score, r.ReachEnd, r.CigarString);
    }
}
