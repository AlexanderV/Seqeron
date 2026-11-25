using System.Text.Json.Serialization;

namespace Seqeron.Mcp.Analysis.Tools;

// B09 A3-8 / A3-12 / A3-13: results of the EMBOSS program ports (etandem / equicktandem, fuzzpro / fuzztran, seqmatchall, palindrome / einverted).

/// <summary>A match of a <c>nucmer_clusters</c> cluster: 1-based starts (the second one in the reverse-complemented query on the reverse strand), the length left after the overlap cut and the overlap.</summary>
public record NucmerClusterMatchItem(long Start1, long Start2, long Length, long Overlap);

/// <summary>A cluster of <c>nucmer_clusters</c>: its strand and matches.</summary>
public record NucmerClusterItem(string Strand, NucmerClusterMatchItem[] Matches);

/// <summary>Result of <c>nucmer_clusters</c>: the clusters of the forward strand, then of the reverse-complement strand.</summary>
public record NucmerClustersResult(NucmerClusterItem[] Clusters);

/// <summary>The text a Mash command writes (<c>mash_bounds</c>, <c>mash_dist</c>, <c>mash_triangle</c>), byte for byte as the program prints it.</summary>
public record MashTextResult(string Text);

/// <summary>A Mash sketch file (<c>mash_sketch</c>, <c>mash_paste</c>): <c>msh:</c> + the base64 of the <c>.msh</c> bytes (accepted as an input by the other Mash tools), their number and the number of sketches in it.</summary>
public record MashSketchFileResult(string Msh, int Bytes, int Sketches);
