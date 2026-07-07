using System.Text.Json.Serialization;
using Seqeron.Mcp.Phylogenetics.Models;

namespace Seqeron.Mcp.Phylogenetics;

[JsonSourceGenerationOptions(WriteIndented = false)]
// Primitives / common containers used as MCP tool inputs.
[JsonSerializable(typeof(string))]
[JsonSerializable(typeof(bool))]
[JsonSerializable(typeof(int))]
[JsonSerializable(typeof(double))]
[JsonSerializable(typeof(string[]))]
[JsonSerializable(typeof(double[]))]
[JsonSerializable(typeof(double[][]))]
[JsonSerializable(typeof(Dictionary<string, string>))]
// Result DTOs (tool outputs).
[JsonSerializable(typeof(PhyloTreeBuildResult))]
[JsonSerializable(typeof(DistanceMatrixResult))]
[JsonSerializable(typeof(PairwiseDistanceResult))]
[JsonSerializable(typeof(ToNewickResult))]
[JsonSerializable(typeof(ParseNewickResult))]
[JsonSerializable(typeof(TreeLeafItem))]
[JsonSerializable(typeof(TreeLeavesResult))]
[JsonSerializable(typeof(TreeLengthResult))]
[JsonSerializable(typeof(TreeDepthResult))]
[JsonSerializable(typeof(RobinsonFouldsDistanceResult))]
[JsonSerializable(typeof(MrcaResult))]
[JsonSerializable(typeof(PatristicDistanceResult))]
[JsonSerializable(typeof(BootstrapSupportItem))]
[JsonSerializable(typeof(BootstrapSupportResult))]
public partial class AppJsonContext : JsonSerializerContext
{
}
