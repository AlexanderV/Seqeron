# Seqeron.Mcp.Phylogenetics

MCP server — **Distance matrices, tree building, Newick I/O, tree statistics and comparison.**

Exposes **13 tools** wrapping the `Seqeron.Genomics` library. Every tool has an
explicit JSON input/output schema, a Schema+Binding test, and per-tool docs under
[`docs/mcp/tools/phylogenetics/`](../../../../docs/mcp/tools/phylogenetics) — see the
campaign ledger [`docs/mcp/MCP_STATUS.md`](../../../../docs/mcp/MCP_STATUS.md).

## Run

```bash
dotnet run --project Seqeron.Mcp.Phylogenetics
```

Register it in any MCP client as a stdio server (`command: dotnet`, `args: ["run","--project","Seqeron.Mcp.Phylogenetics"]`). See [`docs/mcp/README.md`](../../../../docs/mcp/README.md).

## Tools (13)

| Tool | Description |
|------|-------------|
| `bootstrap_support` | Non-parametric bootstrap support for clades of a reference tree. |
| `build_phylogenetic_tree` | Build a phylogenetic tree from a set of named, pre-aligned sequences. |
| `build_tree_from_matrix` | Build a phylogenetic tree directly from a precomputed symmetric distance matrix (UPGMA or NeighborJoining). |
| `distance_matrix` | Compute the symmetric pairwise distance matrix for a list of aligned sequences using the chosen substitution model (PDistance | JukesCant… |
| `mrca` | Most Recent Common Ancestor (MRCA) of two taxa in a rooted tree. |
| `pairwise_distance` | Calculate the evolutionary distance between two aligned sequences under a chosen substitution model: p-distance, Hamming, Jukes-Cantor (J… |
| `parse_newick` | Parse a Newick-format tree string and report a summary: canonical re-serialization (round-trip), taxa list, leaf count, depth, and total… |
| `patristic_distance` | Sum of branch lengths along the unique path between two taxa in a rooted tree (via MRCA). |
| `robinson_foulds_distance` | Robinson-Foulds (symmetric clade-difference) distance between two rooted trees over the same taxon set. |
| `to_newick` | Serialize a phylogenetic tree to Newick format. |
| `tree_depth` | Maximum number of internal-node edges from root to any leaf. |
| `tree_leaves` | Enumerate the leaf (taxon) nodes of a tree, returning each leaf's name and branch length. |
| `tree_length` | Sum of all branch lengths in a tree. |
