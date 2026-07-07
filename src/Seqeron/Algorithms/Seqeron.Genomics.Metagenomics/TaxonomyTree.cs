namespace Seqeron.Genomics.Metagenomics;

/// <summary>
/// A single taxon node in a <see cref="TaxonomyTree"/>.
/// </summary>
/// <param name="Id">Stable integer taxon id (NCBI-style). The root is conventionally id 1.</param>
/// <param name="Name">Human-readable taxon name (e.g. <c>"Escherichia coli"</c>).</param>
/// <param name="Rank">Taxonomic rank label (e.g. <c>"species"</c>, <c>"genus"</c>, <c>"root"</c>).</param>
/// <param name="ParentId">Id of the parent taxon. The root is its own parent.</param>
public readonly record struct TaxonNode(int Id, string Name, string Rank, int ParentId);

/// <summary>
/// An NCBI-style taxonomy tree: a set of <see cref="TaxonNode"/> records connected by
/// parent links, supporting parent-chain walks and lowest-common-ancestor (LCA) queries.
/// </summary>
/// <remarks>
/// <para>
/// This is the data model Kraken builds its database and classification on (Wood &amp; Salzberg,
/// 2014, <em>Genome Biology</em> 15:R46): "a database that contains records consisting of a k-mer
/// and the LCA of all organisms whose genomes contain that k-mer". The same LCA operation is used
/// both when collapsing a k-mer's owning taxa at database-build time and when breaking ties between
/// equally-scoring root-to-leaf paths at classification time.
/// </para>
/// <para>
/// The conventional unclassified / root taxon is id <see cref="RootId"/> (Kraken uses 0 for
/// "unclassified" and 1 for the taxonomy root; here a single root represents both — a read with no
/// hits is reported with <see cref="RootId"/>).
/// </para>
/// </remarks>
public sealed class TaxonomyTree
{
    /// <summary>Conventional id of the taxonomy root (its own parent).</summary>
    public const int RootId = 1;

    private readonly Dictionary<int, TaxonNode> _nodes;

    /// <summary>
    /// Builds a taxonomy tree from a set of nodes. Exactly one node must be the root
    /// (a node that is its own parent); every other node's parent must be present.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="nodes"/> is null.</exception>
    /// <exception cref="ArgumentException">
    /// A duplicate id is supplied, there is not exactly one self-parented root, or a non-root
    /// node references a parent that is not in the set.
    /// </exception>
    public TaxonomyTree(IEnumerable<TaxonNode> nodes)
    {
        ArgumentNullException.ThrowIfNull(nodes);

        _nodes = new Dictionary<int, TaxonNode>();
        foreach (var node in nodes)
        {
            if (!_nodes.TryAdd(node.Id, node))
                throw new ArgumentException($"Duplicate taxon id {node.Id}.", nameof(nodes));
        }

        if (_nodes.Count == 0)
            throw new ArgumentException("A taxonomy tree must contain at least the root node.", nameof(nodes));

        int rootCount = 0;
        int rootId = 0;
        foreach (var node in _nodes.Values)
        {
            if (node.ParentId == node.Id)
            {
                rootCount++;
                rootId = node.Id;
            }
            else if (!_nodes.ContainsKey(node.ParentId))
            {
                throw new ArgumentException(
                    $"Taxon {node.Id} references missing parent {node.ParentId}.", nameof(nodes));
            }
        }

        if (rootCount != 1)
            throw new ArgumentException(
                $"A taxonomy tree must have exactly one self-parented root; found {rootCount}.", nameof(nodes));

        Root = rootId;
    }

    /// <summary>The id of this tree's root taxon.</summary>
    public int Root { get; }

    /// <summary>Number of nodes (taxa) in the tree.</summary>
    public int Count => _nodes.Count;

    /// <summary>Whether a taxon id exists in this tree.</summary>
    public bool Contains(int taxonId) => _nodes.ContainsKey(taxonId);

    /// <summary>Gets the node for a taxon id.</summary>
    /// <exception cref="KeyNotFoundException">The id is not in the tree.</exception>
    public TaxonNode GetNode(int taxonId)
    {
        if (_nodes.TryGetValue(taxonId, out var node))
            return node;
        throw new KeyNotFoundException($"Taxon id {taxonId} is not in the tree.");
    }

    /// <summary>The parent id of a taxon (the root's parent is itself).</summary>
    /// <exception cref="KeyNotFoundException">The id is not in the tree.</exception>
    public int GetParent(int taxonId) => GetNode(taxonId).ParentId;

    /// <summary>
    /// The path from a taxon up to and including the root, ordered taxon-first
    /// (index 0 is <paramref name="taxonId"/>, the last element is <see cref="Root"/>).
    /// </summary>
    /// <exception cref="KeyNotFoundException">The id is not in the tree.</exception>
    public IReadOnlyList<int> GetPathToRoot(int taxonId)
    {
        if (!_nodes.ContainsKey(taxonId))
            throw new KeyNotFoundException($"Taxon id {taxonId} is not in the tree.");

        var path = new List<int>();
        int current = taxonId;
        while (true)
        {
            path.Add(current);
            int parent = _nodes[current].ParentId;
            if (parent == current)
                break; // reached the root
            current = parent;
        }
        return path;
    }

    /// <summary>Depth of a taxon (root has depth 0).</summary>
    public int GetDepth(int taxonId) => GetPathToRoot(taxonId).Count - 1;

    /// <summary>
    /// Whether <paramref name="ancestorId"/> is on the path from <paramref name="taxonId"/> to the
    /// root (a taxon is its own ancestor).
    /// </summary>
    public bool IsAncestorOf(int ancestorId, int taxonId)
    {
        if (!_nodes.ContainsKey(ancestorId) || !_nodes.ContainsKey(taxonId))
            return false;

        int current = taxonId;
        while (true)
        {
            if (current == ancestorId)
                return true;
            int parent = _nodes[current].ParentId;
            if (parent == current)
                return false;
            current = parent;
        }
    }

    /// <summary>
    /// Lowest common ancestor of two taxa: the deepest node that is an ancestor of both.
    /// </summary>
    /// <remarks>
    /// Siblings → their shared parent; an ancestor/descendant pair → the ancestor; a node with
    /// itself → itself; taxa whose only shared ancestor is the root → the root.
    /// </remarks>
    /// <exception cref="KeyNotFoundException">Either id is not in the tree.</exception>
    public int Lca(int a, int b)
    {
        if (!_nodes.ContainsKey(a))
            throw new KeyNotFoundException($"Taxon id {a} is not in the tree.");
        if (!_nodes.ContainsKey(b))
            throw new KeyNotFoundException($"Taxon id {b} is not in the tree.");

        // Collect a's ancestors (including a) into a set, then walk b upward to the first hit.
        var ancestorsOfA = new HashSet<int>();
        int current = a;
        while (true)
        {
            ancestorsOfA.Add(current);
            int parent = _nodes[current].ParentId;
            if (parent == current)
                break;
            current = parent;
        }

        current = b;
        while (true)
        {
            if (ancestorsOfA.Contains(current))
                return current;
            int parent = _nodes[current].ParentId;
            if (parent == current)
                return current; // root: the universal common ancestor
            current = parent;
        }
    }

    /// <summary>
    /// Lowest common ancestor of an arbitrary non-empty set of taxa, computed by folding
    /// <see cref="Lca(int,int)"/> across them. The LCA of a single taxon is that taxon.
    /// </summary>
    /// <exception cref="ArgumentNullException"><paramref name="taxa"/> is null.</exception>
    /// <exception cref="ArgumentException"><paramref name="taxa"/> is empty.</exception>
    /// <exception cref="KeyNotFoundException">A supplied id is not in the tree.</exception>
    public int Lca(IEnumerable<int> taxa)
    {
        ArgumentNullException.ThrowIfNull(taxa);

        int? acc = null;
        foreach (int t in taxa)
            acc = acc is null ? t : Lca(acc.Value, t);

        if (acc is null)
            throw new ArgumentException("LCA requires at least one taxon.", nameof(taxa));
        return acc.Value;
    }
}
