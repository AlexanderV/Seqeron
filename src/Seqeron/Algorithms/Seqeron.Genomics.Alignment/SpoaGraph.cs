using System.Text;

namespace Seqeron.Genomics.Alignment;

/// <summary>
/// The partial order graph of SPOA 4.0.8 (<c>spoa::Graph</c>, <c>src/graph.cpp</c>): nodes with codes, weighted edges labelled by the sequences that use them, aligned-node rings, the topological order,
/// the heaviest-bundle consensus, the multiple sequence alignment, the DOT and the GFA text of the program. Alignments are lists of (node id or -1, sequence position or -1) pairs.
/// </summary>
public sealed class SpoaGraph
{
    /// <summary>A node of the graph.</summary>
    public sealed class Node
    {
        internal Node(uint id, uint code)
        {
            Id = id;
            Code = code;
        }

        /// <summary>The index of the node in the order of creation.</summary>
        public uint Id { get; }

        /// <summary>The code of the node's character.</summary>
        public uint Code { get; }

        /// <summary>The edges that end in this node.</summary>
        public List<Edge> InEdges { get; } = [];

        /// <summary>The edges that start in this node.</summary>
        public List<Edge> OutEdges { get; } = [];

        /// <summary>The nodes aligned to this one (other characters in the same column).</summary>
        public List<Node> AlignedNodes { get; } = [];

        /// <summary>The next node of sequence <paramref name="label"/>, or null.</summary>
        public Node? Successor(uint label)
        {
            foreach (var edge in OutEdges)
            {
                if (edge.Labels.Contains(label))
                    return edge.Head;
            }

            return null;
        }
    }

    /// <summary>An edge of the graph.</summary>
    public sealed class Edge
    {
        internal Edge(Node tail, Node head, uint label, uint weight)
        {
            Tail = tail;
            Head = head;
            Labels = [label];
            Weight = weight;
        }

        /// <summary>The start node.</summary>
        public Node Tail { get; }

        /// <summary>The end node.</summary>
        public Node Head { get; }

        /// <summary>The sequences that use the edge.</summary>
        public List<uint> Labels { get; }

        /// <summary>The total weight (<c>std::int64_t</c>).</summary>
        public long Weight { get; private set; }

        internal void AddSequence(uint label, uint weight)
        {
            Labels.Add(label);
            Weight += weight;
        }
    }

    private uint _numCodes;
    private readonly int[] _coder = Enumerable.Repeat(-1, 256).ToArray();
    private readonly int[] _decoder = Enumerable.Repeat(-1, 256).ToArray();
    private readonly List<Node?> _sequences = [];
    private readonly List<Node> _nodes = [];
    private readonly List<Node> _rankToNode = [];
    private readonly List<Node> _consensus = [];

    /// <summary>The number of different characters seen.</summary>
    public uint NumCodes => _numCodes;

    /// <summary>The nodes in the order of creation.</summary>
    public IReadOnlyList<Node> Nodes => _nodes;

    /// <summary>The nodes in topological order.</summary>
    public IReadOnlyList<Node> RankToNode => _rankToNode;

    /// <summary>The first node of every added sequence (null for an empty alignment of an empty sequence).</summary>
    public IReadOnlyList<Node?> Sequences => _sequences;

    /// <summary>The consensus path of the last <see cref="GenerateConsensus"/> / <see cref="GenerateMultipleSequenceAlignment"/> call.</summary>
    public IReadOnlyList<Node> Consensus => _consensus;

    /// <summary>The character of a code.</summary>
    public char Decoder(uint code) => (char)(byte)_decoder[code];

    private Node AddNode(uint code)
    {
        var node = new Node((uint)_nodes.Count, code);
        _nodes.Add(node);
        return node;
    }

    private void AddEdge(Node tail, Node head, uint weight)
    {
        foreach (var edge in tail.OutEdges)
        {
            if (edge.Head == head)
            {
                edge.AddSequence((uint)_sequences.Count, weight);
                return;
            }
        }

        var created = new Edge(tail, head, (uint)_sequences.Count, weight);
        tail.OutEdges.Add(created);
        head.InEdges.Add(created);
    }

    private Node? AddSequence(string sequence, uint[] weights, uint begin, uint end)
    {
        if (begin == end)
            return null;
        Node? prev = null;
        for (uint i = begin; i < end; i++)
        {
            var curr = AddNode((uint)_coder[sequence[(int)i]]);
            if (prev is not null)
                AddEdge(prev, curr, unchecked(weights[i - 1] + weights[i]));
            prev = curr;
        }

        return _nodes[_nodes.Count - (int)(end - begin)];
    }

    /// <summary>Adds an aligned sequence with the weight 1 per base (<c>Graph::AddAlignment(alignment, sequence)</c>).</summary>
    public void AddAlignment(IReadOnlyList<(int Node, int Position)> alignment, string sequence, uint weight = 1)
    {
        var weights = new uint[sequence.Length];
        Array.Fill(weights, weight);
        AddAlignment(alignment, sequence, weights);
    }

    /// <summary>Adds an aligned sequence with Phred weights (<c>quality[i] - 33</c>, as unsigned 32-bit values).</summary>
    public void AddAlignment(IReadOnlyList<(int Node, int Position)> alignment, string sequence, string quality)
    {
        var weights = new uint[quality.Length];
        for (int i = 0; i < quality.Length; i++)
            weights[i] = unchecked((uint)(quality[i] - 33));
        AddAlignment(alignment, sequence, weights);
    }

    /// <summary>Adds an aligned sequence with explicit weights.</summary>
    /// <exception cref="ArgumentException">The alignment is invalid (<c>std::invalid_argument</c> of SPOA).</exception>
    public void AddAlignment(IReadOnlyList<(int Node, int Position)> alignment, string sequence, uint[] weights)
    {
        uint sequenceLen = (uint)sequence.Length;
        if (sequenceLen == 0)
            return;
        if (sequenceLen != weights.Length)
            throw new ArgumentException("[spoa::Graph::AddAlignment] error: sequence and weights are of unequal size!");
        foreach (char c in sequence)
        {
            if (_coder[c] == -1)
            {
                _coder[c] = (int)_numCodes;
                _decoder[_numCodes++] = c;
            }
        }

        if (alignment.Count == 0)
        {
            _sequences.Add(AddSequence(sequence, weights, 0, sequenceLen));
            TopologicalSort();
            return;
        }

        var valid = new List<uint>();
        foreach (var it in alignment)
        {
            if (it.Position == -1)
                continue;
            if (it.Position < 0 || it.Position >= (int)sequenceLen)
                throw new ArgumentException("[spoa::Graph::AddAlignment] error: invalid alignment");
            valid.Add((uint)it.Position);
        }

        if (valid.Count == 0)
            throw new ArgumentException("[spoa::Graph::AddAlignment] error: missing sequence in alignment");
        var begin = AddSequence(sequence, weights, 0, valid[0]);
        var prev = begin is not null ? _nodes[^1] : null;
        var last = AddSequence(sequence, weights, valid[^1] + 1, sequenceLen);
        foreach (var it in alignment)
        {
            if (it.Position == -1)
                continue;
            uint code = (uint)_coder[sequence[it.Position]];
            Node? curr = null;
            if (it.Node == -1)
            {
                curr = AddNode(code);
            }
            else
            {
                var jt = _nodes[it.Node];
                if (jt.Code == code)
                {
                    curr = jt;
                }
                else
                {
                    foreach (var kt in jt.AlignedNodes)
                    {
                        if (kt.Code == code)
                        {
                            curr = kt;
                            break;
                        }
                    }

                    if (curr is null)
                    {
                        curr = AddNode(code);
                        foreach (var kt in jt.AlignedNodes)
                        {
                            kt.AlignedNodes.Add(curr);
                            curr.AlignedNodes.Add(kt);
                        }

                        jt.AlignedNodes.Add(curr);
                        curr.AlignedNodes.Add(jt);
                    }
                }
            }

            begin ??= curr;
            if (prev is not null)
                AddEdge(prev, curr, unchecked(weights[it.Position - 1] + weights[it.Position]));
            prev = curr;
        }

        if (last is not null)
            AddEdge(prev!, last, unchecked(weights[valid[^1]] + weights[valid[^1] + 1]));
        _sequences.Add(begin);
        TopologicalSort();
    }

    private void TopologicalSort()
    {
        _rankToNode.Clear();
        var marks = new byte[_nodes.Count];
        var ignored = new bool[_nodes.Count];
        var stack = new Stack<Node>();
        foreach (var it in _nodes)
        {
            if (marks[it.Id] != 0)
                continue;
            stack.Push(it);
            while (stack.Count > 0)
            {
                var curr = stack.Peek();
                bool isValid = true;
                if (marks[curr.Id] != 2)
                {
                    foreach (var jt in curr.InEdges)
                    {
                        if (marks[jt.Tail.Id] != 2)
                        {
                            stack.Push(jt.Tail);
                            isValid = false;
                        }
                    }

                    if (!ignored[curr.Id])
                    {
                        foreach (var jt in curr.AlignedNodes)
                        {
                            if (marks[jt.Id] != 2)
                            {
                                stack.Push(jt);
                                ignored[jt.Id] = true;
                                isValid = false;
                            }
                        }
                    }

                    if (isValid)
                    {
                        marks[curr.Id] = 2;
                        if (!ignored[curr.Id])
                        {
                            _rankToNode.Add(curr);
                            foreach (var jt in curr.AlignedNodes)
                                _rankToNode.Add(jt);
                        }
                    }
                    else
                    {
                        marks[curr.Id] = 1;
                    }
                }

                if (isValid)
                    stack.Pop();
            }
        }
    }

    private uint[] InitializeMultipleSequenceAlignment(out uint rowSize)
    {
        var dst = new uint[_nodes.Count];
        uint j = 0;
        int index = 0;
        while (index < _rankToNode.Count)
        {
            var it = _rankToNode[index];
            dst[it.Id] = j;
            index++;
            foreach (var jt in it.AlignedNodes)
            {
                dst[jt.Id] = j;
                index++;
            }

            j++;
        }

        rowSize = j;
        return dst;
    }

    /// <summary>The rows of the multiple sequence alignment, optionally followed by the consensus row.</summary>
    public List<string> GenerateMultipleSequenceAlignment(bool includeConsensus = false)
    {
        var nodeIdToColumn = InitializeMultipleSequenceAlignment(out uint rowSize);
        var dst = new List<string>();
        for (int i = 0; i < _sequences.Count; i++)
        {
            var row = new char[rowSize];
            Array.Fill(row, '-');
            var it = _sequences[i];
            while (it is not null)
            {
                row[nodeIdToColumn[it.Id]] = Decoder(it.Code);
                it = it.Successor((uint)i);
            }

            dst.Add(new string(row));
        }

        if (includeConsensus)
        {
            TraverseHeaviestBundle();
            var row = new char[rowSize];
            Array.Fill(row, '-');
            foreach (var it in _consensus)
                row[nodeIdToColumn[it.Id]] = Decoder(it.Code);
            dst.Add(new string(row));
        }

        return dst;
    }

    /// <summary>The consensus sequence (the heaviest bundle).</summary>
    public string GenerateConsensus()
    {
        TraverseHeaviestBundle();
        var sb = new StringBuilder();
        foreach (var it in _consensus)
            sb.Append(Decoder(it.Code));
        return sb.ToString();
    }

    private void TraverseHeaviestBundle()
    {
        if (_rankToNode.Count == 0)
            return;
        var predecessors = new Node?[_nodes.Count];
        var scores = new long[_nodes.Count];
        Array.Fill(scores, -1L);
        Node? max = null;
        foreach (var it in _rankToNode)
        {
            foreach (var jt in it.InEdges)
            {
                if (scores[it.Id] < jt.Weight || (scores[it.Id] == jt.Weight && scores[predecessors[it.Id]!.Id] <= scores[jt.Tail.Id]))
                {
                    scores[it.Id] = jt.Weight;
                    predecessors[it.Id] = jt.Tail;
                }
            }

            if (predecessors[it.Id] is not null)
                scores[it.Id] += scores[predecessors[it.Id]!.Id];
            if (max is null || scores[max.Id] < scores[it.Id])
                max = it;
        }

        if (max!.OutEdges.Count > 0)
        {
            var nodeIdToRank = new uint[_nodes.Count];
            for (int i = 0; i < _rankToNode.Count; i++)
                nodeIdToRank[_rankToNode[i].Id] = (uint)i;
            while (max.OutEdges.Count > 0)
                max = BranchCompletion(nodeIdToRank[max.Id], scores, predecessors);
        }

        _consensus.Clear();
        while (predecessors[max.Id] is not null)
        {
            _consensus.Add(max);
            max = predecessors[max.Id]!;
        }

        _consensus.Add(max);
        _consensus.Reverse();
    }

    private Node BranchCompletion(uint rank, long[] scores, Node?[] predecessors)
    {
        var start = _rankToNode[(int)rank];
        foreach (var it in start.OutEdges)
        {
            foreach (var jt in it.Head.InEdges)
            {
                if (jt.Tail != start)
                    scores[jt.Tail.Id] = -1;
            }
        }

        Node? max = null;
        for (int i = (int)rank + 1; i < _rankToNode.Count; i++)
        {
            var it = _rankToNode[i];
            scores[it.Id] = -1;
            predecessors[it.Id] = null;
            foreach (var jt in it.InEdges)
            {
                if (scores[jt.Tail.Id] == -1)
                    continue;
                if (scores[it.Id] < jt.Weight || (scores[it.Id] == jt.Weight && scores[predecessors[it.Id]!.Id] <= scores[jt.Tail.Id]))
                {
                    scores[it.Id] = jt.Weight;
                    predecessors[it.Id] = jt.Tail;
                }
            }

            if (predecessors[it.Id] is not null)
                scores[it.Id] += scores[predecessors[it.Id]!.Id];
            if (max is null || scores[max.Id] < scores[it.Id])
                max = it;
        }

        return max!;
    }

    /// <summary><c>Graph::PrintDot</c>: the DOT text (the program writes it to the file of <c>-d</c>).</summary>
    public string PrintDot()
    {
        var consensusRank = new int[_nodes.Count];
        Array.Fill(consensusRank, -1);
        int rank = 0;
        foreach (var it in _consensus)
            consensusRank[it.Id] = rank++;
        var os = new StringBuilder();
        os.Append("digraph ").Append(_sequences.Count).Append(" {\n  graph [rankdir = LR]\n");
        foreach (var it in _nodes)
        {
            os.Append("  ").Append(it.Id).Append("[label = \"").Append(it.Id).Append(" - ").Append(Decoder(it.Code)).Append('"');
            if (consensusRank[it.Id] != -1)
                os.Append(", style = filled, fillcolor = goldenrod1");
            os.Append("]\n");
            foreach (var jt in it.OutEdges)
            {
                os.Append("  ").Append(it.Id).Append(" -> ").Append(jt.Head.Id).Append(" [label = \"").Append(jt.Weight).Append('"');
                if (consensusRank[it.Id] + 1 == consensusRank[jt.Head.Id])
                    os.Append(", color = goldenrod1");
                os.Append("]\n");
            }

            foreach (var jt in it.AlignedNodes)
            {
                if (jt.Id > it.Id)
                    os.Append("  ").Append(it.Id).Append(" -> ").Append(jt.Id).Append(" [style = dotted, arrowhead = none]\n");
            }
        }

        os.Append("}\n");
        return os.ToString();
    }
}
