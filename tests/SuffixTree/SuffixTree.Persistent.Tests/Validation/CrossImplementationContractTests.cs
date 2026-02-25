using System;
using System.Collections.Generic;
using System.Linq;
using NUnit.Framework;
using SuffixTree;
using SuffixTree.Persistent;

namespace SuffixTree.Persistent.Tests.Validation;

[TestFixture]
[Category("Validation")]
public class CrossImplementationContractTests
{
    private static readonly string[] SampleTexts =
    {
        string.Empty,
        "a",
        "banana",
        "mississippi",
        "abracadabra",
        "ACGTACGTACGT",
        "🧬αβγ🧪$",
        "repeat-repeat-repeat"
    };

    [Test]
    [TestCaseSource(nameof(SampleTexts))]
    public void TextAndEmptyParity_WithInMemoryImplementation(string text)
    {
        var reference = global::SuffixTree.SuffixTree.Build(text);
        using var persistentDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable;
        var persistent = (ISuffixTree)persistentDisposable!;

        Assert.Multiple(() =>
        {
            Assert.That(persistent.Text.ToString(), Is.EqualTo(reference.Text.ToString()), "Text mismatch");
            Assert.That(persistent.IsEmpty, Is.EqualTo(reference.IsEmpty), "IsEmpty mismatch");
        });
    }

    [Test]
    [TestCaseSource(nameof(SampleTexts))]
    public void TraverseEventStream_Parity_WithInMemoryImplementation(string text)
    {
        var reference = global::SuffixTree.SuffixTree.Build(text);
        using var persistentDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable;
        var persistent = (ISuffixTree)persistentDisposable!;

        var referenceEvents = new List<string>();
        var persistentEvents = new List<string>();

        reference.Traverse(new EventCollector(referenceEvents));
        persistent.Traverse(new EventCollector(persistentEvents));

        Assert.That(persistentEvents, Is.EqualTo(referenceEvents), "Traverse event stream mismatch");
    }

    [Test]
    [TestCaseSource(nameof(SampleTexts))]
    public void MaxDepth_Parity_WithInMemoryImplementation(string text)
    {
        var reference = global::SuffixTree.SuffixTree.Build(text);
        using var persistentDisposable = PersistentSuffixTreeFactory.Create(new StringTextSource(text)) as IDisposable;
        var persistent = (ISuffixTree)persistentDisposable!;

        Assert.That(
            persistent.MaxDepth,
            Is.EqualTo(reference.MaxDepth),
            "Persistent MaxDepth must match in-memory ISuffixTree contract.");
    }

    private sealed class EventCollector : ISuffixTreeVisitor
    {
        private readonly List<string> _events;

        public EventCollector(List<string> events) => _events = events;

        public void VisitNode(int startIndex, int endIndex, int leafCount, int childCount, int depth)
            => _events.Add($"N({startIndex},{endIndex},{leafCount},{childCount},{depth})");

        public void EnterBranch(int key) => _events.Add($"E({key})");

        public void ExitBranch() => _events.Add("X");
    }
}
