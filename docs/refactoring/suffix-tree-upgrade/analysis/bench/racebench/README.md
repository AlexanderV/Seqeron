# racebench: the Dispose-race client

A process of this client runs many rounds of the race that the lifetime contract of the facade has to survive: readers call a tree while it is disposed after a random delay (ADR-004 section 5a, "1b.8", "1b.9b", "1b.9g").
It is outside the shipped code and not part of `Seqeron.sln`; it uses the public API of `SuffixTree.Esa` only and takes the build under test from a directory, so that one source compares the head with a change:

```
dotnet build RaceBench.csproj -c Release -p:EsaDir=<dir with SuffixTree.Esa.dll, SuffixTree.Core.dll, SuffixTree.Persistent.dll> -o <out dir>
dotnet <out dir>/RaceBench.dll mem  <rounds> <seed>
dotnet <out dir>/RaceBench.dll file <rounds> <seed> <container> <text file>
```

* `mem`: an in-memory tree; four readers walk suffix links (the longest common substring, the maximal exact matches, ranks, counts, a traversal); a thread forces compacting collections and refills the freed memory; the main thread
  disposes after a random 0-2 ms. The defect it finds (1b.9b): `Dispose` that frees the pins of the navigation pack under a reader (the test host dies).
* `file`: a tree over a container (any container of the text file; `Header` and `Structure` levels, `Checked` and `Trusted` access in turn); four readers run every kind of member (count, locate, contains, the longest common
  substring, the maximal exact matches, a traversal, the suffixes, the suffix array, the text); the main thread disposes after a random 0-40 ms. The defect it finds (1b.8): a read of the unmapped view (an access violation).

Every answer is compared with an in-memory tree made by the client. A call must complete with the right answer or throw `ObjectDisposedException`. Output: one `RESULT` line; the exit code is 0 for a clean process, 1 for a wrong answer
or another exception, and a crash of the process (a signal) or a hang (a reader that does not finish within 120 s ends the process with `FailFast`) is a defect of the build. `run_race.sh` runs a number of processes and counts the three outcomes.
