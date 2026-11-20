# Mutants of the audit of 1b.7e that no test kills, in two groups (id -> argument). Read by summarize_1b7e.py and quoted in the ADR subsection.
# EQUIVALENT: the mutant behaves like the original in every state that a caller or a test can reach (the argument says why); DEFENCE: a second line of defence behind a check that is tested at another level, a test hook that no
# fault of the code under test can make observable, or an effect that no assertion of the language runtime can see (a native leak, an fsync): not equivalent, not killed, reported as such. Ids with the prefix d_ are the mutants
# of the audit of 1b.7d that the tests of 1b.7e were to close (SesaGapTests).
EQUIVALENT = {
    'A01': 'a null file falls into the path mode of BuildCore, where MappedOutputFile.Create refuses the null path with the same ArgumentNullException before any output exists; the scratch directories made in between are removed by the failure path',
    'P38': 'the read buffer of the Psi pass is filled by `ReadExactAt` over its whole length before any byte of it is used (zeroing it first changes no byte that is read)',
    'K49': 'every file under ph4 and ph5 is removed by the component that made it (DistScratch, SpillStack) on success and on failure, so both directories are empty when the stage deletes them; the recursive delete only matters for a file that no component accounts for (the failure test with a hook in B2 and B3 leaves nothing)',
    'B169': 'the discard `_ = infoLength;` only silences an unused parameter; removing it changes nothing at run time',
    'B173': 'the deepest repeated substring has the depth of the largest LCP value (LRS_DEPTH = MAX_LCP): the same number in both fields, 0 for a text without a repeat',
    'B187': 'the rename has moved the temporary file, so the cleanup that follows (File.Delete of a path that is gone) deletes nothing; with or without the flag the final state is the same',
    'B189': 'the same as B187: the cleanup deletes the temporary path, which no longer exists after the rename',
    'B241': 'the text is made of integers and enum names only: no culture formats them differently (no group separator, no negative number occurs)',
    'B245': 'the check refuses a big-endian host; the machine is little endian, and no other host is available to the tests (the check is an assertion about the platform, listed under "not verified")',
    'B117': 'BUILD_INFO is a constant text of about 170 bytes with a few numbers: it cannot reach the reserved 4,096 bytes, so the comparison at equality is never evaluated',
    'B113': 'the counts of the BWT are compared with those of PH0 in the line before and are equal when PSI is built',
    'B109': 'the BWT buffers are released before PSI starts and the cap of the windows is the same cap (`WindowCapBytes`): with or without the shared budget object PSI meets the same limit; the parameter is tested at the stage (ThePh6Buffers_AreChargedToTheSharedBudget...)',
    'B39': 'the window of `PackedTextWriterOptions` is read only by the overload that opens its own file; the builder passes a writer of the container, whose window is the one tested (EveryOption_ReachesItsStage)',
    'B37': 'the `Fsync` of `PackedTextWriterOptions` is read only by the overload that completes its own file (PackedTextWriter.cs line 147); the builder passes a writer and completes the container itself',
    'P45': 'a primary row at the first row after a block: the first call scans the whole block and the second scans a negative length (nothing); the next block starts at the primary row and skips it: the same rows are distributed',
    'A09': 'in section mode `owned` is null and `section` the file; in path mode `section` is null and `owned` the file: `section ?? owned` and `owned ?? section` are the same object in both modes',
    'G23': 'Grow requires that no writer is open and is called by the one thread that owns the file (checked by the open-writer test): nothing else takes the gate while it runs',
    'W27': 'row 0 never holds the suffix 0 for n > 0 (SA[0] = n is checked on the same pass), and for n = 0 there is one row: a second zero can only occur at a row above 0, so `primary > 0` and `primary >= 0` agree',
    'C03': 'the hardware path and the table path give the same value (both are tested against the bitwise definition and against each other on every section size); only the speed differs',
    'C15': 'the refusal fires only on a CPU without the instruction; the machine has it, and the path is reached only through `HardwareAccelerated` (a CPU without it takes the table path): the branch cannot be taken here (listed under "not verified": no other CPU)',
    'K01': 'the same check is made by LcpBuild a line later (`ThrowIfNegative(options.SaOffset)`, mutant V01, killed), with the same exception type and before any file of PH4 exists; the scratch directories made in between are removed by the stage',
    'K03': 'a negative cap: the constructor of ArenaBudget refuses it the line after with the same exception type, before any file is created',
    'K35': 'the tagged codec is refused at the entry of the stage (the container has no id for it), so `options.Codec` is always `Packed` where it is read',
    'K25': 'in container mode PH4 creates no file of its own (the placer gives it the sections), so the file options that it would pass to the creation are not read; its scratch files are made by DistScratch, which does not take them',
    'K29': 'the scratch directories ph4 and ph5 are both created and both removed by the stage; PH5 and PH4 use distinct file names, so a shared directory holds the same files (the mutant of PH4 using ph5, K15, is killed by the state of ph4 at B3)',
}
DEFENCE = {
    'P64': 'a native block of 8 x sigma bytes per pass is not freed: as P63',
    'B182': 'the first fsync of the publication: as B183',
    'V11': 'the output writers are created from two ranges of the container and one allocation of the arena that no cap can fail without B2, which holds the peak of PH4, failing first; no other fault is injectable there, so the cleanup of a failed constructor of the output writer is not reached by any test (the cleanup of a failed second writer is, V10)',
    'V15': 'Complete of the base file in path mode only checks that no writer is open (a file is not deleted by Dispose and `fsync: false` flushes nothing): the writers were disposed just before and the owner of the file checks again',
    'W37': 'the same as V15 for the BWT file (the Psi stream opens it by path afterwards and reads it whole)',
    'P63': 'a native block of 8 x sigma bytes per pass is not freed: nothing in the language runtime or in the tests can see a native leak of this size (the arena has guards and counters, these raw blocks do not)',
    'd_Q76': 'ArenaGuards are a test hook that adds guard bands and a poison pattern: a stage that ignores them produces the same bytes and the same counters (as in 1b.7d); they only matter when the stage overruns a buffer, which no injectable fault does',
    'B183': 'the second fsync of the publication has no observable effect on this device, and the counters do not count fsyncs (CountedFile.Flush is not counted): the order of header, fsync, magic, fsync is held by the hook points and tested, the kernel call is not',
    'K17': 'the Threads of PH4 limit the number of B2 tasks; the output and every counter are the same at any thread count (tested), only the time differs',
    'K23': 'ArenaGuards of PH4: as d_Q76',
    'K37': 'ArenaGuards of PH5: as d_Q76',
    'K33': 'PH5 flushes its windows when it is done (`FlushWindows`): the flush is a msync that no read of the file can tell from the kernel writing the pages back',
}

# Work package 1b.9c (2026-11-20): the mutants of the lists that the audit did not run, and the "second line of defence" mutants above, closed by an argued equivalence (a reviewer can check each argument against the code at the
# head named in it). The mutants of the second group that a new test kills are in 1b9c_closed.py; summarize_1b9c.py reads both. Every one of these survived the whole default set of SuffixTree.Esa.Tests (stage 2 of the runner).
EQUIVALENT_1B9C = {
    'B108': 'SesaBuilder.RunPh6: PsiOptions.WindowCapBytes and PsiOptions.Budget carry the same cap (Ph6CapBytes) and PsiStream charges the same bound to both; the shared budget holds nothing at that point (the BWT pass has disposed its arena), so both refuse exactly the same states; the option itself is tested at the stage (Ph6Tests.TheWindowsOfThePsiStream_AreChargedToTheCap_AndGivenBack)',
    'B120': 'SesaBuilder.BuildCore: Complete(fsync: false) of the container only asserts that no writer is open (every writer of PH0-PH6 was disposed by its using before) and that the file is not disposed; it is tested in isolation (MappedOutputFileTests.Complete_NeedsEveryWriterToBeDisposed_AndMayBeRepeated_AndSyncsOrNot); the container is closed, checksummed and renamed right after, so any writer left open would show there',
    'B122': 'SesaBuilder.BuildCore: `text = null` only keeps the `finally` from disposing the reader a second time, and MappedTextReader.Dispose is idempotent (MappedTextReaderTests.AfterDispose_EveryReadThrowsObjectDisposed_NothingIsRead_AndDisposeIsIdempotent)',
    'C08': 'Crc32C.Compute (table path): a text of a multiple of 8 bytes leaves its last 8 bytes to the byte-wise tail loop instead of the word loop; the tail loop is the definition of the recurrence, so the value is the same (the comment of the mutant says it: "same value"); only the speed of 8 bytes differs',
    'K17': 'SectionLcpChildBuilder: LcpBuild reads ExtLcpOptions.Threads only when Workers is 0 (LcpBuild.RunWorkers), and the section builder always passes Workers = p.LcpWorkers, which LcpChildRules.Resolve refuses below 1 (LcpChildBuilder.cs, ThrowIfLessThan(p.LcpWorkers, 1)); the thread count reaches PH4 through p.LcpWorkers, which is clamped to options.Threads and is tested (K18, PH4 gets one worker, is killed)',
    'V15': 'LcpBuild.EmitRankChunks (path mode): Complete(fsync: false) of the base file only asserts that no writer is open, and the writers were disposed by the using of the output writer before; the assertion is tested in isolation (MappedOutputFileTests.Complete_NeedsEveryWriterToBeDisposed_AndMayBeRepeated_AndSyncsOrNot)',
    'W37': 'BwtPass.Build: the same as V15 for the BWT file (the writers are disposed before, and the Psi stream opens the file by path and reads it whole, so a writer left open would not matter to it either)',
    'W08': 'BwtPass.Build: the check "text.Bits > MaxBits" and its mutant "> MaxBits + 1" differ only for a text of exactly 9 bits, and no reader has 9 bits: MappedTextReader.CheckShape accepts 1 to 8 bits (a packed text) or 32 (the reduced string) and refuses every other number (SaGap1b9cTests.ATextInMemory_WithNoLengthOfATextOrNoWidthOfTheCodec_..., which tries 0, -1, 9 and 33)',
    'W14': 'BwtPass.Build: zero: true instead of zero: false only changes the content of the buffer of the suffix array entries before it is read: every entry that the pass uses (saBuf.AsSpan(take)) is written by ReadExactAt, which throws unless it filled the whole span, before the loop reads it, and nothing else reads the buffer',
    'V16': 'LcpBuild.EmitRankChunks (path mode): Complete(fsync: false) of the overflow file only asserts that no writer of the file is open and that it is not disposed, and "using (output)" has disposed both writers when the call is reached; the assertion itself is tested in isolation (MappedOutputFileTests.Complete_NeedsEveryWriterToBeDisposed_AndMayBeRepeated_AndSyncsOrNot); the same argument as L17',
}
