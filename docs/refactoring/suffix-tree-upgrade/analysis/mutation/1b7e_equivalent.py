# Mutants of the audit of 1b.7e that no test kills, in two groups (id -> argument). Read by summarize_1b7e.py and quoted in the ADR subsection.
# EQUIVALENT: the mutant behaves like the original in every state that a caller or a test can reach (the argument says why); DEFENCE: a second line of defence behind a check that is tested at another level, a test hook that no
# fault of the code under test can make observable, or an effect that no assertion of the language runtime can see (a native leak, an fsync): not equivalent, not killed, reported as such. Ids with the prefix d_ are the mutants
# of the audit of 1b.7d that the tests of 1b.7e were to close (SesaGapTests).
EQUIVALENT = {
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
    'W01': 'a negative offset reaches RandomAccess.Read, which refuses it with the same exception type; the BWT file made before is removed by the failure path (tested: no BWT file after any failure)',
    'K35': 'the tagged codec is refused at the entry of the stage (the container has no id for it), so `options.Codec` is always `Packed` where it is read',
    'K25': 'in container mode PH4 creates no file of its own (the placer gives it the sections), so the file options that it would pass to the creation are not read; its scratch files are made by DistScratch, which does not take them',
    'K29': 'the scratch directories ph4 and ph5 are both created and both removed by the stage; PH5 and PH4 use distinct file names, so a shared directory holds the same files (the mutant of PH4 using ph5, K15, is killed by the state of ph4 at B3)',
}
DEFENCE = {
    'V15': 'Complete of the base file in path mode only checks that no writer is open (a file is not deleted by Dispose and `fsync: false` flushes nothing): the writers were disposed just before and the owner of the file checks again',
    'W37': 'the same as V15 for the BWT file (the Psi stream opens it by path afterwards and reads it whole)',
    'P63': 'a native block of 8 x sigma bytes per pass is not freed: nothing in the language runtime or in the tests can see a native leak of this size (the arena has guards and counters, these raw blocks do not)',
    'd_Q76': 'ArenaGuards are a test hook that adds guard bands and a poison pattern: a stage that ignores them produces the same bytes and the same counters (as in 1b.7d); they only matter when the stage overruns a buffer, which no injectable fault does',
    'B183': 'the second fsync of the publication has no observable effect on this device, and the counters do not count fsyncs (CountedFile.Flush is not counted): the order of header, fsync, magic, fsync is held by the hook points and tested, the kernel call is not',
    'K17': 'the Threads of PH4 limit the number of B2 tasks; the output and every counter are the same at any thread count (tested), only the time differs',
    'K23': 'ArenaGuards of PH4: as d_Q76',
    'K37': 'ArenaGuards of PH5: as d_Q76',
    'K33': 'PH5 flushes its windows when it is done (`FlushWindows`): the flush is a msync that no read of the file can tell from the kernel writing the pages back',
    'B53': 'see the rerun: the seed changes the sampling of the splitters of PH2; if the pinned calls do not see it, the output does not depend on it (tested at five seeds)',
}
