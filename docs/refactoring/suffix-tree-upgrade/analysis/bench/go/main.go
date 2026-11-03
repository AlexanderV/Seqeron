// Go standard-library index/suffixarray benchmark harness (stdlib only).
// Same output protocol as the Seqeron harness: "R\top\tmetric\tvalue\tunit\tnotes".
//   gobench all   <text> <qprefix>     build + queries + 1/4-goroutine lookups
//   gobench write <text> <outfile>     build, serialize with Index.Write, report bytes/time
//   gobench load  <file> <qprefix>     suffixarray.Read from file (driver drops caches for cold) + first queries
package main

import (
	"bufio"
	"fmt"
	"index/suffixarray"
	"os"
	"runtime"
	"sort"
	"strconv"
	"strings"
	"sync"
	"time"
)

func emit(op, metric string, v float64, unit, notes string) {
	fmt.Printf("R\t%s\t%s\t%s\t%s\t%s\n", op, metric, strconv.FormatFloat(v, 'g', 17, 64), unit, notes)
}
func emitI(op, metric string, v int64, unit, notes string) {
	fmt.Printf("R\t%s\t%s\t%d\t%s\t%s\n", op, metric, v, unit, notes)
}

func lines(path string) [][]byte {
	b, err := os.ReadFile(path)
	if err != nil {
		panic(err)
	}
	var out [][]byte
	for _, l := range strings.Split(string(b), "\n") {
		if l != "" {
			out = append(out, []byte(l))
		}
	}
	return out
}

func statusKB(key string) int64 {
	b, _ := os.ReadFile("/proc/self/status")
	for _, l := range strings.Split(string(b), "\n") {
		if strings.HasPrefix(l, key) {
			f := strings.Fields(l)
			v, _ := strconv.ParseInt(f[1], 10, 64)
			return v * 1024
		}
	}
	return 0
}

var peakAnon, peakRss int64
var stop = make(chan struct{})

func sampler() {
	go func() {
		for {
			select {
			case <-stop:
				return
			default:
			}
			if v := statusKB("RssAnon:"); v > peakAnon {
				peakAnon = v
			}
			if v := statusKB("VmRSS:"); v > peakRss {
				peakRss = v
			}
			time.Sleep(20 * time.Millisecond)
		}
	}()
}

func procMem(op string) {
	hwm := statusKB("VmHWM:")
	a := statusKB("RssAnon:")
	if peakAnon > a {
		a = peakAnon
	}
	emitI(op, "peak_rss_hwm", hwm, "bytes", "VmHWM at end of process")
	emitI(op, "peak_rss_anon_sampled", a, "bytes", "max RssAnon sampled every 20ms")
	var m runtime.MemStats
	runtime.ReadMemStats(&m)
	emitI(op, "go_sys_bytes_end", int64(m.Sys), "bytes", "runtime.MemStats.Sys")
}

// passes: untimed warm pass, then repeat until >=1s.
func passes(body func() int64) (spp float64, n int, chk int64) {
	chk = body()
	t0 := time.Now()
	for {
		chk = body()
		n++
		if time.Since(t0).Seconds() >= 1.0 {
			break
		}
	}
	return time.Since(t0).Seconds() / float64(n), n, chk
}

func countPass(ix *suffixarray.Index, qs [][]byte) int64 {
	var s int64
	for _, q := range qs {
		s += int64(len(ix.Lookup(q, -1)))
	}
	return s
}
func containsPass(ix *suffixarray.Index, qs [][]byte) int64 {
	var s int64
	for _, q := range qs {
		if len(ix.Lookup(q, 1)) > 0 {
			s++
		}
	}
	return s
}

var hitsOut int64

func locatePass(ix *suffixarray.Index, qs [][]byte) int64 {
	var sum, hits int64
	for _, q := range qs {
		r := ix.Lookup(q, -1)
		hits += int64(len(r))
		for _, p := range r {
			sum += int64(p)
		}
	}
	hitsOut = hits
	return sum
}

func concurrent(ix *suffixarray.Index, qs [][]byte, threads int, op string, pass func(*suffixarray.Index, [][]byte) int64, unit string) {
	// every goroutine runs the FULL query set `reps` times (identical per-thread working set for 1 and 4 threads)
	t0 := time.Now()
	expect := pass(ix, qs)
	single := time.Since(t0).Seconds()
	reps := int(1.0/maxf(single, 1e-6)) + 1
	if reps > 400 {
		reps = 400
	}
	sums := make([]int64, threads)
	var wg, ready sync.WaitGroup
	gate := make(chan struct{})
	wg.Add(threads)
	ready.Add(threads)
	for i := 0; i < threads; i++ {
		go func(k int) {
			defer wg.Done()
			ready.Done()
			<-gate
			var s int64
			for r := 0; r < reps; r++ {
				s = pass(ix, qs)
			}
			sums[k] = s
		}(i)
	}
	ready.Wait()
	t1 := time.Now()
	close(gate)
	wg.Wait()
	el := time.Since(t1).Seconds()
	ok := int64(1)
	for _, s := range sums {
		if s != expect {
			ok = 0
		}
	}
	emit(op, "throughput", float64(threads)*float64(len(qs))*float64(reps)/el, unit, fmt.Sprintf("goroutines=%d reps=%d GOMAXPROCS=%d each goroutine runs the full query set; one shared Index", threads, reps, runtime.GOMAXPROCS(0)))
	emitI(op, "checksum_ok", ok, "bool", "every goroutine's result equals the single-thread result")
}

func maxf(a, b float64) float64 {
	if a > b {
		return a
	}
	return b
}

func querySuite(ix *suffixarray.Index, qp string) {
	present := lines(qp + ".present20.txt")
	absent := lines(qp + ".absent20.txt")
	loc := lines(qp + ".locate8.txt")

	rate := func(op string, qs [][]byte, body func() int64) {
		spp, n, chk := passes(body)
		emit(op, "throughput", float64(len(qs))/spp, "queries/s", fmt.Sprintf("passes=%d", n))
		emit(op, "time_per_pass", spp, "s", fmt.Sprintf("queries/pass=%d", len(qs)))
		emitI(op, "checksum", chk, "count", "sum of results (cross-tool oracle)")
	}
	rate("count_present20", present, func() int64 { return countPass(ix, present) })
	rate("count_absent20", absent, func() int64 { return countPass(ix, absent) })
	rate("contains_present20", present, func() int64 { return containsPass(ix, present) })
	rate("contains_absent20", absent, func() int64 { return containsPass(ix, absent) })
	locRate := func(op string, qs [][]byte) {
		spp, n, chk := passes(func() int64 { return locatePass(ix, qs) })
		hits := hitsOut
		emit(op, "throughput", float64(len(qs))/spp, "patterns/s", fmt.Sprintf("passes=%d", n))
		emit(op, "positions_per_s", float64(hits)/spp, "positions/s", fmt.Sprintf("total_hits_per_pass=%d", hits))
		emitI(op, "total_hits", hits, "count", "occurrences returned per pass (oracle)")
		emitI(op, "checksum", chk, "count", "sum of positions (oracle)")
	}
	locRate("locate8", loc)
	if _, err := os.Stat(qp + ".locate4.txt"); err == nil {
		locRate("locate4", lines(qp+".locate4.txt"))
	}
	// latency
	d := make([]int64, len(present))
	for i, q := range present {
		a := time.Now()
		_ = len(ix.Lookup(q, -1))
		d[i] = time.Since(a).Nanoseconds()
	}
	sort.Slice(d, func(i, j int) bool { return d[i] < d[j] })
	var sm float64
	for _, x := range d {
		sm += float64(x)
	}
	emit("count_present20", "latency_p50", float64(d[len(d)/2]), "ns", "per-query time.Now (includes ~40ns timer overhead)")
	emit("count_present20", "latency_p99", float64(d[int(float64(len(d))*0.99)]), "ns", "")
	emit("count_present20", "latency_mean", sm/float64(len(d)), "ns", "")
	for _, th := range []int{1, 4} {
		concurrent(ix, present, th, fmt.Sprintf("count_present20_%dt", th), countPass, "queries/s")
		concurrent(ix, loc, th, fmt.Sprintf("locate8_%dt", th), locatePass, "patterns/s")
	}
}

type counter struct{ n int64 }

func (c *counter) Write(p []byte) (int, error) { c.n += int64(len(p)); return len(p), nil }

func main() {
	if len(os.Args) < 4 {
		fmt.Fprintln(os.Stderr, "usage: gobench all|write|load ...")
		os.Exit(2)
	}
	sampler()
	switch os.Args[1] {
	case "all":
		data, err := os.ReadFile(os.Args[2])
		if err != nil {
			panic(err)
		}
		qp := os.Args[3]
		runtime.GC()
		var m0, m1 runtime.MemStats
		runtime.ReadMemStats(&m0)
		t0 := time.Now()
		ix := suffixarray.New(data)
		el := time.Since(t0).Seconds()
		hwm := statusKB("VmHWM:")
		emitI("build", "peak_rss_hwm_after_build", hwm, "bytes", "VmHWM right after build, before any query")
		emitI("build", "peak_rss_anon_after_build", maxi(peakAnon, statusKB("RssAnon:")), "bytes", "max RssAnon sampled up to end of build")
		emit("build", "time", el, "s", "suffixarray.New(data) wall time (SA-IS)")
		emit("build", "throughput", float64(len(data))/el, "chars/s", "")
		runtime.GC()
		runtime.ReadMemStats(&m1)
		delta := int64(m1.HeapAlloc) - int64(m0.HeapAlloc)
		emitI("build", "index_bytes_heap_delta", delta, "bytes", "HeapAlloc after-GC delta over Build (suffix array only; Index references the caller's text, 1 B/char, not counted)")
		emit("build", "bytes_per_char", float64(delta)/float64(len(data)), "bytes/char", "heap delta / n")
		emit("build", "bytes_per_char_incl_text", float64(delta+int64(len(data)))/float64(len(data)), "bytes/char", "heap delta + 1 B/char text (comparable to Seqeron file size which embeds text)")
		emitI("build", "total_alloc_during_build", int64(m1.TotalAlloc-m0.TotalAlloc), "bytes", "runtime.MemStats.TotalAlloc delta")
		querySuite(ix, qp)
		procMem("process")
	case "write":
		data, err := os.ReadFile(os.Args[2])
		if err != nil {
			panic(err)
		}
		t0 := time.Now()
		ix := suffixarray.New(data)
		bt := time.Since(t0).Seconds()
		emit("build", "time", bt, "s", "suffixarray.New(data) wall time (SA-IS)")
		cw := &counter{}
		t1 := time.Now()
		if err := ix.Write(cw); err != nil {
			panic(err)
		}
		emit("serialize", "write_time_counting", time.Since(t1).Seconds(), "s", "Index.Write to a byte counter")
		emitI("serialize", "serialized_bytes", cw.n, "bytes", "bytes written by Index.Write (text + varint-style SA)")
		emit("serialize", "serialized_bytes_per_char", float64(cw.n)/float64(len(data)), "bytes/char", "")
		f, _ := os.Create(os.Args[3])
		bw := bufio.NewWriterSize(f, 1<<20)
		t2 := time.Now()
		if err := ix.Write(bw); err != nil {
			panic(err)
		}
		bw.Flush()
		f.Close()
		emit("serialize", "write_time_file", time.Since(t2).Seconds(), "s", "Index.Write to file")
		procMem("process")
	case "load":
		qp := os.Args[3]
		present := lines(qp + ".present20.txt")
		f, err := os.Open(os.Args[2])
		if err != nil {
			panic(err)
		}
		t0 := time.Now()
		ix := suffixarray.New(nil)
		if err := ix.Read(bufio.NewReaderSize(f, 1<<20)); err != nil {
			panic(err)
		}
		emit("load", "time", time.Since(t0).Seconds(), "s", "Index.Read from file (deserialize; NOT mmap, copies into RAM)")
		t1 := time.Now()
		s := countPass(ix, present[:1000])
		emit("load", "first_1000_queries_time", time.Since(t1).Seconds(), "s", "")
		t1 = time.Now()
		s += countPass(ix, present)
		emit("load", "first_100k_queries_time", time.Since(t1).Seconds(), "s", "")
		emitI("load", "checksum", s, "count", "")
		procMem("process")
	}
	close(stop)
}

func maxi(a, b int64) int64 {
	if a > b {
		return a
	}
	return b
}
