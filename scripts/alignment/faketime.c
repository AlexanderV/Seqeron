/* LD_PRELOAD shim: time() returns $FAKE_EPOCH (seconds since 1970, UTC) so that tools printing the run date can be run at chosen dates.
   gettimeofday() (the seed of ajRandomSeed of the EMBOSS library: tv_usec % 9999 + 1) returns $FAKE_EPOCH seconds and $FAKE_USEC microseconds when FAKE_USEC is set; without it the real clock. */
#define _GNU_SOURCE
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <dlfcn.h>
time_t time(time_t *t)
{
    const char *e = getenv("FAKE_EPOCH");
    time_t v = e ? (time_t) atoll(e) : 0;
    if (t) *t = v;
    return v;
}
int gettimeofday(struct timeval *tv, void *tz)
{
    const char *u = getenv("FAKE_USEC");
    if (u)
    {
        const char *e = getenv("FAKE_EPOCH");
        if (tv) { tv->tv_sec = e ? (time_t) atoll(e) : 0; tv->tv_usec = atol(u); }
        return 0;
    }
    int (*real)(struct timeval *, void *) = (int (*)(struct timeval *, void *)) dlsym(RTLD_NEXT, "gettimeofday");
    return real(tv, tz);
}
