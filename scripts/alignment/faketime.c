/* LD_PRELOAD shim: time() returns $FAKE_EPOCH (seconds since 1970, UTC) so that tools printing the run date can be run at chosen dates. */
#include <stdlib.h>
#include <time.h>
time_t time(time_t *t)
{
    const char *e = getenv("FAKE_EPOCH");
    time_t v = e ? (time_t) atoll(e) : 0;
    if (t) *t = v;
    return v;
}
