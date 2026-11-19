#define _GNU_SOURCE
#include <sched.h>
#ifndef SECOM_TEST_CPU_BASE
#define SECOM_TEST_CPU_BASE 24
#endif
/* Test-only, not part of the Fortran runtime or solver build. */
int place(int rank,int lane,int nt) {
    int want=SECOM_TEST_CPU_BASE+2*rank+(nt>1?lane:0);
    cpu_set_t m;
    if(sched_getaffinity(0,sizeof m,&m) || CPU_COUNT(&m)!=1 ||
       !CPU_ISSET(want,&m) || sched_getcpu()!=want) return -1;
    return want;
}
