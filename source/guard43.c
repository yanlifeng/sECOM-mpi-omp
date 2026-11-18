#define _GNU_SOURCE
#include <stdatomic.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#ifdef _OPENMP
#include <omp.h>
#include <pthread.h>
#endif
#ifdef __linux__
#include <sched.h>
#endif

/* Only BARO_SUBCYCLE owns this optional context. Each reservation ends
 * after END PARALLEL, never inside a team. No epoch resets or ABA. */
typedef struct { _Alignas(128) atomic_uint v; } slot;
typedef struct { slot lane[32], captain[4], gate; } domain;
static domain groups[2];
static atomic_flag taken=ATOMIC_FLAG_INIT;
#ifdef _OPENMP
static pthread_once_t once=PTHREAD_ONCE_INIT;
static int enabled, next_call=1;
static void config(void) {
    const char *e=getenv("SECOM_V43_HALO");
    enabled=e && strcmp(e,"1")==0 && atomic_is_lock_free(&groups[0].gate.v);
}
#endif
int claim43_(const int *steps) {
#ifdef _OPENMP
    if (*steps<=0 || *steps>100000 || omp_get_level()!=0 || omp_get_dynamic()) return 0;
    pthread_once(&once,config);
    if (!enabled || atomic_flag_test_and_set_explicit(&taken,memory_order_acquire)) return 0;
    if (next_call>700000000-*steps) {
        atomic_flag_clear_explicit(&taken,memory_order_release);return 0;
    }
    int first=next_call; next_call+=*steps; return first;
#else
    (void)steps; return 0;
#endif
}
void free43_(void) { atomic_flag_clear_explicit(&taken,memory_order_release); }
int place43_(void) {
#ifdef __linux__
    cpu_set_t set;int c=sched_getcpu();
    if(c>=0 && c<CPU_SETSIZE && !sched_getaffinity(0,sizeof set,&set) &&
       CPU_COUNT(&set)==1 && CPU_ISSET(c,&set)) return c;
#endif
    return -1;
}
int map43_(const int *map,const int *size) {
    int n=*size;
    if(n!=32 && n!=64) return 0;
    for(int i=0;i<n;i++) {
        if(map[i]<0) return 0;
        for(int j=0;j<i;j++) if(map[i]==map[j]) return 0;
    }
    return 1;
}
/* One writer per lane. Arrivals broadcast all preceding Fortran writes
 * transitively lane -> captain -> root -> every reader, incl. halos.
 * The next stage cannot finish before its slowest consumer has arrived. */
static void wait_slot(slot *s,unsigned e) {
    unsigned spins=0;struct timespec start,t;
    while(atomic_load_explicit(&s->v,memory_order_acquire)<e) {
        if (!(++spins & 8191)) {
            if(clock_gettime(CLOCK_MONOTONIC,&t)) _Exit(87);
            if(spins==8192) start=t;
            else if(t.tv_sec-start.tv_sec>15) _Exit(86);
        }
    }
}
void localsync43_(const int *tid,const int *ng,const int *epoch) {
    unsigned e=(unsigned)*epoch;
    int n=*ng,t=*tid;
    if(n<1 || n>32 || t<0 || t>=2*n || *epoch<=0) _Exit(84);
    domain *d=groups+t/n;int i=t%n;
    atomic_store_explicit(&d->lane[i].v,e,memory_order_release);
    if(i%8==0) {
        for(int j=i+1;j<n && j<i+8;j++) wait_slot(&d->lane[j],e);
        atomic_store_explicit(&d->captain[i/8].v,e,memory_order_release);
    }
    if(!i) {
        for(int j=1;j<(n+7)/8;j++) wait_slot(&d->captain[j],e);
        atomic_store_explicit(&d->gate.v,e,memory_order_release);
    } else wait_slot(&d->gate,e);
}
