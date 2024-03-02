#include <math.h>
#include <stdio.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "timer.h"

void main() {
    int i;
    int nthreads = 1, thread_num = 0;
    double sum = 0;
    double runtime;
    const long N = 1000000000L;

  
#ifdef _OPENMP
    #pragma omp parallel
    {
    nthreads = omp_get_num_threads();
    thread_num = omp_get_thread_num();
    if (thread_num == 0)
        printf("[OPENMP] Running %d threads\n", nthreads);
    }
#endif

    StartTimer();

    // Parallelize this loop with OpenMP
    #pragma omp parallel for reduction(+:sum)
    for (i = 0; i < N; i++) {
        sum += exp((i % 5) - 2 * (i % 7));
    }
    printf("sum = %g\n", sum);

    runtime = GetTimer();
    printf("total time %f s\n", runtime / 1000);

}