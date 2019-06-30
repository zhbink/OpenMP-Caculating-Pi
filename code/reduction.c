#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <math.h>
typedef long long ll;
#define realPI acos(-1)
const ll num_steps = 1e7;
const int NUM_THREADS = 20;

// 使用reduction子句
int main()
{
    puts("using '#pragma omp parallel for reduction(+:sum)'...");
    printf("PI(%.20f) in step_num=%lld...\n", realPI, num_steps);
    ll i;
    double x, sum = 0.0, step, start_time, end_time, pi = 0.0;
    step = 1.0 / (double)num_steps;
    omp_set_num_threads(NUM_THREADS);
    start_time = omp_get_wtime();
#pragma omp parallel private(i, x)
    {
#pragma omp for reduction(+ \
                          : sum) schedule(static)
        for (i = 0; i < num_steps; i = i + 1)
        {
            x = (i + 0.5) * step;
            sum = sum + 4.0 / (1.0 + x * x);
        }
    }
    pi = step * sum;
    end_time = omp_get_wtime();
    printf("Pi=%.20f  Running time=%f s\n", pi, end_time - start_time);
    return 0;
}


