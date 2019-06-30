#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <math.h>
typedef long long ll;
#define realPI acos(-1)
const ll num_steps = 1e9;
const int NUM_THREADS = 20;

// 使用for制导指令
int main()
{
    puts("using '#pragma omp for'...");
    printf("PI(%.20f) in step_num=%lld...\n", realPI, num_steps);
    ll i;
    double step, start_time, end_time, pi = 0.0, sum[NUM_THREADS];
    step = 1.0 / (double)num_steps;
    omp_set_num_threads(NUM_THREADS);
    start_time = omp_get_wtime();
#pragma omp parallel private(i)
    {
        int id = omp_get_thread_num();
        sum[id] = 0.0;
        double x;
#pragma omp for schedule(static,4)
        for (i = 0; i < num_steps; i++)
        {
            x = (i + 0.5) * step;
            sum[id] += 4.0 / (1.0 + x * x);
        }
    }
    for (i = 0; i < NUM_THREADS; i++)
        pi += sum[i] * step;
    end_time = omp_get_wtime();
    printf("Pi=%.20f  Running time=%f s\n", pi, end_time - start_time);
    return 0;
}
