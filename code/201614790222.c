#include <stdio.h>
#include <time.h>
#include <omp.h>
#include <math.h>
typedef long long ll;
#define realPI acos(-1)
const ll num_steps = 1e9;
const int NUM_THREADS = 20;

//单线程计算Pi
void single()
{
    puts("serial...");
    double step, x, pi, sum = 0.0, timetot;
    step = 1.0 / (double)num_steps;
    clock_t start, end;
    start = clock();
    ll i;
    for (i = 0; i < num_steps; i++)
    {
        x = (i + 0.5) * step;
        sum += 4.0 / (1.0 + x * x);
    }
    pi = step * sum;
    end = clock();
    timetot = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Pi=%.20f  Running time=%f s\n", pi, timetot);
    puts("");
}

// 并行域并行（SPMD并行模式）
// void SPMD()
// {
//     puts("SPMD mode...");
//     ll i;
//     double step, start_time, end_time, pi = 0.0, sum[NUM_THREADS];
//     step = 1.0 / (double)num_steps;
//     omp_set_num_threads(NUM_THREADS);
//     start_time = omp_get_wtime();
// #pragma omp parallel
//     {
//         double x;
//         int id = omp_get_thread_num();
//         for (i = id, sum[id] = 0.0; i < num_steps; i = i + NUM_THREADS)
//         {
//             x = (i + 0.5) * step;
//             sum[id] += 4.0 / (1.0 + x * x);
//         }
//     }
//     for (i = 0; i < NUM_THREADS; i++)
//         pi += sum[i] * step;
//     end_time = omp_get_wtime();
//     printf("Pi=%.20f  Running time=%f s\n", pi, end_time - start_time);
//     puts("");
// }

// 使用for制导指令
void withfor()
{
    puts("with '#pragma omp for'...");
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
#pragma omp for schedule(guided)
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
    puts("");
}

// 使用reduction子句
void reduction()
{
    puts("with '#pragma omp parallel for reduction(+:sum)'...");
    ll i;
    double x, sum = 0.0, step, start_time, end_time, pi = 0.0;
    step = 1.0 / (double)num_steps;
    omp_set_num_threads(NUM_THREADS);
    start_time = omp_get_wtime();
#pragma omp parallel private(i, x)
    {
#pragma omp for reduction(+ \
                          : sum) schedule(guided)
        for (i = 0; i < num_steps; i = i + 1)
        {
            x = (i + 0.5) * step;
            sum = sum + 4.0 / (1.0 + x * x);
        }
    }
    pi = step * sum;
    end_time = omp_get_wtime();
    printf("Pi=%.20f  Running time=%f s\n", pi, end_time - start_time);
    puts("");
}

// 使用critical子句
void critical()
{
    puts("with '#pragma omp critical'...");
    ll i;
    double x, sum = 0.0, step, start_time, end_time, pi = 0.0, aux;
    step = 1.0 / (double)num_steps;
    omp_set_num_threads(NUM_THREADS);
    start_time = omp_get_wtime();
#pragma omp parallel private(i, x, aux) shared(sum)
    {
#pragma omp for schedule(guided)
        for (i = 0; i < num_steps; ++ i)
        {
            x = (i + 0.5) * step;
            aux = 4.0 / (1.0 + x * x);
#pragma omp critical
            sum = sum + aux;
        }
    }
    pi = step * sum;
    end_time = omp_get_wtime();
    printf("Pi=%.20f  Running time=%f s\n", pi, end_time - start_time);
    puts("");
}

int main()
{
    printf("Caculating PI(%.20f) in %lld iterations with %d threads...\n", realPI, num_steps, NUM_THREADS);
    puts("");
    single();
    // SPMD();
    withfor();
    reduction();
    critical();
    return 0;
}
