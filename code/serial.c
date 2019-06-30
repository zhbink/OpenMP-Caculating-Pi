#include <stdio.h>
#include <time.h>
#include <math.h>
typedef long long ll;
#define realPI acos(-1)
const ll num_steps = 1e9;

int main()
{
    printf("PI(%.20f) in step_num=%lld...\n", realPI, num_steps);
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
    return 0;
}
