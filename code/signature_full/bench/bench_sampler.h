#ifndef BENCH_SAMPLER_H
#define BENCH_SAMPLER_H

#include "benchmark.h"

double sample_perturb_bench(timer* t);
double sample_klein_bench(timer* t);
double elliptic_sampler_bench(timer* t);

#endif /* BENCH_SAMPLER_H */

