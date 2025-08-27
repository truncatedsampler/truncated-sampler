#ifndef BENCH_SIGN_H
#define BENCH_SIGN_H

#include "benchmark.h"

double keygen_bench(timer* t);
double sign_bench(timer* t);
double verify_valid_bench(timer* t);
double verify_invalid_bench(timer* t);

#endif /* BENCH_SIGN_H */

