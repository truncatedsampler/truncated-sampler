#include <stdio.h>
#include <time.h>

#include "sign.h"
#include "randombytes.h"
#include "random.h"

#include "bench_sign.h"
#include "bench_sampler.h"

#define BENCH_ITERATIONS 100

int main(void) {
  arith_setup();
  random_init();

  printf("[+] Running ./build/bench with N = %d iterations\n\n", BENCH_ITERATIONS);
  printf("______________________________________________________________________________________________________________________________________________________________\n");

  printf("%-30s\t%26s%4s%27s\t%27s%6s%28s\n", "", "", "Time", "", "", "Cycles", "");
  printf("%-30s\t%57s\t%61s\n", "", "---------------------------------------------------------", "-------------------------------------------------------------");
  printf("%-30s\t%9s\t%9s\t%9s\t%9s\t%13s\t%13s\t%13s\t%13s\n", "Benchmarked Functionality", "mean (ms)", "med (ms)", "min (ms)", "max (ms)", "mean (cycles)", "med (cycles)", "min (cycles)", "max (cycles)");
  printf("______________________________________________________________________________________________________________________________________________________________\n");

  printf("\n");

  benchmark("KEYGEN", BENCH_ITERATIONS, keygen_bench);
  benchmark("SIGN", BENCH_ITERATIONS, sign_bench);
  benchmark("VERIFY (Valid)", BENCH_ITERATIONS, verify_valid_bench);
  benchmark("VERIFY (Invalid)", BENCH_ITERATIONS, verify_invalid_bench);
  printf("\n");
  benchmark("SAMPLE_PERTURB", BENCH_ITERATIONS, sample_perturb_bench);
  benchmark("SAMPLE_KLEIN", BENCH_ITERATIONS, sample_klein_bench);
  benchmark("TRUNCATED_SAMPLER", BENCH_ITERATIONS, truncated_sampler_bench);

  printf("______________________________________________________________________________________________________________________________________________________________\n");


  arith_teardown();
  return 0;
}


