#include "bench_sign.h"

#include "sign.h"
#include "randombytes.h"
#include "random.h"

double keygen_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  keys_init(&pk, &sk);
  keygen(&pk, &sk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  keys_clear(&pk, &sk);
  return time;
}

double sign_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];

  keys_init(&pk, &sk);
  sig_init(&sig);
  keygen(&pk, &sk);
  randombytes(msg, PARAM_N/8);
  randombytes(state, STATE_BYTES);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  sign(&sig, state, &sk, &pk, msg);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!verify(&sig, msg, &pk)) {
    printf("FATAL ERROR: benchmarked signature is not valid\n");
  }
  sig_clear(&sig);
  keys_clear(&pk, &sk);
  return time;
}

double verify_valid_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];

  keys_init(&pk, &sk);
  sig_init(&sig);
  keygen(&pk, &sk);
  randombytes(msg, PARAM_N/8);
  randombytes(state, STATE_BYTES);
  sign(&sig, state, &sk, &pk, msg);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = verify(&sig, msg, &pk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (!is_valid) {
    printf("FATAL ERROR: benchmarked signature is not valid\n");
  }
  sig_clear(&sig);
  keys_clear(&pk, &sk);
  return time;
}

double verify_invalid_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];

  keys_init(&pk, &sk);
  sig_init(&sig);
  keygen(&pk, &sk);
  randombytes(msg, PARAM_N/8);
  randombytes(state, STATE_BYTES);
  sign(&sig, state, &sk, &pk, msg);
  msg[0] ^= 1;
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  int is_valid = verify(&sig, msg, &pk);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  if (is_valid) {
    printf("FATAL ERROR: benchmarked signature is valid for the wrong message\n");
  }
  sig_clear(&sig);
  keys_clear(&pk, &sk);
  return time;
}