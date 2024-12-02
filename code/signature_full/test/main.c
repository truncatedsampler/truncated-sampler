#include <stdio.h>
#include "poly_q_sampling.h"
#include "arith.h"
#include "sign.h"
#include "covariance.h"
#include "randombytes.h"
#include "random.h"

#define NTESTS 1
#define NSUBTESTS 5

static int keygen_test(void)
{
  int rval = 1;
  sk_t sk;
  pk_t pk;
  poly_q_mat_d_d A, RRstar[2][2];

  printf("\nkeygen_test\n");

  keys_init(&pk, &sk);
  poly_q_mat_d_d_init(A);
  poly_q_mat_d_d_init(RRstar[0][0]);
  poly_q_mat_d_d_init(RRstar[0][1]);
  poly_q_mat_d_d_init(RRstar[1][0]);
  poly_q_mat_d_d_init(RRstar[1][1]);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    keygen(&pk, &sk);
    if (sk_sq_spectral_norm(RRstar, sk.R)) {
      printf("keygen generated keys with large spectral norm.\n");
      rval = 0;
      goto keygen_test_cleanup;
    }

    poly_q_mat_d_d_uniform(A, pk.seed, DOMAIN_SEPARATOR_A, 0);
    for (size_t j = 0; j < PARAM_K; j++) {
      poly_q_mat_d_d_mul_mat_d_d(RRstar[0][0], A, sk.R[1][j]);
      poly_q_mat_d_d_add(RRstar[0][0], RRstar[0][0], sk.R[0][j]);
      if (!poly_q_mat_d_d_equal(RRstar[0][0], pk.B[j])) {
        printf("keygen generated keys with B != AR.\n");
        rval = 0;
        goto keygen_test_cleanup;
      }
    }

    poly_q_set_coeff(sk.R[1][0]->rows[0]->entries[0], PARAM_N/2, 3); // changing one coefficient of R
    poly_q_mat_d_d_mul_mat_d_d(RRstar[0][0], A, sk.R[1][0]);
    poly_q_mat_d_d_add(RRstar[0][0], RRstar[0][0], sk.R[0][0]);
    if (poly_q_mat_d_d_equal(RRstar[0][0], pk.B[0])) {
      printf("found R' such that AR' = B.\n");
      rval = 0;
      goto keygen_test_cleanup;
    }

    printf(":");
    fflush(stdout);
  }

keygen_test_cleanup:
  keys_clear(&pk, &sk);
  poly_q_mat_d_d_clear(A);
  poly_q_mat_d_d_clear(RRstar[0][0]);
  poly_q_mat_d_d_clear(RRstar[0][1]);
  poly_q_mat_d_d_clear(RRstar[1][0]);
  poly_q_mat_d_d_clear(RRstar[1][1]);

  return rval;
}

static int sig_test(void)
{
  int rval = 1;
  sk_t sk;
  pk_t pk;
  sig_t sig;
  uint8_t state[STATE_BYTES], msg[PARAM_N/8];
  randombytes(state, STATE_BYTES);

  printf("\nsig_test\n");

  keys_init(&pk, &sk);
  sig_init(&sig);

  for (int i = 0; i < NSUBTESTS; i++)
  {
    keygen(&pk, &sk);
    for (int j = 0; j < NSUBTESTS; j++)
    {
      randombytes(msg, PARAM_N/8);
      sign(&sig, state, &sk, &pk, msg);
      if (!verify(&sig, msg, &pk))
      {
        printf("verify returned zero for a valid signature.\n");
        rval = 0;
        goto sig_test_cleanup;
      }
      msg[0] ^= 1;
      if (verify(&sig, msg, &pk))
      {
        printf("verify returned non-zero for a valid signature on the wrong message.\n");
        rval = 0;
        goto sig_test_cleanup;
      }
      printf(":");
      fflush(stdout);
    }
  }

sig_test_cleanup:
  sig_clear(&sig);
  keys_clear(&pk, &sk);
  return rval;
}

int main(void) {
  int pass = 1;
  arith_setup();
  random_init();
  printf("Hello from the unit tests.\n");
  for (int i = 0; i < NTESTS; i++)
  {
    pass &= keygen_test();
    pass &= sig_test();

    if (!pass)
    {
      printf("FAILED!\n");
      break;
    } else {
      printf(".");
    }
  }
  if (pass)
  {
    printf("\npassed.\n");
  }
  arith_teardown();
  return 0;
}
