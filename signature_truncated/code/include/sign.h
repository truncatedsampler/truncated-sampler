#ifndef SIGN_H
#define SIGN_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

typedef struct {
  poly_q_mat_d_d R[2][PARAM_K_L];
  poly_real_mat_2d_2d S;
} sk_t;

typedef struct {
  poly_q_mat_d_d B[PARAM_K_L];
  uint8_t seed[SEED_BYTES];
} pk_t;

typedef struct {
  poly_q tag;
  poly_q_vec_d v12;
  poly_q_vec_d v2[PARAM_K_L];
  poly_q_vec_k_l v3;
} sig_t;

void keys_init(pk_t *pk, sk_t *sk);
void keys_clear(pk_t *pk, sk_t *sk);
void sig_init(sig_t *sig);
void sig_clear(sig_t *sig);

void keygen(pk_t *pk, sk_t *sk);
void tag_gen(poly_q tag, uint8_t state[STATE_BYTES]);
void sign(sig_t *sig, uint8_t state[STATE_BYTES], const sk_t *sk, const pk_t *pk, const uint8_t msg[PARAM_N/8]);
int verify(const sig_t *sig, const uint8_t msg[PARAM_N/8], const pk_t *pk);

#endif /* SIGN_H */
