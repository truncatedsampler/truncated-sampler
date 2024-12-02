#include "bench_sampler.h"

#include "sign.h"
#include "randombytes.h"
#include "random.h"
#include "arith.h"
#include "sampling.h"
#include "poly_q_sampling.h"

double sample_perturb_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  poly_q tag, tag_star;
  poly_real m_sG_t_tstar;
  poly_q_vec_d pL[PARAM_L], p12, p2[PARAM_K_L];
  uint8_t state[STATE_BYTES];

  keys_init(&pk, &sk);
  poly_q_init(tag);
  poly_q_init(tag_star);
  poly_real_init(m_sG_t_tstar);
  for (size_t i = 0; i < PARAM_L; i++) {
    poly_q_vec_d_init(pL[i]);
  }
  poly_q_vec_d_init(p12);
  for (size_t i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_init(p2[i]);
  }

  keygen(&pk, &sk);
  randombytes(state, STATE_BYTES);
  tag_gen(tag, state);
  poly_q_conjugate(tag_star, tag);
  poly_q_mul(tag, tag, tag_star);
  poly_real_from_poly_q(m_sG_t_tstar, tag);
  poly_real_mul_scalar(m_sG_t_tstar, m_sG_t_tstar, (coeff_real)PARAM_NEG_SGSQ);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  sample_perturb(pL, p12, p2, m_sG_t_tstar, sk.R, sk.S);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  keys_clear(&pk, &sk);
  poly_q_clear(tag);
  poly_q_clear(tag_star);
  poly_real_clear(m_sG_t_tstar);
  for (size_t i = 0; i < PARAM_L; i++) {
    poly_q_vec_d_clear(pL[i]);
  }
  poly_q_vec_d_clear(p12);
  for (size_t i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_clear(p2[i]);
  }
  return time;
}

double sample_klein_bench(timer* t) {
  double time;
  poly_q_vec_d w, z[PARAM_K];
  uint8_t seed[SEED_BYTES];

  poly_q_vec_d_init(w);
  for (size_t i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_init(z[i]);
  }
  poly_q_vec_d_uniform(w, seed, DOMAIN_SEPARATOR_U);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  sample_klein(z,w);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  poly_q_vec_d_clear(w);
  for (size_t i = 0; i < PARAM_K; i++) {
    poly_q_vec_d_clear(z[i]);
  }
  return time;
}

double truncated_sampler_bench(timer* t) {
  double time;
  sk_t sk;
  pk_t pk;
  poly_q tag;
  poly_q_mat_d_d A;
  poly_q_vec_d u, v11, v12, v2[PARAM_K_L];
  uint8_t state[STATE_BYTES];

  keys_init(&pk, &sk);
  poly_q_mat_d_d_init(A);
  poly_q_vec_d_init(u);
  poly_q_vec_d_init(v11);
  poly_q_vec_d_init(v12);
  for (size_t i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_init(v2[i]);
  }
  poly_q_init(tag);

  keygen(&pk, &sk);
  poly_q_mat_d_d_uniform(A, pk.seed, DOMAIN_SEPARATOR_A, 0);
  poly_q_vec_d_uniform(u, pk.seed, DOMAIN_SEPARATOR_U);
  randombytes(state, STATE_BYTES);
  tag_gen(tag, state);
  /* ----------- BEGIN: Code under measurement --------- */
  start_timer(t);
  truncated_sampler(v11, v12, v2, sk.R, A, pk.B, u, tag, sk.S);
  time = stop_timer(t);
  /* ----------- END: Code under measurement ----------- */
  keys_clear(&pk, &sk);
  poly_q_mat_d_d_clear(A);
  poly_q_vec_d_clear(u);
  poly_q_vec_d_clear(v11);
  poly_q_vec_d_clear(v12);
  for (size_t i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_clear(v2[i]);
  }
  poly_q_clear(tag);

  return time;
}
