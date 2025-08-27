#ifndef POLY_Q_SAMPLING_H
#define POLY_Q_SAMPLING_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

void poly_q_mat_d_d_uniform(poly_q_mat_d_d mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, uint8_t offset);
void poly_q_mat_d_k_l_uniform(poly_q_mat_d_k_l mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);
void poly_q_vec_d_uniform(poly_q_vec_d vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator);

void poly_q_mat_d_d_binomial(poly_q_mat_d_d mat, const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator);

void poly_q_binary_fixed_weight(poly_q res, uint8_t state_in[STATE_BYTES]);

#endif /* POLY_Q_SAMPLING_H */
