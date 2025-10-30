#include "poly_q_sampling.h"
#include "fips202.h"
#include "arith.h"
#include "poly_real_vec_2d.h"
#include "random.h"

/*************************************************
* Name:        poly_q_uniform
*
* Description: Sample a uniformly random polynomial modulo
*              PARAM_Q deterministically from a seed.
* 
* Arguments:   - poly_q pout: output uniform polynomial (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - size_t i: index domain separator for XOF
*              - size_t j: index domain separator for XOF
**************************************************/
static void poly_q_uniform(poly_q pout, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, size_t i, size_t j) {
  uint8_t output[SHAKE128_RATE * 2];
  keccak_state state;
  size_t k,cnt,off,bytecnt;
  shake128_init(&state);
  shake128_absorb(&state, seed, SEED_BYTES);
  shake128_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
  shake128_absorb(&state, (const uint8_t*) &i, sizeof(i));
  shake128_absorb(&state, (const uint8_t*) &j, sizeof(j));
  shake128_finalize(&state);
  shake128_squeezeblocks(output, 2, &state);
  bytecnt = 2*SHAKE128_RATE;

  cnt = 0;
  off = 0;
  while (cnt < PARAM_N) {
#if PARAM_Q_BITLEN > 20
#error "PARAM_Q_BITLEN too big for uniform sampling."
#else
    // idea: take 5 bytes, partition into 2 chunks of 20 bits, potentially ignore the MSBs, perform rejection sampling
    if (bytecnt < 5) {
      for (k = 0; k < bytecnt; k++) {
        output[k] = output[off++];
      }
      shake128_squeezeblocks(&output[bytecnt], 1, &state);
      off = 0;
      bytecnt += SHAKE128_RATE;
    }
    uint64_t tmp5byte = output[off] | ((uint64_t)output[off+1] << 8) | ((uint64_t)output[off+2] << 16) | ((uint64_t)output[off+3] << 24) | ((uint64_t)output[off+4] << 32);
    uint64_t tmp = tmp5byte & ((1<<PARAM_Q_BITLEN)-1);
    if (tmp < (coeff_q)PARAM_Q) {
      poly_q_set_coeff(pout, cnt++, tmp);
      if (cnt == PARAM_N) {
        break;
      }
    }
    tmp = (tmp5byte >> 20) & ((1<<PARAM_Q_BITLEN)-1);
    if (tmp < (coeff_q)PARAM_Q) {
      poly_q_set_coeff(pout, cnt++, tmp);
    }
  
    off += 5;
    bytecnt -= 5;
#if PARAM_Q_BITLEN < 19
#warning "PARAM_Q_BITLEN maybe unsuitable for efficient uniform sampling."
#endif
#endif
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D x PARAM_D modulo PARAM_Q deterministically 
*              from a seed.
* 
* Arguments:   - poly_q_mat_d_d mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
*              - uint8_t offset: offset domain separator for XOF
**************************************************/
void poly_q_mat_d_d_uniform(poly_q_mat_d_d mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator, uint8_t offset) {
  size_t i,j;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_q_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j + offset);
    }
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_uniform
*
* Description: Sample a uniformly random polynomial matrix of
*              size PARAM_D x PARAM_K_L modulo PARAM_Q deterministically 
*              from a seed.
* 
* Arguments:   - poly_q_mat_d_k_l mat: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_mat_d_k_l_uniform(poly_q_mat_d_k_l mat, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i,j;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_K_L; j++) {
      poly_q_uniform(mat->rows[i]->entries[j], seed, domain_separator, i, j);
    }
  }
}

/*************************************************
* Name:        poly_q_vec_d_uniform
*
* Description: Sample a uniformly random polynomial vector of
*              size PARAM_D modulo PARAM_Q deterministically 
*              from a seed.
* 
* Arguments:   - poly_q_vec_d vec: output uniform polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_vec_d_uniform(poly_q_vec_d vec, const uint8_t seed[SEED_BYTES], uint32_t domain_separator) {
  size_t i;
  for (i = 0; i < PARAM_D; i++) {
    poly_q_uniform(vec->entries[i], seed, domain_separator, i, 0);
  }
}

/*************************************************
* Name:        poly_q_mat_d_d_binomial
*
* Description: Sample a centered binomial polynomial matrix of
*              size PARAM_D x PARAM_D with binomial parameter 1 
*              deterministically from a seed.
* 
* Arguments:   - poly_q_mat_d_d mat: output binomial polynomial matrix (initialized)
*              - const uint8_t *seed: pointer to byte array containing seed (allocated SEED_BYTES bytes)
*              - uint32_t cnt: repetition domain separator for XOF
*              - uint32_t domain_separator: domain separator for XOF
**************************************************/
void poly_q_mat_d_d_binomial(poly_q_mat_d_d mat, const uint8_t seed[SEED_BYTES], uint32_t cnt, uint32_t domain_separator) {
#if (PARAM_N%64) != 0
#error "PARAM_N must be divisible by 64"
#endif
  uint64_t output[PARAM_N*2/64]; // 2 bits per coefficient
  uint64_t coef_lsb[PARAM_N/64];
  uint64_t coef_sign[PARAM_N/64];
  keccak_state state;
  size_t i,j,k;

  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      shake256_init(&state);
      shake256_absorb(&state, seed, SEED_BYTES);
      shake256_absorb(&state, (const uint8_t*)&cnt, sizeof(uint32_t));
      shake256_absorb(&state, (const uint8_t*)&domain_separator, sizeof(uint32_t));
      shake256_absorb(&state, (const uint8_t*) &i, sizeof(i));
      shake256_absorb(&state, (const uint8_t*) &j, sizeof(j));
      shake256_finalize(&state);
      shake256_squeeze((uint8_t*)output, PARAM_N*2/8, &state);

      for (k = 0; k < PARAM_N/64; k++) {
        coef_lsb[k] = output[2*k] ^ output[2*k+1];
        coef_sign[k] = output[2*k] & output[2*k+1];
      }
      for (k = 0; k < PARAM_N; k++) {
        poly_q_set_coeff(mat->rows[i]->entries[j], k, (int32_t)((coef_lsb[k/64] >> (k%64))&1) + (int32_t)(((coef_sign[k/64] >> ((k%64))) << 1)&2) - 1);
        // we have for sign||lsb either 00 (->-1) or 01 (->0) or 10 (->1), so we reconstruct this and subtract one
      }
    }
  }
}

/*************************************************
* Name:        poly_q_binary_fixed_weight
*
* Description: Sample binary polynomial with fixed Hamming weight PARAM_W
* 
* Arguments:   - poly_q res: output uniform binary polynomial with fixed Hamming weight (initialized)
*              - uint8_t *state_in: state from which to expand the polynomial (allocated STATE_BYTES bytes)
**************************************************/
void poly_q_binary_fixed_weight(poly_q res, uint8_t state_in[STATE_BYTES]) {
  unsigned int i, b, pos = 0;
  uint8_t buf[SHAKE256_RATE];
  keccak_state state;

  shake256_absorb_once(&state, state_in, STATE_BYTES);
  shake256_squeezeblocks(buf, 1, &state);

  for (i = 0; i < PARAM_N; i++) {
    poly_q_set_coeff(res, i, 0);
  }
  for (i = PARAM_N - PARAM_W; i < PARAM_N; i++) {
    do {
      if (pos >= SHAKE256_RATE) {
        shake256_squeezeblocks(buf, 1, &state);
        pos = 0;
      }
      b = buf[pos++];
    } while (b > i);

    poly_q_set_coeff(res, i, poly_q_get_coeff(res, b));
    poly_q_set_coeff(res, b, 1);
  }
}
