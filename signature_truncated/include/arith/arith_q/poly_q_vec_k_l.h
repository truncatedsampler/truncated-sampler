#ifndef POLY_Q_VEC_K_L_H
#define POLY_Q_VEC_K_L_H

#include "poly_q.h"

typedef struct {
	poly_q entries[PARAM_K_L];
} __poly_q_vec_k_l;

typedef __poly_q_vec_k_l poly_q_vec_k_l[1];

void poly_q_vec_k_l_setup(void);
void poly_q_vec_k_l_teardown(void);

void poly_q_vec_k_l_init(poly_q_vec_k_l arg);
void poly_q_vec_k_l_clear(poly_q_vec_k_l arg);
void poly_q_vec_k_l_zero(poly_q_vec_k_l arg);
void poly_q_vec_k_l_set(poly_q_vec_k_l res, const poly_q_vec_k_l arg);
void poly_q_vec_k_l_get_poly(poly_q res, const poly_q_vec_k_l arg, size_t pos);
void poly_q_vec_k_l_set_poly(poly_q_vec_k_l res, const poly_q arg, size_t pos);
void poly_q_vec_k_l_neg(poly_q_vec_k_l res, const poly_q_vec_k_l arg);
void poly_q_vec_k_l_add(poly_q_vec_k_l res, const poly_q_vec_k_l lhs, const poly_q_vec_k_l rhs);
void poly_q_vec_k_l_sub(poly_q_vec_k_l res, const poly_q_vec_k_l lhs, const poly_q_vec_k_l rhs);
void poly_q_vec_k_l_mul_scalar(poly_q_vec_k_l res, const poly_q_vec_k_l arg, coeff_q fac);
void poly_q_vec_k_l_mul_poly(poly_q_vec_k_l res, const poly_q_vec_k_l arg, const poly_q fac);
void poly_q_vec_k_l_mul_inner(poly_q res, const poly_q_vec_k_l lhs, const poly_q_vec_k_l rhs);
uint64_t poly_q_vec_k_l_norm2(const poly_q_vec_k_l arg);
void poly_q_vec_k_l_sample_gaussian_s4(poly_q_vec_k_l res);
int poly_q_vec_k_l_equal(const poly_q_vec_k_l lhs, const poly_q_vec_k_l rhs);
void poly_q_vec_k_l_dump(const poly_q_vec_k_l arg);

#endif /* POLY_Q_VEC_K_L_H */
