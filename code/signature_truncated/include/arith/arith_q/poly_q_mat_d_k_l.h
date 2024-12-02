#ifndef POLY_Q_MAT_D_K_L_H
#define POLY_Q_MAT_D_K_L_H

#include "poly_q.h"
#include "poly_q_vec_k_l.h"
#include "poly_q_vec_d.h"

typedef struct {
    poly_q_vec_k_l rows[PARAM_D];
} __poly_q_mat_d_k_l;

typedef __poly_q_mat_d_k_l poly_q_mat_d_k_l[1];

void poly_q_mat_d_k_l_setup(void);
void poly_q_mat_d_k_l_teardown(void);

void poly_q_mat_d_k_l_init(poly_q_mat_d_k_l res);
void poly_q_mat_d_k_l_clear(poly_q_mat_d_k_l res);
void poly_q_mat_d_k_l_zero(poly_q_mat_d_k_l res);
void poly_q_mat_d_k_l_set(poly_q_mat_d_k_l res, const poly_q_mat_d_k_l arg);
void poly_q_mat_d_k_l_neg(poly_q_mat_d_k_l res, const poly_q_mat_d_k_l arg);
void poly_q_mat_d_k_l_add(poly_q_mat_d_k_l res, const poly_q_mat_d_k_l lhs, const poly_q_mat_d_k_l rhs);
void poly_q_mat_d_k_l_sub(poly_q_mat_d_k_l res, const poly_q_mat_d_k_l lhs, const poly_q_mat_d_k_l rhs);
void poly_q_mat_d_k_l_mul_vec_k_l(poly_q_vec_d res, const poly_q_mat_d_k_l lhs, const poly_q_vec_k_l rhs);
int poly_q_mat_d_k_l_equal(const poly_q_mat_d_k_l lhs, const poly_q_mat_d_k_l rhs);
void poly_q_mat_d_k_l_dump(const poly_q_mat_d_k_l arg);

#endif /* POLY_Q_MAT_D_K_L_H */
