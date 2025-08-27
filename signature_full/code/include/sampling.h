#ifndef SAMPLING_H
#define SAMPLING_H

#include <stdint.h>
#include "params.h"
#include "arith.h"

void sample_perturb(
    poly_q_vec_d p1[2], 
    poly_q_vec_d p2[PARAM_K], 
    const poly_q_mat_d_d R[2][PARAM_K], 
    const poly_real_mat_2d_2d S);
void sample_klein(poly_q_vec_d v[PARAM_K], const poly_q_vec_d w);
void elliptic_sampler(
    poly_q_vec_d v11, 
    poly_q_vec_d v12, 
    poly_q_vec_d v2[PARAM_K], 
    const poly_q_mat_d_d R[2][PARAM_K], 
    const poly_q_mat_d_d A, 
    const poly_q_mat_d_d B[PARAM_K], 
    const poly_q_vec_d u, 
    const poly_q tag, 
    const poly_real_mat_2d_2d S);

#endif /* SAMPLING_H */
