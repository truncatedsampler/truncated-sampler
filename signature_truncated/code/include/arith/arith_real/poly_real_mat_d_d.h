#ifndef POLY_REAL_MAT_D_D_H
#define POLY_REAL_MAT_D_D_H

#include "poly_real.h"
#include "poly_real_vec_d.h"

typedef struct {
    poly_real_vec_d rows[PARAM_D];
} __poly_real_mat_d_d;

typedef __poly_real_mat_d_d poly_real_mat_d_d[1];

void poly_real_mat_d_d_init(poly_real_mat_d_d res);
void poly_real_mat_d_d_clear(poly_real_mat_d_d res);

#endif /* POLY_REAL_MAT_D_D_H */
