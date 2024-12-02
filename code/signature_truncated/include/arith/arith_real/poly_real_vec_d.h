#ifndef POLY_REAL_VEC_D_H
#define POLY_REAL_VEC_D_H

#include "poly_real.h"

typedef struct {
    poly_real entries[PARAM_D];
} __poly_real_vec_d;

typedef __poly_real_vec_d poly_real_vec_d[1];

void poly_real_vec_d_init(poly_real_vec_d res);
void poly_real_vec_d_clear(poly_real_vec_d res);

#endif /* POLY_REAL_VEC_D_H */
