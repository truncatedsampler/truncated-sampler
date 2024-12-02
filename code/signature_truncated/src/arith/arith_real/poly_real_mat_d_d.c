#include "poly_real_mat_d_d.h"

/*************************************************
* Name:        poly_real_mat_d_d_init
*
* Description: Initialize polynomial matrix with PARAM_D x PARAM_D entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_real_mat_d_d res: polynomial matrix to be initialized
**************************************************/
void poly_real_mat_d_d_init(poly_real_mat_d_d res) {
  size_t i,j;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_real_init(res->rows[i]->entries[j]);
    }
  }
}

/*************************************************
* Name:        poly_real_mat_d_d_clear
*
* Description: Clear polynomial matrix with PARAM_D x PARAM_D entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_real_mat_d_d res: polynomial matrix to be cleared
**************************************************/
void poly_real_mat_d_d_clear(poly_real_mat_d_d res) {
  size_t i,j;
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_real_clear(res->rows[i]->entries[j]);
    }
  }
}
