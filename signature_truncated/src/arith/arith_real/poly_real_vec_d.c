#include "poly_real_vec_d.h"

/*************************************************
* Name:        poly_real_vec_d_init
*
* Description: Initialize polynomial vector with PARAM_D entries.
*              This is strictly required before any operations 
*              are done with/on the vector.
* 
* Arguments:   - poly_real_vec_d arg: polynomial vector to be initialized
**************************************************/
void poly_real_vec_d_init(poly_real_vec_d res) {
  size_t i;
  for (i = 0; i < PARAM_D; i++) {
    poly_real_init(res->entries[i]);
  }
}

/*************************************************
* Name:        poly_real_vec_d_clear
*
* Description: Clear polynomial vector with PARAM_D entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial vector must not be used again (unless reinitialized).
* 
* Arguments:   - poly_real_vec_d arg: polynomial vector to be cleared
**************************************************/
void poly_real_vec_d_clear(poly_real_vec_d res) {
  size_t i;
  for (i = 0; i < PARAM_D; i++) {
    poly_real_clear(res->entries[i]);
  }
}
