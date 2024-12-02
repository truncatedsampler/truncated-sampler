#include "arith_q.h"
#include "macros.h"

static poly_q TMP;

/*************************************************
* Name:        poly_q_mat_d_k_l_setup
*
* Description: Initialize and setup the backend for arithmetic 
*              mod PARAM_Q integer matrices with PARAM_D x PARAM_K_L entries. 
*              This is strictly required and must be called once 
*              before any other function from here is used.
**************************************************/
void poly_q_mat_d_k_l_setup(void) {
  poly_q_init(TMP);
}

/*************************************************
* Name:        poly_q_mat_d_k_l_teardown
*
* Description: Clean up and teardown the backend for arithmetic 
*              mod PARAM_Q integer matrices with PARAM_D x PARAM_K_L entries. 
*              This is strictly required and must be called once 
*              at the very end to release any resources.
**************************************************/
void poly_q_mat_d_k_l_teardown(void) {
  poly_q_clear(TMP);
}

/*************************************************
* Name:        poly_q_mat_d_k_l_init
*
* Description: Initialize polynomial matrix with PARAM_D x PARAM_K_L entries.
*              This is strictly required before any operations 
*              are done with/on the matrix.
* 
* Arguments:   - poly_q_mat_d_k_l res: polynomial matrix to be initialized
**************************************************/
void poly_q_mat_d_k_l_init(poly_q_mat_d_k_l res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_l_init(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_clear
*
* Description: Clear polynomial matrix with PARAM_D x PARAM_K_L entries.
*              This is strictly required to avoid memory leaks and the 
*              polynomial matrix must not be used again (unless reinitialized).
* 
* Arguments:   - poly_q_mat_d_k_l res: polynomial matrix to be cleared
**************************************************/
void poly_q_mat_d_k_l_clear(poly_q_mat_d_k_l res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_l_clear(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_zero
*
* Description: Set an initialized polynomial matrix with PARAM_D x PARAM_K_L entries to zero
* 
* Arguments:   - poly_q_mat_d_k_l res: polynomial matrix to be zeroized (initialized)
**************************************************/
void poly_q_mat_d_k_l_zero(poly_q_mat_d_k_l res) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_l_zero(res->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_set
*
* Description: Set a polynomial matrix with PARAM_D x PARAM_K_L entries equal to another polynomial matrix
* 
* Arguments:   - poly_q_mat_d_k_l res: polynomial matrix to be set (initialized)
*              - const poly_q_mat_d_k_l arg: polynomial matrix to be read
**************************************************/
void poly_q_mat_d_k_l_set(poly_q_mat_d_k_l res, const poly_q_mat_d_k_l arg) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_l_set(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_neg
*
* Description: Negate a polynomial matrix with PARAM_D x PARAM_K_L entries
* 
* Arguments:   - poly_q_mat_d_k_l res: polynomial matrix to host the negation (initialized)
*              - const poly_q_mat_d_k_l arg: polynomial matrix to be negated
**************************************************/
void poly_q_mat_d_k_l_neg(poly_q_mat_d_k_l res, const poly_q_mat_d_k_l arg) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_l_neg(res->rows[i], arg->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_add
*
* Description: Add two polynomial matrices with PARAM_D x PARAM_K_L entries
* 
* Arguments:   - poly_q_mat_d_k_l res: polynomial matrix to host the sum (initialized)
*              - const poly_q_mat_d_k_l lhs: first polynomial matrix summand
*              - const poly_q_mat_d_k_l rhs: second polynomial matrix summand
**************************************************/
void poly_q_mat_d_k_l_add(poly_q_mat_d_k_l res, const poly_q_mat_d_k_l lhs, const poly_q_mat_d_k_l rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_l_add(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_sub
*
* Description: Substract two polynomial matrices with PARAM_D x PARAM_K_L entries
* 
* Arguments:   - poly_q_mat_d_k_l res: polynomial matrix to host the difference (initialized)
*              - const poly_q_mat_d_k_l lhs: first polynomial matrix term
*              - const poly_q_mat_d_k_l rhs: second polynomial matrix term
**************************************************/
void poly_q_mat_d_k_l_sub(poly_q_mat_d_k_l res, const poly_q_mat_d_k_l lhs, const poly_q_mat_d_k_l rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_l_sub(res->rows[i], lhs->rows[i], rhs->rows[i]);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_mul_vec_k_l
*
* Description: Product of a polynomial matrix with PARAM_D x PARAM_K_L entries
*              with a polynomial vector with PARAM_K_L entries
* 
* Arguments:   - poly_q_vec_d res: polynomial vector to host the multiplication (initialized)
*              - const poly_q_mat_d_k_l lhs: polynomial matrix to multiply
*              - const poly_q_vec_k_l rhs: polynomial vector to multiply
**************************************************/
void poly_q_mat_d_k_l_mul_vec_k_l(poly_q_vec_d res, const poly_q_mat_d_k_l lhs, const poly_q_vec_k_l rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    poly_q_vec_k_l_mul_inner(TMP, lhs->rows[i], rhs);
    poly_q_vec_d_set_poly(res, TMP, i);
  }
}

/*************************************************
* Name:        poly_q_mat_d_k_l_equal
*
* Description: Equality test between two polynomial matrices with PARAM_D x PARAM_K_L entries
* 
* Arguments:   - const poly_q_mat_d_k_l lhs: first polynomial matrix
*              - const poly_q_mat_d_k_l rhs: second polynomial matrix
* 
* Returns 1 if the polynomial matrices are equal, 0 otherwise
**************************************************/
int poly_q_mat_d_k_l_equal(const poly_q_mat_d_k_l lhs, const poly_q_mat_d_k_l rhs) {
  for (size_t i = 0; i < PARAM_D; ++i) {
    if (!poly_q_vec_k_l_equal(lhs->rows[i], rhs->rows[i])) {
      return 0;
    }
  }
  return 1;
}

/*************************************************
* Name:        poly_q_mat_d_k_l_dump
*
* Description: Print a polynomial matrix with PARAM_D x PARAM_K_L entries
* 
* Arguments:   - const poly_q_mat_d_k_l arg: polynomial matrix to be printed
**************************************************/
void poly_q_mat_d_k_l_dump(const poly_q_mat_d_k_l arg) {
	printf("[");
	for (size_t i = 0; i < PARAM_D - 1; ++i) {
		poly_q_vec_k_l_dump(arg->rows[i]);
		printf(", ");
	}
	poly_q_vec_k_l_dump(arg->rows[PARAM_D - 1]);
	printf("]");
}
