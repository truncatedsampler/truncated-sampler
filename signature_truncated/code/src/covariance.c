#include "covariance.h"
#include "arith.h"

#include <flint/arb.h>
#include <flint/acb_dft.h>

/* Precision for FFT in spectral norm estimation */
#define SPECTRAL_NORM_PRECISION 128

static acb_dft_pre_t pre; // for fft precomputation

/*************************************************
* Name:        fft_precomp_setup
*
* Description: Initialize flint fft precomputation
**************************************************/
void fft_precomp_setup(void) {
  acb_dft_precomp_init(pre, 2*PARAM_N, SPECTRAL_NORM_PRECISION);
}

/*************************************************
* Name:        fft_precomp_teardown
*
* Description: Clear flint fft precomputation
**************************************************/
void fft_precomp_teardown(void) {
  acb_dft_precomp_clear(pre);
}

/*************************************************
* Name:        compute_covariance
*
* Description: Computes Schur complements of covariance for perturbation 
*              sampling from the signer's secret key
*
* Arguments:   - poly_real_mat_2d_2d S: polynomial matrix to host the sampling materials
*              - const poly_q_mat_d_d *RRstar: array of polynomial matrices containing R.R*
**************************************************/
void compute_covariance(poly_real_mat_2d_2d S, const poly_q_mat_d_d RRstar[2][2]) {
  size_t i,j,k;
  poly_real finv, tmppoly;

  // init polynomials
  poly_real_init(finv);
  poly_real_init(tmppoly);

  // convert RRstar to poly_real_mat_2d_2d and write in S and multiply by -1/(1/sG^2 - 1/s4^2)) (RRstar needs to be read in centered representation!)
  for (i = 0; i < 2*PARAM_D; i++) {
    for (j = 0; j < 2*PARAM_D; j++) {
      poly_real_from_poly_q(S->rows[i]->entries[j], RRstar[i/PARAM_D][j/PARAM_D]->rows[i%PARAM_D]->entries[j%PARAM_D]);
      poly_real_mul_scalar(S->rows[i]->entries[j], S->rows[i]->entries[j], PARAM_SGINVSQ_S4INVSQ);
    }
  }

  // Add diag(s1^2, s3^2)
  for (i = 0; i < PARAM_D; i++) {
    poly_real_add_constant(S->rows[i]->entries[i], PARAM_S1SQ);
  }
  for (i = PARAM_D; i < 2*PARAM_D; i++) {
    poly_real_add_constant(S->rows[i]->entries[i], PARAM_S3SQ);
  }

  for (i = 2*PARAM_D-1; i > PARAM_D-1; i--) {
    // i-th row and i-th column as well as everything beyond stands as it is

    // invert f_i
    poly_real_invert(finv, S->rows[i]->entries[i]);

    // multiply finv to si (i-th column until excluding i-th polynomial, which is f_i)
    for (j = 0; j < i; j++) {
      poly_real_mul(S->rows[j]->entries[i], S->rows[j]->entries[i], finv);
    }

    // update S->rows[j]->entries[k] with j=0..i-1, k=0..i-1 by subtracting S->rows[j]->entries[i] times S->rows[i]->entries[k]
    // note that S->rows[j]->entries[i] has finv already multiplied
    for (j = 0; j < i; j++) {
      for (k = 0; k < i; k++) {
        poly_real_mul(tmppoly, S->rows[j]->entries[i], S->rows[i]->entries[k]);
        poly_real_sub(S->rows[j]->entries[k], S->rows[j]->entries[k], tmppoly);
      }
    }
  }
  
  // clean up polynomials
  poly_real_clear(finv);
  poly_real_clear(tmppoly);
}

/*************************************************
* Name:        sk_sq_spectral_norm
*
* Description: Computes an estimate to the square spectral norm of R
*              based on the iterated power method and computes R.R*
*
* Arguments:   - poly_q_mat_d_d *RRstar: array of polynomial matrices hosting R.R*
*              - const poly_q_mat_d_d *R: polynomial matrix containing the secret key R
* 
* Returns a double corresponding to the estimated spectral norm of R
**************************************************/
int sk_sq_spectral_norm(poly_q_mat_d_d RRstar[2][2], const poly_q_mat_d_d R[2][PARAM_K_L]) {
  size_t i,j,k;
  acb_struct num_vec[PARAM_N*2], den_vec[PARAM_N*2];
  poly_q_mat_d_d Rstar[PARAM_K_L][2];
  poly_z_mat_d_d RRstar_z, power_mat, tmp_mat;
  poly_q_mat_d_d tmp_q;
  fmpz_t c;
  arb_t num_embedding, den_embedding;
  double sq_norm_R, test_norm;
  int is_invalid;

  // init flint objects for fft
  for (i = 0; i < PARAM_N*2; i++) {
    acb_init(&num_vec[i]);
    acb_init(&den_vec[i]);
  }
  arb_init(num_embedding);
  arb_init(den_embedding);
  fmpz_init(c);
  
  // init matrices
  poly_z_mat_d_d_init(RRstar_z);
  poly_z_mat_d_d_init(power_mat);
  poly_z_mat_d_d_init(tmp_mat);
  poly_q_mat_d_d_init(tmp_q);
  // init and compute R*
  for (i = 0; i < 2; i++) {
    for (j = 0; j < PARAM_K_L; j++) {
      poly_q_mat_d_d_init(Rstar[j][i]);
      poly_q_mat_d_d_conjugate(Rstar[j][i], R[i][j]);
    }
  }

  // multiply R by R*
  for (i = 0; i < 2; i++) {
    for (j = 0; j < 2; j++) {
      poly_q_mat_d_d_mul_mat_d_d(RRstar[i][j], R[i][0], Rstar[0][j]);
      for (k = 1; k < PARAM_K_L; k++) {
        poly_q_mat_d_d_mul_mat_d_d(tmp_q, R[i][k], Rstar[k][j]);
        poly_q_mat_d_d_add(RRstar[i][j], RRstar[i][j], tmp_q);
      }
    }
  }

  is_invalid = 0;
  for (i = 0; i < 2; i++) {
    // convert RRstar[i][i] to poly_z_mat_d_d and write in RRstar_z[i] (needs to be read in centered representation!)
    poly_z_mat_d_d_from_poly_q_mat_d_d(RRstar_z, RRstar[i][i]);
    poly_z_mat_d_d_set(power_mat, RRstar_z);
    // Iterated power: compute (Ri.Ri*)^32 and store in power_mat and tmp_mat
    for (k = 0; k < PARAM_IT_SPEC_NORM; k++) {
      poly_z_mat_d_d_mul_mat_d_d(tmp_mat, power_mat, power_mat);
      poly_z_mat_d_d_set(power_mat, tmp_mat);
    }

    // multiplying power_mat by Ri.Ri*. power_mat contains (Ri.Ri*)^33 and tmp_mat contains (Ri.Ri*)^32
    poly_z_mat_d_d_mul_mat_d_d(power_mat, tmp_mat, RRstar_z);
    
    // fft of (R1.R1^*)^32[0][0] and (R1.R1^*)^33[0][0]
    for (j = 0; j < PARAM_N; j++) {
      poly_z_get_coeff(c, power_mat->rows[0]->entries[0], j);
      acb_set_fmpz(&num_vec[j], c);
      poly_z_get_coeff(c, tmp_mat->rows[0]->entries[0], j);
      acb_set_fmpz(&den_vec[j], c);
    }
    acb_dft_precomp(num_vec, num_vec, pre, SPECTRAL_NORM_PRECISION);
    acb_dft_precomp(den_vec, den_vec, pre, SPECTRAL_NORM_PRECISION);

    // finding maximal embedding ratio (spectral norm estimate)
    sq_norm_R = 0.0;
    for (j = 0; j < PARAM_N; j++) {
      acb_get_real(num_embedding, &num_vec[2*j+1]);
      acb_get_real(den_embedding, &den_vec[2*j+1]);
      arb_div(num_embedding, num_embedding, den_embedding, SPECTRAL_NORM_PRECISION);
      test_norm = arf_get_d(arb_midref(num_embedding), ARF_RND_DOWN); // round to arb midpoint

      // maximum test
      if (test_norm > sq_norm_R) {
        sq_norm_R = test_norm;
      }
    }
    if (sq_norm_R > PARAM_R_MAX_SQ_SPECTRAL_NORM) {
      is_invalid = 1;
      goto spectral_cleanup;
    }
  }

spectral_cleanup:
  // cleanup flint objects for fft
  for (i = 0; i < 2*PARAM_N; i++) {
    acb_clear(&num_vec[i]);
    acb_clear(&den_vec[i]);
  }
  arb_clear(num_embedding);
  arb_clear(den_embedding);
  fmpz_clear(c);

  // cleanup matrices
  poly_z_mat_d_d_clear(RRstar_z);
  poly_z_mat_d_d_clear(power_mat);
  poly_z_mat_d_d_clear(tmp_mat);
  poly_q_mat_d_d_clear(tmp_q);
  // init and compute R*
  for (i = 0; i < 2; i++) {
    for (j = 0; j < PARAM_K_L; j++) {
      poly_q_mat_d_d_clear(Rstar[j][i]);
    }
  }

  return is_invalid;
}
