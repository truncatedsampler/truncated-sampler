#include "sampling.h"
#include "fips202.h"
#include "arith.h"
#include "poly_real_vec_2d.h"
#include "random.h"

/*************************************************
* Name:        sample_perturb
*
* Description: Gaussian perturbation sampler (SamplePerturb, Algorithm 4.3)
* 
* Arguments:   - poly_q_vec_d *p1: array of polynomial vectors to host top perturbation befor K-projection (initialized)
*              - poly_q_vec_d p12: polynomial vector to host middle perturbation (initialized)
*              - poly_q_vec_d *p2: array of polynomial vectors to host bottom perturbation (initialized)
*              - const poly_q m_sG_t_tstar: polynomial corresponding to -s_G^2.t.t*
*              - const poly_q_mat_d_d *R: array of polynomial matrices, secret key R 
*              - const poly_real_mat_2d_2d *S: array of polynomial matrices, Schur complements for sampling
**************************************************/
void sample_perturb(
    poly_q_vec_d pL[PARAM_L], 
    poly_q_vec_d p12, 
    poly_q_vec_d p2[PARAM_K_L], 
    const poly_real m_sG_t_tstar, 
    const poly_q_mat_d_d R[2][PARAM_K_L], 
    const poly_real_mat_2d_2d S) {
  size_t i,j,k;
  poly_real tmp_sub, tmp_mul;
  poly_real_vec_2d c;
  poly_real_mat_d_d S_d;

  // init vectors and polynomials
  poly_real_init(tmp_sub);
  poly_real_init(tmp_mul);
  poly_real_vec_2d_init(c);
  poly_real_mat_d_d_init(S_d);

  // sample p_{L,2} from Gaussian with variance s_2^2 + m_sG_t_tstar = s_2^2 - s_G^2.t.t*
  poly_real_copy_add_constant(tmp_sub, m_sG_t_tstar, (coeff_real)PARAM_S2SQ);
  poly_real_zero(tmp_mul);
  for (i = 1; i < PARAM_L; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_q_samplefz(pL[i]->entries[j], tmp_sub, tmp_mul);
    }
  }

  // sample p2 from Gaussian distribution with variance sqrt(s_4^2 - s_G^2)
  for (i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_gaussian_sqrt_s4sq_sGsq(p2[i]); // probabilistic
  }

  // pL[0] = R1.p2, p12 = R2.p2  (pL[0] and p12 only used as intermediate variables)
  poly_q_mat_d_d_mul_vec_d(pL[0], R[0][0], p2[0]);
  poly_q_mat_d_d_mul_vec_d(p12, R[1][0], p2[0]);
  for (i = 1; i < PARAM_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(pL[0], R[0][i], p2[i]);
    poly_q_mat_d_d_muladd_vec_d(p12, R[1][i], p2[i]);
  }

  // convert p1 to c, which is a poly_real_vec_2d, and scale with -s_G^2/(s_4^2 - s_G^2)
  for (i = 0; i < PARAM_D; i++) {
    poly_real_from_poly_q(c->entries[i          ], pL[0]->entries[i]);
    poly_real_from_poly_q(c->entries[i + PARAM_D], p12->entries[i]);
    poly_real_mul_scalar(c->entries[i], c->entries[i], PARAM_NEGSGSQ_DIV_S4SQ_SGSQ);
    poly_real_mul_scalar(c->entries[i + PARAM_D], c->entries[i + PARAM_D], PARAM_NEGSGSQ_DIV_S4SQ_SGSQ);
  }

  // sample p12 with pre-computed f_i
  for (i = 2*PARAM_D; i > PARAM_D; i--) {
    // sample p12->entries[i-1-PARAM_D] from a discrete Gaussian distribution with
    // center c->entries[i-1] and variance M_tau(S->rows[i-1]->entries[i-1])
    poly_q_samplefz(p12->entries[i - 1 - PARAM_D], S->rows[i-1]->entries[i-1], c->entries[i-1]); // probabilistic

    // update the entries from c with indices j=0..i-2 by adding S->rows[j]->entries[i-1] (this already contains f_i^-1)
    poly_real_sub_poly_real_poly_q(tmp_sub, p12->entries[i - 1 - PARAM_D], c->entries[i-1]);
    for (j = 0; j < i-1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S->rows[j]->entries[i-1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  // ONLINE PART: sample p_{L,1} with f_i depending on tag
  // setting S_d to F(R1,R2) - s_G^2.t.t*.I_d
  for (i = 0; i < PARAM_D; i++) {
    for (j = 0; j < PARAM_D; j++) {
      poly_real_set(S_d->rows[i]->entries[j], S->rows[i]->entries[j]);
    }
    poly_real_add(S_d->rows[i]->entries[i], S_d->rows[i]->entries[i], m_sG_t_tstar);
  }

  for (i = PARAM_D; i > 0; i--) {
    // sample pL[0]->entries[i-1] from a discrete Gaussian distribution with
    // center c->entries[i-1] and variance M_tau(S->rows[i-1]->entries[i-1])
    poly_q_samplefz(pL[0]->entries[i-1], S_d->rows[i-1]->entries[i-1], c->entries[i-1]); // probabilistic

    // update S_d entries for j=0..i-2, k=0..i-2
    // invert f_i
    poly_real_invert(tmp_sub, S_d->rows[i-1]->entries[i-1]);

    // multiply f_i^{-1} to si (i-th column until excluding i-th polynomial, which is f_i)
    for (j = 0; j < i-1; j++) {
      poly_real_mul(S_d->rows[j]->entries[i-1], S_d->rows[j]->entries[i-1], tmp_sub);
    }

    // update S->rows[j]->entries[k] with j=0..i-2, k=0..i-2 by subtracting S->rows[j]->entries[i-1] times S->rows[i]->entries[k]
    // note that S->rows[j]->entries[i] has finv already multiplied
    for (j = 0; j < i-1; j++) {
      for (k = 0; k < i-1; k++) {
        poly_real_mul(tmp_mul, S_d->rows[j]->entries[i-1], S_d->rows[i-1]->entries[k]);
        poly_real_sub(S_d->rows[j]->entries[k], S_d->rows[j]->entries[k], tmp_mul);
      }
    }

    // update the entries from c with indices j=0..i-2 by adding S_d->rows[j]->entries[i-1] (this already contains f_i^-1)
    poly_real_sub_poly_real_poly_q(tmp_sub, pL[0]->entries[i-1], c->entries[i-1]);
    for (j = 0; j < i-1; j++) {
      poly_real_mul(tmp_mul, tmp_sub, S_d->rows[j]->entries[i-1]);
      poly_real_add(c->entries[j], c->entries[j], tmp_mul);
    }
  }

  //  clean up vectors and polynomials
  poly_real_mat_d_d_clear(S_d);
  poly_real_vec_2d_clear(c);
  poly_real_clear(tmp_mul);
  poly_real_clear(tmp_sub);
}

#if PARAM_K == 5
// Gaussian widths for Klein Sampler
static const double klein_widths[PARAM_K] = {
  3.42997312037747725810L,
  3.44004613446786633446L,
  3.44010564823841225035L,
  3.44010600038776326315L,
  3.45610035816311800261L
};
// Scaled Gram-Schmidt of the basis of the gadget lattice (must be negated)
static const double neg_scaled_gadget_gso[PARAM_K][PARAM_K] = {
  {-0.07647058823529412352L,  -0.00588214820229020913L, -0.00045247284545076348L, -0.00003480560345379704L, -0.00000269779563119137L},
  {0.00588235294117647051L, -0.07646792662977272559L, -0.00588214699085992491L, -0.00045247284489936151L, -0.00003507134320546393L},
  {0.00000000000000000000L, 0.00591695381295464801L,  -0.07646791088117903257L, -0.00588214698369169970L, -0.00045592746167107091L},
  {0.00000000000000000000L, 0.00000000000000000000L,  0.00591715854467268331L,  -0.07646791078799210217L, -0.00592705700172388959L},
  {0.00000000000000000000L, 0.00000000000000000000L,  0.00000000000000000000L,  0.00591715975610271599L,  -0.07705174102241059420L}
};
// Gadget lattice basis
static const int64_t gadget_basis[PARAM_K][PARAM_K] = {
  {13,   0,  0,  0,  4},
  {-1,  13,  0,  0,  4},
  { 0,  -1, 13,  0,  9},
  { 0,   0, -1, 13, 12},
  { 0,   0,  0, -1, 12}
};
#else
#error "Some constants are not generated for the given PARAM_K."
#endif

/*************************************************
* Name:        sample_klein
*
* Description: Gaussian Klein sampler on the gadget lattice for SamplePre
* 
* Arguments:   - poly_q_vec_d *v: array of polynomial vectors to host sample (initialized)
*              - const poly_q_vec_d w: polynomial vector, center
**************************************************/
void sample_klein(poly_q_vec_d v[PARAM_K], const poly_q_vec_d w) {
  int64_t i,j,i_1,i_2,i_2_1,i_2_2;
  double di;
  int64_t zi;
  coeff_q w_coeff, v_coeff;

  for (i = PARAM_N * PARAM_D * PARAM_K; i > 0; i--) {
    i_1   = (i - 1) / (PARAM_N * PARAM_D);
    i_2   = (i - 1) % (PARAM_N * PARAM_D);
    i_2_1 = i_2 / PARAM_N;
    i_2_2 = i_2 % PARAM_N;

    // di = <v - [-w|0|0...], Btilde'[i]>
    w_coeff = poly_q_get_coeff(w->entries[i_2_1], (size_t) i_2_2);
    // coeffs of w do not have to be in centered representation
    if (i == PARAM_N * PARAM_D * PARAM_K) {
      di = neg_scaled_gadget_gso[0][i_1] * ((double) w_coeff);
      // "below" w, everything is zero and does not contribute to the dot product
    } else {
      // coeffs of v have to be in centered representation to compute di
      v_coeff = poly_q_get_coeff_centered(v[0]->entries[i_2_1], (size_t) i_2_2);
      di = neg_scaled_gadget_gso[0][i_1] * ((double) (w_coeff + v_coeff));
      for (j = 1; j < PARAM_K; j++) {
        v_coeff = poly_q_get_coeff_centered(v[j]->entries[i_2_1], (size_t) i_2_2);
        di += neg_scaled_gadget_gso[j][i_1] * ((double) v_coeff);
      }
    }

    // zi following D_{Z, si, di}
    zi = SampleZ(di, klein_widths[i_1]);

    // v = v + zi*B[i]
    for (j = 0; j < PARAM_K; j++) {
      v_coeff = poly_q_get_coeff(v[j]->entries[i_2_1], (size_t) i_2_2);
      v_coeff += zi * gadget_basis[j][i_1];
      poly_q_set_coeff(v[j]->entries[i_2_1], (size_t) i_2_2, v_coeff);
    }
  }
}

/*************************************************
* Name:        truncated_sampler
*
* Description: Elliptic Gaussian preimage sampler (SamplePre)
* 
* Arguments:   - poly_q_vec_d *v1: array of polynomial vectors to host top preimage (initialized)
*              - poly_q_vec_d *v2: array of polynomial vectors to host bottom preimage (initialized)
*              - const poly_q_mat_d_d *R: array of polynomial matrices, secret key R 
*              - const poly_q_mat_d_d A: polynomial matrix, public A'
*              - const poly_q_mat_d_d *B: array of polynomial matrices, public B
*              - const poly_q_vec_d u: polynomial vector, public u
*              - const poly_q tag: polynomial, tag
*              - const poly_real_mat_2d_2d *S: array of polynomial matrices, Schur complements for sampling
**************************************************/
void truncated_sampler(
    poly_q_vec_d v11, 
    poly_q_vec_d v12, 
    poly_q_vec_d v2[PARAM_K_L], 
    const poly_q_mat_d_d R[2][PARAM_K_L], 
    const poly_q_mat_d_d A, 
    const poly_q_mat_d_d B[PARAM_K_L], 
    const poly_q_vec_d u, 
    const poly_q tag, 
    const poly_real_mat_2d_2d S) {
  size_t i;
  poly_q_vec_d pL[PARAM_L];
  poly_q_vec_d p12;
  poly_q_vec_d p2[PARAM_K_L];
  poly_q taginv;
  poly_q_vec_d w, tmp, y[PARAM_K];
  poly_real m_sG_t_tstar;
  int64_t bexpi;

  // init vectors and polynomials
  poly_q_init(taginv);
  poly_real_init(m_sG_t_tstar);
  poly_q_vec_d_init(w);
  poly_q_vec_d_init(tmp);
   for (i = 0; i < PARAM_L; i++) {
    poly_q_vec_d_init(pL[i]);
    poly_q_vec_d_init(y[i]);
  }
  poly_q_vec_d_init(p12);
  for (i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_init(p2[i]);
    poly_q_vec_d_init(y[PARAM_L + i]);
  }

  // compute -s_G^2.t.t* (taginv used as temp variable)
  poly_q_conjugate(taginv, tag);
  poly_q_mul(taginv, taginv, tag);
  poly_real_from_poly_q(m_sG_t_tstar, taginv);
  poly_real_mul_scalar(m_sG_t_tstar, m_sG_t_tstar, (coeff_real)PARAM_NEG_SGSQ);

  // sample_perturb
  sample_perturb(pL, p12, p2, m_sG_t_tstar, R, S);
  
  // w = tag^{-1} * (u - G_L.pL - A.p12 + B.p2) - G_H*p_2
  // compute G_L.pL into pL[0] (stored for later) 
  bexpi = PARAM_B;
  for (i = 1; i < PARAM_L; i++) {
    poly_q_vec_d_mul_scalar(tmp, pL[i], bexpi);
    poly_q_vec_d_add(pL[0], pL[0], tmp);
    bexpi *= PARAM_B;
  }
  // compute w
  poly_q_invert(taginv, tag);
  poly_q_vec_d_sub(w, u, pL[0]);
  poly_q_mat_d_d_mulsub_vec_d(w, A, p12);
  for (i = 0; i < PARAM_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(w, B[i], p2[i]);
  }
  poly_q_vec_d_mul_poly(w, w, taginv);
  for (i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_mul_scalar(tmp, p2[i], bexpi);
    poly_q_vec_d_sub(w, w, tmp);
    bexpi *= PARAM_B;
  }

  // c = G^-1(w) = [w|0|0|...|0], but we only store and pass w
  // sample y from Klein sampler
  sample_klein(y, w);

  // z = c + y; z is represented by y now
  poly_q_vec_d_add(y[0], y[0], w);

  // v11 = G_L.v_L = G_L.pL + t.G_L.z_L + R1.zH
  poly_q_vec_d_mul_poly(y[0], y[0], tag);
  poly_q_vec_d_add(pL[0], pL[0], y[0]);
  bexpi = PARAM_B;
  for (i = 1; i < PARAM_L; i++) {
    poly_q_vec_d_mul_poly(y[i], y[i], tag);
    poly_q_vec_d_mul_scalar(y[i], y[i], bexpi);
    poly_q_vec_d_add(pL[0], pL[0], y[i]);
  }
  poly_q_mat_d_d_mul_vec_d(v11, R[0][0], y[PARAM_L + 0]);
  for (i = 1; i < PARAM_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(v11, R[0][i], y[PARAM_L + i]);
  }
  poly_q_vec_d_add(v11, v11, pL[0]); // pL[0] hosts G_L.pL + t.G_L.z_L

  // v12 = p12 + R2.z_H
  poly_q_mat_d_d_mul_vec_d(v12, R[1][0], y[PARAM_L + 0]);
  for (i = 1; i < PARAM_K_L; i++) {
    poly_q_mat_d_d_muladd_vec_d(v12, R[1][i], y[PARAM_L + i]);
  }
  poly_q_vec_d_add(v12, v12, p12);

  // v2 = p2 + z
  for (i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_add(v2[i], p2[i], y[PARAM_L + i]);
  }

  // clean up vectors and polynomials
  poly_q_clear(taginv);
  poly_real_clear(m_sG_t_tstar);
  poly_q_vec_d_clear(w);
  poly_q_vec_d_clear(tmp);
   for (i = 0; i < PARAM_L; i++) {
    poly_q_vec_d_clear(pL[i]);
    poly_q_vec_d_clear(y[i]);
  }
  poly_q_vec_d_clear(p12);
  for (i = 0; i < PARAM_K_L; i++) {
    poly_q_vec_d_clear(p2[i]);
    poly_q_vec_d_clear(y[PARAM_L + i]);
  }
}