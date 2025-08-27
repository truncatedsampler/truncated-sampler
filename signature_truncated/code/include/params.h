#ifndef PARAMS_H
#define PARAMS_H

/*************************************************
* Domain separators for XOF expansion
**************************************************/
// sep
#define DOMAIN_SEPARATOR_A 0
#define DOMAIN_SEPARATOR_R 1
#define DOMAIN_SEPARATOR_A3 2
#define DOMAIN_SEPARATOR_U 3
#define DOMAIN_SEPARATOR_D 4

/*************************************************
* Signature parameters
**************************************************/
// Ring degree for the signature
#define PARAM_N 256
// Modulus for the signature
#define PARAM_Q 370673L
// Modulus bit-length upper bound for uniform sampling
#define PARAM_Q_BITLEN 19
// Module rank for the signature
#define PARAM_D 4
// Gadget dimension
#define PARAM_K 5
// Number of truncated gadget dimensions
#define PARAM_L 2
// Truncated gadget dimension
#define PARAM_K_L PARAM_K - PARAM_L
// Gadget base
#define PARAM_B 13
// Upper gadget base b^l
#define PARAM_B_POW_L 169
// Number of iterations for the spectral norm estimation
#define PARAM_IT_SPEC_NORM 6
// Hamming weight of the tags
#define PARAM_W 5
// Bound on the square spectral norm of R
#define PARAM_R_MAX_SQ_SPECTRAL_NORM 4909.69550475422420277027
// Gaussian parameter s_4 for v_2 and v_3
#define PARAM_S4 77.68853163270689776709
// Squared Gaussian parameter s_1^2
#define PARAM_S1SQ 29833689.83574411645531654358
// Squared Gaussian parameter s_2^2
#define PARAM_S2SQ 201183.59824153673253022134
// Squared Gaussian parameter s_3^2
#define PARAM_S3SQ 29632506.23750257864594459534
// Gaussian width for p_2 (sqrt(s_4^2 - s_G^2))
#define PARAM_SQRT_S4SQ_SGSQ 63.52563493692179008576
// Negated ratio -1/(1/s_G^2 - 1/s_4^2)
#define PARAM_SGINVSQ_S4INVSQ -2991.20481842988692733343
// Negated ratio -s_G^2/(s_4^2 - s_G^2)
#define PARAM_NEGSGSQ_DIV_S4SQ_SGSQ -0.49560117302052814070
// Negated gadget width -s_G^2
#define PARAM_NEG_SGSQ -2000.00165310704187504598
// Squared verification bound on v_11
#define PARAM_B11SQ 13015501725L
// Squared verification bound on v_12
#define PARAM_B12SQ 6041978260L
// Squared verification bound on v_2
#define PARAM_B2SQ 3365317L
// Squared verification bound on v_3
#define PARAM_B3SQ 954096L

// Length of the public and secret seeds
#define SEED_BYTES 32
// Length of the state
#define STATE_BYTES 64

#endif /* PARAMS_H */
