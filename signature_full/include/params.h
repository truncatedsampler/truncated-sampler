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
// Gadget base
#define PARAM_B 13
// Number of iterations for the spectral norm estimation
#define PARAM_IT_SPEC_NORM 5
// Hamming weight of the tags
#define PARAM_W 5
// Bound on the square spectral norm of R
#define PARAM_R_MAX_SQ_SPECTRAL_NORM 8483.65467712651479814667
// Gaussian parameter s_2 for v_2 and v_3
#define PARAM_S2 63.43242045540068119180
// Squared Gaussian parameter s_1^2
#define PARAM_S1SQ 34135443.48365909606218338013
// Gaussian width for p_2 (sqrt(s_2^2 - s_G^2))
#define PARAM_SQRT_S2SQ_SGSQ 44.98522325968487223236
// Negated ratio -1/(1/s_G^2 - 1/s_2^2)
#define PARAM_SGINVSQ_S2INVSQ -3976.61147401400103262858
// Negated ratio -s_G^2/(s_2^2 - s_G^2)
#define PARAM_NEGSGSQ_DIV_S2SQ_SGSQ -0.98830409356725157366
// Squared verification bound on v_1
#define PARAM_B1SQ 13059023256L
// Squared verification bound on v_2
#define PARAM_B2SQ 3631941L
// Squared verification bound on v_3
#define PARAM_B3SQ 1002342L

// Length of the public and secret seeds
#define SEED_BYTES 32
// Length of the state
#define STATE_BYTES 64

#endif /* PARAMS_H */
