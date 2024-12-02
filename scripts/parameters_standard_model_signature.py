# -*- coding: utf-8 -*-
"""
    Brief: Main parameters classes
"""
#-- Import --#
from math import sqrt, exp, log, log2, floor, ceil, pi, e
from scipy.special import comb
from scipy.optimize import root
from sympy import isprime, prevprime
from estimate_SIS_LWE import *
#-- End Import --#
#-- Global Parameters --#
COST_MODEL      = 'sieving'
QUANTUM         = False
ROUGH           = True
LOG2_EPS        = -40
SPECTRAL_SLACK  = 6
#-- End Global Parameters --#

def c_star(dim:int, sec:float):
    """
    Find the Gaussian tailcut rate c so that the upper bound c*s*sqrt(dim) is
    verified with probability at least 1 - 2^(-sec)
    - input:
        (int)   dim     -- Dimension
        (int)   sec     -- Security parameter
    - output:
        (flt)   c_star  -- Tailcut rate
    """
    f = lambda c: sec + dim * (log2(c * sqrt(2*pi)) + (1/2 - pi*c**2)*log2(e))
    return root(f, 1)['x'][0]

class Full_Gadget_Parameters:
    """
    Main class containing all the parameters for the 
    signature scheme. 
    """

    def __init__(self, target_bitsec:int, n:int, d:int, k:int, q_start:int, no_guessing=False):
        """
        Computing all the signature parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   k               -- Gadget dimension 
            (int)   q_start         -- Starting search modulus (q largest adequately splitted prime below q_start)
            (bool)  no_guessing     -- False if accounting for guessing loss in security estimation
        - output:
            Full_Gadget_Parameters object with the following attributes
                [+] Parameters
                    (int)   security-- Reached bit security
                    (int)   n       -- Ring degree
                    (int)   d       -- Module rank
                    (int)   q       -- Modulus
                    (int)   k       -- Gadget vector dimension
                    (int)   b       -- Gadget base
                    (flt)   s_G     -- Gadget sampling width
                    (flt)   s_1     -- Top preimage sampling width
                    (flt)   s_2     -- Bottom preimage sampling width
                    (int)   w       -- Hamming weight of tags
                    (int)   kappa   -- Number of modulus splitting factor
                    (int)   Q       -- Maximal number of queries
                    (flt)   smoothing_ndk -- Smoothing of Z^(ndk)
                    (flt)   square_bound_sk -- Squared spectral norm bound on sk
                    (flt)   alpha   -- Rejection sampling slack
                    (flt)   M       -- Rejection sampling repetition rate
                    (flt)   B_1     -- Verification bound for v_1
                    (flt)   B_1_s   -- Verification bound for v_1 squared
                    (flt)   B_2     -- Verification bound for v_2
                    (flt)   B_2_s   -- Verification bound for v_2 squared
                    (flt)   B_3     -- Verification bound for v_3
                    (flt)   B_3_s   -- Verification bound for v_3 squared
                [+] Key Recovery Security
                    (int)   mlwe_blocksize  -- Required BKZ blocksize
                    (flt)   mlwe_rhf        -- Root Hermite Factor
                    (int)   mlwe_coresvp_c  -- Classical CoreSVP hardness
                    (int)   mlwe_coresvp_q  -- Quantum CoreSVP hardness
                [+] Forgery Security (for x = I or x = II)
                    (int)   msis_beta_inf_x -- M-SIS Infinity norm bound
                    (flt)   msis_beta_x     -- M-SIS Euclidean norm bound
                    (int)   msis_subdim_x   -- Optimal M-SIS subdimension
                    (flt)   msis_rhf_x      -- Root Hermite Factor
                    (int)   msis_coresvp_c_x-- Classical CoreSVP hardness
                    (int)   msis_coresvp_q_x-- Quantum CoreSVP hardness
                [+] Forgery Security (for x = I or x = II)
                    (flt)   eps_x           -- Hardness bound of type-x forgery
                [+] Sizes
                    (int)   pk_bitsize  -- Bitsize of the public key (B = AR)
                    (int)   sk_bitsize  -- Bitsize of the secret key (R)
                    (int)   sig_bitsize -- Bitsize of the signature (tag + v_{1,2} + v_2 + v_3)
        """
        ### Parameters

        # Degree of the ring R = Z[X]/<X^n + 1>
        self.n = n

        # M-SIS module rank 
        self.d = d

        # Maximal number of signature queries (hardcoded)
        self.Q = 2 ** 32

        # Number of splitting factors for modulus (hardcoded)
        self.kappa = 8

        # Finding the largest prime modulus splitting in kappa factors below q_start
        q = prevprime(q_start + 1)
        while (q%(4*self.kappa) != 2*self.kappa + 1) or (q < sqrt(self.kappa) ** self.kappa):
            q = prevprime(q)
        self.q = q

        # Finding minimum Hamming weight for the tag space
        w = 1
        while comb(self.n, w) < self.Q:
            w += 1
        self.w = w

        # Gadget dimension
        self.k = k

        # Gadget base for the gadget vector g = [1 b b^2 ... b^(k-1)]
        self.b = ceil(q ** (1/self.k))

        # Smoothing parameters
        self.smoothing_ndk = sqrt((log(2 * self.n * self.d * self.k) - LOG2_EPS * log(2)) / pi)

        # Temporary variables: 
        #   - bound on spectral norm of R
        #   - bound on euclidean norm of U*m
        norm_R = 3/4 * (sqrt(2 * self.d * self.n) + sqrt(self.k * self.d * self.n) + SPECTRAL_SLACK)
        norm_sm = self.n * sqrt(self.d)
        self.square_bound_sk = norm_R ** 2

        # Computing Gaussian widths for elliptic sampler
        self.s_G = self.smoothing_ndk * sqrt(self.b ** 2 + 1)
        self.s_1 = sqrt(2) * self.smoothing_ndk * (self.b + 1/self.b) * norm_R
        self.s_2 = sqrt(2) * self.smoothing_ndk * (self.b + 1/self.b) 

        # Computing rejection sampling parameters for security proof
        self.alpha = self.s_1 / norm_sm
        self.M = exp(pi / self.alpha**2)
        if self.M > 2**5:
            print('[Warning] Rejection Sampling incurs more than 5-bit loss')

        ### Security

        # Computing hardness of M-LWE_{n,d,d,q,B_1} for signer key recovery
        Xs = ND.CenteredBinomial(1)
        Xe = ND.CenteredBinomial(1)
        res = estimate_LWE(n = self.n * self.d,
                            q = self.q,
                            Xs = Xs,
                            Xe = Xe,
                            m = self.n * self.d,
                            cost_model = COST_MODEL,
                            rough = ROUGH)
        self.mlwe_blocksize = floor(res[0])
        self.mlwe_rhf = res[1]
        self.mlwe_coresvp_c = res[2]
        self.mlwe_coresvp_q = res[3]

        # Computing M-SIS security for type I+II forgeries
        tail_cut_prob_log = 20
        tail_cut_prob = 2**(-tail_cut_prob_log)
        # Square verification bounds
        self.B_1_s = floor(c_star(2*self.d*self.n, tail_cut_prob_log)**2 * self.s_1**2 * (2*self.d*self.n))
        self.B_1 = sqrt(self.B_1_s)

        self.B_2_s = floor(c_star(self.k*self.d*self.n, tail_cut_prob_log)**2 * self.s_2**2 * (self.k*self.d*self.n))
        self.B_2 = sqrt(self.B_2_s)

        self.B_3_s = floor(c_star(self.k*self.n, tail_cut_prob_log)**2 * self.s_2**2 * (self.k*self.n))
        self.B_3 = sqrt(self.B_3_s)


        # M-SIS hardness
        self.msis_beta_I = sqrt(floor( (self.B_1 + sqrt(self.d*self.n) * self.B_2)**2 + self.B_3_s + self.n + 1 ))
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (2 * self.d + self.k + 2),
                        q = self.q,
                        beta = self.msis_beta_I,
                        cost_model = COST_MODEL)
        self.msis_subdim_I = res[0]
        self.msis_rhf_I = res[1]
        self.msis_bkz_I = res[2]
        self.msis_coresvp_c_I = res[3]
        self.msis_coresvp_q_I = res[4]

        self.msis_beta_II = sqrt(floor( (2*self.B_1 + sqrt(self.d*self.n) * sqrt(4*self.B_2_s + self.n))**2 + 4*self.B_3_s ))
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (2 * self.d + self.k),
                        q = self.q,
                        beta = self.msis_beta_II,
                        cost_model = COST_MODEL)
        self.msis_subdim_II = res[0]
        self.msis_rhf_II = res[1]
        self.msis_bkz_II = res[2]
        self.msis_coresvp_c_II = res[3]
        self.msis_coresvp_q_II = res[4]

        # Computing security for type I+II forgeries
        self.tag_space_size = comb(self.n, self.w)
        # M-LWE and M-SIS loss
        mlwe_hardness_bound = 2 ** (- (self.mlwe_coresvp_q if QUANTUM else self.mlwe_coresvp_c))
        msis_I_hardness_bound = 2 ** (- (self.msis_coresvp_q_I if QUANTUM else self.msis_coresvp_c_I))
        msis_II_hardness_bound = 2 ** (- (self.msis_coresvp_q_II if QUANTUM else self.msis_coresvp_c_II))
        Ek = self.k * mlwe_hardness_bound

        # Trapdoor switching loss
        a = 1 - 1/(2 * target_bitsec)
        o = 2 * target_bitsec
        eps = 2 ** LOG2_EPS
        delta = ((1+eps)/(1-eps))**(4*self.n*self.d + 4) * ((1+eps/(self.n*self.d*self.k))/(1-eps/(self.n*self.d*self.k)))**(2*self.n*self.d*self.k)
        loss_TS = (1 + o*(o-1)/(2*(2-delta)**(o+1)) * (delta-1)**2)**(self.Q/o)

        h = lambda Adv: Ek + loss_TS * (2*Ek + loss_TS * (Ek + Adv)**a)**a # Equation (2)
        C = 2

        # Computing loss through hybrid argument
        self.eps_I = C * msis_I_hardness_bound * (1 if no_guessing else (self.tag_space_size - self.Q))
        self.eps_II = C * msis_II_hardness_bound * (1 if no_guessing else self.Q)
        for _ in range(self.d):
            self.eps_I = h(self.eps_I)
            self.eps_II = h(self.eps_II)
        self.eps_II = mlwe_hardness_bound + 2*self.M*C/(1-tail_cut_prob)**3*(1+eps)/(1-eps)*self.eps_II

        self.security = -log2(2*max(self.eps_I, self.eps_II))
        
        ### Efficiency

        # Size of public key: |pk| = |B| (+ seed)
        self.pk_bitsize = self.d * self.k * self.d * self.n * ceil(log2(self.q))

        # Size of secret key: |sk| = |R| (perturbation sampling material not included)
        self.sk_bitsize = 2 * self.d * self.k * self.d * self.n * ceil(log2(3))

        # Size of signature: |sig| = |tag| + |v_{1,2}| + |v_2| + |v_3|
        self.tag_bitsize    = self.n
        self.v_12_bitsize   = ceil(self.n * self.d * (1/2 + log2(self.s_1)))
        self.v_2_bitsize    = ceil(self.k * self.d * self.n * (1/2 + log2(self.s_2)))
        self.v_3_bitsize    = ceil(self.k * self.n * (1/2 + log2(self.s_2)))
        
        self.sig_bitsize    = self.tag_bitsize + self.v_12_bitsize + self.v_2_bitsize + self.v_3_bitsize       
    
    def __repr__(self):
        """
        Printing a SEP_Parameters object
        """
        tmp = '\n[+] Signature Scheme Parameters\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Ring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Gadget dimension', 'k', self.k)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Gadget base', 'b', self.b)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Top preimage sampling width', 's_1', self.s_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Bottom preimage sampling width', 's_2', self.s_2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Hamming weights of tags', 'w', self.w)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack (security proof)', 'α', self.alpha)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate (security proof)', 'M', self.M)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Maximal number of queries', 'Q', floor(log2(self.Q)))
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Smoothing loss', 'ε', LOG2_EPS)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z^(ndk)', 'η(ndk)', self.smoothing_ndk)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 1', 'B_1', self.B_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 2', 'B_2', self.B_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 3', 'B_3', self.B_3)
        tmp += 100 * '=' + '\n'
        tmp += '\n[+] M-SIS and M-LWE hardness\n'
        tmp += 100 * '=' + '\n'
        tmp += '{:-^100s}'.format('M-SIS Hardness (I)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β (I)', self.msis_beta_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim (I)', self.msis_subdim_I)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0 (I)', self.msis_rhf_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β (I)', self.msis_bkz_I)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec (I)', self.msis_coresvp_c_I)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec (I)', self.msis_coresvp_q_I)
        tmp += '{:-^100s}'.format('M-SIS Hardness (II)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β (II)', self.msis_beta_II)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim (II)', self.msis_subdim_II)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0 (II)', self.msis_rhf_II)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β (II)', self.msis_bkz_II)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec (II)', self.msis_coresvp_c_II)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec (II)', self.msis_coresvp_q_II)
        tmp += '{:-^100s}'.format('M-LWE Hardness') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec', self.mlwe_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.mlwe_coresvp_q)
        tmp += '{:-^100s}'.format('Security') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Type-I Security', 'eps_I', -log2(self.eps_I))
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Type-II Security', 'eps_II', -log2(self.eps_II))
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Forgery Security', 'F-Sec', self.security)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Key Recovery Security', 'KR-Sec', self.mlwe_coresvp_q if QUANTUM else self.mlwe_coresvp_c)
        tmp += '\n[+] Signature Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Public key size (KB)', '|pk|', self.pk_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Secret key size (KB)', '|sk|', self.sk_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Tag size (B)', '|tag|', self.tag_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_{1,2} size (B)', '|v_{1,2}|', self.v_12_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_2 size (B)', '|v_2|', self.v_2_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_3 size (B)', '|v_3|', self.v_3_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Overall Signature size (B)', '|sig|', self.sig_bitsize / 2 ** 3.)
        return tmp


class Truncated_Gadget_Parameters:
    """
    Main class containing all the parameters for the 
    signature scheme. 
    """

    def __init__(self, target_bitsec:int, n:int, d:int, k:int, l:int, q_start:int, no_guessing=False):
        """
        Computing all the signature parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   k               -- Gadget dimension 
            (int)   l               -- Number of dropped gadget entries
            (int)   q_start         -- Starting search modulus (q largest adequately splitted prime below q_start)
        - output:
            SEP_Parameters object with the following attributes
                [+] Parameters
                    (int)   security-- Reached bit security
                    (int)   n       -- Ring degree
                    (int)   d       -- Module rank
                    (int)   q       -- Modulus
                    (int)   k       -- Gadget vector dimension
                    (int)   l       -- Number of dropped gadget entries
                    (int)   b       -- Gadget base
                    (flt)   s_G     -- Gadget sampling width
                    (flt)   s_11    -- Final Gaussian width of v_11
                    (flt)   s_1     -- Preimage sampling width 1
                    (flt)   s_2     -- Preimage sampling width 2
                    (flt)   s_3     -- Preimage sampling width 3
                    (flt)   s_4     -- Preimage sampling width 4
                    (int)   w       -- Hamming weight of tags
                    (int)   kappa   -- Number of modulus splitting factor
                    (int)   Q       -- Maximal number of queries
                    (flt)   smoothing_ndk -- Smoothing of Z^(ndk)
                    (flt)   square_bound_sk -- Squared spectral norm bound on sk
                    (flt)   M_11    -- Rejection sampling repetition rate for v_11
                    (flt)   M_12    -- Rejection sampling repetition rate for v_12
                    (flt)   B_11    -- Verification bound for v_11
                    (flt)   B_11_s  -- Verification bound for v_11 squared
                    (flt)   B_12    -- Verification bound for v_12
                    (flt)   B_12_s  -- Verification bound for v_12 squared
                    (flt)   B_2     -- Verification bound for v_2
                    (flt)   B_2_s   -- Verification bound for v_2 squared
                    (flt)   B_3     -- Verification bound for v_3
                    (flt)   B_3_s   -- Verification bound for v_3 squared
                [+] Key Recovery Security
                    (int)   mlwe_blocksize  -- Required BKZ blocksize
                    (flt)   mlwe_rhf        -- Root Hermite Factor
                    (flt)   mlwe_coresvp_c  -- Classical CoreSVP hardness
                    (flt)   mlwe_coresvp_q  -- Quantum CoreSVP hardness
                [+] M-SIS hardness (for x = I or x = II)
                    (flt)   msis_beta_x     -- M-SIS Euclidean norm bound
                    (int)   msis_subdim_x   -- Optimal M-SIS subdimension
                    (flt)   msis_rhf_x      -- Root Hermite Factor
                    (flt)   msis_coresvp_c_x-- Classical CoreSVP hardness
                    (flt)   msis_coresvp_q_x-- Quantum CoreSVP hardness
                [+] Forgery Security (for x = I or x = II)
                    (flt)   eps_x           -- Hardness bound of type-x forgery
                [+] Sizes
                    (int)   pk_bitsize  -- Bitsize of the public key (B = AR)
                    (int)   sk_bitsize  -- Bitsize of the secret key (R)
                    (int)   sig_bitsize -- Bitsize of the signature (tag + v_{1,2} + v_2 + v_3)
        """
        ### Parameters

        # Degree of the ring R = Z[X]/<X^n + 1>
        self.n = n

        # M-SIS module rank 
        self.d = d

        # Maximal number of signature queries (hardcoded)
        self.Q = 2 ** 32

        # Number of splitting factors for modulus (hardcoded)
        self.kappa = 8

        # Finding the largest prime modulus splitting in kappa factors below q_start
        q = prevprime(q_start + 1)
        while (q%(4*self.kappa) != 2*self.kappa + 1) or (q < sqrt(self.kappa) ** self.kappa):
            q = prevprime(q)
        self.q = q

        # Finding minimum Hamming weight for the tag space
        w = 1
        while comb(self.n, w) < self.Q:
            w += 1
        self.w = w

        # Gadget dimension
        self.k = k

        # Gadget base for the gadget vector g = [1 b b^2 ... b^(k-1)]
        self.b = ceil(q ** (1/self.k))

        # Number of dropped gadget entries 
        self.l = l

        # Smoothing parameters
        self.smoothing_ndk = sqrt((log(2 * self.n * self.d * self.k) - LOG2_EPS * log(2)) / pi)

        # Temporary variables: 
        #   - bound on spectral norm of R
        #   - bound on euclidean norm of U*m
        norm_Ri = 3/4 * (sqrt(self.d * self.n) + sqrt((self.k - self.l) * self.d * self.n) + SPECTRAL_SLACK)
        norm_si_m = self.n * sqrt(self.d / 2) 
        self.square_bound_sk = norm_Ri ** 2

        # Computing Gaussian widths for elliptic sampler
        self.s_G = self.smoothing_ndk * sqrt(self.b ** 2 + 1)
        self.s_1 = self.smoothing_ndk * (self.b + 1/self.b) * sqrt(4*self.w**2 + 3*norm_Ri**2)
        self.s_2 = self.smoothing_ndk * (self.b + 1/self.b) * 2*self.w
        self.s_3 = self.smoothing_ndk * (self.b + 1/self.b) * sqrt(3) * norm_Ri
        self.s_4 = self.smoothing_ndk * (self.b + 1/self.b) * sqrt(3) 

        self.s_11 = self.smoothing_ndk * (self.b + 1/self.b) * sqrt((self.b**(2*self.l) - 1)/(self.b**2 - 1)*4*self.w**2 + 3*norm_Ri**2)

        # Computing rejection sampling parameters for security proof
        self.alpha_11 = self.s_11 / (norm_si_m)
        self.alpha_12 = self.s_3 / (norm_si_m)
        self.M_11 = exp(pi / self.alpha_11 ** 2)
        self.M_12 = exp(pi / self.alpha_12 ** 2)
        if self.M_11 * self.M_12 > 2**5:
            print('[Warning] Rejection Sampling incurs more than 5-bit loss')

        ### Security

        # Computing hardness of M-LWE_{n,d,d,q,B_1} for signer key recovery
        Xs = ND.CenteredBinomial(1)
        Xe = ND.CenteredBinomial(1)
        res = estimate_LWE(n = self.n * self.d,
                            q = self.q,
                            Xs = Xs,
                            Xe = Xe,
                            m = self.n * self.d,
                            cost_model = COST_MODEL,
                            rough = ROUGH)
        self.mlwe_blocksize = floor(res[0])
        self.mlwe_rhf = res[1]
        self.mlwe_coresvp_c = res[2]
        self.mlwe_coresvp_q = res[3]

        # Computing M-SIS security for type I+II forgeries
        tail_cut_prob_log = 20
        tail_cut_prob = 2**(-tail_cut_prob_log)
        # Square verification bounds
        self.B_11_s = floor(c_star(self.d*self.n, tail_cut_prob_log)**2 * self.s_11**2 * (self.d*self.n))
        self.B_11 = sqrt(self.B_11_s)

        self.B_12_s = floor(c_star(self.d*self.n, tail_cut_prob_log)**2 * self.s_3**2 * (self.d*self.n))
        self.B_12 = sqrt(self.B_12_s)

        self.B_2_s = floor(c_star((self.k-self.l)*self.d*self.n, tail_cut_prob_log)**2 * self.s_4**2 * ((self.k-self.l)*self.d*self.n))
        self.B_2 = sqrt(self.B_2_s)

        self.B_3_s = floor(c_star((self.k-self.l)*self.n, tail_cut_prob_log)**2 * self.s_4**2 * ((self.k-self.l)*self.n))
        self.B_3 = sqrt(self.B_3_s)

        # M-SIS hardness
        self.msis_beta_I = sqrt(floor( (sqrt(self.B_11_s + self.B_12_s) + sqrt(self.d * self.n) * self.B_2) ** 2 + self.B_3_s + self.n + 1 ))
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (2 * self.d + self.k - self.l + 2),
                        q = self.q,
                        beta = self.msis_beta_I,
                        cost_model = COST_MODEL)
        self.msis_subdim_I = res[0]
        self.msis_rhf_I = res[1]
        self.msis_bkz_I = res[2]
        self.msis_coresvp_c_I = res[3]
        self.msis_coresvp_q_I = res[4]

        self.msis_beta_II = sqrt(floor( (2 * sqrt(self.B_11_s + self.B_12_s) + sqrt(self.d * self.n) * sqrt(4*self.B_2_s + self.n)) ** 2 + 4*self.B_3_s ))
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (2 * self.d + self.k - self.l),
                        q = self.q,
                        beta = self.msis_beta_II,
                        cost_model = COST_MODEL)
        self.msis_subdim_II = res[0]
        self.msis_rhf_II = res[1]
        self.msis_bkz_II = res[2]
        self.msis_coresvp_c_II = res[3]
        self.msis_coresvp_q_II = res[4]

        # Computing security for type I+II forgeries
        self.tag_space_size = comb(self.n, self.w)
        # M-LWE and M-SIS loss
        mlwe_hardness_bound = 2 ** (- (self.mlwe_coresvp_q if QUANTUM else self.mlwe_coresvp_c))
        msis_I_hardness_bound = 2 ** (- (self.msis_coresvp_q_I if QUANTUM else self.msis_coresvp_c_I))
        msis_II_hardness_bound = 2 ** (- (self.msis_coresvp_q_II if QUANTUM else self.msis_coresvp_c_II))
        Ek = (self.k-self.l) * mlwe_hardness_bound

        # Trapdoor switching loss
        a = 1 - 1/(2 * target_bitsec)
        o = 2 * target_bitsec
        eps = 2 ** LOG2_EPS
        delta = ((1+eps)/(1-eps))**(2*(self.n-1)*self.d*(self.l + 1) + 4*self.d + 4) * ((1+eps/(self.n*self.d*self.k))/(1-eps/(self.n*self.d*self.k)))**(2*self.n*self.d*self.k)
        loss_TS = (1 + o*(o-1)/(2*(2-delta)**(o+1)) * (delta-1)**2)**(self.Q/o)

        h = lambda Adv: Ek + loss_TS * (2*Ek + loss_TS * (Ek + Adv)**a)**a # Equation (2)
        C = 2

        # Computing loss through hybrid argument
        self.eps_I = C * msis_I_hardness_bound * (1 if no_guessing else (self.tag_space_size - self.Q))
        self.eps_II = C * msis_II_hardness_bound * (1 if no_guessing else self.Q)
        for _ in range(self.d):
            self.eps_I = h(self.eps_I)
            self.eps_II = h(self.eps_II)
        self.eps_II = mlwe_hardness_bound + 4*self.M_11*self.M_12*C**2 / (1-tail_cut_prob)**4 *(1+eps)/(1-eps)*self.eps_II

        self.security = -log2(2*max(self.eps_I, self.eps_II))

        ### Efficiency

        # Size of public key: |pk| = |B| (+ seed)
        self.pk_bitsize = self.d * (self.k-self.l) * self.d * self.n * ceil(log2(self.q))

        # Size of secret key: |sk| = |R| (perturbation sampling material not included)
        self.sk_bitsize = 2 * self.d * (self.k-self.l) * self.d * self.n * ceil(log2(3))

        # Size of signature: |sig| = |tag| + |v_{1,2}| + |v_2| + |v_3|
        self.tag_bitsize    = self.n
        self.v_12_bitsize   = ceil(self.n * self.d * (1/2 + log2(self.s_3)))
        self.v_2_bitsize    = ceil((self.k-self.l) * self.d * self.n * (1/2 + log2(self.s_4)))
        self.v_3_bitsize    = ceil((self.k-self.l) * self.n * (1/2 + log2(self.s_4)))
        
        self.sig_bitsize    = self.tag_bitsize + self.v_12_bitsize + self.v_2_bitsize + self.v_3_bitsize       
    
    def __repr__(self):
        """
        Printing a SEP_Parameters object
        """
        tmp = '\n[+] Signature Scheme Parameters\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Ring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Gadget dimension', 'k', self.k)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Gadget base', 'b', self.b)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gadget sampling width', 's_G', self.s_G)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Preimage sampling width 1', 's_1', self.s_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Preimage sampling width 2', 's_2', self.s_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Preimage sampling width 3', 's_3', self.s_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Preimage sampling width 4 (v_2,v_3)', 's_4', self.s_4)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Final preimage sampling width v_11', 's_11', self.s_11)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Hamming weights of tags', 'w', self.w)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack (security proof)', 'α_11', self.alpha_11)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack (security proof)', 'α_12', self.alpha_12)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate (security proof)', 'M_11', self.M_11)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate (security proof)', 'M_12', self.M_12)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Maximal number of queries', 'Q', floor(log2(self.Q)))
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Smoothing loss', 'ε', LOG2_EPS)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z^(ndk)', 'η(ndk)', self.smoothing_ndk)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 11', 'B_11', self.B_11)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 12', 'B_12', self.B_12)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 2', 'B_2', self.B_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 3', 'B_3', self.B_3)
        tmp += '\n[+] M-SIS and M-LWE hardness\n'
        tmp += 100 * '=' + '\n'
        tmp += '{:-^100s}'.format('Type-I Forgeries (M-SIS)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β (I)', self.msis_beta_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim (I)', self.msis_subdim_I)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0 (I)', self.msis_rhf_I)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β (I)', self.msis_bkz_I)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec (I)', self.msis_coresvp_c_I)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec (I)', self.msis_coresvp_q_I)
        tmp += '{:-^100s}'.format('Type-II Forgeries (M-SIS)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β (II)', self.msis_beta_II)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim (II)', self.msis_subdim_II)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0 (II)', self.msis_rhf_II)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β (II)', self.msis_bkz_II)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec (II)', self.msis_coresvp_c_II)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec (II)', self.msis_coresvp_q_II)
        tmp += '{:-^100s}'.format('Signer Key Recovery (M-LWE)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec', self.mlwe_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.mlwe_coresvp_q)
        tmp += '{:-^100s}'.format('Security') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Type-I Security', 'eps_I', -log2(self.eps_I))
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Type-II Security', 'eps_II', -log2(self.eps_II))
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Forgery Security', 'F-Sec', self.security)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Key Recovery Security', 'KR-Sec', self.mlwe_coresvp_q if QUANTUM else self.mlwe_coresvp_c)
        tmp += '\n[+] Signature Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Public key size (KB)', '|pk|', self.pk_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Secret key size (KB)', '|sk|', self.sk_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Tag size (B)', '|tag|', self.tag_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_{1,2} size (B)', '|v_{1,2}|', self.v_12_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_2 size (B)', '|v_2|', self.v_2_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('v_3 size (B)', '|v_3|', self.v_3_bitsize / 2 ** 3.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Overall Signature size (B)', '|sig|', self.sig_bitsize / 2 ** 3.)
        return tmp

sig_full = Full_Gadget_Parameters(target_bitsec=128, n=256, d=4, k=5, q_start=ceil(2 ** 18.5), no_guessing=False)
print(sig_full)
sig_trunc = Truncated_Gadget_Parameters(target_bitsec=128, n=256, d=4, k=5, l=2, q_start=ceil(2 ** 18.5), no_guessing=False)
print(sig_trunc)