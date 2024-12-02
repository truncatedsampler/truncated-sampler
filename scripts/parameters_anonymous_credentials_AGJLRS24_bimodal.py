# -*- coding: utf-8 -*-
"""
    Brief: Main parameters classes
"""
#-- Import --#
from math import sqrt, exp, log, log2, floor, ceil, pi, e
from scipy.special import comb
from scipy.optimize import root
from sympy import isprime, prevprime, nextprime, divisors
from estimate_SIS_LWE import *
#-- End Import --#
#-- Global Parameters --#
COST_MODEL      = 'realistic_sieving'
ROUGH           = True
QUANTUM         = False
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

class SEP_Parameters:
    """
    Main class containing all the parameters for the 
    signature scheme. 
    """

    def __init__(self, target_bitsec:int, n:int, d:int, m:int, q_start:int):
        """
        Computing all the signature parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   m               -- Message dimension (without secret key)
            (int)   q_start         -- Starting search modulus (q largest adequately splitted prime below q_start)
        - output:
            SEP_Parameters object with the following attributes
                [+] Parameters
                    (int)   sec     -- Target bit security or security parameter
                    (int)   n       -- Ring degree
                    (int)   d       -- Module rank
                    (int)   m_s     -- User secret dimension (2d)
                    (int)   m       -- Secret (2d) + Message dimension
                    (int)   q       -- Modulus
                    (int)   k       -- Gadget vector dimension
                    (int)   b       -- Gadget base
                    (flt)   s_G     -- Gadget sampling width
                    (flt)   s_1     -- Top preimage sampling width
                    (flt)   s_2     -- Bottom preimage sampling width
                    (int)   w       -- Hamming weight of tags
                    (int)   kappa   -- Number of modulus splitting factor
                    (int)   Q       -- Maximal number of queries
                    (flt)   r       -- Smoothing of Z
                    (flt)   r_ndk   -- Smoothing of Z^(ndk)
                    (flt)   r_ndk2  -- Smoothing of Z^(nd(2+k))
                    (flt)   alpha   -- Rejection sampling slack
                    (flt)   M       -- Rejection sampling repetition rate
                    (flt)   B_1     -- Verification bound for v_1
                    (flt)   B_1_s   -- Verification bound for v_1 squared
                    (flt)   B_1_prime -- Verification bound for v_1 - r (hiding case)
                    (flt)   B_1_prime_s -- Verification bound for v_1 - r (hiding case) squared
                    (flt)   B_2     -- Verification bound for v_2
                    (flt)   B_2_s   -- Verification bound for v_2 squared
                    (flt)   B_3     -- Verification bound for v_3
                    (flt)   B_3_s   -- Verification bound for v_3 squared
                [+] Required M-SIS security (for x = I or x = II)
                    (int)   req_msis_coresvp_x  -- Required M-SIS CoreSVP hardness
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
                [+] Sizes
                    (int)   pk_bitsize  -- Bitsize of the public key (B = AR)
                    (int)   sk_bitsize  -- Bitsize of the secret key (R)
                    (int)   sig_bitsize -- Bitsize of the signature (tag + v_{1,2} + v_2 + v_3)
        """
        ### Parameters

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R = Z[X]/<X^n + 1>
        self.n = n

        # M-SIS module rank 
        self.d = d

        # User secret key dimensionality (2*d)
        self.m_s = 2 * self.d

        # Message dimensionality including user secret key (2*d)
        self.m = 2 * self.d + m

        # Maximal number of signature queries (hardcoded)
        self.Q = 2 ** 32

        # Number of splitting factors for modulus (hardcoded)
        self.kappa = 2

        # Finding the largest prime modulus splitting in kappa factors below q_start
        q = prevprime(q_start + 1)
        while (q%(4*self.kappa) != 2*self.kappa + 1) or (q < (2*sqrt(self.kappa)) ** self.kappa):
            q = prevprime(q)
        self.q = q

        # Finding minimum Hamming weight for the tag space
        w = 1
        while comb(self.n, w) < self.Q:
            w += 1
        self.w = w

        # Gadget base for the gadget vector g = [1 b b^2 ... b^(k-1)]
        self.b = ceil(q ** (1/5))

        # Gadget dimension
        self.k = ceil(log(self.q) / log(self.b))

        # Smoothing parameters
        self.r = sqrt((log(2) - LOG2_EPS * log(2)) / pi)
        self.r_ndk = sqrt((log(2 * self.n * self.d * self.k) - LOG2_EPS * log(2)) / pi)
        self.r_ndk2 = sqrt((log(2 * self.n * self.d * (2 + self.k)) - LOG2_EPS * log(2)) / pi)

        # Temporary variables: 
        #   - bound on spectral norm of R
        #   - bound on euclidean norm of U*m
        norm_R = 7/10 * (sqrt(2 * self.d * self.n) + sqrt(self.k * self.d * self.n) + SPECTRAL_SLACK)
        self.bound_sk = norm_R
        self.square_bound_sk = norm_R ** 2
        norm_Um = sqrt(self.d * self.n) * sqrt(self.m * self.n)

        # Computing Gaussian widths for elliptic sampler
        self.s_G = self.r_ndk * sqrt(self.b ** 2 + 1)
        s_MP12 = sqrt(2 * self.s_G ** 4 / (self.s_G ** 2 - 1)) * norm_R
        self.s_1 = max( sqrt(pi / log(2)) * (norm_Um + sqrt(2 * self.d * self.n)), s_MP12)
        self.s_2 = sqrt(2 * self.s_G ** 2 + self.r_ndk2 ** 2)

        # Computing rejection sampling parameters for security proof
        self.alpha = self.s_1 / (norm_Um + sqrt(2 * self.d * self.n))
        self.M = exp(pi / self.alpha ** 2)

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
        self.mlwe_S_blocksize = floor(res[0])
        self.mlwe_S_rhf = res[1]
        self.mlwe_S_coresvp_c = res[2]
        self.mlwe_S_coresvp_q = res[3]

        # Computing hardness of M-LWE_{n,d,d,q,U(S_1)} for user key recovery
        Xs = ND.UniformMod(self.q)
        Xe = ND.Uniform(0,1)
        res = estimate_LWE(n = self.n * self.d,
                            q = self.q,
                            Xs = Xs,
                            Xe = Xe,
                            m = 2 * self.n * self.d,
                            cost_model = COST_MODEL,
                            rough = ROUGH)
        self.mlwe_U_blocksize = floor(res[0])
        self.mlwe_U_rhf = res[1]
        self.mlwe_U_coresvp_c = res[2]
        self.mlwe_U_coresvp_q = res[3]

        # Computing M-SIS security for type I+II forgeries
        # Square verification bounds
        self.B_1_s = floor(c_star(2 * self.d * self.n, self.sec + 3) ** 2 * self.s_1 ** 2 * (2 * self.d * self.n))
        self.B_1 = sqrt(self.B_1_s)
        self.B_1_prime_s = floor((self.B_1 + sqrt(2 * self.d * self.n)) ** 2)
        self.B_1_prime = sqrt(self.B_1_prime_s)
        self.B_2_s = floor(c_star(self.k * self.d * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.k * self.d * self.n))
        self.B_2 = sqrt(self.B_2_s)
        self.B_3_s = floor(c_star(self.k * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.k * self.n))
        self.B_3 = sqrt(self.B_3_s)

        self.msis_beta_I = sqrt(floor( (self.B_1_prime + sqrt(self.d * self.n) * self.B_2) ** 2 + self.B_3 ** 2 + self.m * self.n + 1 ))
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (2 * self.d + self.k + self.m + 1),
                        q = self.q,
                        beta = self.msis_beta_I,
                        cost_model = COST_MODEL)
        self.msis_subdim_I = res[0]
        self.msis_rhf_I = res[1]
        self.msis_bkz_I = res[2]
        self.msis_coresvp_c_I = res[3]
        self.msis_coresvp_q_I = res[4]

        self.msis_beta_II = sqrt(floor( (2 * self.B_1_prime + 2 * sqrt(self.d * self.n) * self.B_2 + norm_Um) ** 2 + 4 * self.B_3 ** 2 ))
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

        # M-ISIS hardness
        self.misis_beta = sqrt(floor(2*self.n*self.d))
        res = estimate_SIS(n = self.n * self.d,
                        m = 2 * self.n * self.d,
                        q = self.q,
                        beta = self.misis_beta,
                        cost_model = COST_MODEL)
        self.misis_subdim = 2 * self.n * self.d - res[0]
        self.misis_rhf = res[1]
        self.misis_bkz = res[2]
        self.misis_coresvp_c = res[3]
        self.misis_coresvp_q = res[4]

        # Computing security for type I+II forgeries
        self.tag_space_size = comb(self.n, self.w)
        # M-LWE and M-SIS loss
        mlwe_S_hardness_bound = 2 ** (- (self.mlwe_S_coresvp_q if QUANTUM else self.mlwe_S_coresvp_c))
        msis_I_hardness_bound = 2 ** (- (self.msis_coresvp_q_I if QUANTUM else self.msis_coresvp_c_I))
        msis_II_hardness_bound = 2 ** (- (self.msis_coresvp_q_II if QUANTUM else self.msis_coresvp_c_II))
        Ek = self.k * mlwe_S_hardness_bound
        Em = self.m * mlwe_S_hardness_bound

        # Trapdoor switching loss
        a = 1 - 1/(2 * self.sec)
        o = 2 * self.sec
        eps = 2 ** LOG2_EPS
        delta = ((1+eps)/(1-eps))**(4*self.n*self.d + 4) * ((1+eps/(self.n*self.d*self.k))/(1-eps/(self.n*self.d*self.k)))**(2*self.n*self.d*self.k)
        loss_TS = (1 + o*(o-1)/(2*(2-delta)**(o+1)) * (delta-1)**2)**(self.Q/o)

        h = lambda Adv: Ek + loss_TS * (2*Ek + loss_TS * (Ek + Adv)**a)**a # Equation (2)
        C = 2

        # Computing loss through hybrid argument
        self.eps_I = C * (self.tag_space_size - self.Q) * msis_I_hardness_bound
        self.eps_II = C**2 * self.Q * msis_II_hardness_bound
        for _ in range(self.d):
            self.eps_I = h(self.eps_I)
            self.eps_II = h(self.eps_II)
        self.eps_II = Em + 2*self.M*C*(1+eps)/(1-eps)*self.eps_II

        ### Efficiency

        # Size of public key: |pk| = |B| (+ seed)
        self.pk_bitsize = self.d * self.k * self.d * self.n * ceil(log2(self.q)) + 256

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
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Ring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('User secret dimension (2d)', 'm_s', self.m_s)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Overall message dimension (secret+message)', 'm', self.m)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Gadget dimension', 'k', self.k)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Gadget base', 'b', self.b)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gadget sampling width', 's_G', self.s_G)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Top preimage sampling width', 's_1', self.s_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Bottom preimage sampling width', 's_2', self.s_2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Hamming weights of tags', 'w', self.w)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack (security proof)', 'α', self.alpha)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate (security proof)', 'M', self.M)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Maximal number of queries', 'Q', floor(log2(self.Q)))
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Smoothing loss', 'ε', LOG2_EPS)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z', 'η(1)', self.r)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z^(ndk)', 'η(ndk)', self.r_ndk)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Smoothing of Z^(nd(2+k))', 'η(nd(2+k))', self.r_ndk2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 1', 'B_1', self.B_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Signature verification bound 1 (hiding)', 'B_1\'', self.B_1_prime)
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
        tmp += '{:-^100s}'.format('Impersonation Forgeries (M-ISIS)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-ISIS solution', 'β', self.misis_beta)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-ISIS subdimension', 'sdim', self.misis_subdim)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.misis_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.misis_bkz)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec', self.misis_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.misis_coresvp_q)
        tmp += '{:-^100s}'.format('Signer Key Recovery (M-LWE)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_S_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_S_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec', self.mlwe_S_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.mlwe_S_coresvp_q)
        tmp += '{:-^100s}'.format('User Key Recovery (M-LWE)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_U_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_U_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical CoreSVP hardness', 'CSec', self.mlwe_U_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum CoreSVP hardness', 'QSec', self.mlwe_U_coresvp_q)
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

class Issue_ZKP_Parameters:
    """
    Main class containing all the parameters for the 
    zero-knowledge proof system (Show protocol). 
    """

    def __init__(self, target_bitsec:int, n:int, d:int, m_2:int, q_1_start:int, n_attr:int, gamma:int, D:int, sig, bimodal=True, compression=True, garbage=True):
        """
        Computing all the zero-knowledge argument parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   m_2             -- Commitment randomness dimension
            (int)   q_1_start       -- Starting search modulus (q_1 largest adequately splitted prime below q_1_start)
            (int)   gamma           -- Compression parameter for commitment w
            (int)   D               -- Compression parameter for commitment t_A
            (SEP_Parameters) sig    -- Signature parameters
            (PKE_Parameters) pke    -- Verifiable encryption parameters
        - output:
            Show_ZKP_Parameters object with the following attributes
                [+] Parameters
                    (int)   sec     -- Target bit security or security parameter
                    (int)   n       -- Ring degree
                    (int)   d       -- Module rank
                    (int)   q_1     -- Modulus factor
                    (int)   q_2     -- Second modulus factor (signature modulus)
                    (int)   q_min   -- Smallest modulus factor
                    (int)   q       -- Modulus
                    (int)   l       -- Parameter for soundness amplification
                    (int)   m_1     -- Dimension of witness
                    (int)   m_2     -- Dimension of ABDLOP commitment randomness
                    (int)   xi_s2   -- Infinity norm of ABDLOP commitment randomness
                    (int)   rho     -- Infinity norm of challenges
                    (int)   eta     -- Manhattan-like norm of challenges
                    (int)   challenge_space_size -- Size of challenge space
                    (int)   gamma   -- Compression parameter for commitment w
                    (int)   D       -- Compression parameter for commitment t_A
                    For x = 1, 2, 3
                        (flt)   M_x         -- Rejection sampling repetition rate for z_x
                        (flt)   log_eps_x   -- Rejection sampling (log) loss for z_x
                        (flt)   alpha_x     -- Rejection sampling slack for z_x
                        (flt)   s_x         -- Gaussian width for y_x
                    (flt)   bound_witness -- Euclidean norm bound on the witness
                    (flt)   B_arp     -- Euclidean norm bound from approximate range proofs
                [+] Zero-Knowledge Security
                    (int)   mlwe_blocksize  -- Required BKZ blocksize
                    (flt)   mlwe_rhf        -- Root Hermite Factor
                    (int)   mlwe_coresvp_c  -- Classical hardness bound (-log)
                    (int)   mlwe_coresvp_q  -- Quantum hardness bound (-log)
                [+] Soundness Security
                    (flt)   msis_beta       -- M-SIS Euclidean norm bound
                    (int)   msis_subdim     -- Optimal M-SIS subdimension
                    (flt)   msis_rhf        -- Root Hermite Factor
                    (int)   msis_coresvp_c  -- Classical hardness bound (-log)
                    (int)   msis_coresvp_q  -- Quantum hardness bound (-log)
                    (flt)   soundness_error -- Soundness error
                [+] Sizes
                    (int)   incompressible_bitsize  -- Bitsize of incompressible elements (in R_q)
                    (int)   compressible_bitsize    -- Bitsize of compressible elements (Gaussians)
                    (int)   challenge_bitsize       -- Bitsize of challenge (c)
                    (int)   proof_bitsize           -- Bitsize of proof (π)
        """
        ### Parameters

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R'' = Z[X]/<X^n' + 1>
        self.n = n

        # M-SIS module rank 
        self.d = d

        # Number of disclosed attributes
        self.n_attr = n_attr

        # # Number of splitting factors for modulus (hardcoded)
        self.kappa = 2

        # Finding the largest prime modulus splitting in kappa factors below q_start and such that a divisor of q-1 is close to gamma
        if compression:
            found_q1_gamma = False
            q_1 = q_1_start + 1
            while not(found_q1_gamma):
                q_1 = prevprime(q_1)
                while q_1 % 8 != 5:
                    q_1 = prevprime(q_1)
                q = sig.q * q_1
                divs = divisors(q-1)
                for div in divs:
                    if (gamma <= div <= 5*gamma/4) and (div % 2 == 0): # find even divisor closest to gamma (not larger than 5.gamma/4)
                        self.gamma = div
                        found_q1_gamma = True 
                        break
            # Commitment compression
            self.D = D 
        else:
            self.gamma = 0
            self.D = 0
            q_1 = prevprime(q_1_start + 1)
            while (q_1%(4*self.kappa) != 2*self.kappa + 1) or (q_1 < (2*sqrt(self.kappa)) ** self.kappa):
                q_1 = prevprime(q_1)
        self.q_1 = q_1
        self.q_2 = sig.q
        self.q = self.q_1 * sig.q
        self.q_min = min(self.q_1, sig.q)

        if garbage:
            # Repetition for soundness amplification
            l = ceil(self.sec / log2(self.q_min))
            self.l = ceil(l/2)
        else:
            self.l = ceil(self.sec / log2(self.q_min))

        # Infinity norm bound on the challenges
        self.rho = ceil(1/2 * (2 ** (2*(self.sec + 1)/self.n) - 1))

        # Manhattan-like norm bound on the challenges (hardcoded)
        self.eta = {64:93, 128:42, 256:37, 512:57, 1024:84}[self.n]

        # Size of challenge space
        self.challenge_space_size = (2 * self.rho + 1) ** (self.n // 2) / 2

        # Subring gap
        self.k = sig.n // self.n

        # Witness dimension
        self.m_1 = 2 * sig.d * self.k # commitment randomness r (binary)
        self.m_1 += (sig.m - self.n_attr) * self.k # msg (include usk)
        if bimodal:
            self.m_1 += 1 # bimodal bit

        # Commitment randomness dimension and infinity norm bound (hardcoded)
        self.m_2 = m_2
        self.xi_s2 = 1

        # Bound on Euclidean norm of the witness
        self.bound_witness = sqrt(self.n * self.m_1)

        # Rejection sampling parameters (hardcoded)
        self.M_1 = 2
        self.M_2 = 2
        self.M_3 = 2
        if bimodal:
            self.eps_1 = 0
            self.eps_2 = 0
            self.eps_3 = 0
            self.alpha_1 = sqrt(pi/log(self.M_1))
            self.alpha_2 = sqrt(pi/log(self.M_2))
            self.alpha_3 = sqrt(pi/log(self.M_3))
        else:
            self.eps_1 = 2**(-130)
            self.eps_2 = 2**(-130)
            self.eps_3 = 2**(-130)
            self.alpha_1 = sqrt(pi)/log(self.M_1) * (sqrt(log(1/self.eps_1) + log(self.M_1)) + sqrt(log(1/self.eps_1)))
            self.alpha_2 = sqrt(pi)/log(self.M_2) * (sqrt(log(1/self.eps_2) + log(self.M_2)) + sqrt(log(1/self.eps_2)))
            self.alpha_3 = sqrt(pi)/log(self.M_3) * (sqrt(log(1/self.eps_3) + log(self.M_3)) + sqrt(log(1/self.eps_3)))

        # Gaussian widths
        self.s_1 = self.alpha_1 * self.eta * self.bound_witness
        self.s_2 = self.alpha_2 * self.eta * self.xi_s2 * sqrt(self.n * self.m_2)
        self.s_3 = self.alpha_3 * sqrt(337) * self.bound_witness

        # Checking approximate range proofs bounds
        self.B_256_s = floor(c_star(256, self.sec + 3) ** 2 * self.s_3 ** 2 * 256)
        self.B_256 = sqrt(self.B_256_s)
        self.B_arp = self.B_256 * 2/sqrt(26)
        self.cond_1_bound = 41 * self.n * self.m_1 * self.B_arp # Condition for modular JL bound for ARP bound
        self.cond_2_bound = self.B_arp**2 + self.B_arp * sqrt(self.n*self.m_1) # bound for proving witness binary
        x = (self.q_1*pow(self.q_1, -1, self.q_2) - self.q_2*pow(self.q_2, -1, self.q_1)) % self.q
        if x > self.q/2: x -= self.q
        self.cond_3_bound = abs(x)

        assert (self.q > self.cond_1_bound) and (self.q > self.cond_2_bound) and (self.cond_3_bound > self.B_arp)

        # Square Verification bounds
        self.B_1_s = 4 * floor(c_star(self.m_1 * self.n, self.sec + 3) ** 2 * self.s_1 ** 2 * (self.m_1 * self.n))
        if self.D != 0 and self.gamma != 0:
            self.B_2_s = floor((2 * sqrt(floor(c_star(self.m_2 * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.m_2 * self.n))) + (2**self.D * self.eta + self.gamma)*sqrt(self.n*self.d)) ** 2)
        else:
            self.B_2_s = 4 * floor(c_star(self.m_2 * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.m_2 * self.n))

        ### Security

        # Computing hardness of M-LWE for zero-knowledge
        Xs = ND.CenteredBinomial(self.xi_s2)
        Xe = ND.CenteredBinomial(self.xi_s2)
        res = estimate_LWE(n = self.n * (self.m_2 - (self.d + floor(256/self.n) + self.l + 1)),
                            q = self.q,
                            Xs = Xs,
                            Xe = Xe,
                            m = self.n * self.m_2,
                            cost_model = COST_MODEL,
                            rough = ROUGH)
        self.mlwe_blocksize = floor(res[0])
        self.mlwe_rhf = res[1]
        self.mlwe_coresvp_c = res[2]
        self.mlwe_coresvp_q = res[3]

        # Computing M-SIS security for soundness
        self.msis_beta = 4 * self.eta * sqrt(self.B_1_s + self.B_2_s)
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (self.m_1 + self.m_2),
                        q = self.q,
                        beta = self.msis_beta,
                        cost_model = COST_MODEL)
        self.msis_subdim = res[0]
        self.msis_rhf = res[1]
        self.msis_bkz = res[2]
        self.msis_coresvp_c = res[3]
        self.msis_coresvp_q = res[4]

        self.soundness_error = self.q_min ** (-self.n/self.kappa) + 2 / self.challenge_space_size + \
                                    2 ** (-(self.msis_coresvp_q if QUANTUM else self.msis_coresvp_c))
        if garbage:
            self.soundness_error += self.q_min ** (-2*self.l)
        else:
            self.soundness_error += self.q_min ** (-self.l)

        ### Efficiency

        # CRS size
        ajtai_crs_size = self.d * (self.m_1 + self.m_2) * (self.n * ceil(log2(self.q)))
        bdlop_crs_size = (256 / self.n + self.l + 1 + 1) * self.m_2 * (self.n * ceil(log2(self.q)))
        self.crs_size  = ajtai_crs_size + bdlop_crs_size

        # Proof size
        self.size_z1 = ceil(self.n * self.m_1 * (1/2 + log2(self.s_1)))
        self.size_z2 = ceil(self.n * (self.m_2-(self.d if self.D != 0 else 0)) * (1/2 + log2(self.s_2)))
        self.size_z3 = ceil(256 * (1/2 + log2(self.s_3)))
        self.size_c  = self.n * ceil(log2(2 * self.rho + 1))
        self.size_tA = self.n * self.d * (ceil(log2(self.q)) - self.D)
        self.size_hints = ceil(self.n * self.d * (2.25 if self.D != 0 else 0))
        self.size_tB = self.n * (256 / self.n + self.l) * ceil(log2(self.q)) 
        self.size_h  = self.n * self.l * ceil(log2(self.q))
        self.size_t1 = self.n * ceil(log2(self.q))

        # Size of incompressible elements (those uniform in R_q)
        self.incompressible_bitsize = self.size_tA + self.size_hints + self.size_tB + self.size_h + self.size_t1

        # Size of compressible elements (those Gaussians): 
        self.compressible_bitsize = self.size_z1 + self.size_z2 + self.size_z3

        # Size of proof
        self.proof_bitsize = self.incompressible_bitsize + self.compressible_bitsize + self.size_c
    
    def __repr__(self):
        """
        Printing a Issue_ZKP_Parameters object
        """
        tmp = '\n[+] ISSUANCE ZERO-KNOWLEDGE PROOF PARAMETERS\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Subring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Subring gap factor (n_signature / n_proof)', 'k', self.k)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus factor', 'q_1', self.q_1)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Smallest modulus factor', 'q_min', self.q_min)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Parameter for soundness amplification', 'l', self.l)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Dimension of witness', 'm_1', self.m_1)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Dimension of ABDLOP commitment randomness', 'm_2', self.m_2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Infinity norm of commitment randomness', 'ξ', self.xi_s2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Infinity norm of challenges', 'ρ', self.rho)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Manhattan-like norm of challenges', 'η', self.eta)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Size of challenge space', '|C|', floor(log2(self.challenge_space_size)))        
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Compression parameter 1', 'gamma', self.gamma)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Compression parameter 2', 'D', self.D)        
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 1', 'M_1', self.M_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 2', 'M_2', self.M_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 3', 'M_3', self.M_3)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 1', 'ε_1', 0 if self.eps_1 == 0 else int(log2(self.eps_1)))
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 2', 'ε_2', 0 if self.eps_2 == 0 else int(log2(self.eps_2)))
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 3', 'ε_3', 0 if self.eps_3 == 0 else int(log2(self.eps_3)))
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 1 (security proof)', 'α_1', self.alpha_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 2 (security proof)', 'α_2', self.alpha_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 3 (security proof)', 'α_3', self.alpha_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_1', 'σ_1', self.s_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_2', 'σ_2', self.s_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_3', 'σ_3', self.s_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on the witness', 'B', self.bound_witness)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Norm bound from approximate range proof', 'B_arp', self.B_arp)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 2: 41*n*m_1*B_arp', 'Cond. 1', floor(self.cond_1_bound))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 4: for binary witness', 'Cond. 2', floor(self.cond_2_bound))
        tmp += '| {:60s} | {:^10s} | 2^{:<18.5f} |\n'.format('Soundness error', 'δ_s', log2(self.soundness_error))

        tmp += '\n[+] M-SIS and M-LWE hardness\n'
        tmp += 100 * '=' + '\n'
        tmp += '{:-^100s}'.format('Soundness (M-SIS)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β', self.msis_beta)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim',self.msis_subdim)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.msis_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.msis_bkz)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical hardness bound (-log)', 'CSec', self.msis_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum hardness bound (-log)', 'QSec', self.msis_coresvp_q)
        tmp += '{:-^100s}'.format('Zero-Knowledge (M-LWE)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical hardness bound (-log)', 'CSec', self.mlwe_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum hardness bound (-log)', 'QSec', self.mlwe_coresvp_q)
        tmp += '\n[+] Proof Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Common Random String (KB)', '|crs|', self.crs_size / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Incompressible elements (KB)', '|π_1|', self.incompressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Compressible elements (KB)', '|π_2|', self.compressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Challenge (KB)', '|π_3|', self.size_c / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Total proof bitsize (KB)', '|π|', self.proof_bitsize / 2 ** 13.)
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of z1 (KB)', '|z_1|', self.size_z1 / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of z2 (KB)', '|z_2|', self.size_z2 / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of z3 (KB)', '|z_3|', self.size_z3 / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of c (KB)', '|c|', self.size_c / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of tA (KB)', '|tA|', self.size_tA / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of compression hints (KB)', '|hints|', self.size_hints / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of tB (KB)', '|tB|', self.size_tB / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of h (KB)', '|h|', self.size_h / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of t1 (KB)', '|t1|', self.size_t1 / 2 ** 13.)
        return tmp

class Show_ZKP_Parameters:
    """
    Main class containing all the parameters for the 
    zero-knowledge proof system (Show protocol). 
    """

    def __init__(self, target_bitsec:int, n:int, d:int, m_2:int, q_1_start:int, n_attr:int, gamma:int, D:int, sig, bimodal=True, compression=True, garbage=True):
        """
        Computing all the zero-knowledge argument parameters
        - input:
            (int)   target_bitsec   -- Target bit security or security parameter
            (int)   n               -- Ring degree
            (int)   d               -- Module rank
            (int)   m_2             -- Commitment randomness dimension
            (int)   q_1_start       -- Starting search modulus (q_1 largest adequately splitted prime below q_1_start)
            (int)   n_attr          -- Number of attribute to disclose
            (int)   gamma           -- Compression parameter for commitment w
            (int)   D               -- Compression parameter for commitment t_A
            (SEP_Parameters) sig    -- Signature parameters
        - output:
            Show_ZKP_Parameters object with the following attributes
                [+] Parameters
                    (int)   sec     -- Target bit security or security parameter
                    (int)   n       -- Ring degree
                    (int)   d       -- Module rank
                    (int)   q_1     -- Modulus factor
                    (int)   q_2     -- Second modulus factor (signature modulus)
                    (int)   q_min   -- Smallest modulus factor
                    (int)   q       -- Modulus
                    (int)   l       -- Parameter for soundness amplification
                    (int)   m_1     -- Dimension of witness
                    (int)   m_2     -- Dimension of ABDLOP commitment randomness
                    (int)   xi_s2   -- Infinity norm of ABDLOP commitment randomness
                    (int)   rho     -- Infinity norm of challenges
                    (int)   eta     -- Manhattan-like norm of challenges
                    (int)   challenge_space_size -- Size of challenge space
                    (int)   gamma   -- Compression parameter for commitment w
                    (int)   D       -- Compression parameter for commitment t_A
                    For x = 1, 2, 3
                        (flt)   M_x         -- Rejection sampling repetition rate for z_x
                        (flt)   log_eps_x   -- Rejection sampling (log) loss for z_x
                        (flt)   alpha_x     -- Rejection sampling slack for z_x
                        (flt)   s_x         -- Gaussian width for y_x
                    (flt)   bound_witness -- Euclidean norm bound on the witness
                    (flt)   B_arp     -- Euclidean norm bound from approximate range proofs
                [+] Zero-Knowledge Security
                    (int)   mlwe_blocksize  -- Required BKZ blocksize
                    (flt)   mlwe_rhf        -- Root Hermite Factor
                    (int)   mlwe_coresvp_c  -- Classical hardness bound (-log)
                    (int)   mlwe_coresvp_q  -- Quantum hardness bound (-log)
                [+] Soundness Security
                    (flt)   msis_beta       -- M-SIS Euclidean norm bound
                    (int)   msis_subdim     -- Optimal M-SIS subdimension
                    (flt)   msis_rhf        -- Root Hermite Factor
                    (int)   msis_coresvp_c  -- Classical hardness bound (-log)
                    (int)   msis_coresvp_q  -- Quantum hardness bound (-log)
                    (flt)   soundness_error -- Soundness error
                [+] Sizes
                    (int)   incompressible_bitsize  -- Bitsize of incompressible elements (in R_q)
                    (int)   compressible_bitsize    -- Bitsize of compressible elements (Gaussians)
                    (int)   challenge_bitsize       -- Bitsize of challenge (c)
                    (int)   proof_bitsize           -- Bitsize of proof (π)
        """
        ### Parameters

        # Security parameter
        self.sec = target_bitsec

        # Degree of the ring R'' = Z[X]/<X^n' + 1>
        self.n = n

        # M-SIS module rank 
        self.d = d

        # Number of disclosed attributes
        self.n_attr = n_attr

        # # Number of splitting factors for modulus (hardcoded)
        self.kappa = sig.kappa # 2

        if compression:
            # Finding the largest prime modulus splitting in kappa factors below q_start and such that a divisor of q-1 is close to gamma
            found_q1_gamma = False
            q_1 = q_1_start + 1
            while not(found_q1_gamma):
                q_1 = prevprime(q_1)
                while q_1 % 8 != 5:
                    q_1 = prevprime(q_1)
                q = sig.q * q_1
                divs = divisors(q-1)
                for div in divs:
                    if (gamma <= div <= 5*gamma/4) and (div % 2 == 0): # find even divisor closest to gamma (not larger than 5.gamma/4)
                        self.gamma = div
                        found_q1_gamma = True 
                        break
            # Commitment compression
            self.D = D 
        else:
            self.gamma = 0
            self.D = 0
            q_1 = prevprime(q_1_start + 1)
            while (q_1%(4*self.kappa) != 2*self.kappa + 1) or (q_1 < (2*sqrt(self.kappa)) ** self.kappa):
                q_1 = prevprime(q_1)
        self.q_1 = q_1
        self.q_2 = sig.q
        self.q = self.q_1 * sig.q
        self.q_min = min(self.q_1, sig.q)

        if garbage:
            # Repetition for soundness amplification
            l = ceil(self.sec / log2(self.q_min))
            self.l = ceil(l/2)
        else:
            self.l = ceil(self.sec / log2(self.q_min))

        # Infinity norm bound on the challenges
        self.rho = ceil(1/2 * (2 ** (2*(self.sec + 1)/self.n) - 1))

        # Manhattan-like norm bound on the challenges (hardcoded)
        self.eta = {64:93, 128:42, 256:37, 512:57, 1024:84}[self.n]

        # Size of challenge space
        self.challenge_space_size = (2 * self.rho + 1) ** (self.n // 2) / 2

        # Subring gap
        self.k = sig.n // self.n

        # Witness dimension
        self.m_1 = 2 * sig.d * self.k + 1 # v_1 (+ a_1 for four-square norm proof)
        self.m_1 += sig.d * sig.k * self.k + 1 # v_2 (+ a_2 for four-square norm proof)
        self.m_1 += sig.k * self.k + 1 # v_3 (+ a_3 for four-square norm proof)
        self.m_1 += self.k # tag
        self.m_1 += (sig.m - self.n_attr) * self.k # msg (include usk)
        if bimodal:
            self.m_1 += 1 # bimodal bit

        # Commitment randomness dimension and infinity norm bound (hardcoded)
        self.m_2 = m_2
        self.xi_s2 = 1

        # Bound on Euclidean norm of the witness
        self.bound_witness =  sqrt(sig.B_1_prime_s + sig.B_2_s + sig.B_3_s + sig.w + sig.n * (sig.m - self.n_attr) + (1 if bimodal else 0))

        # Rejection sampling parameters (hardcoded)
        self.M_1 = 2
        self.M_2 = 2
        self.M_3 = 2
        if bimodal:
            self.eps_1 = 0
            self.eps_2 = 0
            self.eps_3 = 0
            self.alpha_1 = sqrt(pi/log(self.M_1))
            self.alpha_2 = sqrt(pi/log(self.M_2))
            self.alpha_3 = sqrt(pi/log(self.M_3))
        else:
            self.eps_1 = 2**(-130)
            self.eps_2 = 2**(-130)
            self.eps_3 = 2**(-130)
            self.alpha_1 = sqrt(pi)/log(self.M_1) * (sqrt(log(1/self.eps_1) + log(self.M_1)) + sqrt(log(1/self.eps_1)))
            self.alpha_2 = sqrt(pi)/log(self.M_2) * (sqrt(log(1/self.eps_2) + log(self.M_2)) + sqrt(log(1/self.eps_2)))
            self.alpha_3 = sqrt(pi)/log(self.M_3) * (sqrt(log(1/self.eps_3) + log(self.M_3)) + sqrt(log(1/self.eps_3)))

        # Gaussian width
        self.s_1 = self.alpha_1 * self.eta * self.bound_witness
        self.s_2 = self.alpha_2 * self.eta * self.xi_s2 * sqrt(self.n * self.m_2)
        self.s_3 = self.alpha_3 * sqrt(337) * self.bound_witness

        # Checking approximate range proofs bounds
        self.B_256_s = floor(c_star(256, self.sec + 3) ** 2 * self.s_3 ** 2 * 256)
        self.B_256 = sqrt(self.B_256_s)
        self.B_arp = self.B_256 * 2/sqrt(26)
        self.cond_1_bound = max(sig.B_1_prime_s, sig.B_2_s, sig.B_3_s, sig.w, 1) # Lower bound -q < -B^2
        self.cond_2_bound = 41 * self.n * self.m_1 * self.B_arp # Condition for modular JL bound for ARP bound
        self.cond_3_bound = self.B_arp**2 - min(sig.B_1_prime_s, sig.B_2_s, sig.B_3_s, sig.w, 1) # Upper bound (2/root(26) * B_e)^2 - B^2 < q
        self.cond_4_bound = sig.w + sqrt(sig.w*self.n*self.k) # bound for proving tag binary
        self.cond_5_bound = self.B_arp**2 + self.B_arp * sqrt(self.n*self.k*(sig.m-self.n_attr)) # bound for msg binary
        x = (self.q_1*pow(self.q_1, -1, self.q_2) - self.q_2*pow(self.q_2, -1, self.q_1)) % self.q
        if x > self.q/2: x -= self.q
        self.cond_6_bound = abs(x)
        assert (self.q > self.cond_1_bound) and (self.q > self.cond_2_bound) and (self.q > self.cond_3_bound) and (self.q > self.cond_4_bound) and (self.q > self.cond_5_bound) and (self.cond_6_bound > self.B_arp)

        # Square Verification bounds
        self.B_1_s = 4 * floor(c_star(self.m_1 * self.n, self.sec + 3) ** 2 * self.s_1 ** 2 * (self.m_1 * self.n))
        if self.D != 0 and self.gamma != 0:
            self.B_2_s = floor((2 * sqrt(floor(c_star(self.m_2 * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.m_2 * self.n))) + (2**self.D * self.eta + self.gamma)*sqrt(self.n*self.d)) ** 2)
        else:
            self.B_2_s = 4 * floor(c_star(self.m_2 * self.n, self.sec + 3) ** 2 * self.s_2 ** 2 * (self.m_2 * self.n))

        ### Security

        # Computing hardness of M-LWE for zero-knowledge
        Xs = ND.CenteredBinomial(self.xi_s2)
        Xe = ND.CenteredBinomial(self.xi_s2)
        res = estimate_LWE(n = self.n * (self.m_2 - (self.d + floor(256/self.n) + self.l + 1)),
                            q = self.q,
                            Xs = Xs,
                            Xe = Xe,
                            m = self.n * self.m_2,
                            cost_model = COST_MODEL,
                            rough = ROUGH)
        self.mlwe_blocksize = floor(res[0])
        self.mlwe_rhf = res[1]
        self.mlwe_coresvp_c = res[2]
        self.mlwe_coresvp_q = res[3]

        # Computing M-SIS security for soundness
        self.msis_beta = 4 * self.eta * sqrt(self.B_1_s + self.B_2_s)
        res = estimate_SIS(n = self.n * self.d,
                        m = self.n * (self.m_1 + self.m_2),
                        q = self.q,
                        beta = self.msis_beta,
                        cost_model = COST_MODEL)
        self.msis_subdim = res[0]
        self.msis_rhf = res[1]
        self.msis_bkz = res[2]
        self.msis_coresvp_c = res[3]
        self.msis_coresvp_q = res[4]

        self.soundness_error = self.q_min ** (-self.n/2) + 2 / self.challenge_space_size + \
                                    2 ** (-(self.msis_coresvp_q if QUANTUM else self.msis_coresvp_c)) 
        if garbage:
            self.soundness_error += self.q_min ** (-2*self.l)
        else:
            self.soundness_error += self.q_min ** (-self.l)

        ### Efficiency

        # CRS size
        ajtai_crs_size = self.d * (self.m_1 + self.m_2) * (self.n * ceil(log2(self.q)))
        bdlop_crs_size = (256 / self.n + self.l + 1 + 1) * self.m_2 * (self.n * ceil(log2(self.q)))
        self.crs_size  = ajtai_crs_size + bdlop_crs_size

        # Proof size
        self.size_z1 = ceil(self.n * self.m_1 * (1/2 + log2(self.s_1)))
        self.size_z2 = ceil(self.n * (self.m_2-(self.d if self.D != 0 else 0)) * (1/2 + log2(self.s_2)))
        self.size_z3 = ceil(256 * (1/2 + log2(self.s_3)))
        self.size_c  = self.n * ceil(log2(2 * self.rho + 1))
        self.size_tA = self.n * self.d * (ceil(log2(self.q)) - self.D)
        self.size_hints = ceil(self.n * self.d * (2.25 if self.D != 0 else 0))
        self.size_tB = self.n * (256 / self.n + self.l) * ceil(log2(self.q))
        self.size_h  = self.n * self.l * ceil(log2(self.q))
        self.size_t1 = self.n * ceil(log2(self.q))

        # Size of incompressible elements (those uniform in R_q)
        self.incompressible_bitsize = self.size_tA + self.size_hints + self.size_tB + self.size_h + self.size_t1

        # Size of compressible elements (those Gaussians): 
        self.compressible_bitsize = self.size_z1 + self.size_z2 + self.size_z3

        # Size of proof
        self.proof_bitsize = self.incompressible_bitsize + self.compressible_bitsize + self.size_c
    
    def __repr__(self):
        """
        Printing a Show_ZKP_Parameters object
        """
        tmp = '\n[+] SHOW ZERO-KNOWLEDGE PROOF PARAMETERS\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Security parameter', 'λ', self.sec)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Subring degree of Z[X]/(X^n + 1)', 'n', self.n)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Subring gap factor (n_signature / n_proof)', 'k', self.k)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Module rank', 'd', self.d)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus factor', 'q_1', self.q_1)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Smallest modulus factor', 'q_min', self.q_min)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Modulus', 'q', self.q)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Parameter for soundness amplification', 'l', self.l)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Dimension of witness', 'm_1', self.m_1)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Dimension of ABDLOP commitment randomness', 'm_2', self.m_2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Infinity norm of commitment randomness', 'ξ', self.xi_s2)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Infinity norm of challenges', 'ρ', self.rho)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Manhattan-like norm of challenges', 'η', self.eta)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Size of challenge space', '|C|', floor(log2(self.challenge_space_size)))        
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Compression parameter 1', 'gamma', self.gamma)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Compression parameter 2', 'D', self.D)        
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 1', 'M_1', self.M_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 2', 'M_2', self.M_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling repetition rate 3', 'M_3', self.M_3)
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 1', 'ε_1', 0 if self.eps_1 == 0 else int(log2(self.eps_1)))
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 2', 'ε_2', 0 if self.eps_2 == 0 else int(log2(self.eps_2)))
        tmp += '| {:60s} | {:^10s} | 2^{:<18d} |\n'.format('Rejection sampling loss 3', 'ε_3', 0 if self.eps_3 == 0 else int(log2(self.eps_3)))
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 1 (security proof)', 'α_1', self.alpha_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 2 (security proof)', 'α_2', self.alpha_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Rejection sampling slack 3 (security proof)', 'α_3', self.alpha_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_1', 'σ_1', self.s_1)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_2', 'σ_2', self.s_2)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Gaussian width for y_3', 'σ_3', self.s_3)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on the witness', 'B', self.bound_witness)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Norm bound from approximate range proof', 'B_arp', self.B_arp)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 1: B_i^2', 'Cond. 1', floor(self.cond_1_bound))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 2: 41*n*m_1*B_arp', 'Cond. 2', floor(self.cond_2_bound))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 3: B_arp^2 - B_i^2', 'Cond. 3', floor(self.cond_3_bound))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 4: for binary tag', 'Cond. 4', floor(self.cond_4_bound))
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('ARP bound condition 4: for binary message', 'Cond. 5', floor(self.cond_5_bound))
        tmp += '| {:60s} | {:^10s} | 2^{:<18.5f} |\n'.format('Soundness error', 'δ_s', log2(self.soundness_error))

        tmp += '\n[+] M-SIS and M-LWE hardness\n'
        tmp += 100 * '=' + '\n'
        tmp += '{:-^100s}'.format('Soundness (M-SIS)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Euclidean norm bound on M-SIS solution', 'β', self.msis_beta)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Optimal M-SIS subdimension', 'sdim',self.msis_subdim)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.msis_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.msis_bkz)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical hardness bound (-log)', 'CSec', self.msis_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum hardness bound (-log)', 'QSec', self.msis_coresvp_q)
        tmp += '{:-^100s}'.format('Zero-Knowledge (M-LWE)') + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20d} |\n'.format('Required BKZ blocksize', 'BKZ-β', self.mlwe_blocksize)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Root Hermite Factor', 'δ_0', self.mlwe_rhf)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Classical hardness bound (-log)', 'CSec', self.mlwe_coresvp_c)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Quantum hardness bound (-log)', 'QSec', self.mlwe_coresvp_q)
        tmp += '\n[+] Proof Estimated Performance (KB)\n'
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Common Random String (KB)', '|crs|', self.crs_size / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Incompressible elements (KB)', '|π_1|', self.incompressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Compressible elements (KB)', '|π_2|', self.compressible_bitsize / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Challenge (KB)', '|π_3|', self.size_c / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Total proof bitsize (KB)', '|π|', self.proof_bitsize / 2 ** 13.)
        tmp += 100 * '=' + '\n'
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of z1 (KB)', '|z_1|', self.size_z1 / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of z2 (KB)', '|z_2|', self.size_z2 / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of z3 (KB)', '|z_3|', self.size_z3 / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of c (KB)', '|c|', self.size_c / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of tA (KB)', '|tA|', self.size_tA / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of compression hints (KB)', '|hints|', self.size_hints / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of tB (KB)', '|tB|', self.size_tB / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of h (KB)', '|h|', self.size_h / 2 ** 13.)
        tmp += '| {:60s} | {:^10s} | {:<20.5f} |\n'.format('Size of t1 (KB)', '|t1|', self.size_t1 / 2 ** 13.)

        return tmp

def estimate_anonymous_credentials(sig_pms, issue_zkp_pms, show_zkp_pms):
    """
    Estimate concrete security and efficiency of the anonymous credentials based on chosen parameters
    - input:
        (SEP_Parameters)        sig_pms         -- Signature parameters
        (Issue_ZKP_Parameters)  issue_zkp_pms   -- Issuance zero-knowledge proof parameters
        (Show_ZKP_Parameters)   show_zkp_pms    -- Show zero-knowledge proof parameters
    """
    # Print all parameters
    print(sig_pms)
    print(issue_zkp_pms)
    print(show_zkp_pms)

    # Compute issuance transcript sizes
    cmt_sz = sig_pms.n*sig_pms.d*ceil(log2(sig_pms.q))
    issproof_sz = issue_zkp_pms.proof_bitsize
    sig_sz = sig_pms.sig_bitsize
    iss_sz = cmt_sz + issproof_sz + sig_sz

    # Compute anonymous credential proof size
    showproof_sz = show_zkp_pms.proof_bitsize

    # Compute blindness security based on security proof 
    e_zk_issue = issue_zkp_pms.eps_1/issue_zkp_pms.M_1 + issue_zkp_pms.eps_2/issue_zkp_pms.M_2 + issue_zkp_pms.eps_3/issue_zkp_pms.M_3 + 2**(-issue_zkp_pms.mlwe_coresvp_c)
    e_zk_show = show_zkp_pms.eps_1/show_zkp_pms.M_1 + show_zkp_pms.eps_2/show_zkp_pms.M_2 + show_zkp_pms.eps_3/show_zkp_pms.M_3 + 2**(-show_zkp_pms.mlwe_coresvp_c)

    e_anonymity = 2 * e_zk_show
    security_anonymity = -log2(e_anonymity)

    # Compute one-more unforgeability security based on security proof 
    e_unforgeability = 3*(e_zk_issue + e_zk_show + issue_zkp_pms.soundness_error + 3*show_zkp_pms.soundness_error + 2**(-sig_pms.mlwe_U_coresvp_c) + sig_pms.tag_space_size*2**(-sig_pms.misis_coresvp_c) + sig_pms.eps_I + sig_pms.eps_II)
    security_unforgeability = -log2(e_unforgeability)

    tmp = '\n[ISSUANCE TRANSCRIPT SIZE]\n'
    tmp += '\n'
    tmp += '    |cmt|   =  %.3f KB        (user)\n' % (cmt_sz/2**13.)
    tmp += '  + |proof| =  %.3f KB       (user)\n' % (issproof_sz/2**13.)
    tmp += '  + |sig|   =  %.3f KB        (signer)\n' % (sig_sz/2**13.)
    tmp += '  _____________________\n'
    tmp += '  =  %.3f KB\n' % (iss_sz/2**13.)
    tmp += '\n'
    tmp += '\n[CREDENTIAL PROOF SIZE]\n'
    tmp += '\n'
    tmp += '    |proof|   =  %.3f KB\n' % (showproof_sz/2**13.)
    tmp += '\n'
    tmp += '[UNFORGEABILITY CONCRETE SECURITY] Achieved: %.3f\n' % (security_unforgeability)
    tmp += '[ANONYMITY CONCRETE SECURITY] Achieved: %.3f\n' % (security_anonymity)

    print(tmp)

ZK_OPTIMIZATIONS = True
sig_pms = SEP_Parameters(target_bitsec=128, n=256, d=4, m=10, q_start=ceil(2 ** 18.7))
issue_zkp_pms = Issue_ZKP_Parameters(target_bitsec=128, n=64, d=20, m_2=58, q_1_start=ceil(2 ** 19), n_attr=0, gamma=2**19, D=11, sig=sig_pms, bimodal=ZK_OPTIMIZATIONS, compression=ZK_OPTIMIZATIONS, garbage=ZK_OPTIMIZATIONS)
show_zkp_pms = Show_ZKP_Parameters(target_bitsec=128, n=64, d=23, m_2=74, q_1_start=ceil(2 ** 39), n_attr=0, gamma=2**30, D=22, sig=sig_pms, bimodal=ZK_OPTIMIZATIONS, compression=ZK_OPTIMIZATIONS, garbage=ZK_OPTIMIZATIONS)

estimate_anonymous_credentials(sig_pms, issue_zkp_pms, show_zkp_pms)