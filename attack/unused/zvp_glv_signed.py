import dcp
import msm
import time
import utils
from sage.all import ZZ, RR, ceil, log, GF
from copy import deepcopy
from utils import GuessVerification


class SignedGuess:
    def __init__(self, k0=None, k1=None):
        self.k0 = [ZZ(b) for b in k0] if k0 is not None else []
        self.k1 = [ZZ(b) for b in k1] if k1 is not None else []

    def add(self, next_guess):
        b0, b1 = next_guess
        self.k0 = [ZZ(b0)] + self.k0
        self.k1 = [ZZ(b1)] + self.k1

    def values(self):
        return msm.from_regular_wnaf(self.k0,1), msm.from_regular_wnaf(self.k1,1)

    def is_positive(self):
        k0, k1 = self.values()
        return k0 >= 0 and k1 >= 0

    def __eq__(self, other):
        return self.k0 == other.k0 and self.k1 == other.k1

    def size(self):
        return len(self.k0)

    def recode(self):
        return msm.from_regular_wnaf(self.k0,1), msm.from_regular_wnaf(self.k1,1)

    def next_guess(self):
        return [[-1, -1], [1, 1], [1, -1], [-1, 1]]

    def __repr__(self):
        return f"Signed(k0={self.k0}, k1={self.k1})"


def zvp_glv_signed(params):
    """ZVP-GLV attack on Shamir's trick"""
    result = {}
    time_zvp = time.time()

    try:
        kipairs, nguesses = zvp_attack_signed_shamir(params)
    except Exception as e:
        params.error_log(e)
        raise e

    result["nguesses"] = nguesses
    result["recovered"] = float(
        max(ZZ(0), RR(2 * len(nguesses) - log(len(kipairs), 2)))
    )
    result["time_zvp"] = float(time.time() - time_zvp)
    result["scalars"] = kipairs
    return result


def verify_guess(s0, s1, next_guess, zvpparams, i):
    g0, g1 = next_guess
    P = None
    glv = zvpparams.glv
    s0, s1 = 2 * s0, 2 * s1
    if (s0, s1) != (0, 0):
        if (g0, g1) == (1, 1):
            dcp_scalars = s1, s1 - s0
            P = dcp.solve_glv_dcp_pari(dcp_scalars, glv, zvpparams.registers)
            if P is None:
                return GuessVerification.UNKNOWN
            P = ZZ(glv.order_field(g0 + g1 * glv.lam) ** (-1)) * P
        if (g0, g1) == (-1, -1):
            dcp_scalars = -s1, s0 - s1
            P = dcp.solve_glv_dcp_pari(dcp_scalars, glv, zvpparams.registers)
            if P is None:
                return GuessVerification.UNKNOWN
            P = ZZ(glv.order_field(g0 + g1 * glv.lam) ** (-1)) * P
        if (g0, g1) == (1, -1):
            dcp_scalars = 2 * s0 - s1, s0 + s1
            P = dcp.solve_glv_multi_dcp_pari(dcp_scalars, 3, glv, zvpparams.registers)
            if P is None:
                return GuessVerification.UNKNOWN
            P = ZZ(glv.order_field(g0 + g1 * glv.lam) ** (-1)) * P
            P = 3 * P
        if (g0, g1) == (-1, 1):
            dcp_scalars = -2 * s0 + s1, -s0 - s1
            P = dcp.solve_glv_multi_dcp_pari(dcp_scalars, 3, glv, zvpparams.registers)
            if P is None:
                return GuessVerification.UNKNOWN
            P = ZZ(glv.order_field(g0 + g1 * glv.lam) ** (-1)) * P
            P = 3 * P
    else:
        return GuessVerification.UNKNOWN
    if msm.signed_shamir_scalar_mul_oracle_positive(P, zvpparams, i, w=1)[1]:
        return GuessVerification.TRUE
    return GuessVerification.FALSE


def zvp_guess_signed_shamir(i, scalar_guesses, zvpparams):
    """ZVP-GLV attack on the i-th upper bits of k0,k1"""
    new_scalar_guesses = []
    measured_guesses = []
    for binguess in scalar_guesses:
        sack0, sack1 = binguess.values()
        for next_guess in binguess.next_guess():
            new_guess = deepcopy(binguess)
            new_guess.add(next_guess)
            oracle_output = verify_guess(sack0, sack1, next_guess, zvpparams, i)
            if oracle_output == GuessVerification.UNKNOWN:
                new_scalar_guesses.append(new_guess)
                continue
            if oracle_output == GuessVerification.TRUE:
                measured_guesses.append(new_guess)
                continue
    if not measured_guesses:
        return new_scalar_guesses
    return measured_guesses


def zvp_attack_signed_shamir(zvpparams):
    """ZVP-GLV attack on MSM"""
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    zvp_target = zvpparams.target_bits

    scalar_guesses = [
        SignedGuess([1], [1]),
        SignedGuess([1], [-1]),
        SignedGuess([-1], [1]),
    ]
    nguesses = [len(scalar_guesses)]
    l = zvpparams.secrets.half_bits
    b0, b1 = msm.regular_wnaf(k0, 1), msm.regular_wnaf(k1, 1)
    for i in range(l - 2, l - 1 - zvp_target, -1):
        scalar_guesses = zvp_guess_signed_shamir(i, scalar_guesses, zvpparams)
        s = SignedGuess(b0[i:], b1[i:])
        assert s in scalar_guesses
        if zvpparams.verbose:
            print("Number of guesses: ", len(scalar_guesses))
        nguesses.append(len(scalar_guesses))
    kipairs = []
    for winguess in scalar_guesses:
        v0, v1 = winguess.k0, winguess.k1
        v00, v01 = msm.fromcutwnaf(v0)
        v10, v11 = msm.fromcutwnaf(v1)
        kipairs.extend([(v00, v10), (v00, v11), (v01, v10), (v00, v11)])
    kipairs = list(set(kipairs))
    return [[int(k0), int(k1)] for k0, k1 in kipairs], [int(g) for g in nguesses]
