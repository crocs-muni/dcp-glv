import dcp
import msm
import time
from sage.all import ZZ, RR, log
from copy import deepcopy
from utils import GuessVerification


class SACGuess:
    def __init__(self, k0=None, k1=None):
        self.k0 = [ZZ(b) for b in k0] if k0 is not None else []
        self.k1 = [ZZ(b) for b in k1] if k1 is not None else []

    def add(self, b0, b1):
        self.k0 = [ZZ(b0)] + self.k0
        self.k1 = [ZZ(b1)] + self.k1

    def values(self):
        return msm.from_bin(self.k0), msm.from_bin(self.k1)

    def is_positive(self):
        k0, k1 = self.values()
        return k0 >= 0 and k1 >= 0

    def __eq__(self, other):
        return self.k0 == other.k0 and self.k1 == other.k1

    def size(self):
        return len(self.k0)

    def __repr__(self):
        return f"SAC(k0={self.k0}, k1={self.k1})"

    def recode(self, l):
        return msm.two_dim_recoding_inverse_upper_bits(l, self.k0, self.k1)

    def next_guess(self):
        return [[1, 1], [1, 0], [-1, -1], [-1, 0]]


def verify_guess(s0, s1, g0, g1, zvpparams, i):
    dcp_scalars = 2 * s0 * g0 + 2 * s1 * g1 - 2 * s0 * g1, 2 * s1 * g0 - 2 * s0 * g1
    P = dcp.solve_glv_dcp_pari(dcp_scalars, zvpparams.glv, zvpparams.registers)
    if P is None:
        return GuessVerification.UNKNOWN
    P = ZZ(zvpparams.glv.order_field(g0 + g1 * zvpparams.glv.lam) ** (-1)) * P
    if msm.scalar_mul_oracle_positive(P, zvpparams, i)[1]:
        return GuessVerification.TRUE
    return GuessVerification.FALSE


def zvp_guess(i, scalar_guesses, zvpparams):
    """ZVP-GLV attack on the i-th upper bits of k0,k1"""
    new_scalar_guesses = []
    measured_guesses = []
    for sacguess in scalar_guesses:
        sack0, sack1 = sacguess.values()
        for next_guess in sacguess.next_guess():
            new_guess = deepcopy(sacguess)
            new_guess.add(*next_guess)
            oracle_output = verify_guess(sack0, sack1, *next_guess, zvpparams, i)
            if oracle_output == GuessVerification.UNKNOWN:
                if new_guess.is_positive():
                    new_scalar_guesses.append(new_guess)
                continue
            if oracle_output == GuessVerification.TRUE:
                measured_guesses.append(new_guess)
    if not measured_guesses:
        return new_scalar_guesses
    return measured_guesses


def zvp_attack(zvpparams):
    """ZVP-GLV attack on MSM"""
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    zvp_target = zvpparams.target_bits
    nguesses = [2]
    l = zvpparams.secrets.half_bits + 1
    b0, b1 = msm.twodim_recoding(msm.to_bin(k0, l), msm.to_bin(k1, l))
    scalar_guesses = [SACGuess([1], [0]), SACGuess([1], [1])]
    for i in range(l - 2, l - 2 - zvp_target + 1, -1):
        scalar_guesses = zvp_guess(i, scalar_guesses, zvpparams)
        assert SACGuess(b0[i:], b1[i:]) in scalar_guesses, (
            b0,
            b1,
            i,
            l,
            zvpparams.secrets.to_dict(),
            scalar_guesses,
        )
        if zvpparams.verbose:
            print("Number of guesses: ", len(scalar_guesses))
        nguesses.append(len(scalar_guesses))
    kipairs = []

    for sacguess in scalar_guesses:
        assert sacguess.size() == zvp_target
        k0bin, k1bin = sacguess.k0, sacguess.k1
        k0bin = [int(b) for b in k0bin]
        k1bin = [int(b) for b in k1bin]
        kipairs.append([k0bin, k1bin])
    return kipairs, [int(g) for g in nguesses]


def zvp_glv_sac(params):
    result = {}
    time_zvp = time.time()
    try:
        kipairs, nguesses = zvp_attack(params)
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
