import dcp
import msm
import time
import utils
from sage.all import ZZ, RR, ceil, log
from copy import deepcopy
from utils import GuessVerification


class BinaryGuess:
    def __init__(self, k0=None, k1=None):
        self.k0 = [ZZ(b) for b in k0] if k0 is not None else []
        self.k1 = [ZZ(b) for b in k1] if k1 is not None else []

    def add(self, next_guess):
        b0, b1 = next_guess
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

    def recode(self):
        return msm.from_bin(self.k0), msm.from_bin(self.k1)

    def next_guess(self):
        return [[0, 0], [0, 1], [1, 0], [1, 1]]


def zvp_glv_binary(params):
    """ZVP-GLV attack on Shamir's trick"""
    result = {}
    time_zvp = time.time()
    kipairs, nguesses = zvp_attack_shamir(params)

    result["nguesses"] = nguesses
    result["recovered"] = float(
        max(ZZ(0), RR(2 * len(nguesses) - log(len(kipairs), 2)))
    )
    result["time_zvp"] = float(time.time() - time_zvp)
    result["scalars"] = kipairs
    return result


def verify_guess(s0, s1, next_guess, zvpparams, i):
    g0, g1 = next_guess
    if g0 == 0 and g1 == 0:
        return GuessVerification.UNKNOWN
    P = None
    if (s0, s1) != (0, 0):
        dcp_scalars = (g0 - g1) * 2 * s0 + g1 * 2 * s1, g0 * 2 * s1 - g1 * 2 * s0
        P = dcp.solve_glv_dcp_pari(dcp_scalars, zvpparams.glv, zvpparams.registers)
    if P is None:
        return GuessVerification.UNKNOWN
    P = ZZ(zvpparams.glv.order_field(g0 + g1 * zvpparams.glv.lam) ** (-1)) * P
    if msm.shamir_scalar_mul_oracle_positive(P, zvpparams, i)[1]:
        return GuessVerification.TRUE
    return GuessVerification.FALSE


def zvp_guess_shamir(i, scalar_guesses, zvpparams):
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


def zvp_attack_shamir(zvpparams):
    """ZVP-GLV attack on MSM"""
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    zvp_target = zvpparams.target_bits

    scalar_guesses = [
        BinaryGuess([0], [0]),
        BinaryGuess([1], [0]),
        BinaryGuess([0], [1]),
        BinaryGuess([1], [1]),
    ]
    nguesses = [len(scalar_guesses)]
    l = zvpparams.secrets.half_bits
    for i in range(l - 2, l - 1 - zvp_target, -1):
        scalar_guesses = zvp_guess_shamir(i, scalar_guesses, zvpparams)
        s = BinaryGuess(msm.to_bin(k0, l)[i:], msm.to_bin(k1, l)[i:])
        assert s in scalar_guesses
        if zvpparams.verbose:
            print("Number of guesses: ", len(scalar_guesses))
        nguesses.append(len(scalar_guesses))
    kipairs = []
    for binguess in scalar_guesses:
        assert binguess.size() == zvp_target
        kipairs.append(binguess.recode())
    kipairs = list(set([(a, b) for a, b in kipairs]))
    return [[int(k0), int(k1)] for k0, k1 in kipairs], [int(g) for g in nguesses]
