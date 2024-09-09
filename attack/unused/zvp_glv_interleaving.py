import dcp
import msm
import utils
import time
from sage.all import ZZ, RR, ceil, log
from copy import deepcopy
from utils import GuessVerification


class WNAFGuess:
    def __init__(self, w, k0=None, k1=None):
        self.w = w
        self.k0 = [ZZ(b) for b in k0] if k0 is not None else []
        self.k1 = [ZZ(b) for b in k1] if k1 is not None else []

    def add(self, b0, b1):
        self.k0 = [ZZ(b0)] + self.k0
        self.k1 = [ZZ(b1)] + self.k1

    def add_s(self, b, s):
        if s == 0:
            self.k0 = [ZZ(b)] + self.k0
        else:
            self.k1 = [ZZ(b)] + self.k1

    def zeros(self):
        c0, c1 = 0, 0
        for b in self.k0:
            if b != 0:
                break
            c0 += 1
        for b in self.k1:
            if b != 0:
                break
            c1 += 1
        return c0, c1

    def values(self):
        return msm.from_wnaf(self.k0), msm.from_wnaf(self.k1)

    def is_positive(self):
        k0, k1 = self.values()
        return k0 >= 0 and k1 >= 0

    def __eq__(self, other):
        return self.k0 == other.k0 and self.k1 == other.k1

    def sizes(self):
        return len(self.k0), len(self.k1)

    def __repr__(self):
        return f"WNAF(k0={self.k0}, k1={self.k1})"

    def next_guess_s(self, s):
        window_values = [2 * i - 1 for i in range(1, 2 ** (self.w - 2) + 1)]
        z = self.zeros()[s]
        scalar = [self.k0, self.k1][s]
        if len(scalar) > z >= self.w - 1:
            return [-v for v in window_values] + [0] + window_values
        if len(scalar) == z:
            return [0] + window_values
        return [0]


def wnaf_values(w):
    values = [2 * i - 1 for i in range(1, 2 ** (w - 2) + 1)]
    return [-v for v in values] + [0] + values


def init_all_windows(w, negative=False):
    windows = [[0]]
    for i in range(1, 2 ** (w - 2) + 1):
        windows.append([2 * i - 1])
    if negative:
        windows.extend([-v for v in windows if v != 0])
    combinations = []
    for w1 in windows:
        for w2 in windows:
            combinations.append(WNAFGuess(w, w1, w2))
    return combinations


def verify_left_guess(s0, s1, left_guess, zvpparams, i, w):
    # if (s0==0 and s1<=2**w-1) or (s1==0 and s0<=2**w-1):
    # return GuessVerification.UNKNOWN
    if left_guess == 0 or (s0 == 0 and s1 == 0):
        return GuessVerification.UNKNOWN
    P = dcp.solve_glv_multi_dcp_pari(
        [2 * s0, 2 * s1], left_guess, zvpparams.glv, zvpparams.registers
    )
    if P is None:
        return GuessVerification.UNKNOWN
    if msm.interleaving_oracle_positive(P, zvpparams, i, w)[1]:
        return GuessVerification.TRUE
    return GuessVerification.FALSE


def verify_right_guess(s0, s1, right_guess, zvpparams, i, w):
    # if (s0==0 and s1<=2**w-1) or (s1==0 and s0<=2**w-1):
    # return GuessVerification.UNKNOWN
    if right_guess == 0 or (s0 == 0 and s1 == 0):
        return GuessVerification.UNKNOWN
    s1 *= 2
    s0, s1 = s1 - s0, -s0
    P = dcp.solve_glv_multi_dcp_pari(
        [s0, s1], right_guess, zvpparams.glv, zvpparams.registers
    )
    if P is None:
        return GuessVerification.UNKNOWN
    P = ZZ(zvpparams.glv.order_field(zvpparams.glv.lam) ** (-1)) * P
    if msm.interleaving_oracle_positive(P, zvpparams, i, w)[1]:
        return GuessVerification.TRUE
    return GuessVerification.FALSE


def zvp_guess(i, scalar_guesses, zvpparams, s):
    new_scalar_guesses = []
    measured_guesses = []
    for winguess in scalar_guesses:
        wg0, wg1 = winguess.values()
        for g in winguess.next_guess_s(s):
            new_guess = deepcopy(winguess)
            new_guess.add_s(g, s)
            if s == 0:
                oracle_output = verify_left_guess(wg0, wg1, g, zvpparams, i, winguess.w)
            if s == 1:
                oracle_output = verify_right_guess(
                    wg0, wg1, g, zvpparams, i, winguess.w
                )
            if oracle_output == GuessVerification.UNKNOWN:
                new_scalar_guesses.append(new_guess)
                # print("tolerated",new_guess)
                continue
            if oracle_output == GuessVerification.TRUE:
                measured_guesses.append(new_guess)
                # print("accepted",new_guess)
                continue
            # print("rejected",new_guess)
    if not measured_guesses:
        return new_scalar_guesses
    return measured_guesses


def zvp_attack(zvpparams, w):
    """ZVP-GLV attack on MSM"""
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    zvp_target = zvpparams.target_bits
    nguesses = []
    l = zvpparams.secrets.half_bits + 1
    b0, b1 = msm.pad_nafs(msm.wnaf(k0, w), msm.wnaf(k1, w), l)
    scalar_guesses = init_all_windows(w)
    for i in range(l - 2, l - 1 - zvp_target, -1):
        scalar_guesses = zvp_guess(
            i, scalar_guesses, zvpparams, 0
        )  # TODO dont recompute division polynomials in pari
        scalar_guesses = zvp_guess(i, scalar_guesses, zvpparams, 1)

        assert WNAFGuess(w, b0[i:], b1[i:]) in scalar_guesses, (
            b0,
            b1,
            i,
            l,
            zvpparams.secrets.to_dict(),
            scalar_guesses,
            WNAFGuess(w, b0[i:], b1[i:]),
        )
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


def zvp_glv_interleaving(params, w):
    result = {}
    time_zvp = time.time()
    try:
        kipairs, nguesses = zvp_attack(params, w)
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


def zvp_glv_interleaving_3(params):
    return zvp_glv_interleaving(params, 3)


def zvp_glv_interleaving_4(params):
    return zvp_glv_interleaving(params, 4)


def zvp_glv_interleaving_5(params):
    return zvp_glv_interleaving(params, 5)
