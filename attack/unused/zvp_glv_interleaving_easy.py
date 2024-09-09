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


def positive_wnaf_values(w):
    values = [2 * i - 1 for i in range(1, 2 ** (w - 2) + 1)]
    return values


def zvp_guess(guess, zvpparams, w):
    w0, w1 = guess
    P = dcp.solve_multi_dcp(w1, w0, zvpparams.glv, zvpparams.registers, mul_lam=True)
    if P is None:
        return []
    result = []
    oracle_output = msm.interleaving_easy_oracle_positive(P, zvpparams, w)[1]
    # assert len(oracle_output)==zvpparams.secrets.half_bits+1
    for v in oracle_output:
        if v:
            result.append(guess)
        else:
            result.append(None)
    # assert len(result)==zvpparams.secrets.half_bits+1
    return result


def zvp_attack(zvpparams, w):
    """ZVP-GLV attack on MSM"""
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    nguesses = []
    l = zvpparams.secrets.half_bits + 1
    b0, b1 = msm.pad_nafs(msm.wnaf(k0, w), msm.wnaf(k1, w), l)
    wnaf_values = positive_wnaf_values(w)
    guesses = [[u, v] for u in wnaf_values for v in wnaf_values]
    position_vectors = []
    for guess in guesses:
        pos_vec = zvp_guess(guess, zvpparams, w)
        if pos_vec:
            position_vectors.append(pos_vec)

    result_vec = [0] * l
    for pos_vec in position_vectors:
        if not pos_vec:
            continue
        for i in range(l):
            if pos_vec[i]:
                result_vec[i] = pos_vec[i]
    ct0, ct1 = 0, 0
    c = 0
    for i in range(l):
        ri = result_vec[i]
        ct0 += b0[i] != 0
        ct1 += b1[i] != 0
        if ri:
            c += 1
            assert (ri[0] == b0[i] or ri[0] == -b0[i]) and (
                ri[1] == b1[i] or ri[1] == -b1[i]
            )
    print(float(c / ct0), float(c / ct1))
    return result_vec


def zvp_glv_interleaving_easy(params, w):
    result = {}
    time_zvp = time.time()
    try:
        kipairs = zvp_attack(params, w)
    except Exception as e:
        params.error_log(e)
        raise e
    result["nguesses"] = 0

    # TODO
    # result["recovered"] = float(max(ZZ(0), RR(2 * len(nguesses) - log(len(kipairs), 2))))
    result["time_zvp"] = float(time.time() - time_zvp)
    result["scalars"] = kipairs
    return result


def zvp_glv_interleaving_easy_3(params):
    return zvp_glv_interleaving_easy(params, 3)


def zvp_glv_interleaving_easy_4(params):
    return zvp_glv_interleaving_easy(params, 4)


def zvp_glv_interleaving_easy_5(params):
    return zvp_glv_interleaving_easy(params, 5)
