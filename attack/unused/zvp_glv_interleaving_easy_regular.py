import dcp
import msm
import utils
import time
from sage.all import ZZ, RR, ceil, log
from copy import deepcopy
from utils import GuessVerification


class RWNAFGuess:
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


    def values(self):
        return msm.from_regular_wnaf(self.k0, self.w), msm.from_regular_wnaf(self.k1, self.w)

    def is_positive(self):
        k0, k1 = self.values()
        return k0 >= 0 and k1 >= 0

    def __eq__(self, other):
        return self.k0 == other.k0 and self.k1 == other.k1

    def sizes(self):
        return len(self.k0), len(self.k1)

    def __repr__(self):
        return f"RWNAF(k0={self.k0}, k1={self.k1})"

    def next_guess_s(self, s):
        window_values = [2 * i - 1 for i in range(1, 2 ** (self.w - 1) + 1)]
        scalar = [self.k0, self.k1][s]
        if scalar>0:
            return [-v for v in window_values] + window_values
        return window_values


def positive_rwnaf_values(w):
    values = [2 * i - 1 for i in range(1, 2 ** (w - 1) + 1)]
    return values

def all_rwnaf_values(w):
    values = positive_rwnaf_values(w)
    return [-v for v in values]+values

def zvp_guess(guess, zvpparams, w):
    w0, w1 = guess
    P = dcp.solve_multi_dcp(w1, w0, zvpparams.glv, zvpparams.registers, mul_lam=True)
    if P is None:
        return []
    result = []
    oracle_output = msm.regular_interleaving_easy_oracle_positive(P, zvpparams, w)[1]
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
    l = ceil(zvpparams.secrets.half_bits/w)
    b0, b1 = msm.regular_wnaf(k0, w), msm.regular_wnaf(k1, w)
    rwnaf_values = all_rwnaf_values(w)
    guesses = [[u, v] for u in rwnaf_values for v in rwnaf_values]
    position_vectors = []
    for guess in guesses:
        pos_vec = zvp_guess(guess, zvpparams, w)
        if pos_vec:
            position_vectors.append(pos_vec)

    result_vec = []
    for i in range(l):
        possible_windows = set([tuple(vec[i]) for vec in position_vectors if vec[i]!=None])
        if not possible_windows:
            result_vec.append([None]) 
        else:
            assert (b0[i],b1[i]) in possible_windows, (possible_windows, b0[i],b1[i])
            result_vec.append(list(possible_windows))
    return result_vec


def zvp_glv_interleaving_easy_regular(params, w):
    result = {}
    time_zvp = time.time()
    try:
        result_vec = zvp_attack(params, w)
    except Exception as e:
        params.error_log(e)
        raise e
    result["nguesses"] = 0

    # TODO this might not be precise
    l = ceil(params.secrets.half_bits/w)
    pos_num = log(len(positive_rwnaf_values(w)),2)
    all_num = log(len(all_rwnaf_values(w)),2)
    poss = 0
    for i,v in enumerate(result_vec):
        if not v:
            poss+=log(len(v),2)
        else:
            if i==0:
                poss+=pos_num
            else:
                poss+=all_num
    result["recovered"] = float(max(ZZ(0), RR(2 * params.secrets.half_bits - poss)))


    result["time_zvp"] = float(time.time() - time_zvp)
    result["scalars"] = result_vec
    return result


def zvp_glv_interleaving_easy_regular_1(params):
    return zvp_glv_interleaving_easy_regular(params, 1)


def zvp_glv_interleaving_easy_regular_2(params):
    return zvp_glv_interleaving_easy_regular(params, 2)


def zvp_glv_interleaving_easy_regular_3(params):
    return zvp_glv_interleaving_easy_regular(params, 3)
