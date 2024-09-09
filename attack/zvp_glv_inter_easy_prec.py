import dcp
import msm
import utils
from sage.all import ZZ, ceil, log
from copy import deepcopy
import json


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
        return msm.from_regular_wnaf(self.k0, self.w), msm.from_regular_wnaf(
            self.k1, self.w
        )

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
        if scalar > 0:
            return [-v for v in window_values] + window_values
        return window_values


def positive_rwnaf_values(w):
    values = [2 * i - 1 for i in range(1, 2 ** (w - 1) + 1)]
    return values


def all_rwnaf_values(w):
    values = positive_rwnaf_values(w)
    return [-v for v in values] + values


def dcp_points_remapping(zvpparams, w):
    zvpparams = deepcopy(zvpparams)
    assert len(zvpparams.registers.polynomials_tuples) == 1
    zvpparams.target_bits = w
    zvpparams.attack = dcp.interleaving_dcp_experiment
    results = utils.load_results(zvpparams, utils.register_submatch, distinct=True)
    all_points = set()
    for result in results:
        if result["result"]["point"]:
            point = zvpparams.glv.curve(*result["result"]["point"])
            all_points.add(point)

    # remapping
    all_pos_scalars = all_rwnaf_values(w)
    point_list = []
    for point in all_points:
        scalars = []
        for w0 in all_pos_scalars:
            for w1 in all_pos_scalars:
                P, Q = w0 * point, w1 * zvpparams.glv.lam * point
                if zvpparams.registers.is_zero(P, Q):
                    scalars.append((int(w0), int(w1)))
        point_list.append([(int(point[0]), int(point[1])), scalars])
    # assert counter_found==len(zvpparams.registers.polynomials_tuples), (counter_found,len(zvpparams.registers.polynomials_tuples),w0,w1)
    to_save = {"register": zvpparams.registers.to_strings(), "points": point_list}

    try:
        with open(f"results/interleaving_secp256k1_remapped_{w}.json") as f:
            data = json.load(f)
    except FileNotFoundError:
        data = []
    data.append(to_save)
    with open(f"results/interleaving_secp256k1_remapped_{w}.json", "w") as f:
        json.dump(data, f)


def zvp_prec(zvpparams, w, dcppoints):
    l = ceil(zvpparams.secrets.half_bits / w)
    awv = all_rwnaf_values(w)
    all_rwnaf_combinations = [(u, v) for u in awv for v in awv]

    result_vec = [set(all_rwnaf_combinations) for _ in range(l)]
    for point, scalars in dcppoints.items():
        oracle_output = msm.regular_interleaving_easy_oracle_positive(
            point, zvpparams, w
        )[1]
        for i, v in enumerate(oracle_output):
            if v:
                result_vec[i] = result_vec[i].intersection(scalars)
            else:
                result_vec[i] = result_vec[i].difference(scalars)

    return result_vec


def zvp_prec_sim(zvpparams, w, dcppoints):
    l = ceil(zvpparams.secrets.half_bits / w)
    awv = all_rwnaf_values(w)
    all_rwnaf_combinations = [(u, v) for u in awv for v in awv]
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    b0, b1 = msm.regular_wnaf(k0, w), msm.regular_wnaf(k1, w)
    result_vec = []
    for w0, w1 in zip(b0, b1):
        possibilities = set(all_rwnaf_combinations)
        for _, scalars in dcppoints.items():
            for reg_scalars in scalars:
                if (w0, w1) in reg_scalars:
                    possibilities = possibilities.intersection(set(reg_scalars))
                else:
                    possibilities = possibilities.difference(set(reg_scalars))
        result_vec.append(possibilities)
    return result_vec


def zvp_attack(zvpparams, w):
    """ZVP-GLV attack on MSM"""
    with open(f"results/interleaving_secp256k1_remapped_{w}.json") as f:
        dcppoints_all = json.load(f)
    dcppoints = {}
    for register_points in dcppoints_all:
        if utils.register_submatch(zvpparams.registers, register_points["register"]):
            for [[x, y], scalars] in register_points["points"]:
                point = zvpparams.glv.curve(x, y)
                if not point in dcppoints:
                    dcppoints[point] = []
                dcppoints[point].append([tuple(s) for s in scalars])
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    b0, b1 = msm.regular_wnaf(k0, w), msm.regular_wnaf(k1, w)
    result_vec = zvp_prec_sim(zvpparams, w, dcppoints)
    for i, v in enumerate(result_vec):
        assert (b0[i], b1[i]) in v
    return result_vec


def zvp_glv_interleaving_easy_regular(params, w):
    result = {}
    try:
        result_vec = zvp_attack(params, w)
    except Exception as e:
        params.error_log(e)
        raise e
    pos = sum(map(lambda x: log(len(x), 2), result_vec))
    result["remaining"] = float(pos)
    return result


def zvp_glv_interleaving_easy_regular_3(params):
    return zvp_glv_interleaving_easy_regular(params, 3)


def zvp_glv_interleaving_easy_regular_4(params):
    return zvp_glv_interleaving_easy_regular(params, 4)


def zvp_glv_interleaving_easy_regular_5(params):
    return zvp_glv_interleaving_easy_regular(params, 5)
