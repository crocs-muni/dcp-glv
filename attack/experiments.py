from sage.all import ZZ, log, sqrt
import zvp_glv_sac
import zvp_glv_inter_easy_prec as inter_easy
import zvp_ltr_signed
from datetime import datetime
import json, os
import numpy
import dcp
import utils
from utils import ZVPparams
from copy import deepcopy

FLOAT_FORMAT = "{0:0.2f}"


def save_results(results, filename):
    if not os.path.exists("results"):
        os.makedirs("results")
    time = datetime.today().strftime("%Y-%m-%d_%H:%M:%S")
    filename = os.path.join("results", f"{filename}_{time}.json")
    with open(
        filename,
        "w",
    ) as f:
        json.dump(results, f)


def dcp_all_secp256k1_interleaving_experiments(zvpparams, w):
    secrets = zvpparams.secrets
    zvpparams.attack = dcp.interleaving_dcp_experiment
    results_list = []
    zvpparams.glv.set_secp256k1()
    wnaf_values = inter_easy.all_rwnaf_values(w)
    for w0 in wnaf_values:
        for w1 in wnaf_values:
            secrets.k = -1
            secrets.k0, secrets.k1 = w0, w1
            zvpparams.do_attack()
            results_list.append(zvpparams.results_to_dict())
    save_results(results_list, f"{zvpparams.filename}_{w}")


def dcp_all_secp256k1_experiments(zvpparams):
    secrets = zvpparams.secrets
    zvpparams.attack = dcp.glv_dcp_experiment
    results_list = []
    lower_bound = 2 ** (zvpparams.target_bits - 1)
    upper_bound = 2 ** (zvpparams.target_bits)
    zvpparams.glv.set_secp256k1()
    for k in range(lower_bound, upper_bound):
        secrets.k = k
        secrets.k0, secrets.k1 = -1, -1
        zvpparams.do_attack()
        results_list.append(zvpparams.results_to_dict())
    save_results(results_list, zvpparams.filename)


def glvdcp_all_secp256k1_experiments(zvpparams):
    secrets = zvpparams.secrets
    zvpparams.attack = dcp.glv_dcp_experiment
    results_list = []
    lower_bound = 2 ** (zvpparams.target_bits - 1)
    upper_bound = 2 ** (zvpparams.target_bits)
    zvpparams.glv.set_secp256k1()
    for k0 in range(lower_bound, upper_bound):
        for k1 in range(lower_bound, upper_bound):
            secrets.k = -1
            secrets.k0, secrets.k1 = k0, k1
            zvpparams.do_attack()
            results_list.append(zvpparams.results_to_dict())
    save_results(results_list, zvpparams.filename)


def zvp_attack_experiments(zvpparams, n_exp, t_bounds):
    for zvp_target in map(ZZ, range(*t_bounds)):
        results_list = []
        zvpparams.target_bits = zvp_target
        for i in range(n_exp):
            if zvpparams.verbose:
                print(f"{i + 1}/{n_exp}")
            zvpparams.generate_secrets()
            zvpparams.do_attack()
            results_list.append(zvpparams.results_to_dict())
        save_results(results_list, zvpparams.filename)


def zvp_attack_experiments_secp256k1(zvpparams, n_exp, t_bounds):
    for zvp_target in map(ZZ, range(*t_bounds)):
        results_list = []
        zvpparams.target_bits = zvp_target
        for i in range(n_exp):
            if zvpparams.verbose:
                print(f"{i + 1}/{n_exp}")
            zvpparams.generate_secp256k1_secrets()
            zvpparams.do_attack()
            results_list.append(zvpparams.results_to_dict())
        save_results(results_list, zvpparams.filename)


def DCP_main():

    """Set up"""
    zvpparams = ZVPparams(bits=256)

    zvpparams.attack = dcp.glv_dcp_experiment
    # zvpparams.attack = dcp.dcp_experiment
    # zvpparams.attack = dcp.interleaving_dcp_experiment
    zvpparams.target_bits = 4

    secp256k1_polynomials = utils.get_secp256k1_polynomials()
    selection = ["f1", "f11", "f6", "f7", "f9"]
    for name in selection:
        zvpparams.registers.add_tuple(*secp256k1_polynomials[name])

    """Run"""
    zvp_attack_experiments_secp256k1(zvpparams, n_exp=10, t_bounds=(4, 5))
    # glvdcp_all_secp256k1_experiments(zvpparams)
    # dcp_all_secp256k1_interleaving_experiments(zvpparams,3)


def interleaving_main():
    zvpparams = ZVPparams(bits=256)
    zvpparams = ZVPparams(bits=256)
    zvpparams.attack = inter_easy.zvp_glv_interleaving_easy_regular_3
    zvpparams.registers = utils.load_efd_secp256k1_registers()["modified:mmadd-2009-bl"]
    zvp_attack_experiments_secp256k1(zvpparams, n_exp=1, t_bounds=(0, 1))


def ZVP_main():

    """Choose bit-length"""
    zvpparams = ZVPparams(bits=256)

    zvpparams.attack = inter_easy.zvp_glv_interleaving_easy_regular_4
    secp256k1_polynomials = utils.get_secp256k1_polynomials()
    selection = ["f1", "f11", "f6", "f7", "f9"]
    for name in selection:
        zvpparams.registers.add_tuple(*secp256k1_polynomials[name])

    zvp_attack_experiments_secp256k1(zvpparams, n_exp=1, t_bounds=(0, 1))


if __name__ == "__main__":

    # DCP_main()
    # ZVP_main()
    # interleaving_main()
    pass
