import dcp
import msm
import time
from sage.all import ZZ, RR, ceil, log


def agreement(scalar_guesses, bits_zvp):
    """ Computes the intersection of our guesses
    currently not using"""
    l = bits_zvp + 1
    result1, result2 = scalar_guesses[0]
    result1bin, result2bin = msm.twodim_recoding(
        msm.to_bin(result1, l), msm.to_bin(result2, l)
    )
    for i in range(l - 2, -1, -1):
        r1, r2 = result1bin[i], result2bin[i]
        agree = True
        for g1, g2 in scalar_guesses[1:]:
            g1bin, g2bin = msm.twodim_recoding(msm.to_bin(g1, l), msm.to_bin(g2, l))
            if g1bin[i] != r1 or g2bin[i] != r2:
                agree = False
                break
        if not agree:
            return result1bin[(i + 1) :], result2bin[(i + 1) :]
    return result1bin, result2bin

def zvp_guess(i, scalar_guesses, params):
    """ZVP-GLV attack on the i-th upper bits of k0,k1"""
    new_scalar_guesses = []

    for scalar_guess_tuple in scalar_guesses:
        scalar_guess = 2 * msm.from_bin(scalar_guess_tuple[0]), 2 * msm.from_bin(
            scalar_guess_tuple[1]
        )
        for guess1 in [0, 1]:
            P = dcp.solve_multi_dcp_pari(scalar_guess, guess1, params)
            if P is None:
                for guess in [[i, i * guess1] for i in [1, -1]]:
                    new_guess = [[guess[j]] + scalar_guess_tuple[j] for j in [0, 1]]
                    if (
                        msm.from_bin(new_guess[0]) >= 0
                        and msm.from_bin(new_guess[1]) >= 0
                    ):
                        new_scalar_guesses.append(new_guess)
                continue
            found = msm.scalar_mul_oracle_positive(
                params["k0"],
                params["k1"],
                P,
                params["lam"],
                params["beta"],
                params["order"],
                i,
            )[1]
            if found:
                return [
                    [[guess[j]] + scalar_guess_tuple[j] for j in [0, 1]]
                    for guess in [[i, i * guess1] for i in [1, -1]]
                ]
    return new_scalar_guesses

def abort_zvp(nguesses, scalar_guesses,params):
    """Abort when the number of scalar_guesses surpasses a threshold (experimentally computed)"""
    bits = ZZ(params["order"]).nbits()
    if bits>100:
        return len(scalar_guesses)>100
    bsgs_timing = {"40": {"2": 1.4586899757385254, "3": 0.6649922847747802, "4": 0.3366295337677002}, "44": {"3": 3.2830397129058837, "4": 1.5843918323516846, "5": 0.7473727703094483}, "48": {"3": 16.571412801742554, "4": 7.442282104492188, "5": 3.481117105484009, "6": 1.4398576736450195}, "52": {"4": 37.54068894386292, "5": 16.539784955978394, "6": 7.099780702590943, "7": 3.4918932914733887}, "56": {"4": 171.46612825393677, "5": 75.4439037322998, "6": 38.669978284835814, "7": 17.450310134887694}, "60": {"5": 390.7318987369537, "6": 168.30803108215332, "7": 83.51950578689575, "8": 37.850191020965575}, "64": {"6": 811.160388469696, "7": 379.64587874412535, "8": 170.56770639419557},"68": {"7": 1112.5570253133774, "8": 479.5018218755722}, "72": {"8": 2355.917799592018, "9": 999.105131149292}}
    msm_dcp_timing = {"40": [0.027345800399780275, 0.028705286979675292, 0.04023177623748779, 0.05399587154388428, 0.12885124683380128], "44": [0.02752363681793213, 0.03117198944091797, 0.03439748287200928, 0.05998327732086182, 0.14451353549957274], "48": [0.02688140869140625, 0.0311568021774292, 0.04145023822784424, 0.07178902626037598, 0.17770805358886718, 0.7723205089569092], "52": [0.030724048614501953, 0.0337996244430542, 0.038432073593139646, 0.06469981670379639, 0.2147299528121948, 0.9621474742889404, 4.545306944847107], "56": [0.028928422927856447, 0.029544734954833986, 0.04648702144622803, 0.07295942306518555, 0.2585240364074707, 1.2503008842468262, 6.229953813552856], "60": [0.02784292697906494, 0.030940818786621093, 0.04206216335296631, 0.07759618759155273, 0.23543293476104737, 1.0590540170669556, 5.0807860612869264, 26.0371036529541], "64": [0.030629730224609374, 0.03234548568725586, 0.053268265724182126, 0.07961657047271728, 0.2805695772171021, 1.262495183944702, 6.156854581832886, 30.771571969985963],"68": [0.03642606735229492, 0.04044771194458008, 0.062442898750305176, 0.11857373714447021, 0.46853954792022706, 2.189353275299072, 10.584056496620178, 51.81921377182007], "72": [0.03867483139038086, 0.03953917026519775, 0.06321074962615966, 0.14159255027770995, 0.5586904764175415, 2.488979959487915, 12.48653199672699, 60.647889852523804, 294.82237355709077]}
    zvp_target = len(nguesses)
    alpha = 0.39
    try:
        bsgs_time = 0.5*len(scalar_guesses)*bsgs_timing[str(bits)][str(zvp_target)]
        next_step_time = msm_dcp_timing[str(bits)][zvp_target]*len(scalar_guesses)*2*((1-alpha)*0.5+alpha)
        expected_guesses = 4*alpha**2*len(scalar_guesses)+2-2*alpha**2
        next_bsgs_time = expected_guesses*bsgs_timing[str(bits)][str(zvp_target+1)]
    except (KeyError, IndexError):
        return False 
    return next_bsgs_time+next_step_time>bsgs_time


def zvp_attack(params, verbose=True):
    """ZVP-GLV attack on MSM"""
    k0, k1 = params["k0"], params["k1"]
    zvp_target = params["zvp_target"]
    nguesses = []
    l = ceil(log(params["order"], 2) / 2) + 1
    b0, b1 = msm.twodim_recoding(msm.to_bin(k0, l), msm.to_bin(k1, l))
    scalar_guesses = [[[ZZ(1)], [ZZ(0)]], [[ZZ(1)], [ZZ(1)]]]
    for i in range(l - 2, l - 2 - zvp_target, -1):
        scalar_guesses = zvp_guess(i, scalar_guesses, params)
        assert [b0[i:], b1[i:]] in scalar_guesses
        if verbose:
            print("Number of guesses: ", len(scalar_guesses))
        nguesses.append(len(scalar_guesses))
        if i>l - 1 - zvp_target and abort_zvp(nguesses,scalar_guesses,params):
            zvp_target = len(nguesses)
            print("aborting")
            break
    kipairs = []
    for s1, s2 in scalar_guesses:
        assert len(s1) == zvp_target + 1
        assert len(s2) == zvp_target + 1
        k1bin, k2bins = msm.two_dim_recoding_inverse_upper_bits(l, s1, s2)
        assert len(k1bin) == zvp_target, f"{k1bin},{zvp_target}"
        for k2bin in k2bins:
            assert len(k2bin) == zvp_target, f"{k2bin},{zvp_target}"
            kipairs.append([msm.from_bin(k1bin), msm.from_bin(k2bin)])
    kipairs = list(set([(a, b) for a, b in kipairs]))
    return [[int(k0), int(k1)] for k0, k1 in kipairs], [int(g) for g in nguesses]


def zvp_attack_timing(params):
    """Measure duration of DCP solver"""
    k0, k1 = params["k0"], params["k1"]
    zvp_target = params["zvp_target"]
    timing = []
    l = ceil(log(params["order"], 2) / 2) + 1
    b0, b1 = msm.twodim_recoding(msm.to_bin(k0, l), msm.to_bin(k1, l))
    for i in range(l - 2, l - 2 - zvp_target, -1):
        t = time.time()
        b0u,b1u =b0[i+1:],b1[i+1:]
        scalar_guess = 2 * msm.from_bin(b0u), 2 * msm.from_bin(b1u)
        dcp.solve_multi_dcp_pari(scalar_guess, 1, params)
        timing.append(time.time()-t)
    return timing 


def zvp_glv_msm(params, verbose=True):
    """ZVP-GLV attack on MSM"""
    result = {"params": params}
    time_zvp = time.time()
    kipairs, result["nguesses"] = zvp_attack(params, verbose)
    result["zvp_target_2"] = len(result["nguesses"]) 
    result["zvp_recovered"] = float(
        max(ZZ(0), RR(2 * len(result["nguesses"]) - log(len(kipairs), 2)))
    )
    result["time_zvp"] = float(time.time() - time_zvp)
    result["kipairs"] = kipairs
    return result


