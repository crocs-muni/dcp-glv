from sage.all import ceil, log, floor, ZZ
import msm
import bsgssolver as bsgs
import time
import zvp_glv_sm, zvp_glv_msm


def twostep_sm(params, verbose=True):
    """2-step attack on SM"""
    curve, lam, order = params["curve"], params["lam"], params["order"]
    k, k0, k1 = params["k"], params["k0"], params["k1"]
    result = zvp_glv_sm.zvp_glv_sm(params, verbose)
    assert k == (k0 + k1 * lam) % order
    lower_bits0 = k0.nbits() - result["zvp_recovered_0"]
    lower_bits1 = k1.nbits() - result["zvp_recovered_1"]
    k0_upper = result["k0_upper"]
    k1_upper = result["k1_upper"]

    # Baby step giant step on lower bits
    P = curve.random_point()
    assert P != curve(0)
    Q = zvp_glv_sm.ltr_oracle(curve, k, P, 0)[0]
    assert Q == k0 * P + k1 * lam * P
    assert (
        k0 // (2**lower_bits0) == k0_upper
    ), f"{bin(k0//(2**lower_bits0))},{bin(k0_upper)}"
    assert k1 // (2**lower_bits1) == k1_upper
    Q2 = Q - 2**lower_bits0 * k0_upper * P - 2**lower_bits1 * k1_upper * lam * P
    ek0_lower, ek1_lower = k0 % (2**lower_bits0), k1 % (2**lower_bits1)
    assert Q2 == ek0_lower * P + ek1_lower * lam * P
    time_bsgs = time.time()
    k0_lower, k1_lower = bsgs.bsgs(
        dimension=2,
        points=[P, Q2],
        endomorphisms=[lam],
        upperbounds=[lower_bits0, lower_bits1],
    )

    # derive the private key
    recovered_k0 = k0_lower + k0_upper * 2**lower_bits0
    recovered_k1 = k1_lower + k1_upper * 2**lower_bits1
    recovered_k = (recovered_k0 + recovered_k1 * lam) % order
    assert k == recovered_k
    if verbose:
        print("Private key:", k, k0, k1)
        print(f"Computed: {recovered_k}, found: {recovered_k==k}")
    result["time_bsgs"] = float(time.time() - time_bsgs)
    return result


def twostep_sm_set(params, verbose=True):
    """2-step attack on SM"""
    curve, lam, order = params["curve"], params["lam"], params["order"]
    k, k0, k1 = params["k"], params["k0"], params["k1"]
    result = zvp_glv_sm.zvp_glv_sm_set(params, verbose)
    assert k == (k0 + k1 * lam) % order
    # Baby step giant step on lower bits
    P = curve.random_point()
    assert P != curve(0)
    Q = zvp_glv_sm.ltr_oracle(curve, k, P, 0)[0]
    assert Q == k0 * P + k1 * lam * P
    lower_bits0 = k0.nbits() - params["zvp_target"]-1
    lower_bits1 = k1.nbits() - params["zvp_target"]-1
    time_bsgs = time.time()
    for k0_upper,k1_upper in result["kipairs"]:
        Q2 = Q - 2**lower_bits0 * k0_upper * P - 2**lower_bits1 * k1_upper * lam * P
        try:
            k0_lower, k1_lower = bsgs.bsgs(
                dimension=2,
                points=[P, Q2],
                endomorphisms=[lam],
                upperbounds=[lower_bits0, lower_bits1],
            )
        except bsgs.NoSolution:
            continue

        # derive the private key
        recovered_k0 = k0_lower + k0_upper * 2**lower_bits0
        recovered_k1 = k1_lower + k1_upper * 2**lower_bits1
        recovered_k = (recovered_k0 + recovered_k1 * lam) % order
        assert k == recovered_k
        if verbose:
            print("Private key:", k, k0, k1)
            print(f"Computed: {recovered_k}, found: {recovered_k==k}")
        result["time_bsgs"] = float(time.time() - time_bsgs)
        return result
    raise Exception



def twostep_msm(params, verbose=True):
    """2-step attack on MSM"""
    curve, lam, order, beta = (
        params["curve"],
        params["lam"],
        params["order"],
        params["beta"],
    )
    k, k0, k1 = params["k"], params["k0"], params["k1"]
    result = zvp_glv_msm.zvp_glv_msm(params, verbose)
    zvp_target = result["zvp_target_2"]
    kipairs = result["kipairs"]
    zvp_recovered = result["zvp_recovered"]
    l = ceil(log(order, 2) / 2) + 1
    if verbose:
        print(f"Recovered bits: {floor(zvp_recovered)}")
    if verbose:
        print(
            f"{bin(k0)[2+zvp_target:]},{bin(k1)[2+zvp_target:]} remains for bsgs with {len(kipairs)} options"
        )

    P = curve.random_point()
    assert P != curve(0)
    Q = msm.scalar_mul_oracle_positive(k0, k1, P, lam, beta, order, 0)[0]
    assert Q == k0 * P + k1 * lam * P
    time_bsgs = time.time()
    for k0_upper, k1_upper in kipairs:
        k1_upper = ZZ(k1_upper)
        lower_bits0 = l - 1 - zvp_target + 1
        lower_bits1 = l - 1 - zvp_target + 1

        # Baby step giant step on lower bits
        Q2 = Q - 2**lower_bits0 * k0_upper * P - 2**lower_bits1 * k1_upper * lam * P
        try:
            k0_lower, k1_lower = bsgs.bsgs(
                dimension=2,
                points=[P, Q2],
                endomorphisms=[lam],
                upperbounds=[lower_bits0, lower_bits1],
            )
        except bsgs.NoSolution:
            continue

        # derive the private key
        if verbose:
            print(bin(k0_lower), bin(k1_upper))
        recovered_k0 = k0_lower + k0_upper * 2**lower_bits0
        recovered_k1 = k1_lower + k1_upper * 2**lower_bits1
        recovered_k = (recovered_k0 + recovered_k1 * lam) % order
        if verbose:
            print("Private key:", k, k0, k1)
            print(f"Computed: {recovered_k}, found: {recovered_k==k}")
        result["time_bsgs"] = float(time.time() - time_bsgs)
        return result
    raise Exception
