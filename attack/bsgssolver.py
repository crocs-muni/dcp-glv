"""
Functions implementing various versions of Baby-step Giant-step
Some functions are calling bsgs binary in pari_tools.
"""

import os
from sage.all import ZZ, randint
import time
import glv


class NoSolution(Exception):
    pass


def bsgs1(P, Q, lower_bits):
    """Simple BSGS for 1-dimension"""
    print(f"BABY STEP GIANT STEP on {lower_bits} bits")
    m1 = lower_bits // 2
    m2 = lower_bits - m1
    # Baby steps
    N = 2**m2
    tbl = {}
    R = P.curve()(0)
    for i in range(N):
        tbl[R] = i
        R += P
    # Giant steps
    P2m = 2**m1 * P
    Q3 = Q
    for j in range(N):
        if Q3 in tbl:
            return [2**m1 * j + tbl[Q3]]
        Q3 -= P2m
    raise NoSolution


def bsgs2(P, Q, lam, lower_bits1, lower_bits2):
    """Simple BSGS for 2-dimensional decomposition k = k1+lam*k2
    lower_bits1,lower_bits2 are upper bounds on bitlengths of k1 and k2"""
    N2 = 2**lower_bits2
    N1 = 2**lower_bits1
    lamP = lam * P
    # Baby steps
    tbl = {}
    R = P.curve()(0)
    for i in range(N2):
        tbl[R] = i
        R += lamP
    # Giant steps
    Q3 = Q
    for j in range(N1):
        if Q3 in tbl:
            return [j, tbl[Q3]]
        Q3 -= P
    raise NoSolution


def bsgs4(P, Q, lam, sigma, lower_bits1, lower_bits2, lower_bits3, lower_bits4):
    """Simple BSGS for 4-dimensional decomposition k = k1+lam*k2+sigma*k3+lam*sigma*k4
    lower_bitsi are upper bounds on bitlengths of ki"""
    N2 = 2**lower_bits2
    N1 = 2**lower_bits1
    N3 = 2**lower_bits3
    N4 = 2**lower_bits4
    lamP = lam * P
    sigmaP = sigma * P
    lamsigmaP = lam * sigma * P
    # Baby steps
    tbl = {}
    R = P.curve()(0)
    for i in range(N3):
        S = P.curve()(0)
        for j in range(N4):
            tbl[R + S] = (i, j)
            S += lamsigmaP
        R += sigmaP
    # Giant steps
    Q3 = Q
    R = P.curve()(0)
    for i in range(N1):
        S = P.curve()(0)
        for j in range(N2):
            Q3 = Q - S - R
            if Q3 in tbl:
                return [i, j, *tbl[Q3]]
            S += lamP
        R += P
    raise NoSolution


def bsgs4_balanced(
    P, Q, lam, sigma, lower_bits1, lower_bits2, lower_bits3, lower_bits4
):
    """BSGS for 4-dimensions with balancing, badly written (made obsolete by pari implementation anyway)"""
    if lower_bits1 + lower_bits2 - lower_bits3 - lower_bits4 in [-1, 0, 1]:
        return bsgs4(
            P, Q, lam, sigma, lower_bits1, lower_bits2, lower_bits3, lower_bits4
        )
    diff = lower_bits1 + lower_bits2 - lower_bits3 - lower_bits4
    if diff > 0:

        shared = diff // 2
        print(
            f"Balancing out to complexity of bsgs: {lower_bits1+lower_bits2-shared}, {lower_bits3+shared+lower_bits4}"
        )
        N2 = 2 ** (lower_bits2 - shared)
        N1 = 2**lower_bits1
        N31 = 2**shared
        N32 = 2**lower_bits3
        N4 = 2**lower_bits4
        lamP = lam * P
        lamN2P = lam * N2 * P
        sigmaP = sigma * P
        lamsigmaP = lam * sigma * P

        tbl = {}
        T = P.curve()(0)
        for k in range(N31):
            R = P.curve()(0)
            for i in range(N32):
                S = P.curve()(0)
                for j in range(N4):
                    tbl[R + S + T] = (i, j, k)
                    S += lamsigmaP
                R += sigmaP
            T += lamN2P

        Q3 = Q
        R = P.curve()(0)
        for i in range(N1):
            S = P.curve()(0)
            for j in range(N2):
                Q3 = Q - S - R
                if Q3 in tbl:
                    e1, e2, e3 = tbl[Q3]
                    return [i, N2 * e3 + j, e1, e2]
                S += lamP
            R += P
    else:

        shared = -diff // 2
        print(
            f"Balancing out to complexity of bsgs: {lower_bits1+lower_bits2+shared}, {lower_bits3-shared+lower_bits4}"
        )
        N2 = 2**lower_bits2
        N1 = 2**lower_bits1
        N23 = 2**shared
        N3 = 2 ** (lower_bits3 - shared)
        N4 = 2**lower_bits4
        lamP = lam * P
        sigmaN3P = sigma * N3 * P
        sigmaP = sigma * P
        lamsigmaP = lam * sigma * P

        tbl = {}

        R = P.curve()(0)
        for i in range(N3):
            S = P.curve()(0)
            for j in range(N4):
                tbl[R + S] = (i, j)
                S += lamsigmaP
            R += sigmaP

        Q3 = Q
        T = P.curve()(0)
        for k in range(N23):
            R = P.curve()(0)
            for i in range(N1):
                S = P.curve()(0)
                for j in range(N2):
                    Q3 = Q - S - R - T
                    if Q3 in tbl:
                        e1, e2 = tbl[Q3]
                        return [i, j, N3 * k + e1, e2]
                    S += lamP
                R += P
            T += sigmaN3P
    return NoSolution


def bsgs2_balanced(P, Q, lam, lower_bits1, lower_bits2):
    """BSGS for 2-dimensions with equalizing"""

    diff = lower_bits1 - lower_bits2
    if diff in [0, 1, -1]:
        return bsgs2(P, Q, lam, lower_bits1, lower_bits2)
    if diff < 0:
        big = lower_bits2
        small = lower_bits1
        shared = -diff // 2
    else:
        big = lower_bits1
        small = lower_bits2
        shared = diff // 2

    print(f"Balancing out to complexity of bsgs: {big-shared}, {shared+small}")
    N1 = 2 ** (big - shared)
    lamP = lam * P
    P2m = N1 * P
    new_tbl = {}
    R = P.curve()(0)
    if diff < 0:
        small_point = P
        shared_point = lam * P2m
        big_point = lamP
    else:
        small_point = lamP
        shared_point = P2m
        big_point = P

    # Baby steps, storing i*shared+j*smallpoint
    for i in range(2**shared):
        S = P.curve()(0)
        for j in range(2**small):
            new_tbl[S + R] = (i, j)
            S += small_point
        R += shared_point

    # Giant steps, computing Q-j*big_point
    Q3 = Q
    for j in range(N1):
        if Q3 in new_tbl:
            e1, e2 = new_tbl[Q3]
            if diff < 0:
                return [e2, e1 * N1 + j]
            else:
                return [j + N1 * e1, e2]
        Q3 -= big_point
    raise NoSolution


def bsgs2_pari(P, Q, lam, l1, l2, balanced=True):
    """BSGS for 2-dimensions with balancing. If you want unbalanced version, change bsgs.c"""
    pari_path = "./pari_tools/bsgs_solver"
    solution_path = "./pari_tools/results/bsgs2_solution"
    E = P.curve()
    with open(solution_path, "w") as f:
        pass
    arguments = f"{pari_path} {E.base_field().order()} {E.a4()} {E.a6()} {P[0]} {P[1]} {Q[0]} {Q[1]} {lam} {l1} {l2} {solution_path}"
    os.system(arguments)
    try:
        with open(solution_path, "r") as f:
            result = f.read()
            assert result != ""
    except (FileNotFoundError, AssertionError):
        with open("error.log", "w") as f:
            f.write(arguments)
        raise NoSolution
    return [ZZ(i) for i in result.split("\n")[:-1]]


def bsgs_timing(bits, lbits, N=10):
    """Computes duration of bsgs"""
    params = glv.gen_glv_curve(bits)
    curve, lam, order = params["curve"], params["lam"], params["order"]
    avg_time = 0
    for _ in range(10):
        k0, k1 = ZZ(randint(2 ** (lbits - 1), 2**lbits)), ZZ(
            randint(2 ** (lbits - 1), 2**lbits)
        )
        k = (k0 + k1 * lam) % order

        P = curve.random_point()
        assert P != curve(0)
        Q = k * P
        time_bsgs = time.time()
        k0c, k1c = bsgs(
            dimension=2,
            points=[P, Q],
            endomorphisms=[lam],
            upperbounds=[lbits, lbits],
        )
        assert k == (k0c + k1c * lam) % order
        avg_time += float(time.time() - time_bsgs)
    return avg_time / N


def bsgs(dimension, points, endomorphisms, upperbounds):
    """General function for Baby step giant step on scalar decomposition"""
    if dimension == 2:
        return bsgs2_pari(*points, *endomorphisms, *upperbounds)
    if dimension == 4:
        return bsgs4_balanced(*points, *endomorphisms, *upperbounds)
