import dcp
import time
from sage.all import PolynomialRing, ZZ, log, RR


def ltr_oracle(curve, k, P, i):
    """LTR multiplier with DCP oracle for f=x1+x2+2"""
    x1, x2, _, _ = PolynomialRing(curve.base_field(), ["x1", "x2", "y1", "y2"]).gens()
    f = x1 + x2 + 2

    if k < 0:
        Q, found = ltr_oracle(-k, P, i)
        return -Q, found
    Q = curve(0)
    found = False
    for counter, bit in enumerate(bin(k)[2:]):
        Q *= 2
        if bit == "1":
            if counter == i:
                found = f(P[0], Q[0], P[1], Q[1]) == 0
            Q += P
    return Q, found


def ltr_oracle_always(curve, k, P, i):
    """LTR-always multiplier with DCP oracle for f=x1+x2+2"""
    x1, x2, _, _ = PolynomialRing(curve.base_field(), ["x1", "x2", "y1", "y2"]).gens()
    f = x1 + x2 + 2

    if k < 0:
        Q, found = ltr_oracle_always(-k, P, i)
        return -Q, found
    Qs = [curve(0), curve(0)]
    found = False
    for counter, bit in enumerate(bin(k)[2:]):
        Qs[0] = 2 * Qs[0]
        if counter == i:
            found = f(P[0], Qs[0][0], P[1], Qs[0][1]) == 0
        Qs[1 - int(bit)] = Qs[0] + P
    return Qs[0], found


def zvp_guess(guesses, k, curve, i):
    """ZVP attack on the i-th upper bit of LTR"""
    new_guesses = []
    for guess in guesses:
        if guess % 2 == 0:
            new_guesses.append(guess)
            continue
        P = dcp.solve_dcp(guess - 1, curve)
        if P is None:
            new_guesses.append(guess)
            continue
        found = ltr_oracle(curve, k, P, i)[1]
        if found:
            return [guess]
        elif guess % 2 == 1 and len(guesses) == 2:
            output = guesses[0] if guesses[0] != guess else guesses[1]
            return [output]
    return new_guesses


def zvp_guess_always(guesses, k, curve, i):
    """ZVP attack on the i-th upper bit of LTR-always"""
    new_guesses = []
    for guess in guesses:
        P = dcp.solve_dcp(2 * guess, curve)
        if P is None:
            new_guesses.append(guess)
            continue
        found = ltr_oracle_always(curve, k, P, i + 1)[1]
        if found:
            return [guess]
        elif len(guesses) == 2:
            output = guesses[0] if guesses[0] != guess else guesses[1]
            return [output]
    return new_guesses


def agreement(scalar_guesses, bits_zvp):
    """Computes the intersection of our guesses"""
    scalar_strings = [format(s, f"0{bits_zvp}b") for s in scalar_guesses]
    s1 = min(scalar_strings)
    s2 = max(scalar_strings)
    for i, c in enumerate(s1):
        if c != s2[i]:
            if s1[:i] == "":
                return -1, i
            return ZZ(s1[:i], 2), i
    return ZZ(s1, 2), len(s1)


def zvp_attack(bits_zvp, k, curve, verbose=True, intersection=True):
    """ZVP attack on k during LTR-always
    (change below to zvp_guess for simple LTR)"""
    scalar_guesses = [ZZ(0), ZZ(1)]
    nguesses = []
    for i in range(bits_zvp):
        scalar_guesses = zvp_guess_always(scalar_guesses, k, curve, i)
        nguesses.append(len(scalar_guesses))
        if verbose:
            print(f"Bit:{i}: {[bin(i) for i in scalar_guesses]}")
        assert ZZ(bin(k)[: (i + 3)]) in scalar_guesses
        new_guesses = []
        for g in scalar_guesses:
            new_guesses.append(g * 2)
            new_guesses.append(g * 2 + 1)
        scalar_guesses = new_guesses
    if not intersection:
        return [int(k_computed) for k_computed in scalar_guesses], [
            int(g) for g in nguesses
        ]
    k_computed, recovered_bits = agreement(scalar_guesses, bits_zvp)
    if k_computed == -1:
        if verbose:
            print("No bits found")
        return int(0), 0, nguesses
    expected = bin(k)[: bits_zvp + 2]
    assert expected.startswith(bin(k_computed)), f"{expected,bin(k_computed)}"
    if verbose:
        print(
            f"Found: {bin(k_computed)} ({recovered_bits} bits) out of {expected} ({bits_zvp} bits)"
        )
    return int(k_computed), int(recovered_bits), [int(g) for g in nguesses]


def zvp_sm(params, verbose=True):
    """Classical ZVP attack"""
    zvp_target = params["zvp_target"]
    curve = params["curve"]
    result = {"params": params}
    time_zvp = time.time()
    result["k_upper"], result["zvp_recovered"], result["nguesses"] = zvp_attack(
        zvp_target, params["k"], curve, verbose
    )
    result["time_zvp"] = float(time.time() - time_zvp)
    return result


def zvp_glv_sm_set(params, verbose=True):
    """ZVP-GLV attack which returns all possible scalars
    and not their common prefix (as opposed to zvp_glv_sm)"""
    zvp_target = params["zvp_target"]
    curve = params["curve"]
    result = {"params": params}
    time_zvp = time.time()
    k0_uppers, result["nguesses_0"] = zvp_attack(
        zvp_target, params["k0"], curve, verbose, intersection=False
    )
    result["zvp_recovered_0"] = float(
        max(ZZ(0), RR(zvp_target - log(len(k0_uppers), 2)))
    )
    k1_uppers, result["nguesses_1"] = zvp_attack(
        zvp_target, params["k1"], curve, verbose, intersection=False
    )
    result["zvp_recovered_1"] = float(
        max(ZZ(0), RR(zvp_target - log(len(k1_uppers), 2)))
    )
    result["kipairs"] = [[k0, k1] for k0 in k0_uppers for k1 in k1_uppers]
    result["time_zvp"] = float(time.time() - time_zvp)
    return result


def zvp_glv_sm(params, verbose=True):
    """ZVP-GLV attack"""
    zvp_target = params["zvp_target"]
    curve = params["curve"]
    result = {"params": params}
    time_zvp = time.time()
    result["k0_upper"], result["zvp_recovered_0"], result["nguesses_0"] = zvp_attack(
        zvp_target, params["k0"], curve, verbose
    )
    result["k1_upper"], result["zvp_recovered_1"], result["nguesses_1"] = zvp_attack(
        zvp_target, params["k1"], curve, verbose
    )
    result["time_zvp"] = float(time.time() - time_zvp)
    return result
