import dcp
import time
import utils
from sage.all import ZZ, log, RR


# def ltr_oracle(curve, k, P, i, registers):
#     """LTR multiplier with DCP oracle for f=x1+x2+2"""

#     if k < 0:
#         Q, found = ltr_oracle(-k, P, i)
#         return -Q, found
#     Q = curve(0)
#     found = False
#     for counter, bit in enumerate(bin(k)[2:]):
#         Q *= 2
#         if bit == "1":
#             if counter == i:
#                 found = registers.is_zero(P,Q)
#             Q += P
#     return Q, found


def ltr_oracle_always(curve, k, P, i, registers):
    """LTR-always multiplier with DCP oracle for f=x1+x2+2"""

    if k < 0:
        Q, found = ltr_oracle_always(-k, P, i)
        return -Q, found
    Qs = [curve(0), curve(0)]
    found = False
    for counter, bit in enumerate(bin(k)[2:]):
        Qs[0] = 2 * Qs[0]
        if counter == i:
            found = found = registers.is_zero(P, Qs[0])
        Qs[1 - int(bit)] = Qs[0] + P
    return Qs[0], found


# def zvp_guess(guesses, k, curve, i):
#     """ZVP attack on the i-th upper bit of LTR"""
#     new_guesses = []
#     for guess in guesses:
#         if guess % 2 == 0:
#             new_guesses.append(guess)
#             continue
#         P = dcp.solve_dcp(guess - 1, {"curve":curve})
#         if P is None:
#             new_guesses.append(guess)
#             continue
#         found = ltr_oracle(curve, k, P, i)[1]
#         if found:
#             return [guess]
#         elif guess % 2 == 1 and len(guesses) == 2:
#             output = guesses[0] if guesses[0] != guess else guesses[1]
#             return [output]
#     return new_guesses


def zvp_guess_always(guesses, zvpparams, i):
    """ZVP attack on the i-th upper bit of LTR-always"""
    new_guesses = []
    curve = zvpparams.glv.curve
    k = zvpparams.secrets.k

    for guess in guesses:
        P = None
        if guess != 0:
            P = dcp.solve_dcp(2 * guess, zvpparams.glv, zvpparams.registers)
        if P is None:
            new_guesses.append(guess)
            continue
        found = ltr_oracle_always(curve, k, P, i + 1, zvpparams.registers)[1]
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


def zvp_attack(zvpparams):
    """ZVP attack on k during LTR-always
    (change below to zvp_guess for simple LTR)"""
    scalar_guesses = [ZZ(0), ZZ(1)]
    secrets = zvpparams.secrets
    nguesses = []
    for i in range(zvpparams.target_bits):
        scalar_guesses = zvp_guess_always(scalar_guesses, zvpparams, i)
        nguesses.append(len(scalar_guesses))
        if zvpparams.verbose:
            print(f"Bit:{i}: {[bin(i) for i in scalar_guesses]}")
        assert ZZ(bin(secrets.k)[: (i + 3)]) in scalar_guesses
        new_guesses = []
        if i < zvpparams.target_bits - 1:
            for g in scalar_guesses:
                new_guesses.append(g * 2)
                new_guesses.append(g * 2 + 1)
            scalar_guesses = new_guesses
    return [int(k_computed) for k_computed in scalar_guesses], [
        int(g) for g in nguesses
    ]


def zvp_ltr(zvpparams):
    """Classical ZVP attack"""
    result = {}
    time_zvp = time.time()
    result["scalars"], result["nguesses"] = zvp_attack(zvpparams)
    result["recovered"] = float(
        max(ZZ(0), RR(len(result["nguesses"]) - log(len(result["scalars"]), 2)))
    )
    result["time_zvp"] = float(time.time() - time_zvp)
    return result
