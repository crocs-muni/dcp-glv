import dcp
import time
from sage.all import ZZ, log, RR
import msm
from copy import deepcopy
from utils import GuessVerification


class SignedGuess:
    def __init__(self, k=None):
        self.k = [ZZ(b) for b in k] if k is not None else []

    def add(self, b):
        self.k = [ZZ(b)] + self.k

    def value(self):
        return msm.from_regular_wnaf(self.k, 1)

    def is_positive(self):
        k = self.value()
        return k >= 0

    def __eq__(self, other):
        return self.k == other.k

    def size(self):
        return len(self.k)

    def __repr__(self):
        return f"Signed(k={self.k})"

    def next_guess(self):
        return [1, -1]


def to_signed(k):
    if k % 2 == 0:
        return msm.regular_wnaf(k - 1, 1), 1
    return msm.regular_wnaf(k, 1), 0


def signed_ltr_oracle(P, zvpparams, i):
    curve = zvpparams.glv.curve
    k = zvpparams.secrets.k

    kbin, ksign = to_signed(k)

    Q = curve(0)
    found = False
    counter = 0
    for b in reversed(kbin):
        Q = 2 * Q
        if counter == i:
            found = found = zvpparams.registers.is_zero(b * P, Q)
        Q = Q + b * P
        counter += 1
    Q += ksign * P
    return Q, found


def verify_guess(s, g, zvpparams, i):
    dcp_scalars = 2 * s, g
    P = dcp.solve_multi_dcp(*dcp_scalars, zvpparams.glv, zvpparams.registers)
    if P is None:
        return GuessVerification.UNKNOWN
    if signed_ltr_oracle(P, zvpparams, i)[1]:
        return GuessVerification.TRUE
    return GuessVerification.FALSE


def zvp_guess(i, scalar_guesses, zvpparams):
    new_scalar_guesses = []
    measured_guesses = []
    for signedguess in scalar_guesses:
        k = signedguess.value()
        for next_guess in signedguess.next_guess():
            new_guess = deepcopy(signedguess)
            new_guess.add(next_guess)
            oracle_output = verify_guess(k, next_guess, zvpparams, i)
            if oracle_output == GuessVerification.UNKNOWN:
                if new_guess.is_positive():
                    new_scalar_guesses.append(new_guess)
                continue
            if oracle_output == GuessVerification.TRUE:
                if new_guess.is_positive():
                    measured_guesses.append(new_guess)
    if not measured_guesses:
        return new_scalar_guesses
    return measured_guesses


def zvp_attack(zvpparams):
    k = zvpparams.secrets.k
    zvp_target = zvpparams.target_bits
    nguesses = [2]
    scalar_guesses = [SignedGuess([1]), SignedGuess([-1])]
    kbin, _ = to_signed(k)  # assuming that the last bit will not be targetted
    for i in range(1, zvp_target):
        scalar_guesses = zvp_guess(i, scalar_guesses, zvpparams)
        assert SignedGuess(kbin[-i - 1 :]) in scalar_guesses, (
            SignedGuess(kbin[-i - 1 :]),
            scalar_guesses,
        )
        if zvpparams.verbose:
            print("Number of guesses: ", len(scalar_guesses))
        nguesses.append(len(scalar_guesses))
    kipairs = []

    for signedguess in scalar_guesses:
        assert signedguess.size() == zvp_target
        kbin = signedguess.k
        kbin = [int(b) for b in kbin]
        kipairs.append(kbin)
    return kipairs, [int(g) for g in nguesses]


def zvp_ltr_signed(zvpparams):
    """Classical ZVP attack"""
    result = {}
    time_zvp = time.time()
    result["scalars"], result["nguesses"] = zvp_attack(zvpparams)
    result["recovered"] = float(
        max(ZZ(0), RR(len(result["nguesses"]) - log(len(result["scalars"]), 2)))
    )
    result["time_zvp"] = float(time.time() - time_zvp)
    return result
