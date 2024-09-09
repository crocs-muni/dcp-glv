"""
Implementation of multi-scalar multiplication as designed in https://eprint.iacr.org/2013/158.pdf for dimension 2
- including symbolic version
- including one leaking zero-registers
"""

from glv import glv_decompose
from sage.all import ceil, log


def to_bin(k, l):
    "Given 0b11,3 outputs [1,1,0], i.e. list of l (padded) bits from LSB"
    kbin = list(map(lambda x: int(x), bin(k)[2:]))
    kbin = list(reversed(kbin))
    kbin.extend([0] * (l - len(kbin)))
    return kbin


def from_bin(kbin):
    "Inverse of to_bin"
    return sum([2**i * x for i, x in enumerate(kbin)])


def zero_allign_bin(k, l):
    "Returns bin padded (with 0s) to l bits"
    kbin = bin(k)[2:]
    return "0b" + (l - len(kbin)) * "0" + kbin


def twodim_recoding(kbin, kpbin):
    """
    Algorithm 1 from https://eprint.iacr.org/2013/158.pdf
    Performs recoding on two bit-lists (outputs of to_bin)
    """
    l = len(kbin)
    bbin = [0] * (l - 1) + [1]
    bpbin = [0] * l
    for i in range(l - 1):
        bbin[i] = 2 * kbin[i + 1] - 1
    for i in range(l):
        bpbin[i] = bbin[i] * kpbin[0]
        kp = from_bin(kpbin)
        kp = kp // 2 - bpbin[i] // 2
        kpbin = to_bin(kp, l)
    return bbin, bpbin


def twodim_recoding2(kbin, kpbin):
    """Just a different implementation twodim_recoding2"""
    l = len(kbin)
    bbin = [0] * (l - 1) + [1]
    bpbin = [0] * l
    for i in range(l - 1):
        bbin[i] = 2 * kbin[i + 1] - 1
    carry = 0
    prev_b = 0
    for i in range(l):
        tmp = (prev_b == -1) + kpbin[i] + carry
        bpbin[i] = prev_b = (tmp % 2) * bbin[i]
        carry = tmp > 1
    return bbin, bpbin


def twodim_recoding_inverse(bbin, bpbin):
    """Inverse of twodim_recoding"""
    l = len(bbin)
    kbin, kpbin = [1] + [0] * (l - 1), [0] * l
    for i in range(1, l):
        kbin[i] = (bbin[i - 1] + 1) // 2
    carry = 0
    prev_b = 0
    for i in range(l):
        kpbin[i] = ((bpbin[i] == bbin[i]) + carry + (prev_b == -1)) % 2
        tmp = (prev_b == -1) + kpbin[i] + carry
        carry = tmp > 1
        prev_b = bpbin[i]
    return kbin, kpbin


def two_dim_recoding_inverse_upper_bits(l, bbin, bpbin):
    """Inverse of twodim_recoding where bbin,bpin are only upper bits of l bit integers"""
    kbin = []
    for i, b in enumerate(bbin[:-1]):
        kbin.append((b + 1) // 2)
    if l == len(bbin):
        kbin = [1] + kbin
    xs = [int(b == -1) for b in bpbin[:-1]]
    if l == len(bpbin):
        xs = [0] + xs
    if xs[-1] == 1 or bpbin[-1] == 0:
        candidates = [[0, 0]]
    else:
        candidates = [[0, 1], [1, 0]]
    for row, x in enumerate(reversed(xs[:-1])):
        new_candidates = []
        for cand in candidates:
            carry = cand[1]
            for i, j in [[0, 0], [0, 1], [1, 0], [1, 1]]:
                if ((i + x + j) > 1) != carry:
                    continue
                if ((i + x + j) * bbin[-(row + 2)]) % 2 != (bpbin[-(row + 2)] % 2):
                    continue
                new_candidates.append([cand[0] * 2 + i, j])
        candidates = new_candidates
    kpbins = [to_bin(i, len(xs)) for i, _ in candidates]
    return kbin, kpbins


def twodim_recoding_extended(kbin, s, kpbin, sp):
    """Algorithm 3 from https://eprint.iacr.org/2013/158.pdf"""
    l = len(kbin)
    bbin = [0] * (l - 1) + [s]
    bpbin = [0] * l
    for i in range(l - 1):
        bbin[i] = s * (2 * kbin[i + 1] - 1)
    for i in range(l):
        bpbin[i] = bbin[i] * kpbin[0]
        kp = from_bin(kpbin)
        kp = kp // 2 - sp * bpbin[i] // 2
        kpbin = to_bin(kp, l)
    return bbin, bpbin


def twodim_recoding_extended_inverse(bbin, bpbin):
    """So far works only for positive (no s,sp in the argument), i.e. the same as nonextended"""
    l = len(bbin)
    kbin, kpbin = [1] + [0] * (l - 1), [0] * l
    for i in range(1, l):
        kbin[i] = (bbin[i - 1] + 1) // 2
    carry = 0
    prev_b = 0
    for i in range(l):
        kpbin[i] = ((bpbin[i] == bbin[i]) + carry + (prev_b == -1)) % 2
        tmp = (prev_b == -1) + kpbin[i] + carry
        carry = tmp > 1
        prev_b = bpbin[i]
    return kbin, kpbin


def scalar_mul(k, P, lamb, beta, n):
    """Algorithm 2 from https://eprint.iacr.org/2013/158.pdf for 2 dimensions, i.e. m = 2"""

    l = ceil(log(n, 2) / 2) + 2

    # Precomputation stage ==========
    table = [P, P.curve()(beta**2 * P[0], -P[1])]
    # P = (x,y), P+lamb*P=(beta^2*x,-y)
    # ===============================

    # Recoding stage ================
    k0, k1 = glv_decompose(k, lamb, n)
    even = k0 % 2
    if even == 0:
        k0 = k0 - 1
    s0, s1 = 1, 1
    if k0 < 0:
        s0 = -1
        k0 = -k0
    if k1 < 0:
        s1 = -1
        k1 = -k1
    k0bin = to_bin(k0, l)
    k1bin = to_bin(k1, l)
    b0bin, b1bin = twodim_recoding_extended(k0bin, s0, k1bin, s1)
    # si = b0bin[i] in the paper
    # Ki = b1bin[i] in the paper
    # ===============================

    # Evaluation stage ==============
    Q = b0bin[l - 1] * table[b1bin[l - 1]]
    for i in range(l - 2, -1, -1):
        Q = 2 * Q
        Q = Q + b0bin[i] * table[abs(b1bin[i])]
    if even == 0:
        Q = Q + P
    # ===============================

    assert Q == k * P, "scalar_mul does not wrong"
    return Q


def scalar_mul_positive(k0, k1, P, lamb, beta, n):
    """Algorithm 2 from https://eprint.iacr.org/2013/158.pdf for 2 dimensions, i.e. m = 2
    Only for positive partial scalars k0,k1
    """

    l = ceil(log(n, 2) / 2) + 1

    # Precomputation stage ==========
    table = [P, P.curve()(beta**2 * P[0], -P[1])]
    # P = (x,y), P+lamb*P=(beta^2*x,-y)
    # ===============================

    # Recoding stage ================
    even = k0 % 2
    if even == 0:
        k0 = k0 - 1
    k0bin = to_bin(k0, l)
    k1bin = to_bin(k1, l)

    b0bin, b1bin = twodim_recoding(k0bin, k1bin)
    # si = b0bin[i] in the paper
    # Ki = b1bin[i] in the paper
    # ===============================

    # Evaluation stage ==============
    Q = b0bin[l - 1] * table[b1bin[l - 1]]
    for i in range(l - 2, -1, -1):
        Q = 2 * Q
        Q = Q + b0bin[i] * table[abs(b1bin[i])]
    if even == 0:
        Q = Q + P
    # ===============================

    assert Q == (k0 + (not even) + k1 * lamb) * P, Q
    return Q


def zip_bits_big_end(k0, k1, pad):
    bits = []
    while k0 > 0 or k1 > 0:
        bits.append([k0 % 2, k1 % 2])
        k0 >>= 1
        k1 >>= 1
    assert pad >= len(bits)
    bits += [[0, 0]] * (pad - len(bits))
    return reversed(bits)


def shamir_scalar_mul_positive(k0, k1, P, lamb, beta, n):
    lamP = P.curve()(beta * (P[0]), P[1])
    lam2P = P.curve()(beta * (lamP[0]), lamP[1])
    table = [[P.curve()(0), lamP], [P, -lam2P]]
    Q = P.curve()(0)
    l = ceil(log(n, 2) / 2) + 1
    for b0, b1 in zip_bits_big_end(k0, k1, l):
        Q *= 2
        Q += table[b0][b1]
    assert Q == (k0 + k1 * lamb) * P, Q
    return Q


def mods(k, w):
    mod = 1 << w
    u = k % mod
    if u >= (mod >> 1):
        u -= mod
    return u


def wnaf(k, w):
    i = 0
    naf = []
    while k >= 1:
        if k % 2 != 0:
            naf.append(mods(k, w))
            k -= naf[-1]
        else:
            naf.append(0)
        k //= 2
        i += 1
    return naf


def from_wnaf(naf):
    v = 0
    for n in reversed(naf):
        v *= 2
        v += n
    return v


def fromcutwnaf(wk):
    for i, n in enumerate(wk):
        if n != 0:
            break
    v = from_wnaf(wk[i:])
    return v, v - 1


def pad_nafs(naf0, naf1, l=None):
    if l is None:
        l = max(len(naf0), len(naf1))
    return naf0 + [0] * (l - len(naf0)), naf1 + [0] * (l - len(naf1))


def interleaving_positive(k0, k1, P, lamb, beta, n, w):
    curve = P.curve()
    naf0, naf1 = wnaf(k0, w), wnaf(k1, w)

    table0 = {2 * i - 1: (2 * i - 1) * P for i in range(1, 2 ** (w - 2) + 1)}
    table0.update({-i: -P for i, P in table0.items()})
    table0[0] = curve(0)
    table1 = {i: lamb * P for i, P in table0.items()}

    Q = curve(0)
    for (n0, n1) in reversed(list(zip(*pad_nafs(naf0, naf1)))):
        Q *= 2
        Q += table0[n0] + table1[n1]
    return Q


def regular_wnaf(k, w):
    i = 0
    m = 2**w
    l = int(k).bit_length()
    assert w in [1,2,3,4,5]  # nothing else tested
    assert k % 2 == 1
    digits = []
    while k > m:
        digits.append((k % (2 * m)) - m)
        k = (k - digits[-1]) // m
        i += 1
    digits.append(k)
    assert len(digits) == ceil(l/w), (len(digits), l)
    return digits


def from_regular_wnaf(rwnaf,w):
    v = 0
    m = 2**w
    for n in reversed(rwnaf):
        v*=m
        v+=n
    return v 

def fromcutregularwnaf(wk, w):
    for i, n in enumerate(wk):
        if n != 0:
            break
    v = from_regular_wnaf(wk[i:], w)
    return v, v - 1


def signed_shamir_scalar_mul_positive(k0, k1, P, lamb, beta, n, w):

    lamP = P.curve()(beta * (P[0]), P[1])
    table = {(i, j): i * P + j * lamP for i, j in [[-1, -1], [-1, 1], [1, -1], [1, 1]]}
    Q = P.curve()(0)
    s0, s1 = regular_wnaf(k0, w), regular_wnaf(k1, w)
    assert len(s0) == len(s1)
    for b0, b1 in reversed(list(zip(s0, s1))):
        Q *= 2
        Q += table[(b0, b1)]
    assert Q == (k0 + k1 * lamb) * P, Q
    return Q


# =========================================================
# =============== Side-channel part begins ================
# =========================================================


def scalar_mul_oracle(k, P, lamb, beta, n, registers):
    """
    scalar_mul with a side-channel oracle
    currently using the positive version
    """
    l = ceil(log(n, 2) / 2) + 2

    # Precomputation stage ==========
    table = [P, P.curve()(beta**2 * P[0], -P[1])]
    # P = (x,y), P+lamb*P=(beta^2*x,-y)
    # ===============================

    # Recoding stage ================
    k0, k1 = glv_decompose(k, lamb, n)
    even = k0 % 2
    if even == 0:
        k0 = k0 - 1
    s0, s1 = 1, 1
    if k0 < 0:
        s0 = -1
        k0 = -k0
    if k1 < 0:
        s1 = -1
        k1 = -k1
    k0bin = to_bin(k0, l)
    k1bin = to_bin(k1, l)

    b0bin, b1bin = twodim_recoding_extended(k0bin, s0, k1bin, s1)
    # si = b0bin[i] in the paper
    # Ki = b1bin[i] in the paper
    # ===============================

    # Evaluation stage ==============
    Q = b0bin[l - 1] * table[b1bin[l - 1]]
    for i in range(l - 2, -1, -1):
        Q = 2 * Q
        Pi = b0bin[i] * table[b1bin[i]]
        Q = Q + Pi
        if registers.is_zero(Pi, Q):
            print(f"LEAK: {i}")
    if even == 0:
        Q = Q + P
    # ===============================

    assert Q == k * P, Q
    return Q


def scalar_mul_oracle_positive(P, zvpparams, iteration):
    """
    scalar_mul_positive with side-channel oracle
    """
    secrets = zvpparams.secrets
    glv = zvpparams.glv
    k = (secrets.k0 + secrets.k1 * glv.lam) % glv.order
    l = secrets.half_bits + 1

    # Precomputation stage ==========
    table = [P, P.curve()(glv.beta**2 * P[0], -P[1])]
    # P = (x,y), P+lamb*P=(beta^2*x,-y)
    # ===============================

    # Recoding stage ================
    even = secrets.k0 % 2
    if even == 0:
        k0 = secrets.k0 - 1
    else:
        k0 = secrets.k0
    k0bin = to_bin(k0, l)
    k1bin = to_bin(secrets.k1, l)

    b0bin, b1bin = twodim_recoding(k0bin, k1bin)
    # si = b0bin[i] in the paper
    # Ki = b1bin[i] in the paper
    # ===============================

    # Evaluation stage ==============
    Q = b0bin[l - 1] * table[b1bin[l - 1]]
    found = False
    for i in range(l - 2, -1, -1):

        Q = 2 * Q
        Pi = b0bin[i] * table[b1bin[i]]
        assert (
            Q
            == 2
            * (from_bin(b0bin[(i + 1) :]) + from_bin(b1bin[(i + 1) :]) * glv.lam)
            * P
        )
        assert Pi[0] == ((glv.beta**2) ** (b1bin[i] != 0)) * (P[0])
        if iteration == i:
            found = zvpparams.registers.is_zero(Pi, Q)
        Q = Q + Pi
    if even == 0:
        Q = Q + P
    # ===============================

    assert (
        Q == (k0 + (not even) + glv.lam * secrets.k1) * P
    ), f"Computed: {Q}, should be {k*P}, for P={P}, lam={glv.lam}, beta={glv.beta}"
    return Q, found


def shamir_scalar_mul_oracle_positive(P, zvpparams, iteration):
    secrets = zvpparams.secrets
    lamP = P.curve()(zvpparams.glv.beta * (P[0]), P[1])
    lam2P = P.curve()(zvpparams.glv.beta * (lamP[0]), lamP[1])
    table = [[P.curve()(0), lamP], [P, -lam2P]]
    Q = P.curve()(0)
    found = False
    l = secrets.half_bits
    i = l - 1
    for [b0, b1] in zip_bits_big_end(secrets.k0, secrets.k1, l):
        Q *= 2
        if iteration == i:
            found = zvpparams.registers.is_zero(table[b0][b1], Q)
        Q += table[b0][b1]
        i -= 1
    assert Q == (secrets.k0 + secrets.k1 * zvpparams.glv.lam) * P, Q
    return Q, found


def signed_shamir_scalar_mul_oracle_positive(P, zvpparams, iteration, w):
    secrets = zvpparams.secrets
    glv = zvpparams.glv
    lamP = P.curve()(glv.beta * (P[0]), P[1])
    table = {(i, j): i * P + j * lamP for i, j in [[-1, -1], [-1, 1], [1, -1], [1, 1]]}
    Q = P.curve()(0)
    s0, s1 = regular_wnaf(secrets.k0, w), regular_wnaf(secrets.k1, w)
    assert len(s0) == len(s1)
    l = ceil(secrets.half_bits/w)
    i = l - 1
    for b0, b1 in reversed(list(zip(s0, s1))):
        Q *= 2
        if iteration == i:
            found = zvpparams.registers.is_zero(table[(b0, b1)], Q)
        Q += table[(b0, b1)]
        i -= 1
    assert Q == (secrets.k0 + secrets.k1 * glv.lam) * P, Q
    return Q, found


def interleaving_oracle_positive(P, zvpparams, iteration, w=3):
    curve = P.curve()
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    lamb = zvpparams.glv.lam
    naf0, naf1 = wnaf(k0, w), wnaf(k1, w)

    table0 = {2 * i - 1: (2 * i - 1) * P for i in range(1, 2 ** (w - 2) + 1)}
    table0.update({-i: -P for i, P in table0.items()})
    table0[0] = curve(0)
    table1 = {i: lamb * P for i, P in table0.items()}

    Q = curve(0)
    l = zvpparams.secrets.half_bits + 1
    i = l - 1
    found = False
    for (n0, n1) in reversed(list(zip(*pad_nafs(naf0, naf1, l)))):
        Q *= 2
        if iteration == i:
            found = zvpparams.registers.is_zero(table0[n0], Q)
        Q += table0[n0]
        if iteration == i:
            found = found or zvpparams.registers.is_zero(table1[n1], Q)
        Q += table1[n1]
        i -= 1
    return Q, found


def regular_interleaving_oracle_positive(P, zvpparams, iteration, w=3):
    curve = P.curve()
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    lamb = zvpparams.glv.lam
    naf0, naf1 = regular_wnaf(k0, w), regular_wnaf(k1, w)

    table0 = {2 * i - 1: (2 * i - 1) * P for i in range(1, 2 ** (w - 1) + 1)}
    table0.update({-i: -P for i, P in table0.items()})
    table0[0] = curve(0)
    table1 = {i: lamb * P for i, P in table0.items()}

    Q = curve(0)
    l = ceil(zvpparams.secrets.half_bits/w)
    assert len(naf0)==l and len(naf1)==l, (l, zvpparams.secrets.half_bits, len(naf0),len(naf1))
    i = l - 1
    found = False
    m = 2**w
    for (n0, n1) in reversed(list(zip(naf0,naf1))):
        Q *= m
        if iteration == i:
            found = zvpparams.registers.is_zero(table0[n0], Q)
        Q += table0[n0]
        if iteration == i:
            found = found or zvpparams.registers.is_zero(table1[n1], Q)
        Q += table1[n1]
        i -= 1
    return Q, found


def interleaving_easy_oracle_positive(P, zvpparams, w=3):
    curve = P.curve()
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    lamb = zvpparams.glv.lam
    naf0, naf1 = wnaf(k0, w), wnaf(k1, w)

    table0 = {2 * i - 1: (2 * i - 1) * P for i in range(1, 2 ** (w - 2) + 1)}
    table0.update({-i: -P for i, P in table0.items()})
    table0[0] = curve(0)
    table1 = {i: lamb * P for i, P in table0.items()}
    l = zvpparams.secrets.half_bits + 1
    Q = curve(0)
    found = []
    for (n0, n1) in reversed(list(zip(*pad_nafs(naf0, naf1, l)))):
        Q *= 2
        found.append(zvpparams.registers.is_zero(table0[n0], table1[n1]))
        Q += table0[n0] + table1[n1]
    return Q, list(reversed(found))


def regular_interleaving_easy_oracle_positive(P, zvpparams, w=3):
    curve = P.curve()
    k0, k1 = zvpparams.secrets.k0, zvpparams.secrets.k1
    lamb = zvpparams.glv.lam
    naf0, naf1 = regular_wnaf(k0, w), regular_wnaf(k1, w)

    table0 = {2 * i - 1: (2 * i - 1) * P for i in range(1, 2 ** (w - 1) + 1)}
    table0.update({-i: -P for i, P in table0.items()})
    table0[0] = curve(0)
    table1 = {i: lamb * P for i, P in table0.items()}
    l = ceil(zvpparams.secrets.half_bits/w)
    assert len(naf0)==l and len(naf1)==l, (l, zvpparams.secrets.half_bits, len(naf0),len(naf1))
    Q = curve(0)
    found = []
    m = 2**w
    for (n0, n1) in reversed(list(zip(naf0, naf1))):
        Q *= m
        found.append(zvpparams.registers.is_zero(table0[n0], table1[n1]))
        Q += table0[n0] + table1[n1]
    return Q, list(reversed(found))