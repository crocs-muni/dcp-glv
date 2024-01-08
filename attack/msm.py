"""
Implementation of multi-scalar multiplication as designed in https://eprint.iacr.org/2013/158.pdf for dimension 2
- including symbolic version
- including one leaking zero-registers
"""

from glv import glv_decompose
from sage.all import PolynomialRing, ceil, log


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


# =========================================================
# ================== Symbolic part begins==================
# =========================================================


def reduce_mod_curve(poly, a, b):
    """Reduces polynomial :poly: modulo curve equation, i.e. replaces y**2 by x**3+ax+b"""
    x, y = poly.parent().gens()
    equation = x**3 + a * x + b
    exponents = poly.exponents()  # list of (i,j) for each term x**i*y**j in poly
    for e_x, e_y in exponents:
        if e_y < 2:
            continue
        coef = poly.coefficient({x: e_x, y: e_y})
        e = e_y % 2
        poly -= coef * x ** (e_x) * y ** (e_y)
        poly += coef * x ** (e_x) * y**e * equation ** ((e_y // 2))
    return poly


def reduce_map(f, a, b):
    """Reduces map on a curve modulo curve equation"""
    y0 = f.parent().gens()[1]
    num = f.numerator()
    den = f.denominator()
    if (den / y0).denominator() == 1:
        num *= y0
        den *= y0
    num = num.numerator()
    den = den.numerator()
    num = reduce_mod_curve(num, a, b)
    den = reduce_mod_curve(den, a, b)
    return num / den


def symbolic_double(P, a, b):
    """Symbolically doubles point P (tuple of rational functions)"""
    x, y = P
    if x == 0 and y == 1:
        return x, y
    s = (3 * x**2 + x.parent()(a)) / (2 * y)
    x3 = s**2 - 2 * x
    y3 = y + s * (x3 - x)
    return (reduce_map(x3, a, b), reduce_map(-y3, a, b))


def symbolic_addition(P, Q, a, b):
    """Symbolically adds points P and Q (tuples of rational functions)"""
    (x1, y1), (x2, y2) = P, Q
    if x1 == x2 and y1 == y2:
        return symbolic_double(P, a, b)
    if x1 == 0 and y1 == 1:
        return Q
    if x2 == 0 and y2 == 1:
        return P
    s = (y1 - y2) / (x1 - x2)
    x3 = s**2 - x1 - x2
    y3 = y1 + s * (x3 - x1)
    return (reduce_map(x3, a, b), reduce_map(-y3, a, b))


def symbolic_scalar_mul(k, P, lamb, beta, n, a, b):
    """Symbolic version of scalar_mul, P is here a tuple of rational functions"""

    l = ceil(log(n, 2) / 2) + 2

    # Precomputation stage ==========
    table = [P, (beta**2 * P[0], -P[1])]
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
    def mul_point(s, point):
        if s == 1:
            return point
        if s == -1:
            return (point[0], -point[1])
        raise Exception(s)

    Q = mul_point(b0bin[l - 1], table[b1bin[l - 1]])
    for i in range(l - 2, -1, -1):
        Q = symbolic_double(Q, a, b)
        Q = symbolic_addition(Q, mul_point(b0bin[i], table[b1bin[i]]), a, b)
    if even == 0:
        Q = symbolic_addition(Q, P, a, b)
    # ===============================

    return Q


def symbolic_scalar_mul_positive(k0, k1, P, lamb, beta, n, a, b):
    """Symbolic version of scalar_mul_positive, P is here a tuple of rational functions"""

    l = ceil(log(n, 2) / 2) + 1

    # Precomputation stage ==========
    table = [P, (beta**2 * P[0], -P[1])]
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
    def mul_point(s, point):
        if s == 1:
            return point
        if s == -1:
            return (point[0], -point[1])
        raise Exception(s)

    Q = mul_point(b0bin[l - 1], table[b1bin[l - 1]])
    for i in range(l - 2, -1, -1):
        Q = symbolic_double(Q, a, b)
        Q = symbolic_addition(Q, mul_point(b0bin[i], table[b1bin[i]]), a, b)
    if even == 0:
        Q = symbolic_addition(Q, P, a, b)
    # ===============================

    return Q


# =========================================================
# =============== Side-channel part begins ================
# =========================================================


def scalar_mul_oracle(k, P, lamb, beta, n, oracle_polynomial):
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
        if oracle_polynomial(Pi[0], Pi[1], Q[0], Q[1]) == 0:
            print(f"LEAK: {i}")
    if even == 0:
        Q = Q + P
    # ===============================

    assert Q == k * P, Q
    return Q


def scalar_mul_oracle_positive(k0, k1, P, lamb, beta, n, iteration):
    """
    scalar_mul_positive with side-channel oracle for f=x1+x2+2
    """
    x1, x2, y1, y2 = PolynomialRing(
        P.curve().base_field(), ["x1", "x2", "y1", "y2"]
    ).gens()
    f = x1 + x2 + 2
    k = (k0 + k1 * lamb) % n
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
    found = False
    for i in range(l - 2, -1, -1):

        Q = 2 * Q
        Pi = b0bin[i] * table[b1bin[i]]
        assert (
            Q
            == 2 * (from_bin(b0bin[(i + 1) :]) + from_bin(b1bin[(i + 1) :]) * lamb) * P
        )
        assert Pi[0] == ((beta**2) ** (b1bin[i] != 0)) * (P[0])
        if iteration == i:
            found = f(Pi[0], Q[0], Pi[1], Q[1]) == 0
        Q = Q + Pi
    if even == 0:
        Q = Q + P
    # ===============================

    assert (
        Q == (k0 + (not even) + lamb * k1) * P
    ), f"Computed: {Q}, should be {k*P}, for P={P}, lam={lamb}, beta={beta}"
    return Q, found
