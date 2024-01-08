"""
Implementation of various functions needed for GLV decomposition
"""

from sage.all import (
    vector,
    EllipticCurve,
    GF,
    PolynomialRing,
    ZZ,
    sqrt,
    QQ,
    matrix,
    random_prime,
)
from random import randint


def extended_euclid(a, b):
    """
    Extended euclidean algorithm

    :param a: integer
    :param b: integer
    :return s1,t1,r1,seq: where a*s1+b*t1=r1=gcd(a,b) and a sequence of s,t,r with a*s+b*t=r (see the algorithm)
    """
    r0, r1, r2 = a, b, None
    s0, s1, s2 = 1, 0, None
    t0, t1, t2 = 0, 1, None
    seq = []
    while True:
        r2, q1 = r0 % r1, r0 // r1
        s2 = s0 - q1 * s1
        t2 = t0 - q1 * t1
        seq.append([r1, s1, t1])
        if r2 == 0:
            break
        s1, s0 = s2, s1
        t1, t0 = t2, t1
        r1, r0 = r2, r1
    assert s1 * a + t1 * b == r1
    return r1, s1, t1, seq[1:]


def ee_basis(n, lam):
    """
     Implementation of extended euclidean as a generator of a sequence of r,s,t
    :param n: integer
    :param lam: integer
    :return r,s,t: where s*n+t*lam = r
    """
    r0, r1, r2 = n, lam, None
    s0, s1, s2 = 1, 0, None
    t0, t1, t2 = 0, 1, None
    seq = []
    while True:
        r2, q1 = r0 % r1, r0 // r1
        s2 = s0 - q1 * s1
        t2 = t0 - q1 * t1
        yield r1, s1, t1
        s1, s0 = s2, s1
        t1, t0 = t2, t1
        r1, r0 = r2, r1


def find_basis(k, n, lam):
    """
    :param k: integer
    :param n: modulus
    :param lam: fast scalar
    :return v1,v2: basis satisfying f(v1)=f(v2)=0
    where f(i,j) = i+lam*j (mod n)
    """
    ee_gen = ee_basis(n, lam)
    for r, s, t in ee_gen:
        if r < sqrt(n):
            break
        r0, s0, t0 = r, s, t
    r1, s1, t1 = next(ee_gen)
    if r0**2 + t0**2 > r1**2 + t1**2:
        return vector((r, -t)), vector((r1, -t1))
    return vector((r, -t)), vector((r0, -t0))


def find_closest(k, vs: list) -> vector:
    """
    :param k: integer
    :param vs: list of vectors
    :return: a closest vector to k in the span of vs
    """
    m = [[QQ(i) for i in v] for v in vs]
    m = matrix(m)
    m = m.transpose()
    sol = m.solve_right(vector([QQ(k)] + [0] * (len(vs) - 1)))
    bs = list(map(lambda x: x.round(), sol))
    return sum(map(lambda x: x[0] * x[1], zip(bs, vs)))


def glv_decompose_simple(k, lam, n) -> vector:
    """
    Uses euclidean extended alg for glv decomposition k=k1+k2*lam (mod n)
    :param k: integer
    :param lam: fast scalar
    :param n: modulus
    :return: vector (k1,k2)
    """
    v1, v2 = find_basis(k, n, lam)
    # assert (v1[0]+v1[1]*lam)%n==0
    # assert (v2[0]+v2[1]*lam)%n==0
    v = find_closest(k, [v1, v2])
    # assert (v[0]+v[1]*lam)%n==0
    short = vector([k, 0]) - v
    assert (short[0] + short[1] * lam) % n == k % n
    return short


def glv_decompose(k, lam, n, dim=2) -> vector:
    """
    Uses LLL for r-dim glv decomposition k = sum k_i*lam^i (mod n)
    :param k: integer
    :param lam: fast scalar
    :param n: modulus
    :param r: dimension of glv decomp
    :return: vector (k0,..,k_{r-1})
    """
    r = dim
    m = [[0] * i + [lam, -1] + (r - 2 - i) * [0] for i in range(r - 2)]
    m += [[0] * i + [n] + (r - i - 1) * [0] for i in range(r)]
    m = matrix([[lam ** (r - 1)] + [0] * (r - 2) + [-1]] + m)
    for u in m:
        assert (
            sum(map(lambda x: x[0] * x[1], zip(list(u), [lam**i for i in range(r)])))
            % n
            == 0
        )
    vs = [v for v in m.LLL() if v != vector([0] * r)]
    for u in vs:
        assert (
            sum(map(lambda x: x[0] * x[1], zip(list(u), [lam**i for i in range(r)])))
            % n
            == 0
        )
    v = find_closest(k, vs)
    assert (
        sum(map(lambda x: x[0] * x[1], zip(list(v), [lam**i for i in range(r)]))) % n
        == 0
    )
    return vector([k] + [0] * (r - 1)) - v


def glv_decompose_two(k, lam, sigma, n):
    """
    Uses LLL for glv decomposition k = k1+k2*lam+k3*sigma+k4*lam*sigma (mod n)
    :param k: integer
    :param lam: fast scalar
    :param sigma: fast scalar
    :param n: modulus
    :return: vector (k1,k2,k3,k4)
    """
    siglam = GF(n)(sigma * lam)
    m = matrix(
        [
            [ZZ(lam), -1, 0, 0],
            [ZZ(sigma), 0, -1, 0],
            [ZZ(siglam), 0, 0, -1],
            [n, 0, 0, 0],
            [0, n, 0, 0],
            [0, 0, n, 0],
            [0, 0, 0, n],
        ]
    )
    for u in m:
        assert (
            sum(
                map(
                    lambda x: x[0] * x[1],
                    zip(list(u), [1, ZZ(lam), ZZ(sigma), ZZ(siglam)]),
                )
            )
            % n
            == 0
        )
    vs = [v for v in m.LLL() if v != vector([0] * 4)]
    for u in vs:
        assert (
            sum(
                map(
                    lambda x: x[0] * x[1],
                    zip(list(u), [1, ZZ(lam), ZZ(sigma), ZZ(siglam)]),
                )
            )
            % n
            == 0
        )
    v = find_closest(k, vs)
    return vector([k, 0, 0, 0]) - v


def find_lambda(curve: EllipticCurve):
    """
    Computes eigenvalue and corresponding x-rational map of endo lambda satisfying lambda^2+lambda+1=0
    :param curve:
    :return lambda, beta: where lambda*(x,y)=(beta*x,y)
    """
    F = curve.base_field()
    n = curve.order()
    p = F.order()
    g = F.multiplicative_generator()
    assert p % 3 == 1
    beta = g ** ((p - 1) // 3)
    z = PolynomialRing(GF(n), "z").gen()
    lam1, lam2 = [ZZ(i[0]) for i in (z**2 + z + 1).roots()]
    while True:
        P = curve.random_point()
        if P != curve(0):
            break
    if curve(beta * P[0], P[1]) == lam1 * P:
        return lam1, beta
    assert curve(beta * P[0], P[1]) == lam2 * P
    return lam2, beta


def find_sigma(curve: EllipticCurve):
    """
    Computes eigenvalue and corresponding x-rational map of endo sigma satisfying sigma^2+3*lambda+3=0
    :param curve:
    :return sigma, sigma_map: where lambda*(x,y)=(sigma_map(x),whatever)
    """
    F = curve.base_field()
    isogs = curve.isogenies_prime_degree(3)
    for isog in isogs:
        if isog.codomain().j_invariant() == curve.j_invariant():
            break
    curve2 = isog.codomain()
    u = (curve.a6() / curve2.a6()).nth_root(3)
    sigma = find_lambda(curve)[0] - 1
    P = curve.random_point()
    sigma_map = PolynomialRing(F, "x").fraction_field()(u * isog.rational_maps()[0])
    assert (sigma * P)[0] == sigma_map(P[0])
    return sigma, sigma_map


def find_curve(bits):
    """
    Finds curve E:y^2=x^3+b with prime order, p=1 (mod 3), p=3 (mod 4)
    :param bits: bitsize of prime p (p will be the smallest s.t. all conditions hold)
    :return curve:
    """
    pi = ZZ(2 ** (bits - 1))
    while True:
        if pi % 3 != 1 or pi % 4 != 3 or not pi.is_prime():
            pi = pi.next_prime()
            continue
        F = GF(pi)
        counter = 10
        for b in F:
            try:
                counter -= 1
                if counter == 0:
                    break
                E = EllipticCurve(F, [0, b])
                n = E.order()
                assert n.is_prime()
                return E
            except:
                continue
        pi = pi.next_prime()


def find_curve_random(bitsize):
    """
    Finds a random curve E:y^2=x^3+b with prime order, p=1 (mod 3), p=3 (mod 4)
    :param bitsize: bitsize of prime p
    :return curve:
    """
    while True:
        pi = random_prime(2**bitsize - 1, True, 2 ** (bitsize - 1))
        if pi % 3 != 1 or pi % 4 != 3 or not pi.is_prime():
            continue
        F = GF(pi)
        counter = 10
        for b in F:
            try:
                counter -= 1
                if counter == 0:
                    break
                E = EllipticCurve(F, [0, b])
                n = E.order()
                assert n.is_prime()
                return E
            except:
                continue


def dcp_instance(bits, k=None):
    """
    Find dcp instance given bits and possibly k.
    :param bits: Size of the base_field prime
    :param k(=None): scalar k (otherwise generated randomly)
    :return E,P,Q,k: such that k*P=Q and P[0]*Q[0]=-1
    """
    Eo = find_curve(bits)
    F = Eo.base_field()
    n = Eo.order()
    if not k:
        k = randint(1, n)
    while True:
        try:
            Po = Eo.random_point()
            Qo = k * Po
            u = (-Po[0] * Qo[0]).nth_root(4)
            assert u != 0
            break
        except:
            continue
    E = EllipticCurve(F, [0, Eo.a6() / u**6])
    P = E(Po[0] / u**2, Po[1] / u**3)
    Q = E(Qo[0] / u**2, Qo[1] / u**3)
    assert k * P == Q
    assert P[0] * Q[0] + 1 == 0
    return E, P, Q, k


def semaev(n, a, b):
    """
    Computes the n-th Semaev polynomial

    :param n: which Semave to compute
    :param a: a coefficient in short weiestrass form
    :param b: b coefficient in short weiestrass form
    :return polynomial:
    """

    if n == 2:
        x1, x2 = PolynomialRing(ZZ, ["x1", "x2"]).gens()
        return x1 - x2
    if n == 3:
        x1, x2, x3 = PolynomialRing(ZZ, ["x1", "x2", "x3"]).gens()
        return (
            (x1 - x2) ** 2 * x3**2
            - 2 * ((x1 + x2) * (x1 * x2 + a) + 2 * b) * x3
            + ((x1 * x2 - a) ** 2 - 4 * b * (x1 + x2))
        )

    Ring = PolynomialRing(ZZ, [f"x{i}" for i in range(1, n + 1)] + ["X"])

    xis = list(Ring.gens()[:-1])
    X = Ring.gens()[-1]

    f = semaev(n - 1, a, b)
    f3 = semaev(3, a, b)

    ff = f(*(xis[:-2] + [X]))
    ff3 = f3(xis[-2], xis[-1], X)

    Ring2 = PolynomialRing(ZZ, xis)

    return Ring2(ff.resultant(ff3, X))


def decomposition_matrix(n, lam, dim):
    """
    Computes the matrix for decomposition

    :param n: modulus
    :param lam: eigenvalue
    :param dim: dimension of decomposition
    :return matrix:
    """
    lam = ZZ(lam)
    m = [i * [0] + [lam, -1] + [0] * (dim - i - 2) for i in range(dim - 1)]
    for i in range(dim):
        m.append(i * [0] + [n] + [0] * (dim - i - 1))
    m = matrix(m).LLL()
    return matrix([v for v in m if v != vector([0] * dim)])


def decompose_precomputed(scalar, lam, m, dim, n):
    """
    Computes decomposition given precomputed matrix

    :param scalar: scalar to decompose
    :param lam: eigenvalue
    :param m: decomposition matrix
    :param dim: dimension of decomposition
    :param n: modulus
    :return vector:
    """
    v = find_closest(scalar, list(m))
    assert (
        sum(map(lambda x: x[0] * x[1], zip(list(v), [lam**i for i in range(dim)])))
        % n
        == 0
    )
    kis = vector([scalar] + [0] * (dim - 1)) - v
    return kis


def gen_glv_curve(bits, random=True):
    """Finds suitable curve for 2-dim decomposition"""
    curve = {}
    if random:
        curve["curve"] = find_curve_random(bits)
    else:
        curve["curve"] = find_curve(bits)
    curve["field"] = curve["curve"].base_field()
    curve["lam"], curve["beta"] = find_lambda(curve["curve"])
    curve["order"] = curve["curve"].order()
    return curve
