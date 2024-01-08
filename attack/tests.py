import bsgssolver
import dcp
import msm
import zvp_glv_sm as zvpsm
import zvp_glv_msm as zvpmsm
from sage.all import ZZ, EllipticCurve, GF, PolynomialRing, ceil, log
from pari_tools.dcp_pari import multidcpsolver_pari2, dcpsolver_pari2


def test_dcp_pari():
    assert multidcpsolver_pari2(1039, 0, 6, 140, 6, 6, 0) == ["47", "594"]
    assert dcpsolver_pari2(1048627, 0, 3, 7) == ["400363"]


def test_bsgs2_balanced():
    p = ZZ(34359738943)
    curve = EllipticCurve(GF(p), [0, 3])
    priv, priv1, priv2, lam, order = (
        ZZ(22160906118),
        ZZ(113686),
        ZZ(68577),
        ZZ(7066965624),
        ZZ(34359695179),
    )
    P = curve(12145174682, 15668259325)
    assert P != curve(0)
    Q = priv * P
    k1_upper, k2_upper = ZZ(0b11011), ZZ(0b1)
    lower_bits1 = priv1.nbits() - k1_upper.nbits()
    lower_bits2 = priv2.nbits() - k2_upper.nbits()
    Q2 = Q - 2**lower_bits1 * k1_upper * P - 2**lower_bits2 * k2_upper * lam * P
    k1_lower, k2_lower = bsgssolver.bsgs2_balanced(P, Q2, lam, lower_bits1, lower_bits2)
    k1 = k1_lower + k1_upper * 2**lower_bits1
    k2 = k2_lower + k2_upper * 2**lower_bits2
    k = (k1 + k2 * lam) % order
    assert k == priv


def test_bsgs2_balanced2():
    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(32), ZZ(8)
    lam, order = ZZ(195), ZZ(1033)
    priv = (k1 + k2 * lam) % order
    P = curve(38, 101)
    Q = priv * P
    assert bsgssolver.bsgs2_balanced(P, Q, lam, k1.nbits(), k2.nbits()) == [32, 8]


def test_bsgs2_balanced3():
    p = ZZ(1099511628079)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    lam = ZZ(533962724094)
    order = ZZ(1099512457333)
    k1, k2 = ZZ(100572), ZZ(8439)
    P = curve.lift_x(F(6))
    priv = (k1 + k2 * lam) % order
    Q = priv * P
    assert bsgssolver.bsgs2_balanced(P, Q, lam, k1.nbits(), k2.nbits()) == [
        100572,
        8439,
    ]


def test_bsgs2_pari():
    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(32), ZZ(8)
    lam, order = ZZ(195), ZZ(1033)
    priv = (k1 + k2 * lam) % order
    P = curve(38, 101)
    Q = priv * P
    assert bsgssolver.bsgs2_pari(P, Q, lam, k1.nbits(), k2.nbits()) == [32, 8]


def test_bsgs2_pari2():
    p = ZZ(1099511628079)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    lam = ZZ(533962724094)  # beta =35271754617
    order = ZZ(1099512457333)
    k1, k2 = ZZ(100572), ZZ(8439)
    P = curve.lift_x(F(6))
    priv = (k1 + k2 * lam) % order
    Q = priv * P
    ck1, ck2 = bsgssolver.bsgs2_pari(P, Q, lam, k1.nbits(), k2.nbits())
    assert [ck1, ck2] == [100572, 8439], f"{ck1},{ck2}"


def test_solve_multi_dcp_pari():
    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(6), ZZ(6)
    lam, order = ZZ(195), ZZ(1033)
    beta = ZZ(140)
    P = curve(47, 353)
    Q = 143 * P
    assert P[0] + Q[0] + 2 == 0
    result = dcp.solve_multi_dcp_pari(
        [k1, k2],
        0,
        {"curve": curve, "order": order, "lam": lam, "beta": beta, "field": F},
    )
    assert int(result[0]) in [594, 47]


def test_solve_multi_dcp():
    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(6), ZZ(6)
    lam, order = ZZ(195), ZZ(1033)
    beta = ZZ(140)
    P = curve(47, 353)
    assert (
        dcp.solve_multi_dcp(
            [k1, k2],
            [1, 0],
            {"curve": curve, "order": order, "lam": lam, "beta": beta, "field": F},
        )[0]
        == P[0]
    )


def test_bin():
    assert msm.to_bin(2, 2) == [0, 1]
    assert msm.from_bin([0, 1, 1]) == 6
    assert 6 == msm.from_bin(msm.to_bin(6, 5))
    assert msm.zero_allign_bin(3, 3) == "0b011"


def test_twodim_recoding():
    l = 5
    kbin, kpbin = msm.to_bin(1, l), msm.to_bin(0, l)
    bbin, bpbin = msm.twodim_recoding(kbin, kpbin)
    bbin2, bpbin2 = msm.twodim_recoding2(kbin, kpbin)
    assert bbin == (l - 1) * [-1] + [1]
    assert bpbin == [0] * l
    assert bbin2 == (l - 1) * [-1] + [1]
    assert bpbin2 == [0] * l
    assert kbin, kpbin == msm.twodim_recoding_inverse(bbin, bpbin)

    kbin, kpbin = [1, 0, 0, 1, 1], [1, 0, 1, 1, 0]
    bbin, bpbin = msm.twodim_recoding(kbin, kpbin)
    bbin2, bpbin2 = msm.twodim_recoding2(kbin, kpbin)
    assert bbin == [-1, -1, 1, 1, 1]
    assert bpbin == [-1, -1, 0, 0, 1]
    assert bbin2 == [-1, -1, 1, 1, 1]
    assert bpbin2 == [-1, -1, 0, 0, 1]
    assert kbin, kpbin == msm.twodim_recoding_inverse(bbin, bpbin)


def test_two_dim_recoding_inverse_upper_bits():
    bbin = [-1, -1, 1, 1, 1][1:]
    bpbin = [-1, -1, 0, 0, 1][1:]
    l = 5
    assert msm.two_dim_recoding_inverse_upper_bits(l, bbin, bpbin) == (
        [0, 1, 1],
        [[0, 1, 0], [1, 1, 0]],
    )
    l = 16
    bbin = [1, 1, 1, -1, 1, 1]
    bpbin = [0, 0, 0, -1, 0, 1]
    assert msm.two_dim_recoding_inverse_upper_bits(l, bbin, bpbin) == (
        [1, 1, 1, 0, 1],
        [[1, 1, 0, 1, 0], [0, 0, 1, 1, 0]],
    )


def test_twodim_recoding_extended():
    # TODO nontrivial example
    l = 5
    bbin, bpbin = msm.twodim_recoding_extended(
        msm.to_bin(1, l), -1, msm.to_bin(0, l), 1
    )
    assert bbin == (l - 1) * [1] + [-1], bbin
    assert bpbin == [0] * l


def test_scalar_mul():
    E = EllipticCurve(GF(31), [0, 3])
    assert msm.scalar_mul(23, E(25, 2), lamb=36, beta=25, n=43) == E(14, 22)
    assert msm.scalar_mul_positive(3, 2, E(25, 2), lamb=36, beta=25, n=43) == E(9, 9)


def test_symbolic_double():
    F = GF(31)
    E = EllipticCurve(F, [0, 3])
    P = E(25, 2)
    R = PolynomialRing(F, ["x0", "y0"]).fraction_field()
    x0, y0 = R.gens()
    assert list(
        map(lambda x: x(P[0], P[1]), msm.symbolic_double((x0, y0), E.a4(), E.a6()))
    ) == list((2 * P)[:2])


def test_symbolic_addition():
    F = GF(31)
    a, b = F(0), F(3)
    E = EllipticCurve(F, [a, b])
    P = E(25, 2)
    R = PolynomialRing(F, ["x0", "y0"]).fraction_field()
    x0, y0 = R.gens()
    Ps = x0, y0
    Qs = msm.symbolic_double(Ps, a, b)
    add = msm.symbolic_addition(Ps, Qs, a, b)
    assert list(map(lambda x: x(P[0], P[1]), add)) == list((3 * P)[:2]), P


def test_symbolic_scalar_mul():
    F = GF(31)
    a, b = F(0), F(3)
    E = EllipticCurve(F, [a, b])
    P = E(25, 2)
    R = PolynomialRing(F, ["x0", "y0"]).fraction_field()
    x0, y0 = R.gens()
    Ps = x0, y0
    Qs = msm.symbolic_scalar_mul(5, Ps, lamb=36, beta=25, n=43, a=a, b=b)
    assert list(map(lambda x: x(P[0], P[1]), Qs)) == list((5 * P)[:2]), P
    Qs = msm.symbolic_scalar_mul_positive(3, 2, Ps, lamb=36, beta=25, n=43, a=a, b=b)
    assert list(map(lambda x: x(P[0], P[1]), Qs)) == list((75 * P)[:2]), P


def test_scalar_mul_oracle_positive():

    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(29), ZZ(26)
    lam, order = ZZ(195), ZZ(1033)
    beta = ZZ(140)
    P = curve(47, 353)
    Q = 143 * P
    Q2 = (6 + 6 * lam) * P
    assert Q2 == Q
    assert P[0] + Q[0] + 2 == 0
    assert msm.scalar_mul_oracle_positive(k1, k2, P, lam, beta, order, 2)[1]


def test_scalar_mul_oracle_positive2():
    p = ZZ(2143)
    F = GF(p)
    curve = EllipticCurve(F, [0, 5])
    k1, k2 = ZZ(0b1010101), ZZ(0b100010)
    lam, order = ZZ(1262), ZZ(2089)
    beta = ZZ(1793)
    l = ceil(log(order, 2) / 2) + 1
    P = curve(1050, 236)
    assert msm.scalar_mul_oracle_positive(k1, k2, P, lam, beta, order, 4)[1]


def test_scalar_mul_oracle_positive3():
    p = ZZ(4294968187)
    F = GF(p)
    curve = EllipticCurve(F, [0, 3])
    k = ZZ(2613661231)
    k1, k2 = ZZ(15788), ZZ(2360)
    lam = ZZ(2168662147)
    beta = ZZ(4037544512)
    order = ZZ(4295070547)
    P = curve(912661074, 590983812)
    assert (k1 + k2 * lam) % order == k
    assert (k1 + k2 * lam) * P == k * P
    assert msm.scalar_mul(k, P, lam, beta, order) == k * P
    assert msm.scalar_mul_positive(k1, k2, P, lam, beta, order) == k * P
    assert msm.scalar_mul_oracle_positive(k1, k2, P, lam, beta, order, 0)[0] == k * P


def test_scalar_mul_oracle_positive4():
    p = ZZ(1099511628079)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k1, k2 = ZZ(0b111001110011100110001), ZZ(0b111111011010100101101)
    lam = ZZ(533962724094)
    beta = ZZ(35271754617)
    order = ZZ(1099512457333)
    P = curve(961630273290, 234600338115)
    k = (k1 + k2 * lam) % order
    assert msm.scalar_mul(k, P, lam, beta, order) == k * P
    assert msm.scalar_mul_positive(k1, k2, P, lam, beta, order) == k * P
    assert msm.scalar_mul_oracle_positive(k1, k2, P, lam, beta, order, 0)[0] == k * P


def test_zpa_attack():
    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k0, k1 = ZZ(29), ZZ(26)
    lam, order = ZZ(195), ZZ(1033)
    beta = ZZ(140)
    P = curve(110, 299)
    assert lam * P == curve(beta * 110, 299)
    l = ceil(log(order, 2) / 2) + 1
    assert k0.nbits() <= l - 1 and k1.nbits() <= l - 1
    params = {
        "curve": curve,
        "field": F,
        "k0": k0,
        "k1": k1,
        "order": order,
        "zvp_target": ZZ(4),
    }
    result = zvpsm.zvp_glv_sm(params, verbose=False)
    assert result["k0_upper"] == params["k0"] // 2
    assert result["k1_upper"] == params["k1"] // 2


def test_zpa_attack3():
    p = ZZ(1039)
    F = GF(p)
    curve = EllipticCurve(F, [0, 6])
    k0, k1 = ZZ(13), ZZ(11)
    order = ZZ(1033)
    params = {
        "curve": curve,
        "field": F,
        "k0": k0,
        "k1": k1,
        "order": order,
        "zvp_target": ZZ(3),
    }
    result = zvpsm.zvp_glv_sm(params, verbose=False)
    assert result["k0_upper"] == 0
    assert result["k1_upper"] == params["k1"] // 2


def test_zpa_attack4():
    p = ZZ(2143)
    F = GF(p)
    curve = EllipticCurve(F, [0, 5])
    k0, k1 = ZZ(0b1010101), ZZ(0b100010)
    lam, order = ZZ(1262), ZZ(2089)
    beta = ZZ(1793)
    params = {
        "curve": curve,
        "field": F,
        "k0": k0,
        "k1": k1,
        "order": order,
        "zvp_target": ZZ(4),
        "lam": lam,
        "beta": beta,
    }
    result = zvpmsm.zvp_glv_msm(params, verbose=False)
    assert [k0 // 8, k1 // 8] in result["kipairs"]


test_zpa_attack4()


def test_ltl_oracle():
    curve = EllipticCurve(GF(101), [1, 1])
    P = curve.lift_x(ZZ(3))
    k = ZZ(5)
    assert zvpsm.ltr_oracle(curve, k, P, 0)[0] == (k * P)
    assert zvpsm.ltr_oracle_always(curve, k, P, 0)[0] == (k * P)
