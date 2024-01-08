"""
DCP solver
"""
import msm
from pari_tools.dcp_pari import dcpsolver_pari2, dcpsolver_pari4, multidcpsolver_pari2


def solve_dcp(scalar, curve):
    """
    Solves DCP for :scalar: and polynomial x1+x2+2
    The polynomial is hardcoded in the pari implementation
    Works for prime fields and quadratic fields
    """
    A, B = curve.a4(), curve.a6()
    F = curve.base_field()
    p = F.characteristic()
    P = None
    if scalar == 0:
        # x2 = 0
        try:
            P = curve.lift_x(F(-2))
        except ValueError as e:
            pass
        return P

    if F.degree() == 1:
        roots = dcpsolver_pari2(p, A, B, scalar)
        for root in [F(int(r)) for r in roots]:
            try:
                P = curve.lift_x(root)
                assert P[0] + (scalar * P)[0] + 2 == 0
            except ValueError as e:
                continue
    elif F.degree() == 2:
        g = F.gen()
        h = g.minpoly()
        a2, a1 = list(A.polynomial())
        b2, b1 = list(B.polynomial())
        p2, p1 = list(h)[:2]
        roots = dcpsolver_pari4(p, p1, p2, a1, a2, b1, b2, scalar)
        for r in roots:
            root = F(int(r[0])) * g + F(int(r[1]))
            try:
                P = curve.lift_x(root)
                assert P[0] + (scalar * P)[0] + 2 == 0
            except (ValueError, AssertionError) as e:
                P = None
                continue
    else:
        raise Exception("Base fields of degree >2 not implemented")
    return P


def solve_multi_dcp_pari(scalar, guess1, params):
    """
    Solves DCP for multi :scalar: and polynomial x1+x2+2
    The polynomial is hardcoded in the pari implementation
    """
    curve, lam, beta = params["curve"], params["lam"], params["beta"]
    A, B = curve.a4(), curve.a6()
    F = curve.base_field()
    p = F.characteristic()
    P = None
    k1, k2 = scalar
    if F.degree() == 1:
        roots = multidcpsolver_pari2(p, A, B, beta, k1, k2, guess1)
        for root in [F(int(r)) for r in roots]:
            try:
                # print(root,"root")
                P = curve.lift_x(root)
                # print(P)
                if (beta ** (2 * (guess1 != 0))) * (P[0]) + ((k1 + k2 * lam) * P)[
                    0
                ] + 2 != 0:
                    P = -P
                assert (beta ** (2 * (guess1 != 0))) * (P[0]) + ((k1 + k2 * lam) * P)[
                    0
                ] + 2 == 0
                break
            except (ValueError, AssertionError) as e:
                P = None
                continue
    return P


def solve_multi_dcp(scalar, guess, params):
    """
    Solves DCP for multi :scalar: and polynomial x1+x2+2
    The polynomial is hardcoded in the sage implementation below
    """
    m1, m2 = scalar
    curve, lam, beta = params["curve"], params["lam"], params["beta"]

    P = None
    if m2 == 0:
        mx = curve.multiplication_by_m(m1)[0]
    else:
        m1map = curve.multiplication_by_m(m1)
        m2map = curve.multiplication_by_m(m2)
        m2map = beta * m2map[0], m2map[1]
        mx, _ = msm.symbolic_addition(m1map, m2map, curve.a4(), curve.a6())
    numx = mx.numerator().univariate_polynomial()
    denx = mx.denominator().univariate_polynomial()
    x = denx.parent().gen()

    if guess[1] == 0:
        # tablepoint = (x,y) or (x,-y)
        f = x * denx + numx + 2 * denx
    else:
        # tablepoint = (beta^2x,-y) or (beta^2x,y)
        f = beta**2 * x * denx + numx + 2 * denx

    for r, _ in f.roots():
        try:
            P = curve.lift_x(r)
            break
        except ValueError as e:
            continue

    return P


def solve_multi_dcp_x(curve, beta, guess, scalar):
    """
    Solves DCP for multi :scalar: and polynomial x1+x2+2
    The polynomial is hardcoded in the Sage implementation below
    """
    m1, m2 = scalar
    m1x = curve.multiplication_by_m(m1, True)
    m2x = curve.multiplication_by_m(m2, True)
    x = m1x.parent().gens()[0]
    if guess[1] == 0:
        X1 = -2 - x
    else:
        X1 = -2 - beta**2 * x
    X2 = m1x
    X3 = beta * m2x
    A, B = curve.a4(), curve.a6()
    f = (
        (X1 - X2) ** 2 * X3**2
        - 2 * ((X1 + X2) * (X1 * X2 + A) + 2 * B) * X3
        + (X1 * X2 - A) ** 2
        - 4 * B * (X1 + X2)
    )
    fnum = f.numerator().univariate_polynomial()
    for r, _ in fnum.roots():
        try:
            P = curve.lift_x(r)
            betaP = curve(beta * P[0], P[1])
            Q = m1 * P + m2 * betaP
            assert P[0] + Q[0] + 2 == 0
            break
        except (ValueError, AssertionError):
            continue

    return P
