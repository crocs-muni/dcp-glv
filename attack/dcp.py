"""
DCP solver
"""
from sage.all import PolynomialRing, ZZ
from pari_tools.dcp_pari import (
    dcpsolver_pari,
    glvdcpsolver_pari,
    multidcpsolver_pari,
    glvdcpmultisolver_pari,
)
from random import randint
import time 


def solve_dcp(scalar, glv, registers, mul_lam=False):
    """
    Solves DCP for :scalar: and polynomial x1+x2+2
    The polynomial is hardcoded in the pari implementation
    Works for prime fields and quadratic fields
    """
    curve = glv.curve
    A, B = curve.a4(), curve.a6()
    F = curve.base_field()
    P = None
    assert scalar != 0
    lam, beta = (1, 1) if not mul_lam else (glv.lam, glv.beta)
    Vpolynomials = [
        prepare_simpleVpolynomial(polynomial, glv, registers, beta)
        for polynomial in registers.polynomials
    ]
    solutions = dcpsolver_pari(glv.p, A, B, scalar, lam, Vpolynomials, registers)
    for x, y in solutions:
        P = curve(F(x), F(y))
        Q = lam * scalar * P
        assert registers.is_zero(P, Q)
    return P


def solve_multi_dcp(scalar, scalar0, glv, registers, mul_lam=False):
    """
    Solves DCP for :scalar: and polynomial x1+x2+2
    The polynomial is hardcoded in the pari implementation
    Works for prime fields and quadratic fields
    """
    curve = glv.curve
    A, B = curve.a4(), curve.a6()
    F = curve.base_field()
    P = None
    assert scalar != 0
    assert scalar0 != 0
    if scalar0 == 1:
        return solve_dcp(scalar, glv, registers, mul_lam)
    lam, beta = (1, 1) if not mul_lam else (glv.lam, glv.beta)
    Vpolynomials = [
        prepare_simpleVpolynomial_multi(polynomial, glv, registers, beta)
        for polynomial in registers.polynomials
    ]
    solutions = multidcpsolver_pari(
        glv.p, A, B, scalar, scalar0, lam, Vpolynomials, registers
    )
    for x, y in solutions:
        P = curve(F(x), F(y))
        Q = lam * scalar * P
        P0 = scalar0 * P
        assert registers.is_zero(P0, Q)
    return P


def prepare_simpleVpolynomial(polynomial, glv, registers, beta=1):
    X1, Y1, X2, Y2, A, B = registers.all_gens

    # we assume that the polynomial is only in X1,X2,A,B
    assert set(polynomial.variables()).issubset(set([X1, X2, A, B]))
    a4, a6 = ZZ(glv.curve.a4()), ZZ(glv.curve.a6())

    x, x1, x2, n, d = PolynomialRing(ZZ, "x,x1,x2,n,d").fraction_field().gens()
    V = x.parent()(polynomial(x1, 1, x2, 1, a4, a6))
    V = V(x, x, beta * n / d, n, d).numerator()

    # clean out the tmp variables
    V = PolynomialRing(ZZ, "x,n,d")(V)
    return V


def prepare_simpleVpolynomial_multi(polynomial, glv, registers, beta=1):
    X1, Y1, X2, Y2, A, B = registers.all_gens

    # we assume that the polynomial is only in X1,X2,A,B
    assert set(polynomial.variables()).issubset(set([X1, X2, A, B]))
    a4, a6 = ZZ(glv.curve.a4()), ZZ(glv.curve.a6())

    n0, d0, x1, x2, n, d = PolynomialRing(ZZ, "n0,d0,x1,x2,n,d").fraction_field().gens()
    V = n0.parent()(polynomial(x1, 1, x2, 1, a4, a6))
    V = V(n0, d0, n0 / d0, beta * n / d, n, d).numerator()
    # clean out the tmp variables
    V = PolynomialRing(ZZ, "n0,d0,n,d")(V)
    return V


def prepare_Vpolynomial(polynomial, registers, glv):
    X1, Y1, X2, Y2, A, B = registers.all_gens

    # we assume that the polynomial is only in X1,X2,A,B
    assert set(polynomial.variables()).issubset(set([X1, X2, A, B]))
    polynomial = polynomial(X1, Y1, X2, Y2, ZZ(glv.curve.a4()), ZZ(glv.curve.a6()))

    # rename the variables to have nested rings
    Qring = PolynomialRing(ZZ, "XQ,XQp")
    XQ, XQp = Qring.gens()
    Pring = PolynomialRing(Qring, "XP")
    XP = Pring.gen()
    f = 0
    for e, c in polynomial.iterator_exp_coeff():
        f += c * (XP ** e[0]) * (XQ ** e[2])

    # Create symmetric polynomial
    revf = 0
    for monom in f.monomials():
        coef = f.monomial_coefficient(monom)
        revf += coef(XQp, XQ) * monom
    symf = revf * f

    # Express XQ,XQp using symmetric polynomials and substitute Seamev coefficients
    x, qe1, qe2, x1, x2, n1, d1, n2, d2 = (
        PolynomialRing(ZZ, "x,qe1,qe2,x1,x2,n1,d1,n2,d2").fraction_field().gens()
    )
    A = ZZ(glv.curve.a4())
    B = ZZ(glv.curve.a6())
    Qprod = XQ * XQp
    Qsum = XQ + XQp
    V = 0
    for monom in symf.monomials():
        coef = symf.monomial_coefficient(monom)
        new_monom = 1
        new_coef = fundamental_sym_thm(qe1, qe2, Qsum, Qprod, XQ, XQp, coef)
        while XP.divides(monom):
            monom //= XP
            new_monom *= x
        V += new_monom * new_coef
    linear = -2 * ((x1 + x2) * (x1 * x2 + A) + 2 * B) / ((x1 - x2) ** 2)
    constant = ((x1 * x2 - A) ** 2 - 4 * B * (x1 + x2)) / ((x1 - x2) ** 2)
    V = V(x, -linear, constant, x1, x2, 1, 1, 1, 1).numerator()
    V = V(x, 1, 1, n1 / d1, glv.beta * n2 / d2, 1, 1, 1, 1).numerator()  # Adding beta

    # clean out the tmp variables
    V = PolynomialRing(ZZ, "x,n1,d1,n2,d2")(V)

    return V




def prepare_Vpolynomial_multi(polynomial, registers, glv):
    X1, Y1, X2, Y2, A, B = registers.all_gens

    # we assume that the polynomial is only in X1,X2,A,B
    assert set(polynomial.variables()).issubset(set([X1, X2, A, B]))
    polynomial = polynomial(X1, Y1, X2, Y2, ZZ(glv.curve.a4()), ZZ(glv.curve.a6()))

    # rename the variables to have nested rings
    Qring = PolynomialRing(ZZ, "XQ,XQp")
    XQ, XQp = Qring.gens()
    Pring = PolynomialRing(Qring, "XP")
    XP = Pring.gen()
    f = 0
    for e, c in polynomial.iterator_exp_coeff():
        f += c * (XP ** e[0]) * (XQ ** e[2])

    # Create symmetric polynomial
    revf = 0
    for monom in f.monomials():
        coef = f.monomial_coefficient(monom)
        revf += coef(XQp, XQ) * monom
    symf = revf * f

    # Express XQ,XQp using symmetric polynomials and substitute Seamev coefficients
    x, n0, d0, qe1, qe2, x1, x2, n1, d1, n2, d2 = (
        PolynomialRing(ZZ, "x,n0,d0,qe1,qe2,x1,x2,n1,d1,n2,d2").fraction_field().gens()
    )
    A = ZZ(glv.curve.a4())
    B = ZZ(glv.curve.a6())
    Qprod = XQ * XQp
    Qsum = XQ + XQp
    V = 0
    for monom in symf.monomials():
        coef = symf.monomial_coefficient(monom)
        new_monom = 1
        new_coef = fundamental_sym_thm(qe1, qe2, Qsum, Qprod, XQ, XQp, coef)
        while XP.divides(monom):
            monom //= XP
            new_monom *= x
        V += new_monom * new_coef
    linear = -2 * ((x1 + x2) * (x1 * x2 + A) + 2 * B) / ((x1 - x2) ** 2)
    constant = ((x1 * x2 - A) ** 2 - 4 * B * (x1 + x2)) / ((x1 - x2) ** 2)
    V = V(x, n0, d0, -linear, constant, x1, x2, 1, 1, 1, 1).numerator()
    V = V(
        n0 / d0, n0, d0, 1, 1, n1 / d1, glv.beta * n2 / d2, 1, 1, 1, 1
    ).numerator()  # Adding beta

    # clean out the tmp variables
    V = PolynomialRing(ZZ, "n0,d0,n1,d1,n2,d2")(V)
    return V


def fundamental_sym_thm(e1, e2, el1, el2, X2, X3, poly):
    elf = 0
    while poly != 0:
        lead_exp = poly.exponents()[0]
        lead_coef = ZZ(poly.coefficient(X2 ** lead_exp[0] * X3 ** lead_exp[1]))
        elf += lead_coef * e2 ** lead_exp[1] * e1 ** (lead_exp[0] - lead_exp[1])
        poly -= lead_coef * el2 ** lead_exp[1] * el1 ** (lead_exp[0] - lead_exp[1])
    return elf


def solve_glv_dcp_pari(glv_scalar, glv, registers):
    """
    Solves DCP for multi :scalar: and polynomial in registers
    """
    A, B = glv.curve.a4(), glv.curve.a6()
    P = None
    k1, k2 = glv_scalar
    if k2 == 0:
        assert k1 != 0
        return solve_dcp(k1, glv, registers)
    if k1 == 0:
        assert k2 != 0
        return solve_dcp(k2, glv, registers, True)

    Vpolynomials = [
        prepare_Vpolynomial(polynomial, registers, glv)
        for polynomial in registers.polynomials
    ]
    solutions = glvdcpsolver_pari(glv.p, A, B, k1, glv.lam, k2, Vpolynomials, registers)
    for x, y in solutions:
        P = glv.curve(glv.field(x), glv.field(y))
        Q = (k1 + k2 * glv.lam) * P
        assert registers.is_zero(P, Q)
        return P
    return P


def solve_glv_multi_dcp_pari(glv_scalar, scalar0, glv, registers):
    # print("solving", glv_scalar,scalar0)
    A, B = glv.curve.a4(), glv.curve.a6()
    P = None
    k1, k2 = glv_scalar
    assert scalar0 != 0
    if scalar0 == 1:
        return solve_glv_dcp_pari(glv_scalar, glv, registers)
    if k2 == 0:
        assert k1 != 0
        return solve_multi_dcp(k1, scalar0, glv, registers)
    if k1 == 0:
        assert k2 != 0
        return solve_multi_dcp(k2, scalar0, glv, registers, True)

    Vpolynomials = [
        prepare_Vpolynomial_multi(polynomial, registers, glv)
        for polynomial in registers.polynomials
    ]
    solutions = glvdcpmultisolver_pari(
        glv.p, A, B, scalar0, k1, glv.lam, k2, Vpolynomials, registers
    )
    for x, y in solutions:
        P = glv.curve(glv.field(x), glv.field(y))
        Q = (k1 + k2 * glv.lam) * P
        P0 = scalar0 * P
        assert registers.is_zero(P0, Q)
        return P
    return P




def dcp_experiment(params):
    result = {}
    scalar = params.secrets.k
    Vpolynomials = [prepare_simpleVpolynomial(polynomial, params.glv, params.registers)
        for polynomial in params.registers.polynomials
    ]
    time_zvp = time.time()
    P = solve_dcp(scalar, params.glv, params.registers)
    result["nguesses"] = int(0)
    result["recovered"] = int(P is not None)
    result["time_zvp"] = float(time.time() - time_zvp)
    if P is not None:
        result["point"] = [int(P[0]),int(P[1])]
    else:
        result["point"] = []
    result["Vpoly"] = [{"degree":V.degree(),"num_monom":len(V.monomials())} for V in Vpolynomials]
    return result

def glv_dcp_experiment(params):
    result = {}
    scalar0 = params.secrets.k0
    scalar1 = params.secrets.k1
    Vpolynomials = [prepare_Vpolynomial(polynomial, params.registers, params.glv)
        for polynomial in params.registers.polynomials
    ]
    time_zvp = time.time()
    P = solve_glv_dcp_pari((scalar0, scalar1), params.glv, params.registers)
    result["nguesses"] = int(0)
    result["recovered"] = int(P is not None)
    result["time_zvp"] = float(time.time() - time_zvp)
    if P is not None:
        result["point"] = [int(P[0]),int(P[1])]
    else:
        result["point"] = []
    result["Vpoly"] = [{"degree":V.degree(),"num_monom":len(V.monomials())} for V in Vpolynomials]
    return result



def interleaving_dcp_experiment(params):
    result = {}
    assert params.secrets.k0 is not None
    scalar0 = params.secrets.k0
    scalar1 = params.secrets.k1
    time_zvp = time.time()
    P = solve_multi_dcp(scalar0, scalar1, params.glv, params.registers, mul_lam=True)
    result["nguesses"] = int(0)
    result["recovered"] = int(P is not None)
    result["time_zvp"] = float(time.time() - time_zvp)
    if P is not None:
        result["point"] = [int(P[0]),int(P[1])]
    else:
        result["point"] = []
    return result  