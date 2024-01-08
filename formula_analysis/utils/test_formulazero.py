from test_formula import give_me_curve
import formula as fm
import operation as op
import formulazero as fm0
from sage.all import PolynomialRing
import json


class NoTestVector(Exception):
    pass


def any_root_multi(f, field):
    for i, g in enumerate(f.parent().gens()):
        if f.variables()[0] == g:
            break
    counter = 0
    while True:
        counter += 1
        if counter == field.order():
            raise NoTestVector()
        s = [field.random_element() for _ in range(len(f.parent().gens()))]
        s[i] = g
        fs = f(*s)
        if len(fs.parent().gens()) > 1:
            fs = fs.univariate_polynomial()
        roots = fs.roots()
        if roots != []:
            solution = s[:i] + [roots[0][0]] + s[i + 1 :]
            break
    return dict(zip(f.parent().gens_dict().keys(), solution))


def rpa_shortw_add(formula, point):
    P = point["X1"], point["Y1"]
    Q = point["X2"], point["Y2"]
    try:
        R = op.Add(formula)(P, Q)
        assert P[0] * Q[0] * P[1] * Q[1] * R[0] * R[1] == 0
    except ZeroDivisionError:
        pass


def rpa_shortw_dbl(formula, point):
    if formula.coordinates == "xz":
        point["Y1"] = 1
    P = point["X1"], point["Y1"]
    try:
        R = op.Dbl(formula)(P)
        if formula.coordinates == "xz":
            assert P[0] * P[1] * R == 0
        else:
            assert P[0] * P[1] * R[0] * R[1] == 0
    except ZeroDivisionError:
        pass


def rpa_shortw_dadd(formula, point):
    P_Q = (point["X1"],)
    P = (point["X2"],)
    Q = (point["X3"],)
    try:
        R = (op.Dadd(formula)(P, Q, P_Q),)
        assert P[0] * Q[0] * R[0] * P_Q[0] == 0
    except ZeroDivisionError:
        pass


def rpa_shortw_ladd(formula, point):
    P_Q = (point["X1"],)
    P = (point["X2"],)
    Q = (point["X3"],)
    try:
        R = op.Ladd(formula)(P, Q, P_Q)
        assert P[0] * Q[0] * R[0] * R[1] * P_Q[0] == 0
    except ZeroDivisionError:
        pass


def rpa_twisted_dbl(formula, point):
    P = point["X1"], point["Y1"]
    try:
        R = op.Dbl(formula)(P)
        assert P[0] * P[1] * R[0] * R[1] == 0
    except ZeroDivisionError:
        pass


def rpa_twisted_add(formula, point):
    P = point["X1"], point["Y1"]
    Q = point["X2"], point["Y2"]
    try:
        R = op.Add(formula)(P, Q)
        assert P[0] * Q[0] * P[1] * Q[1] * R[0] * R[1] == 0
    except ZeroDivisionError:
        pass


def rpa_edwards_add(formula, point):
    P = point["X1"], point["Y1"]
    Q = point["X2"], point["Y2"]
    try:
        R = op.Add(formula)(P, Q)
        assert P[0] * Q[0] * P[1] * Q[1] * R[0] * R[1] == 0
    except ZeroDivisionError:
        pass


def rpa_edwards_dadd(formula, point):
    P_Q = 0, point["Y1"]
    P = 0, point["Y2"]
    Q = 0, point["Y3"]
    try:
        R = (op.Dadd(formula)(P, Q, P_Q),)
        assert P[0] * Q[0] * R[0] * P_Q[0] == 0
    except ZeroDivisionError:
        pass


def rpa_edwards_ladd(formula, point):
    P_Q = 0, point["Y1"]
    P = 0, point["Y2"]
    Q = 0, point["Y3"]
    try:
        R = op.Ladd(formula)(P, Q, P_Q)
        assert P[0] * Q[0] * R[0] * R[1] * P_Q[0] == 0
    except ZeroDivisionError:
        pass


def rpa_edwards_dbl(formula, point):
    if "yz" in formula.coordinates:
        point["X1"] = 1
    P = point["X1"], point["Y1"]
    try:
        R = op.Dbl(formula)(P)
        if "yz" in formula.coordinates:
            assert P[0] * P[1] * R[0] == 0
        else:
            assert P[0] * P[1] * R[0] * R[1] == 0
    except ZeroDivisionError:
        pass


def rpa_montgom_dadd(formula, point):
    P_Q = (point["X1"],)
    P = (point["X2"],)
    Q = (point["X3"],)
    try:
        R = (op.Dadd(formula)(P, Q, P_Q),)
        assert P[0] * Q[0] * R[0] * P_Q[0] == 0
    except ZeroDivisionError:
        pass


def rpa_montgom_ladd(formula, point):
    P_Q = (point["X1"],)
    P = (point["X2"],)
    Q = (point["X3"],)
    try:
        R = op.Ladd(formula)(P, Q, P_Q)
        assert P[0] * Q[0] * R[0] * R[1] * P_Q[0] == 0
    except ZeroDivisionError:
        pass


def rpa_montgom_dbl(formula, point):
    P = (point["X1"],)
    try:
        R = op.Dbl(formula)(P)
        assert P[0] * R == 0
    except ZeroDivisionError:
        pass


def rpa_test(formulas, test_func):
    for name, formula in formulas.items():
        formula.flag_homogeneity()
        formula.flag_multiplication_intermediate_values()
        formula.to_affine()
        field, curve = give_me_curve(formula)
        field_ring = PolynomialRing(field, formula.ring.variable_names())
        formula.set_p(field.order())
        formula.set_curve(curve)
        formula0 = fm0.FormulaZero.from_formula(formula)
        for _, polynomial0 in formula0.coordinate_zvp():
            polynomial = field_ring(polynomial0)
            absub = {
                polynomial.parent().gens_dict()[par]: curve[par]
                for par in formula0.coefficients
            }
            new_ring = PolynomialRing(field, formula.variables)
            try:
                polynomial = new_ring(polynomial.subs(absub))
            except:
                print(
                    polynomial,
                    new_ring,
                    absub,
                    polynomial.parent(),
                    formula0.coefficients,
                )
                raise Exception
            try:
                root = any_root_multi(polynomial, field)
            except NoTestVector:
                print(name, polynomial0, "not tested")
                continue
            test_func(formula, root)


def easy_zpa_operation(path, form_class):

    formulas = fm.load_formulas(path)
    for name, formula in formulas.items():
        formula.flag_homogeneity()
        formula.flag_multiplication_intermediate_values()
        formula.to_affine()
        field, curve = give_me_curve(formula)
        formula.set_p(field.order())
        formula.set_curve(curve)
        formula0 = fm0.FormulaZero.from_formula(formula)
        for poly in formula0.single_point_factors:
            assert len(poly.variables()) <= 2
            if len(poly.variables()) == 2:
                X, Y = map(str, poly.variables())
                assert X[1] == Y[1]
        for poly in formula0.multi_point_factors:
            assert len(poly.variables()) > 1
            if len(poly.variables()) == 2:
                X, Y = map(str, poly.variables())
                assert X[1] != Y[1]


def test_easy_zpa():
    easy_zpa_operation("unrolling/unrolled/edwards_add.json", fm.EdwardsFormula)
    easy_zpa_operation("unrolling/unrolled/edwards_dadd.json", fm.EdwardsFormula)
    easy_zpa_operation("unrolling/unrolled/edwards_ladd.json", fm.EdwardsFormula)

    easy_zpa_operation("unrolling/unrolled/shortw_add.json", fm.ShortWFormula)
    easy_zpa_operation("unrolling/unrolled/shortw_dadd.json", fm.ShortWFormula)
    easy_zpa_operation("unrolling/unrolled/shortw_ladd.json", fm.ShortWFormula)

    easy_zpa_operation("unrolling/unrolled/twisted_add.json", fm.TwistedFormula)

    easy_zpa_operation("unrolling/unrolled/montgom_dadd.json", fm.MontgomFormula)
    easy_zpa_operation("unrolling/unrolled/montgom_ladd.json", fm.MontgomFormula)


def test_rpa():
    rpa_test(fm.load_formulas("unrolling/unrolled/shortw_add.json"), rpa_shortw_add)
    rpa_test(fm.load_formulas("unrolling/unrolled/shortw_dbl.json"), rpa_shortw_dbl)
    rpa_test(fm.load_formulas("unrolling/unrolled/shortw_dadd.json"), rpa_shortw_dadd)
    rpa_test(fm.load_formulas("unrolling/unrolled/shortw_ladd.json"), rpa_shortw_ladd)

    rpa_test(fm.load_formulas("unrolling/unrolled/twisted_add.json"), rpa_twisted_add)
    rpa_test(fm.load_formulas("unrolling/unrolled/twisted_dbl.json"), rpa_twisted_dbl)

    rpa_test(fm.load_formulas("unrolling/unrolled/edwards_add.json"), rpa_edwards_add)
    rpa_test(fm.load_formulas("unrolling/unrolled/edwards_dadd.json"), rpa_edwards_dadd)
    rpa_test(fm.load_formulas("unrolling/unrolled/edwards_dbl.json"), rpa_edwards_dbl)
    rpa_test(fm.load_formulas("unrolling/unrolled/edwards_ladd.json"), rpa_edwards_ladd)

    rpa_test(fm.load_formulas("unrolling/unrolled/montgom_ladd.json"), rpa_montgom_ladd)
    rpa_test(fm.load_formulas("unrolling/unrolled/montgom_dbl.json"), rpa_montgom_dbl)
    rpa_test(fm.load_formulas("unrolling/unrolled/montgom_dadd.json"), rpa_montgom_dadd)
