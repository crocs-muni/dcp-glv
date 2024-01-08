from sage.all import GF, PolynomialRing, ZZ, cartesian_product
import json
from formula import ShortWFormula, MontgomFormula, EdwardsFormula, TwistedFormula
import utils.formula as fm
import utils.formulazero as fm0
import re
import utils.operation as op

p = ZZ(31)
p2 = ZZ(37)
field = GF(p)
field2 = GF(p2)
curves = {
    "twisted-1": {
        "a": field(-1),
        "d": field(3),
        "P": (field(24), field(3)),
        "Q": (field(16), field(6)),
        "P+Q": (field(28), field(23)),
        "2P": (field(15), field(25)),
        "inf": (field(0), field(1)),
    },
    "twisted-1,2": {
        "a": field2(-1),
        "d": field2(3),
        "P": (field2(30), field2(5)),
        "Q": (field2(28), field2(9)),
        "P+Q": (field2(27), field2(18)),
        "2P": (field2(6), field2(0)),
        "inf": (field2(0), field2(1)),
        "y0": (field2(6), field2(0)),
    },
    "edwards": {
        "c": field(2),
        "d": field(3),
        "P": (field(21), field(3)),
        "2P": (field(8), field(15)),
        "Q": (field(17), field(6)),
        "P+Q": (field(17), field(25)),
        "inf": (field(0), field(1)),
        "y0": (field(29), field(0)),
    },
    "edwards-yz": {
        "r": field(5),
        "c": field(1),
        "d": field(25),
        "P": (field(27), field(2)),
        "Q": (field(20), field(7)),
        "P+Q": (field(29), field(4)),
        "inf": (field(0), field(1)),
        "P-Q": (field(7), field(11)),
        "2P": (field(24), field(11)),
        "y0": (field(30), field(0)),
    },
    "shortw": {
        "a": field(24),
        "b": field(10),
        "P": (field(27), field(6)),
        "Q": (field(21), field(17)),
        "P+Q": (field(1), field(29), field(1)),
        "P-Q": (field(2), field(2)),
        "2P": (field(13), field(16)),
    },
    "shortw-1": {
        "a": field(-1),
        "b": field(3),
        "P": (field(14), field(25)),
        "Q": (field(4), field(30)),
        "P+Q": (field(21), field(25), field(1)),
        "2P": (field(28), field(14)),
    },
    "shortw-1,2": {
        "a": field(-1),
        "b": field(4),
        "P": (field(30), field(2)),
        "Q": (field(6), field(20)),
        "P+Q": (field(13), field(24), field(1)),
        "2P": (field(10), field(8)),
        "x0": (field(0), field(2)),
        "y0": (field(5), field(0)),
        "y0-P": (field(3), field(11)),
        "x0-P": (field(17), field(23)),
    },
    "shortw-3": {
        "a": field(-3),
        "b": field(3),
        "P": (field(21), field(5)),
        "Q": (field(22), field(18)),
        "P+Q": (field(2), field(25), field(1)),
        "2P": (field(22), field(18)),
    },
    "shortw-3,2": {
        "a": field(-3),
        "b": field(4),
        "P": (field(4), field(26)),
        "Q": (field(18), field(4)),
        "P+Q": (field(14), field(3), field(1)),
        "2P": (field(20), field(15)),
        "x0": (field(0), field(2)),
        "y0": (field(11), field(0)),
        "y0-P": (field(26), field(24)),
        "x0-P": (field(14), field(3)),
    },
    "shortw-0": {
        "a": field(0),
        "b": field(3),
        "P": (field(25), field(29)),
        "2P": (field(28), field(21)),
    },
    "shortw-0,2": {
        "a": field(0),
        "b": field(4),
        "P": (field(16), field(15)),
        "P-Q": (field(0), field(29), field(1)),
        "P+Q": (field(22), field(9), field(1)),
        "2P": (field(17), field(9)),
        "Q": (field(9), field(12)),
        "x0": (field(0), field(2)),
        "y0": (field(3), field(0)),
        "y0-P": (field(17), field(9)),
        "x0-P": (field(24), field(8)),
    },
    "shortw-a0": {
        "a": field(3),
        "b": field(0),
        "Q": (field(10), field(10)),
        "P+Q": (field(12), field(20)),
        "P": (field(17), field(29)),
        "2P": (field(5), field(27)),
        "y0": (field(0), field(0)),
        "x0": (field(0), field(0)),
        "y0-P": (field(2), field(18)),
        "x0-P": (field(2), field(18)),
    },
    "montgom": {
        "a": field(2),
        "b": field(1),
        "P": (field(2), field(7)),
        "Q": (field(5), field(5)),
        "P+Q": (field(19), field(25)),
        "P-Q": (field(7), field(13)),
        "2P": (field(4), field(21)),
    },
}
p3 = ZZ(100003)
field3 = GF(p3)
p4 = ZZ(107)
field4 = GF(p4)
curves_rpa = {
    "shortw-3": {
        "a": field3(100000),
        "b": field3(4),
        "P": (field3(44862), field3(88218)),
        "x0": (field3(0), field3(2)),
        "y0": (field3(13803), field3(0)),
        "x0-P": (field3(97407), field3(7657)),
        "y0-P": (field3(65540), field3(52830)),
        "x0-2P": (field3(92857), field3(12854)),
        "x0/2": (field3(13509), field3(78751)),
        "x0/2-P": (field3(11219), field3(60224)),
        "y0/2": (field3(16551), field3(18357)),
    },
    "shortw-1": {
        "a": field3(100002),
        "b": field3(4),
        "P": (field3(55494), field3(81655)),
        "x0": (field3(0), field3(2)),
        "y0": (field3(41203), field3(0)),
        "x0-P": (field3(64037), field3(83635)),
        "y0-P": (field3(34199), field3(57024)),
        "y0/2": (field3(79777), field3(14126)),
        "x0/2": (field3(41203), field3(0)),
    },
    "shortw-0": {
        "a": field3(0),
        "b": field3(64),
        "P": (field3(2596), field3(35083)),
        "x0": (field3(0), field3(8)),
        "y0": (field3(28484), field3(0)),
        "x0-P": (field3(90490), field3(50182)),
        "y0-P": (field3(13390), field3(60359)),
    },
    "shortw-a0": {
        "a": field3(3),
        "b": field3(0),
        "P": (field3(37050), field3(48553)),
        "x0": (field3(0), field3(0)),
        "y0": (field3(0), field3(0)),
        "x0-P": (field3(41564), field3(59694)),
        "y0-P": (field3(41564), field3(59694)),
    },
    "shortw-0,2": {
        "a": field4(0),
        "b": field4(4),
        "y0": (field4(71), field4(0)),
        "y0/2": (field4(77), field4(17)),
        "x0": (field4(0), field4(2)),
        "x0/2": (field4(72), field4(6)),
    },
    "shortw-a0,2": {
        "a": field3(1),
        "b": field3(0),
        "y0": (field3(0), field3(0)),
        "y0/2": (field3(100002), field3(44974)),
        "x0": (field3(0), field3(0)),
        "x0/2": (field3(100002), field3(44974)),
    },
}


def give_me_curve(formula):
    """Finds an appropriate curve for the given formula"""
    if formula.form == "edwards":
        if formula.operation in ["ladd", "dadd"]:
            return field, curves["edwards-yz"]
        if formula.coordinates.startswith("yz"):
            return field, curves["edwards-yz"]
        return field, curves["edwards"]
    if formula.form == "twisted":
        return field, curves["twisted-1"]
    if formula.form == "shortw":
        if formula.coordinates.endswith("w12-0"):
            return field, curves["shortw-a0"]
        elif formula.coordinates.endswith("-0"):
            return field, curves["shortw-0,2"]
        elif formula.coordinates.endswith("-1"):
            return field, curves["shortw-1"]
        elif formula.coordinates.endswith("-3"):
            return field, curves["shortw-3"]
        else:
            return field, curves["shortw"]
    assert formula.form == "montgom"
    return field, curves["montgom"]


def twisted_to_affine(formula):
    """Testing function for the reduction of twisted add and dbl"""
    formula.to_affine()
    for c in formula.coefficients:
        assert c in ["a", "d"]
    for c in formula.variables:
        assert c[0] in ["X", "Y"]
        if formula.operation == "dbl":
            assert c[1] in ["1"]
        else:
            assert formula.operation == "add"
            assert c[1] in ["1", "2"]
    for row in formula.formula:
        assert re.match(r"^[XY0-9ad\*\-\+\s\^]*$", str(row.numerator))
        assert re.match(r"^[XY0-9ad\*\-\+\s\^]*$", str(row.denominator))


def twisted_operation(formula):
    formula.to_affine()
    field, curve = give_me_curve(formula)
    formula.set_p(field.order())
    formula.set_a(curve["a"])
    formula.set_d(curve["d"])
    if formula.operation == "add":
        comp_add = op.Add(formula)(curve["P"], curve["Q"])
        assert curve["P+Q"][:2] == comp_add
    elif formula.operation == "dbl":
        comp_dbl = op.Dbl(formula)(curve["P"])
        assert curve["2P"][:2] == comp_dbl, {formula.name}
    else:
        assert False, formula.operation


def twisted_operation_all(path):
    with open(path) as f:
        all_formulas = json.load(f)
    for _, formulas in all_formulas.items():
        for _, formula_dict in formulas.items():
            formula = TwistedFormula(formula_dict)
            twisted_to_affine(formula)
            formula = TwistedFormula(formula_dict)
            twisted_operation(formula)


def edwards_to_affine_operation(path):

    with open(path) as f:
        all_formulas = json.load(f)
    for _, formulas in all_formulas.items():
        for _, formula_dict in formulas.items():
            formula = EdwardsFormula(formula_dict)
            formula.to_affine()
            for c in formula.coefficients:
                assert formula.coordinates.startswith("yz") or c in ["c", "d"]
                assert not formula.coordinates.startswith("yz") or c in ["r"]
            for c in formula.variables:
                assert c[0] in ["X", "Y"]
                if formula.operation == "dbl":
                    assert c[1] in ["1"]
                elif formula.operation == "add":
                    assert c[1] in ["1", "2"]
                elif formula.operation == "dadd":
                    assert c[1] in ["1", "2", "3"]
                else:
                    assert formula.operation == "ladd"
                    assert c[1] in ["1", "2", "3", "4"]
            for row in formula.formula:
                assert re.match(
                    r"^[XY0-9cdr\*\-\+\s\^]*$", str(row.numerator)
                ), row.numerator
                assert re.match(r"^[XY0-9cdr\*\-\+\s\^]*$", str(row.denominator))


def edwards_operation(path):
    with open(path) as f:
        all_formulas = json.load(f)
    for _, formulas in all_formulas.items():
        for name, formula_dict in formulas.items():
            if name in [
                "yzsquared:mladd-2006-g",
                "yzsquared:mladd-2006-g-2",
                "yzsquared:ladd-2006-g",
                "yzsquared:ladd-2006-g-2",
            ]:
                continue
            formula = EdwardsFormula(formula_dict)
            field, curve = give_me_curve(formula)
            formula.to_affine()
            formula.set_p(field.order())

            if formula.operation == "dadd":
                formula.set_r(curve["r"])
                comp_add = op.Dadd(formula)(curve["P"], curve["Q"], curve["P-Q"])
                assert curve["P+Q"][1] == comp_add, name
            elif formula.operation == "ladd":
                formula.set_r(curve["r"])
                comp_add = op.Ladd(formula)(curve["P"], curve["Q"], curve["P-Q"])
                assert curve["P+Q"][1] == comp_add[0], name
            elif formula.operation == "add":
                formula.set_c(curve["c"])
                formula.set_d(curve["d"])
                comp_add = op.Add(formula)(curve["P"], curve["Q"])
                assert curve["P+Q"][:2] == comp_add
            elif formula.operation == "dbl":
                if "r" in curve:
                    formula.set_r(curve["r"])
                formula.set_c(curve["c"])
                formula.set_d(curve["d"])
                if "r" in curve:
                    comp_dbl = op.Dbl(formula)(curve["P"])
                else:
                    comp_dbl = op.Dbl(formula)(curve["P"])
                if formula.coordinates in ["projective", "inverted"]:
                    assert curve["2P"][:2] == comp_dbl, {name}
                else:
                    assert curve["2P"][1] == comp_dbl[0]
            else:
                assert False, formula.operation


def shortw_to_affine_operation(path):

    with open(path) as f:
        all_formulas = json.load(f)
    operation = path.split("_")[-1].split(".")[0]
    for _, formulas in all_formulas.items():
        for formula_name, formula_dict in formulas.items():
            print(formula_name)
            formula = ShortWFormula(formula_dict)
            formula.to_affine()
            for c in formula.coefficients:
                assert c in ["a", "b"]
            for c in formula.variables:
                assert c[0] in ["X", "Y"]
                assert c[0] in ["X", "Y"]
                if operation == "dbl":
                    assert c[1] in ["1"]
                elif operation == "add":
                    assert c[1] in ["1", "2"]
                elif operation == "dadd":
                    assert c[1] in ["1", "2", "3"]
                else:
                    assert operation == "ladd"
                    assert c[1] in ["1", "2", "3", "4"]
            for row in formula.formula:
                assert re.match(
                    r"^[XY0-9ab\*\-\+\s\^]*$", str(row.numerator)
                ), row.numerator
                assert re.match(r"^[XY0-9ab\*\-\+\s\^]*$", str(row.denominator))


def shortw_operation(path):
    with open(path) as f:
        all_formulas = json.load(f)
    for name, formulas in all_formulas.items():
        for name2, formula_dict in formulas.items():
            formula = ShortWFormula(formula_dict)
            formula.to_affine()
            field, curve = give_me_curve(formula)
            formula.set_p(field.order())
            formula.set_a(curve["a"])
            formula.set_b(curve["b"])
            if formula.operation == "add":
                comp_add = op.Add(formula)(curve["P"], curve["Q"])
                assert curve["P+Q"][:2] == comp_add, name2
            elif formula.operation == "dadd":
                comp_add = op.Dadd(formula)(curve["P"], curve["Q"], curve["P-Q"])
                assert curve["P+Q"][0] == comp_add
            elif formula.operation == "ladd":
                comp_add = op.Ladd(formula)(curve["P"], curve["Q"], curve["P-Q"])
                assert curve["P+Q"][0] == comp_add[0]
                assert curve["2P"][0] == comp_add[1], name2
            elif formula.operation == "dbl":
                comp_dbl = op.Dbl(formula)(curve["P"])
                if not isinstance(comp_dbl, tuple):
                    assert curve["2P"][0] == comp_dbl
                else:
                    assert curve["2P"][: len(comp_dbl)] == comp_dbl
            else:
                assert False


def montgom_to_affine_operation(path):
    with open(path) as f:
        all_formulas = json.load(f)
    operation = path.split("_")[-1].split(".")[0]
    for _, formulas in all_formulas.items():
        for _, formula_dict in formulas.items():
            formula = MontgomFormula(formula_dict)
            formula.to_affine()
            for c in formula.coefficients:
                assert c in ["a", "b"]
            for c in formula.variables:
                assert c[0] in ["X", "Y"]
                assert c[0] in ["X", "Y"]
                if operation == "dbl":
                    assert c[1] in ["1"]
                elif operation == "dadd":
                    assert c[1] in ["1", "2", "3"]
                else:
                    assert operation == "ladd"
                    assert c[1] in ["1", "2", "3", "4"]
            for row in formula.formula:
                assert re.match(
                    r"^[XY0-9ab\*\-\+\s\^]*$", str(row.numerator)
                ), row.numerator
                assert re.match(r"^[XY0-9ab\*\-\+\s\^]*$", str(row.denominator))


def montgom_operation(path="unrolling/unrolled/montgom_dadd.json"):
    with open(path) as f:
        all_formulas = json.load(f)
    for _, formulas in all_formulas.items():
        for _, formula_dict in formulas.items():
            formula = MontgomFormula(formula_dict)
            formula.to_affine()
            field, curve = give_me_curve(formula)
            formula.set_p(field.order())
            formula.set_a(curve["a"])
            formula.set_b(curve["b"])
            if formula.operation == "dadd":
                comp_add = op.Dadd(formula)(curve["P"], curve["Q"], curve["P-Q"])
                assert curve["P+Q"][0] == comp_add
            elif formula.operation == "ladd":
                comp_add = op.Ladd(formula)(curve["P"], curve["Q"], curve["P-Q"])
                assert curve["P+Q"][0] == comp_add[0]
                assert curve["2P"][0] == comp_add[1]
            elif formula.operation == "dbl":
                comp_dbl = op.Dbl(formula)(curve["P"])
                assert curve["2P"][0] == comp_dbl
            else:
                assert False


def test_multiplicative_intermediates():
    shortw_add = fm.load_formulas(f"unrolling/unrolled/shortw_add.json")
    f = shortw_add["xyzz:mmadd-2008-s"]
    assert f.row_operation(0) == "-"
    assert f.row_operation(2) == "*"
    assert f.row_intermediate_names(0) == ["X2", "X1"]

    f = shortw_add["jacobian-0:madd-2004-hmv"]
    row = f.formula[4]
    assert f.intermediate_scope(row) == (4, 10)
    assert f.intermediate_operations(row) == set("*")

    f.flag_multiplication_intermediate_values()
    flags = [row.flag for row in f.formula]
    assert flags[:4] == [True, True, False, False]


def test_twisted_operations():
    twisted_operation_all("unrolling/unrolled/twisted_add.json")
    twisted_operation_all("unrolling/unrolled/twisted_dbl.json")


def test_edwards_to_affine():
    edwards_to_affine_operation("unrolling/unrolled/edwards_add.json")
    edwards_to_affine_operation("unrolling/unrolled/edwards_dbl.json")
    edwards_to_affine_operation("unrolling/unrolled/edwards_dadd.json")
    edwards_to_affine_operation("unrolling/unrolled/edwards_ladd.json")


def test_edwards_operations():
    edwards_operation(path="unrolling/unrolled/edwards_ladd.json")
    edwards_operation(path="unrolling/unrolled/edwards_dadd.json")
    edwards_operation(path="unrolling/unrolled/edwards_dbl.json")
    edwards_operation(path="unrolling/unrolled/edwards_add.json")


def test_shortw_to_affine():
    shortw_to_affine_operation("unrolling/unrolled/shortw_add.json")
    shortw_to_affine_operation("unrolling/unrolled/shortw_dbl.json")
    shortw_to_affine_operation("unrolling/unrolled/shortw_dadd.json")
    shortw_to_affine_operation("unrolling/unrolled/shortw_ladd.json")


def test_shortw_operations():
    shortw_operation("unrolling/unrolled/shortw_add.json")
    shortw_operation("unrolling/unrolled/shortw_dbl.json")
    shortw_operation("unrolling/unrolled/shortw_dadd.json")
    shortw_operation("unrolling/unrolled/shortw_ladd.json")


def test_montgom_to_affine():
    montgom_to_affine_operation("unrolling/unrolled/montgom_dbl.json")
    montgom_to_affine_operation("unrolling/unrolled/montgom_dadd.json")
    montgom_to_affine_operation("unrolling/unrolled/montgom_ladd.json")


def test_montgom_operations():
    montgom_operation("unrolling/unrolled/montgom_dbl.json")
    montgom_operation("unrolling/unrolled/montgom_dadd.json")
    montgom_operation("unrolling/unrolled/montgom_ladd.json")
