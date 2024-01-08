import json
import re
from sage.all import PolynomialRing, ZZ, GF


class Formula:
    """
    Class to represent an elliptic curve operation formula in an unrolled form.
    Takes a dictionary on input which must have these keys:
        - 'variables' = list of strings ('X1','Y1',..) used as variables in the formula (not curve coefficients)
        - 'parameters'= list of strings representing curve coefficients (e.g., ['a','b'] for Weierstrass)
        - 'formula': = list of tuples (left, right) where left is a name of the intermediate value and right is an unrolled line.
        - 'original_formula' = original formula in pre-unrolled form
        - 'output_indexes' = dictionary representing output of the formula as {'X3':<index_of_line_in_formula_list>}
        - 'form','coordinates', 'operation' and 'name'
    """

    def __init__(self, dictionary):
        self.name = dictionary["name"]
        self.coordinates = dictionary["coordinates"]
        self.operation = dictionary["operation"]
        self.form = dictionary["form"]
        self.variables = dictionary["variables"]
        self.coefficients = dictionary["parameters"]
        self.original_formula = dictionary["original_formula"]
        self.output_indexes = dictionary["output"]
        self.transformations = dictionary["unified_transformation"]
        self.secondary_output_indexes = dictionary["secondary_output"]
        self.homogeneity_weights = dictionary["homogeneity_weights"]
        self.formula = []
        self.ring = PolynomialRing(ZZ, self.variables + self.coefficients)
        self.curve_field = None
        self.parse_formula(dictionary["formula"])
        self.transformation_methods = {
            "Z_to_1": self.put_Z_1,
            "invert": self.invert_coordinates,
        }

    def full_name(self):
        return f"{self.coordinates}:{self.name}"

    def __repr__(self):
        return f"{self.full_name()} formula for {self.operation} on {self.form} curve in variables {self.variables}"

    def parse_formula(self, formula_strings):
        self.formula = [
            Row(strings, self.ring, i) for i, strings in enumerate(formula_strings)
        ]

    def formula_print(self):
        """pretty print without denominators"""
        for row in self.formula:
            print(f"{row.left} = {row.numerator}")

    def formula_full_print(self):
        """pretty print with denominators"""
        for row in self.formula:
            print(f"{row.left} = {row.numerator/row.denominator}")

    def input_variables_grouped(self):
        input_points = {}
        pattern = re.compile(r"\d")
        for var in self.variables:
            number = int(re.search(pattern, var)[0])
            inpoint = input_points.setdefault(number, [])
            inpoint.append(var)
        return input_points

    def output(self):
        output = {}
        for left, index in self.output_indexes.items():
            row = self.formula[index]
            output[left] = row.numerator / row.denominator
        return output

    def secondary_output(self):
        output = {}
        for left, index in self.secondary_output_indexes.items():
            row = self.formula[index]
            output[left] = row.numerator / row.denominator
        return output

    def is_output(self, row):
        return (row.left[0], row.index) in self.output_indexes.items() or (
            row.left[0],
            row.index,
        ) in self.secondary_output_indexes.items()

    def put_Z_1(self):
        """puts all Z variables to 1"""
        new_formula = []
        self.variables = list(filter(lambda x: "Z" not in str(x), self.variables))
        ring = PolynomialRing(self.ring.base_ring(), self.variables + self.coefficients)
        substitution = {g: (1 if "Z" in str(g) else g) for g in self.ring.gens()}
        for row in self.formula:
            row.substitute(substitution, ring)
        self.ring = ring

    def invert_coordinates(self):
        """Invert coordinates (used for inverted-coordinates formulae), can be used only after put_Z_1"""
        substitution = {
            self.ring.gens_dict()[var]: 1 / self.ring.gens_dict()[var]
            for var in self.variables
        }
        for row in self.formula:
            row.substitute(substitution, self.ring)

    def to_affine(self):
        for transformation in self.transformations:
            self.transformation_methods[transformation]()

    def set_p(self, p):
        self.p = ZZ(p)
        self.curve_field = GF(p)
        self.ring = PolynomialRing(self.curve_field, self.variables + self.coefficients)
        for row in self.formula:
            row.change_ring(self.ring)

    def set_param(self, param, value):
        new_formula = []
        self.coefficients.remove(param)
        new_ring = PolynomialRing(
            self.ring.base_ring(), self.variables + self.coefficients
        )
        g = self.ring.gens_dict()[param]
        for row in self.formula:
            row.substitute({g: value}, new_ring)
        self.ring = new_ring

    def is_constant(self, polynomial):
        if polynomial.is_constant():
            return True
        return set(str(v) for v in polynomial.variables()).issubset(
            set(self.coefficients)
        )

    def flag_homogeneity(self):
        if "mmadd" in self.name:
            return
        homogenity_weights = self.homogeneity_weights
        input_variables_grouped = self.input_variables_grouped()
        if "zadd" in self.name:
            input_variables_grouped = {1: sum(input_variables_grouped.values(), [])}
        ring = self.ring
        for point_index, variables in input_variables_grouped.items():
            weighted_variables = [
                [ring(var), homogenity_weights[strip_numbers(var)]] for var in variables
            ]
            if "madd" in self.name and point_index == 2:
                continue
            for row in self.formula:
                row.flag &= row.is_homogenous(weighted_variables)

    def row_operation(self, index):
        for op in "-+*/":
            if op in self.original_formula[index][1]:
                return op

    def row_intermediate_names(self, index):
        operation = self.row_operation(index)
        right_str = self.original_formula[index][1]
        return right_str.split(operation)

    def intermediate_scope(self, input_row):
        assert self.formula[input_row.index] == input_row
        scope_beg = input_row.index
        scope_end = scope_beg
        for row in self.formula[scope_beg + 1 :]:
            scope_end += 1
            if row.left == input_row.left:
                break
        return scope_beg, scope_end

    def intermediate_operations(self, input_row):
        scope_beg, scope_end = self.intermediate_scope(input_row)
        operations = set()
        for row in self.formula[scope_beg + 1 : scope_end + 1]:
            if not input_row.left in self.row_intermediate_names(row.index):
                continue
            operations.add(self.row_operation(row.index))
        return operations

    def flag_multiplication_intermediate_values(self):
        for row in self.formula:
            if self.is_output(row):
                continue
            operations = self.intermediate_operations(row)
            row.flag &= "*" in operations


def strip_numbers(string: str):
    return "".join(filter(lambda z: not z.isdigit(), string))


class ShortWFormula(Formula):
    def __init__(self, dictionary):
        super().__init__(dictionary)
        for c in self.coefficients:
            assert c in ["a", "b"]
        self.a = None
        if self.coordinates.endswith("-0"):
            self.a = 0
        if self.coordinates.endswith("-3"):
            self.a = 3
        self.b = None
        if self.coordinates == "w12-0":
            self.b = 0
        self.transformation_methods["T_to_a"] = self.remove_T

    def remove_T(self):
        """remove T in modified Jacobian Coordinates"""
        for var in filter(lambda x: x.startswith("T"), self.variables):
            self.variables.remove(var)
            a, T = self.ring.gens_dict()["a"], self.ring.gens_dict()[var]
            self.ring = PolynomialRing(
                self.ring.base_ring(), self.variables + self.coefficients
            )
            for row in self.formula:
                row.substitute({T: a}, self.ring)

    def set_a(self, a):
        self.a = self.curve_field(a)
        if "a" in self.coefficients:
            self.set_param("a", self.a)

    def set_b(self, b):
        self.b = self.curve_field(b)
        if "b" in self.coefficients:
            self.set_param("b", self.b)

    def set_curve(self, params):
        if "p" in params:
            self.set_p(params["p"])
        if "a" in params:
            self.set_a(params["a"])
        if "b" in params:
            self.set_b(params["b"])


class EdwardsFormula(Formula):
    def __init__(self, dictionary):
        super().__init__(dictionary)
        for c in self.coefficients:
            assert c in ["c", "d", "r"]
        self.transformation_methods["xinvert"] = self.xinvert
        self.transformation_methods["xproject"] = self.xproject
        self.transformation_methods["remove_yz"] = self.remove_yz
        self.transformation_methods["remove_yz2"] = self.remove_yzsquared

    def xproject(self):
        """
        Transformation Z1=1,X2=1,Z2=1/x2, Y2=y2/x2 for projective:xmadd-2007-hcd"
        """
        new_formula = []
        variables = self.variables + ["X2"]
        ring_with_x = PolynomialRing(
            self.ring.base_ring(), variables + self.coefficients
        )
        variables = list(filter(lambda x: "Z" not in str(x), variables))
        ring_without_z = PolynomialRing(
            self.ring.base_ring(), variables + self.coefficients
        )
        Z1, Z2, X2, Y2 = (ring_with_x.gens_dict()[g] for g in ["Z1", "Z2", "X2", "Y2"])
        substitute_x = {Z2: 1 / X2, Y2: Y2 / X2}
        substitute_z = {Z1: 1}
        for row in self.formula:
            row.change_ring(ring_with_x)
            row.substitute(substitute_x, ring_with_x)
            row.substitute(substitute_z, ring_without_z)
        self.variables = variables
        self.ring = ring_without_z

    def xinvert(self):
        """
        Transformation Z1=1, Z2=x2,  Y2=x2/y2 for inverted:xmadd-2007-bl"""
        variables = self.variables + ["X2"]
        ring_with_x = PolynomialRing(
            self.ring.base_ring(), variables + self.coefficients
        )
        variables = list(filter(lambda x: "Z" not in str(x), variables))
        ring_without_z = PolynomialRing(
            self.ring.base_ring(), variables + self.coefficients
        )
        Z1, Z2, X2, Y2, X1, Y1 = (
            ring_with_x.gens_dict()[g] for g in ["Z1", "Z2", "X2", "Y2", "X1", "Y1"]
        )
        substitute_x = {Y1: 1 / Y1, X1: 1 / X1, Y2: X2 / Y2}
        substitute_z = {Z1: 1, Z2: X2}
        for row in self.formula:
            row.change_ring(ring_with_x)
            row.substitute(substitute_x, ring_with_x)
            row.substitute(substitute_z, ring_without_z)
        self.variables = variables
        self.ring = ring_without_z

    def remove_yz(self, squared=False):
        """get rid of yz coordinates"""
        self.coefficients, self.variables = ["r"], list(
            filter(lambda x: "Y" in str(x), self.variables)
        )
        Ys = list(self.ring.gens_dict()[g] for g in self.variables)
        r = self.ring.gens_dict()["r"]
        self.ring = PolynomialRing(
            self.ring.base_ring(), self.variables + self.coefficients
        )
        substitute = {Y: r * Y ** (squared + 1) for Y in Ys}
        for row in self.formula:
            row.substitute(substitute, self.ring)

    def remove_yzsquared(self):
        return self.remove_yz(squared=True)

    def set_c(self, a):
        self.a = self.curve_field(a)
        if "c" in self.coefficients:
            self.set_param("c", self.a)

    def set_d(self, b):
        self.b = self.curve_field(b)
        if "d" in self.coefficients:
            self.set_param("d", self.b)

    def set_r(self, r):
        self.r = self.curve_field(r)
        if "r" in self.coefficients:
            self.set_param("r", self.r)

    def set_curve(self, params):
        if "p" in params:
            self.set_p(params["p"])
        if "c" in params:
            self.set_c(params["c"])
        if "d" in params:
            self.set_d(params["d"])
        if "r" in params:
            self.set_r(params["r"])


class MontgomFormula(Formula):
    def __init__(self, dictionary):
        super().__init__(dictionary)
        for c in self.coefficients:
            assert c in ["a", "b"]
        self.a = None
        self.b = None

    def set_a(self, a):
        self.a = self.curve_field(a)
        if "a" in self.coefficients:
            self.set_param("a", self.a)

    def set_b(self, b):
        self.b = self.curve_field(b)
        if "b" in self.coefficients:
            self.set_param("b", self.b)

    def set_curve(self, params):
        if "p" in params:
            self.set_p(params["p"])
        if "a" in params:
            self.set_a(params["a"])
        if "b" in params:
            self.set_b(params["b"])

    def to_affine(self):
        self.put_Z_1()


class TwistedFormula(Formula):
    def __init__(self, dictionary):
        super().__init__(dictionary)
        for c in self.coefficients:
            assert c in ["a", "d"]
        self.transformation_methods["T_to_xy"] = self.remove_extension

    def remove_extension(self):
        """Get rid of extended coordinates"""
        Ts = list(
            self.ring.gens_dict()[g] for g in sorted(self.variables) if "T" in str(g)
        )
        Xs = list(
            self.ring.gens_dict()[g] for g in sorted(self.variables) if "X" in str(g)
        )
        Ys = list(
            self.ring.gens_dict()[g] for g in sorted(self.variables) if "Y" in str(g)
        )
        self.variables = list(filter(lambda x: "T" not in str(x), self.variables))
        self.ring = PolynomialRing(
            self.ring.base_ring(), self.variables + self.coefficients
        )
        substitute = {T: X * Y for (T, (X, Y)) in zip(Ts, zip(Xs, Ys))}
        for row in self.formula:
            row.substitute(substitute, self.ring)

    def set_a(self, a):
        self.a = self.curve_field(a)
        if "a" in self.coefficients:
            self.set_param("a", self.a)

    def set_d(self, d):
        self.d = self.curve_field(d)
        if "d" in self.coefficients:
            self.set_param("d", self.d)

    def set_curve(self, params):
        if "p" in params:
            self.set_p(params["p"])
        if "a" in params:
            self.set_a(params["a"])
        if "d" in params:
            self.set_d(params["d"])


class Row:
    def __init__(
        self, formula_string: tuple[str, str], ring: PolynomialRing, row_index: int
    ):
        self.left, right = formula_string
        self.ring = ring
        self.index = row_index
        self.numerator = self.ring(self.ring.fraction_field()(right).numerator())
        self.denominator = self.ring(self.ring.fraction_field()(right).denominator())
        self.flag = True

    def substitute(self, subs_dict, new_ring):
        subsnum = self.numerator.substitute(subs_dict)
        subsden = self.denominator
        if subsden != 1:
            subsden = self.denominator.substitute(subs_dict)
        new_frac = subsnum / subsden
        self.numerator = new_ring(new_frac.numerator())
        self.denominator = new_ring(new_frac.denominator())
        self.ring = new_ring

    def change_ring(self, new_ring):
        self.numerator = new_ring(self.numerator)
        self.denominator = new_ring(self.denominator)
        self.ring = new_ring

    def is_homogenous(self, weighted_variables):
        monomial_weights_n = compute_monomial_weights(
            self.numerator, weighted_variables
        )
        monomial_weights_d = compute_monomial_weights(
            self.denominator, weighted_variables
        )
        monomial_weights_n = set(monomial_weights_n)
        monomial_weights_d = set(monomial_weights_d)
        return len(monomial_weights_n) <= 1 and len(monomial_weights_d) <= 1


def compute_monomial_weights(polynomial, variable_weights):
    weights = []
    for monomial in polynomial.monomials():
        weights.append(sum(monomial.degree(var) * w for var, w in variable_weights))
    return weights


def instantiate_formula(formula_dict):
    forms = {
        "shortw": ShortWFormula,
        "montgom": MontgomFormula,
        "edwards": EdwardsFormula,
        "twisted": TwistedFormula,
    }
    return forms[formula_dict["form"]](formula_dict)


def load_formulas(path, affine=False):
    """
    Loads all the formulas from prepared json files. See 'unroll.py'
    """
    with open(path) as f:
        formulas_strings = json.load(f)
    all_formulas = {}
    for _, formulas in formulas_strings.items():
        for formula_name, formula_dict in formulas.items():
            formula = instantiate_formula(formula_dict)
            if affine:
                formula.to_affine()
            all_formulas[formula_name] = formula
    return all_formulas
