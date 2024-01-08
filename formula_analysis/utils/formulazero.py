import json
import re
from sage.all import gcd, PolynomialRing, QuotientRing, QQ
import utils.formula as fm


class FormulaZero:

    """
    Class that represents and sorts polynomials in an unrolled formula.
    """

    def __init__(self):
        self.formula = None
        self.full_name = None
        self.coordinates = None
        self.operation = None
        self.coefficients = None
        self.form = None
        self.ring = None
        self.numerators = []
        self.factor_set = set()
        self.factor_list = []
        self.single_point_factors = set()
        self.multi_point_factors = set()
        self.coordinate_factors = set()
        self.unsolvable_factors = set()

    @classmethod
    def from_formula(cls, formula: fm.Formula):
        f0 = cls()
        f0.formula = formula
        f0.full_name = formula.full_name()
        f0.coordinates = formula.coordinates
        f0.operation = formula.operation
        f0.coefficients = formula.coefficients
        f0.form = formula.form
        f0.ring = formula.ring

        f0._parse_numerators(formula)
        f0._create_factor_set(formula)
        f0._sort_numerators(formula)
        return f0

    def to_dict(self):
        result = {}
        result["coordinate_polynomials"] = list(
            (coord, str(f)) for coord, f in self.coordinate_zvp()
        )
        result["hard_polynomials"] = list(str(f) for f in self.hard_zvp())
        result["easy_polynomials"] = list(str(f) for f in self.easy_zvp())
        result["name"] = self.full_name
        result["operation"] = self.operation
        result["form"] = self.form
        result["coordinates"] = self.coordinates
        result["ring_variables"] = self.formula.ring.variable_names()
        result["unsolvable_polynomials"] = list(str(f) for f in self.unsolvable_zvp())
        return result

    @classmethod
    def from_dict(cls, dictionary):
        f0 = cls()
        f0.ring = PolynomialRing(QQ, dictionary["ring_variables"])
        f0.single_point_factors = set(
            f0.ring(f) for f in dictionary["easy_polynomials"]
        )
        f0.multi_point_factors = set(f0.ring(f) for f in dictionary["hard_polynomials"])
        f0.coordinate_factors = set(
            (coord, f0.ring(f)) for coord, f in dictionary["coordinate_polynomials"]
        )
        f0.unsolvable_factors = set(
            f0.ring(f) for f in dictionary["unsolvable_polynomials"]
        )
        f0.full_name = dictionary["name"]
        f0.coordinates = dictionary["coordinates"]
        f0.operation = dictionary["operation"]
        f0.form = dictionary["form"]
        return f0

    def easy_zvp(self):
        return self.single_point_factors

    def hard_zvp(self):
        return self.multi_point_factors

    def unsolvable_zvp(self):
        return self.unsolvable_factors

    def coordinate_zvp(self):
        return self.coordinate_factors

    def zvp_vulnerable(self):
        return self.easy_zvp() != set()

    def zvp_semi_resistant(self):
        return (not self.zvp_vulnerable()) and self.hard_zvp() != set()

    def zvp_resistant(self):
        return not (self.zvp_vulnerable() or self.zvp_semi_resistant())

    def remove_polynomials(self, to_remove_strings: list[str]):
        parsed_to_remove = set(
            normalize_polynomial(self.ring(f)) for f in to_remove_strings
        )
        self.single_point_factors = self.single_point_factors.difference(
            parsed_to_remove
        )
        self.multi_point_factors = self.multi_point_factors.difference(parsed_to_remove)
        self.coordinate_factors = set(
            filter(lambda x: x[1] not in parsed_to_remove, self.coordinate_factors)
        )
        self.unsolvable_factors = parsed_to_remove

    def _parse_numerators(self, formula):
        for row in formula.formula:
            if not row.flag:
                continue
            if row.numerator == 0:
                continue
            numerator = row.numerator / gcd(row.numerator.coefficients())
            if formula.is_constant(numerator):
                continue
            self.numerators.append([row.left, self.ring(numerator)])

    def _create_factor_set(self, formula):
        for [left, right] in self.numerators:
            factors = []
            for f, _ in right.factor():
                if formula.is_constant(f):
                    continue
                if f in self.ring.gens():
                    continue
                f_normalized = normalize_polynomial(f)
                self.factor_set.add(f_normalized)
                factors.append(f_normalized)
            self.factor_list.append([left, factors])

    def _is_single_point_factor(self, polynomial, formula):
        for point_index, variables in formula.input_variables_grouped().items():
            one_point_variables = set(variables).union(self.coefficients)
            if set(str(var) for var in polynomial.variables()).issubset(
                one_point_variables
            ):
                assert self._is_single_index_factor(polynomial)
                return True
        assert not self._is_single_index_factor(polynomial)
        return False

    def _is_single_index_factor(self, polynomial):
        pattern = re.compile(r"\d")
        indexes = set()
        for var in polynomial.variables():
            num = re.search(pattern, str(var))
            if num:
                indexes.add(int(num[0]))
        return len(indexes) == 1

    def _get_coordinate_if_factor(self, polynomial, formula):
        for coordinate, output_function in formula.output().items():
            num = output_function.numerator()
            if divides_mod(self, polynomial, num):
                return coordinate
        return None

    def _sort_numerators(self, formula):
        for polynomial in self.factor_set:
            if coordinate := self._get_coordinate_if_factor(polynomial, formula):
                self.coordinate_factors.add((coordinate, polynomial))
                continue
            if self._is_single_point_factor(polynomial, formula):
                self.single_point_factors.add(polynomial)
                continue
            self.multi_point_factors.add(polynomial)


class FormulaZeroCollection:
    def __init__(self, formula0s: dict[str, FormulaZero]):
        self.formula0s = formula0s

        self.ring = None
        self.create_common_ring()

        self.factor_set = set()
        self.single_point_factors = set()
        self.multi_point_factors = set()
        self.coordinate_factors = set()
        self.sort()

    def update(self):
        self.single_point_factors = set()
        self.multi_point_factors = set()
        self.coordinate_factors = set()
        self.sort()

    def is_empty(self):
        return self.formula0s == {}

    def __getitem__(self, item):
        return self.formula0s[item]

    def create_common_ring(self):
        if self.is_empty():
            return
        base_rings = set()
        variables = set()
        for formula in self:
            variables.update(set(formula.ring.variable_names()))
            base_rings.add(formula.ring.base_ring())
        assert len(base_rings) == 1, base_rings
        base_ring = list(base_rings)[0]
        self.ring = PolynomialRing(base_ring, list(variables))

    def sort(self):
        for formula0 in self:
            for f in formula0.factor_set:
                g = normalize_polynomial(self.ring(f))
                self.factor_set.add(normalize_polynomial(g))
            for coord, f in formula0.coordinate_factors:
                g = normalize_polynomial(self.ring(f))
                self.coordinate_factors.add((coord, normalize_polynomial(g)))
            for f in formula0.multi_point_factors:
                g = normalize_polynomial(self.ring(f))
                self.multi_point_factors.add(normalize_polynomial(g))
            for f in formula0.single_point_factors:
                g = normalize_polynomial(self.ring(f))
                self.single_point_factors.add(normalize_polynomial(g))

    def easy_zvp(self):
        return self.single_point_factors

    def hard_zvp(self):
        return self.multi_point_factors

    def coordinate_zvp(self):
        return self.coordinate_factors

    def zvp_resistant(self):
        collection = {
            formula0.full_name: formula0
            for formula0 in self
            if formula0.zvp_resistant()
        }
        formula0s = FormulaZeroCollection(collection)
        assert formula0s.easy_zvp() == set() and formula0s.hard_zvp() == set()
        return formula0s

    def zvp_semi_resistant(self):
        collection = {
            formula0.full_name: formula0
            for formula0 in self
            if formula0.zvp_semi_resistant()
        }
        formula0s = FormulaZeroCollection(collection)
        assert formula0s.easy_zvp() == set() and (
            formula0s.is_empty() or formula0s.hard_zvp() != set()
        )
        return formula0s

    def zvp_vulnerable(self):
        collection = {
            formula0.full_name: formula0
            for formula0 in self
            if formula0.zvp_vulnerable()
        }
        formula0s = FormulaZeroCollection(collection)
        assert formula0s.is_empty() or formula0s.easy_zvp() != set()
        return formula0s

    def __iter__(self):
        for formula0 in self.formula0s.values():
            yield formula0


def load_formula0s(formulas):
    formula0s = {}
    for name, formula in formulas.items():
        formula.flag_homogeneity()
        #formula.flag_multiplication_intermediate_values()
        formula.to_affine()
        formula0 = FormulaZero.from_formula(formula)
        formula0s[name] = formula0
    return FormulaZeroCollection(formula0s)


def save_zvp_results_per_formula(formula0s, path):
    results = {}
    for formula0 in formula0s:
        results[formula0.full_name] = formula0.to_dict()
    with open(path, "w") as f:
        json.dump(results, f)


def print_classified_formulas(formula0s):

    print("**Resistant:**")
    for formula0 in formula0s.zvp_resistant():
        print(formula0.full_name, end=", ")
    print("\n\n**Semi-resistant:**")
    semi_resistant = formula0s.zvp_semi_resistant()
    hard_zvp_sets = remove_unhashable_duplicities(
        f0.hard_zvp() for f0 in semi_resistant
    )
    for hard_zvp_set in hard_zvp_sets:
        for f0 in semi_resistant:
            if f0.hard_zvp() == hard_zvp_set:
                print(f0.full_name, end=", ")
        print(hard_zvp_set, "\n")
    print("**Vulnerable:**")
    for formula0 in formula0s.zvp_vulnerable():
        print(f"{formula0.full_name}: {formula0s.easy_zvp()}\n")


def print_classified_formulas_dbl(formula0s):

    print("**Resistant:**")
    for formula0 in formula0s.zvp_resistant():
        print(formula0.full_name, end=", ")
    assert formula0s.zvp_semi_resistant().is_empty()
    print("\n**Vulnerable:**")
    for formula0 in formula0s.zvp_vulnerable():
        print(f"{formula0.full_name}: {formula0s.easy_zvp()}\n")


def remove_unhashable_duplicities(list_: list):
    result = []
    for item in list_:
        if item not in result:
            result.append(item)
    return result


def divides_mod(formula0, f, g):
    ring = formula0.ring
    if formula0.formula.form == "shortw" and formula0.formula.operation == "add":
        a, b = ring.gens_dict().get("a", 0), ring.gens_dict().get("b", 0)
        X1, Y1 = ring.gens_dict()["X1"], ring.gens_dict()["Y1"]
        if formula0.formula.operation == "dbl":
            R = QuotientRing(ring, ring.ideal([Y1**2 - X1**3 - a * X1 - b, f]))
            return R(g) == 0
        X2, Y2 = ring.gens_dict()["X2"], ring.gens_dict()["Y2"]
        R = QuotientRing(
            ring,
            ring.ideal(
                [Y2**2 - X2**3 - a * X2 - b, Y1**2 - X1**3 - a * X1 - b, f]
            ),
        )
        g = R(g)
        return R(g) == 0
    return f.divides(g)


def normalize_polynomial(polynomial):
    return polynomial * (-1) ** (polynomial.coefficients()[0] < 0)
