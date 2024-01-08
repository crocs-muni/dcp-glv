from abc import ABC, abstractmethod
import re
from sage.all import PolynomialRing, sqrt
from utils.formulazero import FormulaZero


class Curve(ABC):
    @abstractmethod
    def polynomial_to_univariate(self, polynomial):
        pass

    @abstractmethod
    def solvable_on_curve(self, polynomial):
        pass

    @abstractmethod
    def supports_formula(self, formula: FormulaZero):
        pass

    @abstractmethod
    def has_zero_coordinate(self, coordinate: str):
        pass

    def substitute_coefficients(self, polynomial):
        for param in self.params:
            if gen := polynomial.parent().gens_dict().get(param, None):
                polynomial = polynomial.substitute({gen: self.params[param]})
        return polynomial

    def polynomial_ring_to_base_field(self, polynomial):
        variables = polynomial.parent().variable_names()
        parsed_polynomial = PolynomialRing(self.field, variables)(polynomial)
        return self.substitute_coefficients(parsed_polynomial)


class ShortWeierstrass(Curve):
    def __init__(self, a, b, name):
        self.a = a
        self.b = b
        self.field = a.parent()
        self.name = name
        self.form = "shortw"
        self.params = {"a": self.a, "b": self.b}

    def remove_y(self, polynomial):
        Xs, Ys = get_point_variables(polynomial)
        Ys = [y for y in Ys if y in polynomial.variables()]
        if Ys == []:
            return polynomial
        assert len(Ys) == 1
        Y = Ys[0] 
        X = get_matching_index_variable(Y, Xs)

        if polynomial.degree(Y) == 0:
            return polynomial

        assert polynomial.degree(Y) <= 2
        f0 = polynomial.coefficient({Y: 0})
        f1 = polynomial.coefficient({Y: 1})
        f2 = polynomial.coefficient({Y: 2})
        curve = X**3 + self.a * X + self.b
        return (curve * f2 + f0) ** 2 - curve * f1**2

    def polynomial_to_univariate(self, polynomial):
        return self.remove_y(polynomial)

    def solvable_on_curve(self, polynomial):
        polynomial = self.polynomial_ring_to_base_field(polynomial)
        polynomial = self.polynomial_to_univariate(polynomial)
        if polynomial == 0:
            return False
        if len(polynomial.parent().gens()) > 1:
            polynomial = polynomial.univariate_polynomial()
        for root, ex in polynomial.roots():
            if 0 in set([root**3 + self.a * root + self.b, root]):
                continue
            if (root**3 + self.a * root + self.b).is_square():
                return True
        return False

    def supports_formula(self, formula):
        if formula.form != self.form:
            return False
        coordinates = formula.coordinates
        if coordinates == "w12-0":
            return self.b == 0
        if coordinates[-1].isdigit():
            return self.a == -self.field(coordinates[-1])
        return True

    def has_zero_coordinate(self, coordinate: str):
        if not "X" in coordinate:
            return False
        return self.b.is_square()


class Montgomery(Curve):
    def __init__(self, a, b, name):
        self.a = a
        self.b = b
        self.name = name
        self.field = a.parent()
        self.form = "montgom"
        self.params = {"a": self.a, "b": self.b}

    def remove_y(self, polynomial):
        Xs, Ys = get_point_variables(polynomial)
        Ys = [y for y in Ys if y in polynomial.variables()]
        if Ys == []:
            return polynomial
        assert len(Ys) == 1
        Y = Ys[0]
        X = get_matching_index_variable(Y, Xs)

        if polynomial.degree(Y) == 0:
            return polynomial

        assert polynomial.degree(Y) <= 2
        f0 = polynomial.coefficient({Y: 0})
        f1 = polynomial.coefficient({Y: 1})
        f2 = polynomial.coefficient({Y: 2})
        curve = X**3 + self.a * X**2 + X
        return (curve * f2 + self.b * f0) ** 2 - curve * f1**2 * self.b

    def polynomial_to_univariate(self, polynomial):
        return self.remove_y(polynomial)

    def solvable_on_curve(self, polynomial):
        polynomial = self.polynomial_ring_to_base_field(polynomial)
        polynomial = self.polynomial_to_univariate(polynomial)
        if polynomial == 0:
            return False
        if len(polynomial.parent().gens()) > 1:
            polynomial = polynomial.univariate_polynomial()
        for root, _ in polynomial.roots():
            if 0 in set([root**3 + self.a * root**2 + root, root]):
                continue
            if (
                self.b * root**3 + self.b * self.a * root**2 + self.b * root
            ).is_square():
                return True
        return False

    def supports_formula(self, formula):
        return formula.form == self.form

    def has_zero_coordinate(self, coordinate: str):
        return False


class Edwards(Curve):
    def __init__(self, c, d, name):
        self.c = c
        self.d = d
        self.name = name
        self.field = d.parent()
        self.form = "edwards"
        self.params = {"c": self.c, "d": self.d}
        self.r = sqrt(d) if d.is_square() else None
        if self.r is not None:
            self.params["r"] = self.r

    def remove_x(self, polynomial):
        Xs, Ys = get_point_variables(polynomial)
        Xs = list(filter(lambda x: x in polynomial.variables(), Xs))
        if Xs == []:
            return polynomial
        assert len(Xs) == 1
        X = Xs[0]
        Y = get_matching_index_variable(X, Ys)

        if polynomial.degree(X) == 0:
            return polynomial
        assert polynomial.degree(X) <= 2
        f0 = polynomial.coefficient({X: 0})
        f1 = polynomial.coefficient({X: 1})
        f2 = polynomial.coefficient({X: 2})
        return (
            f0 * (1 - self.d * Y**2 * self.c**2) + f2 * (self.c**2 - Y**2)
        ) ** 2 - f1**2 * (self.c**2 - Y**2) * (1 - self.d * Y**2 * self.c**2)

    def polynomial_to_univariate(self, polynomial):
        return self.remove_x(polynomial)

    def solvable_on_curve(self, polynomial):

        polynomial = self.polynomial_ring_to_base_field(polynomial)
        polynomial = self.polynomial_to_univariate(polynomial)
        if polynomial == 0:
            return False
        if len(polynomial.parent().gens()) > 1:
            polynomial = polynomial.univariate_polynomial()
        for root, _ in polynomial.roots():
            if 0 in set(
                [
                    self.c**2 - root**2,
                    1 - self.c**2 * self.d * root**2,
                    root,
                ]
            ):
                continue
            if (
                (self.c**2 - root**2) / (1 - self.d * root**2 * self.c**2)
            ).is_square():
                return True
        return False

    def supports_formula(self, formula):
        if formula.form != self.form:
            return False
        coordinates = formula.coordinates
        if coordinates.startswith("yz"):
            if self.c == 1 and self.d.is_square():
                self.r = self.d.sqrt()
                return True
            return False
        return True

    def has_zero_coordinate(self, coordinate: str):
        return False


class TwistedEdwards(Curve):
    def __init__(self, a, d, name):
        self.a = a
        self.d = d
        self.name = name
        self.field = d.parent()
        self.form = "twisted"
        self.params = {"a": self.a, "d": self.d}

    def remove_x(self, polynomial):
        Xs, Ys = get_point_variables(polynomial)
        Xs = list(filter(lambda x: x in polynomial.variables(), Xs))
        if Xs == []:
            return polynomial
        assert len(Xs) == 1
        X = Xs[0]
        Y = get_matching_index_variable(X, Ys)

        if polynomial.degree(X) == 0:
            return polynomial
        assert polynomial.degree(X) <= 2
        f0 = polynomial.coefficient({X: 0})
        f1 = polynomial.coefficient({X: 1})
        f2 = polynomial.coefficient({X: 2})
        return (f0 * (self.a - self.d * Y**2) + f2 * (1 - Y**2)) ** 2 - f1**2 * (
            1 - Y**2
        ) * (self.a - self.d * Y**2)

    def polynomial_to_univariate(self, polynomial):
        return self.remove_x(polynomial)

    def solvable_on_curve(self, polynomial):
        polynomial = self.polynomial_ring_to_base_field(polynomial)
        polynomial = self.polynomial_to_univariate(polynomial)
        if polynomial == 0:
            return False
        if len(polynomial.parent().gens()) > 1:
            polynomial = polynomial.univariate_polynomial()
        for root, _ in polynomial.roots():
            if 0 in set([1 - root**2, self.a - self.d * root**2, root]):
                continue
            if ((1 - root**2) / (self.a - self.d * root**2)).is_square():
                return True
        return False

    def supports_formula(self, formula):
        if formula.form != self.form:
            return False
        coordinates = formula.coordinates
        if coordinates[-1].isdigit():
            return self.a == -1
        return True

    def has_zero_coordinate(self, coordinate: str):
        return False


def get_point_variables(polynomial):
    Ys = [v for v in polynomial.parent().gens() if "Y" in str(v)]
    Xs = [v for v in polynomial.parent().gens() if "X" in str(v)]
    return Xs, Ys

def extract_index(variable):
    pattern = re.compile(r"\d")
    return int(re.search(pattern, str(variable))[0])

def get_matching_index_variable(var, candidates):
    for v in candidates:
        if extract_index(v) == extract_index(var):
            return v
        