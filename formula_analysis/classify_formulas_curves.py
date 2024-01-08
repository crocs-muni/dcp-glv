import json, os
from sage.all import GF, PolynomialRing
import pandas, pickle
from copy import deepcopy
from utils.load_curves import load_curves
from utils.curve import Curve
from utils.formulazero import FormulaZero


class ZVPResistanceDatabase:
    def __init__(self, form: str, operation: str):
        self.form = form
        self.operation = operation
        self.curves = []
        self.formulas = []
        self.zvpresistances = []
        self.results_dir = "results"

    def load_curves(self, curves: list[Curve]):
        for curve in curves:
            if curve.form == self.form:
                self.curves.append(curve)

    def load_formulas(self, formulas: list[dict]):
        for _, formula in formulas.items():
            if formula["operation"] == self.operation and formula["form"] == self.form:
                self.formulas.append(FormulaZero.from_dict(formula))

    def load_formulas_from_json(self, path):
        with open(path) as f:
            self.load_formulas(json.load(f))

    def classify_curve_formula(self, curve: Curve, formula: FormulaZero):
        self.zvpresistances.append(Resistance(curve=curve, formula=formula))

    def valid_combinations(self):
        for curve in self.curves:
            for formula in self.formulas:
                if curve.supports_formula(formula):
                    yield curve, formula

    def classify_all(self):
        for curve, formula in self.valid_combinations():
            self.classify_curve_formula(curve, formula)

    def to_pandas_per_polynomial(self):
        dataframe_dict = {
            "curve": [],
            "formula": [],
            "polynomial": [],
            "difficulty": [],
        }
        for curve_formula_resistance in self.zvpresistances:
            for (
                polynomial,
                polynomial_resistance,
            ) in curve_formula_resistance.resistance_per_polynomial():
                dataframe_dict["curve"].append(curve_formula_resistance.curve.name)
                dataframe_dict["formula"].append(
                    curve_formula_resistance.formula.full_name
                )
                dataframe_dict["polynomial"].append(str(polynomial))
                dataframe_dict["difficulty"].append(polynomial_resistance)
        return pandas.DataFrame(dataframe_dict)

    def to_pandas(self):
        dataframe_dict = {"curve": [], "formula": [], "resistance": []}
        for curve_formula_resistance in self.zvpresistances:
            dataframe_dict["curve"].append(curve_formula_resistance.curve.name)
            dataframe_dict["formula"].append(curve_formula_resistance.formula.full_name)
            dataframe_dict["resistance"].append(curve_formula_resistance.resistance())
        return pandas.DataFrame(dataframe_dict)

    def to_pickle(self):
        path = os.path.join(
            "results", f"{self.form}_{self.operation}_polynomial_database.pkl"
        )
        self.to_pandas_per_polynomial().to_pickle(path)

        path = os.path.join(
            "results", f"{self.form}_{self.operation}_resistance_database.pkl"
        )
        self.to_pandas().to_pickle(path)

    def to_html(self):
        path = os.path.join(
            "results", f"{self.form}_{self.operation}_polynomial_database.html"
        )
        self.to_pandas_per_polynomial().to_html(path)

        path = os.path.join(
            "results", f"{self.form}_{self.operation}_resistance_database.html"
        )
        self.to_pandas().to_html(path)

    def save_resistance_changes(self):
        log = ""
        for curve_formula_resistance in self.zvpresistances:
            if not curve_formula_resistance.resistance_changed():
                continue
            log += f"{curve_formula_resistance.curve.name}+{curve_formula_resistance.formula.full_name}"
            log += f":{curve_formula_resistance.resistance_changed()}\n"
        path = os.path.join("results", f"{self.form}_{self.operation}_resistance_changes.txt")
        with open(path, "w") as h:
            h.write(log)


class Resistance:
    def __init__(self, curve: Curve, formula: FormulaZero):
        self.curve = deepcopy(curve)
        self.formula = deepcopy(formula)
        self.hard_flag = "Hard"
        self.easy_flag = "Easy"
        self.unsolvable_flag = "Unsolvable"

        self.resistant_flag = "Resistant"
        self.semi_resistant_flag = "Semi-resistant"
        self.vulnerable_flag = "Vulnerable"

        self.general_formula_resistance = self.resistance()
        self.reclassify_hard()
        self.reclassify_easy()


    def reclassify_easy(self):
        to_remove = []
        for polynomial in self.formula.easy_zvp():
            if not self.formula._is_single_index_factor(polynomial):
                assert False
            if not self.curve.solvable_on_curve(polynomial):
                to_remove.append(str(polynomial))
        self.formula.remove_polynomials(to_remove)

    def reclassify_hard(self):
        for polynomial in self.formula.hard_zvp():
            polynomial = self.curve.polynomial_ring_to_base_field(polynomial)
            if self.formula._is_single_index_factor(polynomial):
                assert (
                    False
                ), "Attempted reclasification of Hard polynomial, Not implemented"

    def resistance(self):
        if self.formula.zvp_vulnerable():
            return self.vulnerable_flag
        if self.formula.zvp_semi_resistant():
            return self.semi_resistant_flag
        if self.formula.zvp_resistant():
            return self.resistant_flag

    def resistance_per_polynomial(self):

        for polynomial in self.formula.hard_zvp():
            yield polynomial, self.hard_flag
        for polynomial in self.formula.easy_zvp():
            yield polynomial, self.easy_flag
        for polynomial in self.formula.unsolvable_zvp():
            yield polynomial, self.unsolvable_flag

    def resistance_changed(self):
        if self.general_formula_resistance != self.resistance():
            return f"{self.general_formula_resistance}>{self.resistance()}"
        return False



def main():
    curves = load_curves()
    for operation in "add", "ladd", "dadd", "dbl":
        for form in "shortw", "edwards", "montgom", "twisted":
            print(f"Classifying {operation} {form} with {len(curves[form])} curves")
            dat = ZVPResistanceDatabase(form, operation)
            dat.load_curves(curves[form])
            try:
                dat.load_formulas_from_json(f"results/{form}_{operation}.json")
            except FileNotFoundError:
                continue
            dat.classify_all()
            dat.to_pickle()
            dat.to_html()
            dat.save_resistance_changes()


if __name__ == "__main__":
    main()
