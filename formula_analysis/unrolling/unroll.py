import os
import re
import json
import sys
from sage.all import PolynomialRing, QQ


def unroll_curve_form():
    """Main Function"""
    fn = FormulaNavigation()
    for form in fn.supported_forms():
        fn.form = form
        if not fn.valid_form_path():
            continue
        unrolled = unroll_form(fn)
        save_formulas_to_json(fn, unrolled)
        save_formula_list(form, unrolled, "formula_list.txt")


class FormulaNavigation:

    """Data structure containing location of formulas to unroll in the efd database"""

    def __init__(self):
        self.path_to_efd = sys.argv[1]
        self.form = None
        self.coordinates = None
        self.operation = None
        self.formula_file = None

        self._supported_forms = ["edwards", "montgom", "shortw", "twisted"]
        self._supported_operations = ["dadd", "dbl", "ladd", "add"]
        self.formulas_to_skip = set(
            [
                ("shortw", "jacobian-0", "dbl-1998-hnm"),
                ("edwards", "projective", "add-20080225-hwcd"),
                ("edwards", "projective", "madd-20080225-hwcd"),
                ("edwards", "projective", "add-2007-bl-4"),
                ("edwards", "projective", "add-20090311-hwcd"),
                ("shortw", "jacobian-0", "add-1998-hnm"),
                ("shortw", "jacobian", "dbl-1998-hnm"),
                ("shortw", "jacobian", "add-1998-hnm"),
                ("shortw", "jacobian-3", "dbl-1998-hnm"),
                ("shortw", "jacobian-3", "add-1998-hnm"),
                ("shortw", "jacobian-3", "dbl-2004-hmv"),
                ("shortw", "jacobian-3", "dbl-1998-hnm-2"),
            ]
        )
        self.operation_dir_names = {
            "dadd": "diffadd",
            "dbl": "doubling",
            "ladd": "ladder",
            "add": "addition",
        }

    def supported_forms(self):
        return self._supported_forms

    def supported_operations(self):
        return self._supported_operations

    def path_to_form(self):
        return os.path.join(self.path_to_efd, self.form)

    def path_to_coordinates(self):
        return os.path.join(self.path_to_form(), self.coordinates)

    def path_to_operation(self):
        return os.path.join(
            self.path_to_coordinates(), self.operation_dir_names[self.operation]
        )

    def path_to_formula_file(self):
        return os.path.join(self.path_to_operation(), self.formula_file)

    def path_to_short_formula_file(self):
        return os.path.join(
            self.path_to_operation(), self.formula_file.split(".op3")[0]
        )

    def path_to_curve_params(self):
        return os.path.join(self.path_to_form(), "coordinates")

    def path_to_formula_variables(self):
        return os.path.join(self.path_to_coordinates(), "variables")

    def formula_name(self):
        return self.formula_file[:-4]

    def formula_full_name(self):
        return f"{self.coordinates}:{self.formula_name()}"

    def supported_coordinates(self):
        return os.listdir(self.path_to_form())

    def supported_formulas(self):
        return os.listdir(self.path_to_operation())

    def valid_form_path(self):
        return os.path.isdir(self.path_to_form())

    def valid_coordinate_path(self):
        return os.path.isdir(self.path_to_coordinates())

    def valid_operation_path(self):
        return os.path.isdir(self.path_to_operation())

    def valid_formula_file(self):
        if not self.formula_file.endswith(".op3"):
            return False
        return (
            not (self.form, self.coordinates, self.formula_name())
            in self.formulas_to_skip
        )
    def path_to_transformation(self):
        filepath = os.path.join(self.path_to_operation(),f"unified_transformation.{self.formula_name()}")
        if os.path.isfile(filepath):
            return filepath
        filepath = os.path.join(self.path_to_coordinates(),f"unified_transformation")
        if os.path.isfile(filepath):
            return filepath

    def path_to_homogeneity_weights(self):
        filepath = os.path.join(self.path_to_coordinates(),f"homogeneity_weights")
        if os.path.isfile(filepath):
            return filepath 



def unroll_form(fn: FormulaNavigation) -> dict:
    """
    Creates dict = {curve_form:unrolled_formulas_in_the_form,..}
    """
    formulas = {}
    for coordinates in fn.supported_coordinates():
        fn.coordinates = coordinates
        if not fn.valid_coordinate_path():
            continue
        formulas[coordinates] = unroll_coordinates(fn)
    return formulas


def unroll_coordinates(fn: FormulaNavigation) -> dict:
    """
    Creates dict = {coordinates:unrolled_formulas_in_the_coordinate_system,..}
    """
    formulas = {}
    for operation in fn.supported_operations():
        fn.operation = operation
        if not fn.valid_operation_path():
            continue
        formulas[operation] = unroll_operation(fn)
    return formulas


def unroll_operation(fn: FormulaNavigation) -> dict:
    """
    Creates dict = {formula_name:unrolled_formula,..}
    """
    formulas = {}
    for formula_file in fn.supported_formulas():
        fn.formula_file = formula_file
        if not fn.valid_formula_file():
            continue
        formula = Formula(fn)
        unroll_formula(formula, fn)
        formulas[formula.formula_full_name] = formula.to_dict()
    return formulas


class Formula:

    """Data structure containing the unrolled formula"""

    def __init__(self, fn: FormulaNavigation):

        with open(fn.path_to_formula_file()) as f:
            lines = f.readlines()
        self.path_to_transformation = fn.path_to_transformation()
        self.path_to_homogeneity_weights = fn.path_to_homogeneity_weights()
        self.lines = lines
        self.unrolled_lines = []
        self.parameters = None
        self.operation = fn.operation
        self.form = fn.form
        self.coordinates = fn.coordinates
        self.output_indexes = {}
        self.secondary_output_indexes = {}
        self.name = fn.formula_name()
        self.formula_full_name = fn.formula_full_name()

    def extract_all_variables(self):
        new_vars = set()
        for left, right in self.unrolled_lines:
            new_vars = new_vars.union(
                set(
                    str(v)
                    for v in right.numerator().variables()
                    + right.denominator().variables()
                )
            )
        return sorted(list(new_vars - set(self.parameters)))

    def save_line_if_output(self, left: str, index: int):
        if self.operation == "ladd":
            if left in ["Z5", "X5", "Y5"]:
                self.output_indexes[left[0]] = index
            if left in ["Z4", "X4", "Y4"]:
                self.secondary_output_indexes[left[0]] = index
            return
        for var_prefix in ["X", "Y", "Z", "ZZ", "ZZZ"]:
            pattern = re.compile(r"^" + var_prefix + r"\d$")
            if re.findall(pattern, left):
                self.output_indexes[var_prefix] = index

    def extract_transformation(self):
        with open(self.path_to_transformation) as f:
            transformation = f.readlines()
        return list(map(lambda x: x.strip(), transformation))

    def extract_homogeneity_weights(self):
        with open(self.path_to_homogeneity_weights) as f:
            weights = f.readlines()
        return {left: int(right) for left,right in map(parse_line_string, weights)}

    def to_dict(self):
        unrolled = {}
        unrolled["parameters"] = self.parameters
        unrolled["variables"] = self.extract_all_variables()
        unrolled["formula"] = [
            [left, polynomial_to_string(right)] for left, right in self.unrolled_lines
        ]
        unrolled["operation"] = self.operation
        unrolled["form"] = self.form
        unrolled["name"] = self.name
        unrolled["output"] = self.output_indexes
        unrolled["secondary_output"] = self.secondary_output_indexes
        unrolled["coordinates"] = self.coordinates
        unrolled["original_formula"] = list(map(parse_line_string, self.lines))
        unrolled["unified_transformation"] = self.extract_transformation()
        unrolled["homogeneity_weights"] = self.extract_homogeneity_weights()
        return unrolled


def unroll_formula(formula: Formula, fn: FormulaNavigation):
    """
    Creates dict = {formula:_,variables:_,parameters:_,operation:_,form:_,name:_}
    """
    substitutions, variable_names, formula.parameters = initial_substitutions(fn)
    ring = PolynomialRing(QQ, variable_names).fraction_field()
    substitutions = {left: ring(right) for left, right in substitutions.items()}
    for index, line in enumerate(formula.lines):
        left, right = parse_line_string(line)
        unrolled_right = unroll_line(right, substitutions, ring)
        substitutions[left] = unrolled_right
        formula.unrolled_lines.append([left, unrolled_right])
        formula.save_line_if_output(left, index)


def unroll_line(rightside_of_line, unrolled_lines, ring):
    variables_to_substitute = extract_variables(rightside_of_line)
    exring = extend_ring(ring, variables_to_substitute)
    casted_right = exring(rightside_of_line)
    substitution = {
        exring(var): exring(unrolled_lines[var]) for var in variables_to_substitute
    }
    unrolled_right = casted_right.substitute(substitution)
    return ring(unrolled_right)


def initial_substitutions(fn: FormulaNavigation) -> tuple[dict, PolynomialRing, list]:
    """
    Transform all conditions on curve_parameters, variables and parameters in assumptions to a dictionary.
    Return the dictionary, polynomial ring and parameters necessary for the substitutions
    """
    curve_parameters = get_symbols(fn.path_to_curve_params(), "parameter")
    variables = get_symbols(fn.path_to_formula_variables(), "variable")
    assumptions = get_symbols(fn.path_to_formula_variables(), "assume")
    parameters = get_symbols(fn.path_to_formula_variables(), "parameter")
    formula_assumptions = get_symbols(fn.path_to_short_formula_file(), "assume")

    substitutions = dict(zip(parameters, parameters))
    substitutions.update(map(parse_line_string, assumptions))

    for cparam in curve_parameters:
        if cparam not in substitutions:
            substitutions[cparam] = cparam
            parameters.append(cparam)

    substitutions.update(map(parse_line_string, formula_assumptions))

    variable_names = [f"{var}{i}" for var in variables for i in range(10)]
    substitutions.update({var: var for var in variable_names})
    variable_names.extend(parameters)

    return substitutions, variable_names, parameters


def parse_line_string(line):
    left, right = map(lambda x: x.strip(), line.split("="))
    return left, right.replace("^", "**")


def polynomial_to_string(polynomial):
    return str(polynomial).replace("^", "**")


def extract_variables(poly_string):
    pattern = re.compile(r"[A-Za-z]+[A-Za-z0-9]*")
    return list(re.findall(pattern, poly_string))


def extend_ring(old_ring, new_variables):
    tmp_vars = list(
        set(old_ring.variable_names() + tuple(new_variables))
    )  # creates a new ring with all accumulated variables and the new ones
    return PolynomialRing(QQ, tmp_vars).fraction_field()


def get_symbols(path: str, name: str) -> list:
    symbols = []
    with open(path) as f:
        lines = f.readlines()
    for line in lines:
        if line.split(" ")[0].strip() == name:
            symbols.append(line.split(name)[1].strip())
    return symbols


def save_formulas_to_json(fn: FormulaNavigation, unrolled: dict):
    """
    Save formulas in form_operation.json for each form/op
    """
    if not os.path.exists("unrolled"):
        print("Creating folder unrolled")
        os.makedirs("unrolled")
    for operation in fn.supported_operations():
        saved_formulas = {}
        for coordinates, unrolled_coordinates in unrolled.items():
            if operation in unrolled_coordinates:
                saved_formulas[coordinates] = unrolled_coordinates[operation]
        with open(os.path.join("unrolled", f"{fn.form}_{operation}.json"), "w") as f:
            print(f"Saving unrolled formulas for {fn.form}, {operation}")
            json.dump(saved_formulas, f, indent=1)


def save_formula_list(form: str, unrolled: dict, path: str):
    """
    Save a list of unrolled formulas
    """
    formula_list = ""
    for _, coordinates in unrolled.items():
        for operation, formulas in coordinates.items():
            formula_list += (
                "\n".join(map(lambda x: f"{form}:{operation}:{x}", formulas.keys()))
                + "\n"
            )
    with open(path, "a") as h:
        h.write(formula_list)


if __name__ == "__main__":
    unroll_curve_form()
