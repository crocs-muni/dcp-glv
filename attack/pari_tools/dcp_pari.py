"""
Functions for solving DCP using pari/gp.
"""
import os
import re


def dcpsolver_pari(p, a, b, k, lam, Vpolynomials, registers):
    """DCP solver for 2-dimensional scalar decomposition.
    Returns a list of x-coordinates."""
    pari_path = "./pari_tools/dcp_solver"
    solution_path = "./pari_tools/results/poly2"
    Vpolynomial_path = "./pari_tools/results/_polynomial"
    fpolynomial_path = "./pari_tools/results/_fpolynomial"
    V_string = polynomials_to_parituples(Vpolynomials, p)
    f_string = polynomials_to_paristring(registers.polynomials_y, p)
    with open(Vpolynomial_path, "w") as f:
        f.write(V_string)
    with open(fpolynomial_path, "w") as f:
        f.write(f_string)

    call_c(
        f"{pari_path} {p} {k} {lam} {a} {b} {Vpolynomial_path} {fpolynomial_path} {solution_path}"
    )
    return read_solution(solution_path)


def multidcpsolver_pari(p, a, b, k, l, lam, Vpolynomials, registers):
    """DCP solver for 2-dimensional scalar decomposition.
    Returns a list of x-coordinates."""
    pari_path = "./pari_tools/multidcp_solver"
    solution_path = "./pari_tools/results/poly2"
    Vpolynomial_path = "./pari_tools/results/_polynomial"
    fpolynomial_path = "./pari_tools/results/_fpolynomial"
    V_string = polynomials_to_parituples(Vpolynomials, p)
    f_string = polynomials_to_paristring(registers.polynomials_y, p)
    with open(Vpolynomial_path, "w") as f:
        f.write(V_string)
    with open(fpolynomial_path, "w") as f:
        f.write(f_string)

    call_c(
        f"{pari_path} {p} {k} {l} {lam} {a} {b} {Vpolynomial_path} {fpolynomial_path} {solution_path}"
    )
    return read_solution(solution_path)


def glvdcpsolver_pari(p, a, b, k1, lam, k2, Vpolynomials, registers):
    """DCP solver for 2-dimensional multiscalar decomposition.
    Returns a list of x-coordinates."""
    pari_path = "./pari_tools/dcp_glv_solver"
    solution_path = "./pari_tools/results/multipoly4"
    Vpolynomial_path = "./pari_tools/results/_Vpolynomial"
    fpolynomial_path = "./pari_tools/results/_fpolynomial"
    V_string = polynomials_to_parituples(Vpolynomials, p)
    f_string = polynomials_to_paristring(registers.polynomials_y, p)
    with open(Vpolynomial_path, "w") as f:
        f.write(V_string)
    with open(fpolynomial_path, "w") as f:
        f.write(f_string)
    call_c(
        f"{pari_path} {p} {a} {b} {k1} {lam} {k2} {Vpolynomial_path} {fpolynomial_path} {solution_path}"
    )
    return read_solution(solution_path)


def glvdcpmultisolver_pari(p, a, b, scalar0, k1, lam, k2, Vpolynomials, registers):
    """DCP solver for 2-dimensional multiscalar decomposition.
    Returns a list of x-coordinates."""
    pari_path = "./pari_tools/multidcp_glv_solver"
    solution_path = "./pari_tools/results/multipoly4"
    Vpolynomial_path = "./pari_tools/results/_Vpolynomial"
    fpolynomial_path = "./pari_tools/results/_fpolynomial"
    V_string = polynomials_to_parituples(Vpolynomials, p)
    f_string = polynomials_to_paristring(registers.polynomials_y, p)
    with open(Vpolynomial_path, "w") as f:
        f.write(V_string)
    with open(fpolynomial_path, "w") as f:
        f.write(f_string)
    call_c(
        f"{pari_path} {p} {a} {b} {scalar0} {k1} {lam} {k2} {Vpolynomial_path} {fpolynomial_path} {solution_path}"
    )
    return read_solution(solution_path)


def polynomials_to_parituples(polynomials, p):
    poly_strings = []
    for polynomial in polynomials:
        assert polynomial.parent().variable_names() in [
            ("x", "n1", "d1", "n2", "d2"),
            ("x", "n", "d"),
            ("n0", "d0", "n1", "d1", "n2", "d2"),
            ("n0", "d0", "n", "d"),
        ], polynomial.parent().variable_names()
        monom_strings = []
        for monomial in polynomial.iterator_exp_coeff():
            monom = [monomial[1]]+list(monomial[0])
            monom_strings.append(str(monom))
        poly_strings.append("["+",".join(monom_strings)+"]")

    return "[" + ",".join(poly_strings) + "]"



def polynomials_to_paristring(polynomials, p):
    poly_strings = []
    for polynomial in polynomials:
        assert polynomial.parent().variable_names() in [
            ("X1", "Y1", "X2", "Y2", "A", "B"),
        ], polynomial.parent().variable_names()
        monom_strings = []
        for monomial in polynomial.monomials():
            coef = polynomial.monomial_coefficient(monomial)
            monom_strings.append(f"Mod({coef},{p})*{monomial}")
        poly_strings.append("+".join(monom_strings))

    return "[" + ",".join(poly_strings) + "]"



def call_c(command):
    os.system(command)


def read_solution(solution_path):
    with open(solution_path, "r") as f:
        result = f.read()
    solution = re.findall(r"Mod\((\d+), \d+\)", result)
    if not solution:
        return []
    x, y = solution
    return [(int(x), int(y))]
