"""
Functions for solving DCP using pari/gp.
"""
import os
import re


def dcpsolver_pari2(p, a, b, k):
    """DCP solver for 2-dimensional scalar decomposition.
    Returns a list of x-coordinates."""
    pari_path = "./pari_tools/dcp_solver"
    solution_path = "./pari_tools/results/poly2"
    os.system(f"{pari_path} {p} {k} {a} {b} {solution_path}")
    with open(solution_path, "r") as f:
        result = f.read()
    return re.findall(r"Mod\((\d+), \d+\)", result)


def dcpsolver_pari4(p, p1, p2, a1, a2, b1, b2, k):
    """DCP solver for 4-dimensional scalar decomposition.
    Returns a list of x-coordinates (tuples)."""
    pari_path = "./dcp_dim4"
    solution_path = "./results/poly4"
    os.system(f"{pari_path} {p} {p1} {p2} {a1} {a2} {b1} {b2} {k} {solution_path}")
    with open(solution_path, "r") as f:
        result = f.read()
    result = result.split("[")[1].split("]")[0].replace(" ", "").split(",")
    if result == [""]:
        return []
    converted = []
    for r in result:
        if "x" in r:
            m = re.search(r"(\d*)\**x\+*(\d*)", r)
            r1, r2 = m.group(1), m.group(2)
            if r1 == "":
                converted.append(("1", r2))
                continue
            if r2 == "":
                converted.append((r1, "0"))
                continue
            converted.append((r1, r2))
        else:
            converted.append(("0", r))

    return converted


def multidcpsolver_pari2(p, a, b, beta, k1, k2, guess):
    """DCP solver for 2-dimensional multiscalar decomposition.
    Returns a list of x-coordinates."""
    pari_path = "./pari_tools/dcp_glv_solver"
    solution_path = "./pari_tools/results/multipoly2"
    os.system(f"{pari_path} {p} {a} {b} {beta} {k1} {k2} {guess} {solution_path}")
    with open(solution_path, "r") as f:
        result = f.read()
    return re.findall(r"Mod\((\d+), \d+\)", result)


def division_polynomial_tuples(p, a, b, m, pari_path, polynomial_path):
    """
    Obsolete function.
    Computes division polynomial using pari/gp.
    Returns the same thing as .division_polynomial in sage in tuple form
    """
    b4 = 2 * a
    b6 = 4 * b
    b8 = -(a**2)
    if m == -1:
        return [("4", "3"), (f"{2*b4}", "1"), (f"{b6}", "0")]
    if m == -2:
        return [
            ("16", "6"),
            (f"{16*b4}", "4"),
            (f"{4*b4**2}", "2"),
            (f"{8*b6}", "3"),
            (f"{4*b4*b6}", "1"),
            (f"{b6**2}", "0"),
        ]

    try:
        open(f"{polynomial_path}/{p}_{a}_{b}_{m}.txt")
    except FileNotFoundError:
        os.system(f"{pari_path} {p} {m} {a} {b} {polynomial_path}")
    with open(f"{polynomial_path}/{p}_{a}_{b}_{m}.txt") as f:
        r = f.read()
    r = r.replace(" ", "")
    result = re.findall(r"Mod\((\d+),\d+\)\*x\^(\d+)", r)
    modsplit = r.split("Mod(")
    try:
        lead = modsplit[0].split("+")[0]
        result = [(lead.split("*")[0], lead.split("^")[1])] + result
    except IndexError:
        pass
    try:
        constant = re.findall(r"Mod\((\d+),\d+\)$", r)[0]
        result.append((constant, "0"))
    except IndexError:
        pass
    try:
        linear = re.findall(r"Mod\((\d+),\d+\)\*x$", r)[0]
        result.append((linear, "1"))
    except IndexError:
        pass
    try:
        linear = re.findall(r"Mod\((\d+),\d+\)\*x\+", r)[0]
        result.append((linear, "1"))
    except IndexError:
        pass
    return result


def division_polynomial(x, p, a, b, m, pari_path, solution_path):
    """
    Obsolete function.
    Computes division polynomial using pari/gp.
    returns the same thing as .division_polynomial in sage"""
    ring = x.parent()
    return sum(
        [
            ring(coef) * x ** int(exp)
            for coef, exp in division_polynomial_tuples(
                p, a, b, m, pari_path, solution_path
            )
        ]
    )


def division_polynomial_0(x, p, a, b, m, pari_path, solution_path):
    """
    Obsolete function.
    Computes division polynomial using pari/gp.
    returns the same thing as .division_polynomial_0 in sage"""
    ring = x.parent()
    f = sum(
        [
            ring(coef) * x ** int(exp)
            for coef, exp in division_polynomial_tuples(
                p, a, b, m, pari_path, solution_path
            )
        ]
    )
    if m % 2 == 0:
        return ring(f / division_polynomial(x, p, a, b, -1))
    return f
