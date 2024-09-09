from sage.all import ZZ, ceil, log, GF, EllipticCurve, PolynomialRing
import glv as glv_module
from random import randint
from enum import Enum
from copy import deepcopy
import os
import json


class GuessVerification(Enum):
    TRUE = 1
    FALSE = 0
    UNKNOWN = 2


class Registers:
    def __init__(self):
        self.ring = PolynomialRing(ZZ, "X1,Y1,X2,Y2,A,B")
        self.gens = self.ring.gens()[:-2]
        self.all_gens = self.ring.gens()
        self.polynomials = []
        self.polynomials_tuples = []

    def empty_out(self):
        self.polynomials = []
        self.polynomials_tuples = []

    def add(self, polynomial):
        self.polynomials.append(self.ring(polynomial))
        self.polynomials_tuples.append((self.ring(polynomial), self.ring(polynomial)))

    def add_tuple(self, xpoly, poly):
        self.polynomials.append(self.ring(xpoly))
        self.polynomials_tuples.append((self.ring(xpoly), self.ring(poly)))

    @property
    def polynomials_y(self):
        return [i[1] for i in self.polynomials_tuples]

    def is_zero_x(self, point1, point2):
        if not self.polynomials:
            raise Exception("No registers set\n.")
        a, b = point1.curve().a4(), point1.curve().a6()
        if point1.curve()(0) == point1 or point1.curve()(0) == point2:
            return False
        for polynomial in self.polynomials:
            # print(polynomial,point1,point2)
            if polynomial(point1[0], point1[1], point2[0], point2[1], a, b) == 0:
                # print(0,polynomial.parent())
                return True
        return False

    def is_zero(self, point1, point2):
        if not self.polynomials:
            raise Exception("No registers set\n.")
        a, b = point1.curve().a4(), point1.curve().a6()
        if point1.curve()(0) == point1 or point1.curve()(0) == point2:
            return False
        for _, polynomial in self.polynomials_tuples:
            # print(polynomial,point1,point2)
            if polynomial(point1[0], point1[1], point2[0], point2[1], a, b) == 0:
                # print(polynomial(point1[0], point1[1], point2[0], point2[1], a, b))
                # print(0,polynomial.parent())
                return True
        return False

    def to_strings(self):
        return [(str(f), str(g)) for f, g in self.polynomials_tuples]


def get_secp256k1_polynomials(extended=False):
    X1, Y1, X2, Y2, A, B = Registers().all_gens
    f1 = X1 + X2, X1 + X2
    # 4

    f2 = (
        -(X1**4)
        + 6 * X1**3 * X2
        - 7 * X1**2 * X2**2
        + 2 * X1 * X2**3
        - X2**4
        + 8 * X1 * B,
        2 * X1**3
        - 4 * X1**2 * X2
        + 2 * X1 * X2**2
        - Y1**2
        + 2 * Y1 * Y2
        - Y2**2,
    )  # = 2*X1*(X1-X2)**2 - (Y1-Y2)**2
    # 2

    f3 = (
        -(X1**4)
        + 4 * X1**3 * X2
        + 6 * X1**2 * X2**2
        + 4 * X1 * X2**3
        - X2**4
        + 24 * X1 * B
        + 24 * X2 * B,
        2 * X1**4
        + 4 * X1**3 * X2
        + 6 * X1**2 * X2**2
        + 4 * X1 * X2**3
        + 2 * X2**4
        - 3 * X1 * Y1**2
        - 3 * Y1**2 * X2
        - 6 * X1 * Y1 * Y2
        - 6 * Y1 * X2 * Y2
        - 3 * X1 * Y2**2
        - 3 * X2 * Y2**2,
    )  # = X1**4+X2**4+(X1+X2)**4- 3*(X1+X2)*(Y1+Y2)**2
    # 2

    f4 = (
        -4 * X1**4
        + 4 * X1**3 * X2
        - 9 * X1**2 * X2**2
        - 4 * X1 * B
        + 4 * X2 * B,
        -(X1**3)
        + 3 * X1**2 * X2
        - 3 * X1 * X2**2
        + X2**3
        - Y1**2
        + 2 * Y1 * Y2
        - Y2**2,
    )  # = (X2-X1)**3 - (Y1-Y2)**2
    # 2

    f5 = (
        -4 * X1**2 * X2**2 - X2**4 + 4 * X1 * B,
        X1**3 - 2 * X1**2 * X2 + X1 * X2**2 - Y1**2 + 2 * Y1 * Y2 - Y2**2,
    )  # = X1*(X1-X2)**2 - (Y1-Y2)**2
    # 2

    f6 = -X1 + X2 + 1, -X1 + X2 + 1
    # 4

    f7 = -X1 + X2 + 2, -X1 + X2 + 2
    # 4

    f8 = -(X1**3) * X2**3 - X1**3 * B - X2**3 * B + 8 * B**2, Y1 * Y2 - 3 * B
    # 4

    f9 = -(X1**3) * X2**3 - X1**3 * B - X2**3 * B + 8 * B**2, Y1 * Y2 + 3 * B
    # 4

    f10 = -(X1**2) * X2**2 + X1 * B + X2 * B, X2 * Y1 + X1 * Y2
    # 2

    f11 = (
        X1**3 * X2**3 - X1**2 * X2**2 + X1**3 * B + X2**3 * B + B**2,
        X1 * X2 + Y1 * Y2,
    )
    # 4

    if not extended:
        return {
            f"f{i+1}": f
            for i, f in enumerate([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11])
        }

    f12 = X1**2 + X1 * X2 + X2**2, Y1 + Y2
    f13 = X1**2 + X1 * X2 + X2**2, Y1 - Y2
    return {
        f"f{i+1}": f
        for i, f in enumerate([f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13])
    }


def load_efd_secp256k1_registers():
    secp256k1_polynomials = get_secp256k1_polynomials()

    registers_dict = {}

    """jacobian-0:madd, jacobian-0:add-1998-cmo, jacobian-0:add-2007-bl, jacobian-0:add-1998-cmo-2, jacobian-0:mmadd-2007-bl, jacobian-0:madd-2007-bl, jacobian-0:madd-2008-g, projective:add-1998-cmo, projective:add-1998-cmo-2, projective:madd-1998-cmo, projective:mmadd-1998-cmo, jacobian:madd, jacobian:add-1998-cmo, jacobian:add-2007-bl, jacobian:add-1998-cmo-2, jacobian:mmadd-2007-bl, jacobian:madd-2007-bl, jacobian:madd-2008-g, xyzz:add-2008-s, xyzz:madd-2008-s, xyzz:mmadd-2008-s, modified:madd-2009-bl, modified:add-1998-cmo-2, modified:add-2009-bl"""
    """represented by jacobian:add-1998-cmo"""
    registers = Registers()
    registers.add_tuple(*secp256k1_polynomials["f4"])
    registers_dict["jacobian:add-1998-cmo"] = registers

    """jacobian-0:zadd-2007-m, jacobian:zadd-2007-m"""
    registers = Registers()
    registers.add_tuple(*secp256k1_polynomials["f5"])
    registers_dict["jacobian:zadd-2007-m"] = registers

    """jacobian-0:madd-2004-hmv, jacobian:madd-2004-hmv"""
    registers = Registers()
    registers.add_tuple(*secp256k1_polynomials["f2"])
    registers_dict["jacobian:madd-2004-hmv"] = registers

    """jacobian-0:add-1986-cc, jacobian-0:add-2001-b, jacobian:add-1986-cc, jacobian:add-2001-b"""
    registers = Registers()
    registers.add_tuple(*secp256k1_polynomials["f1"])
    registers_dict["jacobian:add-1986-cc"] = registers

    """projective:add-2002-bj, projective:add-2007-bl"""
    registers = Registers()
    registers.add_tuple(*secp256k1_polynomials["f1"])
    registers.add_tuple(*secp256k1_polynomials["f3"])
    registers_dict["projective:add-2002-bj"] = registers

    """modified:mmadd-2009-bl"""
    registers = Registers()
    registers.add_tuple(*secp256k1_polynomials["f4"])
    registers.add_tuple(*secp256k1_polynomials["f6"])
    registers.add_tuple(*secp256k1_polynomials["f7"])
    registers_dict["modified:mmadd-2009-bl"] = registers

    """projective:madd-2015-rcb"""
    # f1, f8, f9, f10, f11
    registers = Registers()
    registers.add_tuple(*secp256k1_polynomials["f1"])
    registers.add_tuple(*secp256k1_polynomials["f8"])
    registers.add_tuple(*secp256k1_polynomials["f9"])
    registers.add_tuple(*secp256k1_polynomials["f10"])
    registers.add_tuple(*secp256k1_polynomials["f11"])
    registers_dict["projective:madd-2015-rcb"] = registers

    """our selection"""

    return registers_dict


class GLVCurve:
    def __init__(self):
        self.curve = None
        self.order = None
        self.lam = None
        self.beta = None
        self._order_field = None

    @property
    def field(self):
        return self.curve.base_field()

    @property
    def order_field(self):
        if self._order_field is None:
            self._order_field = GF(self.order)
        return self._order_field

    def generate_random(self, bits):
        self.curve = glv_module.find_curve_random(bits)
        self._compute_properties()

    def set_secp256k1(self):
        self.curve = glv_module.secp256k1["curve"]
        self._compute_properties()

    def generate(self, bits):
        self.curve = glv_module.find_curve(bits)
        self._compute_properties()

    def _compute_properties(self):
        self.lam, self.beta = map(ZZ, glv_module.find_lambda(self.curve))
        self.order = self.curve.order()
        self._order_field = GF(self.order)

    def to_dict(self):
        params = {}
        params["curve"] = [int(self.curve.a4()), int(self.curve.a6())]
        params["field"] = int(self.field.order())
        params["beta"] = int(self.beta)
        params["lambda"] = int(self.lam)
        params["order"] = int(self.order)
        return params

    def from_dict(self, dict):
        self.curve = EllipticCurve(GF(dict["field"]), dict["curve"])
        self.lam = ZZ(dict["lambda"])
        self.beta = ZZ(dict["beta"])
        self.order = ZZ(dict["order"])
        self._order_field = GF(self.order)

    @property
    def p(self):
        return self.field.order()

    @property
    def bits(self):
        return self.p.nbits()


class GLVSecrets:
    def __init__(self, glv=None):
        self.glv = GLVCurve() if glv is None else glv
        self.k = None
        self.k0 = None
        self.k1 = None

    @property
    def half_bits(self):
        return ceil(log(self.glv.order, 2) / 2)

    @property
    def bits(self):
        return self.glv.bits

    def generate_positive(self):
        self.k0 = ZZ(randint(2 ** (self.half_bits - 1), 2**self.half_bits))
        self.k1 = ZZ(randint(2 ** (self.half_bits - 1), 2**self.half_bits))
        self.k = (self.k0 + self.k1 * self.glv.lam) % self.glv.order

    def generate_positive_odd(self):
        self.k0 = (
            ZZ(randint(2 ** (self.half_bits - 2), 2 ** (self.half_bits - 1))) * 2 + 1
        )
        self.k1 = (
            ZZ(randint(2 ** (self.half_bits - 2), 2 ** (self.half_bits - 1))) * 2 + 1
        )
        self.k = (self.k0 + self.k1 * self.glv.lam) % self.glv.order

    def to_dict(self):
        dict = {"k0": int(self.k0), "k1": int(self.k1), "k": int(self.k)}
        dict.update(self.glv.to_dict())
        return dict

    def from_dict(self, dict):
        self.glv.from_dict(dict)
        self.k = ZZ(dict["k"])
        self.k0 = ZZ(dict["k0"])
        self.k1 = ZZ(dict["k1"])


class ZVPparams:
    def __init__(self, bits):
        self.target_bits = 0
        self.bits = bits
        self.secrets = GLVSecrets()

        self.verbose = True
        self.attack = None
        self.results = None
        self.registers = Registers()

    @property
    def glv(self):
        return self.secrets.glv

    def generate_secrets(self):
        self.glv.generate_random(self.bits)
        self.secrets.generate_positive_odd()

    def generate_secp256k1_secrets(self):
        self.glv.set_secp256k1()
        self.secrets.generate_positive_odd()

    def do_attack(self):
        self.results = self.attack(self)

    def results_to_dict(self):
        dict = self.secrets.to_dict()
        dict["target_bits"] = int(self.target_bits)
        dict["result"] = self.results
        dict["registers"] = self.registers.to_strings()
        return dict

    def error_log(self, exception):
        with open("error.log", "w") as f:
            f.write(str(self.secrets.to_dict()))
            f.write("\n\n\n")
            f.write(str(exception))

    @property
    def filename(self):
        if not self.target_bits:
            return f"{self.attack.__name__}_{self.bits}"
        return f"{self.attack.__name__}_{self.bits}_{self.target_bits}"


### Loading results


def extract_result_input(result):
    tmp = deepcopy(result)
    del tmp["result"]
    return tmp


def register_match(registers, registers_str):
    if not len(registers.polynomials_tuples) == len(registers_str):
        return False
    for [f, g] in registers_str:
        if not (registers.ring(f), registers.ring(g)) in registers.polynomials_tuples:
            return False
    return True


def register_submatch(registers, registers_str):
    for [f, g] in registers_str:
        if not (registers.ring(f), registers.ring(g)) in registers.polynomials_tuples:
            return False
    return True


def load_results(zvpparams, register_match_f, distinct):
    path = "results"
    results = []
    assert zvpparams.bits and zvpparams.attack
    for filename in os.listdir(path):
        if not filename.startswith(zvpparams.filename):
            continue
        with open(os.path.join(path, filename)) as handle:
            results.extend(json.load(handle))
    if zvpparams.registers.polynomials:
        results = [
            result
            for result in results
            if register_match_f(zvpparams.registers, result["registers"])
        ]

    if distinct:
        distinct_results = []
        distinct_inputs = []
        for result in results:
            tmp = extract_result_input(result)
            if tmp in distinct_inputs:
                continue
            distinct_results.append(result)
            distinct_inputs.append(tmp)
        results = distinct_results
    return results
