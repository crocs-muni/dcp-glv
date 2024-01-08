import utils.formula as fm
from abc import ABC, abstractmethod


class Operation(ABC):
    def __init__(self, formula: fm.Formula):
        self.formula = formula
        self.output = formula.output()
        self.ring = formula.ring
        self.coordinates = formula.coordinates
        self.curve_field = formula.curve_field


class Add(Operation):
    def __init__(self, formula: fm.Formula):
        super().__init__(formula)
        assert formula.operation == "add"

    def __call__(self, P, Q):
        Xmap, Ymap = self.output["X"], self.output["Y"]

        X1, Y1, X2, Y2 = [
            self.ring.fraction_field().gens_dict().get(i, 1)
            for i in ["X1", "Y1", "X2", "Y2"]
        ]
        Xres = Xmap(X1=P[0], Y1=P[1], X2=Q[0], Y2=Q[1])
        Yres = Ymap(X1=P[0], Y1=P[1], X2=Q[0], Y2=Q[1])

        if self.coordinates.startswith("xyzz"):
            ZZmap, ZZZmap = self.output["ZZ"], self.output["ZZZ"]
            ZZres, ZZZres = ZZmap(X1=P[0], Y1=P[1], X2=Q[0], Y2=Q[1]), ZZZmap(
                X1=P[0], Y1=P[1], X2=Q[0], Y2=Q[1]
            )
            assert ZZres**3 == ZZZres**2
            return Xres / ZZres, Yres / ZZZres
        Zmap = self.output["Z"]
        Zres = Zmap(X1=P[0], Y1=P[1], X2=Q[0], Y2=Q[1])
        if self.coordinates.startswith("jacobian") or self.coordinates.startswith(
            "modified"
        ):
            return Xres / Zres**2, Yres / Zres**3
        if self.coordinates.startswith("projective"):
            return Xres / Zres, Yres / Zres
        if self.coordinates.startswith("w12"):
            return Xres / Zres, Yres / Zres**2
        if self.coordinates.startswith("extended"):
            return Xres / Zres, Yres / Zres
        if self.coordinates.startswith("inverted"):
            return Zres / Xres, Zres / Yres


class Dadd(Operation):
    def __init__(self, formula: fm.Formula):
        super().__init__(formula)
        assert formula.operation == "dadd"

    def __call__(self, P, Q, S):
        if self.coordinates.startswith("xz"):
            Xmap, Zmap = self.output["X"], self.output["Z"]
            X1, X2, X3 = [
                self.ring.fraction_field().gens_dict()[i] for i in ["X1", "X2", "X3"]
            ]
            Xres = Xmap(X1=S[0], X2=P[0], X3=Q[0])
            Zres = Zmap(X1=S[0], X2=P[0], X3=Q[0])
            return Xres / Zres
        if self.coordinates.startswith("yz"):
            Ymap, Zmap = self.output["Y"], self.output["Z"]
            Y1, Y2, Y3 = [
                self.ring.fraction_field().gens_dict()[i] for i in ["Y1", "Y2", "Y3"]
            ]
            Yres = Ymap(Y1=S[1], Y2=P[1], Y3=Q[1])
            Zres = Zmap(Y1=S[1], Y2=P[1], Y3=Q[1])
            if self.coordinates == "yz":
                return Yres / Zres / self.formula.r
            if self.coordinates == "yzsquared":
                return self.curve_field(Yres / Zres / self.formula.r).sqrt()


class Ladd(Operation):
    def __init__(self, formula: fm.Formula):
        super().__init__(formula)
        self.secondary_output = formula.secondary_output()
        assert formula.operation == "ladd"

    def __call__(self, P, Q, S):
        if self.coordinates.startswith("xz"):
            Xmap, Zmap, Xmap2, Zmap2 = (
                self.output["X"],
                self.output["Z"],
                self.secondary_output["X"],
                self.secondary_output["Z"],
            )
            X1, X2, X3 = [
                self.ring.fraction_field().gens_dict()[i] for i in ["X1", "X2", "X3"]
            ]
            Xres = Xmap(X1=S[0], X2=P[0], X3=Q[0])
            Zres = Zmap(X1=S[0], X2=P[0], X3=Q[0])
            X2res = Xmap2(X1=S[0], X2=P[0], X3=Q[0])
            Z2res = Zmap2(X1=S[0], X2=P[0], X3=Q[0])
            return Xres / Zres, X2res / Z2res

        if self.coordinates.startswith("yz"):
            Ymap, Zmap, Y2map, Z2map = (
                self.output["Y"],
                self.output["Z"],
                self.secondary_output["Y"],
                self.secondary_output["Z"],
            )
            Y1, Y2, Y3 = [
                self.ring.fraction_field().gens_dict()[i] for i in ["Y1", "Y2", "Y3"]
            ]
            Yres = Ymap(Y1=S[1], Y2=P[1], Y3=Q[1])
            Zres = Zmap(Y1=S[1], Y2=P[1], Y3=Q[1])
            Y2res = Y2map(Y1=S[1], Y2=P[1], Y3=Q[1])
            Z2res = Z2map(Y1=S[1], Y2=P[1], Y3=Q[1])
            if self.coordinates == "yz":
                return Yres / Zres / self.formula.r, Y2res / Z2res / self.formula.r
            if self.coordinates == "yzsquared":
                return (
                    self.curve_field(Yres / Zres / self.formula.r).sqrt(),
                    self.curve_field(Y2res / Z2res / self.formula.r).sqrt(),
                )


class Dbl(Operation):
    def __init__(self, formula: fm.Formula):
        super().__init__(formula)
        assert formula.operation == "dbl"

    def __call__(self, P):
        if self.coordinates.startswith("yz"):
            Y1 = [self.ring.fraction_field().gens_dict()[i] for i in ["Y1"]]
            Yres = self.output["Y"](Y1=P[1])
            Zres = self.output["Z"](Y1=P[1])
            if self.coordinates == "yz":
                return (Yres / Zres / self.formula.r,)
            if self.coordinates == "yzsquared":
                return ((Yres / Zres / self.formula.r).sqrt(),)

        X1 = self.ring.fraction_field().gens_dict()["X1"]
        if self.coordinates.startswith("xz"):
            return self.output["X"](X1=P[0]) / self.output["Z"](X1=P[0])

        Y1 = [self.ring.fraction_field().gens_dict()[i] for i in ["Y1"]]
        Xres = self.output["X"](X1=P[0], Y1=P[1])
        Yres = self.output["Y"](X1=P[0], Y1=P[1])
        if self.coordinates.startswith("xyzz"):
            ZZmap, ZZZmap = self.output["ZZ"], self.output["ZZZ"]
            ZZres, ZZZres = ZZmap(X1=P[0], Y1=P[1]), ZZZmap(X1=P[0], Y1=P[1])
            assert ZZres**3 == ZZZres**2
            return Xres / ZZres, Yres / ZZZres

        Zres = self.output["Z"](X1=P[0], Y1=P[1])
        if self.coordinates.startswith("jacobian"):
            return Xres / Zres**2, Yres / Zres**3
        if self.coordinates.startswith("modified"):
            return Xres / Zres**2, Yres / Zres**3
        if self.coordinates.startswith("projective"):
            return Xres / Zres, Yres / Zres
        if self.coordinates.startswith("w12"):
            return Xres / Zres, Yres / Zres**2
        if self.coordinates.startswith("inverted"):
            return Zres / Xres, Zres / Yres
        if self.coordinates.startswith("extended"):
            return Xres / Zres, Yres / Zres
