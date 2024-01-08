import urllib.request, json, os
from sage.all import GF
import utils.curve as cv


def load_curves():
    source = "https://dissect.crocs.fi.muni.cz/"
    path_to_curves = "curves.json"

    if not os.path.isfile(path_to_curves):
        download_curves(source)
    curves = load_curves_from_file(path_to_curves)
    sorted_curves = sort_curves_by_form(curves)
    saged_curves = parse_to_sage_field(sorted_curves)
    return saged_curves


def download_curves(source):
    req = urllib.request.Request(f"{source}db/curves", method="GET")
    with urllib.request.urlopen(req) as f:
        curves = json.loads(f.read())["data"]
    curves = filter(
        lambda c: c["field"]["type"] == "Prime" and c["standard"] == True, curves
    )
    curves = list(curves)
    with open("curves.json", "w") as handle:
        json.dump(curves, handle)


def load_curves_from_file(path_to_curves):
    with open(path_to_curves) as handle:
        return json.load(handle)


def sort_curves_by_form(curves):
    filtered_curves = {
        "Weierstrass": [],
        "Montgomery": [],
        "TwistedEdwards": [],
        "Edwards": [],
    }
    for curve in curves:
        filtered_curves[curve["form"]].append(curve)
    return filtered_curves


def parse_weierstrass(curve):
    field = GF(curve["field"]["p"])
    a, b = field(curve["params"]["a"]["raw"]), field(curve["params"]["b"]["raw"])
    return cv.ShortWeierstrass(a=a, b=b, name=curve["name"])


def parse_montgomery(curve):
    field = GF(curve["field"]["p"])
    a, b = field(curve["params"]["a"]["raw"]), field(curve["params"]["b"]["raw"])
    return cv.Montgomery(a=a, b=b, name=curve["name"])


def parse_edwards(curve):
    field = GF(curve["field"]["p"])
    try:
        c = field(curve["params"]["c"]["raw"])
    except KeyError:
        c = field(1)
    d = field(curve["params"]["d"]["raw"])
    return cv.Edwards(c=c, d=d, name=curve["name"])


def parse_twisted(curve):
    field = GF(curve["field"]["p"])
    a, d = field(curve["params"]["a"]["raw"]), field(curve["params"]["d"]["raw"])
    return cv.TwistedEdwards(a=a, d=d, name=curve["name"])


def edwards_to_twisted(edwards):
    return cv.TwistedEdwards(a=edwards.field(1), d=edwards.d, name=edwards.name)


def twisted_to_edwards(twisted):
    return cv.Edwards(c=twisted.field(1), d=twisted.d, name=twisted.name)


def parse_to_sage_field(curves):
    result = {}
    result["shortw"] = list(map(parse_weierstrass, curves["Weierstrass"]))
    result["montgom"] = list(map(parse_montgomery, curves["Montgomery"]))

    twisted_curves = []
    edwards_curves = []
    for curve in curves["Edwards"]:
        edwards = parse_edwards(curve)
        if edwards.c == 1:
            twisted_curves.append(edwards_to_twisted(edwards))
        edwards_curves.append(edwards)

    for curve in curves["TwistedEdwards"]:
        twisted = parse_twisted(curve)
        if twisted.a == 1:
            edwards_curves.append(twisted_to_edwards(twisted))
        twisted_curves.append(twisted)

    result["twisted"] = twisted_curves
    result["edwards"] = edwards_curves
    return result


if __name__ == "__main__":
    load_curves()
