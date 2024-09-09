import numpy
from sage.all import ZZ, log
import zvp_glv_sac
import json
import dcp
from utils import ZVPparams
import utils
from utils import load_results
import matplotlib.pyplot as plt
import itertools
import zvp_glv_inter_easy_prec as inter_easy
import zvp_ltr_signed
from matplotlib.ticker import PercentFormatter


def uncertainity_bits(poly, w):
    with open(f"results/interleaving_secp256k1_remapped_{w}.json") as f:
        points = json.load(f)
    register = utils.Registers()
    register.add_tuple(*poly)
    for register_points in points:
        if utils.register_submatch(register, register_points["register"]):
            avg = 0
            for [_, scalars] in register_points["points"]:
                avg += len(scalars)
            avg /= len(register_points["points"])
    return float(avg)


def zvp_interleaving_heat_row(zvpparams, polynomial_selection, bits_selection):

    poly_labels = {}
    secp256k1_poly = utils.get_secp256k1_polynomials()
    zvpparams.registers.empty_out()
    ring = zvpparams.registers.ring
    for name in polynomial_selection:
        f = secp256k1_poly[name]
        zvpparams.registers.add_tuple(*f)
        poly_labels[f] = name

    row_rec = []
    row_time = []
    for bits in bits_selection:
        zvpparams.target_bits = bits
        results = load_results(zvpparams, utils.register_submatch, distinct=True)
        # collect_scalars
        scalars = {}
        for result in results:
            scalar_tuple = (result["k0"], result["k1"], result["k"])
            if not scalar_tuple in scalars:
                scalars[scalar_tuple] = {}
            registers = result["registers"]
            assert len(registers) == 1
            poly_tuple = (ring(registers[0][0]), ring(registers[0][1]))
            recovered = result["result"]["recovered"]
            time = result["result"]["time_zvp"]
            poly_name = poly_labels[poly_tuple]
            scalars[scalar_tuple][poly_name] = (time, recovered)
        time = 0
        recovered = 0
        for scalar, res in scalars.items():
            assert len(res) == len(polynomial_selection), (res, polynomial_selection)
            for poly in polynomial_selection:
                time += res[poly][0]
                if res[poly][1] == 1:
                    recovered += float(
                        max(
                            ZZ(0),
                            2 * bits
                            - log(uncertainity_bits(secp256k1_poly[poly], bits), 2),
                        )
                    )
                    break
        avg_rec = recovered / len(scalars)
        avg_time = time / len(scalars)
        row_rec.append(avg_rec)
        row_time.append(avg_time)
    return row_time, row_rec


def zvp_interleaving_heat_single():
    zvpparams = ZVPparams(bits=256)
    zvpparams.attack = dcp.interleaving_dcp_experiment
    bit_selection = [3, 4, 5]
    grid_list_t = []
    grid_list_r = []
    xlabels = []

    selections = [
        ["f10", "f11", "f1", "f8", "f9"],
        ["f4", "f6", "f7"],
        ["f3", "f1"],
        ["f4"],
    ]
    selections_tex = [
        "proj:madd-2015-rcb(5)",
        "mod:mmadd-2009-bl(3)",
        "proj:add-2002-bj(2)",
        "jac:add-1998-cmo(1)",
    ]

    grid_list_r_labels = []
    for polynomial_selection in selections:
        row_t, row_r = zvp_interleaving_heat_row(
            zvpparams, polynomial_selection, bit_selection
        )
        row_t = [round(t, 1) for t in row_t]
        row_r_labels = [f"{round(r,1)}b" for r in row_r]
        row_r = [round(100 * r / (bit_selection[i] * 2)) for i, r in enumerate(row_r)]
        grid_list_t.append(row_t)
        grid_list_r.append(row_r)
        grid_list_r_labels.append(row_r_labels)

    xlabels = bit_selection
    ylabels = selections_tex

    xlabels = [f"{x}+{x}b\nw={x}" for x in xlabels]

    figsize, cmap, vmin, vmax, text_color = (8, 3), "inferno", 0, 4500, "white"
    heat(
        grid_list_t,
        xlabels,
        ylabels,
        f"{zvpparams.attack.__name__}",
        figsize,
        cmap,
        vmin,
        vmax,
        text_color,
    )
    figsize, cmap, vmin, vmax, text_color = (8, 3), "viridis", 0, 150, "white"
    heat(
        grid_list_r,
        xlabels,
        ylabels,
        f"{zvpparams.attack.__name__}_rec",
        figsize,
        cmap,
        vmin,
        vmax,
        text_color,
        grid_list_r_labels,
    )


def heat(
    double_list,
    xlabels,
    ylabels,
    filename,
    figsize,
    cmap,
    vmin,
    vmax,
    text_color,
    double_list_labels=None,
):
    nmpy_grid_list = numpy.array(double_list)
    fig, ax = plt.subplots(figsize=figsize)
    im = ax.imshow(
        nmpy_grid_list,
        aspect=0.5,
        extent=[-0.5, len(xlabels) - 0.5, len(ylabels) - 0.5, -0.5],
        vmin=vmin,
        vmax=vmax,
        cmap=cmap,
    )

    ax.set_xticks(numpy.arange(len(xlabels)), labels=xlabels)
    ax.set_yticks(numpy.arange(len(ylabels)), labels=ylabels)
    ax.tick_params(axis="both", which="major", labelsize=15)
    ax.tick_params(axis="both", which="minor", labelsize=15)

    plt.setp(ax.get_xticklabels(), rotation=0, ha="center", rotation_mode="anchor")

    for i in range(len(ylabels)):
        for j in range(len(xlabels)):
            if not double_list_labels:
                time = nmpy_grid_list[i, j]
                unit = "s"
                if time > 500:
                    time = round(time / 60, 1)
                    unit = "m"
                text = ax.text(
                    j,
                    i,
                    f"{time}{unit}",
                    ha="center",
                    va="center",
                    color=text_color,
                    fontsize="xx-large",
                )
            else:
                text = ax.text(
                    j,
                    i,
                    double_list_labels[i][j],
                    ha="center",
                    va="center",
                    color=text_color,
                    fontsize="xx-large",
                )

    fig.tight_layout()
    fig.savefig(f"graphs/{filename}.png", bbox_inches="tight", dpi=300)


def zvp_heat_row(zvpparams, polynomial_selection, bits_selection, value):
    poly_labels = {}
    secp256k1_poly = utils.get_secp256k1_polynomials()
    zvpparams.registers.empty_out()
    for name in polynomial_selection:
        f = secp256k1_poly[name]
        zvpparams.registers.add_tuple(*f)
        poly_labels[f] = name
    results = load_results(zvpparams, utils.register_match, distinct=True)
    tbits_time = {}
    for result in results:
        tbits = result["target_bits"]
        if not tbits in tbits_time:
            tbits_time[tbits] = (0, 0)
        tbits_time[tbits] = (
            tbits_time[tbits][0] + result["result"][value],
            tbits_time[tbits][1] + 1,
        )
    for tbits, v in tbits_time.items():
        tbits_time[tbits] = v[0] / v[1]
    row = []
    for tbits in bits_selection:
        row.append(tbits_time[tbits])
    return row


def zvp_heat(algorithm):
    zvpparams = ZVPparams(bits=256)

    if algorithm == "glv_sac":
        zvpparams.attack = zvp_glv_sac.zvp_glv_sac
        bit_selection = [2, 3, 4, 5]
        bit_potential = [4, 6, 8, 10]
        xlabels = [f"{x}+{x}b" for x in bit_selection]
    else:
        zvpparams.attack = zvp_ltr_signed.zvp_ltr_signed
        bit_selection = [4, 5, 6, 7]
        bit_potential = [4, 5, 6, 7]
        xlabels = [f"{x}b" for x in bit_selection]

    grid_list_time = []
    grid_list_rec = []
    grid_list_rec_labels = []

    selections = [
        ["f10", "f11", "f1", "f8", "f9"],
        ["f4", "f6", "f7"],
        ["f3", "f1"],
        ["f4"],
    ]
    selections_tex = [
        "proj:madd-2015-rcb(5)",
        "mod:mmadd-2009-bl(3)",
        "proj:add-2002-bj(2)",
        "jac:add-1998-cmo(1)",
    ]

    # selections = [["f10","f11","f1","f8","f9"],["f4","f6","f7"],["f4"]]
    # selections_tex = ["proj:madd-2015-rcb","mod:mmadd-2009-bl","jac:add-1998-cmo"]
    for polynomial_selection in selections:
        row_time = zvp_heat_row(
            zvpparams, polynomial_selection, bit_selection, "time_zvp"
        )
        row_time = [round(t, 1) for t in row_time]
        row_rec = zvp_heat_row(
            zvpparams, polynomial_selection, bit_selection, "recovered"
        )
        row_rec_lab = [f"{round(r,1)}b" for r in row_rec]
        row_rec = [100 * r / bit_potential[i] for i, r in enumerate(row_rec)]

        grid_list_time.append(row_time)
        grid_list_rec.append(row_rec)
        grid_list_rec_labels.append(row_rec_lab)
    ylabels = selections_tex

    figsize, cmap, vmin, vmax, text_color = (8, 3), "inferno", 0, 4500, "white"
    heat(
        grid_list_time,
        xlabels,
        ylabels,
        f"{zvpparams.attack.__name__}",
        figsize,
        cmap,
        vmin,
        vmax,
        text_color,
    )
    figsize, cmap, vmin, vmax, text_color = (8, 3), "viridis", 0, 150, "white"
    heat(
        grid_list_rec,
        xlabels,
        ylabels,
        f"{zvpparams.attack.__name__}_rec",
        figsize,
        cmap,
        vmin,
        vmax,
        text_color,
        grid_list_rec_labels,
    )


def heat_map_dcp_timings(algorithm):
    assert algorithm in ["dcp", "dcpglv"]
    zvpparams = ZVPparams(bits=256)
    poly_labels = {}

    selection = ["f2", "f5", "f6", "f9", "f10", "f11"]
    secp256k1_poly = utils.get_secp256k1_polynomials()
    for name in selection:
        f = secp256k1_poly[name]
        zvpparams.registers.add_tuple(*f)
        poly_labels[f] = name

    if algorithm == "dcpglv":
        zvpparams.attack = dcp.glv_dcp_experiment
    else:
        zvpparams.attack = dcp.dcp_experiment

    poly_vector = zvpparams.registers.polynomials_tuples

    poly_tbits_time = {poly_labels[poly]: {} for poly in poly_vector}

    for poly in poly_vector:
        zvpparams.registers.empty_out()
        zvpparams.registers.add_tuple(*poly)
        results = load_results(zvpparams, utils.register_match, distinct=True)
        tbits_time = poly_tbits_time[poly_labels[poly]]
        for result in results:
            tbits = result["target_bits"]
            if not tbits in tbits_time:
                tbits_time[tbits] = (0, 0)
            tbits_time[tbits] = (
                tbits_time[tbits][0] + result["result"]["time_zvp"],
                tbits_time[tbits][1] + 1,
            )
        for tbits, v in tbits_time.items():
            tbits_time[tbits] = float(v[0] / v[1])

    grid_list = []
    xlabels = []
    for _, row in poly_tbits_time.items():
        tbits = sorted(list(row.keys()))
        grid_list.append([round(row[t], 1) for t in tbits])
        xlabels.extend(tbits)
    xlabels = sorted(list(set(xlabels)))
    ylabels = ["$f_2$", "$f_5$", "$f_6$", "$f_9$", "$f_{10}$", "$f_{11}$"]

    if algorithm == "dcpglv":
        xlabels = [f"{x}+{x}b" for x in xlabels]
    else:
        xlabels = [f"{x}b" for x in xlabels]

    figsize, cmap, vmin, vmax, text_color = (10, 4), "inferno", -200, 2000, "white"
    heat(
        grid_list,
        xlabels,
        ylabels,
        f"{zvpparams.attack.__name__}",
        figsize,
        cmap,
        vmin,
        vmax,
        text_color,
    )


def scatter_recovery():
    zvpparams = ZVPparams(bits=256)
    secp256k1_poly = utils.get_secp256k1_polynomials()
    poly_labels = {}
    for name, f in secp256k1_poly.items():
        zvpparams.registers.add_tuple(*f)
        poly_labels[f] = name
    zvpparams.attack = dcp.glv_dcp_experiment
    poly_vector = zvpparams.registers.polynomials_tuples
    results = load_results(zvpparams, utils.register_submatch, distinct=True)
    scalar_dict = {}
    all_poly_tuples = set()
    poly_count = {}
    for result in results:
        registers = result["registers"]
        if len(registers) > 1:
            continue
        poly_tuple = tuple(registers[0])
        all_poly_tuples.add(poly_tuple)
        scalar_tuple = (result["k"], result["k0"], result["k1"])
        if not scalar_tuple in scalar_dict:
            scalar_dict[scalar_tuple] = []

        if not poly_tuple in poly_count:
            poly_count[poly_tuple] = 0
        poly_count[poly_tuple] += 1

        if result["result"]["recovered"] == 1:
            scalar_dict[scalar_tuple].append(poly_tuple)

    # check
    assert len(set(poly_count.values())) == 1

    secp256k1_registers = utils.load_efd_secp256k1_registers()
    xaxis = []
    yaxis = []

    efdxaxis = []
    efdyaxis = []
    ring = poly_vector[0][0].parent()

    for i in range(1, 6):
        for poly_comb in list(itertools.combinations(poly_vector, i)):
            recovery = 0
            poly_set = set(poly_comb)
            for _, poly_tuples_rec in scalar_dict.items():
                poly_set2 = set([(ring(f), ring(g)) for f, g in poly_tuples_rec])
                recovery += len(poly_set2.intersection(poly_set)) > 0
            recovery /= len(scalar_dict)
            recovery *= 100
            xaxis.append(recovery)
            yaxis.append(i)

            for name, registers in secp256k1_registers.items():
                poly_set2 = set(
                    [(ring(f), ring(g)) for f, g in registers.polynomials_tuples]
                )
                if poly_set == poly_set2:
                    efdxaxis.append(recovery)
                    efdyaxis.append(i)
                    break

    xaxis = numpy.array(xaxis)
    yaxis = numpy.array(yaxis)
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.scatter(xaxis, yaxis, marker="o", color="steelblue")
    ax.scatter(efdxaxis, efdyaxis, color="red", marker="o")

    plt.xlabel("DCP solution rate in %", fontsize=12)
    plt.ylabel("Number of IV polynomials", fontsize=12)
    fig.tight_layout()
    fig.savefig(f"graphs/scatter.png", bbox_inches="tight", dpi=300)
    return


def interleaving_easy_success(zvpparams):
    results = utils.load_results(zvpparams, utils.register_match, distinct=True)
    rems = [result["result"]["remaining"] for result in results]
    return rems


def easy_interleaving_bar(window):
    zvpparams = ZVPparams(bits=256)
    # zvpparams.target_bits = 5
    secp256k1_polynomials = utils.get_secp256k1_polynomials()

    attack = {
        5: inter_easy.zvp_glv_interleaving_easy_regular_5,
        3: inter_easy.zvp_glv_interleaving_easy_regular_3,
        4: inter_easy.zvp_glv_interleaving_easy_regular_4,
    }[window]
    zvpparams.attack = attack

    efd = [
        "projective:madd-2015-rcb",
        "modified:mmadd-2009-bl",
        "projective:add-2002-bj",
        "jacobian:add-1998-cmo",
    ]
    short_efd = [
        "madd-2015-rcb(5)",
        "mmadd-2009-bl(3)",
        "add-2002-bj(2)",
        "add-1998-cmo(1)",
    ]
    fig, ax = plt.subplots(figsize=(13, 4.5))
    bins = [i * 10 for i in range(3, 26)]
    plt.xticks(ticks=[i * 10 for i in range(3, 26)])
    all_buckets = []
    for formula in efd:
        zvpparams.registers.empty_out()
        zvpparams.registers = utils.load_efd_secp256k1_registers()[formula]
        rems = interleaving_easy_success(zvpparams)
        # print(formula,sum([r<70 for r in rems]),sum(rems)/len(rems))
        all_buckets.append(rems)

    nonefd = [["f1", "f11", "f6", "f7", "f9"]]
    nonefdlabels = ["$f_1,f_{11},f_6,f_7,f_9$"]
    secp256k1_polynomials = utils.get_secp256k1_polynomials()
    for formula in nonefd:
        zvpparams.registers.empty_out()
        for f in formula:
            zvpparams.registers.add_tuple(*secp256k1_polynomials[f])
        rems = interleaving_easy_success(zvpparams)
        # print(formula,sum([r<70 for r in rems]),sum(rems)/len(rems))
        all_buckets.append(rems)

    colors = ["indianred", "yellowgreen", "steelblue", "hotpink", "gold"]
    # weights = [[1/100 for _ in bucket] for bucket in all_buckets]
    ax.hist(
        all_buckets,
        bins=bins,
        histtype="bar",
        label=short_efd + nonefdlabels,
        color=colors,
        rwidth=1,
    )
    plt.gca().yaxis.set_major_formatter(PercentFormatter())
    ax.legend(fontsize="x-large")
    ax.set_xlabel(
        "Size (in bits) of the space containing the private key", fontsize="large"
    )
    plt.savefig(f"graphs/{attack.__name__}", bbox_inches="tight", dpi=300)


if __name__ == "__main__":
    scatter_recovery()
    heat_map_dcp_timings("dcp")
    heat_map_dcp_timings("dcpglv")
    zvp_heat("glv_sac")
    zvp_heat("ltr_signed")
    zvp_interleaving_heat_single()
    easy_interleaving_bar(3)
    easy_interleaving_bar(4)
    easy_interleaving_bar(5)
    pass
