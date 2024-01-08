import glv
from sage.all import random_prime, GF, PolynomialRing, ZZ, randint, ceil
import matplotlib.pyplot as plt
from pari_tools.dcp_pari import dcpsolver_pari2
import numpy as np
import json

FLOAT_FORMAT = "{0:0.2f}"


def gen_lam(bits, dim):
    n = random_prime(2**bits - 1, True, 2 ** (bits - 1))
    F = GF(n)
    x = PolynomialRing(F, ["x"]).gen()
    while True:
        N_bits = bits // dim
        T_bits = N_bits // 2 + 1
        T = ZZ(randint(2 ** (T_bits - 1), 2**T_bits))
        N = ZZ(randint(2 ** (N_bits - 1), 2**N_bits))
        try:
            lam, _ = (x**2 - T * x + N).roots()[0]
        except:
            continue
        break
    return lam, T, N, n


def exp_max_scalar_size(bits, dim):
    n_exp = 100
    n_exp2 = 100
    avg_max = 0
    for _ in range(n_exp):
        lam, _, _, n = gen_lam(bits, dim)
        m = glv.decomposition_matrix(n, lam, dim)
        for _ in range(n_exp2):
            k = randint(1, n)
            kis = glv.decompose_precomputed(k, lam, m, dim, n)
            avg_max += max([ki.nbits() for ki in kis])
    print(
        f"Bits:{bits}, dim:{dim}, avg max bits: {FLOAT_FORMAT.format(avg_max/(n_exp*n_exp2))}, expected:{FLOAT_FORMAT.format(bits/dim)}"
    )


def shifted_scalar_experiments(bits, dim):
    bits, dim = ZZ(bits), ZZ(dim)
    lam, _, _, n = gen_lam(bits, dim)
    zeroscalars = ceil((dim - 1) / 3)
    n_exp = 10**7
    ks = []
    for _ in range(n_exp):
        k = ZZ(sum([randint(0, n) * lam**i for i in range(dim - 1 - zeroscalars)]))
        ks.append(float(k / n))
    return ks


def plot_distributions(path):
    with open(path) as h:
        dim_ks = json.load(h)
    plt.style.use("seaborn-deep")
    fig, _ = plt.subplots()
    ts = []
    labs = []
    for dim, ks in dim_ks:
        t = np.array(ks, dtype=float)
        ts.append(t)
        labs.append(f"r={dim}")
    fig.set_figwidth(15)
    plt.hist(ts, bins=10, label=labs)
    plt.xlabel("scalars", fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.legend(fontsize=20)
    plt.savefig("text/skewed_histogram.png", bbox_inches="tight")


from scipy.stats import chisquare


def chisquare_test(path):
    with open(path) as h:
        dim_ks = json.load(h)
    nbins = 0
    for dim, ks in dim_ks:
        maxk, mink = max(ks), min(ks)
        nbins = 10
        binsize = (maxk - mink) / nbins
        distribution = {i * binsize: 0 for i in range(1, nbins + 1)}
        for k in ks:
            for i in range(1, nbins + 1):
                if k < i * binsize:
                    distribution[i * binsize] += 1
                    break
        s, p = chisquare(list(distribution.values()))
        print(f"dimension: {dim}, chi-squared test statistic: {s}, p-value: {p}")


def lexperiment():
    bits = ZZ(60)
    dim = ZZ(12)
    lam, _, _, n = gen_lam(bits, dim)
    lam = ZZ(lam)
    m = glv.decomposition_matrix(n, lam, dim)
    zeroscalars = 1
    counter = 0
    n_exp = 1000
    for _ in range(n_exp):
        scalar = ZZ(randint(0, n))
        for l in range(2**2, 2 ** (ceil(bits / dim) + 2)):
            kis = glv.decompose_precomputed(l * scalar, lam, m, dim, n)
            tail = sum([kis[-i] == 0 for i in range(1, zeroscalars + 1)])
            counter += tail == zeroscalars
    print(counter / n_exp)


def count_solution_dcp(curve, scalar):
    A, B = curve.a4(), curve.a6()
    F = curve.base_field()
    p = F.characteristic()
    roots = dcpsolver_pari2(p, A, B, scalar)
    counter = 0
    for root in [F(int(r)) for r in roots]:
        try:
            P = curve.lift_x(root)
            assert P[0] + (scalar * P)[0] + 2 == 0
            counter += 1
        except ValueError as e:
            continue
    return counter


def exp_alpha_probability(bits, k_bits):

    n_exp = 10**2
    root_dist = {}
    for _ in range(n_exp):
        curve = glv.find_curve(bits)
        k = randint(2 ** (k_bits - 1), 2**k_bits)
        nroots = count_solution_dcp(curve, k)
        if not nroots in root_dist:
            root_dist[nroots] = 0
        root_dist[nroots] += 1
    print(root_dist)


if __name__ == "__main__":

    # Experiments for Proposition 1
    # for dim in range(2,8):
    #    exp_max_scalar_size(128,dim)
    # for dim in range(2,8):
    #    exp_max_scalar_size(256,dim)

    # Experiments for Proposition 2
    # path = "results/skewed_scalars.json"
    # chisquare_test(path)
    # with open(path,"w") as h:
    #    json.dump([(dim,shifted_scalar_experiments(128,dim)) for dim in range(3,8)], h)
    # plot_distributions(path)

    # Experiments for Lemma 3
    # for k_bits in range(3,8):
    #    exp_alpha_probability(32,k_bits)
    pass
