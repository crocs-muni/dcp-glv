from sage.all import ZZ, ceil, log, GF, EllipticCurve 
import glv
import zvp_glv_msm
import zvp_glv_sm
import twostep
import bsgssolver as bsgs
from random import randint
from datetime import datetime
import json, os


FLOAT_FORMAT = "{0:0.2f}"

# =========================================================
# ================== Parameter generation =================
# =========================================================


def gen_glv_curve(bits, random=True):
    """Finds suitable curve for 2-dim decomposition"""
    curve = {}
    if random:
        curve["curve"] = glv.find_curve_random(bits)
    else:
        curve["curve"] = glv.find_curve(bits)
    curve["field"] = curve["curve"].base_field()
    curve["lam"], curve["beta"] = glv.find_lambda(curve["curve"])
    curve["order"] = curve["curve"].order()
    return curve


def gen_secret_scalar(bits, lam, order):
    """unused, negative scalars cause problems"""
    k = ZZ(randint(2 ** (bits - 1), 2**bits))
    k0, k1 = glv.glv_decompose_simple(k, lam, order)
    return k, k0, k1


def gen_secret_scalar_positive(order, lam):
    """Generates private key (integer mod order), with positive scalar decomposition"""
    bits = order.nbits()
    k0, k1 = ZZ(randint(2 ** (bits // 2 - 1), 2 ** (bits // 2))), ZZ(
        randint(2 ** (bits // 2 - 1), 2 ** (bits // 2))
    )
    k = (k0 + k1 * lam) % order
    return k, k0, k1


def params_setup(bits):
    params = gen_glv_curve(bits)
    params["k"], params["k0"], params["k1"] = gen_secret_scalar_positive(
        params["order"], params["lam"]
    )
    return params


def gen_secret_scalar_positive_2(order, lam):
    """Generates private key (integer mod order), with positive scalar decomposition"""
    bits = ceil(log(order, 2) / 2)
    k0 = ZZ(randint(2 ** (bits - 2), 2 ** (bits - 1))) * 2 + 1
    k1 = ZZ(randint(2 ** (bits - 1), 2 ** (bits)))
    k = (k0 + k1 * lam) % order
    return k, k0, k1


def params_setup_2(bits):
    params = gen_glv_curve(bits)
    params["k"], params["k0"], params["k1"] = gen_secret_scalar_positive_2(
        params["order"], params["lam"]
    )
    return params

def params_to_dict(params):
    curve = params["curve"]
    field = params["field"]
    params["curve"] = [int(curve.a4()), int(curve.a6())]
    params["field"] = int(field.order())
    for key in "lam", "beta", "order", "k", "k0", "k1":
        params[key] = int(params[key])
    return params 
   
def dict_to_params(params):
    curve = params["curve"]
    field = GF(params["field"])
    params["field"]=field 
    params["curve"] = EllipticCurve(field,[curve[0],curve[1]])
    for key in "lam", "beta", "order", "k", "k0", "k1":
        params[key] = ZZ(params[key])
    return params 


def random_twostep_msm_parameters(bits,N = 10):
    with open(f"results/params_twostep_msm_{bits}.json","w") as h:
        json.dump([params_to_dict(params_setup_2(bits)) for _ in range(N)],h)

def random_twostep_sm_parameters(bits,N = 10):
    with open(f"results/params_twostep_sm_{bits}.json","w") as h:
        json.dump([params_to_dict(params_setup(bits)) for _ in range(N)],h)


# =========================================================
# ====================== Attacks  =========================
# =========================================================


def random_zvp_glv_sm(bits, zvp_target, verbose=True):
    params = params_setup(bits)
    params["zvp_target"] = zvp_target
    return zvp_glv_sm.zvp_glv_sm(params, verbose)


def random_zvp_sm(bits, zvp_target, verbose=True):
    params = params_setup(bits)
    params["zvp_target"] = zvp_target
    return zvp_glv_sm.zvp_sm(params, verbose)


def random_zvp_glv_msm(bits, zvp_target, verbose=True):
    params = params_setup_2(bits)
    params["zvp_target"] = zvp_target
    return zvp_glv_msm.zvp_glv_msm(params, verbose)


def random_twostep_sm(bits, zvp_target, verbose=True):
    params = params_setup(bits)
    params["zvp_target"] = zvp_target
    return twostep.twostep_sm(params, verbose)


def random_twostep_msm(bits, zvp_target, verbose=True):
    params = params_setup_2(bits)
    params["zvp_target"] = zvp_target
    return twostep.twostep_msm(params, verbose)

# =========================================================
# ===================== Experiments =======================
# =========================================================


def save_to_file(results, filename):
    for result in results:
        curve = result["params"]["curve"]
        field = result["params"]["field"]
        result["params"]["curve"] = [int(curve.a4()), int(curve.a6())]
        result["params"]["field"] = int(field.order())
        for key in "lam", "beta", "order", "k", "k0", "k1", "zvp_target":
            result["params"][key] = int(result["params"][key])
    if not os.path.exists("results"):
        os.makedirs("results")
    with open(
        f"results/{filename}_{datetime.today().strftime('%Y-%m-%d_%H:%M')}.json",
        "w",
    ) as f:
        json.dump(results, f)

def zvp_glv_sm_exp(bits, n_exp, c_bounds, verbose=False):
    print("zvp glv sm")
    for zvp_target in map(ZZ, range(*c_bounds)):
        results = []
        for i in range(n_exp):
            print(f"bits:{bits}, zvp_target: {zvp_target}, {i + 1}/{n_exp}")
            results.append(random_zvp_glv_sm(ZZ(bits), zvp_target, verbose))
        save_to_file(results, f"zvp_glv_sm_{zvp_target}")


def zvp_glv_msm_exp(bits, n_exp, c_bounds, verbose=False):
    print("zvp glv msm")
    for zvp_target in map(ZZ, range(*c_bounds)):
        results = []
        for i in range(n_exp):
            print(f"bits:{bits}, zvp_target: {zvp_target}, {i + 1}/{n_exp}")
            results.append(random_zvp_glv_msm(ZZ(bits), zvp_target, verbose))
        save_to_file(results, f"zvp_glv_msm_{zvp_target}")


def zvp_sm_exp(bits, n_exp, c_bounds, verbose=False):
    print("zvp sm")
    for zvp_target in map(ZZ, c_bounds):
        results = []
        for i in range(n_exp):
            print(f"bits:{bits}, zvp_target: {zvp_target}, {i + 1}/{n_exp}")
            results.append(random_zvp_sm(ZZ(bits), zvp_target, verbose))
        save_to_file(results, f"zvp_sm_{zvp_target}")




def twostep_sm_exp_random(params,n_exp,offset=None,verbose=False):
    print("twostep sm random")
    offset = [0] if offset is None else offset
    for bits,zvp_target in params:
        for offset in offset:
            results = []
            for i in range(n_exp):
                print(f"bits: {bits}, zvp target: {zvp_target+offset}, {i + 1}/{n_exp}")
                results.append(random_twostep_sm(ZZ(bits), ZZ(zvp_target+offset), verbose))
            save_to_file(results, f"twostep_sm_{bits}_{zvp_target+offset}")

def twostep_msm_exp_random(params,n_exp,offset=None,verbose=False):
    print("twostep msm random")
    offset = [0] if offset is None else offset
    for bits,zvp_target in params:
        for offset in offset:
            results = []
            for i in range(n_exp):
                print(f"bits: {bits}, zvp target: {zvp_target+offset}, {i + 1}/{n_exp}")
                results.append(random_twostep_msm(ZZ(bits), ZZ(zvp_target+offset), verbose))
            save_to_file(results, f"twostep_msm_{bits}_{zvp_target+offset}")

def twostep_sm_exp(bits,zvp_target,verbose=False):
    filename = f"params_twostep_sm_{bits}.json"
    print(f"twostep sm from {filename}")
    with open("results/"+filename) as f:
        param_list = json.load(f)
    results = []
    for i,params_dict in enumerate(param_list):
        params = dict_to_params(params_dict)
        params["zvp_target"]=zvp_target
        print(f"bits: {bits}, zvp target: {zvp_target}, {i + 1}/{len(param_list)}")
        results.append(twostep.twostep_sm(params,verbose))
    save_to_file(results, f"twostep_sm_{bits}_{zvp_target}")

def twostep_sm_set_exp(bits,zvp_target,verbose=False):
    filename = f"params_twostep_sm_{bits}.json"
    print(f"twostep sm from {filename}")
    with open("results/"+filename) as f:
        param_list = json.load(f)
    results = []
    for i,params_dict in enumerate(param_list):
        params = dict_to_params(params_dict)
        params["zvp_target"]=zvp_target
        print(f"bits: {bits}, zvp target: {zvp_target}, {i + 1}/{len(param_list)}")
        results.append(twostep.twostep_sm_set(params,verbose))
    save_to_file(results, f"twostep_sm_{bits}_{zvp_target}")



def twostep_msm_exp(bits,zvp_target,verbose=False):
    filename = f"params_twostep_msm_{bits}.json"
    print(f"twostep msm from {filename}")
    with open("results/"+filename) as f:
        param_list = json.load(f)
    results = []
    for i,params_dict in enumerate(param_list):
        params = dict_to_params(params_dict)
        params["zvp_target"]=zvp_target
        print(f"bits: {bits}, zvp target: {zvp_target}, {i + 1}/{len(param_list)}")
        results.append(twostep.twostep_msm(params,verbose))
    save_to_file(results, f"twostep_msm_{bits}_{zvp_target}")


def bsgs_avg_timing():
    bit_params = {40:[2,3,4],44:[3,4,5],48:[3,4,5,6],52:[4,5,6,7],56:[4,5,6,7],60:[5,6,7,8],64:[6,7,8],68:[6,7,8],72:[7,8,9]}
    avg_times = {bits:{} for bits in bit_params}
    for bits,lbits_l in bit_params.items():
        for lbits in lbits_l:
            avg_times[bits][lbits] = bsgs.bsgs_timing(bits,bits//2-lbits,N=5)
    with open(f"results/bsbg_timing_{datetime.today().strftime('%Y-%m-%d_%H:%M')}.json","w") as h:
        json.dump(avg_times,h)

def msm_dcp_timing():
    bit_params = {40:5,44:5,48:6,52:7,56:7,60:8,64:8,68:8,72:9}
    timings = {}
    N = 10
    for bits,zvp_t in bit_params.items():
        avg_timing = [0 for _ in range(zvp_t)]
        for _ in range(N):
            params = params_setup_2(bits)
            params["zvp_target"] = zvp_t
            t = zvp_glv_msm.zvp_attack_timing(params)
            for i in range(zvp_t):
                avg_timing[i]+=t[i]
        avg_timing = list(map(lambda x: float(x/N),avg_timing))
        timings[bits]=avg_timing 
    with open(f"results/msm_dcp_timing_{datetime.today().strftime('%Y-%m-%d_%H:%M')}.json","w") as h:
        json.dump(timings,h) 


# =========================================================
# ===================== Printing results ==================
# =========================================================



def zvp_glv_sm_results(results,verbose):
    avg_time_zvp = 0
    avg_recovered = 0
    zvp_target = results[0]['params']['zvp_target']
    for result in results:
        recovered = result["zvp_recovered_0"]+result["zvp_recovered_1"]
        time_zvp = result["time_zvp"]
        avg_time_zvp += time_zvp
        avg_recovered += recovered
        if verbose: print(f'recovered:{recovered} out of {zvp_target*2}')
        if verbose: print(f'time:{time_zvp}')
    print(f"avg bits:{FLOAT_FORMAT.format(avg_recovered/len(results))}, avg time: {FLOAT_FORMAT.format(avg_time_zvp/len(results))}")
        

def zvp_glv_msm_results(results,verbose):
    avg_time_zvp = 0
    avg_recovered = 0
    zvp_target = results[0]['params']['zvp_target']
    for result in results:
        recovered = result["zvp_recovered"]
        time_zvp = result["time_zvp"]
        avg_time_zvp += time_zvp
        avg_recovered += recovered
        if verbose: print(f'recovered:{recovered} out of {zvp_target*2}')
        if verbose: print(f'time:{time_zvp}')
    print(f"avg bits:{FLOAT_FORMAT.format(avg_recovered/len(results))}, avg time: {FLOAT_FORMAT.format(avg_time_zvp/len(results))}")
    

def zvp_sm_results(results,verbose):
    avg_time_zvp = 0
    avg_recovered = 0
    zvp_target = results[0]['params']['zvp_target']
    for result in results:
        recovered = result["zvp_recovered"]
        time_zvp = result["time_zvp"]
        avg_time_zvp += time_zvp
        avg_recovered += recovered
        if verbose: print(f'recovered:{recovered} out of {zvp_target}')
        if verbose: print(f'time:{time_zvp}')
    print(f"avg bits:{FLOAT_FORMAT.format(avg_recovered/len(results))}, avg time: {FLOAT_FORMAT.format(avg_time_zvp/len(results))}")
    

def print_zvp_results(attack,zvp_target, verbose=False):
    path = "results/"
    results = []
    for filename in os.listdir(path):
        if not filename.startswith(f"{attack}_{zvp_target}"):
            continue
        with open(path+filename) as handle:
            results.extend(json.load(handle))
    print(f"{attack}, zvp target: {zvp_target}, results: {len(results)}")
    if attack=="zvp_glv_sm":
        zvp_glv_sm_results(results,verbose)
    elif attack=="zvp_glv_msm":
        zvp_glv_msm_results(results,verbose)
    elif attack=="zvp_sm":
        zvp_sm_results(results,verbose)
        

def print_twostep_results(attack,bits,zvp_target,verbose=False):
    path = "results/"
    results = []
    for filename in os.listdir(path):
        if not filename.startswith(f"{attack}_{bits}_{zvp_target}"):
            continue
        with open(path+filename) as handle:
            results.extend(json.load(handle))
    avg_time_zvp = 0
    avg_time_bsgs = 0
    print(f"{attack},bits={bits}, zvp target: {zvp_target}, results: {len(results)}")
    for i,result in enumerate(results):
        time_zvp = result["time_zvp"]
        time_bsgs = result["time_bsgs"]
        avg_time_zvp += time_zvp
        avg_time_bsgs += time_bsgs
        if verbose:
            print(f'{i}.time_zvp:{time_zvp}',end=", ")
            print(f'time_bsgs:{time_bsgs}')
            if "zvp_target_2" in result and result["zvp_target_2"]!=result["params"]["zvp_target"]:
                print(f'\t aborted from {result["params"]["zvp_target"]} to {result["zvp_target_2"]}')
            
    N = len(results)
    print(f"avg time_zvp:{FLOAT_FORMAT.format(avg_time_zvp/N)}, avg time_bsgs: {FLOAT_FORMAT.format(avg_time_bsgs/N)}, total:{float((avg_time_bsgs+avg_time_zvp)/N)}")

def summary_zvp():
    path = "results/"
    for attack in ["zvp_glv_sm","zvp_glv_msm","zvp_sm"]:
        zvps = set()
        for filename in os.listdir(path):
            if not filename.startswith(f"{attack}"):
                continue
            zvp = filename.split(attack)[1].split("_")[1]
            zvps.add(int(zvp))
        print(attack,"="*50)
        for zvp in sorted(list(zvps)):
            print_zvp_results(attack,zvp)
                
def summary_twostep():
    path = "results/"
    for attack in ["twostep_msm","twostep_sm"]:
        params = set()
        for filename in os.listdir(path):
            if not filename.startswith(f"{attack}"):
                continue
            bits,zvp = filename.split(attack)[1].split("_")[1:3]
            params.add((int(bits),int(zvp)))
        print(attack,"="*50)
        for (bits,zvp) in sorted((list(params))):
            print_twostep_results(attack,bits,zvp)
            print()

def compare_c(attack,bits,c1,c2):
    path = "results/"
    results = []
    for c in c1,c2:
        for filename in os.listdir(path):
            if not filename.startswith(f"{attack}_{bits}_{c}"):
                continue
            with open(f"{path}{filename}") as h:
                results.append(json.load(h))
    for i,(r1,r2) in enumerate(zip(*results)):
        t1 = r1["time_zvp"]+r1["time_bsgs"]
        t2 = r2["time_zvp"]+r2["time_bsgs"]
        if t1>t2:
            print(i,t1,t2,r1["params"]["field"])
 

if __name__ == "__main__":

    # ================= Running experiments =======================
    #zvp_glv_sm_exp(bits=128, n_exp = 10, c_bounds = (4,5), verbose=False)
    #zvp_glv_msm_exp(bits=128, n_exp = 10, c_bounds = (4,5), verbose=False)
    #zvp_sm_exp(bits, n_exp = 10, c_bounds = (4,5), verbose=False)
    #bsgs_avg_timing()
    #msm_dcp_timing()
    #for bits in range(40,73,4):
    #    random_twostep_msm_parameters(bits,N = 30)
    #for bits in range(40,73,4):
    #    random_twostep_sm_parameters(bits,N = 30)
    #twostep_msm_exp(40,4)
    #twostep_sm_exp(40,4)
    
    # ================= Printing results ==========================
    #print_zvp_results("zvp_glv_sm",8)
    #print_zvp_results("zvp_glv_msm",8,verbose=True)
    #print_zvp_results("zvp_sm",8)
    #print_twostep_results("twostep_sm",60,7,True)
    #print_twostep_results("twostep_msm",56,6,True)
    #compare_c("twostep_msm",44,5,4)
    #summary_zvp()
    #summary_twostep()
    pass
