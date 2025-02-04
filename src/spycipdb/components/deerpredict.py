# pylint: skip-file
# flake8: noqa
"""Helper script to run DEERPREdict calculations for the `pre` module."""
import os

import numpy as np
import pandas as pd


try:
    from DEERPREdict.PRE import PREpredict
    from MDAnalysis import Universe
except ModuleNotFoundError:
    # Error message handled in `cli_pre.py`
    pass


default_deerpredict = {
    "atom":"H",
    "temp":298,
    "tau_c":2*1e-9,
    "tau_t":0.5*1e-9,
    "delay":10e-3,
    "r_2":10,
    "wh":750
    }


def calc_pre_predict(pdb, fexp, parameters):
    """Back calculates PRE data intensity ratios based on DEERPREdict."""
    ratios = []
    pdb_name = os.path.splitext(os.path.basename(pdb))[0]
    current_directory = os.getcwd()
    
    exp = pd.read_csv(fexp)
    res1 = exp.res1.values.astype(int)
    unique_res1 = np.unique(res1, return_index=True)[1]
    unique_res1 = res1[np.sort(unique_res1)]
    res2 = exp.res2.values.astype(int)
    res2_splitted = [[] for _ in range(len(unique_res1))]

    for key, value in zip(res1, res2):
        index = np.where(unique_res1 == key)[0][0]
        res2_splitted[index].append(value)

    atom = parameters["atom"]
    temp = parameters["temp"]
    tau_c = parameters["tau_c"]
    tau_t = parameters["tau_t"]
    delay = parameters["delay"]
    r_2 = parameters["r_2"]
    wh = parameters["wh"]
    
    u = Universe(pdb)
    
    for i, res in enumerate(unique_res1):
        pre = PREpredict(
            u,
            residue = res,
            log_file = pdb_name + "_log",
            temperature = temp,
            atom_selection = atom
            )
        pre.run(
            output_prefix = pdb_name,
            tau_c = tau_c,
            tau_t = tau_t,
            delay = delay,
            r_2 = r_2,
            wh = wh
            )
        
        results_dat = f"{current_directory}/{pdb_name}-{res}.dat"
        
        with open(results_dat, 'r') as results_f:
            lines = results_f.readlines()
            lines.pop(0)
            for line in lines:
                splitted = line.split()

                res_num = int(eval(splitted[0]))
                try:
                    ratio = float(eval(splitted[1]))
                except NameError:
                    ratio = float('nan')
                if res_num in res2_splitted[i]:
                    ratios.append(ratio)
        
        try:
            os.remove(results_dat)
            os.remove(f"{current_directory}/{pdb_name}-{res}.pkl")
            os.remove(f"{current_directory}/{pdb_name}-Z-{res}.dat")
            os.remove(f"{current_directory}/{pdb_name}_log")
        except FileNotFoundError:
            pass
    return pdb, ratios
    