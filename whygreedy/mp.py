import os.path
import pandas as pd
import seaborn as sns

import numpy as np
import tqdm

from whygreedy import json_load, Compound, find_lp, find_greedy, find_greedy_old

this_dir = os.path.dirname(os.path.abspath(__file__))
mpdata = os.path.join(this_dir, "../data/mp.json.gz")


def load_mp():
    clean_data = []
    data = json_load(mpdata)
    for compound in data:
        if None not in compound.values():
            clean_data.append(compound)
    discrepancy = 0
    for compound in clean_data:
        if compound['nsites'] != sum(compound['unit_cell_formula'].values()):
            discrepancy += 1
    assert discrepancy == 0
    return clean_data


def find_stable_phases(all_compounds, criteria=50):
    stable_phase = []
    # find all compounds with e_above_hull within 0.05 of 0
    for compound in all_compounds:
        if abs(compound['e_above_hull']) < criteria / 1000:
            stable_phase.append(compound)
    return stable_phase


def mpdata_to_compound(d: dict):
    data = {
        "normalized_formula": d["unit_cell_formula"],
        "formation_energy_per_atom": d["formation_energy_per_atom"],
        "mpid": d["task_id"]
    }
    return Compound.from_dict(data)


def find_pairs(compounds: list[Compound]):
    chemical_system_to_oxides = dict()
    for c in compounds:
        if c.is_oxide:
            key = frozenset(c.elements_exclude_oxygen)
            try:
                chemical_system_to_oxides[key].append(c)
            except KeyError:
                chemical_system_to_oxides[key] = [c]
    pairs = []
    for c in compounds:
        if not c.is_oxide:
            key = frozenset(c.elements)
            try:
                oxides = chemical_system_to_oxides[key]
            except KeyError:
                continue
            pairs.append((c, oxides))
    return pairs


def set_small_to_zeros(a: list[float], eps=1e-5):
    a = np.array(a)
    a[np.abs(a) < eps] = 0
    return a


def load_mp_pairs():
    mp_data = load_mp()
    compounds = find_stable_phases(mp_data, 50)
    compounds = [mpdata_to_compound(c) for c in compounds]
    pairs = find_pairs(compounds)
    return pairs


def compare_greedy_exact(pairs, greedy_type="old"):
    dh_diff = []
    quali_diff = []
    solutions_old = []
    solutions_exact = []
    for original, oxides in tqdm.tqdm(pairs):

        sol_old_min, dh_old_min = None, np.inf
        for i in range(len(oxides)):
            if greedy_type == "old":
                sol_old, dh_old = find_greedy_old(oxides, original, first_choice=i)
            elif greedy_type == "new":
                sol_old, dh_old = find_greedy(oxides, original, first_choice=i, real_greedy=False)
            elif greedy_type == "real":
                sol_old, dh_old = find_greedy(oxides, original, first_choice=i, real_greedy=True)
            else:
                raise NotImplementedError
            if dh_old < dh_old_min:
                dh_old_min = dh_old
                sol_old_min = sol_old

        sol_exact, dh_exact = find_lp(oxides, original)
        dh_diff.append(dh_exact - dh_old_min)

        sol_old_min = set_small_to_zeros(sol_old_min)
        sol_exact = set_small_to_zeros(sol_exact)
        quali_diff.append(not np.allclose(sol_old_min.astype(bool), sol_exact.astype(bool)))

        solutions_old.append(sol_old_min)
        solutions_exact.append(sol_exact)
    dh_diff = set_small_to_zeros(dh_diff)
    return dh_diff, np.array(quali_diff), solutions_old, solutions_exact


def hist_compare(dh_diff, quali_diff, dh_diff_threshold=0.0):
    df = pd.DataFrame()

    mask = np.abs(dh_diff) >= dh_diff_threshold

    df["Delta H difference (eV/atom)"] = dh_diff[mask]
    df["qualitative difference"] = quali_diff[mask]
    histplot = sns.histplot(df, x="Delta H difference (eV/atom)", hue="qualitative difference", hue_order=[True, False], multiple="stack", binwidth=0.1)
    fig = histplot.get_figure()
    return fig


if __name__ == '__main__':
    pairs = load_mp_pairs()
    dh_diff, quali_diff, _, _ = compare_greedy_exact(pairs, greedy_type="old")

    fig = hist_compare(dh_diff, quali_diff, 0)
    fig.savefig("compare_above0.00.png", dpi=300)
    fig.clf()

    fig = hist_compare(dh_diff, quali_diff, 0.05)
    fig.savefig("compare_above0.05.png", dpi=300)
    fig.clf()

    print("total pairs:", len(dh_diff))
    mask = np.abs(dh_diff) > 0.05
    print("dh diff above 0.05:", len(dh_diff[mask]))
    print("qualitative diff:", sum(quali_diff))
    """
    total pairs: 9238
    dh diff above 0.05: 1212
    qualitative diff: 1371
    """