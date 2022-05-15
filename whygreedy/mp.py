import os.path
from collections import OrderedDict

import tqdm

from whygreedy import json_load, Compound

this_dir = os.path.dirname(os.path.abspath(__file__))
mpdata = os.path.join(this_dir, "../data/mp.json.gz")


def load_mp() -> list[dict]:
    clean_data = []
    data = json_load(mpdata)
    neclude = 0
    for compound in data:
        if None not in compound.values():
            clean_data.append(compound)
        else:
            neclude += 1
    print("exclude data as None in required fields: {}".format(neclude))
    discrepancy = 0
    for compound in clean_data:
        if compound['nsites'] != sum(compound['unit_cell_formula'].values()):
            discrepancy += 1
    assert discrepancy == 0
    return clean_data


def find_stable_compounds(all_compounds, criteria=50):
    stable_phase = []
    # find all compounds with e_above_hull within 0.05 eV
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


def find_oxide_pairs_from_compounds(compounds: list[Compound]):
    oxides = []
    non_oxides = []
    chemsys_to_oxide_list = dict()
    for c in compounds:
        if c.is_oxide:
            oxides.append(c)
            if len(c.elements_exclude_oxygen) == 0:
                continue  # ignore pure oxygen
            try:
                chemsys_to_oxide_list[frozenset(c.elements_exclude_oxygen)].append(c)
            except KeyError:
                chemsys_to_oxide_list[frozenset(c.elements_exclude_oxygen)] = [c, ]
        else:
            non_oxides.append(c)

    compound_of_no_oxides = []
    pairs = []
    for non_oxide in tqdm.tqdm(non_oxides):
        oxide_list = []
        non_oxide_element_set = set(non_oxide.elements)
        for element_fset in chemsys_to_oxide_list:
            if non_oxide_element_set.issuperset(element_fset):
                oxide_list += chemsys_to_oxide_list[element_fset]
        if len(oxide_list) == 0:
            compound_of_no_oxides.append(non_oxide)
            continue
        pairs.append((non_oxide, oxide_list))
    print("# of compound that has no oxides", len(compound_of_no_oxides))
    for c in compound_of_no_oxides:
        print(c)
    return pairs


def load_mp_decomposition_pairs():
    """
    mpid -> chemsys -> all possible subset chemsys -> all mpid
    """
    compounds = load_mp()
    compounds = [mpdata_to_compound(c) for c in compounds]

    mpid_to_compound = {c.mpid: c for c in compounds}

    chemsys_to_mpids = OrderedDict()
    mpid_to_chemsys = OrderedDict()
    print("create mpid2chemsys...")
    for c in tqdm.tqdm(compounds):
        try:
            chemsys_to_mpids[frozenset(c.elements)].append(c.mpid)
        except KeyError:
            chemsys_to_mpids[frozenset(c.elements)] = [c.mpid, ]
        mpid_to_chemsys[c.mpid] = frozenset(c.elements)

    chemsys_sets = [cs for cs in chemsys_to_mpids]

    chemsys_to_subsets = OrderedDict()
    print("create chemsys2subsets...")
    for i in tqdm.tqdm(range(len(chemsys_sets))):
        chemsys_i = chemsys_sets[i]
        subsets = []
        for j in range(len(chemsys_sets)):
            chemsys_j = chemsys_sets[j]
            if chemsys_i.issuperset(chemsys_j):
                subsets.append(chemsys_j)
        chemsys_to_subsets[chemsys_i] = subsets

    pairs = []
    print("create pairs...")
    for c in tqdm.tqdm(compounds):
        c_chemsys = mpid_to_chemsys[c.mpid]
        competing_phases = []
        for subset_chemsys in chemsys_to_subsets[c_chemsys]:
            competing_phases += chemsys_to_mpids[subset_chemsys]
        pairs.append((c, [mpid_to_compound[cpid] for cpid in competing_phases if cpid != c.mpid]))
    return pairs


def load_mp_oxidation_pairs():
    mp_data = load_mp()
    compounds = find_stable_compounds(mp_data, 50)
    print("stable compounds:", len(compounds))
    compounds = [mpdata_to_compound(c) for c in compounds]
    pairs = find_oxide_pairs_from_compounds(compounds)
    print("# of pairs loaded:", len(pairs))
    return pairs
