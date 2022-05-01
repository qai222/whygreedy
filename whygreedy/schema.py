import inspect
from itertools import combinations
from typing import Tuple

import numpy as np
from monty.json import MSONable


def normalize_stoi(formula_dictionary: dict[str, int]) -> dict[str, float]:
    natoms = sum(formula_dictionary.values())
    return {k: v / natoms for k, v in formula_dictionary.items()}


def is_close_to_zero(f: float, eps=1e-5):
    return abs(f) < eps


class Compound(MSONable):

    def __init__(
            self, normalized_formula: dict[str, float], formation_energy_per_atom: float,
            mpid: str = None, properties: dict = None
    ):
        self.normalized_formula = normalized_formula
        self.formation_energy_per_atom = formation_energy_per_atom
        self.mpid = mpid
        if properties is None:
            properties = dict()
        self.properties = properties

    def __repr__(self):
        return self.normalized_formula.__repr__()

    @property
    def elements(self):
        return sorted(self.normalized_formula.keys())

    @property
    def elements_exclude_oxygen(self):
        return [e for e in self.elements if e != "O"]

    @property
    def is_oxide(self):
        return "O" in self.elements

    @classmethod
    def from_dict(cls, d: dict):
        keys = [k for k in inspect.signature(Compound.__init__).parameters.keys() if
                k != "self" and not k.startswith("_")]
        data = dict()
        for k in keys:
            if k == "normalized_formula" and k in d:
                data[k] = normalize_stoi(d[k])
            else:
                try:
                    data[k] = d[k]
                except KeyError:
                    continue
        return cls(**data)

    @classmethod
    def random_compound(cls, elements: list[str], seed: int = None):
        rs = np.random.RandomState(seed)
        composition = rs.rand(len(elements)).tolist()
        return cls(
            normalized_formula=normalize_stoi(dict(zip(elements, composition))),
            formation_energy_per_atom=rs.uniform(-1e3, -1e-5),
        )


def is_oxidation_pair(oxide: Compound, original: Compound):
    return set(original.elements).issuperset(set(oxide.elements_exclude_oxygen))


def is_competing_pair(cp: Compound, original: Compound):
    return set(original.elements).issuperset(set(cp.elements))


def compound_subtract(oxide: Compound, original: Compound, for_oxide=True) -> Tuple[float, Compound]:
    # note this will update the original
    if for_oxide:
        assert is_oxidation_pair(oxide, original)
        ratios = [original.normalized_formula[x] / oxide.normalized_formula[x] for x in oxide.elements_exclude_oxygen]
        ratio = min(ratios)
        for e in oxide.elements_exclude_oxygen:
            original.normalized_formula[e] -= ratio * oxide.normalized_formula[e]
        return ratio, original
    else:
        assert is_competing_pair(oxide, original)
        ratios = [original.normalized_formula[x] / oxide.normalized_formula[x] for x in oxide.elements]
        ratio = min(ratios)
        for e in oxide.elements:
            original.normalized_formula[e] -= ratio * oxide.normalized_formula[e]
        return ratio, original


def gen_random_data(elements: list[str], num_oxi_per_chemical_system: int, seed: int) -> Tuple[
    Compound, list[Compound]]:
    element_combinations = []
    for i in range(1, len(elements) + 1):
        element_combinations += list(combinations(elements, i))
    original = Compound.random_compound(elements, seed=seed)
    oxides = []
    for element_combination in element_combinations:
        for _ in range(num_oxi_per_chemical_system):
            oxide = Compound.random_compound(list(element_combination) + ["O", ], seed=len(oxides) + seed)
            oxides.append(oxide)
    return original, oxides
