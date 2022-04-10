from copy import deepcopy
from typing import Tuple

import gurobipy as gp
import numpy as np
from gurobipy import GRB

from whygreedy.Twyman2022ChemMat import find_comp
from whygreedy.schema import Compound, is_close_to_zero, compound_subtract

"""
implement the greedy algorithm proposed by the authors of 10.1021/acs.chemmater.1c02644
note what they proposed is different from what they actually implemented, see the docs of function `find_greedy`

The oxidation enthalpy minimization problem is:
Given an oxygen-free compound [original], find the solution X={x_1, x_2, ...} s.t. the free energy change of reaction
[original] + y[oxygen] --> x1[oxide1] + x2[oxide2] + ... is minimized, where
1. `oxide1, oxide2, ...` are potential oxides
2. y is arbitrary
3. elemental conversation is observed.
"""


def check_solution(solution: list[float], oxides: list[Compound], original: Compound):
    assert len(solution) == len(oxides)
    elements = original.elements
    sum_oxides = dict(zip(elements, np.zeros(len(elements))))
    for x, oxide in zip(solution, oxides):
        for element, v in oxide.normalized_formula.items():
            if element == "O":
                continue
            sum_oxides[element] += v * x
    return all(sum_oxides[e] < original.normalized_formula[e] + 1e-5 for e in elements)


def calculate_ranking_parameter(oxide: Compound, original: Compound) -> float:
    """
    this is the ranking parameter as defined in the paper, a smaller number indicated the oxide is more favored,
    the one with the smallest ranking parameter will be used to "consume" the original compound
    """
    # if the oxide has an element with > 0 composition, and this element is not present or is of 0 composition in the
    # original compound, then the ranking parameter is inf because it is impossible to consume the original compound
    # with this oxide
    # (and the solution for this oxide is 0)
    for e in oxide.elements_exclude_oxygen:
        oxide_composition = oxide.normalized_formula[e]
        original_composition = original.normalized_formula[e]
        if not is_close_to_zero(oxide_composition) and is_close_to_zero(original_composition):
            return np.inf
    p = 0
    for element in oxide.elements:
        if element == "O":
            continue
        c = oxide.normalized_formula[element] / original.normalized_formula[element]
        p += c
    rp = oxide.formation_energy_per_atom / p
    return rp


def find_greedy(
        oxides: list[Compound], original: Compound, first_choice: int = 0, real_greedy=True
) -> Tuple[list[float], float]:
    # placeholder for the solution
    solution = []
    final_enthalpy = 0.0

    # assign an index for oxides
    for ioxi, oxi in enumerate(oxides):
        oxi.properties["index"] = ioxi

    # we will be updating the original compound, better make a deep copy
    updated_original = deepcopy(original)
    sorted_oxides = oxides  # no need for deep copy

    # init the loop and perform the first greedy ranking
    counter = 0
    sorted_oxides = sorted(sorted_oxides, key=lambda x: calculate_ranking_parameter(x, updated_original))

    while len(solution) < len(oxides):
        if real_greedy:
            # greedy means to find the best in each iteration
            # the implementation found on [zenodo](https://zenodo.org/record/5110202#.YlJgpsjMJyg) does not sort the
            # oxide list at every iteration (only at initialization), so it is not a greedy algorithm
            # this becomes even more problematic considering they used different `first_choice`
            sorted_oxides = sorted(sorted_oxides, key=lambda x: calculate_ranking_parameter(x, updated_original))
        # we can force the first choice to be something else, but always choose the best starting the 2nd iteration
        if counter == 0:
            favored_oxide = sorted_oxides[first_choice]
            index_to_pop = first_choice
        else:
            favored_oxide = sorted_oxides[0]
            index_to_pop = 0
        # once the favored oxide is identified, we can calculate the ratio, and subtract it from the original compound
        ratio, updated_original = compound_subtract(favored_oxide, updated_original)
        # remove the favored oxide from ranking
        sorted_oxides.pop(index_to_pop)
        # update solution
        solution.append((favored_oxide.properties["index"], ratio))
        final_enthalpy += ratio * favored_oxide.formation_energy_per_atom
        # update counter before next iteration
        counter += 1
    solution = sorted(solution, key=lambda x: x[0])
    assert len(solution) == len(oxides)
    return [s[1] for s in solution], final_enthalpy - original.formation_energy_per_atom


def find_greedy_old(oxides: list[Compound], original: Compound, first_choice: int = 0) -> Tuple[
    list[float], float]:
    """
    This is just a wrapper for the implementation from 10.1021/acs.chemmater.1c02644
    It is identical to `find_greedy` with `real_greedy` set to False
    """
    stable_oxides = []
    for ioxide, oxide in enumerate(oxides):
        stable_oxide = {
            "nsites": 1,
            # have to make a copy as `find_comp` changes it
            "unit_cell_formula": {k: v for k, v in oxide.normalized_formula.items()},
            "formation_energy_per_atom": oxide.formation_energy_per_atom,
            "index": ioxide,
            "elements": oxide.elements,
        }
        stable_oxides.append(stable_oxide)
    compound_unit_cell = original.normalized_formula

    solution_oxides, final_enthalpy, delta_enthalpy, _, _, _ = find_comp(stable_oxides, compound_unit_cell,
                                                                         original.formation_energy_per_atom, "Oxide",
                                                                         first_choice)
    solution = np.zeros(len(oxides))
    for oxide in solution_oxides:
        solution[oxide["index"]] = oxide["ratio"]
    return solution.tolist(), delta_enthalpy


def find_lp(oxides: list[Compound], original: Compound) -> Tuple[list[float], float]:
    # init gurobi model, suppress output
    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)
        env.setParam('LogToConsole', 0)
        env.start()
        with gp.Model(env=env) as m:

            # add variables
            x = []
            for ioxi, oxi in enumerate(oxides):
                x_i = m.addVar(name=str(ioxi))
                x.append(x_i)

            # add constraints
            c = []
            for e in original.elements:
                element_sum = 0
                for i, x_i in enumerate(x):
                    try:
                        composition = oxides[i].normalized_formula[e]
                    except KeyError:
                        composition = 0
                    element_sum += composition * x_i
                c_e = m.addConstr(element_sum <= original.normalized_formula[e], name=e)
                c.append(c_e)

            # objective function
            objective = 0
            for x_i, oxi in zip(x, oxides):
                objective += x_i * oxi.formation_energy_per_atom

            m.setObjective(objective, GRB.MINIMIZE)
            m.optimize()
            return [v.x for v in m.getVars()], m.objVal - original.formation_energy_per_atom
