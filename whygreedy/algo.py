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
"""


def check_solution(solution: list[float], oxides: list[Compound], original: Compound):
    """
    check if the solution observes elemental conservation

    :param solution: a list of values representing x_i
    :param oxides: a list of compounds representing the decomposition product
    :param original: the reactant of decomposition reaction
    :return: bool
    """
    assert len(solution) == len(oxides)
    elements_in_oxides = []
    for oxi in oxides:
        elements_in_oxides += oxi.elements
    elements_in_oxides = set(elements_in_oxides)
    elements_to_check = sorted(set(original.elements).intersection(elements_in_oxides))
    checks = dict()
    for e in elements_to_check:
        v_original = original.normalized_formula[e]
        v_oxides = 0
        for x, oxide in zip(solution, oxides):
            try:
                v_oxides += x * oxide.normalized_formula[e]
            except KeyError:
                continue
        checks[e] = abs(v_original - v_oxides) < 1e-7
    return all(checks.values())


def calculate_ranking_parameter(oxide: Compound, original: Compound, for_oxide=True) -> float:
    """
    this is the ranking parameter as defined in the paper, a smaller number indicated the oxide is more favored,
    the one with the smallest ranking parameter will be used to "consume" the original compound
    """
    # if the oxide has an element with > 0 composition,
    # and this element is not present or is of 0 composition in the
    # original compound,
    # then the ranking parameter is inf because it is impossible to consume the original compound
    # with this oxide (and the solution for this oxide is 0)
    if for_oxide:
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
    else:
        for e in oxide.elements:
            oxide_composition = oxide.normalized_formula[e]
            original_composition = original.normalized_formula[e]
            if not is_close_to_zero(oxide_composition) and is_close_to_zero(original_composition):
                return np.inf
        p = 0
        for element in oxide.elements:
            c = oxide.normalized_formula[element] / original.normalized_formula[element]
            p += c
        rp = oxide.formation_energy_per_atom / p
        return rp


def find_greedy(
        oxides: list[Compound], original: Compound, first_choice: int = 0, diligent_greedy=True, for_oxide=True,
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
    sorted_oxides = sorted(sorted_oxides,
                           key=lambda x: calculate_ranking_parameter(x, updated_original, for_oxide=for_oxide))

    while len(solution) < len(oxides):
        if diligent_greedy:
            # greedy means to find the best in each iteration
            # the implementation found on [zenodo](https://zenodo.org/record/5110202#.YlJgpsjMJyg) does not sort the
            # oxide list at every iteration (only at initialization), so it is not a greedy algorithm
            # this becomes even more problematic considering they used different `first_choice`
            sorted_oxides = sorted(sorted_oxides,
                                   key=lambda x: calculate_ranking_parameter(x, updated_original, for_oxide=for_oxide))
        # we can force the first choice to be something else, but always choose the best starting the 2nd iteration
        if counter == 0:
            favored_oxide = sorted_oxides[first_choice]
            index_to_pop = first_choice
        else:
            favored_oxide = sorted_oxides[0]
            index_to_pop = 0
        # once the favored oxide is identified, we can calculate the ratio, and subtract it from the original compound
        ratio, updated_original = compound_subtract(favored_oxide, updated_original, for_oxide=for_oxide)
        # remove the favored oxide from ranking
        sorted_oxides.pop(index_to_pop)
        # update solution
        solution.append((favored_oxide.properties["index"], ratio))
        final_enthalpy += ratio * favored_oxide.formation_energy_per_atom
        # if the original has been consumed, fill in solution and break the loop
        if all(is_close_to_zero(v, 1e-7) for v in updated_original.normalized_formula.values()):
            for remaining_oxide in sorted_oxides:
                solution.append((remaining_oxide.properties["index"], 0.0))
            break

        # update counter before next iteration
        counter += 1
    solution = sorted(solution, key=lambda x: x[0])
    assert len(solution) == len(oxides)
    return [s[1] for s in solution], final_enthalpy - original.formation_energy_per_atom


def find_greedy_old(oxides: list[Compound], original: Compound, first_choice: int = 0, for_oxides=True) -> Tuple[
    list[float], float]:
    """
    This is just a wrapper for the implementation from 10.1021/acs.chemmater.1c02644
    It is identical to `find_greedy` with `diligent_greedy` set to False
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

    if for_oxides:
        solution_oxides, final_enthalpy, delta_enthalpy, _, _, _ = find_comp(stable_oxides, compound_unit_cell,
                                                                             original.formation_energy_per_atom,
                                                                             "Oxide", first_choice)
    else:
        solution_oxides, final_enthalpy, delta_enthalpy, _, _, _ = find_comp(stable_oxides, compound_unit_cell,
                                                                             original.formation_energy_per_atom, "non",
                                                                             first_choice)

    solution = np.zeros(len(oxides))
    for oxide in solution_oxides:
        solution[oxide["index"]] = oxide["ratio"]
    return solution.tolist(), delta_enthalpy


def find_lp(oxides: list[Compound], original: Compound) -> Tuple[list[float], float]:
    if len(oxides) == 0:
        return [], - original.formation_energy_per_atom

    elements_in_oxides = []
    for oxi in oxides:
        elements_in_oxides += oxi.elements
    elements_in_oxides = set(elements_in_oxides)
    elements_in_constraints = sorted(set(original.elements).intersection(elements_in_oxides))

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
            for e in elements_in_constraints:
                element_sum = 0
                for i, x_i in enumerate(x):
                    try:
                        composition = oxides[i].normalized_formula[e]
                    except KeyError:
                        continue
                    element_sum += composition * x_i
                c_e = m.addConstr(element_sum == original.normalized_formula[e], name=e)
                c.append(c_e)

            # objective function
            objective = 0
            for x_i, oxi in zip(x, oxides):
                objective += x_i * oxi.formation_energy_per_atom

            m.setObjective(objective, GRB.MINIMIZE)
            m.optimize()
            return [v.x for v in m.getVars()], m.objVal - original.formation_energy_per_atom
