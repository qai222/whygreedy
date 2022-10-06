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


def check_solution(solution: list[float], products: list[Compound], reactant: Compound):
    """
    check if the solution observes elemental conservation

    :param solution: a list of values representing x_i
    :param products: a list of compounds representing the products
    :param reactant: the reactant of decomposition reaction
    :return: bool
    """
    assert len(solution) == len(products)
    elements_in_oxides = []
    for oxi in products:
        elements_in_oxides += oxi.elements
    elements_in_oxides = set(elements_in_oxides)
    elements_to_check = sorted(set(reactant.elements).intersection(elements_in_oxides))
    checks = dict()
    for e in elements_to_check:
        v_original = reactant.normalized_formula[e]
        v_oxides = 0
        for x, oxide in zip(solution, products):
            try:
                v_oxides += x * oxide.normalized_formula[e]
            except KeyError:
                continue
        checks[e] = abs(v_original - v_oxides) < 1e-7
    return all(checks.values())


def calculate_ranking_parameter(product: Compound, reactant: Compound, for_oxide: bool) -> float:
    """
    this is the ranking parameter as defined in the paper, a smaller number indicated the oxide is more favored,
    the one with the smallest ranking parameter will be used to "consume" the original compound
    """
    # if the product has an element with > 0 composition,
    # and this element is not present or is of 0 composition in the reactant,
    # then the ranking parameter is inf because it is impossible to consume the reactant with this product
    # (and x_i for this product is 0)
    if for_oxide:
        assert "O" not in reactant.normalized_formula, "you are calculating ranking param for an oxidation reaction, " \
                                                       "but your reactant has oxygen"
        check_elements = product.elements_exclude_oxygen
    else:
        check_elements = product.elements

    for e in check_elements:
        product_composition = product.normalized_formula[e]
        reactant_composition = reactant.normalized_formula[e]
        if not is_close_to_zero(product_composition) and is_close_to_zero(reactant_composition):
            return np.inf
    p = 0
    for element in check_elements:
        c = product.normalized_formula[element] / reactant.normalized_formula[element]  # the `formation cost`
        p += c
    rp = product.formation_energy_per_atom / p
    return rp


def find_greedy(
        reactant: Compound, products: list[Compound], first_choice: int, diligent_greedy: bool, for_oxide: bool,
) -> Tuple[list[float], float]:
    if len(products) == 0:
        return [], - reactant.formation_energy_per_atom

    # placeholder for the solution
    solution = []
    # sum of formation enthalpies of products
    final_enthalpy = 0.0

    # assign an index for oxides
    for iprod, prod in enumerate(products):
        prod.properties["index"] = iprod

    # we will be updating the original compound, better make a deep copy
    updated_reactant = deepcopy(reactant)
    sorted_products = products  # no need for deep copy

    # init the loop and perform the first greedy ranking
    counter = 0
    sorted_products = sorted(sorted_products,
                             key=lambda x: calculate_ranking_parameter(x, updated_reactant, for_oxide=for_oxide))

    while len(solution) < len(products):
        if diligent_greedy:
            # greedy means to find the best in each iteration
            # the implementation found on [zenodo](https://zenodo.org/record/5110202#.YlJgpsjMJyg) does not sort the
            # product list at every iteration (only at initialization),
            # so strictly speaking it is not a greedy algorithm
            # this becomes even more problematic considering they exhausted all possible `first_choice`
            sorted_products = sorted(sorted_products,
                                     key=lambda x: calculate_ranking_parameter(x, updated_reactant,
                                                                               for_oxide=for_oxide))
        # we can force the first choice to be something else, but always choose the best starting the 2nd iteration
        if counter == 0:
            favored_product = sorted_products[first_choice]
            index_to_pop = first_choice
        else:
            favored_product = sorted_products[0]
            index_to_pop = 0
        # once the favored product is identified, we can calculate the ratio,
        # and subtract it from the reactant
        ratio, updated_reactant = compound_subtract(favored_product, updated_reactant, for_oxide=for_oxide)
        # remove the favored oxide from ranking
        sorted_products.pop(index_to_pop)
        # update solution
        solution.append((favored_product.properties["index"], ratio))
        final_enthalpy += ratio * favored_product.formation_energy_per_atom
        # if the original has been consumed, fill in solution and break the loop
        if all(is_close_to_zero(v, 1e-7) for v in updated_reactant.normalized_formula.values()):
            for remaining_product in sorted_products:
                solution.append((remaining_product.properties["index"], 0.0))
            break

        # update counter before next iteration
        counter += 1
    solution = sorted(solution, key=lambda x: x[0])
    assert len(solution) == len(products)
    return [s[1] for s in solution], final_enthalpy - reactant.formation_energy_per_atom


def find_greedy_old(
        reactant: Compound, products: list[Compound], first_choice: int, for_oxide: bool,
) -> Tuple[list[float], float]:
    """
    This is just a wrapper for the implementation from 10.1021/acs.chemmater.1c02644
    It is identical to `find_greedy` with `diligent_greedy` set to False
    """
    if len(products) == 0:
        return [], - reactant.formation_energy_per_atom

    stable_products = []
    for iproduct, product in enumerate(products):
        stable_product = {
            "nsites": 1,
            # have to make a copy as `find_comp` changes it
            "unit_cell_formula": {k: v for k, v in product.normalized_formula.items()},
            "formation_energy_per_atom": product.formation_energy_per_atom,
            "index": iproduct,
            "elements": product.elements,
        }
        stable_products.append(stable_product)
    compound_unit_cell = reactant.normalized_formula

    if for_oxide:
        solution_oxides, final_enthalpy, delta_enthalpy, _, _, _ = find_comp(stable_products, compound_unit_cell,
                                                                             reactant.formation_energy_per_atom,
                                                                             "Oxide", first_choice)
    else:
        solution_oxides, final_enthalpy, delta_enthalpy, _, _, _ = find_comp(stable_products, compound_unit_cell,
                                                                             reactant.formation_energy_per_atom, "non",
                                                                             first_choice)

    solution = [0.0, ] * len(products)
    for product in solution_oxides:
        solution[product["index"]] = product["ratio"]
    return solution, delta_enthalpy


def find_lp(reactant: Compound, products: list[Compound], ) -> Tuple[list[float], float]:
    if len(products) == 0:
        return [], - reactant.formation_energy_per_atom

    elements_in_products = []
    for product in products:
        elements_in_products += product.elements
    elements_in_products = set(elements_in_products)
    elements_in_constraints = sorted(set(reactant.elements).intersection(elements_in_products))

    # init gurobi model, suppress output
    with gp.Env(empty=True) as env:
        env.setParam('OutputFlag', 0)
        env.setParam('LogToConsole', 0)
        env.start()
        with gp.Model(env=env) as m:

            # add variables
            x = []
            for iproduct, product in enumerate(products):
                x_i = m.addVar(name=str(iproduct))
                x.append(x_i)

            # add stoi constraints
            c = []
            for e in elements_in_constraints:
                element_sum = 0
                for i, x_i in enumerate(x):
                    try:
                        composition = products[i].normalized_formula[e]
                    except KeyError:
                        continue
                    element_sum += composition * x_i
                c_e = m.addConstr(element_sum == reactant.normalized_formula[e], name=e)
                c.append(c_e)
            # add non-negative constraints
            for i, x_i in enumerate(x):
                c_e = m.addConstr(x_i >= 0, name="x_{}".format(i))
                c.append(c_e)

            # objective function
            objective = 0
            for x_i, product in zip(x, products):
                objective += x_i * product.formation_energy_per_atom

            m.setObjective(objective, GRB.MINIMIZE)
            m.optimize()
            return [v.x for v in m.getVars()], m.objVal - reactant.formation_energy_per_atom


def find_greedy_old_first_choices(reactant: Compound, products: list[Compound], for_oxide: bool, firstk: int = None):
    dh_min = np.inf
    sol_min = None
    if firstk == None:
        first_choices = range(len(products))
    else:
        first_choices = range(min([len(products), firstk]))
    for i in first_choices:
        sol, dh = find_greedy_old(reactant, products, first_choice=i, for_oxide=for_oxide)
        # check elemental conservation
        assert check_solution(sol, products, reactant)
        if dh < dh_min:
            dh_min = dh
            sol_min = sol
    return sol_min, dh_min


def find_greedy_first_choices(
        reactant: Compound, products: list[Compound],
        diligent_greedy: bool, for_oxide: bool, firstk: int = None,
):
    dh_min = np.inf
    sol_min = None
    if firstk == None:
        first_choices = range(len(products))
    else:
        first_choices = range(min([len(products), firstk]))
    for i in first_choices:
        sol, dh = find_greedy(reactant, products, first_choice=i, diligent_greedy=diligent_greedy, for_oxide=for_oxide)
        # check elemental conservation
        assert check_solution(sol, products, reactant)
        if dh < dh_min:
            dh_min = dh
            sol_min = sol
    return sol_min, dh_min
