import numpy as np

from whygreedy.algo import check_solution, is_close_to_zero, find_greedy_first_choices, \
    find_greedy_old_first_choices, find_lp
from whygreedy.schema import Compound

"""
functions used in notebooks, define them here so they can be reused/used by pqdm in notebook env
"""


def calculate_diligent_vs_lazy_oxidation(pair: list[Compound, list[Compound]]):
    reactant, products = pair

    sol_old, dh_old = find_greedy_old_first_choices(reactant, products, for_oxide=True)
    sol_lazy, dh_lazy = find_greedy_first_choices(reactant, products, diligent_greedy=False, for_oxide=True)
    sol_diligent, dh_diligent = find_greedy_first_choices(reactant, products, diligent_greedy=True, for_oxide=True)

    # check elemental conservation
    assert check_solution(sol_old, products, reactant)
    assert check_solution(sol_lazy, products, reactant)
    assert check_solution(sol_diligent, products, reactant)

    # confirm we reproduce the old implementation
    assert is_close_to_zero(dh_lazy - dh_old) and np.allclose(sol_lazy, sol_old)

    record = {
        "reactant": reactant.mpid,
        "products": [p.mpid for p in products],
        "sol_old": sol_old,
        "sol_lazy": sol_lazy,
        "sol_diligent": sol_diligent,
        "dh_old": dh_old,
        "dh_lazy": dh_lazy,
        "dh_diligent": dh_diligent,
    }
    return record


def calculate_greedy_vs_lp_oxidation(pair: list[Compound, list[Compound]]):
    reactant, products = pair

    sol_old, dh_old = find_greedy_old_first_choices(reactant, products, for_oxide=True)
    sol_diligent, dh_diligent = find_greedy_first_choices(reactant, products, diligent_greedy=True, for_oxide=True)
    sol_lp, dh_lp = find_lp(reactant, products)

    # check elemental conservation
    assert check_solution(sol_old, products, reactant)
    assert check_solution(sol_diligent, products, reactant)
    assert check_solution(sol_lp, products, reactant)

    # greedy solution should be no better than the exact
    assert dh_lp <= min([dh_diligent, dh_old]) + 1e-7  # floating point error

    record = {
        "sol_old": sol_old,
        "sol_diligent": sol_diligent,
        "sol_lp": sol_lp,
        "dh_old": dh_old,
        "dh_diligent": dh_diligent,
        "dh_lp": dh_lp,
    }
    return record
