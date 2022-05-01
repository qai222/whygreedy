import numpy as np

from whygreedy.algo import find_greedy, find_greedy_old, check_solution, is_close_to_zero, find_lp
from whygreedy.schema import Compound

"""
functions used in notebooks, define them here so they can be used by pqdm in notebook env
"""


def calculate_diligent_vs_lazy_oxides(pair: list[Compound, list[Compound]]):
    original, oxides = pair

    dh_old_min = np.inf
    dh_lazy_min = np.inf
    dh_diligent_min = np.inf

    sol_old_min = None
    sol_lazy_min = None
    sol_diligent_min = None

    for i in range(len(oxides)):
        sol_old, dh_old = find_greedy_old(oxides, original, first_choice=i, for_oxides=True)
        sol_lazy, dh_lazy = find_greedy(oxides, original, first_choice=i, diligent_greedy=False, for_oxide=True)
        sol_diligent, dh_diligent = find_greedy(oxides, original, first_choice=i, diligent_greedy=True, for_oxide=True)

        # check elemental conservation
        assert check_solution(sol_old, oxides, original)
        assert check_solution(sol_lazy, oxides, original)
        assert check_solution(sol_diligent, oxides, original)

        # confirm we reproduce the old implementation
        assert is_close_to_zero(dh_lazy - dh_old) and np.allclose(sol_lazy, sol_old)

        if dh_old < dh_old_min:
            dh_old_min = dh_old
            sol_old_min = sol_old
        if dh_lazy < dh_lazy_min:
            dh_lazy_min = dh_lazy
            sol_lazy_min = sol_lazy
        if dh_diligent < dh_diligent_min:
            dh_diligent_min = dh_diligent
            sol_diligent_min = sol_diligent

    record = {
        "original": original,
        "oxides": oxides,
        "sol_old": sol_old_min,
        "sol_lazy": sol_lazy_min,
        "sol_diligent": sol_diligent_min,
        "dh_old": dh_old_min,
        "dh_lazy": dh_lazy_min,
        "dh_diligent": dh_diligent_min,
    }
    return record


def calculate_greedy_vs_exact_oxides(pair: list[Compound, list[Compound]], return_mpid=False):
    original, oxides = pair

    # greedy solutions
    dh_old_min = np.inf
    dh_diligent_min = np.inf
    sol_old_min = None
    sol_diligent_min = None
    for i in range(len(oxides)):
        sol_old, dh_old = find_greedy_old(oxides, original, first_choice=i)
        sol_diligent, dh_diligent = find_greedy(oxides, original, first_choice=i, diligent_greedy=True)
        if dh_old < dh_old_min:
            dh_old_min = dh_old
            sol_old_min = sol_old
        if dh_diligent < dh_diligent_min:
            dh_diligent_min = dh_diligent
            sol_diligent_min = sol_diligent

    sol_exact, dh_exact = find_lp(oxides, original)

    assert check_solution(sol_old_min, oxides, original)
    assert check_solution(sol_diligent_min, oxides, original)
    assert check_solution(sol_exact, oxides, original)

    # greedy solution should be no better than the exact
    assert dh_exact <= min([dh_diligent_min, dh_old_min]) + 1e-7  # floating point error

    # to reduce size
    if return_mpid:
        original = original.mpid
        oxides = [oxide.mpid for oxide in oxides]

    record = {
        "sol_old": sol_old_min,
        "sol_diligent": sol_diligent_min,
        "sol_exact": sol_exact,
        "dh_old": dh_old_min,
        "dh_diligent": dh_diligent_min,
        "dh_exact": dh_exact,
        "original": original,
        "oxides": oxides,
    }
    return record


def calculate_greedy_vs_exact_cp(pair: list[Compound, list[Compound]], return_mpid=True):
    original, oxides = pair
    # greedy solutions
    dh_old_min = np.inf
    dh_diligent_min = np.inf
    sol_old_min = None
    sol_diligent_min = None
    for i in range(len(oxides)):
        sol_old, dh_old = find_greedy_old(oxides, original, first_choice=i, for_oxides=False)
        sol_diligent, dh_diligent = find_greedy(oxides, original, first_choice=i, diligent_greedy=True, for_oxide=False)
        if dh_old < dh_old_min:
            dh_old_min = dh_old
            sol_old_min = sol_old
        if dh_diligent < dh_diligent_min:
            dh_diligent_min = dh_diligent
            sol_diligent_min = sol_diligent

    sol_exact, dh_exact = find_lp(oxides, original)

    assert check_solution(sol_old_min, oxides, original)
    assert check_solution(sol_diligent_min, oxides, original)
    assert check_solution(sol_exact, oxides, original)

    # greedy solution should be no better than the exact
    assert dh_exact <= min([dh_diligent_min, dh_old_min]) + 1e-7

    # to reduce size
    if return_mpid:
        original = original.mpid
        oxides = [oxide.mpid for oxide in oxides]

    record = {
        "sol_old": sol_old_min,
        "sol_diligent": sol_diligent_min,
        "sol_exact": sol_exact,
        "dh_old": dh_old_min,
        "dh_diligent": dh_diligent_min,
        "dh_exact": dh_exact,
        "original": original,
        "oxides": oxides,
    }
    return record
