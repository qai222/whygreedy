import numpy as np

from whygreedy import Compound, normalize_stoi, find_greedy, find_greedy_old, check_solution, \
    is_close_to_zero, gen_random_data, find_lp


def test_abco4():
    """
    use the test data from `Analyser.py` of the original code
    """
    ABCO4 = {'elements': ['A', 'B', 'C', 'O'], 'formation_energy_per_atom':
        -750, 'nsites': 7, 'unit_cell_formula': {'A': 1, 'B': 1, 'C': 1, 'O': 4}}

    AO = {'elements': ['A', 'O'], 'formation_energy_per_atom': -100,
          'nsites': 8, 'unit_cell_formula': {'A': 4, 'O': 4}}

    BO2 = {'elements': ['B', 'O'], 'formation_energy_per_atom': -100,
           'nsites': 6, 'unit_cell_formula': {'B': 2, 'O': 4}}

    C2O = {'elements': ['C', 'O'], 'formation_energy_per_atom': -300,
           'nsites': 24, 'unit_cell_formula': {'C': 16, 'O': 8}}

    A2BO6 = {'elements': ['A', 'B', 'O'], 'formation_energy_per_atom': -380,
             'nsites': 9, 'unit_cell_formula': {'A': 2, 'B': 1, 'O': 6}}

    A2CO4 = {'elements': ['A', 'C', 'O'], 'formation_energy_per_atom': -620,
             'nsites': 63, 'unit_cell_formula': {'A': 18, 'C': 9, 'O': 36}}

    original = {'A': 4, 'B': 8, 'C': 100}
    oxides = [ABCO4, AO, BO2, C2O, A2BO6, A2CO4]
    original = Compound(normalized_formula=normalize_stoi(original), formation_energy_per_atom=-400)
    oxides = [Compound(normalized_formula=normalize_stoi(oxi['unit_cell_formula']),
                       formation_energy_per_atom=oxi['formation_energy_per_atom']) for oxi in oxides]

    for first_choice in range(len(oxides)):
        sol_new, dh_new = find_greedy(oxides, original, first_choice=first_choice, real_greedy=False)
        check_solution(sol_new, oxides, original)
        sol_old, dh_old = find_greedy_old(oxides, original, first_choice=first_choice)
        check_solution(sol_old, oxides, original)
        assert np.allclose(sol_new, sol_old) and is_close_to_zero(dh_old - dh_new), \
            "New implementation is different from the old one!!"


def test_random(n=1000):
    print("\t".join(["seed[0-{}]".format(n), "old", "new", "real", "exact"]))
    n_way_worse_than_exact = 0
    for seed in range(n):
        original, oxides = gen_random_data(["A", "B", "C", "D"], 5, seed)

        sol_old, dh_old = find_greedy_old(oxides, original, first_choice=0)
        sol_new, dh_new = find_greedy(oxides, original, first_choice=0, real_greedy=False)
        sol_real, dh_real = find_greedy(oxides, original, first_choice=0, real_greedy=True)
        sol_exact, dh_exact = find_lp(oxides, original)

        assert check_solution(sol_old, oxides, original)
        assert check_solution(sol_new, oxides, original)
        assert check_solution(sol_real, oxides, original)
        assert check_solution(sol_exact, oxides, original)

        # the new implementation should reproduce old results
        assert np.allclose(sol_old, sol_new) and is_close_to_zero(dh_old - dh_new)
        # printout inconsistent results
        if not is_close_to_zero(dh_old - dh_real):
            print("{}\t{:.3f}\t{:.3f}\t{:.3f}\t{:.3f}".format(seed, dh_old, dh_new, dh_real, dh_exact))

        # greedy solution should be no better than the exact
        assert dh_exact <= min([dh_real, dh_old, dh_new]) + 1e-5
        if dh_exact < 0 and dh_exact < min([dh_real, dh_old, dh_new]) * 1.05:
            n_way_worse_than_exact += 1
    print("# of old 5% worse than exact: {}".format(n_way_worse_than_exact))


if __name__ == '__main__':
    test_abco4()
    test_random(1000)
    """
    seed[0-1000]	old	new	real	exact
    20	-5808.957	-5808.957	-5874.024	-5874.024
    447	-2069.512	-2069.512	-2080.216	-2722.686
    553	-3611.081	-3611.081	-3599.704	-3611.081
    555	-3155.043	-3155.043	-3324.340	-3324.340
    914	-2097.586	-2097.586	-2212.100	-2212.100
    915	-2156.731	-2156.731	-2090.781	-2156.731
    # of old 5% worse than exact: 116
    """
