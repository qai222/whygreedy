"""
Created on Thu Apr  4 18:07:03 2019

@author: NickT
"""


def find_comp(stable_oxides, compound_unit_cell, compound_formE, condition, n):
    """
    Finds complementary oxide or competing phases group and associated total
    heat of oxidation

    args:
        stable_oxides - list of dictionaries of stable oxides or competing
        phases with lower formation energy than original material
        compound_unit_cell - dict of elements in unit cell of original compound
        ompound_formE - formation energy of original compound
        condition - string dictating whether it is for comp oxide or comp
        competing phases
        n - forced first list selection (n = 0, 1, 2, 3, ...)

    output:
        tuple: (list of dicitionatries of predicted materials,
        combined formation energy of these materials (with appropriate ratios),
        combined formation energy of these materials (with appropriate ratios)
        minus original form E,
        number of compounds in complementary group,
        whether the grouping fuction finished early (boolean),
        n - forced first list selection (n = 0, 1, 2, 3, ...))

    notes:
        intersect_rank: used to find limiting element by finding ratio of
        normalised stochiometry between original material and oxide

    """
    result = []
    FinishEarly = False
    # what if positive formE

    orig_natoms = sum(compound_unit_cell.values())
    normalised_unit_cell = dict((a, b / orig_natoms) for a, b in compound_unit_cell.items())  # normalise stoichiometry

    for oxide in stable_oxides:
        oxide['el_weight'] = dict(
            (a, b / oxide['nsites']) for a, b in oxide['unit_cell_formula'].items())  # normalise stoichiometry

        if condition == 'Oxide':
            del oxide['el_weight']['O']

        oxide['ranker'] = dict(
            (a, b / normalised_unit_cell[a]) for a, b in oxide['el_weight'].items())  # find greedy ranking parameter

        # define how much material is used in total per unit cell
        oxide['ranking_no'] = sum(oxide['ranker'].values())

    # Order by energy per unit used up
    sort_oxides = sorted(stable_oxides, key=lambda oxide: (oxide['formation_energy_per_atom'] / oxide['ranking_no']))

    sort_oxides1 = sort_oxides[:]

    total_formE = 0
    counter = 0
    # if all atoms in unit cell not yet accounted for
    while sum(normalised_unit_cell.values()) != 0 and sort_oxides1 != []:
        if counter == 0:
            oxide = sort_oxides1[n]
        else:
            oxide = sort_oxides1[0]  # to allow forced initial choice

        # elements in oxide and orig. material
        intersection = list(set(oxide['elements']).intersection(normalised_unit_cell.keys()))
        # shouldnt we remove O from intersection???
        if len(intersection) == 0:
            print(compound_unit_cell)
            print(oxide['unit_cell_formula'])
            print(oxide['nsites'])
        intersect_rank = {}

        for element in intersection:
            # same as 1/ranker values
            intersect_rank[element] = normalised_unit_cell[element] / (
                    oxide['unit_cell_formula'][element] / oxide['nsites'])

        # find limiting element
        limiting_element = min(intersect_rank, key=intersect_rank.get)
        ratio = intersect_rank[limiting_element]  # (value)
        oxide['ratio'] = ratio  # For PBR calculation
        used_up_elements = []
        for element in intersection:

            normalised_unit_cell[element] = normalised_unit_cell[element] - (
                    ratio * oxide['unit_cell_formula'][element] / oxide['nsites'])

            # inequality because of != 0 problem
            if abs(normalised_unit_cell[element]) < 1e-7:  # was 1e-4
                used_up_elements.append(element)

        result.append(oxide)
        sort_oxides1.remove(oxide)
        total_formE += oxide['formation_energy_per_atom'] * ratio

        # remove oxides in list which arent useful (dont have new elements)
        sort_oxides1 = [oxide for oxide in sort_oxides1 if \
                        len(set(oxide['elements']).intersection(used_up_elements)) == 0]
        counter += 1

    # inequality because of != 0 problem
    if len(sort_oxides1) == 0 and abs(sum(normalised_unit_cell.values())) > 0.0001:
        FinishEarly = True

    return (result, total_formE, total_formE - compound_formE, len(result), FinishEarly, n)
