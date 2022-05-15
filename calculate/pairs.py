from whygreedy import load_mp_oxidation_pairs, load_mp_decomposition_pairs, file_exists, pkl_load, pkl_dump

# a `pair` is a tuple of (reactant, product list)
# each pair correspond to a reaction, based on which the reaction enthalpy minimization is performed

mp_oxidation_pairs_pkl = "mp_oxidation_pairs.pkl"  # the pairs for oxidation reactions from stable, oxygen-free compounds
mp_decomp_pairs_pkl = "mp_decomp_pairs.pkl"  # the pairs for vacuum decomposition

if __name__ == '__main__':

    for pairs_pkl, load_pairs_function in zip(
            [mp_oxidation_pairs_pkl, mp_decomp_pairs_pkl],
            [load_mp_oxidation_pairs, load_mp_decomposition_pairs],
    ):
        if file_exists(pairs_pkl):
            print("Found file: {}".format(pairs_pkl))
            pairs = pkl_load(pairs_pkl)
        else:
            print("File not found: {}".format(pairs_pkl))
            print("loading with: {}".format(load_pairs_function.__name__))
            pairs = load_pairs_function()
            pkl_dump(pairs, pairs_pkl)
        print("# of pairs: {}".format(len(pairs)))
