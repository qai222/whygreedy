"""
combine results for decompositions/oxidations
"""

from whygreedy import pkl_load, pkl_dump

if __name__ == '__main__':
    # pairs to construct reactions
    mp_oxidation_pairs = pkl_load("mp_oxidation_pairs.pkl")

    # precomputed data for oxidation
    mp_oxidation_records_lazy_first3 = pkl_load("mp_oxidation_records_lazy_first3.pkl")
    mp_oxidation_records_lp = pkl_load("mp_oxidation_records_lp.pkl")

    # combine results for oxidations
    mp_oxidation_records = []
    for i in range(len(mp_oxidation_pairs)):
        record = {
            "sol_lazy_f3": mp_oxidation_records_lazy_first3[i]["sol"],
            "sol_lp": mp_oxidation_records_lp[i]["sol"],
            "dh_lazy_f3": mp_oxidation_records_lazy_first3[i]["dh"],
            "dh_lp": mp_oxidation_records_lp[i]["dh"],
        }
        mp_oxidation_records.append(record)
    assert len(mp_oxidation_records) == len(mp_oxidation_pairs)
    pkl_dump(mp_oxidation_records, "mp_oxidation_records.pkl")
