"""
combine results for decompositions/oxidations
"""

from whygreedy import pkl_load, pkl_dump

if __name__ == '__main__':
    # pairs to construct reactions
    mp_oxidation_pairs = pkl_load("mp_oxidation_pairs.pkl")
    mp_decomposition_pairs = pkl_load("mp_decomp_pairs.pkl")

    # precomputed data for oxidation
    mp_oxidation_records_lazy = pkl_load("mp_oxidation_records_lazy.pkl")
    mp_oxidation_records_lazy_first3 = pkl_load("mp_oxidation_records_lazy_first3.pkl")
    mp_oxidation_records_diligent = pkl_load("mp_oxidation_records_diligent.pkl")
    mp_oxidation_records_diligent_first3 = pkl_load("mp_oxidation_records_diligent_first3.pkl")
    mp_oxidation_records_lp = pkl_load("mp_oxidation_records_lp.pkl")

    # precomputed data for decomposition
    mp_decomp_records_lazy = pkl_load("mp_decomp_records_lazy.pkl")
    mp_decomp_records_lazy_first3 = pkl_load("mp_decomp_records_lazy_first3.pkl")
    mp_decomp_records_diligent = pkl_load("mp_decomp_records_diligent.pkl")
    mp_decomp_records_diligent_first3 = pkl_load("mp_decomp_records_diligent_first3.pkl")
    mp_decomp_records_lp = pkl_load("mp_decomp_records_lp.pkl")
    mp_decomp_records_pmgehull = pkl_load("mp_decomp_records_pmgehull.pkl")
    
    # combine results for decompositions
    mp_decomp_records = []
    for i in range(len(mp_decomposition_pairs)):
        record = {
            "sol_lazy_all": mp_decomp_records_lazy[i]["sol"],
            "sol_diligent_all": mp_decomp_records_diligent[i]["sol"],
            "sol_lazy_f3": mp_decomp_records_lazy_first3[i]["sol"],
            "sol_diligent_f3": mp_decomp_records_diligent_first3[i]["sol"],
            "sol_lp": mp_decomp_records_lp[i]["sol"],
            "sol_pmgehull": mp_decomp_records_pmgehull[i]["sol"],
            "dh_lazy_all": mp_decomp_records_lazy[i]["dh"],
            "dh_diligent_all": mp_decomp_records_diligent[i]["dh"],
            "dh_lazy_f3": mp_decomp_records_lazy_first3[i]["dh"],
            "dh_diligent_f3": mp_decomp_records_diligent_first3[i]["dh"],
            "dh_lp": mp_decomp_records_lp[i]["dh"],
            "dh_pmgehull": mp_decomp_records_pmgehull[i]["dh"],
        }
        mp_decomp_records.append(record)
    assert len(mp_decomp_records) == len(mp_decomposition_pairs)
    pkl_dump(mp_decomp_records, "mp_decomposition_records.pkl")
    
    # combine results for oxidations
    mp_oxidation_records = []
    for i in range(len(mp_oxidation_pairs)):
        record = {
            "sol_lazy_all": mp_oxidation_records_lazy[i]["sol"],
            "sol_diligent_all": mp_oxidation_records_diligent[i]["sol"],
            "sol_lazy_f3": mp_oxidation_records_lazy_first3[i]["sol"],
            "sol_diligent_f3": mp_oxidation_records_diligent_first3[i]["sol"],
            "sol_lp": mp_oxidation_records_lp[i]["sol"],
            "dh_lazy_all": mp_oxidation_records_lazy[i]["dh"],
            "dh_diligent_all": mp_oxidation_records_diligent[i]["dh"],
            "dh_lazy_f3": mp_oxidation_records_lazy_first3[i]["dh"],
            "dh_diligent_f3": mp_oxidation_records_diligent_first3[i]["dh"],
            "dh_lp": mp_oxidation_records_lp[i]["dh"],
        }
        mp_oxidation_records.append(record)
    assert len(mp_oxidation_records) == len(mp_oxidation_pairs)
    pkl_dump(mp_oxidation_records, "mp_oxidation_records.pkl")

    