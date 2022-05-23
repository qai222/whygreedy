# lazy decomposition
python calculate.py --records_pkl mp_decomp_records_lazy.pkl --pairs_pkl mp_decomp_pairs.pkl --reaction_type decomposition --method lazy --parallel
# CRITICAL:root:time cost: 211545.9852 s

# lazy decomposition first 3
python calculate.py --records_pkl mp_decomp_records_lazy_first3.pkl --pairs_pkl mp_decomp_pairs.pkl --reaction_type decomposition --method lazy --firstk 3
# CRITICAL:root:time cost: 543.6980 s

# diligent decomposition
python calculate.py --records_pkl mp_decomp_records_diligent.pkl --pairs_pkl mp_decomp_pairs.pkl --reaction_type decomposition --method diligent --parallel
# CRITICAL:root:time cost: 252894.2312 s

# diligent decomposition first 3
python calculate.py --records_pkl mp_decomp_records_diligent_first3.pkl --pairs_pkl mp_decomp_pairs.pkl --reaction_type decomposition --method diligent --firstk 3
# CRITICAL:root:time cost: 560.9483 s

# lp decomposition
python calculate.py --records_pkl mp_decomp_records_lp.pkl --pairs_pkl mp_decomp_pairs.pkl --reaction_type decomposition --method lp
# CRITICAL:root:time cost: 545.8920 s

# pmg ehull
python calculate.py --records_pkl mp_decomp_records_pmgehull.pkl --pairs_pkl mp_decomp_pairs.pkl --reaction_type decomposition --method pmg


# lazy oxidation
python calculate.py --records_pkl mp_oxidation_records_lazy.pkl --pairs_pkl mp_oxidation_pairs.pkl --reaction_type oxidation --method lazy
# CRITICAL:root:time cost: 792.6788 s

# lazy oxidation first 3
python calculate.py --records_pkl mp_oxidation_records_lazy_first3.pkl --pairs_pkl mp_oxidation_pairs.pkl --reaction_type oxidation --method lazy --firstk 3
# CRITICAL:root:time cost: 19.7088 s

# diligent oxidation
python calculate.py --records_pkl mp_oxidation_records_diligent.pkl --pairs_pkl mp_oxidation_pairs.pkl --reaction_type oxidation --method diligent
# CRITICAL:root:time cost: 1365.1205 s

# diligent oxidation first 3
python calculate.py --records_pkl mp_oxidation_records_diligent_first3.pkl --pairs_pkl mp_oxidation_pairs.pkl --reaction_type oxidation --method diligent --firstk 3
# CRITICAL:root:time cost: 33.6773 s

# lp oxidation
python calculate.py --records_pkl mp_oxidation_records_lp.pkl --pairs_pkl mp_oxidation_pairs.pkl --reaction_type oxidation --method lp
# CRITICAL:root:time cost: 51.7470 s
