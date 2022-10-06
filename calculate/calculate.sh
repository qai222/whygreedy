# lazy oxidation first 3
python calculate.py --records_pkl mp_oxidation_records_lazy_first3.pkl --pairs_pkl mp_oxidation_pairs.pkl --reaction_type oxidation --method lazy --firstk 3
# CRITICAL:root:time cost: 19.7088 s

# lp oxidation
python calculate.py --records_pkl mp_oxidation_records_lp.pkl --pairs_pkl mp_oxidation_pairs.pkl --reaction_type oxidation --method lp
# CRITICAL:root:time cost: 51.7470 s
