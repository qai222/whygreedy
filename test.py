import random

import numpy as np
import pytest

from whygreedy import pkl_load, json_load, find_lp, find_greedy_first_choices


class TestChemmat:

    @pytest.fixture
    def oxidation_pairs(self):
        return pkl_load("data/mp_oxidation_pairs.pkl")

    @pytest.fixture
    def oxidation_records(self):
        return pkl_load("data/mp_oxidation_records.pkl")

    @pytest.fixture
    def mp_data(self):
        return json_load("data/mp.json.gz")

    def test_dataload(self, oxidation_pairs, oxidation_records, mp_data):
        assert len(mp_data) == 126335
        assert len(oxidation_pairs) == 39634
        assert len(oxidation_records) == len(oxidation_pairs)

    def test_lp(self, oxidation_pairs, oxidation_records):
        random.seed(42)
        for i in random.sample(range(len(oxidation_pairs)), 100):
            sol, dh = find_lp(*oxidation_pairs[i])
            assert np.allclose(sol, oxidation_records[i]['sol_lp'])
            assert np.allclose(dh, oxidation_records[i]['dh_lp'])

    def test_lazyf3(self, oxidation_pairs, oxidation_records):
        random.seed(42)
        for i in random.sample(range(len(oxidation_pairs)), 100):
            sol, dh = find_greedy_first_choices(*oxidation_pairs[i], for_oxide=True, diligent_greedy=False, firstk=3)
            assert np.allclose(sol, oxidation_records[i]['sol_lazy_f3'])
            assert np.allclose(dh, oxidation_records[i]['dh_lazy_f3'])
