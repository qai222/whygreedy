from typing import Tuple, Callable

from pqdm.processes import pqdm
from tqdm import tqdm

from whygreedy.schema import Compound


class Calculator:
    def __init__(self, pairs: list[Tuple[Compound, list[Compound]]], name: str,
                 cal_function: Callable, cal_function_kwargs: dict):
        self.pairs = pairs
        self.name = name
        self.cal_function = cal_function
        self.cal_function_kwargs = cal_function_kwargs

    def cal_serial(self, k: int = None):
        if k is None:
            pairs = self.pairs
        else:
            pairs = self.pairs[:k]
        records = []
        for p in tqdm(pairs):
            record = self.cal_one(p)
            records.append(record)
        return records

    def cal_one(self, p):
        reactant, products = p
        sol, dh = self.cal_function(reactant=reactant, products=products, **self.cal_function_kwargs)
        return dict(
            sol=sol, dh=dh,
            reactant=reactant.mpid,
            products=[prod.mpid for prod in products]
        )

    def cal_parallel(self, n_jobs=8):
        pairs = self.pairs
        records = pqdm(pairs, self.cal_one, n_jobs=n_jobs)
        return records
