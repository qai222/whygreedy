import gzip
import json
import os
import pickle
import time
from pathlib import Path
from typing import Union

import numpy as np
from monty.json import MontyDecoder, MontyEncoder

file_type = Union[Path, str]


def json_dump(o, fn: file_type, compress=True) -> None:
    if compress:
        with gzip.open(fn, 'wt', encoding='UTF-8') as zipfile:
            json.dump(o, zipfile, cls=MontyEncoder)
    else:
        with open(fn, "w") as f:
            json.dump(o, f, cls=MontyEncoder)


def json_load(fn: file_type) -> dict:
    try:
        with gzip.open(fn, 'rt', encoding='UTF-8') as zipfile:
            d = json.load(zipfile, cls=MontyDecoder)
    except:
        with open(fn, "r") as f:
            d = json.load(f, cls=MontyDecoder)
    return d


def pkl_dump(o, fn: file_type) -> None:
    ts1 = time.perf_counter()
    with open(fn, "wb") as f:
        pickle.dump(o, f)
    ts2 = time.perf_counter()
    print("dumped {} in: {:.4f} s".format(os.path.basename(fn), ts2 - ts1))


def pkl_load(fn: file_type):
    ts1 = time.perf_counter()
    with open(fn, "rb") as f:
        d = pickle.load(f)
    ts2 = time.perf_counter()
    print("loaded {} in: {:.4f} s".format(os.path.basename(fn), ts2 - ts1))
    return d


def file_exists(fn: file_type):
    return os.path.isfile(fn) and os.path.getsize(fn) > 0


def set_small_to_zeros(a: list[float], eps=1e-5):
    a = np.array(a)
    a[np.abs(a) < eps] = 0
    return a
