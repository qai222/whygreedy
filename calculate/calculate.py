import argparse
import inspect
import logging
import os
import time

from whygreedy import pkl_load, pkl_dump
from whygreedy.algo import find_greedy_first_choices, find_greedy_old_first_choices, find_lp, find_pmgehull
from whygreedy.calculator import Calculator
from whygreedy.utils import file_type, file_exists


def get_kwargs():
    frame = inspect.currentframe().f_back
    keys, _, _, values = inspect.getargvalues(frame)
    kwargs = {}
    for key in keys:
        if key != 'self':
            kwargs[key] = values[key]
    return kwargs


def compute(
        method: str, records_pkl: file_type,
        pairs_pkl: file_type, firstk: int or None,
        reaction_type: str, parallel: bool
):
    if not file_exists(pairs_pkl):
        raise FileNotFoundError("pairs file not found!")
    logging.info("loading pairs file: {}".format(pairs_pkl))
    pairs = pkl_load(pairs_pkl)

    if file_exists(records_pkl):
        logging.info("found records file: {}".format(records_pkl))
        logging.info("will not compute anything, just sanity check")
        records = pkl_load(records_pkl)
        if len(records) != len(pairs):
            logging.critical("records has length: {}".format(len(records)))
            logging.critical("but pairs has length: {}".format(len(pairs)))
        if not all(isinstance(d, dict) for d in records):
            logging.critical("some records are not dictionary!")
        return records

    name = str(get_kwargs())

    cal_function_kwargs = {}
    if reaction_type == "oxidation":
        cal_function_kwargs["for_oxide"] = True
    elif reaction_type == "decomposition":
        cal_function_kwargs["for_oxide"] = False
    else:
        raise ValueError("reaction_type is: {}".format(reaction_type))

    if method == "lazy":
        cal_function = find_greedy_old_first_choices
        cal_function_kwargs["firstk"] = firstk
    elif method == "diligent":
        cal_function = find_greedy_first_choices
        cal_function_kwargs["diligent_greedy"] = True
        cal_function_kwargs["firstk"] = firstk
    elif method == "lp":
        cal_function = find_lp
        cal_function_kwargs = {}  # this does not take any kwarg
    elif method == "pmg":
        cal_function = find_pmgehull
        cal_function_kwargs = {}  # this does not take any kwarg
    else:
        raise ValueError("method is: {}".format(method))

    if method == "pmg" and reaction_type == "oxidation":
        raise ValueError("this cannot be done: method=={}, reaction_type=={}".format(method, reaction_type))

    calculator = Calculator(pairs=pairs, name=name, cal_function=cal_function, cal_function_kwargs=cal_function_kwargs)
    ts1 = time.perf_counter()
    if parallel:
        records = calculator.cal_parallel(n_jobs=os.cpu_count())
    else:
        records = calculator.cal_serial()
    pkl_dump(records, records_pkl)
    ts2 = time.perf_counter()
    logging.critical("time cost: {:.4f} s".format(ts2 - ts1))
    return records


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Compute reaction enthalpies.')
    parser.add_argument('--records_pkl', dest='records_pkl', metavar='records_pkl', type=str, nargs='?',
                        help='pkl filename for resulting records', default="mp_decomp_records_lp.pkl")
    parser.add_argument('--pairs_pkl', dest='pairs_pkl', metavar='pairs_pkl', type=str, nargs='?',
                        help='existing pkl file for `pairs` describing reactions', default="mp_decomp_pairs.pkl")
    parser.add_argument('--reaction_type', dest='reaction_type', metavar='reaction_type', type=str, nargs='?',
                        help='oxidation, decomposition', default='decomposition')
    parser.add_argument('--method', dest='method', type=str, nargs='?',
                        help='method for minimizing delta H', default='lp',
                        choices=['lazy', 'diligent', 'lp', 'pmg'])
    parser.add_argument('--firstk', dest='firstk', type=int, nargs='?',
                        help='how many different first choices to try in a greedy algorithm, default all choices',
                        default=None, )
    parser.add_argument('--parallel', action='store_true')

    args = parser.parse_args()
    logging.warning("arguments: {}".format(vars(args)))
    records = compute(
        method=args.method,
        records_pkl=args.records_pkl,
        pairs_pkl=args.pairs_pkl,
        firstk=args.firstk,
        reaction_type=args.reaction_type,
        parallel=args.parallel,
    )
