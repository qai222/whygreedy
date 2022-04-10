import gzip
import json
from pathlib import Path
from typing import Union

file_type = Union[Path, str]


def json_load(jsonfilename: file_type):
    with gzip.open(jsonfilename, 'rt', encoding='UTF-8') as zipfile:
        my_object = json.load(zipfile)
    return my_object


def json_dump(data, jsonfilename: file_type):
    with gzip.open(jsonfilename, 'wt', encoding='UTF-8') as zipfile:
        json.dump(data, zipfile)
