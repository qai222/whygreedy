# Linear programming for environmental stability of crystals

In a recent paper by Twyman et al. 
[(Chem. Mater. 2022, 34, 2545-2552)](https://pubs.acs.org/doi/abs/10.1021/acs.chemmater.1c02644), 
the environmental stability of crystals is 
investigated by determining the most favored oxidation reaction (lowest reaction enthalpy).
This is done using a greedy heuristic algorithm.

With two notebooks in this repository, we show that:
1. [greedy_versions.ipynb](greedy_versions.ipynb): The original implementation, 
kindly provided by Twyman et al. 
through [zenodo](https://zenodo.org/record/5110202#.YlJgpsjMJyg), 
is a "lazy" version of the greedy algorithm as described in their manuscript.
2. [whygreedy.ipynb](whygreedy.ipynb) The environmental stability problem can be solved exactly using linear programming, 
which gives solutions always better than or equivalent to the approximations by greedy algorithms.

## Usage
The content of both notebooks should be self-explanatory. 
To reproduce the notebooks, please follow these instructions:
1. Clone this repository and add the root folder to `$PYTHONPATH$`
2. Create a virtual environment and install dependencies: `pip install -r requirements.txt`
3. Extract the content of [data/data.7z](data/data.7z) to the [data/](data/) folder. 
   - latest `7z` can be downloaded [here](https://www.7-zip.org/download.html)
4. Start a Jupyter server and run the notebooks.

### precomputed results
The `*.pkl` files in the extracted content are precomputed results
to save computation time in notebooks. To reproduce them:
1. download materials project as `mp.json.gz` using `pymatgen` as described in [downloader.py](data/downloader.py)
2. extract reactions from `mp.json.gz` using [pairs.py](calculate/pairs.py)
3. calculate reaction enthalpies with [calculate.py](calculate/calculate.py), 
commands can be found in [calculate.sh](calculate/calculate.sh), and results will be saved as `*_records_*.pkl`.
4. [combine.py](calculate/combine.py) combines `*_records_*.pkl` to `mp_oxidation_records.pkl` and
`mp_decomposition_records.pkl` used in notebooks.
