from pymatgen.ext.matproj import MPRester

from whygreedy.utils import json_dump

if __name__ == '__main__':
    mat_api_key = 'Your MP API key'

    mpr = MPRester(mat_api_key)

    all_compounds = mpr.query({}, properties=["task_id", "pretty_formula", 'e_above_hull',
                                              'elements', 'volume', 'formation_energy_per_atom', 'band_gap',
                                              'nsites', 'unit_cell_formula', 'energy_per_atom'])

    json_dump(all_compounds, "mp.json.gz")
