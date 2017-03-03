"""
Functions to run checks on data.
"""

import numpy as np

def generate_checks_dict(data, code, test=""):
    """
    Generate a possibly nested dictionary with key-value pairs,
    where the key denotes a check and the value is always true
    if the check succeeded.

    Warning:
    When using NumPy, make sure to convert to native types before (even floats and bools),
    otherwise JSON serialization will fail.
    """


    if code == 'CP2K':
        checks = {
            'no_warnings': data['warnings_count'] == 0,
            'converged': not any('SCF run NOT converged' in w for w in data['warnings']),
            }

        if test.startswith('deltatest'):
            charges = [e['charge'] for e in data['mulliken_population_analysis']['per-atom']]

            checks['deltatest'] = {
                # For the mono-crystals in the deltatest we expect all charges to be equal,
                # check this by looking at the standard deviation
                'equal_charges_on_atoms': float(np.std(charges)) < 1e-3
                }

        return checks

    raise NotImplementedError
