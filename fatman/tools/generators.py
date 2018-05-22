
import copy
from io import StringIO, BytesIO
from collections import OrderedDict

import numpy as np
from ase import io as ase_io

import click

from . import mergedicts
from .cp2k import dict2cp2k


def cp2k_array_merge_strategy(key, larr, rarr):
    """
    Merge strategy for merging arrays (currently only the kinds) of the CP2K config.
    Doesn't allow nuking elements in the list for now.
    """

    if key != "kind":
        raise RuntimeError("cp2k array merger only knows how to merge the kind section")

    merged = copy.deepcopy(larr)

    # check all elements in the source list
    for rentry in rarr:
        try:
            # figure out whether it should be overwritten
            idx, lentry = [(i, e) for i, e in enumerate(merged) if e['_'] == rentry['_']][0]
            # and use mergedicts to generate a merged dict and add it to the output list
            merged[idx] = dict(mergedicts(lentry, rentry))
        except IndexError:
            # if not found, simply add the entry
            merged.append(rentry)

    return merged


def generate_CP2K_inputs(settings, basis_sets, pseudos, struct, tagline, overrides=None):
    """
    Generate the inputs for CP2K based on the given data for CP2K

    Args:
        settings: CP2K-specific input settings dictionary
        basis_sets: List of tuples (type, id, element, family, basis)
        pseudos: List of tuples (id, element, family, ncore_el, pseudo)
        struct: A Python ASE atoms structure
        tagline: Comment line to add to generated files
        overrides: Input settings to be merged after autogenerating, just before generating the actual file

    Returns:
        a dictionary of (filename, bytebuf) objects
    """

    inputs = {}

    generated_input = {
        'global': {
            'project': "fatman.calc",
            },
        'force_eval': {
            'dft': {
                'basis_set_file_name': "./BASIS_SETS",
                'potential_file_name': "./POTENTIALS",
                'poisson': {
                    'periodic': None,
                    },
                },
            'subsys': {
                'cell': {  # filled out below
                    'a': None,  # filled out below
                    'b': None,
                    'c': None,
                    'periodic': None,
                    },
                'topology': {
                    'coord_file': "./struct.xyz",
                    'coord_file_format': 'XYZ',
                    },
                'kind': [],  # filled out below
                },
            },
        }

    periodic = "".join(axis*bool(enabled) for axis, enabled in zip("XYZ", struct.get_pbc()))
    if periodic == "":
        periodic = "NONE"

    generated_input['force_eval']['dft']['poisson']['periodic'] = periodic

    cell = {
        'a': (('[angstrom]',) + tuple(struct.get_cell()[0, :])),
        'b': (('[angstrom]',) + tuple(struct.get_cell()[1, :])),
        'c': (('[angstrom]',) + tuple(struct.get_cell()[2, :])),
        'periodic': periodic,
        }
    generated_input['force_eval']['subsys']['cell'] = cell

    if 'key_value_pairs' in struct.info:
        # it seems Python ASE is unable to handle nested dicts in
        # the Atoms.info attribute when writing XYZ, even though it
        # creates it in the first place
        # see https://gitlab.com/ase/ase/issues/60
        struct.info = dict(mergedicts(
            {k: v for k, v in struct.info.items() if k != 'key_value_pairs'},
            struct.info['key_value_pairs']))

    inputs['BASIS_SETS'] = BytesIO()
    inputs['BASIS_SETS'].write("# BASIS_SETS: {}\n".format(tagline).encode('utf-8'))

    # for the basis sets we have to be able to lookup the entry by element
    kind = {s: {'_': s, 'element': s, 'basis_set': [], 'potential': None} for s in struct.get_chemical_symbols()}

    for btype, _, element, family, _ in basis_sets:
        kind[element]['basis_set'].append(('ORB' if btype == 'default' else btype.upper(), family))

    # to write the basis set we can drop the type and therefore avoid
    # writing a basis twice in case we use the same basis for different types
    for basis_set in set(b[1:] for b in basis_sets):
        inputs['BASIS_SETS'].write(("# Basis Set ID {0}\n"
                                    "{1} {2}\n"
                                    "{3}\n")
                                   .format(*basis_set)
                                   .encode('utf-8'))

    inputs['BASIS_SETS'].seek(0)

    inputs['POTENTIALS'] = BytesIO()
    inputs['POTENTIALS'].write("# POTENTIALS: {}\n".format(tagline).encode('utf-8'))
    for pseudo in pseudos:
        kind[pseudo[1]]['potential'] = ("{2}-q{3}".format(*pseudo))
        # the format is checked when creating the Calculation
        inputs['POTENTIALS'].write(("# Pseudopotential ID {0}\n"
                                    "{1} {2}-q{3} {2}\n"
                                    "{4}\n")
                                   .format(*pseudo)
                                   .encode('utf-8'))
    inputs['POTENTIALS'].seek(0)

    # some older structures contain additional settings in the key_value_pairs dict,
    # add them to the generate settings
    if 'key_value_pairs' in struct.info:
        if 'multiplicity' in struct.info['key_value_pairs']:
            generated_input['force_eval']['dft']['multiplicity'] = struct.info['key_value_pairs']['multiplicity']

        if 'charge' in struct.info['key_value_pairs']:
            generated_input['force_eval']['dft']['charge'] = struct.info['key_value_pairs']['charge']

    # if we have any initial magnetic moment defined in the structure,
    # we need unrestricted KS + setting the magnetic moment and calculate the multiplicity
    if struct.get_initial_magnetic_moments().any():
        generated_input['force_eval']['dft']['uks'] = True
        # calculate the total magnetic moment, if not already set (for example by the structure above)
        if 'multiplicity' not in generated_input['force_eval']['dft']:
            # the total magnetization gives the difference in # of electrons between α and β spin, so 0.5*2*tmom=tmom
            generated_input['force_eval']['dft']['multiplicity'] = int(sum(struct.get_initial_magnetic_moments())) + 1

        # enumerate the different element+magmoms
        symmagmom = list(zip(struct.get_chemical_symbols(), struct.get_initial_magnetic_moments()))
        # and create a mapper (can not use set() here since it doesn't preserve the order)
        symmagmom2key = {(sym, magmom): '{}{}'.format(sym, num)
                         for num, (sym, magmom) in enumerate(OrderedDict.fromkeys(symmagmom), 1)}

        # replace the generated kind list by one containing the MAGNETIZATION
        mkind = {}
        for (sym, magmom), key in symmagmom2key.items():
            mkind[key] = kind[sym].copy()
            mkind[key]['_'] = key
            mkind[key]['magnetization'] = magmom
        kind = mkind

        # and add a new column to the atoms struct containing those labels
        struct.new_array('cp2k_labels', np.array([symmagmom2key[sm] for sm in symmagmom]))

    elif int(sum(struct.get_atomic_numbers() + struct.get_initial_charges())) % 2:
        generated_input['force_eval']['dft']['uks'] = True  # force enable LSD for systems with odd number of electrons

    else:
        # if no magnetic moments are required, simply copy the chemical symbols
        struct.new_array('cp2k_labels', np.array(struct.get_chemical_symbols()))

    # if charges are defined, set the charges keyword to be the total charge if nothing is set yet
    if struct.get_initial_charges().any() and 'charge' not in generated_input['force_eval']['dft']:
        generated_input['force_eval']['dft']['charge'] = int(sum(struct.get_initial_charges()))

    # in the CP2K Python dict struct the kinds are stored as list
    generated_input['force_eval']['subsys']['kind'] = list(kind.values())

    # merge the provided settings over the enerated input, giving the user the possibility
    # to override even auto-generated values
    combined_input = dict(mergedicts(generated_input, settings, cp2k_array_merge_strategy))

    # make some last adjustments which depend on a merged input structure

    try:
        # if scf is itself a dict we get a reference here
        scf = combined_input['force_eval']['dft']['scf']
        syms = struct.get_chemical_symbols()
        if 'smear' in scf.keys() and 'added_mos' not in scf.keys():
            n_mos = 0
            # when calculating the number of MOs on the other hand, we only want the default (for CP2K the "ORB")
            # type of basis sets since we don't want to count the AUX/RI/.. sets as well
            for _, _, element, _, basis in [b for b in basis_sets if b[0] == 'default']:
                basis_lines = basis.split('\n')
                nsets = int(basis_lines[0])
                lineno = 1  # start at the first set
                for _ in range(nsets):
                    # the number of MOs depends on the basis set
                    econfig_string = basis_lines[lineno]
                    econfig = [int(n) for n in econfig_string.split()]
                    # sum over (the number of m's per l quantum number times
                    # the number of functions per m) times the number of atoms of this kind:
                    n_mos += syms.count(element)*np.dot([2*l+1 for l in range(econfig[1], econfig[2]+1)], econfig[4:])
                    lineno += int(econfig[3]) + 1  # skip the block of coefficients and go to the next set

            scf['added_mos'] = max(int(0.3*n_mos), 1)  # at least one MO must be added
    except KeyError:
        pass

    # enforce same periodicity for CELL_REF (if specified) as for the cell itself
    try:
        combined_input['force_eval']['subsys']['cell']['cell_ref']['periodic'] = \
                combined_input['force_eval']['subsys']['cell']['periodic']
    except KeyError:
        pass

    # merge any additional settings on top if not None or empty
    if overrides:
        combined_input = dict(mergedicts(combined_input, overrides, cp2k_array_merge_strategy))

    inputs['calc.inp'] = BytesIO()
    inputs['calc.inp'].write("# calc.inp: {}\n".format(tagline).encode('utf-8'))
    dict2cp2k(combined_input, inputs['calc.inp'], parameters=struct.info)
    inputs['calc.inp'].seek(0)

    # write the structure after the input since we relabel the atoms for magmoms

    # we can't use TextIOWrapper here since this will close the underlying BytesIO
    # on destruction, resulting in a close BytesIO when leaving the scope
    stringbuf = StringIO()
    ase_io.write(stringbuf, struct, format='xyz', columns=['cp2k_labels', 'positions'])
    inputs['struct.xyz'] = BytesIO(stringbuf.getvalue().encode("utf-8"))
    inputs['struct.xyz'].seek(0)

    return inputs


def test():
    from . import json2atoms

    settings = {
        "force_eval": {
            "method": "Quickstep",
            "dft": {
                "multiplicity": 3,
                "scf": {
                    "eps_scf": 1e-08,
                    "smear": {
                        "method": "FERMI_DIRAC",
                        "_": True,
                        "electronic_temperature": "[K] 300"
                        },
                    "mixing": {
                        "method": "BROYDEN_MIXING",
                        "alpha": 0.4
                        }
                    },
                "xc": {"xc_functional": {"_": "PBE"}},
                "print": {"overlap_condition": {"1-norm": True, "_": "ON", "diagonalization": True}},
                "qs": {"method": "GPW", "extrapolation": "USE_GUESS"},
                "mgrid": {"cutoff": 1000, "rel_cutoff": 100},
                "kpoints": {
                    "full_grid": True,
                    "symmetry": False,
                    "parallel_group_size": -1,
                    "scheme": "MONKHORST-PACK {kpoints[0]} {kpoints[1]} {kpoints[2]}"
                    }
                },
            "subsys": {
                "kind": [
                    {"_": "O1", "magnetization": 2},
                    {"_": "O2", "magnetization": -2},
                    ]
                },
            },
        "global": {"run_type": "ENERGY", "print_level": "MEDIUM"}
        }
    basis_sets = [('default', 1, 'O', 'DZVP-GTH', ''' 1
2 0 2 7 2 2 1
    12.015954705512 -0.060190841200  0.065738617900  0.036543638800 -0.034210557400  0.014807054400
     5.108150287385 -0.129597923300  0.110885902200  0.120927648700 -0.120619770900  0.068186159300
     2.048398039874  0.118175889400 -0.053732406400  0.251093670300 -0.213719464600  0.290576499200
     0.832381575582  0.462964485000 -0.572670666200  0.352639910300 -0.473674858400  1.063344189500
     0.352316246455  0.450353782600  0.186760006700  0.294708645200  0.484848376400  0.307656114200
     0.142977330880  0.092715833600  0.387201458600  0.173039869300  0.717465919700  0.318346834400
     0.046760918300 -0.000255945800  0.003825849600  0.009726110600  0.032498979400 -0.005771736600
     ''')]
    pseudos = [(1, 'O', 'GTH-PBE', 6, ''' 2    4
     0.24455430    2   -16.66721480     2.48731132
    2
     0.22095592    1    18.33745811
     0.21133247    0
''')]
    struct = json2atoms('''
        {
            "numbers": [8, 8, 8, 8],
            "positions": [[4.27100416, 2.48306235e-17, 0.615915813], [1.80729466, 1.44465234e-16, 3.58341475], [2.13518916, 2.13841, 0.615915813], [3.94310966, 2.13841, 3.58341475]],
            "initial_magmoms": [1.5, 1.5, -1.5, -1.5],
            "cell": [[4.27163, 0.0, 0.0], [2.61879696e-16, 4.27682, 0.0], [1.80666882, 1.69295858e-16, 4.19933056]],
            "pbc": [true, true, true],
            "key_value_pairs": {"dataset": "deltatest", "identifier": "deltatest_O_1.00", "kpoints": [26, 24, 24]},
            "unique_id": "d54e984f10fd34539f0a8d1d4f89a7f4"
        }
        ''')
    tagline = "Generated by generator test"

    overrides = None

    inputs = generate_CP2K_inputs(settings, basis_sets, pseudos, struct, tagline, overrides=overrides)

    for filename, content in inputs.items():
        print("=== BEGIN: {}".format(filename))
        print(content.read().decode("utf-8"))
        print("=== END: {}\n".format(filename))


@click.command()
@click.argument("settings", type=click.File("r"))
@click.argument("struct", type=click.File("r"))
@click.argument("basisset", type=str)
@click.argument("pseudo", type=str)
def cli(settings, struct, basisset, pseudo):

    inputs = generate_CP2K_inputs(settings, basisset, pseudo, struct, tagline, overrides=overrides)

    for filename, content in inputs.items():
        print("=== BEGIN: {}".format(filename))
        print(content.read().decode("utf-8"))
        print("=== END: {}\n".format(filename))


if __name__ == '__main__':
    cli()
