"""
Functions to parse inputs & outputs
"""

import regex as re


class OutputParseError(Exception):
    """
    Thrown when we were unable to parse the specified output
    """
    pass


CP2K_GW_MATCH = re.compile(r'''
# anchor to indicate beginning of the HOMO/LUMO GW energy table
^[ \t]* MO [ \t]* E_SCF [ \t]* Sigc [ \t]* Sigc_fit [ \t]* Sigx-vxc [ \t]* Z [ \t]* E_GW \n

(
  ^
  [ \t]+ \d+
  [ \t]+ \(\ occ\ \)
  [ \t]+ (?P<HOMO_energy_SCF>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ ([\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ ([\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ ([\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ ([\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ (?P<HOMO_energy_GW>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  \n
)+  # the regex engine will only keep the last capture groups, effectively selecting the LUMO
(
  ^
  [ \t]+ \d+
  [ \t]+ \(\ vir\ \)
  [ \t]+ (?P<LUMO_energy_SCF>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ ([\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ ([\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ ([\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ ([\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]+ (?P<LUMO_energy_GW>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  \n
)  # do not match more since we want only the HOMO
''', re.MULTILINE | re.VERBOSE)


CP2K_MULLIKEN_MATCH = re.compile(r'''
# anchor to indicate beginning of the Mulliken Population Analysis
^[ \t]* Mulliken\ Population\ Analysis [ \t]* \n
 [ \t]* \n
 [ \t]* \#  [\w \t]+\n  # match the header
(
  ^
  [ \t]* (?P<atom>\d+)
  [ \t]* (?P<element>\w+)
  [ \t]* (?P<kind>\d+)
  [ \t]* (?P<population>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]* (?P<charge>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  \n
)+
(
  ^
  [ \t]* \#\ Total\ charge
  [ \t]* (?P<total_population>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
  [ \t]* (?P<total_charge>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
)
''', re.VERSION1 | re.MULTILINE | re.VERBOSE)


def parse_cp2k_output(fhandle):
    """
    Parse output from CP2K.

    Some basic info must be present in the input file.
    Others are purely optional
    """

    content = fhandle.read()

    def key_value_match(key, value_type=str):
        """
        Key-Value matching for specified strings on the content,
        with type conversion.
        """

        if value_type == float:
            value = r'[\+\-]?((\d*[\.]\d+)|(\d+[\.]?\d*))([Ee][\+\-]?\d+)?'
        elif value_type == int:
            value = r'[\+\-]?\d+'
        else:
            value = r'.+'

        match = re.search(r'^[ \t]*{key}[ \t]+(?P<value>{value})$'.format(
            key=re.escape(key),
            value=value  # match any number
            ), content, re.MULTILINE)

        if not match:
            return None

        return value_type(match.group('value'))

    data = {
        'version': key_value_match('CP2K| source code revision number:'),
        'mpiranks': key_value_match('GLOBAL| Total number of message passing processes', int),
        'threads': key_value_match('GLOBAL| Number of threads for this process', int),
        'username': key_value_match('**    ****   ******    PROGRAM STARTED BY'),
        'total_energy': key_value_match('ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):', float),
        'warnings_count': key_value_match('The number of warnings for this run is :', int),
        # the definition for the following is based on src/common/{termination,cp_error_handling}.F
        'warnings': re.findall(r'^[ \t]*\*{3} WARNING (in .+ :: .+) \*{3}', content, re.MULTILINE),
        }

    # only when doing calculations with kpoints
    nkpoints = key_value_match('BRILLOUIN| List of Kpoints [2 Pi/Bohr]', int)
    if nkpoints is not None:
        data['nkpoints'] = nkpoints

    match = CP2K_GW_MATCH.search(content)
    if match:
        data.update(match.groupdict())

    match = CP2K_MULLIKEN_MATCH.search(content)
    if match:
        captures = match.capturesdict()
        per_atom = []

        for idx in range(len(captures['atom'])):
            per_atom.append({
                'element': captures['element'][idx],
                'kind': int(captures['kind'][idx]),
                'population': float(captures['population'][idx]),
                'charge': float(captures['charge'][idx]),
                })

        data['mulliken_population_analysis'] = {
            'per-atom': per_atom,
            'total': {
                'population': float(match.group('total_population')),
                'charge': float(match.group('total_charge')),
                },
            }

    return data


def get_data_from_output(fhandle, code):
    """
    Parse code output data.
    Currently only implemented for QE & CP2K
    """

    if code == 'CP2K':
        return parse_cp2k_output(fhandle)

    elif code == 'espresso':
        data = {}

        for line in fhandle:
            if 'Program PWSCF v' in line:
                data['version'] = line.split()[2]
            if 'Number of MPI processes:' in line:
                data['mpiranks'] = int(line.split()[-1])
            if 'Threads/MPI process:' in line:
                data['threads'] = int(line.split()[-1])
            if 'number of k points=' in line:
                data['nkpoints'] = int(line.split()[4])
            if 'total cpu time spent up to now is' in line:
                data['runtime'] = float(line.split()[-2])
            if '!    total energy' in line:
                # extract the total energy and convert from Ry to eV
                data['total_energy'] = float(line.split()[-2])*13.605697827758654

        return data

    raise OutputParseError("Unknown code: %s" % code)
