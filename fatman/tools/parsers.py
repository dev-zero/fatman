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


CP2K_CONDITION_NUMBER_MATCH = re.compile(r'''
# anchor to indicate beginning the overlap matrix condition number section
^[ \t]* OVERLAP\ MATRIX\ CONDITION\ NUMBER\ AT\ GAMMA\ POINT [ \t]* \n
 [ \t]* 1-Norm\ Condition\ Number\ \(Estimate\) [ \t]* \n
 [ \t]* CN\ :\ \|A\|\*\|A\^-1\|:
   [ \t]* (?P<norm1_estimate_A>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* \*
   [ \t]* (?P<norm1_estimate_Ainv>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* =
   [ \t]* (?P<norm1_estimate>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* Log\(1-CN\):
   [ \t]* (?P<norm1_estimate_log>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* \n

 [ \t]* 1-Norm\ and\ 2-Norm\ Condition\ Numbers\ using\ Diagonalization [ \t]* \n

 [ \t]* CN\ :\ \|A\|\*\|A\^-1\|:
   [ \t]* (?P<norm1_diag_A>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* \*
   [ \t]* (?P<norm1_diag_Ainv>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* =
   [ \t]* (?P<norm1_diag>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* Log\(1-CN\):
   [ \t]* (?P<norm1_diag_log>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* \n

 [ \t]* CN\ :\ max/min\ ev:
   [ \t]* (?P<norm2_diag_max_ev>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* /
   [ \t]* (?P<norm2_diag_min_ev>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* =
   [ \t]* (?P<norm2_diag>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* Log\(2-CN\):
   [ \t]* (?P<norm2_diag_log>[\+\-]?(\d*[\.]\d+|\d+[\.]?\d*)([Ee][\+\-]?\d+)?)
   [ \t]* \n
''', re.MULTILINE | re.VERBOSE)

# the definition for the following is based on src/common/{termination,cp_error_handling}.F
CP2K_WARNINGS_MATCH = re.compile(r'''
^[ \t] \*{3} [ \t]+ WARNING [ \t]+ (in [ \t]+ .+ [ \t]+ :: [ \t]+ .+) [ \t]+ \*{3} [ \t]* \n
(?:
  ^[ \t] \*{3} [ \t]+ (.+?) [ \t]+ \*{3} \n
)*
''', re.MULTILINE | re.VERBOSE)


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
            value = r'[\+\-]?\d+' # for ints we don't allow scientific notation since python doesn't accept it
        else:
            value = r'.+' # this will match everthing except a newline

        # value now contains the regex corresponding to a type

        # complete the parametrized regex and search for the given key
        # and specified value format, capture the value
        match = re.search(r'^[ \t]*{key}[ \t]+(?P<value>{value})$'.format(
            key=re.escape(key), # make sure we match key as it is (no regex interpretation)
            value=value
            ), content, re.MULTILINE)

        if not match:
            return None

        # convert the captured value string to the requested python type
        return value_type(match.group('value'))

    data = {
        'version': key_value_match('CP2K| source code revision number:'),
        'mpiranks': key_value_match('GLOBAL| Total number of message passing processes', int),
        'threads': key_value_match('GLOBAL| Number of threads for this process', int),
        'username': key_value_match('**    ****   ******    PROGRAM STARTED BY'),
        'total_energy': key_value_match('ENERGY| Total FORCE_EVAL ( QS ) energy (a.u.):', float),
        'warnings_count': key_value_match('The number of warnings for this run is :', int),
        'warnings': [" ".join(warning_match.groups()) for warning_match in CP2K_WARNINGS_MATCH.finditer(content)],
        }

    # only when doing calculations with kpoints
    nkpoints = key_value_match('BRILLOUIN| List of Kpoints [2 Pi/Bohr]', int)
    if nkpoints is not None:
        data['nkpoints'] = nkpoints

    match = CP2K_GW_MATCH.search(content)
    if match:
        data.update(match.groupdict())

    match = CP2K_MULLIKEN_MATCH.search(content)
    # for this one we needed the extended regex library https://pypi.python.org/pypi/regex
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

    match = CP2K_CONDITION_NUMBER_MATCH.search(content)
    if match:
        captures = match.groupdict()
        data['overlap_matrix_condition_number'] = {
            '1-norm (estimate)': {
              '|A|': float(captures['norm1_estimate_A']),
              '|A^-1|': float(captures['norm1_estimate_Ainv']),
              'CN': float(captures['norm1_estimate']),
              'Log(CN)': float(captures['norm1_estimate_log']),
              },

            '1-norm (using diagonalization)': {
              '|A|': float(captures['norm1_diag_A']),
              '|A^-1|': float(captures['norm1_diag_Ainv']),
              'CN': float(captures['norm1_diag']),
              'Log(CN)': float(captures['norm1_diag_log']),
              },

            '2-norm (using diagonalization)': {
              'max EV': float(captures['norm2_diag_max_ev']),
              'min EV': float(captures['norm2_diag_min_ev']),
              'CN': float(captures['norm2_diag']),
              'Log(CN)': float(captures['norm2_diag_log']),
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
