"""
Functions related to CP2K
"""

from io import BufferedIOBase


def dict2line_iter(nested, ilevel=0):
    """
    Iterator to convert a nested python dict to a CP2K input file.
    """

    indent = ' '*ilevel

    def _keyfunc(keyval):
        """Custom sorting function to get stable CP2K output"""
        key, val = keyval

        # ensure subsection generating entrys are sorted after keywords
        if isinstance(val, dict) or isinstance(val, list):
            # .. by prefixing them with the last possible character before sorting
            return chr(0x10ffff) + key.lower()

        return key.lower()

    for key, val in sorted(nested.items(), key=_keyfunc):
        if isinstance(val, dict):
            if '_' in val:
                yield "{}&{} {}".format(indent, key.upper(), val.pop('_'))
            else:
                yield "{}&{}".format(indent, key.upper())

            for line in dict2line_iter(val, ilevel + 3):
                yield line

            yield "{}&END {}".format(indent, key.upper())

        elif isinstance(val, list):
            for listitem in val:
                for line in dict2line_iter({key: listitem}, ilevel):
                    yield line

        elif isinstance(val, tuple):
            yield "{}{} {}".format(indent, key.upper(),
                                   ' '.join(str(v) for v in val))

        elif isinstance(val, bool):
            yield "{}{} {}".format(
                indent, key.upper(),
                '.TRUE.' if val else '.FALSE.')

        else:
            yield "{}{} {}".format(indent, key.upper(), val)


def dict2cp2k(data, output=None, parameters={}):  # pylint: disable=locally-disabled, dangerous-default-value
    """
    Convert and write a nested python dict to a CP2K input file.

    Some of this code is either taken from AiiDA or heavily
    inspired from there.

    Writes to a file if a handle or filename is given
    or returns the generated file as a string if not.

    The parameters are passed to a .format() call on the final input string.
    """

    if output:
        if isinstance(output, str):
            with open(output, 'w') as fhandle:
                fhandle.write("\n"
                              .join(dict2line_iter(data))
                              .format(**parameters))
        elif isinstance(output, BufferedIOBase):
            # bytes-like object, need to encode
            output.write("\n"
                         .join(dict2line_iter(data))
                         .format(**parameters)
                         .encode('utf-8'))
        else:
            output.write("\n".join(dict2line_iter(data)))
    else:
        return ("\n"
                .join(dict2line_iter(data))
                .format(**parameters))
