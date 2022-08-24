"""
Parsing routines for different data structures.

All functions in this module work with Python3 native
datastructures.

Inspired from:
https://github.com/julie-forman-kay-lab/IDPConformerGenerator/blob/3aef6b085ec09eeebc5812639a5eb6832c0215cd/src/idpconfgen/libs/libparse.py
"""
import ast


def values_to_dict(values):
    """
    Generalization of converting parameters to dict.

    Adapted from:
    https://github.com/joaomcteixeira/taurenmd/blob/6bf4cf5f01df206e9663bd2552343fe397ae8b8f/src/taurenmd/libs/libcli.py#L94-L138

    Parameters
    ----------
    values : string
        List of values with the format "par1=1 par2='string' par3=[1,2,3]

    Returns
    -------
    param_dict : dictionary
        Converted string above to dictionary with `=` denoting linkage
        E.g. {'par1': 1, 'par2':'string', 'par3': [1,2,3]}
    """
    bool_value = {
        'true': True,
        'false': False,
        }

    param_dict = {}
    for kv in values:
        # print(param_dict, kv)
        try:
            k, v = kv.split('=')
        except ValueError:
            param_dict[kv] = True
        else:
            if ',' in v:
                vs = v.split(',')
                try:
                    param_dict[k] = tuple(ast.literal_eval(i) for i in vs)
                except (ValueError, TypeError, SyntaxError):
                    param_dict[k] = tuple(i for i in vs)
            else:
                try:
                    param_dict[k] = ast.literal_eval(v)
                except (ValueError, TypeError):  # is string or list
                    param_dict[k] = bool_value.get(v.lower(), v)
                except (SyntaxError):
                    param_dict[k] = v

    return param_dict
