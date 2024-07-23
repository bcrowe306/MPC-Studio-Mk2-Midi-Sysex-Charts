from typing import Callable
from inspect import signature

_d= '"'*3


def getfuncmetastr(f: Callable):
    """Gets Function Metadata in Python Code Format

    Args:
        f: function

    Returns:
        string function signature and docstring
    """
    return f'def {f.__name__}{signature(f)}:\n    {_d}{f.__doc__}{_d}\n'
