"""    File utilities module
 -----------------------------------
| Copyright 2022 by Joel C. Alcarez |
| [joelalcarez1975@gmail.com]       |
|-----------------------------------|
|    We make absolutely no warranty |
| of any kind, expressed or implied |
|-----------------------------------|
|       The primary author and any  |
| any subsequent code contributors  |
| shall not be liable in any event  |
| for  incidental or consequential  |
| damages  in connection with,  or  |
| arising out from the use of this  |
| code in current form or with any  |
| modifications.                    |
 -----------------------------------
"""

from os.path import isfile
from typing import Callable
from .messages import sysmsg
from functools import wraps


def checklink(func: Callable):
    """Decorator to test if the first
        parameter in a function is
        a valid file

    Args:
        function

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        if isfile(args[0]):
            return(func(*args, **kwargs))
        else:
            print(sysmsg['filenotexist'])
    return(callf)


def checklinks(func: Callable):
    """Decorator to test if the two
        parameters in a function
        are valid files

    Args:
        function

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        if isfile(args[0]) and \
           isfile(args[1]):
            return(func(*args, **kwargs))
        else:
            print(sysmsg['filenotexist'])
    return(callf)
