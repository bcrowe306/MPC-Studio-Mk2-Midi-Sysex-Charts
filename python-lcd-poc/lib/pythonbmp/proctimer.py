"""        Timer module
 -----------------------------------
| Copyright 2022 by Joel C. Alcarez |
| [joelalcarez1975@gmail.com]       |
|-----------------------------------|
|    We make absolutely no warranty |
| of any kind, expressed or implied |
|-----------------------------------|
|   Contact primary author          |
|   if you plan to use this         |
|   in a commercial product at      |
|   joelalcarez1975@gmail.com       |
 -----------------------------------
"""

from time import process_time_ns
from functools import wraps


def elaspedtimeinseconds(inittime):
    """Get elasped time in seconds

    Args:
        inittime: int start time in
                  nanoseconds

    Returns:
        int elapsed time in seconds
    """
    return (process_time_ns() - inittime) / \
            1000000000


def hhmmsselaspedtime(inittime: int
                           ) -> str:
    """Get elasped time

    Args:
        inittime: int start time in
                  nanoseconds

    Returns:
      formatted string for elapsed time
    """
    secs, ns = divmod(
                (process_time_ns() - inittime),
                1000000000)
    mins, secs = divmod(secs, 60)
    hrs, mins = divmod(mins, 60)
    return (
        ((f'{str(hrs).zfill(2)}:' + str(mins).zfill(2)) + ':')
        + str(secs).zfill(2)
        + '.'
    ) + str(ns)


def functimer(func):
    """Function timer Decorator

    Args:
        function

    Returns:
        caller function

    """
    @wraps(func)
    def callf(*args, **kwargs):
        print((f'Applying {func.__name__}' + ' please wait...'))
        inittime = process_time_ns()
        r = func(*args, **kwargs)
        print("Done in: " + \
               hhmmsselaspedtime(inittime))
        return(r)
    return(callf)
