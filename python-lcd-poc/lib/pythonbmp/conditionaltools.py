"""  Conditional functions module
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

def iif(boolcond: bool,
         trueval: any,
        falseval: any) -> any:
    """Returns trueval if
        boolcond is true
        else return falseval

    Args:
        boolcond: an expression that
                  evaluates as either
                  True or False
        trueval : value to return if
                  boolcond evaluates
                  to True
        falseval: value to return if
                  boolcond evaluates
                  to False

    Returns:
        a value depending on boolcond
    """
    return trueval if boolcond else falseval


def swapif(val1: any, val2: any,
       boolcond: bool):
    """Swaps val1 and val2 if
        boolcond is true

    Args:
        boolcond  : an expression that
                    evaluates as either
                    True or False
        val1, val2: values to swap if
                    boolcond is True

    Returns:
        values depending on boolcond
    """
    if boolcond:
        val1, val2 = val2, val1
    return val1, val2
