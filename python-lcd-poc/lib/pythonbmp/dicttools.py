""" Dictionaries and lists module
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


def dict2descorderlist(
        d: dict) -> list:
    """Creates a sorted list in
        decending order from
        a dictionary with counts

    Args:
        dict: histogram or
              frequency count

    Returns:
        list
    """
    l = [[v, k] for k, v in d.items()]
    l.sort()
    l.reverse()
    return l