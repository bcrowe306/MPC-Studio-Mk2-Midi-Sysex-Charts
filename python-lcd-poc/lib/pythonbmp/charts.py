"""
 Chart or graph specific numeric mod
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

from numbers import Number


def getdatalisttotal(dlist: list[Number]
                            ) -> Number:
    """Returns the total of a
        list of ints or floats

    Args:
        dlist: list of ints or floats

    Returns:
        float or int
    """
    return sum(d[0] for d in dlist)


def genpiechartdata(dlist: list):
    """Preprocess data to make
        it suitable for a pie chart

    Args:
        dlist: [[20, c['red']],
                [30, c['brightyellow']],
                ...]

    Returns:
        list and large value (if any)
    """
    sa = 0
    tot = getdatalisttotal(dlist)
    alist = []
    big = -1
    for d in dlist:
        p = d[0] / tot
        ea = sa + p * 360
        p *= 100
        alist.append([sa, ea, d[1], d[0], p])
        if p >= 50:
            big = dlist.index(d)
        sa = ea
    return alist, big
