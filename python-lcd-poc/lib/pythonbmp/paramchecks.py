"""
Function decorators for parameter check
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

from functools import wraps
from .messages import sysmsg
from .bmpconstants import(
    bmpcolorbits, bmpx, bmpy)

from .inttools import readint

from .primitives2D import(
    entirecircleisinboundary)


f = lambda bmp, x, y: (x < readint(bmpx,4,bmp) and \
                       y < readint(bmpy,4,bmp)) and \
                      (x > -1 and y > -1)


def intcircleparam(func):
    """Decorator to test if the
        2nd, 3rd, 4th parameters
        in a function that renders
        circle are ints

    Args:
        function(bmp:array,x,y,r....)

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        if (type(args[1]) == int and \
            type(args[2]) == int) and \
            type(args[3]) == int:
            return(func(*args, **kwargs))
        else:
            print(sysmsg['inttypereq'])
    return(callf)


def intcircleparam24bitonly(func):
    """Decorator to test if 2nd, 3rd,
        4th parameters in a function
        that renders circle are ints
        and restrict the use of this
        function to only 24-bit or
        RGB bitmaps (1st parameter)

    Args:
        function(bmp:array,x,y,r....)

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        if args[0][bmpcolorbits] != 24:
            print(sysmsg['not24bit'])
        elif (type(args[1]) == int and \
                type(args[2]) == int) and \
                type(args[3]) == int:
            return(func(*args, **kwargs))
        else:
            print(sysmsg['inttypereq'])
    return(callf)


def func24bitonly(func):
    """Decorator to restrict the
        use of this function to only
        24-bit or RGB bitmaps
        (1st parameter)

    Args:
        function(bmp:array,...)

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        if args[0][bmpcolorbits] != 24:
            print(sysmsg['not24bit'])
        else:
            return(func(*args, **kwargs))
    return(callf)


def func24bitonlyandentirerectinboundary(func):
    """Decorator to restrict the
        use of this function to only
        24-bit or RGB bitmaps
        (1st parameter) and ensure that
        the 2nd, 3rd, 4th and 5th
        parameters are ints whose
        values when interpreted as
        x and y coordinates lay
        within the RGB bitmap.

    Args:
        function

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        bmp = args[0]
        x1 = args[1]
        y1 = args[2]
        x2 = args[3]
        y2 = args[4]

        if bmp[bmpcolorbits] != 24:
            print(sysmsg['not24bit'])
        elif (type(x1) == int and type(x2) == int) and \
               (type(y1) == int and type(y2) == int):

            if f(bmp, x1, y1) and f(bmp, x2, y2):
                return(func(*args, **kwargs))
            else:
                print(sysmsg['regionoutofbounds'])
        else:
            print(sysmsg['inttypereq'])
    return(callf)


def func24bitonlyandentirecircleinboundary(func):
    """Decorator to restrict the
        use of this function to only
        24-bit bitmaps (1st parameter)
        and ensure that the 2nd, 3rd,
        4th parameters are ints whose
        values when interpreted as
        x, y and radius of a circle
        lay within the RGB bitmap.

    Args:
        function(bmp:array,x:int,y:int,r:int...)

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        bmp = args[0]
        x = args[1]
        y = args[2]
        r = args[3]
        if bmp[bmpcolorbits] != 24:
            print(sysmsg['not24bit'])
        elif (type(x) == int and type(y) == int) \
                               and type(r) == int:
            if entirecircleisinboundary(
                x, y, -1, readint(bmpx,4,bmp),
                      -1, readint(bmpy,4,bmp), r):
                return(func(*args, **kwargs))
            else:
                print(sysmsg['regionoutofbounds'])
        else:
            print(sysmsg['inttypereq'])
    return(callf)


def func8and24bitonlyandentirecircleinboundary(func):
    """Decorator to restrict the
        use of this function to only
        24-bit or 8-bit bitmaps
        (1st parameter) and ensure
        that the 2nd, 3rd, 4th
        parameters are ints whose
        values when interpreted as
        x, y and radius of a circle
        lay within the RGB bitmap.

    Args:
        function(bmp:array,x:int,y:int,r:int...)

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        bmp = args[0]
        x = args[1]
        y = args[2]
        r = args[3]
        if bmp[bmpcolorbits] not in [24, 8]:
            print(sysmsg['not24or8bit'])
        elif (type(x) == int and type(y) == int) \
                               and type(r) == int:
            if entirecircleisinboundary(
                x, y, -1, readint(bmpx, 4, bmp),
                      -1, readint(bmpy, 4, bmp), r):
                return(func(*args, **kwargs))
            else:
                print(sysmsg['regionoutofbounds'])
        else:
            print(sysmsg['inttypereq'])
    return(callf)


def func8and24bitonly(func):
    """Decorator to restrict the
        use of this function to
        only 24-bit or 8-bit bitmaps
        (1st parameter)

    Args:
        function(bmp:array,...)

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        if args[0][bmpcolorbits] not in [24, 8]:
            print(sysmsg['not24or8bit'])
        else:
            return(func(*args, **kwargs))
    return(callf)


def func4and8bitonly(func):
    """Decorator to restrict the
        use of this function to
        only 4-bit or 8-bit bitmaps
        (1st parameter)

    Args:
        function(bmp:array,...)

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        if args[0][bmpcolorbits] not in [4, 8]:
            print(sysmsg['not4or8bit'])
        else:
            return(func(*args, **kwargs))
    return(callf)


def func8and24bitonlyandentirerectinboundary(func):
    """Decorator to restrict the
        use of this functiom to only
        24 bit or 8 bit bitmaps
        (1st parameter) and ensure
        that the 2nd, 3rd, 4th and
        5th parameters are ints whose
        values when interpreted as
        x and y coordinates of a
        rectangle lay within
        the RGB bitmap.

    Args:
        function

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        bmp = args[0]
        x1 = args[1]
        y1 = args[2]
        x2 = args[3]
        y2 = args[4]
        if bmp[bmpcolorbits] not in [24, 8]:
            print(sysmsg['not24or8bit'])
        elif (type(x1) == int and type(x2) == int) and \
               (type(y1) == int and type(y2) == int):
            if f(bmp, x1, y1) and f(bmp, x2, y2):
                return(func(*args, **kwargs))
            else:
                print(sysmsg['regionoutofbounds'])
        else:
            print(sysmsg['inttypereq'])
    return(callf)


def entirerectinboundary(func):
    """Decorator to ensure that the
        2nd, 3rd, 4th and 5th
        parameters are ints whose
        values when interpreted as
        x and y coordinates of a
        rectangle lay within the
        bitmap.

    Args:
        function

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        bmp = args[0]
        x1 = args[1]
        y1 = args[2]
        x2 = args[3]
        y2 = args[4]
        if (type(x1) == int and type(x2) == int) and \
           (type(y1) == int and type(y2) == int):
            if not (f(bmp, x1, y1) and f(bmp, x2, y2)):
                print(sysmsg['regionoutofbounds'])
            else:
                return(func(*args, **kwargs))
        else:
            print(sysmsg['inttypereq'])
    return(callf)


def entirecircleinboundary(func):
    """Decorator to ensure that the 2nd,
        3rd, 4th parameters are ints
        whose values when interpreted
        as the centerpoint x, y
        and radius r of a circle
        lay within the RGB bitmap.

    Args:
        function(bmp:array,x:int,y:int,r:int...)

    Returns:
        caller function
    """
    @wraps(func)
    def callf(*args, **kwargs):
        bmp = args[0]
        x = args[1]
        y = args[2]
        r = args[3]
        if (type(x) == int and type(y) == int) \
                           and type(r) == int:
            if entirecircleisinboundary(
                    x, y, -1, readint(bmpx, 4, bmp),
                          -1, readint(bmpx, 4, bmp), r):
                return(func(*args, **kwargs))
            else:
                print(sysmsg['regionoutofbounds'])
        else:
            print(sysmsg['inttypereq'])
    return(callf)
