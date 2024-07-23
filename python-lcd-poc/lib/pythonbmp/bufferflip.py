"""
Bit rotation and buffer flipping mod
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

from array import array
from .mathlib import rotatebits
from .colors import RGB2BGRbuf


def rotatebitsinbuf(
        buf: array) -> array:
    """Does a bit rotate to the bytes
        in an unsigned byte array

    Args:
        buf: unsigned byte array

    Returns:
        unsigned byte array
    """
    return array('B', [rotatebits(b)
                       for b in buf])


def flipnibbleinbuf(
        buf: array) -> array:
    """Flips a 4-bit image buffer

    Args:
        buf: unsigned byte array

    Returns:
        unsigned byte array
    """
    return array('B', [(b >> 4) + ((b % 16) << 4)
                    for b in buf])


def flip24bitbuf(buf: array) -> array:
    """Flips a 24-bit buffer

    Args:
        buf: unsigned byte array

    Returns:
        unsigned byte array
    """
    buf.reverse()
    RGB2BGRbuf(buf)
    return array('B',buf)


def flipbuf(buf: array, bits: int
            ) -> array:
    """Flips/rotates bits in buffer

    Args:
        buf : unsigned byte array
        bits: (1, 4, 8, 24) bits

    Returns:
        unsigned byte array
    """
    if bits == 24:
        buf = flip24bitbuf(buf)
    else:
        buf.reverse()
    if bits == 4:
        buf = flipnibbleinbuf(buf)
    if bits == 1:
        buf = rotatebitsinbuf(buf)
    return buf