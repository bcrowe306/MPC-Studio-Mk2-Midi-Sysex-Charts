"""
 Buffer resize, pack and unpack mod
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
from .mathlib import(
    addvectinlist,
    intscalarmulvect,
    enumbits
    )

from .buffersplit import(
    altsplitbuf3way,
    altsplitbufnway
    )

from .colors import makeBGRbuf


def adjustbufsize(
      bufsize: int, bits: int) -> int:
    """Adjust buffer size to
        account for bit depth

    Args:
        bufsize: initial estimate
                 of buffer size
        bits   : bit depth of bitmap
                 (1, 4, 8, 24) bits

    Returns:
        An adjusted int value of the
        buffer size
    """
    if bits == 1:
        bufsize >>= 3
    elif bits == 24:
        bufsize *= 3
    elif bits == 4:
        bufsize >>= 1
    return bufsize


def adjustxbufsize(bmp: array,
    x1: int, x2: int) -> int:
    """Returns a 32 -bit padded
        int buffer size for a
        buffer with the bit depth
        stored in byte number 28
        of the bitmap bmp and an
        int length of x2 - x1 + 1

    Args:
        bmp   : unsigned byte array
                with bmp format
        x1, x2: int params for
                a x coord line
                in the bitmap

    Returns:
        int adjusted buffer size
    """
    return adjustbufsize(x2 - x1 + 1,
                            bmp[28])


def packbitlisttobuf(blist: list[int]
                       ) -> list[int]:
    """Packs literal list of ones and
        zeros to a list of bytes

    Args:
        blist: literal list
               of ones and zeros

    Returns:
        list
    """
    retval = []
    j = len(blist) + 1
    i = 1
    b = 0
    while i < j:
        m = i % 8
        b += blist[i - 1] << (7 - m)
        if m == 0 and i > 1:
            retval += [b]
            b = 0
        i += 1
    return retval


def resizebitpattenNtimesbigger(
           byteval: int, n: int):
    """Resize byte n times bigger
        bit wise

    Args:
        buf: unsigned byte
        n  : buffer multiplier

    Returns:
        list of ones and zeroes
    """
    retval = []
    for bit in enumbits(byteval):
        retval += [bit] * n
    return retval


def resize1bitbufNtimesbigger(
        buf: array, n: int):
    """Resize a 1-bit buffer
        n times bigger

    Args:
        buf: unsigned byte array
        n  : buffer multiplier

    Returns:
        unsigned byte array
    """
    retval = []
    for b in buf:
        retval += packbitlisttobuf(
          resizebitpattenNtimesbigger(
                        b, n))
    return array('B', retval)


def unpack4bitbuf(buf: list[int]
                   ) ->list[int]:
    """Unpacks a 4-bit buffer
        into a list

    Args:
        buf: unsigned byte array

    Returns:
        list
    """
    retval = []
    for b in buf:
        retval += [b >> 4, b & 0xf]
    return retval


def unpack4bitbufresizeNtimesbigger(
        buf: list[int], n: int
        ) -> list[int]:
    """Unpacks a 4-bit buffer into a
        list and repeats 4-bit units
        to resize the buffer int n
        times bigger

    Args:
        buf: unsigned byte array
        n  : unsigned int multiplier
             to resize buffer

    Returns:
        list
    """
    retval = []
    for b in buf:
        retval += [b >> 4] * n
        retval += [b & 0xf] * n
    return retval

def pack4bitbuf(
        unpackedbuf: list[int]
                ) -> list[int]:
    """Packs an unpacked 4-bit buffer
        into a list

    Args:
        buf: unsigned byte array
             or list

    Returns:
        list
    """
    retval = []
    j = len(unpackedbuf) - 1
    i = 0
    while i < j:
        retval += \
            [(unpackedbuf[i] << 4) + \
              unpackedbuf[i + 1]]
        i += 2
    return retval


def resize4bitbufNtimesbigger(
        buf: array, n: int
        ) -> array:
    """Resize a 4-bit buffer n times
        bigger

    Args:
        buf: An unsigned byte array
        n  : buffer multiplier

    Returns:
        unsigned byte array
    """
    return array('B', pack4bitbuf(
      unpack4bitbufresizeNtimesbigger(
          buf, n)))


def resize8bitbufNtimesbigger(
        buf: array, n: int
        ) -> array:
    """Resize a 8-bit buffer n times
        bigger

    Args:
        buf: unsigned byte array
        n  : buffer multiplier

    Returns:
        unsigned byte array
    """
    retval = []
    for b in buf:
        retval += [b] * n
    return array('B', retval)


def resizesmaller24bitbuf(buf: array
                          ) -> array:
    """Resize a 24-bit buffer
        n times smaller

    Args:
        buf: unsigned byte array
        n  : buffer multiplier

    Returns:
        unsigned byte array
    """
    n = len(buf)
    f = 1 / n**2
    b = altsplitbuf3way(
          addvectinlist(buf))
    return  makeBGRbuf(
            *[intscalarmulvect(
                 addvectinlist(
               altsplitbufnway(i, n)), f)
                           for i in b])


def resize24bitbufNtimesbigger(
            buf: array, n: int):
    """Resize a 24-bit buffer n times
        bigger

    Args:
        buf: unsigned byte array
        n  : buffer multiplier

    Returns:
        unsigned byte array
    """
    buf = altsplitbuf3way(buf)
    return makeBGRbuf(
           *[resize8bitbufNtimesbigger(
                   i, n) for i in buf])


def resizebufNtimesbigger(buf: array,
         n: int, bits: int) -> array:
    """Resize a buffer n times bigger
        given a particular bit depth n

    Args:
        buf : array to resize
        n   : resize factor
        bits: bit depth of
              color info
              (1, 4, 8, 24) bits

    Returns:
        list
    """
    f={24: resize24bitbufNtimesbigger,
        8: resize8bitbufNtimesbigger,
        4: resize4bitbufNtimesbigger,
        1: resize1bitbufNtimesbigger}[bits]
    return f(buf, n)
