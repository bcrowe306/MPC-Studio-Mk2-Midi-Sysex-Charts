"""Integer buffer functions module
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


def writeint(offset: int, cnt: int,
            arr: array, value: int):
    """Writes an integer value to an
        unsigned byte array

    Args:
        offset: uint starting offset in
                buffer or array
                to write to
        cnt   : uint length of int data
                to write
        arr   : unsigned byte array
                to write int data in
        value : value of uint data to
                write in buffer or
                array

    Returns:
        byref unsigned byte array
    """
    j = cnt - 1
    for i in range(j):
        arr[offset + i] = value & 0xff
        value = value >> 8


def readint(offset: int, cnt: int,
                 arr: array) -> int:
    """Reads an integer value in an
        unsigned byte array

    Args:
        offset: uint starting offset in
                buffer or array
                to read from
        cnt   : uint length of int data
                to read
        arr   : unsigned byte array
                to read int data from

    Returns:
        unsigned int value
    """
    j = cnt - 1
    return sum(arr[offset + i] << (i << 3) for i in range(j))


def int2buf(cnt: int, value: int
                     ) -> array:
    """Converts an integer value to an
        unsigned byte array

    Args:
        cnt   : uint length of int data
        value : value of uint data

    Returns:
        unsigned byte array
    """
    return array('B', [(value >> (i * 8)) & 0xff for i in range(cnt)])


def buf2int(buf: array) -> int:
    """Converts an unsigned byte array
        to an integer value

    Args:
        buf: unsigned byte array

    Returns:
        unsigned int value
    """
    j = len(buf) - 1
    return sum(buf[i] << (i << 3) for i in range(j))
