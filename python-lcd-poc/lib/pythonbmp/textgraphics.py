"""  Text based graphics module
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


def computetextrowandcolcnt(text: str) -> list[int, int]:
    """Returns row and col counts of a multiline text

    Args:
        single of multiline text (string)

    Returns:
        [row: int, col: int]
    """
    l = text.replace('\t', '    ')
    l = l.split('\n') \
        if text.find('\n') > -1 \
        else [l]
    return [len(max(l, key = len)), len(l)]


def plotbitsastext(bits: int):
    """Outputs the bits of byte to
        console

    Args:
        bits: byte value

    Returns:
        space for 0
        *     for 1
    """
    mask = 128
    while mask > 0:
        if (mask & bits) == 0:
            print(' ', end='')
        else:
            print('*', end='')
        mask >>= 1


def plot8bitpatternastext(
        bitpattern: list[int],
        onechar: str,
        zerochar: str):
    """Outputs the bits of a list
        of bytes to console

    Args:
        bitpattern: list of bytes
        onechar   : char to display
                    if bit is 1
        zerochar : char to display
                    if bit is 0

    Returns:
        console output
    """
    s = ""
    for bits in bitpattern:
        mask = 128
        while mask > 0:
            s += onechar if mask & bits > 0 else zerochar
            mask >>= 1
        s += '\n'
    return s
