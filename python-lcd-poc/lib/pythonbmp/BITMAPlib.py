""" Main library include this to use
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
|-----------------------------------|
|   Contact primary author          |
|   if you plan to use this         |
|   in a commercial product at      |
|   joelalcarez1975@gmail.com       |
|-----------------------------------|
|   Educational or hobby use is     |
|   highly encouraged...            |
|   have fun coding !               |
|-----------------------------------|
|   This graphics library outputs   |
|   to a bitmap file.               |
 -----------------------------------
"""

from . import shims
from array import array
from math import(
    sin,
    cos,
    radians,
    pi
    )

from os import (
    TMP_MAX,
    environ,
    getcwd,
    mkdir
    )

from os.path import pathsep
from pathlib import Path
from random import random
from typing import Callable
from numbers import Number
from .proctimer import functimer

from .bmpconstants import(
    bmpheaderid,
    bmpfilesize,
    bmphdrsize,
    bmpcolorbits,
    bmpx,
    bmpy,
    bmppal,
    bmpheadersizedict
    )

from .mathlib import(
    addvect,
    addvectinlist,
    andvect,
    anglebetween2Dlines,
    centerpoint,
    computerotvec,
    cosaffin,
    clamp,
    delta,
    dist,
    enumbits,
    getcomplexdomainbounds,
    getdomainbounds,
    histogram2D,
    histogram2Dcomplex,
    initmatrix,
    intscalarmulvect,
    isinrange,
    lerp,
    invlerp,
    LSMslope,
    LSMYint,
    mirror,
    polar2rectcoord2D,
    range2baseanddelta,
    ravel2D,
    rotatebits,
    roundvect,
    setmax,
    setmin,
    setminmax,
    scalarmulvect,
    smoothstep,
    smootherstep,
    smootheststep,
    smoothsteplerp,
    smoothersteplerp,
    smootheststeplerp,
    invsmoothstep,
    subvect,
    swapxy,
    trans,
    vmag,
    xorvect
    )

from .primitives2D import(
    arcvert,
    bsplinevert,
    circleinvolutevert,
    cornuspiralvert,
    eggcurvevert,
    ellipsevert,
    entirecircleisinboundary,
    entireellipseisinboundary,
    epicycloidvert,
    flowervert,
    gearcurvevert,
    heartcurvevert,
    horizontalvert,
    hypotrochoidvert,
    isinrectbnd,
    iterbeziercurve,
    iterbspline,
    itercircle,
    itercircleinvolute,
    itercirclepart,
    itercirclepartlineedge,
    itercirclepartvertlineedge,
    itercornuspiral,
    iterdrawvec,
    itereggcurve,
    iterellipse,
    iterellipsepart,
    iterellipserot,
    iterepicycloid,
    iterflower,
    itergearcurve,
    itergetneighbors,
    iterheartcurve,
    iterhypotrochoid,
    iterline,
    iterlissajouscurve,
    iterparallelogram,
    iterspirograph,
    itersquircle,
    itersuperellipse,
    lissajouscurvevert,
    listinrecbnd,
    rectboundarycoords,
    recvert,
    regpolygonvert,
    sortrecpoints,
    spiralcontrolpointsvert,
    spirographvert,
    squirclevert,
    superellipsevert,
    verticalvert
    )

from .solids3D import(
    conevertandsurface,
    cubevert,
    cylindervertandsurface,
    decahedvertandsurface,
    fillpolydata,
    gensides,
    getshapesidedict,
    hexahedravert,
    icosahedvertandsurface,
    octahedravert,
    perspective,
    polyboundary,
    rotvec3D,
    spherevertandsurface,
    surfplot3Dvertandsurface,
    tetrahedravert
    )

from .colors import(
    applybrightnessadjtoBGRbuf,
    applycolorfiltertoBGRbuf,
    applygammaBGRbuf,
    applymonochromefiltertoBGRbuf,
    applythresholdadjtoBGRbuf,
    bmpstdpal,
    bmpvalidcolorbits,
    brightnessadjust,
    colorfilter,
    colorfiltertoBGRbuf,
    colormix,
    dichromaticpal,
    gammaBGRbuf,
    gammacorrect,
    getdefaultbitpal,
    getdefaultlumrange,
    interpolateRGB,
    int2BGRarr,
    int2RGB,
    int2RGBarr,
    int2RGBlist,
    invertbitsinbuffer,
    isvalidcolorbit,
    makeBGRbuf,
    matchRGBtopal,
    monochrome,
    monochromefiltertoBGRbuf,
    monochromepal,
    monoshiftablepal,
    monoinverseshiftablepal,
    monoinverseshiftableBGRpal,
    monoinverseshiftableBGRApal,
    probplotRGBto1bit,
    RGB2BGRAlist,
    RGB2BGRbuf,
    RGB2HSL,
    RGB2int,
    RGBfactors2RGB,
    RGBfactorstoBaseandRange,
    RGB2BGRarr,
    thresholdadjust
    )

from .paramchecks import(
    entirecircleinboundary,
    entirerectinboundary,
    func24bitonly,
    func24bitonlyandentirecircleinboundary,
    func24bitonlyandentirerectinboundary,
    func8and24bitonly,
    func4and8bitonly,
    func8and24bitonlyandentirecircleinboundary,
    func8and24bitonlyandentirerectinboundary,
    intcircleparam,
    intcircleparam24bitonly
    )

from .bufferflip import(
    flipnibbleinbuf,
    rotatebitsinbuf
    )

from .chartools import(
    char2int,
    enumletters,
    enumreverseletters
    )

from .conditionaltools import iif, swapif

from .dicttools import dict2descorderlist

from .fileutils import checklink, checklinks

from .fonts import(
    font8x8,
    font8x8thin,
    font8x14,
    font8x14thin,
    font8x8sanserif,
    font8x14sanserif,
    getcharfont)

from .fractals import(
    cliffordattractorlist,
    fractaldomainparamdict,
    fractaldreamattractorlist,
    funcparamdict,
    getIFSparams,
    gumowskimiraattractorlist,
    hilbertvert,
    hopalongattractorlist,
    ikedaattractorlist,
    iterastroid,
    iterbarnsleytree,
    itercosjulia,
    iterIFS,
    iterlamdbafractal,
    iterlemniscate,
    itermarekdragon,
    itermulti2ndtetrationbiomorph,
    itermultibiomorph,
    itermultibiomorphvariant,
    itermultibiomorphphasevariant,
    itermultibrot,
    itermulticircle,
    itermulticorn,
    itermulticosbiomorph,
    itermulticoshbiomorph,
    itermultiexpbiomorph,
    itermultihyperbola,
    itermultijulia,
    itermultisinbiomorph,
    itermultisinhbiomorph,
    itermultisuperellipse,
    itermultitanbiomorph,
    itermultitanhbiomorph,
    itermultizconjugatebiomorph,
    iternewtonsfractal,
    iterngonfractal,
    itersinjulia,
    iterspiraljulia,
    itertetration,
    itertricorn,
    iterxordivfractal,
    iterxorfractal,
    kochcurvevert,
    kochsnowflakevert,
    nattractorlist,
    peterdejongattractorlist,
    symmetriciconattractorlist,
    duffingattractorlist
    )

from .inttools import(
    readint,
    writeint,
    int2buf,
    buf2int
    )

from .messages import sysmsg

from .funcmeta import getfuncmetastr

from .bufresize import(
    adjustxbufsize,
    resizebufNtimesbigger,
    resizesmaller24bitbuf
    )

from .charts import(
    getdatalisttotal,
    genpiechartdata
    )

from .textgraphics import(
    computetextrowandcolcnt,
    plotbitsastext,
    plot8bitpatternastext
    )

from .colordict import(
    getcolorname2RGBdict,
    getcolorname2HSLdict,
    getRGBfactors
    )

from .X11colordict import(
    getX11colorname2RGBdict,
    getX11colorname2HSLdict,
    getX11RGBfactors
    )

from .XKCDcolordict import(
    getXKCDcolorname2RGBdict,
    getXKCDcolorname2HSLdict,
    getXKCDRGBfactors
    )

def _setx(bmp: array, xmax: int):
    """Sets the x value stored in the windows bmp header

    Args:
        bmp: unsigned byte array
             with bmp format
        int: value of x dimension

    Returns:
        byref modified unsigned byte array
    """
    writeint(bmpx, 4, bmp, xmax)


def getmaxx(bmp: array) -> int:
    """Gets the x value stored in the windows bmp header

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        int value of x bmp dimension
    """
    return readint(bmpx, 4, bmp)


def _sety(bmp: array, ymax: int):
    """Sets the y value stored in the windows bmp header

    Args:
        bmp: unsigned byte array
             with bmp format
        int: value of y dimension

    Returns:
        byref modified unsigned byte array
    """
    writeint(bmpy, 4, bmp, ymax)


def getmaxy(bmp: array) -> int:
    """Gets the y value stored in the windows bmp header

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        int value of y bmp dimension
    """
    return readint(bmpy, 4, bmp)


def getmaxxy(bmp: array) -> tuple:
    """Gets the max x and y values stored in the bmp header

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        tuple (x:int,y:int)
    """
    return (readint(bmpx, 4, bmp),
            readint(bmpy, 4, bmp))


def bottomrightcoord(
        bmp: array) -> tuple:
    """Gets the maximum bottom right coordinates of a bmp

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        tuple (x:int,y:int)
    """
    return (readint(bmpx, 4, bmp) - 1,
            readint(bmpy, 4, bmp) - 1)


def centercoord(bmp: array) -> tuple:
    """Gets the central coordinates of a BMP

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        tuple (x:int,y:int)
    """
    return ((readint(bmpx, 4, bmp) - 1) >> 1,
            (readint(bmpy, 4, bmp) - 1) >> 1)


def isinBMPrectbnd(bmp: array,
        x: int, y: int) -> bool:
    """Checks if (x,y) coordinates are within the BMP

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int value
              of location
              in x-axis and y-axis

    Returns:
        True if within bounds
        False if out of bounds
    """
    return (x < readint(bmpx, 4, bmp) and
            y < readint(bmpy, 4, bmp)) and \
           (x > -1 and y > -1)


def listinBMPrecbnd(bmp: array,
        xylist: list) -> bool:
    """Checks if a list of (x, y) coordinates are within the BMP

    Args:
        bmp   : unsigned byte array
                with bmp format
        xylist: list of (x, y)
                coordinates
                to be checked

    Returns:
        True if within bounds
        False if out of bounds
    """
    for (x, y) in xylist:
        if isinBMPrectbnd(
            bmp, x, y) == False:
            break
    return True


def _getclrbits(bmp: array) -> int:
    """Get the bit depth of BMP

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        int value of bit depth
        (1, 4, 8, 24) bits
    """
    return bmp[bmpcolorbits]


def _xchrcnt(bmp: array) -> int:
    """Get the chars or bytes in row of a BMP

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        count of bytes or
        chars in a row (x dim)
    """
    return _xbytes(
                readint(bmpx, 4, bmp),
                    bmp[bmpcolorbits])


def _setflsz(bmp: array, size: int):
    """Set the file size of a BMP

    Args:
        bmp  : unsigned byte array
               with bmp format
        size : unsigned int value
               of size of
               the bmp file

    Returns:
        byref modified unsigned byte array
        with new file size
    """
    writeint(bmpfilesize, 8, bmp, size)


def _flsz(bmp: array) -> int:
    """Get the file size of a win bmp
        from its bmp header

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        unsigned int value of size
        of the bmp header
    """
    return readint(bmpfilesize, 8, bmp)


def _sethdsz(bmp: array, hdsize: int):
    """Set the header size of a win bmp

    Args:
        bmp   : unsigned byte array
                with bmp format
        hdsize: unsigned int value
                of size of the
                bmp header

    Returns:
        byref modified byte array
        with new header size
    """
    writeint(bmphdrsize, 4, bmp, hdsize)


def _hdsz(bmp: array) -> int:
    """Get the header size of a BMP

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        int value of header size
    """
    return readint(bmphdrsize, 4, bmp)


def getmaxcolors(bmp: array) -> int:
    """Get the maximum number of colors supported by a BMP

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        int value
    """
    return 1 << bmp[bmpcolorbits]


def _24bmof(bmp: array,
        x: int, y: int) -> int:
    """Get the offset in a byte array with RGB data given x and y

    Args:
        bmp: unsigned byte array
             with bmp format
        x,y: unsigned int value
             of location in
             x-axis and y-axis

    Returns:
        int value of offset to that data in byte array
    """
    return (x * 3) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
        _xbytes(readint(bmpx, 4, bmp), 24))


def _24bmofhd(bmp: array,
        x: int, y: int) -> int:
    """Get the offset in a byte array
    with RGB data given x and y with bmp header

    Args:
        bmp:  unsigned byte array
              with bmp format
        x, y: unsigned int value
              of location in
              x-axis and y-axis

    Returns:
        int value of offset to that data in byte array
    """
    return (x * 3) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
        _xbytes(readint(bmpx, 4, bmp), 24 )) + 54


def _8bmof(bmp: array,
        x: int, y: int) -> int:
    """Get the offset in a byte array
    with 8 bit color data given x and y data

    Args:
        bmp:  unsigned byte array
              with bmp format
        x, y: unsigned int value
              of location in
              x-axis and y-axis

    Returns:
        int value of offset to
        that data in byte array
    """
    return x + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
            _xbytes(readint(bmpx, 4, bmp), 8))


def _8bmofhd(bmp: array,
        x: int, y: int) -> int:
    """Get the offset in a byte array
    with 8 bit color data given
    x and y with adjustments
    made for a bmp header

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int value
              of location in
              x-axis and y-axis

    Returns:
        int value of offset to
        that data in byte array
    """
    return x + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
             _xbytes(readint(bmpx, 4, bmp), 8)) + \
                 1078


def _4bmof(bmp: array,
        x: int, y: int) -> int:
    """Get the offset in a byte array
    with 4 bit color data given x and y data

    Args:
        bmp:  unsigned byte array
              with bmp format
        x, y: unsigned int value
              of location in
              x-axis and y-axis

    Returns:
        int value of offset to that data in byte array
    """
    return (x >> 1) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
             _xbytes(readint(bmpx, 4, bmp), 4))


def _1bmof(bmp: array,
        x: int, y: int) -> int:
    """Get the offset in a byte array
    with 1 bit color data given x and y data

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int value
              of location in
              x-axis and y-axis

    Returns:
        int value of offset to that data in byte array
    """
    return (x >> 3) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
            _xbytes(readint(bmpx, 4, bmp), 1))


def _4bmofhd(bmp: array,
        x: int, y: int) -> int:
    """Get the offset in a byte array
    with 4 bit color data given
    x and y with adjustments
    made due to a header

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int value
              of location in
              x-axis and y-axis

    Returns:
        int value of offset to that data in byte array
    """
    return (x >> 1) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
            _xbytes(readint(bmpx, 4, bmp), 4)) + \
                118


def _1bmofhd(bmp: array,
        x: int, y: int) -> int:
    """Get the offset in a byte array
        with 1 bit color data given
        x and y with adjustments
        made for a bmp header

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int value
              of location in
              x-axis and y-axis

    Returns:
        int value of offset to that data in byte array
    """
    return (x >> 3) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
            _xbytes(readint(bmpx, 4, bmp), 1)) + \
                62


def _getBMoffhdfunc(bmp: array):
    """Returns the correct function
    to use in computing offsets
    with headers in a given bmp

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        function to compute offsets with headers
    """
    return {24: _24bmofhd,
             8: _8bmofhd,
             4: _4bmofhd,
             1: _1bmofhd}[bmp[bmpcolorbits]]


def _getBMofffunc(bmp: array):
    """Returns the correct function
    to use in computing the offset
    in a given bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        function to compute offsets
    """
    return {24: _24bmof,
             8: _8bmof,
             4: _4bmof,
             1: _1bmof}[bmp[bmpcolorbits]]


def _BMoffset(bmp: array,
        x: int, y: int) -> int:
    """Returns the offset given a bmp and (x, y) coordinates

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int location
              in x-axis and y-axis

    Returns:
        unsigned int offset to data in buffer
    """
    return _getBMofffunc(bmp)(bmp, x, y)


def _BMoffsethd(bmp: array,
        x: int, y: int) -> int:
    """Returns the offset given a bmp
    and (x, y) coordinates with header considered

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int location
              in x-axis and y-axis

    Returns:
        unsigned int offset to data in buffer
    """
    return _getBMoffhdfunc(bmp)(bmp, x, y)


def getmaxxyandbits(bmp: array) -> tuple:
    """Returns bitmap metadata
       (x-dimension, y-dimension, bit depth)

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        (x-dimension, y-dimension, bit depth)
    """
    return (((readint(bmpx, 4, bmp),
              readint(bmpy, 4, bmp)),
              bmp[bmpcolorbits]))


def _xbytes(x: int, bits: int) -> int:
    """Get the number of bytes in
    a bmp row given x dimension and bit depth

    Args:
        x   : unsigned int value of
              x-dimension
        bits: unsigned int value of
              bit depth (1, 4, 8, 24)

    Returns:
        int value of number of bytes in a row
    """
    if bits <= 8:
        xperbyte = 8 // bits
        xbytes, rem = divmod(x, xperbyte)
        if rem > 0:
            xbytes += 1
    if bits == 24:
        xbytes = x * 3
    rem = xbytes & 3
    if rem > 0:
        xbytes = xbytes + (4 - rem)
    return xbytes


def _pdbytes(x: int, bits: int) -> int:
    """Get the number of bytes used to pad for
    32-bit alignment given x dimension and bit depth

    Args:
        x   : unsigned int value of
              x-dimension
        bits: unsigned int value of
              bit depth (1, 4, 8, 24)

    Returns:
        int value of number of pad bytes
    """
    if bits <= 8:
        xperbyte = 8 // bits
        xbytes, rem = divmod(x, xperbyte)
        if rem > 0:
            xbytes += 1
    if bits == 24:
        xbytes = x * 3
    rem = xbytes & 3
    return iif(rem > 0, 4 - rem, 0)


def _getbmflsz(x: int, y: int,
            bits: int) -> int:
    """Computes bitmap file size

    Args:
        x, y: unsigned int value
              of x and y dimensions
        bits: bit depth (1, 4, 8, 24)

    Returns:
        int value of file size
    """
    return _xbytes(x, bits) * y + \
           bmpheadersizedict[bits]


def _bmmeta(x: int, y: int,
        bits: int) -> tuple:
    """Computes bitmap meta data

    Args:
        x, y : unsigned int
               values of the
               x and y dimension
        bits : bit depth (1, 4, 8, 24)

    Returns:
        unsigned int values for
        (filesize, headersize,
        xdim, ydim, bitdepth)
    """
    return (_getbmflsz(x, y, bits),
            bmpheadersizedict[bits], x, y, bits)


def isdefaultpal(bmp: array) -> bool:
    """Checks if bitmap has a default RGB color palette

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        True if default
        False if not default
    """
    return getdefaultbitpal(_getclrbits(bmp)) == \
                getallRGBpal(bmp)


def getBMPimgbytes(bmp: array) -> list:
    """Gets the raw image buffer of a bmp without the header

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        list of unsigned bytes
    """
    return bmp[_hdsz(bmp): _flsz(bmp)]


def setBMPimgbytes(
        bmp: array, buf: array):
    """Sets the raw image buffer of a bitmap

    Args:
        bmp: unsigned byte array
             with bmp format
        buf: array of unsigned bytes

    Returns:
        byref modified  unsigned byte array
    """
    bmp[_hdsz(bmp): _flsz(bmp)] = buf


def setbmppal(bmp: array,
        pallist: list):
    """Sets the RGB palette of a bitmap

    Args:
        bmp    : unsigned byte array
                 with bmp format
        pallist: [(r: byte,
                   g: byte,
                   b: byte), ...]

    Returns:
        byref modified unsigned byte array
    """
    if len(pallist) == getmaxcolors(bmp):
        s = bmppal
        for (r, g, b) in pallist:
            bmp[s: s + 3] = array('B', [b, g, r])
            s += 4


def setbmppalwithBGRlist(bmp: array,
        bgrpallist: list):
    """Sets the RGB palette of a bitmap
    with a BGR list

    Args:
        bmp    : unsigned byte array
                 with bmp format
        bgrpallist: [(b: byte,
                      g: byte,
                      r: byte), ...]

    Returns:
        byref modified unsigned byte array
    """

    if len(bgrpallist) == getmaxcolors(bmp):
        s = bmppal
        for p in bgrpallist:
            bmp[s: s + 3] = array('B', p)
            s += 4


def setbmppalwithBGRAlist(bmp: array,
        bgrapallist: list):
    """Sets the RGB palette of a bitmap
    with a BGRA list

    Args:
        bmp    : unsigned byte array
                 with bmp format
        bgrpallist: [(b: byte,
                      g: byte,
                      r: byte,
                      a: byte), ...]

    Returns:
        byref modified unsigned byte array
    """
    if len(bgrapallist) == getmaxcolors(bmp):
        p = ravel2D(bgrapallist)
        bmp[bmppal: bmppal + len(p)] = array('B', p)


def getallRGBpal(
        bmp: array
        ) -> list[list[int, int, int]]:
    """Gets the RGB palette of a bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        [(r: byte, g: byte, b: byte), ...]
    """
    colors = getmaxcolors(bmp)
    return [getRGBpal(bmp, c) for c in range(colors)]


def getRGBpal(bmp: array,
        c: int) -> list[int, int, int]:
    """Gets the [R, G, B] values
        of color c in a bitmap

    Args:
        bmp: unsigned byte array
             with bmp format
        c  : unsigned int color

    Returns:
        [R: byte, G: byte, B:byte]
    """
    i = bmppal + (c << 2)
    return [bmp[i + 2], bmp[i + 1], bmp[i]]


def setRGBpal(bmp: array, c: int,
        r: int, g: int, b: int):
    """Sets the r,g,b values of color c in a bitmap

    Args:
        bmp    : unsigned byte array
                 with bmp format
        r, g, b: unsigned byte values
                 for red, green and
                 blue
        c      : unsigned int color

    Returns:
        byref modified unsigned byte array
    """
    s = bmppal + (c << 2)
    bmp[s: s + 3] = array('B', [b, g, r])


@func4and8bitonly
def palshift(bmp: array, shift: int = 1):
    """Returns a palette shifted bitmap

    Args:
        bmp    : unsigned byte array
                 with bmp format
    Returns:
        byref modified unsigned byte array
    """
    m = getmaxcolors(bmp)
    shift = abs(shift) % m
    e = (m - shift) * 4 + bmppal
    shift *= 4
    p = bmppal + shift
    q = e + shift
    bmp[e: q], bmp[bmppal: e] = \
    bmp[bmppal: p], bmp[p: q]


@func4and8bitonly
def itercolorcycle(bmp: array):
    """Yields a palette shifted bitmaps
    for a full color cycle
    (4 and 8 bit BMP only)

    Args:
        bmp    : unsigned byte array
                 with bmp format
    Returns:
        byref modified unsigned byte array
    """
    m = getmaxcolors(bmp) - 1
    yield bmp
    for _ in range(m):
        palshift(bmp)
        yield bmp


def colorcycle2dir(ExistingBMPfile: str,
        NewDir: str, overwritedir: bool = False):
    """Saves to a directory
    palette shifted bitmaps
    for a full color cycle
    (4 and 8 bit BMP only)

    Args:
        bmp    : unsigned byte array
                 with bmp format
    Returns:
        Bitmap files in specified directory
    """
    _usefndirsv(ExistingBMPfile, NewDir, itercolorcycle, overwritedir)


def colorhistogram(bmp: array) -> list:
    """Creates a color histogram

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        list sorted in descending order of color frequencies
    """
    d = {}
    for v in iterimagecolor(bmp,
                sysmsg['colorhist'],
                '*', sysmsg['done']):
        c = v[1]
        if c not in d:
            d.setdefault(c,1)
        else:
            d[c] += 1
    return  dict2descorderlist(d)


def makenewpalfromcolorhist(
        chist: list, colors: int,
        similaritythreshold: float) -> list:
    """Creates a new palatte based on a color histogram

    Args:
        chist              : list sorted in
                             descending order
                             of color
                             frequencies
        colors             : maximum colors
                             of new palette
        similaritythreshold: controls how
                             close palette
                             entries can be

    Returns:
        unsigned byte array with bmp format
    """
    newpal = [[0, 0, 0]]
    palcnt = 1
    colorcount = len(chist)
    for i in range(colorcount):
        rgb = int2RGBlist(chist[i][1])
        addcl = True
        for palentry in newpal:
            if dist(palentry, rgb) < \
                similaritythreshold:
                addcl = False
        if addcl:
            print(palcnt, '=' ,rgb)
            newpal.append(rgb)
            palcnt += 1
        if palcnt == colors:
            break
    return newpal


def copyBMPhdr(bmp: array) -> array:
    """Copies the bitmap header of an in-memory bmp
    to a new unsigned byte array

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        unsigned byte array with bmp format
    """
    hdsz = _hdsz(bmp)
    newbmp = array('B', [0] * _flsz(bmp))
    newbmp[0: hdsz] = bmp[:hdsz]
    return newbmp


def copyRGBpal(sourceBMP: array,
        destBMP: array):
    """Copies the RGB palette info from
    a source unsigned byte array to
    a destination unsigned byte array

    Args:
        sourceBMP and destBMP are both
        unsigned byte arrays with bmp format

    Returns:
        byref modified destBMP
        (unsigned byte array)
    """
    l = _hdsz(sourceBMP)
    destBMP[bmppal: l] = \
  sourceBMP[bmppal: l]


def _setmeta(bmpmeta: list) -> array:
    """Creates a new bitmap
        with the properties set
        by bmpmeta to
        a new unsigned byte array

    Args:
        bmpmeta: [filesize, hdrsize,
                  x, y, bits]

    Returns:
        unsigned byte array with bmp format
    """
    bmp = array('B')
    [filesize, hdrsize, x, y, bits] = \
        bmpmeta
    bmp.frombytes(b'\x00' * (filesize))
    bmp[0: 2] = bmpheaderid
    _setflsz(bmp, filesize)
    _sethdsz(bmp, hdrsize)
    _setx(bmp, x)
    _sety(bmp, y)
    bmp[bmpcolorbits] = bits
    bmp[14] = 40 #required
    if bits < 24:
        setbmppal(bmp,
        getdefaultbitpal(bits))
    return bmp


def setnewpalfromsourcebmp(
        sourcebmp: array, newbmp: array,
        similaritythreshold: float) -> list:
    """Copies the RGB palette info from
    a source unsigned byte array
    to a destination unsigned byte array
    (source and destination
    can have different bit depths)

    Args:
        sourceBMP, newBMP  : unsigned
                             byte arrays
                             with
                             bmp format
        similaritythreshold: how close
                             can a color
                             in a palette
                             entry be to
                             another color

    Returns:
        byRef modified newBMP
        (unsigned byte array) and a
        list of new palette entries
        based on source bitmap
    """
    newpal = makenewpalfromcolorhist(
        colorhistogram(sourcebmp),
        getmaxcolors(newbmp),
        similaritythreshold)
    setbmppal(newbmp, newpal)
    return newpal


def RGBpalbrightnessadjust(bmp: array,
        percentadj: float)-> list:
    """Copies the RGB palette info from
    a source unsigned byte array to
    a destination unsigned byte array

    Args:
        bmp       : unsigned byte array
                    with bmp format
        percentadj: signed float
                    brightness adj in %

    Returns:
        list of modified RGB values
    """
    return [brightnessadjust(c,percentadj)
                for c in getallRGBpal(bmp)]


def setBMP2monochrome(bmp: array,
        RGBfactors: list[float, float, float]) -> list:
    """Sets a bitmap to use a monochrome palette

    Args:
        bmp       : unsigned byte array
                    with bmp format
        RGBfactors: (r, g, b)
                    all values range
                    from 0.0 to 1.0

    Returns:
        list of modified RGB values
        byref modified byte array
    """
    newpal = monochromepal(
           _getclrbits(bmp), RGBfactors)
    setbmppal(bmp, newpal)
    return newpal


def setBMP2invshiftedmonochrome(bmp: array,
        RGBfactors: list[float, float, float],
        mult: int = 20, shift: int = 0) -> list:
    """Sets a bitmap to use a monochrome palette

    Args:
        bmp       : unsigned byte array
                    with bmp format
        RGBfactors: (r, g, b)
                    all values range
                    from 0.0 to 1.0

    Returns:
        list of modified BGRA values
        byref modified byte array
    """
    newpal = monoinverseshiftableBGRApal(
           _getclrbits(bmp), RGBfactors, mult, shift)
    setbmppalwithBGRAlist(bmp, newpal)
    return newpal


def setBMP2dichromaticpal(bmp: array,
        c1: int, c2: int,
        ) -> list:
    """Sets a bitmap to use a dichromatic palette

    Args:
        bmp    : unsigned byte array
                 with bmp format
        c1, c2 : color as int packed rgb triplets

    Returns:
        list of modified BGRA values
        byref modified byte array
    """
    newpal = RGB2BGRAlist(dichromaticpal(c1, c2, getmaxcolors(bmp), 3))
    setbmppalwithBGRAlist(bmp, newpal)
    return newpal


def newBMP(x: int, y: int,
        colorbits: int) -> array:
    """Creates a new in-memory bitmap

    Args:
        x, y     : unsigned int values
                   of x and y dims
        colorbits: bit depth
                   (1, 4, 8, 24) bits

    Returns:
        unsigned byte array with bitmap layout
    """
    return _setmeta(
           _bmmeta(x, y, colorbits))


def CopyBMPxydim2newBMP(bmp: array,
        newbits: int) -> array:
    """Creates a new bitmap with the same dimensions as bmp

    Args:
        bmp    : unsigned byte array
                 with bmp format
        newbits: bit depth
                 (1, 4, 8, 24)

    Returns:
        unsigned byte array with bitmap layout
    """
    return newBMP(getmaxx(bmp),
                  getmaxy(bmp),
                  newbits)


@checklink
def loadBMP(filename: str) -> array:
    """Load bitmap to a byte array
       (uncompressed bitmap only)

    Args:
        filename: full path to
                  the file to be loaded

    Returns:
        byte array with bmp file contents
    """
    a = array('B')
    with open(filename, 'rb') as f:
        hd = f.read(2)
        if hd != b'BM':
            print(sysmsg['notBMP'])
        else:
            fsize = char2int(f.read(8))
            if fsize > 54:
                f.seek(0)
                a.frombytes(f.read(fsize))
            else:
                print(sysmsg['notBMP'])
        f.close()
    return a


def ensurePaintWrapperInPath():
    """
	Ensure that the parent directory containing
    the 'mspaint' wrapper is in the
    PATH environment var.
    """
    # Directory containing `mspaint`
    d: Path = Path(__file__).parent.parent
    d_abs = str(d.absolute())
    entries = environ.get("PATH", "").split(pathsep)
    if d_abs not in entries:
        entries.insert(0, d_abs)
        environ["PATH"] = pathsep.join(entries)


def saveBMP(filename: str, bmp: array):
    """Saves bitmap to file

    Args:
        filename: full path to
                  the file to be saved
        bmp     : unsigned byte array
                  with the layout of
                  a bitmap file
    Returns:
        A Bitmap File
    """
    ensurePaintWrapperInPath()
    with open(filename, 'wb') as f:
        f.write(bmp)
        f.close()


def BMPbitBLTput(bmp: array,
        offset: int, arraybuf: array):
    """Sets offset in array to arraybuf

    Args:
        bmp     : unsigned byte array
                  with bmp format
        offset  : unsigned int
                  address in buffer
        arraybuf: unsigned byte array

    Returns:
        byref modified unsigned byte array
    """
    hdrsize = _hdsz(bmp)
    bufsize = len(arraybuf)
    startoff = hdrsize + offset
    endoff = startoff + bufsize
    if offset >= 0 and endoff <= _flsz(bmp):
        bmp[startoff: endoff] = arraybuf
    else:
        print (sysmsg['invalidoffset'])


def BMPbitBLTget(bmp: array,
        offset: int,
        bufsize: int) -> array:
    """Gets [offset: offset + bufsize] to a new array

    Args:
        bmp    : unsigned byte array
                 with bmp format
        offset : unsigned int
                 address in buffer
        bufsize: unsigned int
                 size of buffer

    Returns:
        unsigned byte array
    """
    hdrsize = _hdsz(bmp)
    retval = array('B', [])
    startoff = hdrsize + offset
    endoff = startoff + bufsize
    if (offset >= 0 and bufsize > 0 ) and \
            endoff <= _flsz(bmp):
        retval = bmp[startoff: endoff]
    else:
        print(sysmsg['invalidoffset'])
    return retval


def vertBMPbitBLTget(bmp: array,
        x: int, y1: int,
                y2: int) -> array:
    """Gets vertical slice to a new array

    Args:
        bmp   : unsigned byte array
                with bmp format
        x     : unsigned int x
        y1, y2  and y coordinates

    Returns:
        unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    r = _xchrcnt(bmp)
    m = getmaxxy(bmp)
    c = _getBMoffhdfunc(bmp)
    if isinrange(x, m[0], -1):
        y1, y2 = swapif(y1 ,y2, y1 > y2)
        y1 = setmin(y1, 0)
        y2 = setmax(y2, m[1] - 1)
        s = c(bmp, x, y2)
        e = c(bmp, x, y1) + r
        if bits == 24:
            return makeBGRbuf(bmp[s    : e - 2: r],
                              bmp[s + 1: e - 1: r],
                              bmp[s + 2: e    : r])
        else:
            return bmp[s: e: r]
    else:
        print(sysmsg['lineoutofbnd'])


def _fnwithpar2vertslice(bmp: array,
        x: int, y1: int, y2: int,
        func: Callable, funcparam):
    """Apply a user defined function to vertical slices

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y1, y2 : unsigned int
                    x and y coordinates
        func      : user defined function
        funcparam : parameters of
                    the user defined
                    function

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    r = _xchrcnt(bmp)
    (mx, my) = getmaxxy(bmp)
    c = _getBMoffhdfunc(bmp)
    if isinrange(x, mx, -1):
        y1, y2 = swapif(y1, y2, y1 > y2)
        y1 = setmin(y1, 0)
        y2 = setmax(y2, my - 1)
        s = c(bmp, x, y2)
        e = c(bmp, x, y1) + r
        if bits == 24:
            e2 = e - 2
            bmp[s: e2: r] = \
                func(bmp[s: e2: r],
                        funcparam)
            s1 = s + 1
            e1 = e - 1
            bmp[s1: e1: r] = \
                func(bmp[s1: e1: r],
                         funcparam)
            s2 = s + 2
            bmp[s2: e : r] = \
                func(bmp[s2: e: r],
                        funcparam)
        else:
            bmp[s: e: r] = \
                func(bmp[s: e: r],
                        funcparam)
    else:
        print(sysmsg['lineoutofbnd'])


def plotRGBxybit(bmp: array,
        x: int, y: int, rgb: list):
    """Sets pixel at (x, y) in a bitmap to color [R, G, B]

    Args:
        bmp: unsigned byte array
             with bmp format
        x,y: unsigned int locations
             in x and y
        rgb: color defined by
             [R: byte, G: byte, B: byte]

    Returns:
        byref modified unsigned byte array
    """
    if isinBMPrectbnd(bmp, x, y):
        if bmp[bmpcolorbits] == 24:
            offset = _24bmofhd(bmp, x, y)
            (r, g, b) = rgb
            bmp[offset: offset + 3] = \
                array('B', [b, g, r])
        else:
            plotxybit(bmp,
                x, y, matchRGBtopal(rgb,
                        getallRGBpal(bmp)))


def plotxybit(bmp: array,
        x: int, y: int, c: int):
    """Sets pixel at (x, y) in a bitmap to color c

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int
              locations in x and y
        c   : unsigned int color

    Returns:
        byref modified unsigned byte array
    """
    if not isinBMPrectbnd(bmp, x, y):
        return
    bits = bmp[bmpcolorbits]
    if bits == 24:
        offset = (x * 3) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
        _xbytes(readint(bmpx, 4, bmp), 24 )) + 54
        bmp[offset: offset + 3] = \
            array('B', [c & 0xff,
                      (c >> 8) & 0xff,
                       c >> 16])
    elif bits == 8:
        bmp[x + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
             _xbytes(readint(bmpx, 4, bmp), 8)) + \
                 1078] = c & 0xff
    elif bits == 4:
        offset = (x >> 1) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
            _xbytes(readint(bmpx, 4, bmp), 4)) + \
                118
        c &= 0xf
        if x & 1 == 1:
            bmp[offset] = \
                (bmp[offset] & 0xf0) + c
        else:
            bmp[offset] = \
                (c << 4) + (bmp[offset] & 0xf)
    elif bits == 1:
        offset = (x >> 3) + \
        ((readint(bmpy, 4, bmp) - y - 1) * \
            _xbytes(readint(bmpx, 4, bmp), 1)) + \
                62
        b = bmp[offset]
        mask = 1 << (7 - (x % 8))
        c = 1 if c > 0 else 0
        if c == 1:
           b |= mask
        elif (b & mask) == mask \
        and c == 0:
           b = b ^ mask
        bmp[offset] = b


def getxybit(bmp: array,
        x: int, y: int) -> int:
    """Gets color of pixel at (x, y) in a BMP

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int
              locations in x and y

    Returns:
        unsigned int color value
    """
    retval = 0
    if isinBMPrectbnd(bmp, x, y):
        offset = _BMoffsethd(bmp, x, y)
        bits = bmp[bmpcolorbits]
        if bits == 1:
            mask = 7 - (x % 8)
            retval = \
                (bmp[offset] & (1 << mask)) >> mask
        elif bits == 4:
            retval = (bmp[offset] & 0xf) if x & 1 == 1 else (bmp[offset] & 0xf0) >> 4
        elif bits == 8:
            retval = bmp[offset]
        elif bits == 24:
            retval = \
                RGB2int(bmp[offset + 2],
                        bmp[offset + 1],
                        bmp[offset])
    else:
        retval = -1
    return retval


def getRGBxybitvec(bmp: array,
        v: list) -> list:
    """Gets the RGB of a pixel at (x,y) in a BMP

    Args:
        bmp: unsigned byte array
             with bmp format
        v  : pixel location
             (x:int, y:int)

    Returns:
        [R: byte, G: byte, B: byte]
    """
    return getRGBxybit(bmp, v[0], v[1])


def getRGBxybit(bmp: array,
        x: int, y: int) -> list[int,
                                int,
                                int]:
    """Gets [R:byte, G:byte, B:byte]
    of pixel at (x, y) in a bitmap

    Args:
        bmp : unsigned byte array
              with bmp format
        x, y: unsigned int
              locations in x and y

    Returns:
        [R: byte, G: byte, B: byte]
    """
    retval = []
    if isinBMPrectbnd(bmp, x, y):
        if _getclrbits(bmp) == 24:
            i = _24bmofhd(bmp, x, y)
            retval = [bmp[i + 2],
                      bmp[i + 1],
                      bmp[i    ]]
        else:
            retval = getRGBpal(bmp,
                        getxybit(bmp, x, y))
    return retval


def getxybitvec(bmp: array,
        v: list) -> int:
    """Gets color of pixel at (x, y)

    Args:
        bmp: unsigned byte array
             with bmp format
        v  : (x: int, y: int)
             pixel coordinates

    Returns:
        unsigned int color value
    """
    return getxybit(bmp, v[0], v[1])


def intplotvecxypoint(bmp: array,
        v: list[int, int], c: int):
    """Sets the color of a pixel at (x, y)

    Args:
        bmp: unsigned byte array
             with bmp format
        v  : (x:int, y:int)
             pixel coordinates
        c  : unsigned int
             color value

    Returns:
        byref modified unsigned byte array
    """
    plotxybit(bmp, v[0], v[1], c)


def plotvecxypoint(bmp: array,
        v: list, c: int):
    """Sets the color of a pixel
        at (x, y)

    Args:
        bmp: unsigned byte array
             with bmp format
        v  : (x:float, y:float)
                    or
             (x:int, y:int)
        c  : unsigned int
             color value

    Returns:
        byref modified
        unsigned byte array
    """
    (x, y) = roundvect(v)
    plotxybit(bmp, x, y, c)


def plotRGBxybitvec(bmp: array,
        v: list, rgb: list):
    """Sets [R, G, B] of pixel at (x, y)

    Args:
        bmp: unsigned byte array
             with bmp format
        v  : (x: float, y: float)
        rgb: [R: byte,
              G: byte,
              B: byte]

    Returns:
        byref modified unsigned byte array
    """
    plotRGBxybit(bmp, v[0], v[1], rgb)


def plotxypointlist(bmp: array,
        vlist: list, penradius: int,
        color: int):
    """Draws a circle or a point
    depending on the penradius
    with a given color for
    all points in a point list

    Args:
        bmp      : unsigned byte array
                   with bmp format
        vlist    : [(x: uint, y: uint) ,...]
                   list of points
        penradius: radius of the pen
                   (in pixels)
        color    : color of the pen

    Returns:
        byref modified unsigned byte array
    """
    for v in vlist:
        roundpen(bmp, v,
            penradius, color)


def roundpen(bmp: array, point: list,
        penradius: int, color: int):
    """Draws a circle or a point
    depending on the penradius
    with a given color

    Args:
        bmp      : unsigned byte array
                   with bmp format
        point    : (x:uint,y:uint)
                   centerpoint
        penradius: radius of the
                   pen in pixels
        color    : color of the pen

    Returns:
        byref modified unsigned byte array
    """
    (x, y) = point
    if penradius <= 1:
        plotxybit(bmp, x, y, color)
    else:
        circle(bmp, x, y, penradius,
                        color, True)


def swapcolors(bmp: array,
        p1: list, p2: list):
    """Swaps the colors of two points in a BMP

    Args:
        bmp   : unsigned byte array
                with bmp format
        p1, p2: endpoints of the
                line(x: uint, y: uint)

    Returns:
        byref modified unsigned byte array
    """
    c = getxybitvec(bmp, p1)
    intplotvecxypoint(bmp, p1,
        getxybitvec(bmp, p2))
    intplotvecxypoint(bmp, p2, c)


def line(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int, color: int):
    """Draw a Line in a bitmap

    Args:
        bmp   : unsigned byte array
                with bmp format
        x1, y1: endpoint 1
        x2, y2: endpoint 2
        color : color of the line

    Returns:
        byref modified
        unsigned byte array
    """
    (mx, my) = bottomrightcoord(bmp)
    x1 = setminmax(x1 ,0, mx)
    x2 = setminmax(x2, 0, mx)
    y1 = setminmax(y1, 0, my)
    y2 = setminmax(y2, 0, my)
    if x1 == x2:
        vertline(bmp, x1, y1, y2, color)
    elif y1 == y2:
        horiline(bmp, y1, x1, x2, color)
    else:
        bits = bmp[bmpcolorbits]
        p1 = (x1, y1)
        p2 = (x2, y2)
        c = _getBMoffhdfunc(bmp)
        if bits == 24:
            buf = int2BGRarr(color)
            for (x, y) in iterline(p1, p2):
                s = c(bmp, x, y)
                bmp[s: s + 3] = buf
        elif bits ==  8:
            color &= 0xff
            for (x, y) in iterline(p1, p2):
                bmp[c(bmp, x, y)] = \
                    color
        elif bits == 4:
            if color > 15:
                color &= 0xf
            for (x, y) in iterline(p1, p2):
                s = c(bmp, x, y)
                if x & 1 == 1:
                    bmp[s] = \
                        (bmp[s] & 0xf0) + color
                else:
                    bmp[s] = (color << 4) + \
                             (bmp[s] & 0xf)
        elif bits == 1:
            if color > 1:
                color &= 0x1
            for (x, y) in iterline(p1, p2):
                plotxybit(bmp, x, y, color)


def horiline(bmp: array, y: int,
        x1: int, x2: int, color: int):
    """Draw a Horizontal Line

    Args:
        bmp   : unsigned byte array
                with bmp format
        y     : constant y value
                of the line
        x1    : starts at x1
        x2    : ends at x2
        color : color of the line

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    (mx, my) = getmaxxy(bmp)
    if isinrange(y, my, -1):
        x1, x2 = swapif(x1, x2, x1 > x2)
        x1 = setmin(x1, 0)
        x2 = setmax(x2, mx - 1)
        dx = x2 - x1 + 1
        s = _BMoffsethd(bmp, x1, y)
        if bits == 24:
            e = s + (dx * 3)
            (r, g, b)= int2RGBlist(color)
            bmp[s: e] = array('B', [b, g, r] * dx)
        elif bits == 8:
            e = s + dx
            color &= 0xff
            bmp[s:e] = \
                array('B', [color] * dx)
        elif bits == 4:
            dx >>= 1
            e = s + dx
            color &= 0xf
            c2 = (color << 4) + color
            c1 = (bmp[s] & 0xf0) if x1 & 1 == 1 else -1
            bmp[s: e] = array('B', [c2] * dx)
            if c1 != -1:
                bmp[s] = c1 + (bmp[s] & 0xf)
            if x2 & 1 == 1:
                plotxybit(bmp,
                    x2 - 1, y, color)
        elif bits == 1:
            x2 += 1
            for x in range(x1, x2):
                plotxybit(bmp,
                    x, y, color)
    else:
        print(sysmsg['lineoutofbnd'])


def vertline(bmp: array, x: int,
        y1: int, y2: int, color: int):
    """Draw a Vertical Line

    Args:
        bmp  : unsigned byte array
               with bmp format
        x    : constant x value
               of the line
        y1   : starts at y1
        y2   : ends at y2
        color: color of the line

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    if bits not in [24, 8]:
        y1 , y2 = \
            swapif(y1, y2, y1 > y2)
        y2 += 1
        for y in range(y1, y2):
            plotxybit(bmp, x, y, color)
    else:
        r = _xchrcnt(bmp)
        (mx, my) = getmaxxy(bmp)
        c = _getBMoffhdfunc(bmp)
        if isinrange(x, mx, -1):
            y1, y2 = \
                swapif(y1, y2, y1 > y2)
            y1 = setmin(y1, 0)
            y2 = setmax(y2, my - 1)
            dy = y2 - y1 + 1
            rgb = int2RGBlist(color)
            s = c(bmp, x, y2)
            e = c(bmp, x, y1) + r
            if bits == 24:
                bmp[s: e - 2: r] = \
                    array('B', [rgb[2]] * dy)
                bmp[s + 1: e - 1: r] = \
                    array('B', [rgb[1]] * dy)
                bmp[s + 2: e: r] = \
                    array('B', [rgb[0]] * dy)
            elif bits==8:
                bmp[s: e: r] = \
                    array('B', [color % 256] * dy)
        else:
            print(sysmsg['lineoutofbnd'])


def fillbackgroundwithgrad(bmp: array,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float],
        direction: int):
    """Fills entire bitmap with a linear gradient

    Args:
        bmp       : unsigned byte array
                    with bmp format
        lumrange  : [byte,byte] that
                    define the range
                    the of gradient
        RGBfactors: [r,g,b] each item
                    in list are unsigned
                    floats from 0 to 1
        direction : 0 - vertical
                    1 - horizontal

    Returns:
        byref modified unsigned byte array
    """
    filledgradrect(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1,
        lumrange, RGBfactors,
        direction)


@func8and24bitonlyandentirerectinboundary
def filledgradrect(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float],
        direction: int,
        fn: Callable = lerp):
    """Draw a filled rectangle with a linear gradient

    Args:
        bmp            : unsigned
                         byte array
                         with bmp format
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        lumrange       : [byte, byte]
                         defines
                         the range
                         of the
                         gradient
        RGBfactors     : [r, g, b] items
                         in list are
                         unsigned floats
                         from 0.0 to 1.0
        direction      : 0 - vertical
                         1 - horizontal
        fn             : optional gradient function

    Returns:
        byref modified
        unsigned byte array
    """
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    dx = x2 - x1 + 1
    dy = y2 - y1 + 1
    (r1, g1, b1) = scalarmulvect(RGBfactors, lumrange[0])
    (r2, g2, b2) = scalarmulvect(RGBfactors, lumrange[1])
    if direction == 0:
        xlim = x2 + 1
        for x in range(x1, xlim):
            c = RGB2int(*interpolateRGB(
             fn, r1, r2, g1, g2, b1, b2,
             setminmax(x / dx, 0.0, 1.0)))
            if bmp[bmpcolorbits] != 24:
                c = matchRGBtopal(
                        int2RGBarr(c),
                        getallRGBpal(bmp))
            vertline(bmp, x, y1, y2, c)
    else:
        ylim = y2 + 1
        for y in range(y1, ylim):
            c = RGB2int(*interpolateRGB(
             fn, r1, r2, g1, g2, b1, b2,
             setminmax(y / dy, 0.0, 1.0)))
            if bmp[bmpcolorbits] != 24:
                c = matchRGBtopal(
                        int2RGBarr(c),
                        getallRGBpal(bmp))
            horiline(bmp, y, x1, x2, c)



@func24bitonlyandentirerectinboundary
def filleddichromaticgradrect(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c1: int,
        c2: int,
        direction: int,
        fn: Callable = lerp):
    """Draw a filled rectangle with a dichromatic linear gradient

    Args:
        bmp            : unsigned
                         byte array
                         with bmp format
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        c1, c2         : int packed RGB
        direction      : 0 - vertical
                         1 - horizontal

    Returns:
        byref modified
        unsigned byte array
    """
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    dx = x2 - x1 + 1
    dy = y2 - y1 + 1
    r1, g1, b1 = int2RGB(c1)
    r2, g2, b2 = int2RGB(c2)
    if direction == 0:
        xlim = x2 + 1
        for x in range(x1, xlim):
            vertline(bmp, x, y1, y2,
             RGB2int(*interpolateRGB(
             fn, r1, r2, g1, g2, b1, b2,
             setminmax(x / dx, 0.0, 1.0))))
    else:
        ylim = y2 + 1
        for y in range(y1, ylim):
            horiline(bmp, y, x1, x2,
             RGB2int(*interpolateRGB(
             fn, r1, r2, g1, g2, b1, b2,
             setminmax(y / dy, 0.0, 1.0))))


def fillbackgroundwithdichromaticgrad(bmp: array,
        lumrange: list[int, int],
        c1: int,
        c2: int,
        direction: int):
    """Fills entire bitmap with a dichromatic linear gradient

    Args:
        bmp       : unsigned byte array
                    with bmp format
        lumrange  : [byte,byte] that
                    define the range
                    the of gradient
        RGBfactors: [r,g,b] each item
                    in list are unsigned
                    floats from 0 to 1
        direction : 0 - vertical
                    1 - horizontal

    Returns:
        byref modified unsigned byte array
    """
    filleddichromaticgradrect(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1,
        c1, c2,
        direction)


@entirerectinboundary
def itercopyrect(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int) -> array:
    """Scan a rectangular area and yield scan lines

    Args:
        bmp            : unsigned
                         byte array
                         with bmp format
        x1, y1, x2, y2 : defines the
                         rectangle

    Yields:
        unsigned byte array
        scanlines of the area
    """
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    bufsize = adjustxbufsize(bmp, x1, x2)
    r = _xchrcnt(bmp)
    offset = _BMoffset(bmp, x1, y2)
    lim = _BMoffset(bmp, x2, y1) + r
    while (offset + bufsize) <= lim:
        yield BMPbitBLTget(bmp,
               offset, bufsize)
        offset += r


def intlinevec(bmp: array,
        u: list, v: list,
        color: int):
    """Draw a line in a bitmap

    Args:
        bmp   : unsigned byte array
                with bmp format
        u, v  : (x: int, y: int) points
                that defines the line
        color : color of the line

    Returns:
        byref modified unsigned byte array
    """
    line(bmp, u[0], u[1],
              v[0], v[1], color)


def linevec(bmp: array,
        u: list, v: list,
        color: int):
    """Draw a line in a bitmap

    Args:
        bmp  : unsigned byte array
               with bmp format
        u, v : (x:float,y:float)
               the endpoints
               of the line
        color: the color of the line

    Returns:
        byref modified unsigned byte array
    """
    intlinevec(bmp, roundvect(u),
                    roundvect(v), color)


def filledparallelogram(bmp: array,
        p1: list, p2: list, p3: list,
        color: int):
    """Creates a filled parallelogram defined by 3 points

    Args:
        bmp       : unsigned byte array
                    with bmp format
        p1, p2, p3:(x: float, y: float)
                    points that define
                    a parallelogram
        color     : color of the
                    filled parallelgram

    Returns:
        byref unsigned modified byte array
    """
    for (u, v) in iterparallelogram(
                p1, p2, p3):
        linevec(bmp, u, v, color)


def drawvec(bmp: array,
        u: list, v: list,
        headsize: int,
        color: int,
        penradius: int = None):
    """Draws a vector (line segment with arrow head)

    Args:
        bmp      : unsigned byte array
                   with bmp format
        u        : (x: float, y: float)
                    point 1 origin
        v        : (x: float, y: float)
                   point 2 has arrow
        headsize : size of the arrow
                   0 for default size
        color    : color of the vector
        penradius: optional parameter
                   for thick arrow

    Returns:
        byref modified unsigned byte array
    """
    for (a, b) in iterdrawvec(u, v, headsize):
        if type(penradius) == int:
            thickroundline(bmp, a, b, penradius, color)
        else:
            linevec(bmp, a, b, color)


def thickroundline(bmp: array,
        p1: list, p2: list,
        penradius: int, color: int):
    """Draw a Thick Rounded Line

    Args:
        bmp       : unsigned byte array
                    with bmp format
        p1, p2    : (x, y) endpoints
                    of the line
        penradius : radius of pen
                    in pixels
        color     : color of the line

    Returns:
        byref modified unsigned byte array
    """
    for (x, y) in iterline(p1, p2):
        circle(bmp, x, y,
        penradius, color, True)


def gradthickroundline(bmp: array,
        p1: list, p2: list,
        penradius: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float]):
    """Draw a Thick Rounded Line with a Gradient

    Args:
        bmp      : unsigned byte array
                   with bmp format
        p1, p2   : (x, y) endpoints
                    of the line
        penradius: radius of pen
        lumrange : list of two
                   byte values
                   [gradstart,gradend]
                   that define the
                   luminosity gradient
        RGBfactors: [r, g, b]
                    each item
                    in list is an
                    unsigned float
                    with a range
                    from 0.0 to 1.0

    Returns:
        byref modified
        unsigned byte array
    """
    lum1, lumrang = range2baseanddelta(lumrange)
    for i in range(penradius, 0, -1):
        c = colormix(int(lum1 +
            (lumrang * i / penradius)),
            	RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                    int2RGBarr(c),
                    getallRGBpal(bmp))
        thickroundline(bmp, p1, p2, i, c)


@intcircleparam24bitonly
def _usenopar24btfn2circreg(bmp: array,
        x: int, y: int, r: int,
        func: Callable):
    """Apply a no parameter function
    to a circular region with
    center at (x, y) and a radius r
    in a 24-bit bitmap

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of the circular area
        func   : function to apply

    Returns:
        byref modified unsigned byte array
    """
    c = _getBMoffhdfunc(bmp)
    if entirecircleisinboundary(
            x, y, - 1,
            getmaxx(bmp), -1,
            getmaxy(bmp), r):
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            s1 = c(bmp, x1, y1)
            e1 = c(bmp, x2, y1)
            s2 = c(bmp, x1, y2)
            e2 = c(bmp, x2, y2)
            bmp[s1: e1] = func(bmp[s1: e1])
            if y1 != y2:
                bmp[s2: e2] = func(bmp[s2: e2])
    else:
        xmax = getmaxx(bmp)
        ymax = getmaxy(bmp)
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            x1 = setmin(x1, 0)
            x2 = setmax(x2, xmax - 1)
            if isinrange(y2, ymax, -1):
                s2 = c(bmp, x1 ,y2),
                e2 = c(bmp, x2, y2)
                bmp[s2: e2] = func(bmp[s2: e2])
            if isinrange(y1, ymax, -1) and y2 != y1:
                s1 = c(bmp, x1, y1)
                e1 = c(bmp, x2, y1)
                bmp[s1: e1] = func(bmp[s1: e1])


@intcircleparam24bitonly
def _use24btfn2circreg(bmp: array,
        x: int, y: int, r: int,
        func: Callable, funcparam):
    """Apply function func to a
    circular region with center at
    (x, y) and a radius r that is
    within a 24-bit bitmap

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y, r   : center (x, y)
                    and radius r
                    of the circular area
        func      : function to apply
        funcparam : parameters of the
                    function

    Returns:
        byref modified
        unsigned byte array
    """
    c = _getBMoffhdfunc(bmp)
    if entirecircleisinboundary(
            x, y, -1,
            getmaxx(bmp), -1,
            getmaxy(bmp), r):
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            s1 = c(bmp, x1, y1)
            e1 = c(bmp, x2, y1)
            s2 = c(bmp, x1, y2)
            e2 = c(bmp, x2, y2)
            bmp[s1: e1] = \
                func(bmp[s1: e1], funcparam)
            if y2 != y1:
                bmp[s2: e2] = \
                    func(bmp[s2: e2], funcparam)
    else:
        xmax = getmaxx(bmp)
        ymax = getmaxy(bmp)
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            x1 = setmin(x1, 0)
            x2 = setmax(x2, xmax - 1)
            if isinrange(y2, ymax, -1):
                s2 = c(bmp, x1, y2)
                e2 = c(bmp, x2, y2)
                bmp[s2: e2] = \
                    func(bmp[s2: e2], funcparam)
            if isinrange(y1, ymax, -1) and y2 != y1:
                s1 = c(bmp, x1, y1)
                e1 = c(bmp, x2, y1)
                bmp[s1: e1] = \
                    func(bmp[s1: e1], funcparam)


@entirecircleinboundary
def copycircregion2buf(bmp: array,
        x: int, y: int,
        r: int) -> list:
    """Copies a circular region to a
    buffer which is defined by
    centerpoint at (x, y) and radius r

    Args:
        bmp     : unsigned byte array
                  with bmp format
        x, y, r : center (x, y)
                  and radius r
                  of the circular area

    Returns:
        list with buffer of circular region
    """
    copybuf = [_getclrbits(bmp), r]
    c = _getBMoffhdfunc(bmp)
    for (v0, v1) in itercirclepartlineedge(r):
        x1, x2 = mirror(x, v0)
        y1, y2 = mirror(y, v1)
        copybuf += [[bmp[c(bmp, x1, y1): c(bmp, x2, y1)],
                     bmp[c(bmp, x1, y2): c(bmp, x2, y2)]]]
    return copybuf


def pastecirularbuf(bmp: array,
        x: int, y: int, circbuf: list):
    """Paste a circular buffer with a
    given radius to a centerpoint at (x, y)

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y   : center of circular
                 region
        circbuf: list generated by
                 copycircregion2buf

    Returns:
        byref modified unsigned byte array
    """
    if circbuf is None:
        print(sysmsg['invalidbuf'])

    else:
        r = circbuf[1]
        c = _getBMoffhdfunc(bmp)
        if entirecircleisinboundary(
                x, y, -1,
                getmaxx(bmp), -1,
                getmaxy(bmp), r):
            if _getclrbits(bmp) == circbuf[0]:
                for i, (v0, v1) in enumerate(itercirclepartlineedge(r), start=2):
                    x1, x2 = mirror(x, v0)
                    y1, y2 = mirror(y, v1)
                    bmp[c(bmp, x1, y1): c(bmp, x2, y1)] = \
                        circbuf[i][0]
                    bmp[c(bmp, x1, y2): c(bmp, x2, y2)] = \
                        circbuf[i][1]
            else:
                raise(sysmsg['bitsnotequal'])
        else:
            print(sysmsg['regionoutofbounds'])


def copycircregion(bmp: array,
        x: int, y: int, r: int,
        newxy: list):
    """Copy a circular buffer with center at (x, y)
    and a radius r to a new area with centerpoint at
    newxy [x, y]

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of the
                 circular region
        newxy :  center of
                 circular area
                 to paste
                 the buffer into

    Returns:
        byref modified unsigned byte array
    """
    pastecirularbuf(
      bmp, newxy[0], newxy[1],
      copycircregion2buf(bmp, x, y, r))


@intcircleparam
def _usenoparfn2circreg(bmp: array,
        x: int, y: int, r: int,
        func: Callable):
    """Apply a no parameter function
    to a circular region with a
    center at (x,y) and a radius r

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x,y)
                 and radius r
                 of the
                 circular area
        func   : function
                 to apply
                 to the area

    Returns:
        byref modified unsigned byte array
    """
    c = _getBMoffhdfunc(bmp)
    if entirecircleisinboundary(
            x, y, -1, getmaxx(bmp),
                  -1, getmaxy(bmp), r):
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            s1 = c(bmp, x1, y1)
            e1 = c(bmp, x2, y1)
            s2 = c(bmp, x1, y2)
            e2 = c(bmp, x2, y2)
            bmp[s1: e1] = func(bmp[s1: e1])
            if y1 != y2:
                bmp[s2: e2] = func(bmp[s2: e2])
    else:
        xmax = getmaxx(bmp)
        ymax = getmaxy(bmp)
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            x1 = setmin(x1, 0)
            x2 = setmax(x2, xmax - 1)
            if isinrange(y2, ymax, -1):
                s2 = c(bmp, x1, y2)
                e2 = c(bmp, x2, y2)
                bmp[s2: e2] = func(bmp[s2: e2])
            if isinrange(y1, ymax, -1) and y2 != y1:
                s1 = c(bmp, x1, y1)
                e1 = c(bmp, x2, y1)
                bmp[s1: e1] = func(bmp[s1: e1])


@entirecircleinboundary
def flipXYcircregion(bmp: array,
        x: int, y: int, r: int):
    """Flip the X and Y coordinates of a circular area

    Args:
        bmp     : unsigned byte array
                  with bmp format
        x, y, r : center (x,y) and
                  radius r of region

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    c = _getBMoffhdfunc(bmp)
    if bits not in [24, 8]:
        n = flipXY(bmp)
        for (v0, v1) in itercirclepartlineedge(r):
            x3, x4 = mirror(y, v0)
            y3, y4 = mirror(x, v1)
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            bmp[c(bmp, x1, y1): c(bmp, x2, y1)] = \
              n[c(n, x3, y3): c(n, x4, y3)]
            bmp[c(bmp, x1, y2): c(bmp, x2, y2)] = \
              n[c(n, x3, y4): c(n, x4, y4)]
    else:
        buf = []
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            x2 += 1
            x3, x4 = mirror(x, v1)
            y3, y4 = mirror(y, v0)
            buf += [[x1, y1, x2, y2,
                     vertBMPbitBLTget(
                      bmp, x3, y3, y4),
                     vertBMPbitBLTget(
                      bmp, x4, y3, y4)]]
        for b in buf:
            x1, y1, x2, y2 = b[0: 4]
            bmp[c(bmp, x1, y1): c(bmp, x2, y1)] = b[4]
            bmp[c(bmp, x1, y2): c(bmp, x2, y2)] = b[5]


def fliphoricircregion(bmp: array,
        x: int, y: int, r: int):
    """Flips horizontally a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of circular region

    Returns:
        byref modified
        unsigned byte array
    """
    horitransformincircregion(
        bmp, x, y, r, 'F')


@func8and24bitonlyandentirecircleinboundary
def horitransformincircregion(
        bmp: array,
        x: int, y: int, r: int,
        trans: str):
    """Horizontal transform to a circular area

    Args:
        bmp  :   unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r of region
        trans:   single letter
                 transform code
                 'L' -> mirror left
                 'R' -> mirror right
                 'F' -> flip

    Returns:
        byref modified
        unsigned byte array
    """
    def _24F():
        bmp[s1: e1 - 2: k], bmp[s2: e2 -2: k], bmp[s1 + 1: e1-1: k], bmp[s2 + 1: e2-1: k], bmp[s1 + 2: e1: k], bmp[s2 + 2: e2: k] = \
        bmp[s2: e2 - 2: k], bmp[s1: e1 - 2: k], bmp[s2 + 1: e2 - 1: k], bmp[s1 + 1: e1 - 1: k], bmp[s2 + 2: e2: k], bmp[s1 + 2: e1: k]

    def _8F():
        bmp[s1: e1: k], bmp[s2: e2: k] = \
        bmp[s2: e2: k], bmp[s1: e1: k]

    def _24R():
        bmp[s1: e1 - 2: k], bmp[s1 + 1: e1 -1: k], bmp[s1 + 2: e1: k] = \
        bmp[s2: e2 - 2: k], bmp[s2 + 1: e2 - 1: k], bmp[s2 + 2: e2: k]

    def _8R():
        bmp[s1: e1: k] = bmp[s2: e2: k]

    def _24L():
        bmp[s2: e2 - 2: k], bmp[s2 + 1: e2 - 1: k], bmp[s2 + 2: e2: k] = \
        bmp[s1: e1 - 2: k], bmp[s1 + 1: e1 - 1: k], bmp[s1 + 2: e1: k]

    def _8L():
        bmp[s2: e2: k] = bmp[s1: e1: k]

    c = _getBMoffhdfunc(bmp)

    f ={8: {'F': _8F,
            'R': _8R,
            'L': _8L
             },
        24: {'F': _24F,
             'R': _24R,
             'L': _24L
             }
        }[bmp[bmpcolorbits]][trans]

    for v in itercirclepartvertlineedge(r):
        x1, x2 = mirror(x, v[0])
        y1, y2 = mirror(y, v[1])
        k = _xchrcnt(bmp)
        s1 = c(bmp, x1, y2)
        e1 = c(bmp, x1, y1) + k
        s2 = c(bmp, x2, y2)
        e2 = c(bmp, x2, y1) + k
        f()


def mirrorleftincircregion(
        bmp: array,
        x: int, y: int, r: int):
    """Mirrors the top-left of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y) and
                 radius r of region

    Returns:
        byref modified unsigned byte array
    """
    horitransformincircregion(
        bmp, x, y, r, 'L')


def mirrorrightincircregion(
        bmp: array,
        x: int, y: int, r: int):
    """Mirrors the right-half of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of a region

    Returns:
        byref modified unsigned byte array
    """
    horitransformincircregion(
        bmp, x, y, r, 'R')


def flipvertcircregion(bmp: array,
        x: int, y: int, r: int):
    """Vertical Flip of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of a region

    Returns:
        byref modified unsigned byte array
    """
    verttransformincircregion(
        bmp, x, y, r, 'F')


@entirecircleinboundary
def verttransformincircregion(
        bmp: array,
        x: int, y: int, r: int,
        trans: str):
    """Applies a vertical transform to
    a circular region with a center
    at (x, y) and radius r

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of region
        trans  :single letter
                transform code
                'T' mirror top
                'B' mirror bottom
                'F' flip

    Returns:
        byref modified unsigned byte array
    """
    def _T():
        bmp[s2: e2] = bmp[s1: e1]

    def _B():
        bmp[s1: e1] = bmp[s2: e2]

    def _F():
        bmp[s1: e1], bmp[s2: e2] = \
        bmp[s2: e2], bmp[s1: e1]

    f = {'T': _T,
         'B': _B,
         'F' : _F}[trans]

    c=_getBMoffhdfunc(bmp)
    for (v0, v1) in itercirclepartlineedge(r):
        x1, x2 = mirror(x, v0)
        y1, y2 = mirror(y, v1)
        s1 = c(bmp, x1, y1)
        e1 = c(bmp, x2, y1)
        s2 = c(bmp, x1, y2)
        e2 = c(bmp, x2, y2)
        f()


def mirrortopincircregion(
        bmp: array,
        x: int, y: int, r: int):
    """Mirror the top-half of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r of area

    Returns:
        byref modified unsigned byte array
    """
    verttransformincircregion(
        bmp, x, y, r, 'T')


def mirrorbottomincircregion(
        bmp: array,
        x: int, y: int, r: int):
    """Mirror the bottom-half of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of region

    Returns:
        byref modified unsigned byte array
    """
    verttransformincircregion(
        bmp, x, y, r, 'B')


def mirrortopleftincircregion(
        bmp: array,
        x: int, y: int, r: int):
    """Mirror the top-left of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of region

    Returns:
        byref modified unsigned byte array
    """
    mirrorleftincircregion(bmp, x, y, r)
    mirrortopincircregion(bmp, x, y, r)


def mirrortoprightincircregion(
        bmp: array,
        x: int, y: int, r: int):
    """Mirror the top-right of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y) and
                 radius r of region

    Returns:
        byref modified unsigned byte array
    """
    mirrorrightincircregion(bmp, x, y, r)
    mirrortopincircregion(bmp, x, y, r)


def mirrorbottomleftincircregion(
        bmp: array,
        x: int, y: int, r: int):
    """Mirror the bottom-left of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x,y) and
                 radius r of region

    Returns:
        byref modified unsigned byte array
    """
    mirrorleftincircregion(bmp, x, y, r)
    mirrorbottomincircregion(bmp, x, y, r)


def mirrorbottomrightincircregion(
        bmp: array,
        x: int, y: int, r: int):
    """Mirror the bottom-right of a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x,y) and
                 radius r of region

    Returns:
        byref modified unsigned byte array
    """
    mirrorrightincircregion(bmp, x, y, r)
    mirrorbottomincircregion(bmp, x, y, r)


@func24bitonlyandentirecircleinboundary
def vertbrightnessgrad2circregion(
        bmp: array,
        x: int, y: int, r: int,
        lumrange: list[int, int]):
    """Vertical brightness gradient adjustment to a circular area

    Args:
        bmp     : unsigned byte array
                  with bmp format
        x, y, r : center (x,y)
                  and radius r
        lumrange: [byte,byte]
                  that define
                  the range of
                  luminosity

    Returns:
        byref modified unsigned byte array
    """
    c = _getBMoffhdfunc(bmp)
    f =  applybrightnessadjtoBGRbuf
    l = lumrange[0]
    dl = delta(lumrange) / (2 * r)
    b = y - r
    for (v0, v1) in itercirclepartlineedge(r):
        x1, x2 = mirror(x, v0)
        y1, y2 = mirror(y, v1)
        l1 = l + (y1 - b) * dl
        l2 = l + (y2 - b) * dl
        s1 = c(bmp, x1, y1)
        e1 = c(bmp, x2, y1)
        s2 = c(bmp, x1, y2)
        e2 = c(bmp, x2, y2)
        bmp[s1: e1], bmp[s2: e2] = \
        f(bmp[s1: e1], l1), f(bmp[s2: e2], l2)


@func24bitonlyandentirecircleinboundary
def horibrightnessgrad2circregion(
        bmp: array,
        x: int, y: int, r: int,
        lumrange: list[int, int]):
    """Horizontal brightness gradient adjustment to a circular area

    Args:
        bmp     : unsigned byte array
                  with bmp format
        x, y, r : center (x, y)
                  and radius r
                  of a circular area
        lumrange: [byte, byte] the
                  luminosity range

    Returns:
        byref modified unsigned byte array
    """
    f =  applybrightnessadjtoBGRbuf
    l = lumrange[0]
    dl = delta(lumrange) / (2 * r)
    b = x - r
    for (v0, v1)  in itercirclepartvertlineedge(r):
        x1, x2 = mirror(x, v0)
        y1, y2 = mirror(y, v1)
        _fnwithpar2vertslice(bmp,
            x1, y1, y2,
            f, l + (x1 - b) * dl)
        if x1 != x2:
            _fnwithpar2vertslice(bmp,
                x2, y1, y2,
                f, l + (x2 - b) * dl)

@intcircleparam
def outlinecircregion(bmp: array,
        x: int, y: int, r: int):
    """Outlines a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of the region

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    c = _getBMoffhdfunc(bmp)
    if entirecircleisinboundary(
            x, y,
            -1, getmaxx(bmp),
            -1, getmaxy(bmp), r):
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            s1 = c(bmp, x1, y1)
            e1 = c(bmp, x2, y1)
            s2 = c(bmp, x1, y2)
            e2 = c(bmp, x2, y2)
            if bits == 24:
                bmp[s1: e1] = \
                array('B', xorvect(
                      bmp[s1: e1],
                      bmp[s1 + 3: e1 + 3]))
                if y1 != y2:
                    bmp[s2: e2] = \
                    array('B', xorvect(
                          bmp[s2: e2],
                          bmp[s2 + 3: e2 + 3]))
            else:
                bmp[s1: e1] = \
                array('B', xorvect(
                      bmp[s1: e1],
                      bmp[s1 + 1: e1 + 1]))
                if y1 != y2:
                    bmp[s2: e2] = \
                    array('B', xorvect(
                          bmp[s2: e2],
                          bmp[s2 + 1: e2 + 1]))
    else:
        (xmax, ymax) = getmaxxy(bmp)
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            x1 = setmin(x1, 0)
            x2 = setmax(x2, xmax - 1)
            if isinrange(y2, ymax, -1):
                s2 = c(bmp, x1, y2)
                e2 = c(bmp, x2, y2)
                if bits == 24:
                    bmp[s2: e2] = \
                    array('B',
                     xorvect(bmp[s2: e2],
                     bmp[s2 + 3: e2 + 3]))
                else:
                    bmp[s2: e2] = \
                    array('B',
                      xorvect(bmp[s2: e2],
                      bmp[s2 + 1: e2 + 1]))
            if isinrange(y1, ymax, -1) and y2 != y1:
                s1 = c(bmp, x1, y1)
                e1 = c(bmp, x2, y1)
                if bits == 24:
                    bmp[s1: e1] = \
                        array('B',
                          xorvect(bmp[s1: e1],
                          bmp[s1 + 3: e1 + 3]))
                else:
                    bmp[s1: e1] = \
                        array('B', xorvect(bmp[s1: e1],
                                   bmp[s1 + 1: e1 + 1]))


def monocircle(bmp: array,
        x: int, y: int, r: int):
    """Monochrome filter to a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of the region

    Returns:
        byref modified unsigned byte array
    """
    _usenopar24btfn2circreg(bmp, x, y, r,
        monochromefiltertoBGRbuf)


def colorfiltercircregion(bmp: array,
        x: int, y: int, r: int,
        rgbfactors: list[float, float, float]):
    """Color Filter to a circular area

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y, r   : center (x, y)
                    and radius r
                    of the region
        rgbfactors: [r, g, b] range of
                    r, g and b are from
                    0.0 min to 1.0 max

    Returns:
        byref modified unsigned byte array
    """
    _use24btfn2circreg(bmp, x, y, r,
        colorfiltertoBGRbuf,
        rgbfactors)


def gammacorrectcircregion(
        bmp: array,
        x: int, y: int, r: int,
        gamma: float):
    """Gamma correction to a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of the region
        gamma  : float value of the
                 gamma adjustment

    Returns:
        byref modified unsigned byte array
    """
    _use24btfn2circreg(bmp, x, y, r,
        gammaBGRbuf, gamma)


def brightnessadjcircregion(
        bmp: array,
        x: int, y: int, r: int,
        percentadj: float):
    """Brightness adjustment to a circular area

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y, r   : center (x, y)
                    and radius r
                    of the region
        percentadj: float percent of the
                    brightness adjustment

    Returns:
        byref modified unsigned byte array
    """
    _use24btfn2circreg(bmp, x, y, r,
         applybrightnessadjtoBGRbuf, percentadj)


def invertbitsincircregion(bmp: array,
        x: int, y: int, r: int):
    """Inverts the bits in a circular area

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x,y)
                 and radius r
                 of the region

    Returns:
        byref modified unsigned byte array
    """
    _usenoparfn2circreg(
        bmp, x, y, r,  invertbitsinbuffer)


def circlevec(bmp: array,
        v: list, r: int, color: int,
        isfilled: bool = None):
    """Draws a circle

    Args:
        bmp     : unsigned byte array
                  with bmp format
        v       : (x, y) center of
                  the circular region
        r       : radius of the
                  circular region
        color   : color of the circle
        isfilled: toggles if the
                  circle is filled

    Returns:
        byref modified unsigned byte array
    """
    (x, y) = roundvect(v)
    circle(bmp, x, y, r, color, isfilled)


def filledcircle(bmp: array,
        x: int, y: int, r: int,
        color: int):
    """Draws a Filled Circle

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of filled circle
        color  : color of the
                 filled circle

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    if bits < 8:
        for (v0, v1) in itercirclepartlineedge(r):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            horiline(bmp, y1, x1, x2, color)
            horiline(bmp, y2, x1, x2, color)
    else:
        (xmax, ymax) = getmaxxy(bmp)
        xmax -= 1
        if bits == 24:
            for (v0, v1) in itercirclepartlineedge(r):
                x1, x2 = mirror(x, v0)
                y1, y2 = mirror(y, v1)
                x1 = setmin(x1, 0)
                x2 = setmax(x2, xmax)
                dx = x2 - x1 + 1
                rgb = int2RGBlist(color)
                colorbuf = array('B', [rgb[2],
                                       rgb[1],
                                       rgb[0]] * dx)
                lbuf = dx * 3
                if isinrange(y2, ymax, -1):
                    s = _24bmofhd(bmp, x1, y2)
                    bmp[s: s + lbuf] = colorbuf
                if isinrange(y1, ymax, -1):
                    s = _24bmofhd(bmp, x1, y1)
                    bmp[s: s+ lbuf] = colorbuf
        elif bits == 8:
            for (v0, v1) in itercirclepartlineedge(r):
                x1, x2 = mirror(x, v0)
                y1, y2 = mirror(y, v1)
                x1 = setmin(x1, 0)
                x2 = setmax(x2, xmax)
                dx = x2 - x1 + 1
                colorbuf = \
                    array('B', [color & 0xff] * dx)
                if isinrange(y2, ymax, -1):
                    s = _8bmofhd(bmp, x1, y2)
                    bmp[s: s + dx] = colorbuf
                if isinrange(y1, ymax, -1):
                    s = _8bmofhd(bmp, x1, y1)
                    bmp[s: s + dx] = colorbuf


def circle(bmp: array,
        x: int, y: int, r: int,
        color: int,
        isfilled: bool = None):
    """Draws a Circle

    Args:
        bmp     : unsigned byte array
                  with bmp format
        x, y, r : center (x,y)
                  and radius r
                  of the circle
        color   : color of the circle
        isfilled: toggles if the
                  circle is filled
                  True -> filled circle
                  False -> circle outline

    Returns:
        byref modified unsigned byte array
    """
    if isfilled:
        filledcircle(bmp, x, y, r,
                            color)
    else:
        (mx, my) = getmaxxy(bmp)
        bits = bmp[bmpcolorbits]
        c = _getBMoffhdfunc(bmp)
        dobndcheck = \
            not entirecircleisinboundary(
                    x, y, -1, mx,
                          -1, my, r)
        if bits == 24:
            color = int2BGRarr(color)
            if dobndcheck:
                for (px, py) in itercircle(x, y, r):
                    if isinBMPrectbnd(
                            bmp, px, py):
                        s = c(bmp, px, py)
                        bmp[s: s + 3] = color
            else:
                for (px, py) in itercirclepart(r):
                    x1, x2 = mirror(x, px)
                    y1, y2 = mirror(y, py)
                    s1 = c(bmp, x1, y1)
                    s2 = c(bmp, x2, y2)
                    s3 = c(bmp, x1, y2)
                    s4 = c(bmp, x2, y1)
                    bmp[s1: s1 + 3] = \
                    bmp[s2: s2 + 3] = \
                    bmp[s3: s3 + 3] = \
                    bmp[s4: s4 + 3] = color
        elif bits == 8:
            if dobndcheck:
                for (px, py) in itercircle(x, y, r):
                    if isinBMPrectbnd(bmp, px, py):
                        bmp[c(bmp, px, py)] = color
            else:
                for (px, py) in itercirclepart(r):
                    x1, x2 = mirror(x, px)
                    y1, y2 = mirror(y, py)
                    bmp[c(bmp, x1, y1)] = \
                    bmp[c(bmp, x2, y2)] = \
                    bmp[c(bmp, x1, y2)] = \
                    bmp[c(bmp, x2, y1)] = color
        else:
            for (px, py) in itercircle(x, y, r):
                plotxybit(bmp, px, py, color)



def squircle(bmp: array,
        x: int, y: int, r: int,
        color: int):
    """Draws a Squircle an intermediate
    between a square and a circle.

    The word "squircle" is a portmanteau
    of the words "square" and "circle".

    Squircles have been applied in design
    and optics.

    Args:
        bmp     : unsigned byte array
                  with bmp format
        x, y, r : center (x,y)
                  and radius r
                  of the circle
        color   : color of the circle
        isfilled: toggles if the
                  circle is filled
                  True -> filled circle
                  False -> circle outline

    Returns:
        byref modified unsigned byte array
    """
    plotpoly(bmp, squirclevert(x, y, r), color)


def thickcircle(bmp: array,
        x: int, y: int, r: int,
        penradius: int, color: int):
    """Draws a Thick Circle

    Args:
        bmp      : unsigned byte array
                   with bmp format
        x, y, r  : center (x, y)
                   and radius r
        penradius: radius of round pen
        color    : color of the circle

    Returns:
        byref modified unsigned byte array
    """
    for (p0, p1) in itercircle(x, y, r):
        circle(bmp, p0, p1,
        penradius, color, True)


def gradthickcircle(bmp: array,
        x: int, y: int, r: int,
        penradius: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float]):
    """Thick Circle with a Gradient

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y, r   : center (x, y)
                    and radius r
        penradius : radius of the
                    round pen
        lumrange  : [byte,byte] range
                    of the luminosity
                    gradient
        rgbfactors: [r,g,b] range are
                    from 0.0 to 1.0

    Returns:
        byref modified unsigned byte array
    """
    lum1, lumrang = range2baseanddelta(lumrange)
    for i in range(penradius, 0, -1):
        c = colormix(int(
                lum1 + (lumrang * i / penradius)), RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                    int2RGBarr(c), getallRGBpal(bmp))
        thickcircle(bmp, x, y, r, i, c)


def gradcircle(bmp: array,
        x: int, y: int, r: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float],
        fn = lerp):
    """Filled Circle with Gradient

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y, r   : center (x, y)
                    and radius r
        lumrange  : [byte,byte] range of
                    the gradient
        rgbfactors: [r, g, b] range are
                    from 0.0 to 1.0

    Returns:
        byref modified unsigned byte array
    """
    (a, b) = lumrange
    for i in range(r - 1, 0, -1):
        c = colormix(
                int(fn(a, b, i / r)),
                RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                    int2RGBarr(c),
                    getallRGBpal(bmp))
        thickcircle(bmp, x, y, i, 2, c)


@func24bitonly
def graddichromiccircle(bmp: array,
        x: int, y: int, r: int,
        c1: int, c2: int, fn: Callable = lerp):
    """Filled Circle with dichromic Gradient

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y, r   : center (x, y)
                    and radius r
        c1, c2    : int RGB packed color

    Returns:
        byref modified unsigned byte array
    """
    r1, g1, b1 = int2RGB(c1)
    r2, g2, b2 = int2RGB(c2)
    for i in range(r - 1, 0, -1):
        thickcircle(bmp, x, y, i, 2,
        RGB2int(*interpolateRGB(fn,
        r1, r2, g1, g2, b1, b2, i / r)))


def orb(bmp: array,
        x: int, y: int, r: int,
        RGBfactors: list[float, float, float]):
    """Draw a Glowy Orb

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y, r   : center (x, y)
                    and radius r
        rgbfactors: [r, g, b] range are
                    from 0.0 to 1.0

    Returns:
        byref modified unsigned byte array
    """
    j = r >> 1
    c = colormix(255, RGBfactors)
    if bmp[bmpcolorbits] != 24:
        c = matchRGBtopal(
                int2RGBarr(c),
                getallRGBpal(bmp))
    filledcircle(bmp, x, y, j, c)
    for i in range(r - 1, j, -1):
        c = colormix(
                int(255 * (r - i) / j),
                RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                    int2RGBarr(c),
                    getallRGBpal(bmp))
        thickcircle(bmp, x, y, i, 2, c)


def thickellipserot(bmp: array,
        x: int, y: int,
        b: int, a: int, degrot: float,
        penradius: int, color: int):
    """Draws a Thick Ellipse

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : center of ellipse
        b, a      : major and minor axis
        degrot    : rotation of
                    the ellipse in degrees
        penradius : thickness of the pen
        color     : color of the ellipse

    Returns:
        byref modified
        unsigned byte array
    """
    for (p0, p1) in iterellipserot(x, y,
                b, a, degrot):
        circle(bmp, p0, p1,
            penradius, color, True)


def gradthickellipserot(bmp: array,
        x: int, y: int, b: int, a: int,
        degrot: float, penradius: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float], fn = lerp):
    """Thick Ellipse with a Gradient fill

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : center of ellipse
        b, a      : major amd minor axis
        degrot    : rotation of
                    the ellipse in degrees
        penradius : defines the
                    thickness of the pen
        lumrange  : [byte:byte] sets
                    the range of the
                    luminosity gradient
        rgbfactors: [r, g, b] range are
                    from 0.0 min to 1.0 max
        fn        : optional gradient function

    Returns:
        byref modified unsigned byte array
    """
    (l1, l2) = lumrange
    for i in range(penradius, 0, -1):
        c = colormix(int(lerp(l1, l2, i / penradius)),
                RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                    int2RGBarr(c),
                    getallRGBpal(bmp))
        thickellipserot(bmp,
            x, y, b, a, degrot, i, c)


def filledellipse(bmp: array,
        x: int, y: int, b: int, a: int,
        color: int):
    """Filled Ellipse

    Args:
        bmp  : unsigned byte array
               with bmp format
        x, y : center of ellipse
        b, a : major amd minor axis
        color: color of the ellipse

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    if bits < 8:
        for (v0, v1) in iterellipsepart(b, a):
            x1, x2 = mirror(x, v0)
            y1, y2 = mirror(y, v1)
            horiline(
                bmp, y1, x1, x2, color)
            horiline(
                bmp, y2, x1, x2, color)
    else:
        (mx, my) = getmaxxy(bmp)
        if bits == 24:
            for (v0, v1) in iterellipsepart(b, a):
                x1, x2 = mirror(x, v0)
                y1, y2 = mirror(y, v1)
                x1 = setmin(x1, 0)
                x2 = setmax(x2, mx - 1)
                dx = x2 - x1 + 1
                rgb = int2RGBlist(color)
                colorbuf = \
                 array('B', [rgb[2],
                             rgb[1],
                             rgb[0]] * dx)
                lbuf = dx * 3
                if isinrange(y2, my, -1):
                    s = _24bmofhd(bmp, x1, y2)
                    bmp[s: s + lbuf] = colorbuf
                if isinrange(y1, my, -1):
                    s = _24bmofhd(bmp, x1, y1)
                    bmp[s: s + lbuf] = colorbuf
        elif bits == 8:
            for (v0, v1) in iterellipsepart(b, a):
                x1, x2 = mirror(x, v0)
                y1, y2 = mirror(y, v1)
                x1 = setmin(x1, 0)
                x2 = setmax(x2, mx - 1)
                dx = x2 - x1 + 1
                colorbuf = \
                    array('B', [color & 0xff] * dx)
                if isinrange(y2, my, -1):
                    s = _8bmofhd(bmp, x1, y2)
                    bmp[s: s + dx] = colorbuf
                if isinrange(y1, my, -1):
                    s = _8bmofhd(bmp, x1, y1)
                    bmp[s: s + dx] = colorbuf


def ellipse(bmp: array,
        x: int, y: int, b: int, a: int,
        color: int,
        isfilled: bool = None):
    """Draws an Ellipse

    Args:
        bmp     : unsigned byte array
                  with bmp format
        x, y    : center of ellipse
        b, a    : major amd minor axis
        color   : color of the ellipse
        isfilled: True -> filled
                          ellipse
                  False-> one pixel
                          thick
                          ellipse

    Returns:
        byref modified unsigned byte array
    """
    if isfilled:
        filledellipse(
            bmp, x, y, b, a, color)
    else:
        (mx, my) = getmaxxy(bmp)
        bits = bmp[bmpcolorbits]
        c = _getBMoffhdfunc(bmp)
        dobndcheck = not entireellipseisinboundary(
                            x, y, -1, mx,
                                  -1, my, b, a)
        if bits == 24:
            color = int2BGRarr(color)
            for (px, py) in iterellipse(x, y, b, a):
                if dobndcheck:
                    if isinBMPrectbnd(bmp, px, py):
                        s = c(bmp, px, py)
                        bmp[s: s + 3] = color
                else:
                    s = c(bmp, px, py)
                    bmp[s: s + 3] = color
        elif bits == 8:
            for (px, py) in iterellipse(x, y, b, a):
                if dobndcheck:
                    if isinBMPrectbnd(bmp, px, py):
                        bmp[c(bmp, px, py)] = color
                else:
                    bmp[c(bmp, px, py)] = color
        else:
            for (px, py) in iterellipse(x, y, b, a):
                plotxybit(bmp, px, py, color)


def gradellipse(bmp: array,
        x: int, y: int, b: int, a: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float],
        fn :Callable = lerp):
    """Ellipical gradient

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : center of ellipse
        b, a      : major amd minor axis
        lumrange  : [byte:byte] controls
                    the range of the
                    luminosity gradient
        rgbfactors: [r, g, b] range
                    are from 0.0 to 1.0
        fn        : optional gradent function

    Returns:
        byref modified unsigned byte array
    """
    (l1, l2) = lumrange
    r = max(a, b)
    a -= r
    b -= r
    for i in range(r, 0, -1):
        c = colormix(
                int(fn(l1, l2, i / r)),
                RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                    int2RGBarr(c),
                    getallRGBpal(bmp))
        ellipse(bmp, x, y,
        b + i, a + i, c, True)


@intcircleparam
def drawarc(bmp: array,
        x: int, y: int, r: int,
        startdegangle: float,
        enddegangle: float,
        color: int, showoutline: bool,
        fillcolor: int, isfilled: bool):
    """Draws an Arc

    Args:
        bmp        : unsigned
                     byte array
                     with bmp format
        x, y, r    : defines a circle
                     that contains
                     the arc
        color      : color of arc
        showoutline: True -> draw arc
                             outline
        fillcolor  : color of arc
                     filling
        isfilled   : True -> set area
                             inside the
                             arc to
                             fillcolor

    Returns:
        byref modified unsigned byte array
    """
    av = arcvert(x, y, r,
         startdegangle, enddegangle)
    for p in av:
        plotxybit(bmp, p[0], p[1],
            color)
    if isfilled:
        fillboundary(bmp,
            fillpolydata(av,
                getmaxx(bmp),
                getmaxy(bmp)),
                fillcolor)


def rectangle(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int, color: int):
    """Draws a Rectangle

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: the rectangular
                        region
        color         : color of the
                        rectangle

    Returns:
        byref modified unsigned byte array
    """
    plotpoly(bmp,
        recvert(x1, y1, x2, y2), color)


@entirerectinboundary
def filledrect(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int, color: int):
    """Draws a Filled Rectangle

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region
        color         : color of the
                        rectangle

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    x1, y1, x2, y2 = sortrecpoints(
                        x1, y1, x2, y2)
    if bits not in [8, 24]:
        y2 += 1
        for y in range(y1, y2):
            horiline(bmp, y, x1, x2,
                     color)
    else:
        dx = x2 - x1 + 1
        r = _xchrcnt(bmp)
        rgb = int2RGB(color)
        c = _getBMoffhdfunc(bmp)
        if bits == 24:
            buf = array('B', [rgb[2],
                              rgb[1],
                              rgb[0]] * dx)
        elif bits == 8:
            buf = array('B', [matchRGBtopal(rgb,
                              getallRGBpal(bmp))] * dx)
        offset = c(bmp, x1, y2)
        lim = c(bmp, x2, y1)
        bufsize = len(buf)
        fsize = _flsz(bmp)
        while (offset <= lim) and \
              ((offset + bufsize) <= fsize):
            bmp[offset: offset + bufsize] = buf
            offset += r


def beziercurve(bmp: array,
        pntlist: list[list[Number, Number]],
        penradius: int, color: int):
    """Draws a Bezier Curve

    Args:
        bmp      : unsigned byte array
                   with bmp format
        pntlist  : [(x,y)...]
                   the control points
        penradius: radius of pen
        color    : color of the
                   bezier curve

    Returns:
        byref modified unsigned byte array
    """
    for v in iterbeziercurve(pntlist):
        roundpen(bmp, v, penradius, color)


def bspline(bmp: array, pntlist: list,
        penradius: int, color: int,
        isclosed: bool,
        curveback: bool):
    """Draws a Bspline

    Args:
        bmp      : unsigned byte array
                   with bmp format
        pntlist  : [(x, y)...]
                   the control points
        penradius: radius of pen
        color    : color of curve
        isclosed : True means the
                   curve is closed
        curveback: True means
                   extra computation
                   for curve to loop
                   back on itself

    Returns:
        byref modified unsigned byte array
    """
    for v in iterbspline(pntlist,
             isclosed, curveback):
        roundpen(bmp, v, penradius, color)


def plotrotated8bitpatternwithfn(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int,
        fn: Callable):
    """Draws a 8-bit pattern with
    the bits rotated with a function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : where to draw
                    the pattern
        bitpattern: list of bytes
                    that make a pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of pattern

    Returns:
        byref modified unsigned byte array
    """
    inc = scale - 1 - pixspace
    for bits in bitpattern:
        ox = x
        mask = 128
        bits = rotatebits(bits)
        while mask > 0:
            if (mask & bits) > 0:
                if scale == 1 or inc <= 0:
                    plotxybit(bmp, x, y, color)
                else:
                    fn(bmp, x, y, inc, color)
            mask >>= 1
            x += scale
        y += scale
        x = ox


def plotrotateditalic8bitpatternwithfn(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int,
        fn: Callable):
    """Bit rotated Italic 8-bit pattern
    with a function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the
                    pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of pattern

    Returns:
        byref modified unsigned byte array
    """
    inc = scale - 1 - pixspace
    i = scale >> 1
    h = len(bitpattern)
    if scale > 1:
        x += (i * h)
    else:
        x += (h >> 1)
    if type(color) == int:
        for bits in bitpattern:
            ox = x
            mask = 128
            bits = rotatebits(bits)
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, color)
                    else:
                        fn(bmp, x, y, inc, color)
                mask >>= 1
                x += scale
            y += scale
            x = ox - i
            if y % 2 == 0 and scale == 1:
                x -= 1
    elif type(color) == list or type(color) == tuple:
        colorcnt = len(color)
        colorind = 0
        for bits in bitpattern:
            ox = x
            mask = 128
            bits = rotatebits(bits)
            cl = color[colorind]
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, cl)
                    else:
                        fn(bmp, x, y, inc, cl)
                mask >>= 1
                x += scale
            y += scale
            x = ox - i
            colorind += 1
            if colorind == colorcnt:
                colorind = 0
            if y % 2 == 0 and scale == 1:
                x -= 1


def plot8bitpatternwithfn(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int,
        fn: Callable):
    """8-bit pattern with a function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the
                    pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of pattern

    Returns:
        byref modified unsigned byte array
    """
    inc = scale - 1 - pixspace
    if type(color) == int:
        for bits in bitpattern:
            ox = x
            mask = 128
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, color)
                    else:
                        fn(bmp, x, y, inc, color)
                mask >>= 1
                x += scale
            y += scale
            x = ox
    elif type(color) == list or type(color) == tuple:
        colorcount = len(color)
        colorind = 0
        for bits in bitpattern:
            ox = x
            mask = 128
            cl = color[colorind]
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, cl)
                    else:
                        fn(bmp, x, y, inc, cl)
                mask >>= 1
                x += scale
            y += scale
            colorind += 1
            if colorind == colorcount:
                colorind = 0
            x = ox


def plotitalic8bitpatternwithfn(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int,
        fn: Callable):
    """Italic 8-bit pattern with a function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the
                    pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of pattern

    Returns:
        byref modified unsigned byte array
    """
    inc = scale - 1 - pixspace
    i = scale >> 1
    h = len(bitpattern)
    if scale > 1:
        x += (i * h)
    else:
        x += (h >> 1)
    if type(color) == int:
        for bits in bitpattern:
            ox = x
            mask = 128
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, color)
                    else:
                        fn(bmp, x, y, inc, color)
                mask >>= 1
                x += scale
            y += scale
            x = ox - i
            if y % 2 == 0 and scale == 1:
                x -= 1
    elif type(color) == list or type(color) == tuple:
        colorind = 0
        colorcnt = len(color)
        for bits in bitpattern:
            ox = x
            mask = 128
            cl = color[colorind]
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, cl)
                    else:
                        fn(bmp, x, y, inc, cl)
                mask >>= 1
                x += scale
            y += scale
            x = ox - i
            colorind += 1
            if colorind == colorcnt:
                colorind = 0
            if y % 2 == 0 and scale == 1:
                x -= 1


def plotitalic8bitpatternupsidedownwithfn(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int,
        fn: Callable):
    """Italic 8-bit pattern upsidedown with a function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern
        fn        : function

    Returns:
        byref modified unsigned byte array
    """
    inc = scale - 1 - pixspace
    i = len(bitpattern) - 1
    j = scale >> 1
    if type(color) == int:
        while i > -1:
            bits = bitpattern[i]
            ox = x
            mask = 128
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y,
                            color)
                    else:
                        fn(bmp,x, y, inc,
                            color)
                mask >>= 1
                x += scale
            y += scale
            x = ox - j
            i -= 1
            if y % 2 == 0 and scale == 1:
                x -= 1
    elif type(color) == list or type(color) == tuple:
        colorcnt = len(color)
        colorind = 0
        while i > -1:
            bits = bitpattern[i]
            ox = x
            mask = 128
            cl = color[colorind]
            while mask > 0:
                if (mask & bits)>0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, cl)
                    else:
                        fn(bmp,x, y, inc, cl)
                mask >>= 1
                x += scale
            y += scale
            x = ox - j
            i -= 1
            colorind += 1
            if colorcnt == colorind:
                colorind = 0
            if y % 2 == 0 and scale == 1:
                x -= 1


def plot8bitpattern(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit pattern

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plot8bitpatternwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        lambda bmp, x, y, inc, color: \
            filledrect(bmp, x, y, \
                x + inc, y + inc, color))


def plotitalic8bitpattern(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit italic pattern

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the
                    pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit
        color     : color of pattern

    Returns:
        byref modified unsigned byte array
    """
    plotitalic8bitpatternwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        lambda bmp, x, y, inc, color: \
            filledrect(bmp, x, y, \
                x + inc, y + inc, color))


def plotcircinsqr(bmp, x, y, d, color):
    """Draws a circle in an
        invisible square with side
        equal to the circle's diameter
        and positioned by (x, y)
        at  the upper left
        of the bounding square

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        d         : diameter of the
                    circle
     Returns:
        byref modified unsigned byte array
    """
    d >>= 1
    filledcircle(bmp, x + d, y + d, d,
                                color)


def plot8bitpatternasdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit pattern as circles

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the
                    pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plot8bitpatternwithfn(bmp, x, y,
        bitpattern, scale, pixspace,
        color, plotcircinsqr)


def plotitalic8bitpatternasdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit italic pattern as dots

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plotitalic8bitpatternwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        plotcircinsqr)


def plotitalic8bitpatternupdsidedownasdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit italic pattern upsidedown as dots

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plotitalic8bitpatternupsidedownwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        plotcircinsqr)


def plotreverseditalic8bitpatternasdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit reversed italic pattern as dots

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plotrotateditalic8bitpatternwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        plotcircinsqr)


def plotreverseditalic8bitpattern(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit reversed italic pattern

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plotrotateditalic8bitpatternwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        lambda bmp, x, y, inc, color: \
            filledrect(bmp, x, y, \
                x + inc, y + inc, color))


def plotupsidedownitalic8bitpattern(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit italic pattern upsidedown

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plotitalic8bitpatternupsidedownwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        lambda bmp, x, y, inc, color: \
            filledrect(bmp, x, y, \
                x + inc, y + inc, color))


def plotupsidedownitalic8bitpatternasdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit italic pattern upsidedown as dots

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to draw
                    the pattern
        bitpattern: list of bytes
                    that makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plotitalic8bitpatternupsidedownwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        plotcircinsqr)


def plotrotated8bitpattern(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit pattern with the bits rotated

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : where to draw
                    the pattern
        bitpattern: list of bytes
                    that make a pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the
                    pattern

    Returns:
        byref modified unsigned byte array
    """
    plotrotated8bitpatternwithfn(bmp,
        x, y, bitpattern, scale,
        pixspace, color,
        lambda bmp, x, y, inc, color: \
            filledrect(bmp, x, y, \
                x + inc, y + inc, color))


def plotrotated8bitpatternasdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit pattern with the bits rotated

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : where to draw
                    the pattern
        bitpattern: list of bytes
                    that make a pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the
                    pattern

    Returns:
        byref modified unsigned byte array
    """
    plotrotated8bitpatternwithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        plotcircinsqr)


def plot8bitpatternupsidedownwithfn(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int,
        fn: Callable):
    """Draws a 8-bit pattern upsidedown with a function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern
        fn        : function

    Returns:
        byref modified unsigned byte array
    """
    inc = scale - 1 - pixspace
    i = len(bitpattern) - 1
    if type(color) == int:
        while i > -1:
            bits = bitpattern[i]
            ox = x
            mask = 128
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y,
                            color)
                    else:
                        fn(bmp,x, y, inc,
                            color)
                mask >>= 1
                x += scale
            y += scale
            x = ox
            i -= 1
    elif type(color) == list or type(color) == tuple:
        colorcnt = len(color)
        colorind = 0
        while i > -1:
            bits = bitpattern[i]
            ox = x
            mask = 128
            cl = color[colorind]
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, cl)
                    else:
                        fn(bmp,x, y, inc, cl)
                mask >>= 1
                x += scale
            y += scale
            x = ox
            i -= 1
            colorind += 1
            if colorind == colorcnt:
                colorind = 0


def plot8bitpatternupsidedown(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit pattern upsidedown

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the
                    pattern

    Returns:
        byref modified unsigned byte array
    """
    plot8bitpatternupsidedownwithfn(
        bmp, x, y, bitpattern, scale,
        pixspace, color,
        lambda bmp, x, y, inc, color:
        filledrect(bmp, x, y,
            x + inc, y + inc, color))


def plot8bitpatternupsidedownasdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit pattern upsidedown with dots

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : controls how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plot8bitpatternupsidedownwithfn(
        bmp, x, y, bitpattern, scale,
        pixspace, color, plotcircinsqr)


def plot8bitpatternsidewaywithfn(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int,
        fn: Callable):
    """Draws a 8-bit pattern sideways with a function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    inc = scale - 1 - pixspace
    if type(color) == int:
        for bits in bitpattern:
            oy = y
            mask = 128
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y,
                                      color)
                    else:
                        fn(bmp, x, y, inc,
                                    color)
                mask >>= 1
                y -= scale
            x += scale
            y = oy
    elif type(color) == list or type(color) == tuple:
        colorcnt = len(color)
        colorind = 0
        for bits in bitpattern:
            oy = y
            mask = 128
            cl = color[colorind]
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, cl)
                    else:
                        fn(bmp, x, y, inc, cl)
                mask >>= 1
                y -= scale
                colorind += 1
                if colorcnt == colorind:
                    colorind = 0
            x += scale
            y = oy


def plotitalic8bitpatternsidewaywithfn(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int,
        fn: Callable):
    """Draws an Italic 8-bit pattern sideways with a function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    inc = scale - 1 - pixspace
    i = scale >> 1
    if type(color) == int:
        for bits in bitpattern:
            oy = y
            mask = 128
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y,
                                        color)
                    else:
                        fn(bmp, x, y, inc,
                                        color)
                mask >>= 1
                y -= scale
            x += scale
            y = oy + i
            if x % 2 == 0 and scale == 1:
                y -= 1
    elif type(color) == list or type(color) == tuple:
        colorcnt = len(color)
        colorind = 0
        for bits in bitpattern:
            oy = y
            mask = 128
            cl = color[colorind]
            while mask > 0:
                if (mask & bits) > 0:
                    if scale == 1 or inc <= 0:
                        plotxybit(bmp, x, y, cl)
                    else:
                        fn(bmp, x, y, inc, cl)
                mask >>= 1
                y -= scale
            x += scale
            y = oy + i
            colorind += 1
            if colorind == colorcnt:
                colorind = 0
            if x % 2 == 0 and scale == 1:
                y -= 1

def plot8bitpatternsideway(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit pattern sideways

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plot8bitpatternsidewaywithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        lambda bmp, x, y, inc, color:
        filledrect(bmp, x, y,
            x + inc, y + inc, color))


def plotitalic8bitpatternsideway(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws an italic 8-bit pattern sideways

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plotitalic8bitpatternsidewaywithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        lambda bmp, x, y, inc, color:
        filledrect(bmp, x, y,
            x + inc, y + inc, color))


def plot8bitpatternsidewaywithdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws a 8-bit pattern sideways with dots

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plot8bitpatternsidewaywithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        plotcircinsqr)


def plotitalic8bitpatternsidewayasdots(
        bmp: array, x: int, y: int,
        bitpattern: list, scale: int,
        pixspace: int, color: int):
    """Draws an italic 8-bit pattern sideways as dots

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : sets where to
                    draw the pattern
        bitpattern: list of bytes that
                    makes the pattern
        scale     : control how big
                    the pattern is
        pixspace  : space between
                    each bit in pixels
        color     : color of the pattern

    Returns:
        byref modified unsigned byte array
    """
    plotitalic8bitpatternsidewaywithfn(
        bmp, x, y, bitpattern,
        scale, pixspace, color,
        plotcircinsqr)


def plotstringfunc(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        orderfunc: Callable,
        fontrenderfunc: Callable):
    """Draws a string using a function

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        orderfunc       : function that
                          enumerates
                          each char
                          in the
                          input string
        fontrenderfunc  : function that
                          renders the font

    Returns:
        byref modified unsigned byte array
    """
    if spacebetweenchar == 0:
        spacebetweenchar = 1
    ox = x
    xstep = (scale << 3) + \
                spacebetweenchar
    ystep = fontbuf[0] * scale + \
                spacebetweenchar
    for c in orderfunc(str2plot):
        if c == '\n':
            y += ystep
            x = ox
        elif c == '\t':
            x += xstep << 2
        else:
            fontrenderfunc(bmp, x, y,
                getcharfont(fontbuf, c),
                scale, pixspace, color)
            x += xstep


def plotstringfunc2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        orderfunc: Callable,
        fontrenderfunc: Callable,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a string using a function to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        orderfunc       : function that
                          enumerates
                          each char
                          in the
                          input string
        fontrenderfunc  : function that
                          renders the font
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    if spacebetweenchar == 0:
        spacebetweenchar = 1
    x = y = scale
    ox = x
    xstep = (scale << 3) + \
                spacebetweenchar
    h = fontbuf[0]
    ystep = h * scale + \
                spacebetweenchar
    (row, col) = computetextrowandcolcnt(str2plot)
    pad = (scale << 1)
    mx = row * xstep + pad
    my = col * ystep + pad
    if fontrenderfunc.__name__.find("italic") > 0:
        if scale > 1:
            mx += (h * (scale >> 1))
        else:
            mx += (h >> 1)
    bmp = newBMP(mx, my, bitdepth)
    if backgroundcolor > 0:
        mx -= 1
        my -= 1
        filledrect(bmp, 0, 0, mx, my, backgroundcolor)
    for c in orderfunc(str2plot):
        if c == '\n':
            y += ystep
            x = ox
        elif c == '\t':
            x += xstep << 2
        else:
            fontrenderfunc(bmp, x, y,
                getcharfont(fontbuf, c),
                scale, pixspace, color)
            x += xstep
    saveBMP(file, bmp)


def plotstringverticalwithfn2file(
        file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        orderfunc: Callable,
        fontrenderfunc: Callable,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a string vertically using a function
    to a new file


    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        orderfunc       : function that
                          enumerates
                          each char
                          in the
                          input string
        fontrenderfunc  : function that
                          renders the font
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    if spacebetweenchar == 0:
        spacebetweenchar = 1
    h = fontbuf[0]
    xstep = (scale << 3) + \
                spacebetweenchar
    ystep = h * scale + \
                spacebetweenchar
    x = y = scale
    oy = y
    (row, col) = computetextrowandcolcnt(str2plot)
    pad = scale << 1
    mx = col * xstep + pad
    my = row * ystep + pad
    if fontrenderfunc.__name__.find("italic") > 0:
        if scale > 1:
            mx += (h * (scale >> 1))
        else:
            mx += (h >> 1)
    bmp = newBMP(mx, my, bitdepth)
    if backgroundcolor > 0:
        mx -= 1
        my -= 1
        filledrect(bmp, 0, 0, mx, my, backgroundcolor)
    for c in orderfunc(str2plot):
        if c == '\n':
            x += xstep
            y = oy
        elif c == '\t':
            y += ystep << 2
        else:
            fontrenderfunc(bmp, x, y,
               getcharfont(fontbuf,c),
               scale, pixspace, color)
            y += ystep
    saveBMP(file, bmp)



def plotstring(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the
                          string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
                          (can be an int
                               or a
                          list of ints)
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumletters,
        plot8bitpattern)


def plotstring2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a string to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumletters,
        plot8bitpattern,
        backgroundcolor,
        bitdepth)


def plotverticalstring2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a vertical string to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringverticalwithfn2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumletters,
        plot8bitpattern,
        backgroundcolor,
        bitdepth)


def plotitalicstring(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string as Italic

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the
                          string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumletters,
        plotitalic8bitpattern)


def plotitalicstring2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws an italic string to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumletters,
        plotitalic8bitpattern,
        backgroundcolor,
        bitdepth)


def plotverticalitalicstring2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a vertical italic string to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringverticalwithfn2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumletters,
        plotitalic8bitpattern,
        backgroundcolor,
        bitdepth)


def plotreverseditalicstring2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a reversed italic string to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumreverseletters,
        plotreverseditalic8bitpattern,
        backgroundcolor,
        bitdepth)


def plotverticalstringasdots2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a vertical string as Dots to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringverticalwithfn2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumletters,
        plot8bitpatternasdots,
        backgroundcolor,
        bitdepth)


def plotverticalitalicstringasdots2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a vertical italic string as Dots to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringverticalwithfn2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumletters,
        plotitalic8bitpatternasdots,
        backgroundcolor,
        bitdepth)


def plotstringasdots(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string as Dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumletters,
        plot8bitpatternasdots)


def plotstringasdots2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a string as Dots to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumletters,
        plot8bitpatternasdots,
        backgroundcolor,
        bitdepth)


def plotitalicstringasdots2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws an italic string as Dots to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumletters,
        plotitalic8bitpatternasdots,
        backgroundcolor,
        bitdepth)


def plotreversedstring2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a reversed string to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumreverseletters,
        plotrotated8bitpattern,
        backgroundcolor,
        bitdepth)


def plotreversedstringasdots2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws a reversed string as Dots to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumreverseletters,
        plotrotated8bitpatternasdots,
        backgroundcolor,
        bitdepth)


def plotitalicstringasdots(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string as Italic dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the
                          string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumletters,
        plotitalic8bitpatternasdots)


def plotreverseditalicstringasdots(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a Reversed String as Italic dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the
                          string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumreverseletters,
        plotreverseditalic8bitpatternasdots)


def plotreverseditalicstring(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a reversed string as Italic

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumreverseletters,
        plotreverseditalic8bitpattern)


def plotstringupsidedown(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string upsidedown

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar, color,
        fontbuf, enumreverseletters,
        plot8bitpatternupsidedown)


def plotupsidedownstring2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws an upsidedown string to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumreverseletters,
        plot8bitpatternupsidedown,
        backgroundcolor,
        bitdepth)


def plotupsidedownstringasdots2file(file: str,
        str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        backgroundcolor: int = 0,
        bitdepth: int = 24):
    """Draws an upsidedown string as Dots to a file

    Args:
        file            : new file
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)
        backgroundcolor : optional background color
                          (default = 0)
        bitdepth        : optional bitdepth
                          of bitmap default = 24
                          (1, 4, 8, 24) bits

    Returns:
        new bitmap file
    """
    plotstringfunc2file(
        file, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf,
        enumreverseletters,
        plot8bitpatternupsidedownasdots,
        backgroundcolor,
        bitdepth)


def plotupsidedownitalicstring(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draw an italic string upsidedown

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the
                          string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumreverseletters,
        plotupsidedownitalic8bitpattern)


def plotupsidedownitalicstringasdots(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draw an italic string upsidedown as dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the
                          string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumreverseletters,
        plotupsidedownitalic8bitpatternasdots)


def plotstringupsidedownasdots(
        bmp: array, x: int, y: int,
        str2plot: str, scale: int,
        pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string upsidedown as dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar, color,
        fontbuf, enumletters,
        plot8bitpatternupsidedownasdots)


def plotreversestring(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string reversed

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar,
        color, fontbuf, enumreverseletters,
        plotrotated8bitpattern)


def plotreversestringasdots(
        bmp: array, x: int, y: int,
        str2plot: str, scale: int,
        pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string reversed with dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringfunc(bmp, x, y, str2plot,
        scale, pixspace,
        spacebetweenchar, color,
        fontbuf, enumreverseletters,
        plotrotated8bitpatternasdots)


def plotstringsidewayfn(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        fn: Callable):
    """Draws a string sideways with a function

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    if spacebetweenchar == 0:
        spacebetweenchar = 1
    oy = y
    xstep = (scale << 3 ) + spacebetweenchar
    ystep = fontbuf[0] * scale + spacebetweenchar
    for c in enumletters(str2plot):
        if c == '\n':
            x += ystep #we swap x and y since sideways
            y = oy
        elif c == '\t':
            y -= xstep << 2 #we swap x and y since sideways
        else:
            fn(bmp, x, y,
               getcharfont(fontbuf, c),
               scale, pixspace, color)
            y -= xstep


def plotstringsideway(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string sideways

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringsidewayfn(bmp, x, y,
        str2plot, scale, pixspace,
        spacebetweenchar, color, fontbuf,
        plot8bitpatternsideway)


def plotitalicstringsidewayasdots(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws an italic string sideways as dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringsidewayfn(bmp, x, y,
        str2plot, scale, pixspace,
        spacebetweenchar, color, fontbuf,
        plotitalic8bitpatternsidewayasdots)


def plotitalicstringsideway(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws an Italic String Sideways

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
                          can be an int
                              or a
                          list of ints
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringsidewayfn(bmp, x, y,
        str2plot, scale, pixspace,
        spacebetweenchar, color, fontbuf,
        plotitalic8bitpatternsideway)


def plotstringsidewayasdots(
        bmp: array, x: int, y: int,
        str2plot: str, scale: int,
        pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string sideways as dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringsidewayfn(bmp, x, y,
        str2plot, scale, pixspace,
        spacebetweenchar, color, fontbuf,
        plot8bitpatternsidewaywithdots)


def plotstringverticalwithfn(
        bmp: array, x: int, y: int,
        str2plot: str, scale: int,
        pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        fn: Callable):
    """Draws a string vertically using a function

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                          (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    if spacebetweenchar == 0:
        spacebetweenchar = 1
    oy = y
    xstep = (scale << 3) + \
                spacebetweenchar
    ystep = fontbuf[0] * scale + \
                spacebetweenchar
    for c in enumletters(str2plot):
        if c == '\n':
            x += xstep
            y = oy
        elif c == '\t':
            y += ystep << 2
        else:
            fn(bmp, x, y,
               getcharfont(fontbuf,c),
               scale, pixspace, color)
            y += ystep


def plotstringvertical(bmp: array,
        x: int, y: int, str2plot: str,
        scale: int, pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string vertically

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringverticalwithfn(bmp,
        x, y, str2plot, scale,
        pixspace, spacebetweenchar,
        color, fontbuf, plot8bitpattern)


def plotstringverticalasdots(
        bmp: array, x: int, y: int,
        str2plot: str, scale: int,
        pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws a string vertically with dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringverticalwithfn(bmp,
        x, y, str2plot, scale,
        pixspace, spacebetweenchar,
        color, fontbuf,
        plot8bitpatternasdots)


def plotitalicstringverticalasdots(
        bmp: array, x: int, y: int,
        str2plot: str, scale: int,
        pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws an italic string vertically with dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringverticalwithfn(bmp,
        x, y, str2plot, scale,
        pixspace, spacebetweenchar,
        color, fontbuf,
        plotitalic8bitpatternasdots)


def plotitalicstringvertical(
        bmp: array, x: int, y: int,
        str2plot: str, scale: int,
        pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list):
    """Draws an italic string vertically with dots

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y            : sets where to
                          draw the string
        str2plot        : string to draw
        scale           : control how big
                          the font is
        pixspace        : space between
                          each bit
        spacebetweenchar: space between
                          the characters
        color           : color of
                          the font
        fontbuf         : the font
                         (see fonts.py)

    Returns:
        byref modified unsigned byte array
    """
    plotstringverticalwithfn(bmp,
        x, y, str2plot, scale,
        pixspace, spacebetweenchar,
        color, fontbuf,
        plotitalic8bitpattern)


def fillboundary(bmp: array,
        bndfilldic: dict, color: int):
    """Draws lines in a boundary to fill it

    Args:
        bmp       : unsigned byte array
                    with bmp format
        bndfilldic: boundary dictionary
        color     : color of fill

    Returns:
        byref modified unsigned byte array
    """
    for y in bndfilldic:
        b =  bndfilldic[y]
        yint = len(b)
        if yint == 1:
            plotxybit(bmp, b[0], y, color)
        else:
            for j in range(1, yint):
                horiline(bmp, y, b[j - 1], b[j], color)


def plotpolyfill(bmp: array,
        vertlist: list[list[Number, Number]],
        color: int):
    """Draws a filled polygon with a given color

    Args:
        bmp     : unsigned byte array
                  with bmp format
        vertlist: [(x, y), ...]
                  list of vertices
        color   : color of the
                  filled polygon

    Returns:
        byref modified unsigned byte array
    """
    fillboundary(bmp,
        fillpolydata(polyboundary(vertlist),
            getmaxx(bmp), getmaxy(bmp)), color)


def thickplotpoly(bmp: array,
        vertlist: list[list[Number, Number]],
        penradius: int, color: int):
    """Draws a polygon of a given color and thickness

    Args:
        bmp      : unsigned byte array
                   with bmp format
        vertlist : [(x, y)...]
                   list of vertices
        penradius: radius of pen
        color    : color of the polygon

    Returns:
        byref modified unsigned byte array
    """
    vertcount = len(vertlist)
    for i in range(vertcount):
        if i > 0:
            thickroundline(bmp,
                roundvect(vertlist[i - 1]),
                roundvect(vertlist[i]),
                penradius, color)
    thickroundline(bmp,
        roundvect(vertlist[0]),
        roundvect(vertlist[vertcount - 1]),
        penradius, color)


def gradthickplotpoly(bmp: array,
        vertlist: list[list[Number, Number]],
        penradius: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float],
        fn: Callable = lerp):
    """Draws a polygon of a given gradient and thickness

    Args:
        bmp       : unsigned byte array
                    with bmp format
        vertlist  : [(x,y)...] the
                    list of vertices
        penradius : radius of pen
        lumrange  : [byte,byte] range
                    of the gradient
        RGBfactors: [r, g, b] value
                    range from
                    0.0 to 1.0
        fn        : optional gradient function

    Returns:
        byref modified unsigned byte array
    """
    (l1, l2) = lumrange
    for i in range(penradius, 0, -1):
        c = colormix(int(
                fn(l1, l2, i / penradius)),
                RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                int2RGBarr(c),
                getallRGBpal(bmp))
        thickplotpoly(bmp, vertlist, i, c)


def gradplotlines(bmp: array,
        vertlist: list[list[Number, Number]] ,
        penradius: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float],
        fn : Callable = lerp):
    """Draws connected lines of a given gradient and thickness

    Args:
        bmp       : unsigned byte array
                    with bmp format
        vertlist  : [(x,y)...] the
                    list of vertices
        penradius : radius of pen
        lumrange  : [byte,byte] range
                    of the gradient
        RGBfactors: [r, g, b] value
                    range from
                    0.0 to 1.0
        fn        : optional gradient function

    Returns:
        byref modified unsigned byte array
    """
    (l1, l2) = lumrange
    for i in range(penradius, 0, -1):
        c = colormix(int(
                fn(l1, l2, i / penradius)),
                RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                int2RGBarr(c),
                getallRGBpal(bmp))
        plotlines(bmp, vertlist, c, i)



def plotlines(bmp: array,
        vertlist: list[list[Number, Number]],
        color: int,
        penradius: int = 1
        ):
    """Draws connected lines defined by a list of vertices

    Args:
        bmp      : unsigned byte array
                   with bmp format
        vertlist : [(x:uint,y:uint),...]
                   list of vertices
        color    : color of the lines
        penradius: optional parameter
                   for thick line


    Returns:
        byref modified unsigned byte array
    """
    vertcount = len(vertlist)
    if penradius <= 1:
        for i in range(vertcount):
            if i > 0:
                linevec(bmp,
                    vertlist[i - 1],
                    vertlist[i], color)
    elif penradius >= 2:
        for i in range(vertcount):
            if i > 0:
                thickroundline(bmp,
                    roundvect(vertlist[i - 1]),
                    roundvect(vertlist[i]),
                    penradius, color)


def plotpoly(bmp: array,
        vertlist: list[list[Number, Number]],
        color: int):
    """Draws a polygon defined by a list of vertices

    Args:
        bmp     : unsigned byte array
                  with bmp format
        vertlist: [(x:uint,y:uint),...]
                  list of vertices
        color   : color of the lines

    Returns:
        byref modified unsigned byte array
    """
    plotlines(bmp, vertlist, color)
    linevec(bmp, vertlist[0], vertlist[-1], color)


def plotpolylist(bmp: array,
        polylist: list[list[Number, Number]],
        color: int):
    """Draws a list of polygons of a given color

    Args:
        bmp      : unsigned byte array
                   with bmp format
        polytlist: [[(x:uint,y:uint),
                    ...],...]
                   list of polygons
        color    : color of the lines

    Returns:
        byref modified unsigned byte array
    """
    for poly in polylist:
        plotpoly(bmp, poly, color)


def plotpolyfillist(
        bmp: array,
        sides: list[list[list[list]],
               list[list[float, float, float]]],
        RGBfactors: list[float, float]):
    """3D polygon rendering function

    Args:
        bmp       : unsigned byte array
                    with bmp format
        sides     : list of polygons
                    and normals
        RGBfactors: [r, g, b]
                    r, g, b are
                    float values
                    from 0.0 to 1.0

    Returns:
        byref modified unsigned byte array
    """
    (polylist, normlist) = sides
    for i, poly in enumerate(polylist):
        c = (
            colormix(
                int(cosaffin(normlist[i], [0, 0, 1]) * 128) + 127, RGBfactors
            )
            if normlist[i] != [0, 0, 0]
            else colormix(127, RGBfactors)
        )

        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                    int2RGBarr(c),
                    getallRGBpal(bmp))
        plotpolyfill(bmp, poly, c)


def plot3d(bmp: array,
        sides: list[list, list],
        issolid: bool,
        RGBfactors: list[float, float],
        showoutline: bool,
        outlinecolor: int):
    """The 3D rendering function

    Args:
        bmp         : unsigned
                      byte array
                      with bmp format
        sides       : list of polygons
                      and normals
        isolid      : toggles solid render
        RGBfactors  : [r,g,b] r, g, b
                      range in value
                      from 0.0 to 1.0
        showoutine  : toggles the
                      polygon outline
        outlinecolor: color of the
                      polygon outline

    Returns:
        byref modified unsigned byte array
    """
    if issolid:
        plotpolyfillist(
            bmp, sides, RGBfactors)
    if showoutline:
        plotpolylist(
            bmp, sides[0], outlinecolor)


def plot3Dsolid(bmp: array,
        vertandsides: list[list, list],
        issolid: bool,
        RGBfactors: list[float, float, float],
        showoutline: bool,
        outlinecolor: int,
        rotvect: list[float, float, float],
        transvect3D: list[float, float, float],
        d: int,
        transvect: list[int, int]):
    """3D solid rendering function

    Args:
        bmp         : unsigned
                      byte array
                      with bmp format
        sides       : list of polygons
                      and normals
        isolid      : toggles the
                      solid render
        RGBfactors  : [r,g,b] r, g, b
                      range in value
                      from 0.0 to 1.0
        showoutine  : toggles the
                      polygon outline
        outlinecolor: color of the
                      polygon outline
        rotvect     : rotation vector
        transvect3D : 3D translation
                      vector
        d           : distance of the
                      observer from
                      the screen
        transvect   : 2D translation
                      vector for
                      screen position

    Returns:
        byref modified unsigned byte array
    """
    plot3d(bmp, gensides(perspective(
           vertandsides[0], rotvect,
           transvect3D, d),
        transvect, vertandsides[1]),
        issolid, RGBfactors,
        showoutline,
        outlinecolor)


def gradvert(bmp: array,
        vertlist: list[list[int, int]],
        penradius: int,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float],
        fn: Callable = lerp):
    """List of 2d vertices as spheres of a given color

    Args:
        bmp       : unsigned byte array
                    with bmp format
        vertlist  : [(x, y), ...]
                    list of vertices
        penradius : radius of
                    the spheres
        lumrange  : [byte, byte] sets
                    the luminosity
                    gradient
        RGBfactors: (r, g, b)
                    values range from
                    min 0.0 to 1.0 max
        fn        : optional gradient function

    Returns:
        byref modified
        unsigned byte array
    """
    (l1, l2) = lumrange
    for i in range(penradius, 0, -1):
        c = colormix(
            fn(l1, l2, i / penradius),
                RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = matchRGBtopal(
                    int2RGBarr(c),
                    getallRGBpal(bmp))
        for point in vertlist:
            roundpen(bmp, point, i, c)


@entirerectinboundary
def xygrid(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        xysteps: list[int, int],
        color: int):
    """Draws a grid

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: sets limits
                        of the grid
        xysteps       : [x, y] sets the
                        increments
        color         : sets the color
                        of the grid

    Returns:
        byref modified unsigned byte array
    """
    x1, y1, x2, y2 = sortrecpoints(
                        x1, y1, x2, y2)
    xstep, ystep = xysteps
    for x in range(x1, x2, xstep):
        vertline(bmp, x, y1, y2, color)
    for y in range(y1, y2, ystep):
        horiline(bmp, y, x1, x2, color)


def xygridvec(bmp: array,
        u: list[int, int],
        v: list[int, int],
        steps: list[int, int],
        gridcolor: int):
    """Grid using (x, y) point pairs u and v

    Args:
        bmp  : unsigned byte array
               with bmp format
        u, v : (x, y) sets limits
                of the  grid
        steps: (x, y) -> sets the
               increments for x and y
        color: sets the color
               of the grid

    Returns:
        byref modified unsigned byte array
    """
    xygrid(bmp, u[0], u[1],
                v[0], v[1],
                steps, gridcolor)


def numbervert(bmp: array,
        vlist: list[list[int, int]],
        xadj: int, yadj: int,
        scale: int,
        valstart: Number,
        valstep: Number,
        pixspace: int,
        spacebetweenchar: int,
        color: int, fontbuf: list,
        suppresszero: bool,
        suppresslastnum: bool,
        rightjustify: bool):
    plot = False
    maxv = len(vlist) - 1
    for v in vlist:
        i = vlist.index(v)
        plot = True if i > 0 else not suppresszero
        if i == maxv and suppresslastnum:
            plot = False
        stval = str(valstart + i * valstep)
        rjust = 0
        if rightjustify:
            rjust = (len(stval) - 1) << 3
        if plot:
            plotstring(bmp,
                v[0] + xadj - rjust,
                v[1] + yadj,
                stval, scale, pixspace,
                spacebetweenchar,
                color, fontbuf)


def vertlinevert(bmp: array,
        vlist: list[list[int, int]],
        linelen: int,
        yadj: int, color: int):
    """Vertical line marks at vertices in vlist

    Args:
        bmp    : unsigned byte array
                 with bmp format
        vlist  : [(x,y),...] the
                 list of vertices
        linelen: lenght of the
                 vertical lines
        yadj   : sets an adjustment
                 for y coordinates
        color  : color of the line

    Returns:
        byref modified unsigned byte array
    """
    for (x, y) in vlist:
        vertline(bmp, x, y, y + linelen + yadj, color)


def horilinevert(bmp: array,
        vlist: list[list[int, int]],
        linelen: int,
        xadj: int, color: int):
    """Horizontal line marks at vertices in vlist

    Args:
        bmp    : unsigned byte array
                 with bmp format
        vlist  : [(x, y), ...]
                 the list of
                 2D vertices
        linelen: length of the
                 vertical lines
        xadj   : sets an adjustment
                 for x coordinates
        color  : color of the line

    Returns:
        byref modified unsigned byte array
    """
    for (x, y) in vlist:
        horiline(bmp, y, x, x + linelen + xadj, color)


def XYaxis(bmp: array,
        origin: list[int, int],
        steps: list[int, int],
        xylimits: list[int, int],
        xyvalstarts: list[Number, Number],
        xysteps: list[Number, Number],
        color: int, textcolor: int,
        showgrid: bool,
        gridcolor: int):
    """XY axis with tick marks and numbers

    Args:
        bmp      : unsigned byte array
                   with bmp format
        origin   : (x, y) on screen
                    origin point
                    of the axis
        steps    : (x, y) steps between
                   tick marks onscreen
        xylimits : (x, y) sets where the
                   graph ends onscreen
        xyvalstarts: (x, y) sets the
                      start point of
                      x and y number
                      lines
        xysteps   : (x, y) sets the
                    number increment
                    along the x and y
                    numberlines
        color    : color of the lines
        textcolor: color of the
                   numberline text
        showgrid : True -> display
                           gridline
                   False -> no grid
        gridcolor: color of the grid

    Returns:
        byref modified
        unsigned byte array
    """
    hvert = horizontalvert(origin[1],
                origin[0], xylimits[0],
                steps[0])
    vvert = verticalvert(origin[0],
                origin[1], xylimits[1],
                -steps[1])
    if showgrid:
        xygridvec(bmp, origin, xylimits,
            steps, gridcolor)
    drawvec(bmp, origin, [origin[0],
            xylimits[1]], 10, color)
    drawvec(bmp, origin,
      [xylimits[0], origin[1]], 10, color)
    vertlinevert(bmp, hvert, 5, -2, color)
    horilinevert(bmp, vvert, -3, 0, color)
    font = font8x14
    numbervert(bmp, vvert, -15, -4, 1,
        xyvalstarts[1],
        xysteps[1],
        0, 0, textcolor, font,
        False, False, True)
    numbervert(bmp, hvert, -4, 7, 1,
        xyvalstarts[0],
        xysteps[0],
        0, 0, textcolor, font,
        False, False, False)
    xvalmax = xyvalstarts[0] + \
            (len(hvert) - 1) * xysteps[0]
    yvalmax = xyvalstarts[1] + \
            (len(vvert) - 1) * xysteps[1]
    return (origin, steps, xylimits,
            xyvalstarts, xysteps,
            (xvalmax, yvalmax))


def userdef2Dcooordsys2screenxy(
        x: int, y: int,
        lstcooordinfo: list):
    """2D coordinate trans from user to screen

    Args:
        x, y         : user coordinates
        lstcooordinfo: info on how to
                       transform the
                       2D coordinate
                       system
                       [origin,
                        steps,
                        xylimits,
                        xyvalstarts,
                        xysteps]
                        all
                        (x: int,
                         y: int) pairs

    Returns:
        [x: int, y: int] screen coordinates

    """
    (xo, yo), (xs, ys), (xlim, ylim), (xvst, yvst), (xvstp, yvstp) = lstcooordinfo[0:5]
    x = xo + ((x - xvst) / xvstp) * xs
    y = yo - ((y - yvst) / yvstp) * ys
    if x > xlim or y < ylim:
        x = -1
        y = -1
        print(sysmsg['regionoutofbounds'])
    return [x, y]


def XYscatterplot(
        bmp: array,
        XYdata: list,
        XYcoordinfo: list,
        showLinearRegLine: bool,
        reglinecolor: int):
    """Create a XY scatterplot

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        XYData          : [[x,y,
                          radius
                          (max radius is 5),
                          isfilled],...]
        lstcooordinfo   : info on how to
                          transform the
                          coordinate system
                          [origin,
                           -> origin point
                              on screen
                          steps,
                           -> on screen steps
                          xylimits,
                           -> x and y
                              clipping limit
                          xyvalstarts,
                           -> number line
                              start values
                          xysteps
                           -> increments for
                              the number lines
                            ]
                          * (x:int,y:int) pairs
        howLinearRegLine: True -> display linear
                                  regression line
        reglinecolor    : color of linear
                          regression line

    Returns:
        byref modified unsigned byte array
    """
    for v in XYdata:
        r = v[2]
        w = userdef2Dcooordsys2screenxy(
                v[0], v[1], XYcoordinfo)
        if r>1:
            circlevec(bmp, w,
                setmax(r, 5), v[3], v[4])
        else:
            plotvecxypoint(bmp, w, v[3])
    if showLinearRegLine:
        m = LSMslope(XYdata)
        b = LSMYint(XYdata)
        xmin = XYcoordinfo[3][0]
        xmax = XYcoordinfo[5][0]
        ymin = xmin * m + b
        ymax = xmax * m + b
        u = userdef2Dcooordsys2screenxy(
                xmin, ymin, XYcoordinfo)
        w = userdef2Dcooordsys2screenxy(
                xmax, ymax, XYcoordinfo)
        if w == [-1, -1]:
            ymax = XYcoordinfo[5][1]
            xmax = (ymax - b) / m
            w = userdef2Dcooordsys2screenxy(
                    xmax, ymax, XYcoordinfo)
        linevec(bmp, u, w, reglinecolor)


def getneighborcolorlist(
        bmp: array, v: list):
    return [getRGBxybitvec(bmp, u)
            for u in itergetneighbors(
                v, getmaxx(bmp), getmaxy(bmp), True)]


def iterimagedgevert(
        bmp: array,
        similaritythreshold: float):
    """Find edges in an image

    Args:
        bmp                : unsigned
                             byte array
                             with
                             bmp format
        similaritythreshold: how close
                             to the color
                             before we
                             yield it

    Yields:
        (x: int, y: int)
    """
    (mx, my) = getmaxxy(bmp)
    for (v0, v1) in iterimageRGB(bmp,
                sysmsg['edgedetect'],
                '*' , sysmsg['done']):
        for u in itergetneighbors(
                    v0, mx, my, False):
            if dist(
                    getRGBxybitvec(bmp, u), v1) > \
                        similaritythreshold:
                yield u
                break


def iterimageregionvertbyRGB(
        bmp: array,
        rgb: list[int, int, int],
        similaritythreshold: int):
    """RGB Color selection by color similarity

    Args:
        bmp                : unsigned
                             byte array
                             with bmp
                             format
        rgb                : (r: byte,
                              g: byte,
                              b: byte)
        similaritythreshold: how close
                             to the color
                             before we
                             yield it

    Yields:
        ((x: int, y: int), (r: byte, g: byte, b: byte))
    """
    for v in iterimageRGB(bmp,
                sysmsg['edgedetect'],
                '*', sysmsg['done']):
        if dist(rgb, v[1]) < similaritythreshold:
            yield v[0]


def getimageregionbyRGB(
        bmp: array,
        rgb: list[int, int, int],
        similaritythreshold: int):
    """Select a region by color

    Args:
        bmp                 :unsigned
                             byte array
                             with bmp
                             format
        rgb                 :(r: byte,
                              g: byte,
                              b: byte)
        similaritythreshold: controls
                             the edge
                             detection
                             sensitivity

    Returns:
        list of vertices
    """
    return list(iterimageregionvertbyRGB(bmp, rgb, similaritythreshold))


def getimagedgevert(bmp: array,
        similaritythreshold: int):
    """Find edges in an image

    Args:
        bmp                : unsigned
                             byte array
                             with bmp
                             format
        similaritythreshold: how close
                             to the color
                             before we
                             yield it

    Yields:
        [(x: int, y: int),...]
    """
    return list(iterimagedgevert(bmp, similaritythreshold))


def plotimgedges(bmp: array,
        similaritythreshold: int,
        edgeradius: int,
        edgecolor: int):
    """Draw edges

    Args:
        bmp                : unsigned
                             byte array
                             with bmp
                             format
        similaritythreshold: controls
                             the edge
                             detection
                             sensitivity
        edgeradius         : radius and
        edgecolor            color of
                             the pen
                             used to
                             draw the
                             edges

    Returns:
        byref modified unsigned byte array
    """
    plotxypointlist(bmp,
        getimagedgevert(bmp,
            similaritythreshold),
        edgeradius, edgecolor)


def getBGRpalbuf(bmp: array):
    """Gets bitmap palette as stored in the bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        unsigned byte array (BGR)
    """
    return bmp[bmppal: _hdsz(bmp)]


def convertbufto24bit(buf: array,
        bgrpalbuf: array,
        bits: int) -> array:
    """Converts 1, 4 and 8-bit buffers to a BGR buffer

    Args:
        buf      : unsigned byte array
        bgrpalbuf: BGR palette info
        bits     : color depth
                   (1, 4, 8)

    Returns:
        unsigned byte array
    """
    retval = []
    for b in buf:
        if bits == 8:
            s = b << 2
            retval += bgrpalbuf[s: s + 3]
        elif bits == 4:
            s1, s2 = divmod(b, 16)
            s1 <<= 2
            s2 <<= 2
            retval += bgrpalbuf[s1: s1 + 3] + \
                      bgrpalbuf[s2: s2 + 3]
        elif bits == 1:
            for i in enumbits(b):
                s = i << 2
                retval += bgrpalbuf[s: s + 3]
    return array('B', retval)


def upgradeto24bitimage(bmp: array):
    """Upgrade an image to 24-bits

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    if bits == 24:
        print(sysmsg['is24bit'])
        nbmp = bmp
    else:
        nbmp = CopyBMPxydim2newBMP(
                           bmp, 24)
        bgrpal = getBGRpalbuf(bmp)
        mx = getmaxx(bmp) - 1
        my = getmaxy(bmp) - 1
        offset = _BMoffset(nbmp, 0, my)
        r = _xchrcnt(nbmp)
        for buf in itercopyrect(bmp,
                        0, 0, mx, my):
            BMPbitBLTput(nbmp, offset,
                convertbufto24bit(buf,
                    bgrpal, bits))
            offset += r
    return nbmp


def iterimageRGB(bmp: array,
        waitmsg: str, rowprocind: str,
        finishmsg: str):
    """Yields (r, g, b) information for the entire bitmap

    Args:
        bmp       : unsigned byte array
                    with bmp format
        waitmsg   : what to display
                    in terminal at
                    process start
        rowprocind: char to display as
                    a row is processed
        finishmsg : what to display
                    in terminal at
                    process end

    Yields:
        ((x: int, y: int), (r: byte, g: byte, b: byte))
    """
    if waitmsg != '':
        print(waitmsg)
    y = getmaxy(bmp) - 1
    offset = 0
    b = getBMPimgbytes(bmp)
    maxoffset = len(b)
    x = 0
    mx = getmaxx(bmp)
    bits = bmp[bmpcolorbits]
    padbytes = _pdbytes(mx, bits)
    if bits < 24:
        p = getallRGBpal(bmp)
        doff = 1
    else:
        doff = 3
    while offset < maxoffset:
        if bits == 1:
            c = b[offset]
            for bit in range(7, -1, -1):
                if x < mx:
                    yield ((x, y), p[(c & 1 << bit) >> bit])
                x += 1
        elif bits == 4:
            c0, c1 = divmod(b[offset], 16)
            if x < mx:
                yield ((x, y), p[c0])
            x += 1
            if x < mx:
                yield ((x, y), p[c1])
            x += 1
        elif bits == 8:
            if x < mx:
                yield ((x, y), p[b[offset]])
            x += 1
        elif bits == 24:
            if x < mx:
                yield ((x, y), (b[offset + 2],
                                b[offset + 1],
                                b[offset]))
            x += 1
        if x == mx:
            x = 0
            y -= 1
            offset += padbytes
            if rowprocind != '':
                print(rowprocind, end = '')
        offset += doff
    print('\n')
    if finishmsg != '':
        print(finishmsg)


def iterimagecolor(bmp: array,
        waitmsg: str, rowprocind: str,
        finishmsg: str):
    """Yields color information for entire bitmap

    Args:
        bmp       : unsigned byte array
                    with bmp format
        waitmsg   : what to display
                    in the terminal
                    when process starts
        rowprocind: string to print in
                    the terminal as a
                    row is processed as
                    a process indicator
        finishmsg : what to display
                    in the terminal
                    when process ends

    Yields:
        ((x: int, y: int), color: int)
    """
    if waitmsg != '':
        print(waitmsg)
    y = getmaxy(bmp) - 1
    offset = 0
    b = getBMPimgbytes(bmp)
    maxoffset = len(b)
    x = 0
    mx = getmaxx(bmp)
    bits = bmp[bmpcolorbits]
    pb = _pdbytes(mx, bits)
    doff = 3 if bits == 24 else 1
    while offset < maxoffset:
        if bits == 1:
            c = b[offset]
            for bit in range(7, -1, -1):
                if x < mx:
                    yield ((x, y), (c & 1 << bit) >> bit)
                x += 1
        elif bits == 4:
            c0, c1 = divmod(b[offset], 16)
            if x < mx:
                yield ((x, y), c0)
            x += 1
            if x < mx:
                yield ((x, y), c1)
            x += 1
        elif bits == 8:
            if x < mx:
                yield ((x, y), b[offset])
            x += 1
        elif bits == 24:
            if x < mx:
                yield ((x, y), (b[offset + 2] << 16) + (b[offset + 1] << 8) + b[offset])
            x += 1
        if x >= mx:
            x = 0
            y -= 1
            offset += pb
            if rowprocind != '':
                print(rowprocind, end='')
        offset += doff
    print('\n')
    if finishmsg != '':
        print(finishmsg)


@entirerectinboundary
def copyrect(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int) -> array:
    """Copy a rectangular region to a custom buffer

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        custom unsigned byte array
    """
    retval = array('B', [bmp[bmpcolorbits]])
    x1, y1, x2, y2 = \
       sortrecpoints(x1, y1, x2, y2)
    retval += int2buf(2, x2 - x1 + 2)
    retval += int2buf(2, y2 - y1 + 1)
    retval += int2buf(2,
                adjustxbufsize(bmp, x1, x2))
    for buf in itercopyrect(
                    bmp, x1, y1, x2, y2):
        retval += buf
    return retval


def pasterect(bmp: array, buf: array,
                   x1: int, y1: int):
    """Paste a rectangular area from a buffer to a bmp

    Args:
        bmp   : unsigned byte array
                with bmp format
        buf   : rectangular
                image buffer
        x1, y1: point to paste
                the buffer

    Returns:
        byref modified unsigned byte array
    """
    if bmp[bmpcolorbits] != buf[0]:
        print(sysmsg['bitsnotequal'])
    else:
        x2 = x1 + buf2int(buf[1: 3])
        y2 = y1 + buf2int(buf[3: 5])
        r = _xchrcnt(bmp)
        if listinBMPrecbnd(bmp, ((x1, y1),
                                (x2, y2))):
            offset = _BMoffset(bmp, x1, y2)
            br = buf2int(buf[5: 7])
            startoff = 7
            endoff = br + 7
            bufsize = len(buf)
            while startoff < bufsize:
                BMPbitBLTput(bmp,
                    offset,
                    buf[startoff: endoff])
                offset += r
                startoff = endoff
                endoff += br
        else:
            print(sysmsg['regionoutofbounds'])


def convertselection2BMP(buf: array):
    """Converts custom unsigned byte array to bmp format

    Args:
        buf: unsigned byte array

    Returns:
        unsigned byte array with bmp format
    """
    bmp = -1
    bits = buf[0]
    if not isvalidcolorbit(bits):
        print (sysmsg['invalidbuf'])
    else:
        bmp = newBMP(
                buf2int(buf[1: 3]) + 1,
                buf2int(buf[3: 5]) + 1,
                bits)
        pasterect(bmp, buf, 0, 0)
    return bmp


def invertimagebits(bmp: array):
    """Inverts the bits in a bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified unsigned byte array
    """
    offset = _hdsz(bmp)
    maxoffset = _flsz(bmp)
    while offset < maxoffset:
        bmp[offset] ^= 255
        offset += 1


def erasealternatehorizontallines(
        bmp: array,
        int_eraseeverynline: int,
        int_eraseNthline: int,
        bytepat: int):
    """Erase every nth line

    Args:
        bmp                : unsigned
                             byte array
                             with
                             bmp format
        int_eraseeverynline: erase every
                             nth line
                             in the
                             region
        int_eraseNthline   : control
                             which line
                             every
                             n lines
                             to erase
        bytepat            : byte pattern
                             to overwrite
                             the erased
                             lines

    Returns:
        byref modified unsigned byte array
    """
    bufsize = _xchrcnt(bmp)
    s1 = _hdsz(bmp)
    s2 = _flsz(bmp) - bufsize
    bytepat &= 0xff
    blank = array('B', [bytepat] * bufsize)
    i = 1
    while s2 > s1:
        if i % int_eraseeverynline == \
               int_eraseNthline:
            bmp[s2: s2 + bufsize] = blank
        s2 -= bufsize
        i += 1


def eraseeverynthhorizontalline(
        bmp: array, n: int):
    """Erase every nth horizontal line

    Args:
        bmp: unsigned byte array
             with bmp format
        n  : erase every nth line

    Returns:
        byref modified unsigned byte array
    """
    erasealternatehorizontallines(
        bmp, n, 0, 0)


@entirecircleinboundary
def erasealternatehorizontallinesincircregion(
        bmp: array,
        x: int, y: int, r: int,
        int_eraseeverynline: int,
        int_eraseNthline: int,
        bytepat: int):
    """Erase every nth line
        in a circular region

    Args:
        bmp                : unsigned
                             byte array
                             with bmp
                             format
        x, y, r            : (x, y)
                             centerpoint
                             and radius r
                             of the
                             circular area
        int_eraseeverynline: erase every
                             nth line
                             in the
                             circular
                             region
        int_eraseNthline   : control which
                             line every
                             n lines
                             to erase
        bytepat            : pattern to
                             overwrite
                             the erased
                             lines

    Returns:
        byref modified unsigned byte array
    """
    c = _getBMoffhdfunc(bmp)
    bytepat &= 0xff
    for v in itercirclepartlineedge(r):
        x1, x2 = mirror(x, v[0])
        y1, y2 = mirror(y, v[1])
        s1 = c(bmp, x1, y1)
        e1 = c(bmp, x2, y1)
        s2 = c(bmp, x1, y2)
        e2 = c(bmp, x2, y2)
        if y1 % int_eraseeverynline == \
                int_eraseNthline:
            bmp[s1: e1] = \
            array('B', [bytepat] * (e1 - s1))
        if y2 % int_eraseeverynline == \
                int_eraseNthline:
            bmp[s2: e2] = \
            array('B', [bytepat] * (e2 - s2))


def eraseeverynthhorizontallineinccircregion(
        bmp: array,
        x: int, y: int, r: int, n: int):
    """Erase every nth horizontal line in a circular region

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
                 of the circular area
        n      : erase every nth line

    Returns:
        byref modified unsigned byte array
    """
    erasealternatehorizontallinesincircregion(
        bmp, x, y, r, n, 0, 0)


@entirerectinboundary
def erasealternatehorizontallinesinregion(
        bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        int_eraseeverynline: int,
        int_eraseNthline: int,
        bytepat: int):
    """Erase every nth line in a rectangular region

    Args:
        bmp                : unsigned
                             byte array
                             with bmp
                             format
        x1, y1, x2, y2     : ints that
                             defines the
                             rectangular
                             region
        int_eraseeverynline: erase every
                             nth line
                             in the region
        int_eraseNthline   : control which
                             line every
                             n lines
                             to erase
        bytepat            : pattern to
                             overwrite
                             the erased
                             lines

    Returns:
        byref modified unsigned byte array
    """
    bytepat &= 0xff
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    bufsize = adjustxbufsize(bmp, x1, x2)
    r = _xchrcnt(bmp)
    f = _getBMoffhdfunc(bmp)
    s1 = f(bmp, x1, y2)
    s2 = f(bmp, x1, y1)
    blank = array('B', [bytepat] * bufsize)
    i = 1
    while s2 > s1:
        if i % int_eraseeverynline == \
               int_eraseNthline:
            bmp[s2: s2 + bufsize] = blank
        s2 -= r
        i += 1


def eraseeverynthhorilineinregion(
        bmp: array,
        x1: int, y1: int,
        x2: int, y2: int, n: int):
    """Erase every nth line in a
        rectangular region

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region
        n             : erase every
                        nth line in the
                        rectangular area

    Returns:
        byref modified
        unsigned byte array
    """
    erasealternatehorizontallinesinregion(
        bmp, x1, y1, x2, y2, n, 0, 0)


def verttrans(bmp: array, trans: str):
    """Do vertical image transforms

    Args:
        bmp : unsigned byte array
              with bmp format
        tran: single letter
              transform code
              'T' - mirror top-half
              'B' - mirror bottom-half
              'F' - flip

    Returns:
        byref modified
        unsigned byte array
    """
    def _F():
        bmp[s1: e1], bmp[s2: e2] = \
        bmp[s2: e2], bmp[s1: e1]

    def _T():
        bmp[s1: s1 + bufsize] = \
        bmp[s2: s2 + bufsize]

    def _B():
        bmp[s2: s2 + bufsize] = \
        bmp[s1: s1 + bufsize]

    f = {'F': _F,
         'T': _T,
         'B': _B}[trans]
    bufsize = _xchrcnt(bmp)
    s1 = _hdsz(bmp)
    s2 = _flsz(bmp) - bufsize
    while s1 < s2:
        e1 = s1 + bufsize
        e2 = s2 + bufsize
        f()
        s1 += bufsize
        s2 -= bufsize


def flipvertical(bmp: array):
    """Does an vertical flip of a bmp

    Args:
        bmp: unsigned byte array
        with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    verttrans(bmp, 'F')


def mirrortop(bmp: array):
    """Mirrors the top-half of a bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    verttrans(bmp, 'T')


def mirrorbottom(bmp: array):
    """Mirrors the bottom-half of a bmp

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    verttrans(bmp, 'B')


@entirerectinboundary
def verttransregion(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int, trans: str):
    """Do vertical image transforms
        in a rectangular region

    Args:
        bmp            : unsigned
                         byte array
                         with bmp format
        x1, y1, x2, y2 : ints that
                         defines the
                         rectangular
                         region
        trans          : single letter
                         transform code
                'T' - mirror top half
                'B' - mirror bottom half
                'F' - flip

    Returns:
        byref modified
        unsigned byte array
    """
    def _F():
        bmp[s1: e1], bmp[s2: e2] = \
        bmp[s2: e2], bmp[s1: e1]

    def _T():
        bmp[s1: s1 + bufsize] = \
        bmp[s2: s2 + bufsize]

    def _B():
        bmp[s2: s2 + bufsize] = \
        bmp[s1: s1 + bufsize]

    f = {'F': _F,
         'T': _T,
         'B': _B}[trans]
    x1, y1, x2, y2 = sortrecpoints(
                        x1, y1, x2, y2)
    bufsize = adjustxbufsize(bmp, x1, x2)
    r = _xchrcnt(bmp)
    c = _getBMoffhdfunc(bmp)
    s1 = c(bmp, x1, y2)
    s2 = c(bmp, x1, y1)
    while s1 < s2:
        e1 = s1 + bufsize
        e2 = s2 + bufsize
        f()
        s1 += r
        s2 -= r


def flipverticalregion(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int):
    """Flips vertical a rectangular area

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    verttransregion(bmp, x1, y1,
                         x2, y2, 'F')


def mirrorbottominregion(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirror the bottom-half
        of a rectangular region

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    verttransregion(bmp, x1, y1,
                         x2, y2, 'B')


def mirrortopinregion(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirror the top-half of
        a rectangular region

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region

    Returns:
        byref modified
        unsigned byte array
    """
    verttransregion(bmp, x1, y1,
                         x2, y2, 'T')

@entirerectinboundary
def fliphorzontalpixelbased(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Flips horizontal
        a rectangular region
        using pixel addressing
        (slightly slow)

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region

    Returns:
        byref modified
        unsigned byte array

    """
    m = x1 + ((x2 - x1) >> 1)
    while x1 <= m:
        y = y1
        while y <= y2:
            swapcolors(bmp, (x1, y),
                            (x2, y))
            y += 1
        x1 += 1
        x2 -= 1


@entirerectinboundary
def fliphverticalalpixelbased(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Flips vertical
        a rectangular region
        using pixel addressing
        (slightly slow)

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: the rectangular
                        region

    Returns:
        byref modified
        unsigned byte array
    """
    m = y1 + ((y2 - y1) >> 1)
    while y1 <= m:
        x = x1
        while x <= x2:
            swapcolors(bmp, (x, y1),
                            (x, y2))
            x += 1
        y1 += 1
        y2 -= 1


@entirerectinboundary
def horizontalbulkswap(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        swapfunc: Callable):
    """Applies function swapfunc
        to a rectangular area

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: the rectangular
                        region

    Returns:
        byref modified
        unsigned byte array
    """
    dx = {24: 1,
           8: 1,
           4: 2,
           1: 8}[bmp[bmpcolorbits]]
    c = _getBMoffhdfunc(bmp)
    r = _xchrcnt(bmp)
    y1, y2 = swapif(y1, y2, y1 > y2)
    x1, x2 = swapif(x1, x2, x1 > x2)
    while x1 < x2:
        swapfunc(bmp,
            c(bmp, x1, y2), c(bmp, x1, y1) + r,
            c(bmp, x2, y2), c(bmp, x2, y1) + r, r)
        x1 += dx
        x2 -= dx


def fliphorizontalregion(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int):
    """Does a horizontal flip
        of a rectangular area

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region

    Returns:
        byref modified
        unsigned byte array
    """
    def _24(bmp, s1, e1, s2, e2, r):
        bmp[s1: e1 - 2: r], bmp[s2: e2 - 2: r], bmp[s1 + 1: e1 - 1: r], bmp[s2 + 1: e2 - 1: r], bmp[s1 + 2: e1: r], bmp[s2 + 2: e2: r] = \
        bmp[s2: e2 - 2: r], bmp[s1: e1 - 2: r], bmp[s2 + 1: e2 - 1: r], bmp[s1 + 1: e1 - 1: r], bmp[s2 + 2: e2: r], bmp[s1 + 2: e1: r]

    def _8(bmp, s1, e1, s2, e2, r):
        bmp[s1: e1: r], bmp[s2: e2: r]= \
        bmp[s2: e2: r], bmp[s1: e1: r]

    def _4(bmp, s1, e1, s2, e2, r):
        bmp[s1: e1: r], bmp[s2: e2: r] = \
         flipnibbleinbuf(bmp[s2: e2: r]),  flipnibbleinbuf(bmp[s1: e1: r])

    def _1(bmp, s1, e1, s2, e2, r):
        bmp[s1: e1: r], bmp[s2: e2: r] = \
        rotatebitsinbuf(bmp[s2: e2: r]), rotatebitsinbuf(bmp[s1: e1: r])

    horizontalbulkswap(
        bmp,
        x1, y1, x2, y2,
        {24: _24,
          8: _8,
          4: _4,
          1: _1}[bmp[bmpcolorbits]])


def mirrorleftinregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Mirrors the left-half
        of a rectangular area

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    def _24(bmp, s1, e1, s2, e2, r):
        bmp[s2: e2 - 2: r], bmp[s2 + 1: e2 - 1: r], bmp[s2 + 2: e2: r]= \
        bmp[s1: e1 - 2: r], bmp[s1 + 1: e1 - 1: r], bmp[s1 + 2: e1: r]

    def _8(bmp, s1, e1, s2, e2, r):
        bmp[s2: e2: r] = bmp[s1: e1: r]

    def _4(bmp, s1, e1, s2, e2, r):
        bmp[s2: e2: r] = \
             flipnibbleinbuf(bmp[s1: e1: r])

    def _1(bmp, s1, e1, s2, e2, r):
        bmp[s2: e2: r] = \
            rotatebitsinbuf(bmp[s1: e1: r])

    horizontalbulkswap(bmp,
        x1, y1, x2, y2,
        {24: _24,
          8: _8,
          4: _4,
          1: _1}[bmp[bmpcolorbits]])


def mirrorrightinregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Mirrors the right-half of
        a rectangular area in a bitmap

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    def _24(bmp, s1, e1, s2, e2, r):
        bmp[s1: e1 - 2: r], bmp[s1 + 1: e1 - 1: r], bmp[s1 + 2: e1: r]=\
        bmp[s2: e2 - 2: r], bmp[s2 + 1: e2 - 1: r], bmp[s2 + 2: e2: r]

    def _8(bmp, s1, e1, s2, e2, r):
        bmp[s1: e1: r] = bmp[s2: e2: r]

    def _4(bmp, s1, e1, s2, e2, r):
        bmp[s1: e1: r] = \
             flipnibbleinbuf(bmp[s2: e2: r])

    def _1(bmp, s1, e1, s2, e2, r):
        bmp[s1: e1: r] = \
            rotatebitsinbuf(bmp[s2: e2: r])

    horizontalbulkswap(
        bmp, x1, y1, x2, y2,
        {24: _24,
          8: _8,
          4: _4,
          1: _1}[bmp[bmpcolorbits]])


def mirrorleft(bmp: array):
    """Mirrors the left-half of an
        in-memory bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    mirrorleftinregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1)


def mirrorright(bmp: array):
    """Mirrors the right-half of an
        in-memory bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    mirrorrightinregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1)


def mirrortopleftinregion(
        bmp: array, x1: int, y1: int,
                   x2: int, y2: int):
    """Mirrors the top-left of a
        rectangular region defined by
        (x1, y1) and (x2, y2)

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    mirrorleftinregion(
        bmp, x1, y1, x2, y2)
    mirrortopinregion(
        bmp, x1, y1, x2, y2)


def mirrortoprightinregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Mirrors the top-right of a
        rectangular region defined by
        (x1, y1) and (x2, y2)

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    mirrorrightinregion(
        bmp, x1, y1, x2, y2)
    mirrortopinregion(
        bmp, x1, y1, x2, y2)


def mirrorbottomleftinregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Mirrors the bottom-left of a
        rectangular region defined
        by (x1, y1) and (x2, y2)

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    mirrorleftinregion(
        bmp, x1, y1, x2, y2)
    mirrorbottominregion(
        bmp, x1, y1, x2, y2)


def mirrorbottomrightinregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Mirrors the bottom-right of a
        rectangular region defined by
        (x1, y1) and (x2, y2)

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: the rectangular
                        region

    Returns:
        byref modified
        unsigned byte array
    """
    mirrorrightinregion(
        bmp, x1, y1, x2, y2)
    mirrorbottominregion(
        bmp, x1, y1, x2, y2)


def mirrortopleft(bmp):
    """Mirrors the top-left part
        of an in-memory bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    mirrorleftinregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
        (getmaxy(bmp) - 1)//2)
    mirrortop(bmp)


def mirrortopright(bmp: array):
    """Mirrors the top-right part
        of an in-memory bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    mirrorrightinregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
       (getmaxy(bmp) - 1) // 2)
    mirrortop(bmp)


def mirrorbottomleft(bmp: array):
    """Mirrors the bottom-left part
        of an in-memory bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    ymax = getmaxy(bmp) - 1
    mirrorleftinregion(bmp, 0,
        ymax // 2, getmaxx(bmp) - 1,
        ymax)
    mirrorbottom(bmp)


def mirrorbottomright(bmp: array):
    """Mirrors the bottom right part
        of an in-memory bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    ymax = getmaxy(bmp) - 1
    mirrorrightinregion(bmp, 0,
        ymax // 2, getmaxx(bmp) - 1,
        ymax)
    mirrorbottom(bmp)


def fliphorizontal(bmp: array):
    """Does a horizontal flip of an
        in-memory bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    fliphorizontalregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1)


def flipXY(bmp: array):
    """Flips the x and y coordinates of
        an in-memory bitmap for a
        90 degree rotation

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    mx = getmaxy(bmp)
    my = getmaxx(bmp)
    bits = bmp[bmpcolorbits]
    nbmp = newBMP(mx, my, bits)
    if bits in [8, 24]:
        r = _xchrcnt(nbmp)
        mx -= 1
        if bits < 24:
            copyRGBpal(bmp, nbmp)
        offset = 0
        for y in range(my):
            BMPbitBLTput(nbmp,
                offset,
                array('B', vertBMPbitBLTget(bmp, y, 0, mx)))
            offset += r
    else:
        copyRGBpal(bmp, nbmp)
        for v in iterimagecolor(bmp,
                    sysmsg['flipXY'], '*',
                    sysmsg['done']):
            plotxybit(nbmp, v[0][1],
                            v[0][0],
                            v[1])
    return nbmp


@entirerectinboundary
def itergetcolorfromrectregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Yields color info of
        a rectangular area

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangle

    Yields:
        ((x: int, y: int), color: int)
        for all points in area
    """
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    x, y = x1, y1
    while y <= y2:
        x = x1
        while x <= x2:
            yield ((x, y), getxybit(
                            bmp, x, y))
            x += 1
        y += 1


@entirerectinboundary
def crop(bmp: array, x1: int, y1: int,
                     x2: int, y2: int):
    """Crops the image to a rectangular
        region defined by (x1, y1)
                      and (x2, y2)

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: the rectangular
                        region

    Returns:
        unsigned byte array
        with bitmap layout
    """
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    bits = bmp[bmpcolorbits]
    nbmp = newBMP(x2 - x1 + 1,
                  y2 - y1 + 1, bits)
    if bits < 8:
        copyRGBpal(bmp, nbmp)
        for v in itergetcolorfromrectregion(
                    bmp, x1, y1, x2, y2):
            intplotvecxypoint(nbmp,
                subvect(v[0], (x1, y1)), v[1])
    else:
        offset = 0
        r = _xchrcnt(nbmp)
        for buf in itercopyrect(
                        bmp, x1, y1, x2, y2):
            BMPbitBLTput(nbmp, offset, buf)
            offset += r
    return nbmp


@entirerectinboundary
def invertregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Inverts the bits in a
        rectangular region defined
        by (x1,y1) and (x2,y2)

    Args:
        bmp          :  unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: the rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    bits = bmp[bmpcolorbits]
    if bits < 8:
        c = getmaxcolors(bmp) - 1
        for v in itergetcolorfromrectregion(
                    bmp, x1, y1, x2, y2):
            intplotvecxypoint(bmp,
                v[0], getxybitvec(
                        bmp, v[0]) ^ c)
    else:
        offset = iif(bits == 24,
                    _24bmof(bmp, x1, y2),
                    _8bmof(bmp, x1, y2))
        r = _xchrcnt(bmp)
        for buf in itercopyrect(
                        bmp, x1, y1, x2, y2):
            BMPbitBLTput(bmp,
                offset,  invertbitsinbuffer(buf))
            offset += r


def monofilterto24bitregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int):
    """Applies a monochrome filter
        to a rectangular area
        defined by (x1, y1) and
        (x2, y2) in a 24 bit bitmap

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: the rectangle

    Returns:
        byref modified
        unsigned byte array
    """
    _usebyrefnopar24bitfn2reg(
        bmp, x1, y1, x2, y2,
        applymonochromefiltertoBGRbuf)


def monofilterto24bitimage(bmp: array):
    """Applies a mono filter
        to a 24 bit in-memory bitmap

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified
        unsigned byte array
    """
    monofilterto24bitregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1)


def horizontalbrightnessgradto24bitimage(
        bmp: array,
        lumrange: list[int, int]):
    """Applies a horizontal brightness
        gradient to a 24-bit bitmap

    Args:
        bmp     : unsigned byte array
                  with bmp format
        lumrange: [byte,byte]
                  the range of the
                  luminosity gradient

    Returns:
        byref modified
        unsigned byte array
    """
    horizontalbrightnessgradto24bitregion(
        bmp, 0, 0, getmaxx(bmp) - 1,
        getmaxy(bmp) - 1, lumrange)


def resizeNtimesbigger(
        bmp: array, n: int):
    """Resize an in-memory bmp
        n times bigger

    Args:
        buf : array to resize
        n   : resize factor
        bits: bit depth of
              the color info
              (1, 4, 8, 24)

    Returns:
        unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    (mx, my) = getmaxxy(bmp)
    nx = mx * n
    ny = my * n
    nbmp = newBMP(nx, ny, bits)
    mx -= 1
    my -= 1
    ny -= 1
    if bits < 24:
        copyRGBpal(bmp, nbmp)
    offset = _BMoffset(nbmp, 0, ny)
    r = _xchrcnt(nbmp)
    for buf in itercopyrect(bmp, 0, 0,
                        mx, my):
        nbuf = resizebufNtimesbigger(buf, n, bits)
        for _ in range(n):
            BMPbitBLTput(nbmp, offset, nbuf)
            offset += r
    return nbmp


def colorfilterto24bitregion(
        bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        rgbfactors: list[float,
                         float,
                         float]):
    """Applies a color filter to
        a rectangular area defined
        by (x1, y1) and (x2, y2)
        in a 24-bit bitmap

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region
        rgbfactors    : color filter
                        r, g and b
                        range from
                        0.0 to 1.0

    Returns:
        byref modified
        unsigned byte array

    """
    _usebyref24btfn2reg(
        bmp, x1, y1, x2, y2,
        applycolorfiltertoBGRbuf, rgbfactors)


def colorfilterto24bitimage(
        bmp: array,
        rgbfactors: list[float,
                         float,
                         float]):
    """Applies a color filter
        to a whole image in
        an in-memory 24 bit bitmap

    Args:
        bmp       : unsigned byte array
                    with bmp format
        rgbfactors: color filter
                    r,g and b values
                    are from 0.0 to 1.0

    Returns:
        byref modified
        unsigned byte array
    """
    colorfilterto24bitregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1, rgbfactors)


def brightnesseadjto24bitregion(
        bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        percentadj: float):
    """Brightness adjustment to a rectangular area

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines
                        the rectangle
        percentadj    : float percentage
                        brightness adjust
                        can be positive
                        or negative

    Returns:
        byref modified unsigned byte array
    """
    _use24bitfn2reg(bmp, x1, y1, x2, y2,
         applybrightnessadjtoBGRbuf, percentadj)


def thresholdadjto24bitregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int,
        lumrange: list[int, int]):
    """Threshold adjustment to a rectangular area

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: ints that defines
                        the rectangular
                        region
        lumrange      : (byte:byte)
                        threshold
                        adjustment
                        luminosity range

    Returns:
        byref modified unsigned byte array
    """
    _use24bitfn2reg(bmp, x1, y1, x2, y2,
        applythresholdadjtoBGRbuf, lumrange)


def thresholdadjcircregion(
        bmp: array,
        x: int, y: int, r: int,
        lumrange: list[int, int]):
    """Threshold adjustment to a circular area

    Args:
        bmp      : unsigned byte array
                   with bmp format
        x, y, r  : centerpoint (x, y)
                   and radius r
        lumrange : (byte: byte)
                   threshold adjustment
                   luminosity range

    Returns:
        byref modified unsigned byte array
    """
    _use24btfn2circreg(bmp, x, y, r,
        applythresholdadjtoBGRbuf, lumrange)


def brightnesseadjto24bitimage(
        bmp: array, percentadj: float):
    """Brightness adjustment to a whole BMP

    Args:
        bmp       : unsigned byte array
                    with bmp format
        percentadj: float percentage
                    brightness
                    adjustment
                    can be positive
                    or negative

    Returns:
        byref modified unsigned byte array
    """
    brightnesseadjto24bitregion(
        bmp, 0, 0, getmaxx(bmp) - 1,
        getmaxy(bmp) - 1, percentadj)


def thresholdadjto24bitimage(
        bmp: array,
        lumrange: list[int, int]):
    """Threshold adjustment to a whole BMP

    Args:
        bmp     : unsigned byte array
                  with bmp format
        lumrange: (byte: byte)
                  threshold adjustment
                  luminosity range

    Returns:
        byref modified unsigned byte array
    """
    thresholdadjto24bitregion(
        bmp, 0, 0, getmaxx(bmp) - 1,
        getmaxy(bmp) - 1, lumrange)


def verticalbrightnessgradto24bitimage(
        bmp: array,
        lumrange: list[int, int]):
    """Applies a vertical brightness gradient

    Args:
        bmp     : unsigned byte array
                  with bmp format
        lumrange: (byte: byte) the
                  brightness gradient

    Returns:
        byref modified unsigned byte array
    """
    verticalbrightnessgradto24bitregion(
        bmp, 0, 0, getmaxx(bmp) - 1,
        getmaxy(bmp) - 1, lumrange)


def plotmultifractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        func: Callable,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi Fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : power to raise z to
        func          : fractal function
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    maxcolors = getmaxcolors(bmp)
    mcolor = maxcolors - 1
    bits = bmp[bmpcolorbits]
    if bits == 24:
        for (x, y, c) in func(x1, y1, x2, y2, d,
                          domain, maxiter):
            plotxybit(bmp, x, y,
            colormix(((255 - c) * 20) % 256,
                        RGBfactors))
    elif bits == 8 and len(RGBfactors) == 3:
        for (x, y, c) in func(x1, y1, x2, y2, d,
                          domain, maxiter):
            bmp[_8bmofhd(bmp, x, y)] = c
    elif bits == 4 and len(RGBfactors) == 3:
        for (x, y, c) in func(x1, y1, x2, y2, d,
                          domain, maxiter):
            plotxybit(bmp, x, y, c)
    elif bits == 1:
        for (x, y, c) in func(x1, y1, x2, y2, d,
                          domain, maxiter):
            plotxybit(bmp, x, y, c % maxcolors)
    else:
        for (x, y, c) in func(x1, y1, x2, y2, d,
                          domain, maxiter):
            plotxybit(bmp, x, y, (mcolor - c) % maxcolors)


def plotmultifractalcomplexpar(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        func: Callable,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multifractal with a complex parameter

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        func          : fractal function
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    maxcolors = getmaxcolors(bmp)
    mcolor = maxcolors - 1
    bits = bmp[bmpcolorbits]
    if bits == 24:
        for (x, y, v) in func(x1, y1, x2, y2, c, d,
                           domain, maxiter):
            plotxybit(bmp, x, y,
            colormix(((255 - v) * 20) % 256,
                        RGBfactors))
    elif bits == 8 and len(RGBfactors) == 3:
        for (x, y, v) in func(x1, y1, x2, y2, c, d,
                          domain, maxiter):
            bmp[_8bmofhd(bmp, x, y)] = v
    elif bits == 4 and len(RGBfactors) == 3:
        for (x, y, v) in func(x1, y1, x2, y2, c, d,
                          domain, maxiter):
            plotxybit(bmp, x, y, v)
    elif bits == 1:
        for (x, y, v) in func(x1, y1, x2, y2, c, d,
                           domain, maxiter):
            plotxybit(bmp, x, y, v % maxcolors)
    else:
        for (x, y, v) in func(x1, y1, x2, y2, c, d,
                           domain, maxiter):
            plotxybit(bmp, x, y, (mcolor - v) % maxcolors)


def newtonsfractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: list[Callable, Callable],
        domain: list[float, float, float, float],
        RGBfactorslist: list[list[float, float, float]],
        maxiter: int) -> list:
    """Draw Newtons Fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : function pair
                        (func, derivative func)
        domain        : coordinates in real
                        and imaginary plane
        rgbfactorslist: [[r, g, b],...] values
                        range from 0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
        list of roots
    """
    tol = 10e-3
    maxcolors = getmaxcolors(bmp)
    mcolor = maxcolors - 1
    roots = []
    for (x, y, r) in iternewtonsfractal(x1, y1,
                     x2, y2, d, domain, maxiter):
        rt, c = r
        if rt != None:
            flag = False
            for tst_rt in roots:
                if abs(tst_rt - rt) < tol:
                    rt = tst_rt
                    flag = True
                    break
            if not flag:
                roots.append(rt)
            if bmp[bmpcolorbits] == 24:
                c = colormix(((255 - c) * 20) % 256,
                    RGBfactorslist[roots.index(rt)])
            else:
                c = mcolor - c % maxcolors
            plotxybit(bmp, x, y, c)
    return roots


def mandelbrot(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        domain    : list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Mandelbrot set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, 2,
        itermultibrot, domain,
        RGBfactors, maxiter)


def multibrot(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multibrot set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, d,
        itermultibrot, domain,
        RGBfactors, maxiter)


def multicircle(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multicircle fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, d,
        itermulticircle, domain,
        RGBfactors, maxiter)


def multihyperbola(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multihyperbola fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, d,
        itermultihyperbola, domain,
        RGBfactors, maxiter)


def multisuperellipse(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multisuperellipse fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, d,
        itermultisuperellipse, domain,
        RGBfactors, maxiter)


def astroidfractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        r: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw an Astroid fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        r             : radius
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, r,
        iterastroid, domain,
        RGBfactors, maxiter)


def lemniscatefractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        r: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Lemniscate fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        r             : radius
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, r,
        iterlemniscate, domain,
        RGBfactors, maxiter)


def xorfractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Xor fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : int modulo
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, d,
        iterxorfractal, domain,
        RGBfactors, maxiter)


def xordivfractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Xor int div fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : int div
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, d,
        iterxordivfractal, domain,
        RGBfactors, maxiter)


def tetrationfractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Tetration Fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : threshold
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, d,
        itertetration, domain,
        RGBfactors, maxiter)


def sinjulia(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Sin(z) Julia Set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, c,
        itersinjulia, domain,
        RGBfactors, maxiter)


def cosjulia(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Cos(z) Julia Set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, c,
        itercosjulia, domain,
        RGBfactors, maxiter)


def spiraljulia(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Spiral Julia Set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, c,
        iterspiraljulia, domain,
        RGBfactors, maxiter)


def lambdafractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Lambda Fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(bmp, x1, y1, x2, y2, c,
        iterlamdbafractal, domain,
        RGBfactors, maxiter)


@functimer
def savelambdafractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Lambda fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file,
        x, y,
        lambdafractal,
        c,
        domain,
        rgbfactors,
        bitdepth,
        maxiter)


@functimer
def savetetrationfractal2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Tetration Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : threshold
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        tetrationfractal,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savexorfractal2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 24,
        maxiter: int = 255):
    """Saves a Xor Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : int modulo
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        xorfractal,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savexordivfractal2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 24,
        maxiter: int = 255):
    """Saves a Xor int div Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : int div
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        xordivfractal,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savesinjulia2file(
        file: str,
        x: int, y: int,
        c: complex,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Sin(z) Julia Set to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        sinjulia,
        c, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savespiraljulia2file(
        file: str,
        x: int, y: int,
        c: complex,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Spiral Julia Set to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        spiraljulia,
        c, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savecosjulia2file(
        file: str,
        x: int, y: int,
        c: complex,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Cos(z) Julia Set to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        cosjulia,
        c, domain,
        rgbfactors, bitdepth,
        maxiter)


def julia(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Julia set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, 2,
        itermultijulia, domain,
        RGBfactors, maxiter)


def multijulia(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multijulia set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultijulia, domain,
        RGBfactors, maxiter)


def ngonfractal(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: float,
        n: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a N-gon fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : float constant
        n             : float exponent (sides)
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, n,
        iterngonfractal, domain,
        RGBfactors, maxiter)


def multibiomorphvariant(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi Biomorph Variant fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultibiomorphvariant, domain,
        RGBfactors, maxiter)


def multibiomorphphasevariant(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi Biomorph Phase Variant fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultibiomorphphasevariant, domain,
        RGBfactors, maxiter)


def multibiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultibiomorph, domain,
        RGBfactors, maxiter)


def multisinbiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi sin(z) Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultisinbiomorph, domain,
        RGBfactors, maxiter)


def multicosbiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi cos(z) Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermulticosbiomorph, domain,
        RGBfactors, maxiter)


def multicoshbiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi cosh(z) Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermulticoshbiomorph, domain,
        RGBfactors, maxiter)


def multisinhbiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi sinh(z) Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultisinhbiomorph, domain,
        RGBfactors, maxiter)


def multiexpbiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi exp(z) Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultiexpbiomorph, domain,
        RGBfactors, maxiter)


def multi2ndtetrationbiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi 2nd Tetration Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermulti2ndtetrationbiomorph,
        domain, RGBfactors, maxiter)


def multizconjugatebiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi z conjugate Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultizconjugatebiomorph, domain,
        RGBfactors, maxiter)


def multitanbiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi tan(z) Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultitanbiomorph, domain,
        RGBfactors, maxiter)


def multitanhbiomorph(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multi tanh(z) Biomorph fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        c             : complex number
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractalcomplexpar(
        bmp, x1, y1, x2, y2, c, d,
        itermultitanhbiomorph, domain,
        RGBfactors, maxiter)


def tricorn(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Tricorn set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(
        bmp, x1, y1, x2, y2, 2,
        itermulticorn, domain,
        RGBfactors, maxiter)


def multicorn(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Multicorn set

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : power to raise z to
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(
        bmp, x1, y1, x2, y2, d,
        itermulticorn, domain,
        RGBfactors, maxiter)


def barnsleytree(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: complex,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Barnsley Tree Fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : complex number
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(
        bmp, x1, y1, x2, y2, d,
        iterbarnsleytree, domain,
        RGBfactors, maxiter)


def marekdragon(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        d: float,
        domain: list[float, float, float, float],
        RGBfactors: list[float, float, float],
        maxiter: int):
    """Draw a Marek Dragon Fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular area
                        to draw in
        d             : irrational number
        domain        : coordinates in real
                        and imaginary plane
        rgbfactors    : [r, g, b] values
                        range from
                        0.0 to 1.0
        maxiter       : when to break
                        color compute

    Returns:
        byref modified unsigned byte array
    """
    plotmultifractal(
        bmp, x1, y1, x2, y2, d,
        itermarekdragon, domain,
        RGBfactors, maxiter)


def IFS(bmp: array,
        IFStransparam: tuple,
        x1: int, y1: int,
        x2: int, y2: int,
        xscale: int, yscale: int,
        xoffset: int, yoffset: int,
        color: int, maxiter: int):
    """Draw an Interated Function System Fractal

    Args:
        bmp            : unsigned
                         byte array
                         with bmp format
        IFStransparam  : see fractals.py
        x1, y1, x2, y2 : rectangular
                         region
                         to draw in
        xscale,yscale  : scaling factors
        xoffset,yoffset: used to move
                         the fractal
        color          : color of fractal
        maxiter        : when to break
                         color compute

    Returns:
        byref modified unsigned byte array
    """
    for (px, py) in iterIFS(IFStransparam,
        x1, y1, x2, y2, xscale, yscale,
        xoffset, yoffset, maxiter):
        plotxybit(bmp, px, py, color)


def plotflower(bmp: array,
        cx: int, cy: int, r: int,
        petals: float, angrot: float,
        lumrange: list[int, int],
        RGBfactors: list[float, float, float]):
    """Draw a flower

    Args:
        bmp       : unsigned byte array
                    with bmp format
        cx, cy, r : center (cx,cy)
                    and radius r
        petals    : number of petals
        angrot    : angle of rotation
        lumrange  : (byte:byte) range
                    of brightness for
                    the gradient
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0

    Returns:
        byref modified unsigned byte array
    """
    maxcolors = getmaxcolors(bmp)
    mcolor = maxcolors - 1
    lum1, dlum = range2baseanddelta(lumrange)
    for (x, y) in iterflower(cx, cy, r, petals, angrot):
        c = colormix(setmax(abs(int(lum1 + dlum * (dist([x, y], [cx, cy]) / r))),255), RGBfactors)
        if bmp[bmpcolorbits] != 24:
            c = mcolor - c % maxcolors
        plotxybit(bmp, x, y, c)


def newBMPwithPAL(x, y, bitdepth, rgbfactors):
    """Create a new in memory BMP with palette
    based on rgbfactors and bitdepth

    Args:
        x       : width of bitmap
        y       : height of bitmap
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
    """
    bmp = newBMP(x, y, bitdepth)
    if bitdepth < 24:
        cp = len(rgbfactors)
        if cp == 3:
            setBMP2invshiftedmonochrome(bmp, rgbfactors)
        elif cp == 2:
            setBMP2dichromaticpal(bmp, rgbfactors[0], rgbfactors[1])
    return bmp


def savefractal2file(
        file: str,
        x: int, y: int,
        f: Callable,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        f       : fractal function
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    bmp = newBMPwithPAL(x, y, bitdepth, rgbfactors)
    f(bmp, 0, 0, x, y,
    domain, rgbfactors, maxiter)
    saveBMP(file, bmp)


def savemultifractal2file(
        file: str,
        x: int, y: int,
        f: Callable,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 24,
        maxiter: int = 255):
    """Saves a Multi Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        f       : fractal function
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    bmp = newBMPwithPAL(x, y, bitdepth, rgbfactors)
    f(bmp, 0, 0, x, y, d,
    domain, rgbfactors, maxiter)
    saveBMP(file, bmp)


@functimer
def savemandelbrotfractal2file(
        file: str,
        x: int, y: int,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Mandelbrot Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractal2file(
        file, x, y,
        mandelbrot,
        domain, rgbfactors,
        bitdepth, maxiter)


@functimer
def savemultibrotfractal2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Multibrot Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        multibrot,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savemulticirclefractal2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Multicircle Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        multicircle,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savemultihyperbolafractal2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Multihyperbola Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        multihyperbola,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savemultisuperellipsefractal2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 24,
        maxiter: int = 255):
    """Saves a Multisuperellipse Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        multisuperellipse,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def saveastroidfractal2file(
        file: str,
        x: int, y: int,
        r: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves an Astroid Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        r       : radius
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        astroidfractal,
        r, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savelemniscatefractal2file(
        file: str,
        x: int, y: int,
        r: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Lemniscate Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        r       : radius
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        lemniscatefractal,
        r, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savetricornfractal2file(
        file: str,
        x: int, y: int,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Tricorn Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractal2file(
        file, x, y,
        tricorn,
        domain, rgbfactors,
        bitdepth, maxiter)


@functimer
def savemulticornfractal2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Multicorn Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        multicorn,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savebarnsleytreefractal2file(
        file: str,
        x: int, y: int,
        d: complex,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Barnsley Tree Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : complex number
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        barnsleytree,
        d, domain,
        rgbfactors, bitdepth,
        maxiter)


@functimer
def savemarekdragon2file(
        file: str,
        x: int, y: int,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Marek Dragon Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : irrational number
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savemultifractal2file(
        file, x, y,
        marekdragon,
        d, domain,
        rgbfactors,
        bitdepth,
        maxiter)


def savefractalwithparam2file(
        file: str,
        x: int, y: int,
        f: Callable,
        p: any,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 24,
        maxiter: int = 255):
    """Saves a Fractal with a parameter to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        f       : fractal function
        p       : any
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    bmp = newBMP(x, y, bitdepth)
    if bitdepth < 24:
        cp = len(rgbfactors)
        if   cp == 3:
            setBMP2invshiftedmonochrome(bmp, rgbfactors)
        elif cp == 2:
            setBMP2dichromaticpal(bmp, rgbfactors[0],
                                       rgbfactors[1])
    retval = f(bmp, 0, 0, x, y, p,
    domain, rgbfactors, maxiter)
    saveBMP(file, bmp)
    return retval


def savefractalwith2param2file(
        file: str,
        x: int, y: int,
        f: Callable,
        p1: any,
        p2: any,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 24,
        maxiter: int = 255):
    """Saves a Fractal with two parameters to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        f       : fractal function
        p1, p1  : any (param1, param2)
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    bmp = newBMPwithPAL(x, y, bitdepth, rgbfactors)
    retval = f(bmp, 0, 0, x, y, p1, p2,
    domain, rgbfactors, maxiter)
    saveBMP(file, bmp)
    return retval


@functimer
def savejuliafractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Julia Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwithparam2file(file, x, y,
        julia,
        c, domain, rgbfactors, bitdepth,
        maxiter)


@functimer
def savemultijuliafractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Multi Julia Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multijulia,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savengonfractal2file(
        file: str,
        x: int, y: int,
        c: float,
        n: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a n gon Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : float constant
        n       : float exponent (sides)
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    ngonfractal,
    c, n, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultibiomorphvariantfractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Multi Biomorph Variant Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multibiomorphvariant,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultibiomorphphasevariantfractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 8,
        maxiter: int = 255):
    """Saves a Multi Biomorph Phase Variant Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multibiomorphphasevariant,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultibiomorphfractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multibiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultisinbiomorphfractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi sin(z) Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multisinbiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemulticosbiomorphfractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi cos(z) Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multicosbiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemulticoshbiomorphfractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi cosh(z) Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multicoshbiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultisinhbiomorphfractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi sinh(z) Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multisinhbiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultiexpbiomorphfractal2file(
        file: str,
        x: int, y: int,
        c: complex,
        d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi exp(z) Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multiexpbiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemulti2ndtetrationbiomorphfractal2file(
        file: str, x: int, y: int,
        c: complex, d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi 2nd Tetration Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multi2ndtetrationbiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultizconjugatebiomorphfractal2file(
        file: str, x: int, y: int,
        c: complex, d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi z conjugate Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multizconjugatebiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultitanbiomorphfractal2file(
        file: str, x: int, y: int,
        c: complex, d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi tan(z) Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multitanbiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savemultitanhbiomorphfractal2file(
        file: str, x: int, y: int,
        c: complex, d: float,
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 4,
        maxiter: int = 255):
    """Saves a Multi tanh(z) Biomorph Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        c       : complex number
        d       : power to raise z to
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                        or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    savefractalwith2param2file(file, x, y,
    multitanhbiomorph,
    c, d, domain, rgbfactors, bitdepth,
    maxiter)


@functimer
def savenewtonsfractal2file(
        file: str, x: int, y: int,
        d: list[Callable, Callable],
        domain: list[float, float, float, float],
        rgbfactors: list[float, float, float],
        bitdepth: int = 24,
        maxiter: int = 255):
    """Saves a Newtons Fractal to a file

    Args:
        file    : full path to new file
        x       : width of bitmap
        y       : height of bitmap
        d       : [function and,
                  function derivative pair]
        domain  : location in real and
                  imaginary plane
                  (minreal, maxreal,
                   minimag, maximag)
        rgbfactors: [r, g, b] values
                    all range from
                    0.0 to 1.0
                         or
                    [int, int] values
                    packed RGB
        bitdepth: optional parameter
                  for bit depth
                  (1, 4, 8, 24) bits
        maxiter : optional parameter
                  to set maximum iteration

    Returns:
        a bitmap file
    """
    return savefractalwithparam2file(
        file, x, y,
        newtonsfractal,
        d, domain, rgbfactors,
        bitdepth, maxiter)


@functimer
def savehilbertcurve2file(
        file: str, x: int, y: int,
        order: int,
        bitdepth: int = 1,
        color: int = 1,
        backgroundcolor: int = 0,
        penradius: int = 1,
        orientation: int = 0
        ):
    """Saves a Hilbert Curve to a file

    Args:
        file     : full path to new file
        x        : width of bitmap
        y        : height of bitmap
        order    : order of the curve
        bitdepth : optional bit depth
                   (1, 4, 8, 24)
                   default = 1
        color    : optional color
                   default = 1
        backgroundcolor: optional backgroundcolor
                         default = 0
        penradius: optional penradius
                   default = 1
        orientation: can be 0 or 1
                     default = 1

    Returns:
        a bitmap file
    """
    bmp = newBMP(x, y, bitdepth)
    x -= 1
    y -= 1
    if backgroundcolor > 0:
        filledrect(bmp, 0, 0, x, y, backgroundcolor)
    o =  (0, x, y, 0) if orientation == 1 \
    else (x, 0, 0, y)
    l = []
    hilbertvert(l, (0, 0), o, order)
    plotlines(bmp, l, color, penradius)
    saveBMP(file, bmp)


@functimer
def savekochsnowflake2file(
        file: str, r: int, order: int,
        angle: float = 0,
        bitdepth: int = 1,
        color: int = 1,
        backgroundcolor: int = 0,
        penradius: int = 1,
        ):
    """Saves a Koch snowflake to a file

    Args:
        file     : full path to new file
        r        : radius of the snowflake
        order    : order of the curve
        angle    : rotation of snowflake
                   around center in degrees
        bitdepth : optional bit depth
                   (1, 4, 8, 24)
                   default = 1
        color    : optional color
                   default = 1
        backgroundcolor: optional backgroundcolor
                         default = 0
        penradius: optional penradius
                   default = 1

    Returns:
        a bitmap file
    """
    x = y = (r << 1) + (penradius << 1) + 2
    bmp = newBMP(x, y, bitdepth)
    x -= 1
    y -= 1
    if backgroundcolor > 0:
        filledrect(bmp, 0, 0, x, y, backgroundcolor)
    cx = x >> 1
    cy = y >> 1
    l = kochsnowflakevert(cx, cy, r, angle, order)
    if penradius > 1:
        thickplotpoly(bmp, l, penradius, color)
    else:
        plotpoly(bmp, l, color)
    saveBMP(file, bmp)


def plotfilledflower(bmp: array,
        cx: int, cy: int, r: int,
        petals: float, angrot: float,
        lumrange: list[int, int],
        RGBfactors: list[float,
                         float,
                         float]):
    """Draw a filled flower

    Args:
        bmp       : unsigned byte array
                    with bmp format
        cx, cy, r : center (cx, cy)
                    and radius r
        petals    : number of petals
        angrot    : angle of rotation
        lumrange  : (byte:byte) range
                    of brightness
        rgbfactors: [r, g, b] values
                    of r, g and b
                    range from
                    0.0 to 1.0

    Returns:
        byref modified unsigned byte array

    """
    for nr in range (r, 2, -1):
        plotflower(bmp, cx, cy, nr,
            petals, angrot,
            lumrange, RGBfactors)


def plotbmpastext(bmp: array):
    """Plot a bitmap as text
    (cannot output 24-bit bmp)

    Args:
        bmp : unsigned byte array
              with bmp format

    Returns:
        console text output
        for debug and ascii art

    """
    bits = _getclrbits(bmp)
    my = getmaxy(bmp) - 1
    r = _xchrcnt(bmp)
    offset = _hdsz(bmp)
    for y in range(my, 0, -1):
        for x in range(r):
            i = offset + r * y + x
            if bits == 1:
                plotbitsastext(bmp[i])
            if bits == 4:
                c0, c1 = divmod(bmp[i], 16)
                print(chr(97 + c0) +
                      chr(97 + c1), end='')
            if bits == 8:
                print(chr(bmp[i]), end='')
        print()


def piechart(bmp: array,
        x: int, y: int, r: int,
        dataandcolorlist: list):
    """Draw a piechart

    Args:
        bmp             : unsigned
                          byte array
                          with bmp format
        x, y, r         : center (x, y)
                          and radius r
        dataandcolorlist: stuff to plot
                          + color

    Returns:
        byref modified unsigned byte array

    """
    alist, big = genpiechartdata(dataandcolorlist)
    if big > -1: #for speed more computions in drawarc
            circle(bmp, x, y, r,
                alist[big][2], True)
    for a in alist:
        if a[4] < 50:
            drawarc(bmp, x, y, r,
                a[0], a[1], a[2],
                True, a[2], True)
    return [alist, big]


@func24bitonlyandentirerectinboundary
def _usebyrefnopar24bitfn2reg(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int,
        func: Callable):
    """Apply a func to a rectangular area in a 24-bit bitmap

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: the rectangular
                        area
        func          : user defined
                        function

    Returns:
        byref modified unsigned byte array

    """
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    offset = _24bmof(bmp, x1, y2)
    r = _xchrcnt(bmp)
    for buf in itercopyrect(
                    bmp, x1, y1, x2, y2):
        func(buf)
        BMPbitBLTput(bmp, offset, buf)
        offset += r


@func24bitonlyandentirerectinboundary
def _usebyref24btfn2reg(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        func: Callable, funcparam):
    """Apply byref func(funcparam) to a rectangular area
    in a 24-bit bitmap

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular area
        func          : user defined
                        function
        funcparam     : parameters of
                        the function

    Returns:
        byref modified unsigned byte array
    """
    x1, y1, x2, y2 = sortrecpoints(
                        x1, y1, x2, y2)
    offset = _24bmof(bmp, x1, y2)
    r = _xchrcnt(bmp)
    for buf in itercopyrect(
                    bmp, x1, y1, x2, y2):
        func(buf, funcparam)
        BMPbitBLTput(bmp,offset,buf)
        offset += r


@func24bitonlyandentirerectinboundary
def _use24bitfn2reg(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        func: Callable, funcparam):
    """Apply func(funcparam) to a rectangular area
    in a 24-bit bitmap

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular area
        func          : user defined
                        function
        funcparam     : parameters of
                        the function

    Returns:
        byref modified unsigned byte array
    """
    x1, y1, x2, y2 = sortrecpoints(
                        x1, y1, x2, y2)
    offset = _24bmof(bmp, x1, y2)
    r = _xchrcnt(bmp)
    for buf in itercopyrect(
                    bmp, x1, y1, x2, y2):
        BMPbitBLTput(bmp, offset,
            func(buf, funcparam))
        offset += r


@func24bitonlyandentirerectinboundary
def verticalbrightnessgradto24bitregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int,
        lumrange: list[int, int]):
    """Apply a vertical brightness gradient
    to a rectangular area in a 24-bit bitmap

    Args:
        bmp          :  unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region
        lumrange      : (byte:byte)
                        brightness
                        gradient

    Returns:
        byref modified unsigned byte array
    """
    x1, y1, x2, y2 = \
        sortrecpoints(x1, y1, x2, y2)
    offset = _24bmof(bmp, x1, y2)
    r = _xchrcnt(bmp)
    lum = lumrange[1]
    dlum = (lumrange[0] - lum) / (y2 - y1)
    for buf in itercopyrect(bmp, x1, y1, x2, y2):
        BMPbitBLTput(bmp, offset,
             applybrightnessadjtoBGRbuf(buf, lum))
        offset += r
        lum += dlum


@func24bitonlyandentirerectinboundary
def horizontalbrightnessgradto24bitregion(
        bmp: array, x1: int, y1: int,
                    x2: int, y2: int,
        lumrange: list[int, int]):
    """Apply a horizontal brightness gradient
    to a rectangular area in a 24-bit bitmap

    Args:
        bmp            : unsigned
                         byte array
                         with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region
        lumrange      : (byte: byte)
                        defines
                        the brightness
                        gradient

    Returns:
        byref modified unsigned byte array
    """
    r = _xchrcnt(bmp)
    c = _getBMoffhdfunc(bmp)
    x1, y1, x2, y2 = sortrecpoints(x1, y1, x2, y2)
    xlim = x2 + 1
    f =  applybrightnessadjtoBGRbuf
    l = lumrange[0]
    dl = delta(lumrange) / (x2 - x1)
    for x in range(x1, xlim):
        s = c(bmp, x, y2)
        e = c(bmp, x, y1)
        bmp[s: e - 2: r], bmp[s + 1: e - 1: r], bmp[s + 2: e: r] = \
        f(bmp[s: e - 2: r], l), f(bmp[s + 1: e -1 : r], l), f(bmp[s + 2: e: r], l)
        l += dl


@func24bitonlyandentirecircleinboundary
def magnifyNtimescircregion(bmp: array,
        x: int, y: int, r: int, n: int):
    """Magnify a circular region in a bitmap file by int n

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x,y)
                 and radius r
        n      : int magnification
                 factor

    Returns:
        byref modified unsigned byte array
    """
    nx = x * n
    ny = y * n
    b = resizeNtimesbigger(bmp, n)
    c = _BMoffsethd
    for (v0, v1) in itercirclepartlineedge(r):
        x1, x2 = mirror(x, v0)
        y1, y2 = mirror(y, v1)
        x3, x4 = mirror(nx, v0)
        y3, y4 = mirror(ny, v1)
        bmp[c(bmp, x1, y1): c(bmp, x2, y1)], bmp[c(bmp, x1, y2): c(bmp, x2, y2)] = \
        b[c(b, x3, y3): c(b, x4, y3)], b[c(b, x3, y4): c(b, x4, y4)]


@func24bitonlyandentirecircleinboundary
def pixelizenxncircregion(bmp: array,
        x: int, y: int, r: int, n: int):
    """Pixelize a circular region in a BMP by n

    Args:
        bmp    : unsigned byte array
                 with bmp format
        x, y, r: center (x, y)
                 and radius r
        n      : integer pixellation
                 dimension n by n

    Returns:
        byref modified unsigned byte array
    """
    b = pixelizenxn(bmp, n)
    c = _BMoffsethd
    for (v0, v1) in itercirclepartlineedge(r):
        x1, x2 = mirror(x, v0)
        y1, y2 = mirror(y, v1)
        bmp[c(bmp, x1, y1): c(bmp, x2, y1)], bmp[c(bmp, x1, y2): c(bmp,x2,y2)]= \
        b[c(b, x1, y1): c(b, x2, y1)], b[c(b, x1, y2): c(b, x2, y2)]


@func24bitonly
def resizeNtimessmaller(
        bmp: array, n: int) -> array:
    """Resize a whole image int n times smaller

    Args:
        bmp: unsigned byte array
             with bmp format
        n  : int resize factor

    Returns:
        byref modified unsigned byte array
    """
    bits = bmp[bmpcolorbits]
    nbmp = -1
    (mx, my) = getmaxxy(bmp)
    nx = mx // n
    ny = my // n
    nbmp = newBMP(nx, ny, bits)
    mx -= 1
    my -= 1
    ny -= 1
    offset = _24bmof(nbmp, 0, ny)
    r = _xchrcnt(nbmp)
    i = 1
    bufl = []
    y = 0
    for buf in itercopyrect(
                    bmp, 0, 0, mx, my):
        j = i % n
        bufl += [buf]
        if j == 0:
            BMPbitBLTput(nbmp, offset,
            resizesmaller24bitbuf(bufl))
            offset += r
            bufl = []
            y += 1
            if y > ny:
                break
        i += 1
    return nbmp


def pixelizenxn(bmp: array,
        n: int) -> array:
    """Pixelize a whole image with n by n areas
    in which colors are averaged

    Args:
        bmp: unsigned byte array
             with bmp format
        n  : size of pixel blur

    Returns:
        byref modified unsigned byte array
    """
    return resizeNtimesbigger(
           resizeNtimessmaller(bmp, n), n)


def adjustcolordicttopal(
        bmp: array, colordict: dict):
    """Adjust a color dictionary to match
    as closely as possible a bitmap palette

    Args:
        bmp      : unsigned byte array
                   with bmp format
        colordict: dictionary of colors

    Returns:
        byref modified dictionary of colors
    """
    if _getclrbits(bmp) < 24:
        for color in colordict:
            colordict[color] = \
                matchRGBtopal(
                    int2RGBarr(colordict[color]),
                    getallRGBpal(bmp))


def gammaadjto24bitregion(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        gamma: float):
    """Gamma correction to a rectangular area in a 24-bit BMP

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines the
                        rectangular
                        region
        gamma         : gamma correction

    Returns:
        byref modified unsigned byte array
    """
    _usebyref24btfn2reg(
        bmp, x1, y1, x2, y2,
        applygammaBGRbuf, gamma)


def gammaadjto24bitimage(
        bmp: array, gamma: float):
    """Gamma correction to an in-memory 24-bit BMP

    Args:
        bmp  : unsigned byte array
               with bmp format
        gamma: gamma correction

    Returns:
        byref modified unsigned byte array
    """
    gammaadjto24bitregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1, gamma)


@entirerectinboundary
def _cmpimglines(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int,
        func: Callable):
    offset = _BMoffset(bmp, x1, y2)
    r = _xchrcnt(bmp)
    oldbuf = []
    for buf in itercopyrect(
                    bmp, x1, y1, x2, y2):
        if oldbuf != []:
            BMPbitBLTput(bmp,
                offset, array('B', func(buf, oldbuf)))
            offset += r
        oldbuf = buf
    BMPbitBLTput(bmp, offset,
        array('B', func(buf, oldbuf)))


def outlineregion(bmp: array,
        x1: int, y1: int,
        x2: int, y2: int):
    """Outines a rectangular region in a BMP

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: defines
                        the rectangular
                        region

    Returns:
        byref modified unsigned byte array
    """
    _cmpimglines(bmp, x1, y1,
                      x2, y2, xorvect)


def outline(bmp: array):
    """Applies an Outline Filter

    Args:
        bmp: unsigned byte array
             with bmp format

    Returns:
        byref modified unsigned byte array
    """
    outlineregion(bmp, 0, 0,
        getmaxx(bmp) - 1,
        getmaxy(bmp) - 1)


def plotattractorcomplex(
        v: list[complex],
        x: int, y: int, bits: int):
    """Draws an Attractor from a list of Complex Numbers

    Args:
        v   : list of complex numbers
        x, y: int dimensions of bmp
        bits: int bit depth

    Returns:
        byref unsigned byte array
    """
    h = histogram2Dcomplex(v, x, y)
    fh = ravel2D(h)
    minval = min(fh)
    bmp = newBMP(x, y, bits)
    f = (getmaxcolors(bmp) - 1)/ (max(fh) - minval)
    for q in range(y):
        for p in range(x):
            plotxybit(bmp, p, q, int((h[q][p] - minval) * f))
    return bmp


def plotattractor(
        v: list[list[float, float]],
        x: int, y: int, bits: int):
    """Draws an Attractor from a list of x, y pairs

    Args:
        v   : list of [x, y] coordinate pairs
        x, y: int dimensions of bmp
        bits: int bit depth

    Returns:
        byref unsigned byte array
    """
    h = histogram2D(v, x, y)
    fh = ravel2D(h)
    minval = min(fh)
    bmp = newBMP(x, y, bits)
    f = (getmaxcolors(bmp) - 1)/ (max(fh) - minval)
    for q in range(y):
        for p in range(x):
            plotxybit(bmp, p, q, int((h[q][p] - minval) * f))
    return bmp


def ikedaattractor(
        x: int, y: int, bits: int,
        a: float, b: float,
        k: float, p: float, n: int):
    """Draws an Ikeda Attractor

    Args:
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, k, p: float coefficients
        n: number of terms to compute

    Returns:
        byref unsigned byte array
    """
    return plotattractorcomplex(ikedaattractorlist(a, b, k, p, n), x, y, bits)


def nattractor(
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, d: float, n: int):
    """Draws a N Attractor

    Args:
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, c, d: float coefficients
        n: number of terms to compute

    Returns:
        byref unsigned byte array
    """
    return plotattractor(nattractorlist(a, b, c, d, n), x, y, bits)


def peterdejongattractor(
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, d: float, n: int):
    """Draws a Peter de Jong Attractor

    Args:
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, c, d: float coefficients
        n: number of terms to compute

    Returns:
        byref unsigned byte array
    """
    return plotattractor(peterdejongattractorlist(a, b, c, d, n), x, y, bits)


def duffingattractor(
        x: int, y: int, bits: int,
        a: float, b: float,
        w: float, dt: float, n: int):
    """Draws a Peter de Jong Attractor

    Args:
        x, y    : int dimensions of bmp
        bits    : int bit depth
        a, b, w : float coefficients
        dt      : time interval
        n: number of terms to compute

    Returns:
        byref unsigned byte array
    """
    return plotattractor(duffingattractorlist(0, 0, 0, dt, a, b, w, n), x, y, bits)


def hopalongattractor(
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, n: int):
    """Draws a Hopalong Attractor

    Args:
        x, y   : int dimensions of bmp
        bits   : int bit depth
        a, b, c: float coefficients
        n: number of terms to compute

    Returns:
        byref unsigned byte array
    """
    return plotattractor(hopalongattractorlist(a, b, c, n), x, y, bits)


def symmetriciconattractor(
        x: int, y: int, bits: int,
        ox: float, oy: float,
        a: float, b : float, g: float,
        o: float, l: float,
        d: int, n: int):
    """Draws a Symmetric Icon Attractor

    Args:
        x, y   : int dimensions of bmp
        bits   : int bit depth
        ox, oy: float starting coordinates
        a, b, g, o, l: float coefficients
        d: int degree
        n: number of terms to compute

    Returns:
        byref unsigned byte array
    """
    return plotattractorcomplex(symmetriciconattractorlist(ox, oy, a, b, g, o, l, d, n), x, y, bits)


def cliffordattractor(
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, d: float, n: int):
    """Draws a Clifford Attractor

    Args:
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, c, d: float coefficients
        n: number of terms to compute

    Returns:
        byref unsigned byte array
    """
    return plotattractor(cliffordattractorlist(a, b, c, d, n), x, y, bits)


def fractaldreamattractor(
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, d: float, n: int):
    """Draws a Fractal Dream Attractor

    Args:
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, c, d: float coefficients
        n: number of terms to compute

        a and b are floating point values
        between -3 and +3

        c and d are floating point values
        between -0.5 and +1.5

    Returns:
        byref unsigned byte array
    """
    return plotattractor(fractaldreamattractorlist(a, b, c, d, n), x, y, bits)


def gumowskimiraattractor(
        x: int, y: int, bits: int,
        ox: float, oy: float,
        a: float, b: float, n: int):
    """Draws a Gumowski-Mira Attractor

    Args:
        x, y: int dimensions of bmp
        bits: int bit depth
        ox, oy: starting coordinates
        a, b: float coefficients
        n: number of terms to compute

    ox and oy -> any floating point
    value between -20 and +20

    The a and b -> any floating point
    value between -1 and +1.

    Returns:
        byref unsigned byte array
    """
    return plotattractor(gumowskimiraattractorlist(ox, oy, a, b, n), x, y, bits)


@functimer
def saveikedaattractor2file(
        file: str,
        x: int, y: int, bits: int,
        a: float, b: float,
        k: float, p: float, n: int):
    """Draws an Ikeda Attractor to file

    Args:
        file      : full path to new file
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, k, p: float coefficients
        n: number of terms to compute

    Returns:
        a bitmap file
    """
    saveBMP(file, ikedaattractor(x, y, bits, a, b, k, p, n))


@functimer
def savegumowskimiraattractor2file(
        file: str,
        x: int, y: int, bits: int,
        ox: float, oy: float,
        a: float, b: float, n: int):
    """Draws a Gumowski-Mira attractor to file

    Args:
        file: full path to new file
        x, y: int dimensions of bmp
        bits: int bit depth
        ox, oy: starting coordinates
        a, b: float coefficients
        n: number of terms to compute

    ox and oy -> any floating point
    value between -20 and +20

    The a and b -> any floating point
    value between -1 and +1.

    Returns:
        a bitmap file
    """
    saveBMP(file, gumowskimiraattractor(x, y, bits, ox, oy, a, b, n))


@functimer
def savenattractor2file(
        file: str,
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, d: float, n: int):
    """Draws a N Attractor to file

    Args:
        file      : full path to new file
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, c, d: float coefficients
        n: number of terms to compute

    Returns:
        a bitmap file
    """
    saveBMP(file, nattractor(x, y, bits, a, b, c, d, n))


@functimer
def savepeterdejongattractor2file(
        file: str,
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, d: float, n: int):
    """Draws a Peter de Jong Attractor to file

    Args:
        file      : full path to new file
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, c, d: float coefficients
        n: number of terms to compute

    Returns:
        a bitmap file
    """
    saveBMP(file, peterdejongattractor(x, y, bits, a, b, c, d, n))


@functimer
def saveduffingattractor2file(
        file: str,
        x: int, y: int, bits: int,
        a: float, b: float,
        w: float, dt: float, n: int):
    """Draws a Duffing Attractor to file

    Args:
        file   : full path to new file
        x, y   : int dimensions of bmp
        bits   : int bit depth
        a, b, c: float coefficients
        dt     : time interval
        n: number of terms to compute


    Returns:
        a bitmap file
    """
    saveBMP(file, duffingattractor(x, y, bits, a, b, w, dt, n))


@functimer
def savehopalongattractor2file(
        file: str,
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, n: int):
    """Draws a Hopalong Attractor to file

    Args:
        x, y   : int dimensions of bmp
        bits   : int bit depth
        a, b, c: float coefficients
        n: number of terms to compute

    Returns:
        a bitmap file
    """
    saveBMP(file, hopalongattractor(x, y, bits, a, b, c, n))


@functimer
def savesymmetriciconattractor2file(
        file: str,
        x: int, y: int, bits: int,
        ox: float, oy: float,
        a: float, b : float, g: float,
        o: float, l: float,
        d: int, n: int):
    """Draws a Symmetric Icon Attractor to file

    Args:
        x, y   : int dimensions of bmp
        bits   : int bit depth
        ox, oy: float starting coordinates
        a, b, g, o, l: float coefficients
        d: int degree
        n: number of terms to compute

    Returns:
        a bitmap file
    """
    saveBMP(file, symmetriciconattractor(x, y, bits, ox, oy, a, b, g, o, l, d, n))


@functimer
def savecliffordattractor2file(
        file: str,
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, d: float, n: int):
    """Draws a Clifford Attractor to file

    Args:
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, c, d: float coefficients
        n: number of terms to compute

    Returns:
        a bitmap file
    """
    saveBMP(file, cliffordattractor(x, y, bits, a, b, c, d, n))


@functimer
def savefractaldreamattractor2file(
        file: str,
        x: int, y: int, bits: int,
        a: float, b: float,
        c: float, d: float, n: int):
    """Draws a Fractal Dream Attractor to file

    Args:
        x, y      : int dimensions of bmp
        bits      : int bit depth
        a, b, c, d: float coefficients
        n: number of terms to compute

        a and b are floating point values
        between -3 and +3

        c and d are floating point values
        between -0.5 and +1.5

    Returns:
        a bitmap file
    """
    saveBMP(file, fractaldreamattractor(x, y, bits, a, b, c, d, n))


@intcircleparam
def sphere(bmp: array,
        x: int, y: int, r: int,
        rgbfactors: list[float, float, float],
        fn: Callable = lerp):
    """Draws a Rendered Sphere

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : center of sphere
                    in the image
        r         : radius of sphere
                    in pixels
        rgbfactors: (r,g,b) r, g and b
                    values range from
                    0.0 to 1.0
                        or
                    (int, int)
                    packed RGB
        fn: optional gradient function

    Returns:
        byref modified unsigned byte array
    """
    cp = len(rgbfactors)
    if cp == 3:
        gradcircle(bmp, x, y, r,
            [255, 0], rgbfactors, fn)
    if cp == 2:
        graddichromiccircle(bmp, x, y, r,
        rgbfactors[0], rgbfactors[1], fn)


@intcircleparam
def thickencirclearea(bmp: array,
        x: int, y: int, r: int,
        rgbfactors: list[float, float, float]):
    """Encircle area with a gradient

    Args:
        bmp       : unsigned byte array
                    with bmp format
        x, y      : center of circle
        r         : radius of circle
        rgbfactors: (r,g,b) r, g and b
                    values are 0 to 1
                    unsigned floats

    Returns:
        byref modified unsigned byte array
    """
    gradthickcircle(bmp, x, y, r,
        8, [255,0], rgbfactors)


@entirerectinboundary
def fern(bmp: array, x1: int, y1: int,
         x2: int, y2: int, color: int):
    """Draws an IFS fern fractal

    Args:
        bmp           : unsigned
                        byte array
                        with bmp format
        x1, y1, x2, y2: rectangular
                        area to draw
                        the fern in
        color         : color of the
                        fern fractal

    Returns:
        byref modified unsigned byte array
    """
    y = abs(y2 - y1) // 10
    IFS(bmp, getIFSparams()['fern'],
        x1, y1, x2, y2, y, y,
        abs(x2 - x1) // 2, 0,
        color, 100000)


@checklink
def _usebyreffnwithpar2regnsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        func: Callable, *param):
    """Apply a byref function with parameters
    to a rectangular area and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         area
        func           : user defined
                         function
        *param  : function parameters

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    func(bmp, x1, y1, x2, y2, *param)
    saveBMP(NewBMPfile, bmp)
    print(sysmsg['savedareafunc'] %
    (func.__name__, x1, y1, x2, y2,
     ExistingBMPfile, NewBMPfile))


@checklink
def _use24btbyrefclrfn2regnsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        func: Callable, funcparam):
    """Apply a 24-bit byref
    color modification function
    to a rectangular area and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2:  the rectangular
                         area
        func           : user defined
                         function
        funcparam      : function
                         parameters

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        print(sysmsg['not24bit'])
    else:
        func(bmp, x1, y1, x2, y2,
        funcparam)
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savedareafunc'] %
        (func.__name__, x1, y1, x2, y2,
          ExistingBMPfile, NewBMPfile))


@checklink
def _usebyref24btfn2regnsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        func: Callable):
    """Apply a 24-bit byref function
    with no parameters to a
    rectangular area and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        func           : user defined
                         function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        print(sysmsg['not24bit'])
    else:
        func(bmp, x1, y1, x2, y2)
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savedareafunc'] %
        (func.__name__, x1, y1, x2, y2,
        ExistingBMPfile, NewBMPfile))


@checklink
def _usebyreffn2regnsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        func: Callable):
    """Apply a byref function
    with no parameters to a
    rectangular region and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        func           : user defined
                         function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    func(bmp, x1, y1, x2, y2)
    saveBMP(NewBMPfile, bmp)
    print(sysmsg['savedareafunc'] %
     (func.__name__, x1, y1, x2, y2,
      ExistingBMPfile, NewBMPfile))


@checklink
def _usefn2regsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        func: Callable):
    """Apply a function with no parameters to a
    rectangular area and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        func           : user defined
                         function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    bmp = func(bmp, x1, y1, x2, y2)
    if bmp != None:
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savedareafunc'] %
        (func.__name__, x1, y1, x2, y2,
          ExistingBMPfile, NewBMPfile))


@checklink
def _usebyreffnsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable):
    """Apply a by-ref function with no parameters and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        func           : user defined
                         function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    func(bmp)
    saveBMP(NewBMPfile, bmp)
    print(sysmsg['savefunc'] %
        (func.__name__ ,
        ExistingBMPfile, NewBMPfile))


@checklink
def _usebyreffnwithparnsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable,
        funcparam: any):
    """Apply a by-ref function with parameter and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        func           : user defined
                         function
        funcparam      : function parameter

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    func(bmp, funcparam)
    saveBMP(NewBMPfile, bmp)
    print(sysmsg['savesingleparamfunc'] %
    (func.__name__, str(funcparam),
    ExistingBMPfile, NewBMPfile))


@checklink
def _usefnsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable):
    """Apply a function with no parameters and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        func           : user defined
                         function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    saveBMP(NewBMPfile, func(bmp))
    print(sysmsg['savefunc'] %
    (func.__name__ ,
    ExistingBMPfile, NewBMPfile))


@checklink
def _usefndirsv(
        ExistingBMPfile: str,
        NewDir: str,
        func: Callable,
        overwritedir: bool = False):
    """Apply a function with no parameters and save
    images to a directory

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        func           : user defined
                         function

    Returns:
        new bitmap files in a new directory
    """
    obmp = loadBMP(ExistingBMPfile)
    cwd = getcwd()
    try:
        mkdir(NewDir)
    except FileExistsError:
        if not overwritedir:
            raise FileExistsError(sysmsg['direxist'])
    for (i, bmp) in enumerate(func(obmp)):
        NewBMPfile = f'{cwd}/{NewDir}/{i}.bmp'
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savefunc'] %
        (func.__name__ ,
        ExistingBMPfile, NewBMPfile))


@checklink
def _use24btfnwithparnsv(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable, funcparam):
    """Apply a 24-bit only function with parameters and save

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        func           : user defined
                         function
        funcparam      : parameters of
                         the function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        print(sysmsg['not24bit'])
    else:
        saveBMP(NewBMPfile,
            func(bmp, funcparam))
        print(sysmsg['savesingleparamfunc'] %
        (func.__name__, str(funcparam),
          ExistingBMPfile, NewBMPfile))


@checklink
def _usefn2circreg(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable,
        x: int, y: int, r: int):
    """Apply a user provided function (no parameters)
    to a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r
        func           : user defined
                         function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    func(bmp, x, y, r)
    saveBMP(NewBMPfile, bmp)
    print(sysmsg['savecircfunc'] %
         (func.__name__ , x, y, r,
         ExistingBMPfile, NewBMPfile))


@checklink
def _usefnwithpar2circreg(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable,
        x: int, y: int, r: int,
        funcparam):
    """Apply a user provided function with parameters
    to a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x,y)
                         and radius r
        func           : user defined
                         function
        funcparam      : parameters of
                         the function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    func(bmp, x, y, r, funcparam)
    saveBMP(NewBMPfile, bmp)
    print(sysmsg['savecircfuncwithparam'] %
            (func.__name__ , x, y, r,
            funcparam, ExistingBMPfile,
            NewBMPfile))


@checklink
def _use24btclrfntocircregion(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable,
        x: int, y: int, r: int):
    """Apply a no parameter color adjustment function
    to a circular area (24-bit only)

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r
        func           : user defined
                         function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        print(sysmsg['not24bit'])
    else:
        func(bmp, x, y, r)
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savecircfunc'] %
            (func.__name__, x, y, r,
            ExistingBMPfile, NewBMPfile))


@checklink
def _use24btclrfnwithpar2circreg(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable,
        x: int, y: int, r: int,
        *funcparam):
    """Apply a user provided color adjustment function
    to a circular area (24-bit only)

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x,y)
                         and radius r
        func           : user defined
                         function
        *funcparam     : parameters of
                         the function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        print(sysmsg['not24bit'])
    else:
        func(bmp, x, y, r, *funcparam)
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savecircfuncwithparam'] %
            (func.__name__, x, y, r, *funcparam,
            ExistingBMPfile, NewBMPfile))


def _use24btclrfnwithmultpar2circreg(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable,
        x: int, y: int, r: int,
        funcparam,
        *otherparam):
    """Apply a user provided color adjustment function
    to a circular area (24-bit only)

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x,y)
                         and radius r
        func           : user defined
                         function
        funcparam      : parameters of
                         the function
        *otherparam

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        print(sysmsg['not24bit'])
    else:
        func(bmp, x, y, r, funcparam, *otherparam)
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savecircfuncwithparam'] %
            (func.__name__, x, y, r, funcparam,
            ExistingBMPfile, NewBMPfile))


@checklink
def _useclradjfn(ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable, funcparam):
    """Apply a user provided color adjustment function
    to an existing bitmap

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        func           : user defined
                         function
        funcparam      : parameters of
                         the function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        setbmppal(bmp,[func(c, funcparam)
              for c in getallRGBpal(bmp)])
    else:
        if func.__name__ == 'colorfilter':
            colorfilterto24bitimage(
                bmp, funcparam)
        elif func.__name__ == 'brightnessadjust':
            brightnesseadjto24bitimage(
                bmp, funcparam)
        elif func.__name__ == 'gammacorrect':
            gammaadjto24bitimage(
                bmp, funcparam)
        elif func.__name__ == 'thresholdadjust':
            thresholdadjto24bitimage(
                bmp, funcparam)
        else:
            for v in iterimageRGB(bmp,
                        sysmsg['coloradj'], '*',
                        sysmsg['done']):
                plotRGBxybitvec(
                    bmp, v[0],
                    func(v[1],funcparam))
    saveBMP(NewBMPfile,bmp)
    print(sysmsg['savesingleparamfunc'] %
         (func.__name__, str(funcparam),
          ExistingBMPfile, NewBMPfile))


@checklink
def _use24btclrfn(ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable, funcparam):
    """Apply a user provided color adjustment function
    to a 24-bit bitmap

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save the
                         changes in
        func           : user defined
                         function
        funcparam      : parameters of
                         the function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        print(sysmsg['not24bit'])
    else:
        func(bmp, funcparam)
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savesingleparamfunc'] %
        (func.__name__, str(funcparam),
          ExistingBMPfile, NewBMPfile))


@checklink
def _usenoparclradjfn(
        ExistingBMPfile: str,
        NewBMPfile: str,
        func: Callable):
    """Apply a user provided no parameter color
        adjustment function to an existing bitmap

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        func           : user defined
                         function

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    if bmp[bmpcolorbits] != 24:
        setbmppal(bmp,
            [func(c) for c in getallRGBpal(bmp)])
    else:
        if func.__name__ == 'monochrome':
            monofilterto24bitimage(bmp)
        else:
            for v in iterimageRGB(
                        bmp,
                        sysmsg['coloradj'],
                        '*',
                        sysmsg['done']):
                plotRGBxybitvec(bmp, v[0],
                                func(v[1]))
    saveBMP(NewBMPfile, bmp)
    print(sysmsg['savenoparamfunc'] %
         (func.__name__, ExistingBMPfile,
                             NewBMPfile))


@checklink
def cropBMPandsave(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Crops and saves a rectangular area to a BMP

    Args:
        ExistingBMPfile: Whole path
                         to an
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x1, y1, x2, y2 : the rectagular
                         region

    Returns:
        new bitmap file
    """
    _usefn2regsv(ExistingBMPfile,
        NewBMPfile, x1, y1, x2, y2,
        crop)


@checklink
def cropBMPandsaveusingrectbnd(
        ExistingBMPfile: str,
        NewBMPfile: str,
        rectbnd: list):
    """Crops and saves a rectangular area to a BMP

    Args:
        ExistingBMPfile: Whole path
                         to existing file
        NewBMPfile     : New file to
                         save changes in
        rectbnd        : list defining
                         a rectangular
                         region
                         [(x1, y1),
                          (x2, y2),
                          (x3, y3),
                          (x4, y4)]

    Returns:
        new bitmap file
    """
    _usefn2regsv(ExistingBMPfile,
        NewBMPfile, rectbnd[0][0],
        rectbnd[0][1],
        rectbnd[1][0] + 1,
        rectbnd[1][1] + 1, crop)


@checklinks
def imagecomp(
        inputfile1: str,
        inputfile2: str,
        diff_file: str,
        func: Callable):
    """Perform a bitwise comparison of two bitmap files
    with the same x and y dimensions and bit depth
    using a user defined bitwise comparator function

    Args:
        Inputfile1: path to bmp file
        Inputfile2: path to bmp file
        diff_file : New file to save
                    comparison in
        func      : User provided
                    bitwise function

    Returns:
        new bitmap file
    """
    bmp1 = loadBMP(inputfile1)
    bmp2 = loadBMP(inputfile2)
    s1 = getmaxxyandbits(bmp1)
    s2 = getmaxxyandbits(bmp2)
    if s1 != s2:
        print(sysmsg['cantcomparefiles'] % (s1, s2))
    else:
        bits = s1[1]
        nbmp=CopyBMPxydim2newBMP(bmp1, bits)
        if bits < 24:
            pal1 = getallRGBpal(bmp1)
            pal2 = getallRGBpal(bmp2)
            if pal1 != pal2:
                print(sysmsg['diffpal'])
            copyRGBpal(bmp1, nbmp)
        setBMPimgbytes(nbmp, array('B',
             func(getBMPimgbytes(bmp1),
                  getBMPimgbytes(bmp2))))
        saveBMP(diff_file, nbmp)
        print(sysmsg['savdifffile'] % diff_file)


@checklink
@functimer
def reduce24bitimagebits(
        Existing24BMPfile: str,
        NewBMPfile: str, newbits: int,
        similaritythreshold: float,
        usemonopal: bool,
        RGBfactors: list[float, float, float] = None):
    """Reduce bits used to encode color in a 24-bit BMP

    Args:
        ExistingBMPfile    : Whole path
                             to existing file
        NewBMPfile         : New file to
                             save changes in
        newbits            : can be 1, 4
                             or 8 bits
        similaritythreshold: how close can
                             a color be to
                             another color
        usemonopal         : True -> image
                             will be mono
        RGBfactors         : (r: float,
                              b: float,
                              g: float)
                             values range
                             from 0 to 1
                             used only if
                             usemonopal
                             is True
    Returns:
        new bitmap file

    """
    sbmp = loadBMP(Existing24BMPfile)
    if sbmp[bmpcolorbits] != 24:
        print(sysmsg['not24bit'])
    else:
        bmp = CopyBMPxydim2newBMP(
                sbmp, newbits)
        if newbits > 1:
            if usemonopal:
                newpal = setBMP2monochrome(
                            bmp, RGBfactors)
            else:
                newpal = setnewpalfromsourcebmp(
                            sbmp, bmp,
                            similaritythreshold)
        for v in iterimageRGB(sbmp,
                    sysmsg['colorquant'],
                    '*',
                    sysmsg['done']):
                if newbits == 1:
                    c = probplotRGBto1bit(v[1], 2)
                else:
                    c = matchRGBtopal(v[1], newpal)
                intplotvecxypoint(bmp, v[0], c)
        saveBMP(NewBMPfile, bmp)
        print(sysmsg['savemod'] % (
                 Existing24BMPfile,
                        NewBMPfile))


@checklink
@functimer
def imgregionbyRGB2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        edgeradius: int,
        edgecolor: int,
        rgb: list[int, int, int],
        similaritythreshold: float,
        showedgeonly:bool):
    """Gets an image region by rgb in a bitmap file

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        edgeradius     : radius of pen
                         for the edge
        edgecolor      : color of edge
        rgb            : (r: byte,
                          g: byte,
                          b: byte)
                         color to select
        similaritythreshold: how close
                             the color
                             is to rgb
                             before
                             selection
        showedgeonly: True-> only edges
                 False-> edge and image

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    edge = getimageregionbyRGB(
            bmp, rgb, similaritythreshold)
    if showedgeonly:
        bmp = copyBMPhdr(bmp)
    plotxypointlist(bmp,
        edge, edgeradius, edgecolor)
    saveBMP(NewBMPfile,bmp)
    print(sysmsg['savemod'] % (
         ExistingBMPfile, NewBMPfile))


@functimer
def invertbits2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Inverts the bits in a bmp file

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, invertimagebits)


@functimer
def flipvertical2file(
    ExistingBMPfile: str,
    NewBMPfile: str):
    """Flips a bitmap file vertically

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, flipvertical)


@functimer
def mirrortop2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Mirrors the top-half of a bitmap

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, mirrortop)


@functimer
def mirrortopleft2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Mirrors the top-left of a bitmap

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, mirrortopleft)


@functimer
def mirrortopright2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Mirrors the top-right of a
        bitmap file

    Args:
        ExistingBMPfile: Whole path
                         to existing
                         file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, mirrortopright)


@functimer
def mirrorbottomleft2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Mirrors the bottom-left of a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, mirrorbottomleft)


@functimer
def mirrorbottomright2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Mirrors the bottom-right of a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, mirrorbottomright)


@functimer
def mirrorbottom2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Mirrors the bottom half of a BMP

    Args:
        ExistingBMPfile: Whole path
                         to an existing
                         file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, mirrorbottom)


@functimer
def mirrorleft2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Mirrors the left-half of a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, mirrorleft)


@functimer
def mirrorright2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Mirrors the right half of a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, mirrorright)


@functimer
def fliphorizontal2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Flips horizontally the image in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
        NewBMPfile, fliphorizontal)


@functimer
def flipXY2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Flips the x and y coordinates to rotate by 90 degrees

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usefnsv(ExistingBMPfile,
        NewBMPfile, flipXY)


@functimer
def palshift2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        shift: int = 1):
    """Shifts the palette of a 4 or 8 bit bitmap

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        shift          : int shift value
                         (default = 1)

    Returns:
        new bitmap file
    """
    _usebyreffnwithparnsv(ExistingBMPfile,
        NewBMPfile, palshift, shift)


@functimer
def flipverticalregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Flips vertically a rectangular area in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        flipverticalregion)


@functimer
def fliphorizontalregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Flips horizontally a rectangular area in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        fliphorizontalregion)


@functimer
def mirrorleftinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirrors the left-half region in a rectangular area in a bmp

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        mirrorleftinregion)


@functimer
def mirrorrightinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirrors the right-half region in a rectangular area in a bmp

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        mirrorrightinregion)


@functimer
def mirrortopinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirrors the top-half region in a rectangular area in a bmp

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         region

    Returns:
        new bitmap file

    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        mirrortopinregion)


@functimer
def mirrorbottominregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirrors the bottom-half region in a rectangular area in a bmp

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        mirrorbottominregion)


@functimer
def mirrortopleftinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirrors the top-left region in a rectangular area in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        mirrortopleftinregion)


@functimer
def mirrortoprightinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirrors the top-right region in a rectangular area in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        mirrortoprightinregion)


@functimer
def mirrorbottomleftinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirrors the bottom-left region in a rectangular area in a bmp

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        mirrorbottomleftinregion)


@functimer
def mirrorbottomrightinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Mirrors the bottom-right region in a rectangular area in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2  :the rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        mirrorbottomrightinregion)


@functimer
def invertregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Inverts the bits in a rectangular area in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : the rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2, invertregion)


@functimer
def autocropimg2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        similaritythreshold: float):
    """Perform an auto crop to the image in a bitmap file

    Args:
        ExistingBMPfile    : Whole path to
                             existing file
        NewBMPfile         : New file to
                             save changes in
        similaritythreshold: used to tune
                             autocrop

    Returns:
        new bitmap file
    """
    bmp = loadBMP(ExistingBMPfile)
    cropBMPandsaveusingrectbnd(
        ExistingBMPfile, NewBMPfile,
        rectboundarycoords(
            getimagedgevert(bmp,
                similaritythreshold)))


@functimer
def adjustbrightness2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        percentadj: float):
    """Apply a brightness adjustment to the image in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        percentadj     : can be a
                         positive or
                         negative value
                         (signed float)

    Returns:
        new bitmap file
    """
    _useclradjfn(ExistingBMPfile,
        NewBMPfile, brightnessadjust,
        percentadj)


@functimer
def thresholdadjust2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        lumrange: list[int, int]):
    """Apply a threshold adjustment

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        lumrange       : (byte:byte)
                         threshold
                         to apply

    Returns:
        new bitmap file
    """
    _useclradjfn(ExistingBMPfile,
        NewBMPfile, thresholdadjust,
        lumrange)


@functimer
def adjustbrightnessinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        percentadj: float):
    """Brightness adjustment to rectangular area in a 24-bit BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        percentadj     : can be a
                         positive or
                         negative
                         adjustment
                        (signed float)

    Returns:
        new bitmap file
    """
    _use24btbyrefclrfn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        brightnesseadjto24bitregion,
        percentadj)


@functimer
def adjustthresholdinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        lumrange: list):
    """Threshold adjust in a rectangular area in a 24-bit BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        lumrange       : (byte:byte)
                         threshold to
                         apply

    Returns:
        new bitmap file
    """
    _use24btbyrefclrfn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        thresholdadjto24bitregion,
        lumrange)


@functimer
def colorfilter2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        rgbfactors: list[float, float, float]):
    """Applies Color Filter rgbfactors to a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        rgbfactors     : (r,g,b) r, g, b
                         values range
                         from 0.0 to 1.0

    Returns:
        new bitmap file
    """
    _useclradjfn(ExistingBMPfile,
        NewBMPfile, colorfilter,
        rgbfactors)


@functimer
def monochrome2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Applies a monochrome filter to a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usenoparclradjfn(ExistingBMPfile,
               NewBMPfile, monochrome)


@functimer
def colorfilterinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        rgbfactors: list[float, float, float]):
    """Color Filter to a Rectangular Area in a 24-bit BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        rgbfactors     : (r,g,b)
                         values range
                         from
                         0.0 to 1.0

    Returns:
        new bitmap file
    """
    _use24btbyrefclrfn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        colorfilterto24bitregion,
        rgbfactors)


@functimer
def monofilterinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Monochrome filter to rectangular area in a 24-bit BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyref24btfn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        monofilterto24bitregion)


@functimer
def pixelizenxntofile(
        ExistingBMPfile: str,
        NewBMPfile: str, n: int):
    """Pixellate a bitmap file with n by n pixel areas

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _use24btfnwithparnsv(
        ExistingBMPfile, NewBMPfile,
        pixelizenxn, n)


@functimer
def resizeNtimessmaller2file(
        ExistingBMPfile: str,
        NewBMPfile: str, n: int):
    """Resize a bitmap file n times smaller

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _use24btfnwithparnsv(
        ExistingBMPfile, NewBMPfile,
        resizeNtimessmaller, n)


@functimer
def resizeNtimesbigger2file(
        ExistingBMPfile: str,
        NewBMPfile: str, n: int):
    """Resize a bitmap file n times bigger

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _use24btfnwithparnsv(
        ExistingBMPfile, NewBMPfile,
        resizeNtimesbigger, n)


@functimer
def upgradeto24bitimage2file(
        ExistingBMPfile: str,
        NewBMPfile: str):
    """Upgrades a bitmap file to 24-bits

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usefnsv(ExistingBMPfile,
      NewBMPfile, upgradeto24bitimage)


@functimer
def gammaadj2file(
        ExistingBMPfile: str,
        NewBMPfile: str, gamma: float):
    """Applies a Gamma Correction

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        gamma          : gamma
                         correction

    Returns:
        new bitmap file
    """
    _useclradjfn(ExistingBMPfile,
      NewBMPfile, gammacorrect, gamma)


@functimer
def gammaadjtoregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        gamma: float):
    """Gamma Correction to a Rectangular Region

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        gamma          : gamma
                         correction

    Returns:
        new bitmap file
    """
    _use24btbyrefclrfn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        gammaadjto24bitregion, gamma)


@functimer
def eraseeverynthhorilineinregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int, n: int):
    """Erase every nth line in a rectangular region

    Args:
        ExistingBMPfile: Whole path to
                         an existing
                         file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        n              : erase every
                         nth line

    Returns:
        new bitmap file
    """
    _usebyreffnwithpar2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        eraseeverynthhorilineinregion, n)


@functimer
def rectangle2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        color: int):
    """Draws a Rectangle

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular region
        color          : color of rectangle

    Returns:
        new bitmap file
    """
    _usebyreffnwithpar2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2, rectangle,
        color)


@functimer
def filledrect2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        color: int):
    """Draws a Filled Rectangle

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular region
        color          : color of rectangle

    Returns:
        new bitmap file
    """
    _usebyreffnwithpar2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2, filledrect,
        color)


@functimer
def filledgradrect2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        lumrange: list[int, int],
        rgbfactors: list[float, float, float],
        direction: int,
        fn: Callable = lerp):
    """Draws a Filled Gradient Rectangle to file

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular region
        lumrange       : [byte, byte]
                         luminosity range
        rgbfactors     : [r, g, b] values
                         range from
                         0.0 to 1.0
        direction      : 0 - vertical
                         1 - horizontal
        fn             : optional gradient function

    Returns:
        new bitmap file
    """
    _usebyreffnwithpar2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2, filledgradrect,
        lumrange, rgbfactors, direction, fn)


@functimer
def filleddichromaticgradrect2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        c1: int, c2: int,
        direction: int,
        fn: Callable = lerp):
    """Draws a Filled Dichromatic Gradient Rectangle to file

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular region
        c1, c2         : int packed RGB
        direction      : 0 - vertical
                         1 - horizontal

    Returns:
        new bitmap file
    """
    _usebyreffnwithpar2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2, filleddichromaticgradrect,
        c1, c2, direction, fn)


@functimer
def fern2file(ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        color: int):
    """Fern Fractal in a bounding rectangle

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular region
        color          : color of the
                         fern fractal

    Returns:
        new bitmap file
    """
    _usebyreffnwithpar2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2, fern, color)


@functimer
def eraseeverynthhoriline2file(
        ExistingBMPfile : str,
        NewBMPfile: str, n: int):
    """Erase every nth line

    Args:
        ExistingBMPfile: Whole path
                         to existing file
        NewBMPfile     : New file to
                         save changes in
        n              : erase every
                         nth line

    Returns:
        new bitmap file
    """
    _usebyreffnwithparnsv(
        ExistingBMPfile, NewBMPfile,
        eraseeverynthhorizontalline, n)


@functimer
def outlineregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int):
    """Outline filter to rectangular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x1, y1, x2, y2 : defines the
                         rectangular
                         region

    Returns:
        new bitmap file
    """
    _usebyreffn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2, outlineregion)


@functimer
def outline2file(ExistingBMPfile: str,
                 NewBMPfile: str):
    """Applies an outline filter

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in

    Returns:
        new bitmap file
    """
    _usebyreffnsv(ExistingBMPfile,
             NewBMPfile, outline)


@functimer
def imagediff(inputfile1: str,
     inputfile2: str, diff_file: str):
    """Compares 2 files and saves diff to a bitmap file

    Args:
        inputfile1: Whole paths
        inputfile2 to existing files

        diff_file: New file
                   to store diff

    Returns:
        new bitmap file
    """
    imagecomp(inputfile1, inputfile2,
        diff_file, xorvect)


@functimer
def showsimilarparts(
        inputfile1: str,
        inputfile2: str,
        diff_file: str):
    """Compares 2 files and saves the similar parts to a BMP

    Args:
        inputfile1: Whole paths
        inputfile2  to existing files

        diff_file : New file to
                    store similar parts

    Returns:
        new bitmap file
    """
    imagecomp(inputfile1, inputfile2,
        diff_file, andvect)


@functimer
def monochromecircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Monochrome filter to a circular region

    Args:
        ExistingBMPfile: Whole path
                         to an existing
                         file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _use24btclrfntocircregion(
        ExistingBMPfile, NewBMPfile,
        monocircle, x, y, r)


@functimer
def invertbitsincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Inverts bits in a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        invertbitsincircregion, x, y, r)


@functimer
def colorfiltercircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        rgbfactors: list[float, float, float]):
    """Applies a color filter to a circular area in a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
        rgbfactors     : (r, g, b) values
                         range from 0 to 1

    Returns:
        new bitmap file
    """
    _use24btclrfnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        colorfiltercircregion,
        x, y, r, rgbfactors)


@functimer
def thresholdadjcircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        lumrange: list[int, int]):
    """Threshold adjustment to a circular region in a 24-bit BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
        lumrange       : (byte:byte)
                         threshold range

    Returns:
        new bitmap file
    """
    _use24btclrfnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        thresholdadjcircregion,
        x, y, r, lumrange)


@functimer
def gammacorrectcircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        gamma: float):
    """Gamma Adjust to a circular area in a BMP

    Args:
        ExistingBMPfile: Whole path
                         to existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
                         of a circular
                         region
        gamma          : gamma
                         adjustment

    Returns:
        new bitmap file
    """
    _use24btclrfnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        gammacorrectcircregion,
        x, y, r, gamma)


@functimer
def sphere2file(ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        rgbfactors: list[float, float, float],
        fn: Callable = lerp):
    """Renders a sphere

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r
        rgbfactors     : (r, g, b)
                         values
                         range from
                         0.0 to 1.0
                             or
                         (int, int)
                         packed RGB
        fn: optional gradient function


    Returns:
        new bitmap file
    """
    _use24btclrfnwithmultpar2circreg(
        ExistingBMPfile, NewBMPfile,
        sphere, x, y, r, rgbfactors, fn)


@functimer
def filledcircle2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        color: int):
    """Draws a Filled Circle

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
        color          : color of the
                         circular region

    Returns:
        new bitmap file
    """
    _usefnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        filledcircle, x, y, r, color)


@functimer
def circle2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        color: int):
    """Draws a Circle

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r
        color          : color of circle

    Returns:
        new bitmap file
    """
    _usefnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        circle, x, y, r, color)


@functimer
def thickencirclearea2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        rgbfactors: list[float, float, float]):
    """"Encircle area with a gradient and save to a file

    Args:
        bmp       : unsigned byte array
                    with bmp format

        x, y      : center of circle
        r         : radius of circle
        rgbfactors: (r, g, b) values
                    are from 0.0 to 1.0

    Returns:
        new bitmap file
"""
    _use24btclrfnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        thickencirclearea, x, y, r,
        rgbfactors)


@functimer
def brightnessadjcircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        percentadj: float):
    """Brightness gradient to a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
                         of a circular
                         region
                         in which we
                         apply a
        percentadj     : brightness
                         adjustment
                         that can be
                         positive
                         or negative

    Returns:
        new bitmap file
    """
    _use24btclrfnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        brightnessadjcircregion,
        x, y, r, percentadj)


@functimer
def vertbrightnessgrad2circregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        lumrange: list[int, int]):
    """Vertical brightness gradient to a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
                         of a circular
                         region
        lumrange       : brightness
                         gradient
                         (byte: byte)
                         adjust

    Returns:
        new bitmap file
    """
    _use24btclrfnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        vertbrightnessgrad2circregion,
        x, y, r, lumrange)


@functimer
def horibrightnessgrad2circregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        lumrange: list[int, int]):
    """Horizontal brightness gradient to a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
                         of a circular
                         region in which
                         we apply a
        lumrange       : brightness gradient
                         (byte, byte) adjust

    Returns:
        new bitmap file
    """
    _use24btclrfnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        horibrightnessgrad2circregion,
        x, y, r, lumrange)


@functimer
def flipvertcircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Vertical Flip of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        flipvertcircregion, x, y, r)


@functimer
def eraseeverynthhorilineinccircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        n: int):
    """Erase every nth horzontal line in a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
                         in which we
        n              : omit every
                         nth line

    Returns:
        new bitmap file
    """
    _usefnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        eraseeverynthhorizontallineinccircregion,
        x, y, r, n)


@functimer
def mirrortopincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Mirror the top-half of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        mirrortopincircregion, x, y, r)


@functimer
def mirrorbottomincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Mirror the bottom-half of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        mirrorbottomincircregion,
        x, y, r)


@functimer
def mirrorleftincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Mirror the left-half of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        mirrorleftincircregion,
        x, y, r)


@functimer
def mirrorrightincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Mirror the right-half of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        mirrorrightincircregion,
        x, y, r)


@functimer
def mirrortopleftincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Mirror the top-left of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        mirrortopleftincircregion,
        x, y, r)


@functimer
def mirrorbottomleftincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Mirror the bottom-left of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        mirrorbottomleftincircregion,
        x, y, r)


@functimer
def mirrortoprightincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Mirror the top-right of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        mirrortoprightincircregion,
        x, y, r)


@functimer
def mirrorbottomrightincircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Mirror the bottom-right of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        mirrorbottomrightincircregion,
        x, y, r)


@functimer
def fliphoricircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Horizontal Flip of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x,y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        fliphoricircregion, x, y, r)


@functimer
def outlinecircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Outlines area in a circular region

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        outlinecircregion, x, y, r)


@functimer
def flipXYcircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int):
    """Flips the x and y coordinates of a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r

    Returns:
        new bitmap file
    """
    _usefn2circreg(
        ExistingBMPfile, NewBMPfile,
        flipXYcircregion, x, y, r)


@functimer
def magnifyNtimescircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        intmagfactor: int):
    """Magnify a circular region by an integer factor n times

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
        intmagfactor   : int magnification
                         factor

    Returns:
        new bitmap file
    """
    _usefnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        magnifyNtimescircregion,
        x, y, r, intmagfactor)


@functimer
def pixelizenxncircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        intpixsize: int):
    """Apply a Pixel Blur in a circular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        x, y, r        : center (x, y)
                         and radius r
        intpixsize     : n by n
                         pixel blur size

    Returns:
        new bitmap file
    """
    _usefnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        pixelizenxncircregion,
        x, y, r, intpixsize)


@functimer
def copycircregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x: int, y: int, r: int,
        newxycenterpoint: list):
    """Copy/Paste of a circular area in a BMP

    Args:
        ExistingBMPfile : Whole path
                          to existing file
        NewBMPfile      : New file to
                          save changes to
        x, y, r         : center (x,y)
                          and radius r
        newxycenterpoint: (x:int,y:int)
                          where to paste

    Returns:
        new bitmap file
    """
    _usefnwithpar2circreg(
        ExistingBMPfile, NewBMPfile,
        copycircregion, x, y, r,
        newxycenterpoint)


@functimer
def horizontalbrightnessgrad2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        lumrange: list[int, int]):
    """Horizontal brightness gradient to a BMP

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes to
        lumrange       : (byte:byte)
                         defines the
                         brightness
                         gradient

    Returns:
        new bitmap file
    """
    _use24btclrfn(
        ExistingBMPfile, NewBMPfile,
        horizontalbrightnessgradto24bitimage,
        lumrange)


@functimer
def horizontalbrightnessgradregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        lumrange: list[int, int]):
    """Horizontal brightness gradient to a rectangular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2 ,y2 : defines the
                         rectangular
                         region
        lumrange       : (byte:byte)
                         brightness
                         gradient

    Returns:
        new bitmap file
    """
    _use24btbyrefclrfn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        horizontalbrightnessgradto24bitregion,
        lumrange)


@functimer
def verticalbrightnessgrad2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        lumrange: list[int, int]):
    """Applies a Vertical brightness gradient

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        lumrange       : (byte:byte)
                         defines the
                         brightness
                         gradient

    Returns:
        new bitmap file
    """
    _use24btclrfn(
        ExistingBMPfile, NewBMPfile,
        verticalbrightnessgradto24bitimage,
        lumrange)


@functimer
def verticalbrightnessgradregion2file(
        ExistingBMPfile: str,
        NewBMPfile: str,
        x1: int, y1: int,
        x2: int, y2: int,
        lumrange: list[int, int]):
    """Vertical brightness gradient to a rectangular area

    Args:
        ExistingBMPfile: Whole path to
                         existing file
        NewBMPfile     : New file to
                         save changes in
        x1, y1, x2, y2 : defines the
                         rectangular
                         region
        lumrange       : (byte:byte)
                         defines the
                         brightness
                         gradient

    Returns:
        new bitmap file
    """
    _use24btbyrefclrfn2regnsv(
        ExistingBMPfile, NewBMPfile,
        x1, y1, x2, y2,
        verticalbrightnessgradto24bitregion,
        lumrange)
