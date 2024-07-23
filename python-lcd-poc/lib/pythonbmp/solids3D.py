"""   3D math and solids module
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

from math import(
    sqrt,
    radians
    )

from numbers import Number
from typing import Callable

from .mathlib import(
    adddimz,
    addvect,
    computerotvec,
    cylindrical2rectcoord3D,
    dist,
    getnormvec,
    roundvect,
    spherical2rectcoord3D,
    subvect,
    trans
    )

from .primitives2D import(
    floatregpolygonvert,
    iterline,
    rectboundarycoords,
    regpolygonvert
    )

from .messages import sysmsg


def getshapesidedict() -> dict:
    """Returns a dictionary of side
        or polygon definitions for
        simple solids

    Returns:
        dictionary of side or polygon
        definitions for basic solids
    """
    return {"tetrahedra":((2, 1, 0),
                          (2, 3, 1),
                          (0, 3, 2),
                          (3, 0, 1)),
                  "cube":((1, 2, 3, 0),
                          (5, 6, 7, 4),
                          (0, 3, 6, 5),
                          (4, 7, 2, 1),
                          (4, 1, 0, 5),
                          (2, 7, 6, 3)),
             "hexahedra":((2, 3, 1),
                          (0, 3, 2),
                          (3, 0, 1),
                          (1, 4, 2),
                          (2, 4, 0),
                          (1, 0, 4)),
             "octahedra":((1, 2, 0),
                          (4, 1, 0),
                          (3, 4, 0),
                          (2, 3, 0),
                          (2, 1, 5),
                          (1, 4, 5),
                          (4, 3, 5),
                          (3, 2, 5))}


def tetrahedravert(x: float
                ) -> list[list[float,
                               float,
                               float]]:
    """Returns a list of vertices
        for a tetrahedron

    Args:
        x: length of a side

    Returns:
        list (x: float,
              y: float,
              z: float)
    """
    x_sqr = x**2
    halfx= x / 2
    return [[0, 0, 0],
            [halfx, sqrt(x_sqr / 2), 0],
            [x, 0, 0],
            [halfx, halfx, sqrt(3 / 8 * x_sqr)]]


def cubevert(x: float
            ) -> list[list[float,
                           float,
                           float]]:
    """Returns a list of vertices
        for a cube

    Args:
        x: length of a side

    Returns:
        list (x: float,
              y: float,
              z: float)
    """
    return [[0, 0, 0],
            [0, x, 0],
            [x, x, 0],
            [x, 0, 0],
            [0, x, x],
            [0, 0, x],
            [x, 0, x],
            [x, x, x]]


def hexahedravert(x: float
            ) -> list[list[float,
                           float,
                           float]]:
    """Returns a list of vertices
        for a hexahedron

    Args:
        x: length of a side

    Returns:
        list (x: float,
              y: float,
              z: float)
    """
    x_sqr = x**2
    halfx= x / 2
    z = sqrt(3 / 8 * x_sqr)
    return [[0, 0, 0],
            [halfx, sqrt(x_sqr / 2), 0],
            [x, 0, 0],
            [halfx, halfx, z],
            [halfx, halfx, -z]]


def octahedravert(x: float
        ) -> list[list[float,
                       float,
                       float]]:
    """Returns a list of vertices
        for an octrahedron

    Args:
        x: length of a side

    Returns:
        list (x: float,
              y: float,
              z: float)
    """
    halfx = x / 2
    return [[halfx, halfx, halfx],
            [0, 0, 0],
            [x, 0, 0],
            [x, x, 0],
            [0, x, 0],
            [halfx, halfx, -halfx]]


def decahedvertandsurface(
        x: float
        ) -> list[list[float,
                       float,
                       float]]:
    """Returns a list of vertices
        for a decahedron

    Args:
        x: min radius of sphere
           that can hold
           the decahedron

    Returns:
        list (x: float,
              y: float,
              z: float)
    """
    pts = regpolygonvert(0, 0, x, 5, 0)
    z = sqrt((dist(pts[0], pts[1]) ** 2 - x**2))

    return [[[0, 0, -z]] + \
            adddimz(pts, 0) + \
            [[0, 0, z]],
            ((1, 2, 0),
             (5, 1, 0),
             (3, 4, 0),
             (2, 3, 0),
             (4, 5, 0),
             (2, 1, 6),
             (1, 5, 6),
             (4, 3, 6),
             (3, 2, 6),
             (5, 4, 6))]


#don't edit this it took much
# computation to make
def icosahedvertandsurface(
        x: float) -> list[list[list[float,
                                    float,
                                    float]],
                          tuple]:
    """Returns a list of vertices
        and surfaces for an icosahedron

    Args:
        x: min radius of sphere
           that can hold
           the icosahedron

    Returns:
        list (x: float,
              y: float,
              z: float)
    """
    pts = floatregpolygonvert(
                0, 0, x, 5, 0)
    pts1 = floatregpolygonvert(
                0, 0, x, 5, 36)
    z = sqrt((dist(pts[0], pts[1]) ** 2 - x**2))
    z1 = 2 * x - z
    return [[[0, 0, -z]] + \
             adddimz(pts, 0) + \
             adddimz(pts1, z1 - z) + \
            [[0, 0, z1]],
           ((1, 2, 0),
            (5, 1, 0),
            (3, 4, 0),
            (2, 3, 0),
            (4, 5, 0),
            (2, 1, 6),
            (1, 5, 10),
            (4, 3, 8),
            (3, 2, 7),
            (5, 4, 9),
            (6, 7, 2),
            (7, 8, 3),
            (8, 9, 4),
            (9, 10, 5),
            (10, 6, 1),
            (7, 6, 11),
            (9, 8, 11),
            (8, 7, 11),
            (10, 9, 11),
            (6, 10, 11))]


def rotvec3D(roll: float,
            pitch: float,
              yaw: float) -> tuple:
    """Returns a 3D rotation vector

    Args:
        All input arguements are in
        degrees (roll, pitch, yaw)

    Returns:
        tuple ((float, float),
               (float, float),
               (float, float))
    """
    return (computerotvec(roll),
            computerotvec(pitch),
            computerotvec(yaw))


#translated from C code by Roger Stevens
def perspective(
        vlist: list[list[Number,
                         Number,
                         Number]],
        rotvec: list[list[float, float],
                     list[float, float],
                     list[float, float]],
        dispvec: list[Number,
                      Number,
                      Number],
        d: float) -> tuple:
    """Projects 3D points to 2D and
        apply rotation and translation
        vectors

    Args:
        vlist  : list of 3D vertices
        rotvec : 3D rotation vector
        dispvec: 3D translation vector
        d      : Distance of observer
                 from the screen

    Returns:
        tuple (list, list)
    """
    projvlist = []
    rotvlist = []
    ((sroll, croll),
    (spitch, cpitch),
    (syaw, cyaw)) = rotvec
    (dx, dy, dz) = dispvec
    for p in vlist:
        (px, py, pz) = p
        x1 = -cyaw * px - syaw * pz
        y1 = croll * py - sroll * x1
        z1 = -syaw * px + cyaw * pz
        x = croll * x1 + sroll * py
        y = spitch * z1 + cpitch * y1
        z = cpitch * z1 - spitch * y1
        x += dx
        y += dy
        z += dz
        rotvlist.append([x, y, z])
        projvlist.append([-d * x / z,
                          -d * y / z])
    return (rotvlist, projvlist)


#may be slow if polygon goes offscreen
def fillpolydata(
        polybnd: list[list[int, int]],
        xlim: int,
        ylim: int) -> list:
    """Generates a list of x values per
        y values for filling polygon
        boundaries

    Args:
        polybnd : list of 2D vertices
                  list[list[x: int,
                            y: int]]
                  that forms a
                  complete polygon
                  boundary (see
                  def polyboundary)
        xlim    : Screen limit x dim
        ylim    : Screen limit x dim

    Returns:
        A dictionary with y values
        as key and list of x values
        per key (if within bounds)
        or an empty list if out of
        bounds
    """
    filld = {}
    ((minx, miny), (maxx, maxy)) = \
    rectboundarycoords(polybnd)
    maxx += 1
    maxy += 1
    if (minx >= 0 and miny >= 0) and \
       (maxx <= xlim and maxy <= ylim):
        for y in range(miny, maxy):
            filld.setdefault(y, [])
            for x in range(minx, maxx):
                if [x, y] in polybnd:
                    filld[y] += [x]
    else:
        print(sysmsg['regionoutofbounds'])
        filld = []
    return filld


def polyboundary(
        vertlist: list[list[int, int]]
             ) -> list[list[int, int]]:
    """Generates a polygon boundary
        from a list of 2D vertices

    Args:
        polybnd : list of 2D vertices
                  list[list[x: int,
                            y: int]]

    Returns:
        list[list[x: int, y: int]]
        A list of vertices that traces
        the boundaries of the polygon
    """
    px = []
    vertcount = len(vertlist)
    for i in range(vertcount):
        if i>0:
            v1 = roundvect(
                    vertlist[i - 1])
            v2 = roundvect(vertlist[i])
            for p in iterline(v1, v2):
                if p not in px:
                    px.append(p)
    v2 = roundvect(vertlist[0])
    v1 = roundvect(
            vertlist[vertcount - 1])
    for p in iterline(v1, v2):
        if p not in px:
            px.append(p)
    return px


def gensides(
        pointlists: list[list, list],
        transvect: list[float,
                        float,
                        float],
        sides: list[list[float]]
                          ) -> tuple:
    """Generates a list of polygons
        and normals from 3D polygon
        and side data to a list of
        sides and normals with
        with the hidden surfaces
        removed that is suitable
        for rendering on a 2D surface
        and also applies a
        3D translation vector for
        positioning
    """
    (plist, slist) = pointlists
    polylist = []
    normlist = []
    for sidepts in sides:
        (s0, s1, s2) = sidepts[0:3]
        u = plist[s0]
        v = plist[s1]
        w = plist[s2]
        if surfacetest(u, v, w) <= 0:
            polylist.append(
                trans(
                    [slist[i]
                     for i in sidepts], transvect))
            normlist.append(getnormvec(u, v, w))
    return (polylist, normlist)


def spherevert(vcen: list[float,
                          float,
                          float],
                   r: float,
        deganglestep: float) -> list:
    """Returns a list of sparse vertices
        for a sphere

    Args:
        vcen       : [x: float, center
                      y: float, of the
                      z: float] sphere
        r           : spherical radius
        deganglestep: angle step between
        vertices that controls how sparse
        the list will be

    Returns:
        list (x: float,
              y: float,
              z: float)
    """
    plist = []
    for theta in range(0, 360,
                deganglestep):
        for phi in range(0, 180,
                  deganglestep):
            p = addvect(vcen,
                   spherical2rectcoord3D(
                       [r, radians(theta),
                           radians(phi)]))
            if p not in plist:
                plist.append(p)
    p = [0, 0, -r]
    if p not in plist:
        plist.append(p)
    return plist


#we have a 3D object and we take z slices
def zlevelcoords(verlist: list
                )     -> list[list, list]:
    """Takes a list of 3D vertices
        and makes slices along the
        z-axis

    Args:
        verlist: list of 3D vertices

    Returns:
        (list of z-axis slices,
         dict of z-axis slices
         with 3D vertices)
    """
    zlist = []
    zord = {}
    for p in verlist:
        pind = [verlist.index(p)]
        z = p[2]
        if z not in zord:
            zord.setdefault(z, pind)
            zlist.append(z)
        else:
            zord[z] += pind
    return (zlist, zord)


def genspheresurfaces(
        zlevelcoord:list[list, list]
                     ) -> list:
    """Generates the tiled surfaces
        of a sphere

    Args:
        zlevelcoord:
        (list of z-axis slices,
         dict of z-axis slices
         with 3D vertices)

    Returns:
        list of polygons to render
        the surface of a sphere
    """
    (zl, vl) = zlevelcoord
    surf = []
    levels = len(zl) - 1
    northpole = vl[zl[0]][0]
    southpole = vl[zl[levels]][0]
    npts = vl[zl[1]]
    spts = vl[zl[levels - 1]]
    lpts = len(npts)
    maxp = lpts - 1
    i = 0
    for i in range(lpts):
        j = 0 if i == maxp else i + 1
        surf += [[northpole, npts[j], npts[i]],
                 [southpole, spts[i], spts[j]]]
    if levels > 2:
        mxlv = levels - 1
        for k in range(1, mxlv):
            pts = vl[zl[k]]
            adjpts = vl[zl[k + 1]]
            lpts = len(pts)
            maxadjpts = len(adjpts)
            maxp = lpts - 1
            i = 0
            for i in range(lpts):
                j = 0 if i == maxp else i + 1
                surf += [[adjpts[j % maxadjpts],
                          adjpts[i % maxadjpts],
                          pts[i], pts[j]]]
    return surf


def spherevertandsurface(
        vcen: list[float, float, float],
        r: float,
        deganglestep: float)-> tuple:
    """Returns a list of sparse
        vertices and tiled surfaces
        for a sphere

    Args:
        vcen       : [x: float, center
                      y: float, of the
                      z: float] sphere
        r           : spherical radius
        deganglestep: angle step between
        vertices that controls how sparse
        the list will be

    Returns:
        list of vertices and surfaces
        for plot3Dsolid()
        see Hello_DiscoBall.py
        and Hello_Globe.py
    """
    s = spherevert(vcen, r, deganglestep)
    return (s, genspheresurfaces(zlevelcoords(s)))


def cylindervertandsurface(
        vcen: list[float, float, float],
        r: float,
        zlen: float,
        deganglestep: float) -> tuple:
    """Returns a list of sparse vertices
       and tiled surfaces for a cylinder

    Args:
        vcen       : [x: float, center
                      y: float, of the
                      z: float] sphere
        r           : radius
        zlen        : height of the
                      cylinder
        deganglestep: angle step between
        vertices that controls how sparse
        the list will be

    Returns:
        list of vertices and surfaces
        for plot3Dsolid()
        see Hello_Coin.py
    """
    z = zlen / 2
    i = 0
    plist = []
    top = []
    bottom = []
    side = []
    maxang = 360 + (deganglestep << 1)
    for theta in range(0, maxang, deganglestep):
        plist.extend(
            (
                addvect(
                    vcen, cylindrical2rectcoord3D([r, radians(theta), -z])
                ),
                addvect(vcen, cylindrical2rectcoord3D([r, radians(theta), z])),
            )
        )

        top.append(i)
        bottom.append(i + 1)
        if i > 2:
            j = i - 2
            side.append([i, i + 1, j + 1, j])
        i += 2
    top.reverse()
    return (plist, [top] + side + [bottom])


def conevertandsurface(
        vcen: list[float, float, float],
        r: float,
        zlen: float,
        deganglestep: float) -> tuple:
    """Returns a list of sparse
        vertices and tiled surfaces
        for a cone

    Args:
        vcen       : [x: float, center
                      y: float, of the
                      z: float] sphere
        r           : radius of
                      conical base
        zlen        : height of
                      the cone
        deganglestep: angle step between
        vertices that controls how sparse
        the list will be

    Returns:
        list of vertices and surfaces
        for plot3Dsolid()
        see Hello_Cone.py
    """
    z = zlen / 2
    bottom = []
    side = []
    maxang = 360 + (deganglestep << 1)
    plist = [subvect(vcen, [0, 0, z])]
    for i, theta in enumerate(range(0, maxang, deganglestep), start=1):
        plist.append(addvect(vcen,
                cylindrical2rectcoord3D(
                    [r, radians(theta), z])))
        bottom.append(i)
        if i > 2:
            side.append([0, i, i - 1])
    return (plist, side + [bottom])


def surfplot3Dvertandsurface(
        x1: int, y1: int,
        x2: int, y2: int,
        step: int, fnxy: Callable
        ) -> tuple:
    """Does a 3D surface plot of a
        function z = fnxy(x, y)

    Args:
        x1, y1, x2, y2: set drawing area
        fnxy          : fnxy(x, y)
                        Callable
                        (lambda or fn)

    Returns:
        list of vertices and surfaces
        for plot3Dsolid()
        see Hello_3D_surfaceplot.py
    """
    vlist = []
    surf = []
    for y in range(y1, y2, step):
        for x in range(x1, x2, step):
            z = fnxy(x, y)
            vlist.append([x, y, z])
    dx = abs(x2 - x1) // step
    dx1 = dx - 1
    vl = len(vlist)
    for v in vlist:
        i = vlist.index(v)
        idx = i + dx
        idx1 = idx + 1
        if (vl - idx1) >= 0 and (i % dx) < dx1:
             surf.append([idx, idx1, i + 1, i])
    return (vlist, surf)


def surfacetest(
        p1: list[float, float, float],
        p2: list[float, float, float],
        p3: list[float, float, float]
               ) -> float:
    """Hidden surface removal test
        returns a positive number
        if hidden and negative
        number if not

    Args:
        p1, p2, p3 : [x: float,
                      y: float,
                      z: float]
        the first 3 vertices of a
        surface or polygon that
        can be triangle or square

    Returns:
        float -> a positive number
        if hidden and negative
        number if not
    """
    (p1x, p1y, p1z) = p1
    (p2x, p2y, p2z) = p2
    (p3x, p3y, p3z) = p3
    return p1x * (p3y * p2z - p2y * p3z) - \
           p2x * (p3y * p1z - p1y * p3z) - \
           p3x * (p1y * p2z - p2y * p1z)
