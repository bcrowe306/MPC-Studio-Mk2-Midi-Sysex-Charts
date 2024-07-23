"""
Math library for vectors, statistics
basic linear algebra, bit operations
and other maths
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
from cmath import exp
from math import(
    asin,
    acos,
    atan,
    cos,
    cosh,
    degrees,
    pi,
    tau,
    radians,
    sin,
    sinh,
    sqrt
    )

from numbers import Complex, Number
from random import randint
from functools import reduce
from typing import Callable
from .conditionaltools import iif


def setmaxvec(
        vlist: list[Number],
        maxval: Number) -> list[Number]:
    """Sets the maximum value in a list of numbers to maxval

    Args:
        vlist : list of Numbers
        maxval: upper limit for
                the list of
                Numbers

    Returns:
        list of Numbers
    """
    return [setmax(v, maxval)
               for v in vlist]


def setminmaxvec(
        vlist: list[Number],
        minval: Number,
        maxval: Number) -> list[Number]:
    """Sets the minimum and maximum values in a list of numbers
    (a vector) to minval and maxval respectively

    Args:
        vlist : list of Numbers
        minval: lower limit for
                the list of
                Numbers
        maxval: upper limit for
                the list of
                Numbers

    Returns:
        list of Numbers
    """
    return [setminmax(v, minval, maxval)
                  for v in vlist]


def intsetminmaxvec(
        vlist: list[Number],
        minval: int,
        maxval: int) -> list[int]:
    """Sets the minimum and maximum values in a list of numbers
    to minval and maxval respectively and return int values

    Args:
        vlist : list of Numbers
        minval: lower limit for
                the list of
                Numbers
        maxval: upper limit for
                the list of
                Numbers

    Returns:
        list of ints
    """
    return [intsetminmax(v, minval, maxval)
                for v in vlist]


def range2baseanddelta(
        lst_range: list[int, int]):
    """Gets the base and range values in a list of numbers

    Args:
        lst_range: list[min: int,
                        max: int]

    Returns:
        minimum value and delta of min and max value
        in lst_range
    """
    return lst_range[0], delta(lst_range)


def roundvect(v: list[Number]
            ) -> list[int]:
    """Rounds off the components of a vector
       (list of floats -> list of ints)

    Args:
        v: list of floats

    Returns:
        list of ints
    """
    return [round(n) for n in v]


def roundvectlist(
        vlist: list[list[Number]]
          ) -> list[list[int]]:
    """Rounds the components of a vector in a list within a list
        list of (list of floats) -> list of (list of ints)

    Args:
        vlist: list[list[floats]

    Returns:
        list[list[ints]
    """
    return [roundvect(v) for v in vlist]


def addvect(u: list[Number],
            v: list[Number]
          ) -> list[Number]:
    """Add vectors u and v by adding their components

    Args:
        u, v: list of ints or floats

    Returns:
        list of ints or floats
    """
    return [i + j for i, j in zip(u, v)]


def trans(vlist: list[list[Number]],
              u: list[Number]
            ) -> list[list[Number]]:
    """Translates list of vectors by adding vector u
    to all vectors in the list of vectors

    Args:
        vlist: list of vectors
        u    : translation vector

    Returns:
        list of vectors
    """
    return [addvect(v, u) for v in vlist]


def subvect(u: list[Number],
            v: list[Number]
          ) -> list[Number]:
    """Subtracts vectors u and v by subtracting their components

    Args:
        u, v: list of ints or floats

    Returns:
        list of ints or floats
    """
    return [i - j
            for i, j in zip(u, v)]


def mulvect(u: list[Number],
            v: list[Number]
          ) -> list[Number]:
    """Gets the inner product of vector u and vector v

    Args:
        u, v: list of ints or floats

    Returns:
        list of ints or floats
    """
    return [i * j for i, j in zip(u, v)]


def divvect(u: list[Number],
            v: list[Number]
          ) -> list[float]:
    return [i / j for i, j in zip(u, v)]


def scalarmulvect(
            v: list[Number],
         scalarval: Number
          ) -> list[Number]:
    """Scales a vector by multiplying a scalar value (float)
    to all components of the vector or a list of numbers

    Args:
        v        : the vector or
                   a list of
                   ints or floats
        scalarval: scalar value
                   (float or int)

    Returns:
        list of ints or floats
    """
    return [s * scalarval for s in v]


def scalardivvect(
            v: list[Number],
         scalarval: Number
          ) -> list[Number]:
    """Scales a vector by dividing with a scalar value (float)
    all the components of the vector or a list of numbers

    Args:
        v        : the vector or
                   a list of
                   ints or floats
        scalarval: scalar value
                   (float or int)

    Returns:
        list of ints or floats
    """
    return [s / scalarval for s in v]


def intscalarmulvect(vec: list[Number],
                    scalarval: Number
                    ) -> list[int]:
    """Scales a vector by multiplying a scalar value (float)
    to all components of the vector or a list of numbers
    then rounds off values in the list to a whole number

    Args:
        v        : the vector or
                   a list of
                   ints or floats
        scalarval: scalar value
                   (float or int)

    Returns:
        list of ints
    """
    return [round(s * scalarval)
            for s in vec]


def mean(v: list[Number]) -> float:
    """Gets the average of a list of numbers

    Args:
        v: the vector or a list of
           ints or floats

    Returns:
        float
    """
    return sum(v) / len(v)


def meanlist(
      vlist: list[list[float]]
        ) -> list[float]:
    """Gets the averages of a list of list of numbers

    Args:
        vlist: the vector or a list
               of ints or floats

    Returns:
        list[float, ...]
    """
    return [mean(v) for v in vlist]


def pivotlist(
        vlist:list[list[any]]
         ) -> list[list[any]]:
    j = len(vlist[0])
    """Pivot table with a list of list of anything

    Args:
        vlist: a list of lists

    Returns:
        list of lists
    """
    return [[v[i] for v in vlist] for i in range(j)]


def variance(v: list[Number]
           ) -> list[float]:
    """Computes the variance of the elements in a list of numbers

    Args:
        vlist: a list of
               ints or floats


    Returns:
        list of floats
    """
    ave = mean(v)
    return [n - ave for n in v]


def isinrange(value: Number,
          highlimit: Number,
           lowlimit: Number) -> bool:
    """Checks is value is within high and low limits

    Args:
        value    : numeric variable
                   to check
        highlimit: upper limit of
                   the variable
        lowlimit : lower limit of
                   the variable

    Returns:
        True if within bounds
    """
    return (value > lowlimit) and \
           (value < highlimit)


def setmax(val: Number,
        maxval: Number) -> Number:
    """Set the value of val to maxval if val > maxval

    Args:
        val   : numeric variable
        maxval: upper limit of variable

    Returns:
        Number
    """
    return maxval if val > maxval else val


def setmin(val: Number,
        minval: Number) -> Number:
    """Set the value of val to minval if val < minval

    Args:
        val   : numeric variable
        minval: lower limit of variable

    Returns:
        Number
    """
    return minval if val < minval else val


def setminmax(val: Number,
           minval: Number,
           maxval: Number) -> Number:
    """Set the value of val to minval if val < minval
    or the value of val to maxval if val > maxval

    Args:
        val   : numeric variable
        minval: lower limit of variable
        maxval: upper limit of variable

    Returns:
        Number
    """
    if val > maxval:
        val = maxval
    if val < minval:
        val = minval
    return val


def intsetminmax(val: Number,
              minval: int,
              maxval: int) -> int:
    """Set the value of val to minval if val < minval
    or the value of val to maxval if val > maxval
    and return an int

    Args:
        val   : numeric variable
        minval: lower limit of variable
        maxval: upper limit of variable

    Returns:
        int
    """
    val = round(val)
    val = min(val, maxval)
    val = max(val, minval)
    return val


def sign(val: Number) -> int:
    """Returns an int depending on the sign of val
    Args:
        val   : a Number

    Returns:
        Postive value ->  1
        zero          ->  0
        Negative value-> -1
    """
    retval = 0
    if val > 0:
        return 1
    elif val == 0:
        return 0
    else:
        return -1


def LSMslope(XYdata: list) -> float:
    """Slope of a line obtained by Least Squares Method

    Args:
        XYdata: list of vectors
        first two values in the
        list must be [[x, y, ..,], ...]

    Returns:
        float slope of line
    """
    XY , N = pivotlist(XYdata), len(XYdata)
    X =  XY[0]
    Y =  XY[1]
    EX = sum(X)
    EY = sum(Y)
    EXY = sum(mulvect(X,Y))
    EsqX = sum(mulvect(X,X))
    return ((N * EXY) - (EX * EY)) / \
           ((N * EsqX) - (EX ** 2))


def LSMYint(XYdata: list) -> float:
    """Returns the y-intercept of a line obtained
    by Least Squares Method

    Args:
        XYdata: list of vectors
        first two values in the
        list must be [[x, y, ..,], ...]

    Returns:
        float y-intercept of line
    """
    XY = pivotlist(XYdata)
    meanX = mean(XY[0])
    meanY = mean(XY[1])
    return meanY - LSMslope(XYdata) * meanX


def PearsonsR(XYdata: list) -> float:
    """Pearsons r is a measure of linear
correlation between two sets of data.

It is the ratio between the covariance of
two variables and the product of their
standard deviations; thus it is essentially
a normalized measurement of the covariance,
such that the result always has a value
between -1 and 1.

A value of +1 implies that all data points lie
on a line for which Y increases as X increases,
and vice versa for -1.

A value of 0 implies that there is no linear
dependency between the variables.

    Args:
        XYdata: list[list[x: Number,
                          y: Number,
                               ...]]

    Returns:
        a value between -1 and 1
"""
    XY = pivotlist(XYdata)
    N = len(XYdata)
    X = XY[0]
    Y = XY[1]
    EX = sum(X)
    EY = sum(Y)
    EXY = sum(mulvect(X, Y))
    EsqX = sum(mulvect(X, X))
    EsqY = sum(mulvect(Y, Y))
    return ((N * EXY) - (EX * EY)) / \
        sqrt(abs(((N * EsqX) - (EX ** 2))) * \
            abs(((N * EsqY) - (EY ** 2))))


def Rsquare(XYdata: list) -> float:
    """Squared Pearson correlation coefficient
not to be confused with the coefficient of determination (R²)

In linear least squares multiple regression with an estimated
intercept term, R2 equals the square of the Pearson correlation
coefficient between the observed y and modeled (predicted)
data values of the dependent variable.

In a linear least squares regression with an intercept term
and a single explanator, this is also equal to the squared
Pearson correlation coefficient ofthe dependent variable y
and explanatory variable x.

     Args:
        XYdata: list[list[x: Number,
                          y: Number,
                               ...]]

    Returns:
        a value between 0 and 1
    """
    return PearsonsR(XYdata) ** 2


def TTest(XYdata: list) -> float:
    r = PearsonsR(XYdata)
    return r * sqrt((len(XYdata) - 2) / (1 - r ** 2))


def StdDev(v: list[Number]) -> float:
    """Computes the standard deviation of a list of numbers

    Args:
        v:  list of numbers

    Returns:
        float
    """
    vr = variance(v)
    return sqrt(sum(mulvect(vr, vr)) / len(v))


def slope(u: list[float, float],
          v: list[float, float]) -> float:
    """Computes the slope of a line
    using points u(x1, y1) and v(x2, y2)

    Args:
        u, v:  line endpoints (x, y)

    Returns:
        float
    """
    return (v[1] - u[1]) / \
           (v[0] - u[0])


def coefvar(v: list[Number]) -> float:
    """Coefficient of variation (CV) is a statistical measure
    of the relative dispersion of data points in a data series
    around the mean.

    It is also known as (RSD) or relative standard deviation
    and is defined as the ratio of the standard deviation
    to the mean

    Args:
        v: list of ints or floats

    Returns:
        float
    """
    return StdDev(v) / mean(v)


def vectiszero(v: list[Number]) -> bool:
    """Checks if a vector is all zero

    Args:
        v: list of ints or floats

    Returns:
        True  :if list is all zero
        False :if there is any non zero
               in the list
    """
    return all(u == 0 for u in v)


def isorthogonal(
        u: list[Number],
        v: list[Number]) -> bool:
    """Checks if vector u and vector v are orthogonal to each other

    Args:
        u, v: list of ints or floats

    Returns:
        True  : vectors are orthogonal
                since the inner product
                of u and v is a
                zero vector or a
                null vector that
                has a zero magnitude
                and no direction

        False : inner product of
                u and v is a
                non zero vector
                so they are
                nonorthogonal vectors
    """
    return vectiszero(mulvect(u, v))


def crossprod3d(
        u: list[float, float, float],
        v: list[float, float, float]
      ) -> list[float, float, float]:
    """Compute the cross product of 3D vectors u and v.

    The cross product is perpendicular to both u and v

    Args:
        u, v: list of ints
              or floats
              [x, y, z]

    Returns:
        [x, y, z] : after u x v
                    where x is
                    crossprod
    """
    x = y = z = 0
    if len(u) == 3 and len(v) == 3:
        x = u[1] * v[2] - u[2] * v[1]
        y = u[2] * v[0] - u[0] * v[2]
        z = u[0] * v[1] - u[1] * v[0]
    return [x, y, z]


def getnormvec(
        p1: list[float, float, float],
        p2: list[float, float, float],
        p3: list[float, float, float]) -> list:
    """Gets the normal to a plane given 3 vertices (x, y, z)
    that define the surface in 3D spaceusing a cross product
    between two vectors a and b:

    The normal vector is orthogonal to vectors a and b.
    (see below)

        a = (p2 - p1)
        b = (p3 - p1)
        n = (a x b) where x is crossprod
        n = ((p2 - p1) x (p3 - p1))
    Returns:
        A Normal Vector to the surface
    """
    return crossprod3d(subvect(p2, p1),
                       subvect(p3, p1))


def dotprod(u: list[Number],
            v: list[Number]) -> float:
    """Dot product is an algebraic operation that takes two
    equal-length sequences of numbers and returns a float.

    It is the sum of the products of the corresponding entries
    of the two sequences of numbers

    Args:
        u, v : list of ints or floats

    Returns:
        float
    """
    return sum(mulvect(u, v))


def vmag(v: list[float]) -> float:
    """Compute the Magnitude or length of a vector v
    of arbitrary dimension n equal to len(v)

    Args:
        v: list of ints or floats

    Returns:
        float
    """
    return sqrt(sum(i * i for i in v))


def dist(u: list[float],
             v: list[float]) -> float:
    """Compute the Distance or length of a vector v
    of arbitrary dimension n in a n-dimensional
    Euclidean space where u and v are both vectors
    with n components

    Args:
        u, v: list of ints or floats

    Returns:
        float
    """
    return vmag(subvect(u, v))


def distancetable(vertlist: list) -> list:
    """Compute the Distances between the vectices in a list

    Args:
        vertlist: list of vertices

    Returns:
        list of vertices and distances
    """
    dlist = []
    for v in vertlist:
        dlist.extend(
            [vertlist.index(v), vertlist.index(u), dist(u, v)]
            for u in vertlist
            if u != v
        )

    return dlist


def countdist(distlist: list) -> dict:
    d = {}
    for i in distlist:
        dist = i[2]
        if dist not in d:
            d.setdefault(dist,1)
        else:
            d[dist] += 1
    return d


def det3D(
    a1: list[Number, Number, Number],
    a2: list[Number, Number, Number],
    a3: list[Number, Number, Number]
         ) -> float:
    """Compute the Determinant of a 3D matrix

    Args:
        a1, a2, a3 : list of ints
                     or floats
                     [x, y, z]

    Returns:
        float
    """
    return a1[0] * a2[1] * a3[2] + \
           a1[1] * a2[2] * a3[0] + \
           a1[2] * a2[0] * a3[1] - \
           a1[0] * a2[2] * a3[1] - \
           a1[1] * a2[0] * a3[2] - \
           a1[2] * a2[1] * a2[0]


def cosaffin(u: list[Number],
             v: list[Number]) -> float:
    """Compute Cosine Similarity or Cosine Affinity

    Args:
        u, v : list of ints or floats

    Returns:
        float value
            proportional vectors = 1
            orthogonal vectors = 0
            opposite vectors = -1
            and values in between
    """
    return dotprod(u, v) / \
         (vmag(u) * vmag(v))


def dircos(v: list[Number]
         ) -> list[float]:
    """Direction cosine of the vector computed by dividing
    the corresponding coordinates of a vector by the vector length.

    The unit vector coordinate is equal to thedirection cosine.
    One such property of the direction cosine is that the
    addition of the squares of the direction cosines is
    equivalent to one

    Args:
        v : list of ints or floats
            (vector components)

    Returns:
        Direction Cosine list
        list[float,...]
    """
    mag = vmag(v)
    return [i / mag for i in v]


def diracos(dcos: list[float]
              )-> list[float]:
    """Takes the Direction cosine of the vector and returns
    a list of angles in radians

    Args:
        v : Direction Cosine
            list[float]

    Returns:
        list of angles in radians per vector component
        list[float,...]

    """
    return [acos(d) for d in dcos]


def dirdeg(raddir: list[float]
              ) -> list[float]:
    """Takes the Directional angles of the vector in radians
    and returns a list of angles in degrees

    Args:
        v : list of angles
            per vector component
            in radians

    Returns:
        list of angles in degrees per vector component
        list[float,...]
    """
    return [degrees(d) for d in raddir]


def rect2sphericalcoord3D(
        v: list[Number, Number, Number]
      ) -> list[float, float, float]:
    """3D coordinate transform from rectangular to spherical

        p = The length of the hypotenuse
            or the magnitude of the
            vector

        theta = is the angle between the
                positive x-axis and p
                (azimuth)

        phi = is the angle between the
              positive z-axis and p
              (colatitude)

    Args:
        vspherecoord: [p, theta, phi]
                      spherical
                      coordinates

    Returns:
        [p: float,
        theta: float,
        phi: float]
    """
    p = vmag(v)
    return [p, atan(v[1] / v[0]),
               acos(v[2] / p)]


def spherical2rectcoord3D(
    vspherecoord: list[float, float, float]
             ) -> list[float, float, float]:
    """3D coordinate transform from spherical to rectangular

        p = The length of the hypotenuse
            or the magnitude of the
            vector

        theta = is the angle between the
                positive x-axis and p
                (azimuth)

        phi = is the angle between the
              positive z-axis and p
              (colatitude)

    Args:
        vspherecoord: [p, theta, phi]
                      spherical
                      coordinates

    Returns:
        [x:float, y: float, z: float]
    """
    (p, theta, phi) = vspherecoord
    sinphi = sin(phi)
    return [p * sinphi * cos(theta),
            p * sinphi * sin(theta),
            p * cos(phi)]


def rect2cylindricalcoord3D(
        v: list[Number, Number, Number]
      ) -> list[float, float, float]:
    """Converts rectangular to cylindrical coordinates
    with origin at (0, 0)

    Args:
        v: a vector (x: Number,
                     y: Number,
                     z: Number)

    Returns:
        [magnitude: float,
             angle: float,
             z: Number]
    """
    return [vmag(v),
            atan(v[1] / v[0]),
            v[2]]


def cylindrical2rectcoord3D(
        vcylindcoord: list[float,
                           float,
                           float]
        ) -> list[float, float, float]:
    """Converts from cylindrical coordinates with
    origin at (0, 0) to 3D rectangular coordinates

    Args:
        vcylindcoord:(r: float,
                  theta: float,
                      z: float)

    Returns:
        [x: float, y: float, z: float]
    """
    (r, theta, z) = vcylindcoord
    return [r * cos(theta),
            r * sin(theta), z]


def polar2rectcoord2D(
        vpolarcoord: list[float,
                            float]
                  ) -> list[float,
                            float]:
    """Converts from polar coordinates with
    origin at (0, 0) to 2D rectangular coordinates

    Args:
        vcylindcoord:(r: float,
                  theta: float)

    Returns:
        [x: float,
         y: float]
    """
    (r, theta) = vpolarcoord
    return [r * cos(theta),
            r * sin(theta)]


def rect2polarcoord2D(
        v: list[Number, Number]
      ) -> list[float, float]:
    """Converts rectangular to polar coordinates with
    origin at (0,0)

    Args:
        v: a vector (x: Number,
                     y: Number)

    Returns:
        [magnitude: float, angle: float]
    """
    return [vmag(v), atan(v[1] / v[0])]


def polarcoordangle2D(v:list) -> float:
    """Gets the polar coordinate angle from a vector v(x, y)
    with the origin set at (0, 0)

    Args:
        v: a vector (x: Number, y: Number)

    Returns:
        float angle in radians
    """
    a = acos(cosaffin(v, [0, -1]))
    if v[0] < 0:
        a = 2 * pi - a
    return a


def rect2polarcoord2Dwithcenter(
        vcen: list[Number, Number],
        vpnt: list[Number, Number]
         ) -> list[float, float]:
    """Converts rectangular to polar coordinates with
    origin at vcen(x, y)

    Args:
        vcen: origin (x: Number, y: Number)

        vpnt: a vector (x: Number, y: Number)

    Returns:
        [magnitude: float, angle: float]
    """
    v = subvect(vpnt, vcen)
    return [vmag(v),
            polarcoordangle2D(v)]


def computerotvec(degrot: float
                ) -> list[float,
                          float]:
    """2D rotation vector based on angle degrot in degrees

    Args:
        degrot: rotation in degrees

    Returns:
        [float, float] rotation vector
    """
    a = radians(degrot)
    return (sin(a), cos(a))


def rotvec2D(v: list[Number, Number],
        rotvec: list[float, float]
           ) -> list[float, float]:
    """2D rotation vector based on angle degrot in degrees and v(x, y)

    Args:
        v     : centerpoint
        degrot: rotation in degrees

    Returns:
        [float, float] rotation vector
    """
    (x, y) = v
    (a, b) = rotvec
    return [x * b - y * a,
            x * a + y * b]


def mirrorx(p: list[Number, Number], x: Number
          ) -> list[Number, Number]:
    """Mirrors the x value at point p(ox, y) to return
    two points at p1(ox - x, y) and p2(ox + x, y)

    Args:
        p(ox, y): original point
        x       : value to mirror

    Returns:
        Two 2D points [x1, y], [x2, y]
    """
    (a, y) = p
    return [a - x, y], [a + x, y]


def mirrory(p: list[Number, Number], y: Number
         ) -> list[Number, Number] :
    """Mirrors the y value at point p(x, oy) to return
    two points at p1(x, oy - y) and p2(x, oy + y)

    Args:
        p(x, oy): original point
        y       : value to mirror

    Returns:
        Two 2D points [x, y1], [x, y2]
    """
    (x, b) = p
    return [x, b - y], [x, b + y]


def mirrorvec(vcen: list[Number],
                 v: list[Number]
               ) -> list[Number]:
    """Mirrors a vector v of arbitrary dimension n
    in a n-dimensional Euclidean space where
    v is the vector to mirror and vcen is the origin point
    both with n components

    Args:
        v, vcen: list of numbers
                 with equal lengths

    Returns:
        list of two vectors
        [list[Number], list[Number]]
    """
    return [subvect(vcen, v),
            addvect(vcen, v)]


def mirror(pt: float, delta: float):
    """Mirrors a value in a numberline

    Args:
        pt   : real value in numberline
        delta: value to mirror

    Returns:
        pt - delta, pt + delta
    """
    return pt - delta, pt + delta


def randomvect(minrnd: int,
               maxrnd: int
              ) -> list[int, int, int]:
    """Creates a random 3D vector

    Args:
        minrnd: minimum int random value
                of 3D vector component
        maxrnd: maximum int random value
                of 3D vector component
    Returns:
        list[x: int, y: int , z: int]
    """
    return [randint(minrnd, maxrnd),
            randint(minrnd, maxrnd),
            randint(minrnd, maxrnd)]


def addrndtovert(vertlist: list,
                   minrnd: int,
                   maxrnd: int) -> list:
    """Adds random int values to a list of vertices

    Args:
        vertlist: list of vertices
        minrnd  : minimum random value
        maxrnd  : maximum random value
    Returns:
        list of vertices
    """
    return [addvect(pt, randomvect(minrnd, maxrnd))
                for pt in vertlist]


def adddimz(
         vlist2D: list[Number, Number],
         value: Number
             ) -> list[Number, Number,
                               Number]:
    """Add the z dimension in a list of (x, y) vertices

    Args:
        vlist2D: list of 2D vertices
                 [(x, y), ....]
        value  : constant z value

    Returns:
        list of 3D vertices
        [(x, y, z), ....]
    """
    return [[v[0], v[1], value]
              for v in vlist2D]


def anglebetween2Dlines(
        u: list[Number, Number],
        v: list[Number, Number]
        ) -> float:
    """Compute the angle between two lines of vectors

    Args:
        u, v: list[Number, Number]

    Returns:
        float angle in radians
    """
    return (
        atan(slope(u, v))
        if u[0] != v[0]
        else iif(u[0] < v[0], 1.5707963267948966, 4.71238898038469)
    )


def rotatebits(bits: int) -> int:
    """Rotates the bits in a byte

    Args:
        bits: the int 8 bits to rotate

    Returns:
        int value of rotated 8 bits
    """
    return sum(
        ((bits & (1 << bit)) >> bit) << (7 - bit) for bit in range(7, 0, -1)
    )


def mirror1stquad(x: Number, y: Number,
        v: list[Number, Number]
      ) -> list[list[Number, Number],
                list[Number, Number],
                list[Number, Number],
                list[Number, Number]]:
    """Mirrors a 2D vector in the 1st quadrant to create
    four symmetrical 2D vectors

    Args:
        x: x-axis of symmetry
        y: y-axis of symmetry
        v: vector (x:Number, y:Number)
           to mirror along x and y

    Returns:
        four symmetrical vectors (x, y)
    """
    xmin, xmax = mirror(x, v[0])
    ymin, ymax = mirror(y, v[1])
    return [[xmin, ymax], [xmax, ymax],
            [xmin, ymin], [xmax, ymin]]


def xorvect(u: list[int],
            v: list[int]) -> list[int]:
    """Applies a xor operation of between the elements of
    two lists of ints

    Args:
        v      : list[int]
        bitmask: int

    Returns:
        list[int]
    """
    return [i ^ j for i, j in zip(u, v)]


def andvect(u: list[int],
            v: list[int]) -> list[int]:
    """Bitwise and operation of between
    the elements of two lists of ints

    Args:
        v      : list[int]
        bitmask: int

    Returns:
        list[int]
    """
    return [i & j for i, j in zip(u, v)]


def bitmaskvect(v: list[int],
        bitmask: int) -> list[int]:
    """Applies a bitmask to a list of ints

    Args:
        v      : list[int]
        bitmask: int

    Returns:
        list[int]
    """
    return [b & bitmask for b in v]


def orvect(u: list[int],
           v: list[int]) -> list[int]:
    """Bitwise or operation between
    the elements of two lists of ints

    Args:
        v      : list[int]
        bitmask: int

    Returns:
        list[int]
    """
    return [i | j for i, j in zip(u, v)]


def addvectinlist(
        vlist: list[list[Number]]
               ) -> Number:
    """Gets the sum of the vectors in a list of vectors

    Args:
        vlist: list of vectors

    Returns:
        list or vector
    """
    return reduce(addvect,vlist)


def addvectpairlist(
    vpair: list[list[Number],
                list[Number]]
     )  -> list[Number]:
    """Adds pairwise the numbers in a list of
    two lists of numbers and returns a single
    list of numbers

    Args:
        vpair: list of two lists[Number]

    Returns:
        list[Number]
    """
    return addvect(vpair[0], vpair[1])


def addvecttripletlist(
    vtriplet: list[list[Number],
                   list[Number],
                   list[Number]]
         ) -> list[Number]:
    """Adds pairwise the numbers in a list of three lists
    of numbers and returns a single list of numbers

    Args:
        vtriplet: list[list[Number],
                       list[Number],
                       list[Number]]

    Returns:
        list[Number]
    """
    return addvect(
           addvect(vtriplet[0],
                   vtriplet[1]),
                   vtriplet[2])


def addvectlist(
        vlist1: list[list[Number]],
        vlist2: list[list[Number]]
           ) -> list[list[Number]]:
    """Adds pairwise the vectors (list of numbers) in a list
    then return a list of vectors

    Args:
        vlist1: First list of vectors
                list[lists[Number]]
        vlist2: Second list of vectors
                list[lists[Number]]

    Returns:
        A list of vectors
        list[list[Number]]
    """
    return [addvect(u, v)
            for u, v in zip(vlist1, vlist2)]


def swapxy(v:list) -> list:
    """Swaps the first two values in a list

    Args:
        list[x, y]

    Returns:
        list[y, x]
    """
    return [v[1], v[0]]


def centerpoint(x1: int, y1: int,
                x2: int, y2: int):
    """Returns the centerpoint x, y in rectangular area

    Args:
        x1, y1 : defines the
        x2, y2   rectangular area

    Returns:
        x: int, y: int centerpoint
    """
    return ((x2 - x1) >> 1) + x1, ((y2 - y1) >> 1) + y1


def enumbits(byteval: int):
    """Yields the 8 bits in a byte

    Args:
        byteval: int value
                 from 0 to 255

    Yields:
        Eight bits that is either
        int 0 or int 1
    """
    for bit in range(7, -1, -1):
        yield  ((byteval & (1 << bit)) >> bit)


def delta(v: list[Number, Number]):
    """Gets the delta from a list of two numbers

    Args:
        v: list[Number, Number]

    Returns:
        delta of numbers
    """
    return (v[1] - v[0])


def newtonmethod(x0: Number,
            f: Callable,
       fprime: Callable,
      maxiter: int = 40,
          tol: float = 1.e-8) -> tuple:
    """The Newton method for finding roots

    Args:
        x0     : initial guess
        f      : function
        fprime : derivative of f
        maxiter: maximum iterations
        tol    : tolerance

    Returns:
        Returns the root and iter if found or
        None if no convergence was reached

    """
    null = (None, None)
    x = x0
    for i in range(maxiter):
        fx = f(x)
        # Checks for Overflow / Underflow
        #if abs(fx) > 10e10 or abs(fx) < 10e-14:
        #    return null
        fpx = fprime(x)
        if fpx != 0:
            dx = fx / fpx
        else:
            return null
        if abs(dx) < tol:
            return (x, i)
        x -= dx
    return null


def tetration(x: Number, n: Number) -> Number:
    """ Tetration by recursion

    Args
        x: base
        n: exponent

    Returns:
        number
    """
    if n == 0:
        return 1
    return x ** tetration(x, n-1)


fftfn = lambda a: exp(a * -1j)
ifftfn = lambda a: exp(a * 1j)


def applyfft(fn : Callable, v: list[complex], tmp: list = []):
    """
    A recursive implementation of
    the 1D Cooley-Tukey FFT, the
    input should have a length of
    power of 2.

    Args:
       fn: function
        v: vector or list of complex numbers
      tmp: optional buffer (list)

    Returns:
        vector or list of complex numbers
    """
    n = len(v)
    if n > 1:
        l = n >> 1
        ve = v[::2]
        vo = v[1::2]
        applyfft(fn, ve, v)
        applyfft(fn, vo, v)
        for m in range(l):
            w = fn(tau * m / n) * vo[m]
            vm = ve[m]
            v[m] = vm + w
            v[m + l] = vm - w


def fft(v: list[complex]):
    """
    FFT

    A recursive implementation of
    the 1D Cooley-Tukey FFT, the
    input should have a length of
    power of 2.

    Args:
        v: vector or list of complex number

    Returns:
        vector or list of complex numbers
    """
    applyfft(fftfn, v)


def ifft(v: list[complex]):
    """
    Inverse FFT

    A recursive implementation of
    the 1D Cooley-Tukey FFT, the
    input should have a length of
    power of 2.

    Args:
        v: vector or list of complex numbers

    Returns:
        vector or list of complex numbers
    """
    applyfft(ifftfn, v)


def normvec(v: list[Number]) -> list[float]:
    """Normalize a vector

    To normalize a vector is to take a vector of
    any length and, keeping it pointing in the
    same direction, change its length to 1,
    turning it into what is called a unit vector.

    It describes a vector's direction without
    regard to its length.

    Args
        v: a list of Numbers

    Returns:
        a unit vector (list of floats)
    """
    return scalardivvect(v, vmag(v))


def getcomplexdomainbounds(v: list[Complex]) -> list[float, float, float, float]:
    """Get the bounds of a list of complex numbers
    in the complex plane

    Args
        v: a list of complex numbers

    Returns:
        (float: minreal, float: maxreal, float: minimag, float: maximag)
    """
    p = [z.real for z in v]
    q = [z.imag for z in v]
    return (min(p), max(p), min(q), max(q))


def getdomainbounds(v: list[list[float, float]]) -> list[float, float, float, float]:
    """Get the bounds of a list of xy coordinate pairs

    Args
        v: a list of [x, y]

    Returns:
        (float: minreal, float: maxreal, float: minimag, float: maximag)
    """
    p = [z[0] for z in v]
    q = [z[1] for z in v]
    return (min(p), max(p), min(q), max(q))


def initmatrix(i: int, j: int, initval: any = 0):
    """Creates a matrix with  i rows and j columns

    Args
        i: int rows of the matrix
        j: int cols of the matrix

    Returns:
        [[c: int:, c1: int, ...], ...]
    """
    return  [[initval for i in range(i)] for j in range(j)]


def histogram2D(
        v: list[float, float],
        xbins: int,
        ybins: int) -> list[list[int]]:
    """Computes the histogram of a list of x, y points

    Args
        v: a list of xy pairs

    Returns:
        [[c: int:, c1: int, ...], ...]
    """
    (xmin, xmax, ymin, ymax) = getdomainbounds(v)
    xbinsize = (xmax - xmin) / (xbins - 1)
    ybinsize = (ymax - ymin) / (ybins - 1)
    h = initmatrix(xbins, ybins)
    for z in v:
        h[int((z[1] - ymin) / ybinsize)][int((z[0] - xmin) / xbinsize)] += 1
    return h


def histogram2Dcomplex(
        v: list[Complex],
        pbins: int,
        qbins: int) -> list[list[int]]:
    """Computes the histogram of a list of points
    in the complex plane centered at 0, 0

    Args
        v: a list of complex numbers

    Returns:
        [[c: int:, c1: int, ...], ...]
    """
    (pmin, pmax, qmin, qmax) = getcomplexdomainbounds(v)
    pbinsize = (pmax - pmin) / (pbins - 1)
    qbinsize = (qmax - qmin) / (qbins - 1)
    h = initmatrix(pbins, qbins)
    for z in v:
        h[int((z.imag - qmin)/ qbinsize)][int((z.real - pmin)/ pbinsize)] += 1
    return h


def ravel2D(v: list [list[any]]) -> list[any]:
    """Flatten a nested list

    Args
        v: [[x: any:, x1: any, ...], ...]

    Returns:
        [x: any:, x1: any, ...]
    """
    return [i for sv in v for i in sv]


def lerp(a: float,
         b: float,
         f: float) -> float:
    """
    Calculates a number between two numbers at a specific increment.

    Args:
        a, b: float values
        f   : the amount to interpolate between the two values
              where 0.0 equal to the first point,
              0.1 is very near the first point
              0.9 is very near the second point
              etc

    Returns:
        a float value between a and b
    """
    return (1.0 - f) * a + f * b


def invlerp(a: float,
            b: float,
            v: float) -> float:
    """
    Returns a number between 0 and 1
    for a value (v) given that is
    between values a and b


    Args
        a, b: float values
        v   : number between a and b

    Returns:
        a float value between 0.0 and 1.0
    """
    d = b - a
    if d != 0:
        v = (v - a) / d
    else:
        v = 1
    return v


def clamp(x: float,
          a: float,
          b: float) -> float:
    """
    Returns a number between a and b
    for a value x given that is
    between values a and b


    Args
        x, a, b,: float values


    Returns:
        a float value between a and b
    """
    if x < a:
        x = a
    if x > b:
        x = b
    return x


def smoothstep(a: float,
               b: float,
               v: float) -> float:
    """
    Sigmoid-like interpolation

    Returns a number between 0.0 and 1.0
    for a value (v) given that is
    between values a and b


    Args
        a, b: float values
        v   : number between a and b

    Returns:
        a float value between 0.0 and 1.0
    """
    d = b - a
    if d != 0:
        v = (v - a) / d
    else:
        v = 1
    return v * v * (3 - 2 * v)


def smootherstep(a: float,
                 b: float,
                 v: float) -> float:
    """
    Sigmoid-like interpolation

    Returns a number between 0.0 and 1.0
    for a value (v) given that is
    between values a and b


    Args
        a, b: float values
        v   : number between a and b

    Returns:
        a float value between 0.0 and 1.0
    """
    d = b - a
    if d != 0:
        v = (v - a) / d
    else:
        v = 1
    return v * v * v * (v * (v * 6 - 15) + 10)


def smootheststep(a: float,
                 b: float,
                 v: float) -> float:
    """
    Sigmoid-like interpolation

    Returns a number between 0.0 and 1.0
    for a value (v) given that is
    between values a and b


    Args
        a, b: float values
        v   : number between a and b

    Returns:
        a float value between 0.0 and 1.0
    """
    d = b - a
    if d != 0:
        v = (v - a) / d
    else:
        v = 1
    return v * v * v * v * (v * (v * (-20 * v + 70) - 84) + 35)


def invsmoothstep(x: float):
    """Inverse Sigmoid-like interpolation

    Returns a number between 0.0 and 1.0
    for a value x given that is
    between 0.0 and 1.0

    Args
        x: float value between 0.0 and 1.0

    Returns:
        a float value between 0.0 and 1.0
    """
    return 0.5 - sin(asin(1.0 - 2.0 * x) / 3.0)


def smoothsteplerp(a: float,
         b: float,
         f: float) -> float:
    """
    Calculates a number between two numbers at a specific increment.
    with sigmoid like smoothing

    Args:
        a, b: float values
        f   : the amount to interpolate between the two values
              where 0.0 equal to the first point,
              0.1 is very near the first point
              0.9 is very near the second point
              etc

    Returns:
        a float value between a and b
    """
    return lerp(a, b, smoothstep(a, b, lerp(a, b, f)))


def smoothersteplerp(a: float,
         b: float,
         f: float) -> float:
    """
    Calculates a number between two numbers at a specific increment.
    with sigmoid like smoothing

    Args:
        a, b: float values
        f   : the amount to interpolate between the two values
              where 0.0 equal to the first point,
              0.1 is very near the first point
              0.9 is very near the second point
              etc

    Returns:
        a float value between a and b
    """
    return lerp(a, b, smootherstep(a, b, lerp(a, b, f)))


def smootheststeplerp(a: float,
         b: float,
         f: float) -> float:
    """
    Calculates a number between two numbers at a specific increment.
    with sigmoid like smoothing

    Args:
        a, b: float values
        f   : the amount to interpolate between the two values
              where 0.0 equal to the first point,
              0.1 is very near the first point
              0.9 is very near the second point
              etc

    Returns:
        a float value between a and b
    """
    return lerp(a, b, smootheststep(a, b, lerp(a, b, f)))