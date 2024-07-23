"""
 *Do not edit this file it contains
  bitmap file Address locations ...*
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

bmpheaderid = array('B', [66, 77])
bmpfilesize = 2
bmphdrsize = 10
bmpcolorbits = 28
bmpx = 18
bmpy = 22
bmppal = 54
bmpheadersizedict = \
    {1: 62,
     4: 118,
     8: 1078,
     24: 54}
