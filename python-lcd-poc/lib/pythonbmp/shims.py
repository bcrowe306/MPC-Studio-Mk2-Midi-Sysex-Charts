"""
conftest.py - builtin classes generics shim for python < 3.9
"""
from sys import version_info

if version_info < (3, 9):
    # For each of the following types, add a
    # new class function called `__class_getitem__`.
    # This is called when the type is subscripted,
    # such as in:
    #
    #    list[str]
    #    dict[int, list[int]]
    #
    # This will make these typing expressions work
    # seamlessly on python versions >= 3.7 and < 3.9
    # (support was added, but not used out-of-the-box),
    # in python 3.7).
    import gc, types, typing
    from typing import _GenericAlias as GenericAlias
    for t in (list, dict, set, tuple, frozenset):
      r = gc.get_referents(t.__dict__)[0]
      r.update({
        "__class_getitem__": classmethod(GenericAlias),
      })

