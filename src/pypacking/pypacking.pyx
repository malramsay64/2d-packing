# This is a wrapper around the 2d-packing library for python

cdef extern from "wallpaper.h":
    ctypedef struct shapetype:
        pass

cdef class Shape:

    def __init__(self, name: str, points: np.ndarray):
        self.name = name
        self.resolution = len(points)


cpdef struct shapetype define_shape():
