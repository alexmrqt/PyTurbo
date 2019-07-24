# distutils: language = c++

# Copyright 2019 Free Software Foundation, Inc.
#
# This is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3, or (at your option)
# any later version.
#
# This software is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this software; see the file COPYING.  If not, write to
# the Free Software Foundation, Inc., 51 Franklin Street,
# Boston, MA 02110-1301, USA.


from libcpp.vector cimport vector

cdef extern from "viterbi.cc":
    pass

cdef extern from "viterbi.h":
    cppclass viterbi:
        viterbi(int, int, int, vector[int], vector[int]) except +
        void viterbi_algorithm(int K, int S0, int, const float*, unsigned char*)
        int get_I()
        int get_S()
        int get_O()

cdef extern from "log_bcjr.cc":
    pass

cdef extern from "log_bcjr.h":
    cppclass log_bcjr:
        log_bcjr(int, int, int, vector[int], vector[int]) except +
        void log_bcjr_algorithm(vector[float], vector[float], vector[float], vector[float])
        int get_I()
        int get_S()
        int get_O()

import numpy

cdef class PyViterbi:
    cdef int I, S, O
    cdef viterbi* cpp_viterbi

    def __cinit__(self, int I, int S, int O, vector[int] NS, vector[int] OS):
        self.cpp_viterbi = new viterbi(I, S, O, NS, OS)
        self.I = self.cpp_viterbi.get_I()
        self.S = self.cpp_viterbi.get_S()
        self.O = self.cpp_viterbi.get_O()

    def __dealloc__(self):
        del self.cpp_viterbi

    def viterbi_algorithm(self, S0, SK, float[::1] _in):
        cdef int K = _in.shape[0]/self.O
        cdef unsigned char[::1] _out = numpy.zeros(K, dtype=numpy.uint8)

        self.cpp_viterbi.viterbi_algorithm(K, S0, SK, &_in[0], &_out[0])

        return numpy.asarray(_out)

cdef class PyLogBCJR:
    cdef int I, S, O
    cdef log_bcjr* cpp_log_bcjr

    def __cinit__(self, int I, int S, int O, vector[int] NS, vector[int] OS):
        self.cpp_log_bcjr= new log_bcjr(I, S, O, NS, OS)
        self.I = self.cpp_log_bcjr.get_I()
        self.S = self.cpp_log_bcjr.get_S()
        self.O = self.cpp_log_bcjr.get_O()

    def __dealloc__(self):
        del self.cpp_log_bcjr

    def viterbi_algorithm(self, S0, SK, float[::1] _in):
        cdef int K = _in.shape[0]/self.O
        cdef vector[float] _in_vec = _in.tolist()
        cdef vector[float] _out = numpy.zeros(self.S*self.I*K, dtype=numpy.float32)

        self.cpp_viterbi.viterbi_algorithm(K, S0, SK, _in_vec, _out)

        return numpy.asarray(_out)
