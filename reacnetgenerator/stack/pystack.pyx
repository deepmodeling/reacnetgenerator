# distutils: language = c++
# distutils: sources = c_stack.cpp
# cython: language_level=3
"""This file is copied from https://zhuanlan.zhihu.com/p/38212302"""

from cpython.ref cimport PyObject,Py_INCREF,Py_DECREF

cdef extern from 'c_stack.h':
    cdef cppclass C_Stack:
        PyObject* peek();

        void push(PyObject* val);

        PyObject* pop();

class StackEmpty(Exception):
    pass

cdef class Stack:
    cdef C_Stack _c_stack

    cpdef object peek(self):
        cdef PyObject* val
        val=self._c_stack.peek()
        if val==NULL:
            raise StackEmpty
        return <object>val

    cpdef object push(self,object val):
        Py_INCREF(val);
        self._c_stack.push(<PyObject*>val);
        return None

    cpdef object pop(self):
        cdef PyObject* val
        val=self._c_stack.pop()
        if val==NULL:
            raise StackEmpty
        cdef object rv=<object>val;
        Py_DECREF(rv)
        return rv