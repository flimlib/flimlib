/*-
 * #%L
 * FLIMLib package for exponential curve fitting of fluorescence lifetime data.
 * %%
 * Copyright (C) 2010 - 2022 University of Oxford and Board of Regents of the
 * University of Wisconsin-Madison.
 * %%
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/gpl-3.0.html>.
 * #L%
 */
#define PY_SSIZE_T_CLEAN
#include <Python.h>

/*
 * This is an empty Python module that allows us to use setuptools' build
 * system to build flimlib as a Python extension, thereby automatically taking
 * care of the platform and ABI tags in the DLL and wheel filenames.
 * Actual access to fliblib only happens through ctypes.
 */

static PyModuleDef module_def = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_flimlib",
    .m_doc = "Dummy module for flimlib C library",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit__flimlib(void)
{
    PyObject *m;

    m = PyModule_Create(&module_def);
    if (m == NULL)
        return NULL;

    return m;
}
