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