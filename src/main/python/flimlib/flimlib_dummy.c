#define PY_SSIZE_T_CLEAN
#include <Python.h>

static PyModuleDef dummy_module = {
    PyModuleDef_HEAD_INIT,
    .m_name = "_flimlib_dummy",
    .m_doc = "dummy module for flimlib",
    .m_size = -1,
};

PyMODINIT_FUNC
PyInit__flimlib_dummy(void)
{
    PyObject *m;

    m = PyModule_Create(&dummy_module);
    if (m == NULL)
        return NULL;

    return m;
}