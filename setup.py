import glob
import os.path
import setuptools
import sys


# normpath cleans up slashes on Windows
flimlib_c_sources = [os.path.normpath(p) for p in
                     glob.glob('src/main/c/**/*.c', recursive=True)]
flimlib_cpp_sources = [os.path.normpath(p) for p in
                       glob.glob('src/main/c/**/*.cpp', recursive=True)]

flimlib_msvc_def = 'src/main/c/flimlib.def'

module_source = 'src/main/python/flimlib/flimlib_dummy.c'

c_sources = flimlib_c_sources + flimlib_cpp_sources + [module_source]

# Technically we should make conditional on compiler, not OS, but we only
# support MSVC for now.
link_args = ['/DEF:' + flimlib_msvc_def] if sys.platform == 'win32' else []


flimlib_ext = setuptools.Extension(
    'flimlib._flimlib',
    sources=c_sources,
    extra_link_args=link_args,
)

setuptools.setup(
    install_requires=["numpy>=1.12.0"],
    ext_modules=[flimlib_ext],
)
