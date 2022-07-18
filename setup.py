import glob
import os.path
import re
import sys
from distutils.command.build_ext import build_ext
from distutils.unixccompiler import UnixCCompiler
from xml.etree import ElementTree

import setuptools
from wheel.bdist_wheel import bdist_wheel


def read_pom_version():
    pom = ElementTree.parse("pom.xml")
    pom_root = pom.getroot()
    m = re.match(r"\{(.*)\}project", pom_root.tag)
    xmlns = m.group(1)
    assert xmlns.startswith("http://maven.apache.org/POM/")
    namespaces = {"pom": xmlns}
    pom_version = pom_root.findall("pom:version", namespaces)
    assert len(pom_version) == 1
    return pom_version[0].text.strip()


def get_pep440_version():
    pom_version = read_pom_version()
    m = re.fullmatch(r"([0-9]+)\.([0-9]+)\.([0-9]+)(-SNAPSHOT)?", pom_version)
    assert m, "pom version must match major.minor.patch[-SNAPSHOT]"
    major, minor, patch = m.group(1, 2, 3)
    is_snapshot = m.group(4) is not None
    suffix = ".dev0" if is_snapshot else ""
    return f"{major}.{minor}.{patch}{suffix}"


# Workaround for setuptools/distutils limitations described in:
# https://github.com/pypa/setuptools/issues/1192 and
# https://github.com/pypa/setuptools/issues/1732
# In our case, we need -std=c++14 (or c++11) for macOS clang++, but setuptools
# uses 'clang $CFLAGS' for both C and C++ sources. Passing a C++ std to the C
# compiler is an error, so we need a mechanism to pass a flag only when
# compiling C++ sources. (Credit to https://stackoverflow.com/a/65865696 for
# this idea.)

unix_cxx_flags = "-std=c++14"  # We could also read os.environ.get('CXXFLAGS')


class UnixCCxxCompiler(UnixCCompiler):
    def _compile(self, obj, src, ext, cc_args, extra_postargs, pp_opts):
        is_cxx_source = os.path.splitext(src)[-1] in (".cpp", ".cxx", ".cc")
        my_cc_args = [unix_cxx_flags] + cc_args if is_cxx_source else cc_args
        super()._compile(obj, src, ext, my_cc_args, extra_postargs, pp_opts)


class build_ext_c_cxx(build_ext):
    def build_extension(self, ext):
        if self.compiler.compiler_type == "unix":
            compiler = UnixCCxxCompiler()
            for attr, value in self.compiler.__dict__.items():
                setattr(compiler, attr, value)
            self.compiler = compiler
        super().build_extension(ext)


# normpath cleans up slashes on Windows
flimlib_c_sources = [
    os.path.normpath(p) for p in glob.glob("src/main/c/**/*.c", recursive=True)
]
flimlib_cpp_sources = [
    os.path.normpath(p)
    for p in glob.glob("src/main/c/**/*.cpp", recursive=True)
]

flimlib_msvc_def = "src/main/c/flimlib.def"

module_source = "src/main/python/flimlib/flimlib_dummy.c"

c_sources = flimlib_c_sources + flimlib_cpp_sources + [module_source]

# Technically we should make conditional on compiler, not OS, but we only
# support MSVC for now.
link_args = ["/DEF:" + flimlib_msvc_def] if sys.platform == "win32" else []


flimlib_ext = setuptools.Extension(
    "flimlib._flimlib",
    sources=c_sources,
    extra_link_args=link_args,
)


class bdist_wheel_abi3(bdist_wheel):
    # See https://github.com/joerick/python-abi3-package-sample
    def get_tag(self):
        python, abi, plat = super().get_tag()
        if python.startswith("cp"):
            return "cp36", "abi3", plat
        return python, abi, plat


setuptools.setup(
    version=get_pep440_version(),
    install_requires=["numpy>=1.12.0"],
    ext_modules=[flimlib_ext],
    py_limited_api=True,
    define_macros=[
        ("Py_LIMITED_API", "0x03060000"),
    ],
    cmdclass={
        "build_ext": build_ext_c_cxx,
        "bdist_wheel": bdist_wheel_abi3,
    },
)
