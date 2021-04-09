import setuptools
import glob
import os
import sys
import distutils.command.build
import distutils.spawn


class build_clib(distutils.cmd.Command):
    user_options = []

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        if sys.platform == 'win32':
            try:
                mvn_home = os.environ['MAVEN_HOME']
            except KeyError:
                raise Exception('MAVEN_HOME not set')
            if not mvn_home:
                raise Exception('MAVEN_HOME not set')
            maven = os.path.join(mvn_home, 'mvn.cmd')

        else:
            maven = distutils.spawn.find_executable('mvn')
            if maven is None:
                raise FileNotFoundError("unable to find mvn")
        distutils.spawn.spawn([maven])


setuptools.setup(
    install_requires=["numpy>=1.12.0"],
    cmdclass={'build_clib': build_clib},
    data_files=[(os.path.dirname(os.path.realpath(__file__)), [os.path.join('target', 'natives', 'flimlib.dll')])])
