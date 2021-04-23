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
        old_dir = os.getcwd()
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        try:
            distutils.spawn.spawn([maven])
        finally:
            os.chdir(old_dir)
    def should_run_mvn(self):
        return True

distutils.command.build.build.sub_commands.append(('build_clib',build_clib.should_run_mvn))

setuptools.setup(
    install_requires=["numpy>=1.12.0"],
    cmdclass={'build_clib': build_clib},
    # Empty string means target directory
    data_files=[("", [os.path.join('target', 'natives', 'flimlib.dll')])],
    )
