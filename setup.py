import distutils.command.build
import distutils.spawn
import glob
import os
import sys

import setuptools
import setuptools.command.build_py
import setuptools.command.build_ext


class build_mvn(distutils.cmd.Command):
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

class build_ext(setuptools.command.build_ext.build_ext):
    def run(self):
        self.run_command('build_mvn')
        return super().run()

_flimlib_dummy = setuptools.Extension(
    'flimlib._flimlib_dummy',
    sources=[
        'src/main/python/flimlib/flimlib_dummy.c',
    ]
)

distutils.command.build.build.sub_commands.append(
    ('build_mvn', build_mvn.should_run_mvn))

setuptools.setup(
    install_requires=["numpy>=1.12.0"],
    cmdclass={'build_mvn': build_mvn, 'build_ext': build_ext},
    # Empty string means target directory
    #data_files=[(os.path.join('lib', 'site-packages', 'flimlib'),
    #             [os.path.join('target', 'natives', 'flimlib.dll')])],
    data_files=[('flimlib',
                 [os.path.join('target', 'natives', 'flimlib.dll')])],
    include_package_data=True,
    ext_modules=[_flimlib_dummy]
)
