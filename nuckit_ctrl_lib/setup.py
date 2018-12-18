from setuptools import setup

setup(name = 'nuckit_ctrl_lib',
      packages = [ 'nuckit_ctrl_lib', 'nuckit_ctrl_lib.scripts' ],
      entry_points = { 'console_scripts' : [
          'nuc = nuckit_ctrl_lib.scripts.command:main'
          ] }
     )
