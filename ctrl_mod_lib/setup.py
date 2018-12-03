from setuptools import setup

setup(name = 'ctrl_mod_lib',
      packages = [ 'ctrl_mod_lib', 'ctrl_mod_lib.scripts' ],
      entry_points = { 'console_scripts' : [
          'nuc = ctrl_mod_lib.scripts.command:main'
          ] }
     )
