from setuptools import setup

setup(
    name = 'nuckit_ctrl',
    packages = [ 'nuckit_ctrl', 'nuckit_ctrl.scripts'  ],
    entry_points = { 'console_scripts' : [
        'nuc = nuckit_ctrl.scripts.command:main'
        ] },
)
