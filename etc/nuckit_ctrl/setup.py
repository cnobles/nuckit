from os import getenv, getcwd
from subprocess import run, PIPE
from setuptools import setup, find_packages

def get_nuckit_version(with_hash = False):
    nuckit_version_path = getenv("NUCKIT_DIR", getcwd()) + "/.version"
    nuckit_version = open(nuckit_version_path, "r").readlines()[0].rstrip()
    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
      )
    commit_str = commit_hash.stdout.decode('utf-8').rstrip()
    if with_hash:
        return nuckit_version + "+" + commit_str
    else:
        return nuckit_version

setup(
    name = 'nuckit_ctrl',
    version = get_nuckit_version(),
    packages = find_packages(),
    entry_points = { 'console_scripts' : [
        'nuc = nuckit_ctrl.scripts.command:main',
        'nuckit = nuckit_ctrl.scripts.command:main'
    ] }
)
