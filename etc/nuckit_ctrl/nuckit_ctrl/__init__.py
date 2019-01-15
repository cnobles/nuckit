__author__ = "Christopher Nobles"
__license__ = "GPL3+"

from os import getenv, getcwd, path
from subprocess import run, PIPE

def get_nuckit_version(with_hash = False):
    nuckit_version_path = getenv("NUCKIT_DIR", "None")
    
    if nuckit_version_path in ["None"]:
        raise SystemExit(
          print("  NUCKIT_DIR cannot be found as an environmental variable.\n"
                "  Check to make sure your NucKit environment is active,   \n"
                "  you may need to restart your environment, update, or    \n"
                "  reinstall nuckit with the install.sh script.")
        )
    else:
        nuckit_version_path = nuckit_version_path + "/.version"
    
    if not path.exists(nuckit_version_path):
        raise SystemExit(
          print("  NucKit version cannot be located. Check environmental\n"
                "  variables, such as NUCKIT_DIR, otherwise you may want\n"
                "  to restart your environment, update, or reinstall    \n"
                "  NucKit using the install.sh script.")
        )
    
    nuckit_version = open(nuckit_version_path, "r").readlines()[0].rstrip()
    
    commit_hash = run(
      ["git", "rev-parse", "--short", "HEAD"], stdout=PIPE
    )
    
    commit_str = commit_hash.stdout.decode('utf-8').rstrip()
    
    if with_hash:
        return nuckit_version + "+" + commit_str
    else:
        return nuckit_version

__version__ = get_nuckit_version( with_hash = True )
