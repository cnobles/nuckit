import os
import sys
import argparse
import subprocess

from pathlib import Path
from nuckit_ctrl import __version__

def main(argv = sys.argv):

    usage_str = "%(prog)s [-h/--help, -v/--version] <subcommand>"
    description_str = (
        "NucKit - Nucleotide sequence analysis took kit!\n"
        "\n"
        "subcommands:\n"
        "  demulti   \tDemultiplex Illumina sequence files given various barcoding schemes.\n"
        "  trim      \tTrim leading and/or over-reading sequences from sequence files.\n"
        "  filt      \tFilter sequences in one or more sequence files by various criteria.\n"
        "  consol    \tConsolidate reads to unique sequences and generate a key file.\n"
        "  couple    \tCouple independent BLAT alignments from paired-end sequences.\n"
    ).format(version=__version__)

    parser = argparse.ArgumentParser(
        prog = "nuckit",
        usage = usage_str,
        description = description_str,
        epilog = "For more help, see the docs at http://nuckit.readthedocs.io.",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        add_help = False
    )

    parser.add_argument(
        "command", help = argparse.SUPPRESS, nargs = "?"
    )

    parser.add_argument(
        "subcommand", help = argparse.SUPPRESS, default = 'None', nargs = "?"
    )

    parser.add_argument(
        "--nuckit_dir", default = os.getenv("NUCKIT_DIR", os.getcwd()),
        help = "Path to NucKit installation. (default: %(default)s)"
    )

    parser.add_argument(
        "-v", "--version", action = "version",
        version = "%(prog)s {}".format(__version__)
    )

    args, remaining = parser.parse_known_args(argv)

    sub_cmds = ["demulti", "trim", "filt", "consol", "couple"]
    if not args.subcommand in sub_cmds:
        parser.print_help()
        if not args.subcommand in ['None']:
            sys.stderr.write("Unrecognized subcommand, '{}'.\n".format(
                args.subcommand
            ))
        sys.exit(1)

    r_script = Path(args.nuckit_dir + "/scripts/" + args.subcommand + ".R")
    if not r_script.is_file():
        sys.stderr.write(
            "Error: Could not find a {0} in directory '{1}'\n".format(
                (args.subcommand + ".R"), args.nuckit_dir)
        )
        sys.exit(1)

    r_comps = ["Rscript", str(r_script)] + remaining

    cmd = subprocess.run(r_comps)

    sys.exit(cmd.returncode)
