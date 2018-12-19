import os
import sys
import argparse
import subprocess

def main(argv = sys.argv):

    usage_str = "%(prog)s [-h/--help, -v/--version] <subcommand>"
    description_str = (
        "subcommands:\n"
        "  demulti   \tDemultiplex Illumina sequence files given various barcoding schemes.\n"
        "  trim      \tTrim leading and/or over-reading sequences from sequence files.\n"
        "  filt      \tFilter sequences in one or more sequence files by various criteria.\n"
        "  consol    \tConsolidate reads to unique sequences and generate a key file.\n"
        "  couple    \tCouple independent BLAT alignments from paired-end sequences.\n"
    )

    parser = argparse.ArgumentParser(
        prog = "nuckit",
        usage = usage_str,
        description = description_str,
        epilog = "For more help, see the docs at http://nuckit.readthedocs.io.",
        formatter_class = argparse.RawDescriptionHelpFormatter,
        add_help = False
    )

    parser.add_argument("command", help = argparse.SUPPRESS, nargs = "?")

    parser.add_argument(
        "--nuckit_dir", default = os.getenv("NUCKIT_DIR", os.getcwd()),
        help = "Path to nuckit installation. (default: %(default)s)'"
    )

#    parser.add_argument(
#        "-v", "--version", action = "version",
#        version = "%(prog)s v{}".format(nuckit_ctrl.__version__)
#    )

    args, remaining = parser.parse_known_args(argv)

    sub_cmds = ["demulti", "trim", "filt", "consol", "couple"]
    if not args.command in sub_cmds:
        parser.print_help()
        sys.stderr.write("Unrecognized subcommand, '{}'.\n".format(
            args.command
        ))
        sys.exit(1)

    r_script = Path(args.nuckit_dir) + "/scripts/" + args.command + ".R"
    if not r_script.exists():
        sys.stderr.write(
            "Error: Could not find a {0} in directory '{1}'\n".format(
                (args.command + ".R"), args.nuckit_dir)
        )
        sys.exit(1)

    r_comps = ["Rscript", r_script] + remaining

#    cmd = subprocess.run(r_comps, shell = True)

#    sys.exit(cmd.returncode)
    print(' '.join(r_comps))
    sys.exit(0)
