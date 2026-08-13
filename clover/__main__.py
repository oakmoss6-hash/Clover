"""Command-line entry point for Clover."""

import sys

from clover.load_config import HELP_TEXT
from clover.main import run_cli


def main():
    if any(arg in {"-h", "--help"} for arg in sys.argv[1:]):
        print(HELP_TEXT)
        return

    run_cli()


if __name__ == "__main__":
    main()
