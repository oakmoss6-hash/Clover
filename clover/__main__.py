"""Command-line entry point for Clover."""

import runpy


if __name__ == "__main__":
    runpy.run_module(
        "clover.main",
        run_name="__main__",
    )
