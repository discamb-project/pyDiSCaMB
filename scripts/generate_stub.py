# Script to generate a stub file for the C++ code

from pathlib import Path
import subprocess


def main():
    # Check dependency

    import pydiscamb
    import mypy

    root = Path(__file__).parent.parent

    subprocess.call(
        f"stubgen -m pydiscamb._cpp_module --include-docstrings -o {root}".split()
    )


if __name__ == "__main__":
    main()
