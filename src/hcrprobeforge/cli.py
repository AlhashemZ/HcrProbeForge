"""Console entry point for the HCRProbeForge CLI."""

from __future__ import annotations

import sys
from typing import Sequence

from . import core, references


def main(argv: Sequence[str] | None = None) -> int:
    """Execute the compatibility CLI with an optional argument list."""
    values = list(argv) if argv is not None else sys.argv[1:]
    if values and values[0] == "index":
        return references.index_main(values[1:])
    if values and values[0] == "species":
        return references.species_main(values[1:])
    return core.main(values)


if __name__ == "__main__":
    raise SystemExit(main())
