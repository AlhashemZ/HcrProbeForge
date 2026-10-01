"""Installable HCRProbeForge package.

The scientific implementation lives in :mod:`hcrprobeforge.core`.  Keeping
that module intact and routing every front end through ``core.main`` preserves
the established HCRProbeDesign, Primer3, adaptive-curation, and reporting
behavior while adding installable entry points.
"""

from .core import SCRIPT_BUILD, __version__

__all__ = ["SCRIPT_BUILD", "__version__"]
