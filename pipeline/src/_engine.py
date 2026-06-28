"""
_engine.py — bridge to the verified physics single source of truth.

The whole `src/` pipeline delegates its physics to `Skills/sei_physics.py`, so the
regenerated tables/figures/CSV can never diverge from the theory or the skills.
This module only resolves the import path and re-exports the physics module `P`.
"""

from __future__ import annotations

import sys
from pathlib import Path


def _find_skills() -> Path:
    """Search upward for the sibling Skills/ folder (location-independent)."""
    for up in Path(__file__).resolve().parents:
        cand = up / "Skills" / "sei_physics.py"
        if cand.exists():
            return cand.parent
    raise ImportError("Could not locate Skills/sei_physics.py above this file.")


_SKILLS = _find_skills()
if str(_SKILLS) not in sys.path:
    sys.path.insert(0, str(_SKILLS))

import sei_physics as P  # noqa: E402,F401  (re-exported)
