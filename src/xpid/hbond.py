"""Compatibility wrapper for :mod:`xpid.hbond_annotations`."""

from .hbond_annotations import *  # noqa: F401,F403
from .hbond_annotations import (
    annotate_cone_hbond_context as annotate_hbond_competition,
    annotate_explicit_hbond_context as annotate_rigid_hbond,
)
