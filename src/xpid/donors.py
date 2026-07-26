"""Compatibility wrapper for :mod:`xpid.rotatable_groups`."""

from .rotatable_groups import (  # noqa: F401
    RotatableGroupDefinition as DonorDefinition,
    RotatableGroupKind as DonorKind,
    ROTATABLE_GROUPS as ROTATABLE_DONORS,
    find_parent_atom as parent_atom,
    get_rotatable_group as get_definition,
    is_rotatable,
)
