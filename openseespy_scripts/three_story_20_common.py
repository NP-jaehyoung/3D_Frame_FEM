"""Shared constants and small helpers for 20-story 3D OpenSeesPy examples."""

from __future__ import annotations

import numpy as np


NUM_STORIES = 20
FLOOR_HEIGHT = 4.0
SPAN_X = 4.0
SPAN_Z = 4.0
NUM_BAYS_X = 3
NUM_BAYS_Z = 3

POINT_LOAD_VALUE = -3000.0
POINT_LOAD_DIR = "UX"

FLOOR_LOAD_PRESSURE_KPA = -5.0
FLOOR_LOAD_DIR = "UY"
INCLUDE_BASE_FLOOR_LOAD = False

E_DEFAULT = 2.0e7
NU_DEFAULT = 0.167
WIDTH_DEFAULT = 1.0
DEPTH_DEFAULT = 1.0

DT_DEFAULT = 0.01
TIME_END_DEFAULT = 20.0
AG_PEAK_G_DEFAULT = 0.30
AG_FREQ_DEFAULT = 1.8
AG_CENTER_DEFAULT = 5.0
AG_WIDTH_DEFAULT = 1.1
DAMPING_RATIO_DEFAULT = 0.02
DEAD_LOAD_G_DEFAULT = 9.81

INCLUDE_DEAD_LOAD_AS_MASS = True

STATIC_PREFIX = "three_story_20_static"
STATIC_PREFIX_OPENSEES = "three_story_20_opensees"
DYNAMIC_PREFIX = "three_story_20_dynamic"
DYNAMIC_PREFIX_OPENSEES = "three_story_20_opensees"


def level_count(num_stories: int = NUM_STORIES) -> int:
    """Return number of levels including ground level."""
    if num_stories < 0:
        raise ValueError("num_stories must be zero or positive.")
    return num_stories + 1


def nodes_per_level(
    num_bays_x: int = NUM_BAYS_X,
    num_bays_z: int = NUM_BAYS_Z,
) -> int:
    """Return node count per floor level for the rectangular bay grid."""
    if num_bays_x < 1 or num_bays_z < 1:
        raise ValueError("num_bays_x and num_bays_z must be positive.")
    return (num_bays_x + 1) * (num_bays_z + 1)


def floor_node_id(
    level_index: int,
    ix: int,
    iz: int,
    num_bays_x: int = NUM_BAYS_X,
    num_bays_z: int = NUM_BAYS_Z,
) -> int:
    """Return 1-based node id for a zero-based level/grid location."""
    if level_index < 0:
        raise ValueError("level_index must be zero or positive.")
    if not (0 <= ix <= num_bays_x and 0 <= iz <= num_bays_z):
        raise ValueError("Grid index is outside the bay range.")
    return level_index * nodes_per_level(num_bays_x, num_bays_z) + iz * (num_bays_x + 1) + ix + 1


def top_node_id(
    num_levels: int,
    num_bays_x: int = NUM_BAYS_X,
    num_bays_z: int = NUM_BAYS_Z,
) -> int:
    """Return node id of the roof corner at max X and Z=0."""
    if num_levels < 1:
        raise ValueError("num_levels must be at least 1.")
    return floor_node_id(num_levels - 1, num_bays_x, 0, num_bays_x, num_bays_z)


def floor_loads(
    num_levels: int,
    span_x: float = SPAN_X,
    span_z: float = SPAN_Z,
    num_bays_x: int = NUM_BAYS_X,
    num_bays_z: int = NUM_BAYS_Z,
    pressure_kpa: float = FLOOR_LOAD_PRESSURE_KPA,
    include_base: bool = INCLUDE_BASE_FLOOR_LOAD,
) -> np.ndarray:
    """Return per-level distributed floor load (total force by level, kN)."""
    floor_area = (num_bays_x * span_x) * (num_bays_z * span_z)
    loads = np.full(num_levels, pressure_kpa * floor_area, dtype=float)
    if not include_base and num_levels > 0:
        loads[0] = 0.0
    return loads


def section_properties(
    width: float = WIDTH_DEFAULT,
    depth: float = DEPTH_DEFAULT,
    E: float = E_DEFAULT,
    nu: float = NU_DEFAULT,
):
    """Return A, G, J, Iy, Iz for the sample square/rectangular section."""
    if width <= 0 or depth <= 0:
        raise ValueError("width and depth must be positive.")
    A = width * depth
    Iy = width * depth**3 / 12.0
    Iz = depth * width**3 / 12.0
    if width >= depth:
        ar = depth / width
        J = depth * width**3 * (1 / 3 - 0.21 * ar + 0.063 * ar**5)
    else:
        ar = width / depth
        J = width * depth**3 * (1 / 3 - 0.21 * ar + 0.063 * ar**5)
    if J <= 0.0:
        raise ValueError("Invalid torsional constant J.")
    G = E / (2.0 * (1.0 + nu))
    return A, G, J, Iy, Iz

