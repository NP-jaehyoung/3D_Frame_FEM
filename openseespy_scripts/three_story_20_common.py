"""Shared constants and small helpers for 20-story 3D OpenSeesPy examples."""

from __future__ import annotations

import numpy as np


NUM_STORIES = 20
FLOOR_HEIGHT = 4.0
SPAN_X = 4.0
SPAN_Z = 4.0

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


def top_node_id(num_levels: int) -> int:
    """Return node id of top node (corner 4 on top level)."""
    return 4 * num_levels


def floor_loads(
    num_levels: int,
    span_x: float = SPAN_X,
    span_z: float = SPAN_Z,
    pressure_kpa: float = FLOOR_LOAD_PRESSURE_KPA,
    include_base: bool = INCLUDE_BASE_FLOOR_LOAD,
) -> np.ndarray:
    """Return per-level distributed floor load (total force by level, kN)."""
    loads = np.full(num_levels, pressure_kpa * span_x * span_z, dtype=float)
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

