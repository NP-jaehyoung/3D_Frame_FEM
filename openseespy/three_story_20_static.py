"""OpenSeesPy Static Analysis for 20-story 3D rigid-diaphragm frame.

- 20 floors, 4 nodes per floor (4m x 4m square plan, 4m story height)
- 6 DOF node (UX,UY,UZ,RX,RY,RZ)
- Rigid diaphragm enforced per floor
- Loads:
    - Gravity: Distributed floor load (5 kN/m^2)
    - Lateral: Point load at top (-3000 kN UX)
- Output:
    - Static Displacements
    - Base Reactions
    - Natural Frequencies (Modal Analysis)
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from three_story_20_common import (
    FLOOR_HEIGHT,
    DEAD_LOAD_G_DEFAULT,
    FLOOR_LOAD_DIR,
    FLOOR_LOAD_PRESSURE_KPA,
    INCLUDE_BASE_FLOOR_LOAD,
    INCLUDE_DEAD_LOAD_AS_MASS,
    NUM_STORIES,
    POINT_LOAD_DIR,
    POINT_LOAD_VALUE,
    SPAN_X,
    SPAN_Z,
    floor_loads,
    level_count,
    section_properties,
    top_node_id,
)

try:
    import openseespy.opensees as ops
except ImportError as err:
    raise RuntimeError(
        "OpenSeesPy is required: pip install openseespy"
    ) from err


def _build_geometry(num_floors: int, floor_height: float, span_x: float, span_z: float):
    """Build square 4-node floor mesh and element connectivity."""
    node_id = 1
    node_xy = []
    floor_nodes = []

    for i_floor in range(num_floors):
        y = floor_height * i_floor
        floor = []
        for x, z in ((0.0, 0.0), (0.0, span_z), (span_x, span_z), (span_x, 0.0)):
            ops.node(node_id, x, y, z)
            floor.append(node_id)
            node_xy.append((node_id, x, y, z))
            node_id += 1
        floor_nodes.append(floor)

    ele_id = 1
    for i_floor in range(num_floors - 1):
        base = 4 * i_floor
        # columns
        for a in range(4):
            n1 = base + a + 1
            n2 = n1 + 4
            ops.element("elasticBeamColumn", ele_id, n1, n2,
                        A, E, G, J, Iy, Iz, 1)
            ele_id += 1
        # floor beams (closed square) at TOP of the story
        # This will create beams at levels 1 to num_floors-1 (i.e. Floor 1 to Roof)
        # Base (Floor 0) will have no beams.
        base_top = base + 4
        n1 = base_top + 1
        n2 = base_top + 2
        n3 = base_top + 3
        n4 = base_top + 4
        # (n1,n2)=Z, (n2,n3)=X, (n3,n4)=Z, (n4,n1)=X
        # Use transf 2 for Z, 1 for X
        beam_pairs = [(n1, n2, 2), (n2, n3, 1), (n3, n4, 2), (n4, n1, 1)]
        for (i, j, transf) in beam_pairs:
            ops.element("elasticBeamColumn", ele_id, i, j,
                        A, E, G, J, Iy, Iz, transf)
            ele_id += 1

    return node_xy, floor_nodes


def _apply_rigid_diaphragms(num_floors: int, floor_nodes):
    for i_floor in range(1, num_floors):
        # Floor normal is global Y axis (2), so rigid diaphragm keeps X/Y/Z rigidly coupled.
        # This captures in-plane rigid-body deformation including torsion-like in-plane rotation.
        floor = floor_nodes[i_floor]
        master = floor[0]
        slave_nodes = floor[1:]
        if slave_nodes:
            ops.rigidDiaphragm(2, master, *slave_nodes)


def _add_nodal_loads(num_floors: int, floor_nodes, distributed_floor_loads, distributed_floor_dir, point_load_node, point_load_value, point_load_dir):
    component_map = {"UX": 1, "UY": 2, "UZ": 3, "RX": 4, "RY": 5, "RZ": 6}

    # point load (single node)
    p_dir_idx = component_map[point_load_dir]
    point_load = np.zeros(6)
    point_load[p_dir_idx - 1] = point_load_value
    ops.load(point_load_node, *point_load)

    # distributed floor load: split by floor nodes
    f_dir_idx = component_map[distributed_floor_dir]
    if distributed_floor_loads is not None:
        for i_floor, load in enumerate(distributed_floor_loads, start=1):
            if abs(load) <= 1e-12:
                continue
            nodes = floor_nodes[i_floor - 1]
            if not nodes:
                continue
            share = load / len(nodes)
            for n in nodes:
                vals = np.zeros(6)
                vals[f_dir_idx - 1] = share
                ops.load(n, *vals)


def _dof_to_name(dof: int) -> str:
    return ["UX", "UY", "UZ", "RX", "RY", "RZ"][dof - 1]


def _write_static_results(out_dir: Path, top_node: int, n_floors: int, floor_nodes):
    # nodal displacement table
    num_nodes = 4 * n_floors
    with open(out_dir / "three_story_20_static_displacements.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Node", "DOF", "Component", "Displacement"])
        for n in range(1, num_nodes + 1):
            for dof in range(1, 7):
                writer.writerow([n, dof, _dof_to_name(dof), ops.nodeDisp(n, dof)])

    # support reactions
    with open(out_dir / "three_story_20_static_reactions.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["SupportNode", "DOF", "Component", "Reaction"])
        for n in range(1, 5):
            for dof in range(1, 7):
                writer.writerow([n, dof, _dof_to_name(dof), ops.nodeReaction(n, dof)])

    # top-node displacement summary
    top_disp = [ops.nodeDisp(top_node, 1), ops.nodeDisp(top_node, 2), ops.nodeDisp(top_node, 3)]
    print(f"\nTop node {top_node} displacement: UX={top_disp[0]:.6e}, UY={top_disp[1]:.6e}, UZ={top_disp[2]:.6e}")


def run_static_analysis():
    # -----------------------------
    # model settings
    # -----------------------------
    num_stories = NUM_STORIES
    num_floors = level_count(num_stories)
    floor_height = FLOOR_HEIGHT
    span_x = SPAN_X
    span_z = SPAN_Z

    # material + section
    global E, G, J, A, Iy, Iz
    A, G, J, Iy, Iz = section_properties()

    # loads
    point_load_node = top_node_id(num_floors)
    point_load_value = POINT_LOAD_VALUE
    point_load_dir = POINT_LOAD_DIR
    
    # 5 kN/m^2 distributed load (total load per floor = pressure * area)
    # Applied in negative Y direction (gravity)
    # Floor 0 (Base) has NO load by default.
    distributed_floor_loads = floor_loads(
        num_floors,
        span_x=span_x,
        span_z=span_z,
        pressure_kpa=FLOOR_LOAD_PRESSURE_KPA,
        include_base=INCLUDE_BASE_FLOOR_LOAD,
    )
    distributed_floor_dir = FLOOR_LOAD_DIR
    
    include_dead_load_as_mass = INCLUDE_DEAD_LOAD_AS_MASS
    floor_dead_load_kilonewton = distributed_floor_loads

    # -----------------------------
    # model build
    # -----------------------------
    script_dir = Path(__file__).resolve().parent
    out_dir = script_dir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    ops.wipe()
    ops.model("Basic", "-ndm", 3, "-ndf", 6)
    
    # geom transformation for 3D linear elements
    # 1: vecxz along global Z (for columns and X-beams)
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    # 2: vecxz along global X (for Z-beams)
    ops.geomTransf("Linear", 2, 1.0, 0.0, 0.0)

    nodes, floor_nodes = _build_geometry(num_floors, floor_height, span_x, span_z)
    if len(nodes) != 4 * num_floors:
        raise RuntimeError("Node count mismatch.")

    # base support
    for n in floor_nodes[0]:
        ops.fix(n, 1, 1, 1, 1, 1, 1)

    _apply_rigid_diaphragms(num_floors, floor_nodes)

    # Apply static loads (Gravity + Lateral)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)

    _add_nodal_loads(
        num_floors,
        floor_nodes,
        distributed_floor_loads,
        distributed_floor_dir,
        point_load_node,
        point_load_value,
        point_load_dir,
    )
    
    def _add_dead_load_mass(num_floors, floor_nodes, floor_dead_load_kilonewton, dead_load_accel=DEAD_LOAD_G_DEFAULT):
        if floor_dead_load_kilonewton is None: return
        for i_floor, floor_mass in enumerate(abs(floor_dead_load_kilonewton) / dead_load_accel):
            nodes = floor_nodes[i_floor]
            if not nodes: continue
            node_count = len(nodes)
            weight = np.ones(node_count) / node_count
            for n, w in zip(nodes, weight):
                # mass is scalar or vector? openseespy takes args for mass(node, m1, m2, m3, r1, r2, r3)
                m = floor_mass * w
                current = ops.nodeMass(n)
                if current is None: current = [0]*6
                # we add to translational mass (X,Y,Z)
                # Note: rigid diaphragm handles coupling but typically we put mass on nodes
                # floor_mass is total floor mass. each node gets share.
                new_mass = [current[0] + m, current[1] + m, current[2] + m, 0.0, 0.0, 0.0]
                ops.mass(n, *new_mass)

    if include_dead_load_as_mass:
        _add_dead_load_mass(num_floors, floor_nodes, floor_dead_load_kilonewton)

    # -----------------------------
    # static + modal
    # -----------------------------
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-6, 40)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ok = ops.analyze(1)
    if ok != 0:
        raise RuntimeError("Static analysis did not converge.")

    _write_static_results(out_dir, point_load_node, num_floors, floor_nodes)

    # modal
    eig = np.array(ops.eigen(10))
    eig = eig[eig > 1e-12]
    omegas = np.sqrt(eig)
    freqs = omegas / (2.0 * math.pi)
    with open(out_dir / "three_story_20_static_modal.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Mode", "Frequency_Hz", "Omega_rad_s"])
        for i, (f_hz, om) in enumerate(zip(freqs, omegas), start=1):
            writer.writerow([i, f_hz, om])

    print("Natural frequencies (Hz):")
    for i, f_hz in enumerate(freqs, start=1):
        print(f"Mode {i:2d}: {f_hz:.8e}")

    print(f"\nStatic Results written to: {out_dir}")


if __name__ == "__main__":
    run_static_analysis()
