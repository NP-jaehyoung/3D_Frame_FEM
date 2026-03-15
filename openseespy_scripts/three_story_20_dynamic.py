"""OpenSeesPy Dynamic Analysis for 20-story 3D rigid-diaphragm frame.

- 20 floors, 4 nodes per floor (4m x 4m square plan, 4m story height)
- 6 DOF node (UX,UY,UZ,RX,RY,RZ)
- Rigid diaphragm enforced per floor
- Loads:
    - Gravity: Distributed floor load (5 kN/m^2) -> Converted to Mass
    - Dynamic: Uniform Ground Exception (Sine Pulse)
- Lateral Point Load is EXCLUDED.
- Output:
    - Top Node Displacement Time History
    - Inter-story Drift Time History
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

from three_story_20_common import (
    AG_CENTER_DEFAULT,
    AG_FREQ_DEFAULT,
    AG_WIDTH_DEFAULT,
    AG_PEAK_G_DEFAULT,
    DEAD_LOAD_G_DEFAULT,
    DAMPING_RATIO_DEFAULT,
    DT_DEFAULT,
    FLOOR_HEIGHT,
    FLOOR_LOAD_DIR,
    FLOOR_LOAD_PRESSURE_KPA,
    INCLUDE_BASE_FLOOR_LOAD,
    INCLUDE_DEAD_LOAD_AS_MASS,
    NUM_STORIES,
    SPAN_X,
    SPAN_Z,
    TIME_END_DEFAULT,
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


def _add_gravity_loads(num_floors: int, floor_nodes, distributed_floor_loads, distributed_floor_dir):
    component_map = {"UX": 1, "UY": 2, "UZ": 3, "RX": 4, "RY": 5, "RZ": 6}
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


def _run_dynamic(n_steps: int, dt: float, floor_nodes, top_node: int, out_dir: Path):
    # record responses during dynamic
    time = np.zeros(n_steps)
    top_ux = np.zeros(n_steps)
    top_uy = np.zeros(n_steps)
    top_uz = np.zeros(n_steps)
    floor_count = len(floor_nodes)

    # helper to sample inter-story drift in UX for each floor
    levels = floor_nodes
    drift = np.zeros((n_steps, floor_count - 1))

    print("Running dynamic analysis...")
    for i in range(n_steps):
        time[i] = ops.getTime()
        top_ux[i] = ops.nodeDisp(top_node, 1)
        top_uy[i] = ops.nodeDisp(top_node, 2)
        top_uz[i] = ops.nodeDisp(top_node, 3)

        for i_floor in range(1, floor_count):
            top_nodes = levels[i_floor]
            bot_nodes = levels[i_floor - 1]
            mean_ux_top = np.mean([ops.nodeDisp(n, 1) for n in top_nodes])
            mean_ux_bot = np.mean([ops.nodeDisp(n, 1) for n in bot_nodes])
            drift[i, i_floor - 1] = mean_ux_top - mean_ux_bot

        if i < n_steps - 1:
            ok = ops.analyze(1, dt)
            if ok != 0:
                raise RuntimeError(f"Transient analysis failed at step {i}, t={ops.getTime():.3f}")
                
    print("Dynamic analysis completed.")

    with open(out_dir / "three_story_20_dynamic_top.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Time_s", "TopUx_m", "TopUy_m", "TopUz_m"])
        for t, ux, uy, uz in zip(time, top_ux, top_uy, top_uz):
            writer.writerow([t, ux, uy, uz])

    with open(out_dir / "three_story_20_dynamic_story_drift.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Time_s"] + [f"Story{i+1}" for i in range(floor_count - 1)])
        for i in range(n_steps):
            writer.writerow([time[i], *drift[i].tolist()])


def run_dynamic_analysis():
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
    
    # 5 kN/m^2 distributed load (total load per floor = pressure * area)
    # Applied in negative Y direction (gravity)
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

    # dynamic setup
    dt = DT_DEFAULT
    time_end = TIME_END_DEFAULT
    ag_peak = AG_PEAK_G_DEFAULT * DEAD_LOAD_G_DEFAULT
    ag_freq = AG_FREQ_DEFAULT
    ag_center = AG_CENTER_DEFAULT
    ag_width = AG_WIDTH_DEFAULT
    time_steps = np.arange(0.0, time_end + 0.5 * dt, dt)
    ground_acc = ag_peak * np.exp(-((time_steps - ag_center) / ag_width) ** 2) * np.sin(2 * np.pi * ag_freq * time_steps)
    n_steps = len(time_steps)
    damping_ratio = DAMPING_RATIO_DEFAULT

    # -----------------------------
    # model build
    # -----------------------------
    script_dir = Path(__file__).resolve().parent
    out_dir = script_dir / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    ops.wipe()
    ops.model("Basic", "-ndm", 3, "-ndf", 6)
    
    # geom transformation for 3D linear elements
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.geomTransf("Linear", 2, 1.0, 0.0, 0.0)

    nodes, floor_nodes = _build_geometry(num_floors, floor_height, span_x, span_z)
    
    # base support
    for n in floor_nodes[0]:
        ops.fix(n, 1, 1, 1, 1, 1, 1)

    _apply_rigid_diaphragms(num_floors, floor_nodes)

    # 1. Apply Gravity Loads (Pattern 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    
    _add_gravity_loads(
        num_floors,
        floor_nodes,
        distributed_floor_loads,
        distributed_floor_dir
    )
    
    # Define Mass
    if include_dead_load_as_mass:
        if floor_dead_load_kilonewton is not None:
             for i_floor, floor_mass in enumerate(abs(floor_dead_load_kilonewton) / DEAD_LOAD_G_DEFAULT):
                nodes = floor_nodes[i_floor]
                if not nodes: continue
                node_count = len(nodes)
                weight = np.ones(node_count) / node_count
                for n, w in zip(nodes, weight):
                    m = floor_mass * w
                    current = ops.nodeMass(n)
                    if current is None: current = [0]*6
                    new_mass = [current[0] + m, current[1] + m, current[2] + m, 0.0, 0.0, 0.0]
                    ops.mass(n, *new_mass)

    # -----------------------------
    # Gravity Analysis & Fix
    # -----------------------------
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-6, 40)
    ops.algorithm("Newton")
    ops.integrator("LoadControl", 1.0)
    ops.analysis("Static")
    ops.analyze(1)
    
    # Fix gravity loads
    ops.loadConst("-time", 0.0)
    
    # -----------------------------
    # Eigen (needed for Damping)
    # -----------------------------
    eig = np.array(ops.eigen(10))
    eig = eig[eig > 1e-12]
    freqs = np.sqrt(eig) / (2.0 * math.pi)
    
    if len(freqs) >= 2:
        w1, w2 = freqs[0] * 2.0 * math.pi, freqs[1] * 2.0 * math.pi
        alpha_m = 2.0 * damping_ratio * w1 * w2 / (w1 + w2)
        beta_k = 2.0 * damping_ratio / (w1 + w2)
    elif len(freqs) >= 1:
        w1 = freqs[0] * 2.0 * math.pi
        alpha_m = 0.0
        beta_k = 2.0 * damping_ratio / w1
    else:
        alpha_m = 0.0
        beta_k = 0.0
        
    ops.rayleigh(alpha_m, beta_k, 0.0, 0.0)

    # -----------------------------
    # Dynamic Analysis
    # -----------------------------
    # Pattern 2: Ground Motion
    ops.timeSeries("Path", 2, "-dt", dt, "-values", *ground_acc.tolist(), "-factor", 1.0)
    ops.pattern("UniformExcitation", 2, 1, "-accel", 2) # Direction 1 (UX)
    
    ops.wipeAnalysis()
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-6, 40)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    
    _run_dynamic(n_steps, dt, floor_nodes, point_load_node, out_dir)

    print(f"\nDynamic Results written to: {out_dir}")


if __name__ == "__main__":
    run_dynamic_analysis()
