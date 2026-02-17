"""OpenSeesPy example for 20-story 3D rigid-diaphragm frame.

Comparable setup to ThreeDimFrame_20story.m:
- 20 floors, 4 nodes per floor (4m x 4m square plan, 4m story height)
- 6 DOF node (UX,UY,UZ,RX,RY,RZ)
- Rigid diaphragm enforced per floor
- Separate control of:
    - point load (single node)
    - distributed floor load (per-floor vector, distributed to floor nodes)
- static response
- natural frequencies
- dynamic response by prescribed ground acceleration
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np

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
        # floor beams (closed square)
        n1 = base + 1
        n2 = base + 2
        n3 = base + 3
        n4 = base + 4
        for (i, j) in ((n1, n2), (n2, n3), (n3, n4), (n4, n1)):
            ops.element("elasticBeamColumn", ele_id, i, j,
                        A, E, G, J, Iy, Iz, 1)
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


def _add_dead_load_mass(num_floors: int, floor_nodes, floor_dead_load_kilonewton, dead_load_distribution="AREA", dead_load_accel=9.81):
    if floor_dead_load_kilonewton is None:
        return
    if np.isscalar(floor_dead_load_kilonewton):
        floor_dead_load_kilonewton = np.full(num_floors, float(floor_dead_load_kilonewton))
    else:
        floor_dead_load_kilonewton = np.asarray(floor_dead_load_kilonewton, dtype=float)
    if len(floor_dead_load_kilonewton) != num_floors:
        raise ValueError("floorDeadLoadKiloNewton must be scalar or length numFloors.")

    for i_floor, floor_mass in enumerate(abs(floor_dead_load_kilonewton) / dead_load_accel):
        nodes = floor_nodes[i_floor]
        if not nodes:
            continue
        node_count = len(nodes)
        if node_count == 0:
            continue
        if dead_load_distribution.upper() == "UNIFORM":
            weight = np.ones(node_count) / node_count
        elif dead_load_distribution.upper() == "AREA":
            # fallback: same as uniform for rectangle geometry
            weight = np.ones(node_count) / node_count
        else:
            raise ValueError("dead_load_distribution must be UNIFORM or AREA.")

        for n, w in zip(nodes, weight):
            m = floor_mass * w
            current = ops.nodeMass(n)
            if current is None:
                current = [0, 0, 0, 0, 0, 0]
            new_mass = [current[0] + m, current[1] + m, current[2] + m, 0.0, 0.0, 0.0]
            ops.mass(n, *new_mass)


def _dof_to_name(dof: int) -> str:
    return ["UX", "UY", "UZ", "RX", "RY", "RZ"][dof - 1]


def _write_static_results(out_dir: Path, top_node: int, n_floors: int, floor_nodes):
    # nodal displacement table
    num_nodes = 4 * n_floors
    with open(out_dir / "three_story_20_opensees_static_displacements.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Node", "DOF", "Component", "Displacement"])
        for n in range(1, num_nodes + 1):
            for dof in range(1, 7):
                writer.writerow([n, dof, _dof_to_name(dof), ops.nodeDisp(n, dof)])

    # support reactions
    with open(out_dir / "three_story_20_opensees_reactions.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["SupportNode", "DOF", "Component", "Reaction"])
        for n in range(1, 5):
            for dof in range(1, 7):
                writer.writerow([n, dof, _dof_to_name(dof), ops.nodeReaction(n, dof)])

    # top-node displacement summary
    top_disp = [ops.nodeDisp(top_node, 1), ops.nodeDisp(top_node, 2), ops.nodeDisp(top_node, 3)]
    print(f"\nTop node {top_node} displacement: UX={top_disp[0]:.6e}, UY={top_disp[1]:.6e}, UZ={top_disp[2]:.6e}")

    # applied loads summary (non-zero)
    with open(out_dir / "three_story_20_opensees_nonzero_loads.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Node", "DOF", "Component", "Load"])
        for n in range(1, num_nodes + 1):
            for dof in range(1, 7):
                val = ops.getLoad(n, dof)
                if abs(val) > 1e-12:
                    writer.writerow([n, dof, _dof_to_name(dof), val])


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

    with open(out_dir / "three_story_20_opensees_dynamic_top.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Time_s", "TopUx_m", "TopUy_m", "TopUz_m"])
        for t, ux, uy, uz in zip(time, top_ux, top_uy, top_uz):
            writer.writerow([t, ux, uy, uz])

    with open(out_dir / "three_story_20_opensees_dynamic_story_drift.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Time_s"] + [f"Story{i+1}" for i in range(floor_count - 1)])
        for i in range(n_steps):
            writer.writerow([time[i], *drift[i].tolist()])


def run_example():
    # -----------------------------
    # model settings
    # -----------------------------
    num_floors = 20
    floor_height = 4.0
    span_x = 4.0
    span_z = 4.0

    # material + section
    global E, G, J, A, Iy, Iz
    E = 2.0e7
    rho = 7850
    width = 1.0
    depth = 1.0
    A = width * depth
    nu = 0.167
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

    # loads
    point_load_node = 4 * num_floors
    point_load_value = -3000.0
    point_load_dir = "UX"
    distributed_floor_loads = np.zeros(num_floors)
    distributed_floor_dir = "UY"

    include_dead_load_as_mass = False
    floor_dead_load_kilonewton = np.zeros(num_floors)

    # dynamic setup
    dt = 0.01
    time_end = 20.0
    ag_peak = 0.30 * 9.81
    ag_freq = 1.8
    ag_center = 5.0
    ag_width = 1.1
    time_steps = np.arange(0.0, time_end + 0.5 * dt, dt)
    ground_acc = ag_peak * np.exp(-((time_steps - ag_center) / ag_width) ** 2) * np.sin(2 * np.pi * ag_freq * time_steps)
    n_steps = len(time_steps)
    damping_ratio = 0.02

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

    nodes, floor_nodes = _build_geometry(num_floors, floor_height, span_x, span_z)
    if len(nodes) != 4 * num_floors:
        raise RuntimeError("Node count mismatch.")

    # base support
    for n in floor_nodes[0]:
        ops.fix(n, 1, 1, 1, 1, 1, 1)

    _apply_rigid_diaphragms(num_floors, floor_nodes)

    _add_nodal_loads(
        num_floors,
        floor_nodes,
        distributed_floor_loads,
        distributed_floor_dir,
        point_load_node,
        point_load_value,
        point_load_dir,
    )

    if include_dead_load_as_mass:
        _add_dead_load_mass(
            num_floors,
            floor_nodes,
            floor_dead_load_kilonewton,
            dead_load_distribution="AREA",
            dead_load_accel=9.81,
        )

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
    with open(out_dir / "three_story_20_opensees_modal.csv", "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Mode", "Frequency_Hz", "Omega_rad_s"])
        for i, (f_hz, om) in enumerate(zip(freqs, omegas), start=1):
            writer.writerow([i, f_hz, om])

    print("Natural frequencies (Hz):")
    for i, f_hz in enumerate(freqs, start=1):
        print(f"Mode {i:2d}: {f_hz:.8e}")

    if len(freqs) >= 2:
        w1, w2 = freqs[0] * 2.0 * math.pi, freqs[1] * 2.0 * math.pi
        alpha_m = 2.0 * damping_ratio * w1 * w2 / (w1 + w2)
        beta_k = 2.0 * damping_ratio / (w1 + w2)
    elif len(freqs) == 1:
        w1 = freqs[0] * 2.0 * math.pi
        alpha_m = 0.0
        beta_k = 2.0 * damping_ratio / w1
    else:
        alpha_m = 0.0
        beta_k = 0.0

    ops.rayleigh(alpha_m, beta_k, 0.0, 0.0)

    # -----------------------------
    # dynamic (ground acceleration only)
    # -----------------------------
    ops.timeSeries("Path", 1, "-dt", dt, "-values", *ground_acc.tolist(), "-factor", 1.0)
    ops.pattern("UniformExcitation", 1, 1, 1)

    # keep existing static loads and apply inertial excitation
    ops.loadConst("-time", 0.0)

    ops.wipeAnalysis()
    ops.system("BandGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.test("NormDispIncr", 1e-6, 40)
    ops.algorithm("Newton")
    ops.integrator("Newmark", 0.5, 0.25)
    ops.analysis("Transient")
    ops.setTime(0.0)

    _run_dynamic(n_steps, dt, floor_nodes, point_load_node, out_dir)

    print(f"\nResults written to: {out_dir}")


if __name__ == "__main__":
    run_example()
