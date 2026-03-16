"""Visualize the structure using OpenSeesPy's native plotting tools.

This script rebuilds the model (geometry only) and uses `openseespy.postprocessing.Get_Rendering.plot_model`
to generate a 3D view.
"""

import sys
from pathlib import Path

import sys
import struct
import platform

print("=== Python Environment Check ===")
print("Version      :", sys.version)
print("Executable   :", sys.executable)
print("Architecture :", struct.calcsize("P") * 8, "bit")
print("Platform     :", platform.platform())
print("=================================")
import openseespy.opensees as ops
print("OpenSeesPy import OK")

# Add script directory to sys.path to allow importing modules if needed,
# though we will just copy the model building logic for standalone execution.
script_dir = Path(__file__).resolve().parent
sys.path.append(str(script_dir))

#import openseespy.opensees as ops
import vfo.vfo as vfo
import pyvista as pv
from three_story_20_common import E_DEFAULT, FLOOR_HEIGHT, NUM_BAYS_X, NUM_BAYS_Z, NUM_STORIES, SPAN_X, SPAN_Z, level_count, section_properties

def build_model_for_vis():
    # -----------------------------
    # model settings (Same as three_story_20_static.py)
    # -----------------------------
    num_floors = level_count(NUM_STORIES)
    floor_height = FLOOR_HEIGHT
    span_x = SPAN_X
    span_z = SPAN_Z
    num_bays_x = NUM_BAYS_X
    num_bays_z = NUM_BAYS_Z

    # material + section
    E = E_DEFAULT
    A, G, J, Iy, Iz = section_properties(E=E)

    ops.wipe()
    ops.model("Basic", "-ndm", 3, "-ndf", 6)
    
    ops.geomTransf("Linear", 1, 0.0, 0.0, 1.0)
    ops.geomTransf("Linear", 2, 1.0, 0.0, 0.0)

    # -----------------------------
    # Geometry Build
    # -----------------------------
    node_id = 1
    floor_nodes = []

    for i_floor in range(num_floors):
        y = floor_height * i_floor
        floor = []
        for iz in range(num_bays_z + 1):
            for ix in range(num_bays_x + 1):
                x = ix * span_x
                z = iz * span_z
                ops.node(node_id, x, y, z)
                floor.append(node_id)
                node_id += 1
        floor_nodes.append(floor)

    ele_id = 1
    for i_floor in range(num_floors - 1):
        floor_bottom = floor_nodes[i_floor]
        floor_top = floor_nodes[i_floor + 1]

        for local_idx, n1 in enumerate(floor_bottom):
            n2 = floor_top[local_idx]
            ops.element("elasticBeamColumn", ele_id, n1, n2, A, E, G, J, Iy, Iz, 1)
            ele_id += 1

        for iz in range(num_bays_z + 1):
            row = iz * (num_bays_x + 1)
            for ix in range(num_bays_x):
                n1 = floor_top[row + ix]
                n2 = floor_top[row + ix + 1]
                ops.element("elasticBeamColumn", ele_id, n1, n2, A, E, G, J, Iy, Iz, 1)
                ele_id += 1

        for iz in range(num_bays_z):
            row = iz * (num_bays_x + 1)
            next_row = (iz + 1) * (num_bays_x + 1)
            for ix in range(num_bays_x + 1):
                n1 = floor_top[row + ix]
                n2 = floor_top[next_row + ix]
                ops.element("elasticBeamColumn", ele_id, n1, n2, A, E, G, J, Iy, Iz, 2)
                ele_id += 1
            
    # base support
    for n in floor_nodes[0]:
        ops.fix(n, 1, 1, 1, 1, 1, 1)

def run_visualization():
    build_model_for_vis()
    print("Model built. Generating vfo visualization...")
    
    # Configure pyvista for off-screen rendering to save image
    pv.OFF_SCREEN = True
    
    try:
        # vfo.plot_model(model_name="3D_Frame")
        output_image = "structure_view_vfo"
        vfo.plot_model(filename=output_image)
        print(f"Visualization saved to {output_image}.png")
    except Exception as e:
        print(f"Error during visualization: {e}")

if __name__ == "__main__":
    run_visualization()
