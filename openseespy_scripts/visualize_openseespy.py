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

def build_model_for_vis():
    # -----------------------------
    # model settings (Same as three_story_20_static.py)
    # -----------------------------
    num_floors = 21 # 21 Levels (0 to 20)
    floor_height = 4.0
    span_x = 4.0
    span_z = 4.0

    # material + section
    E = 2.0e7
    G = E / (2.0 * (1.0 + 0.167))
    J = 1.0 # arbitrary non-zero
    A = 1.0
    Iy = 1.0
    Iz = 1.0

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
        for x, z in ((0.0, 0.0), (0.0, span_z), (span_x, span_z), (span_x, 0.0)):
            ops.node(node_id, x, y, z)
            floor.append(node_id)
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
        
        # floor beams at TOP of the story (using corrected indices)
        base_top = base + 4
        n1 = base_top + 1
        n2 = base_top + 2
        n3 = base_top + 3
        n4 = base_top + 4
        
        beam_pairs = [(n1, n2, 2), (n2, n3, 1), (n3, n4, 2), (n4, n1, 1)]
        for (i, j, transf) in beam_pairs:
            ops.element("elasticBeamColumn", ele_id, i, j,
                        A, E, G, J, Iy, Iz, transf)
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
