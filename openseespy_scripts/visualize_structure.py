"""Visualize the 20-story multi-bay 3D frame structure.

This script generates a 3D plot of the structure:
- Frame (Columns/Beams): Black lines
- Slabs (Floors): Green semi-transparent surfaces
- Supports: Marked at the base
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from three_story_20_common import FLOOR_HEIGHT, NUM_BAYS_X, NUM_BAYS_Z, SPAN_X, SPAN_Z

def visualize_structure(
    num_floors=21,
    floor_height=FLOOR_HEIGHT,
    span_x=SPAN_X,
    span_z=SPAN_Z,
    num_bays_x=NUM_BAYS_X,
    num_bays_z=NUM_BAYS_Z,
):
    fig = plt.figure(figsize=(10, 15))
    ax = fig.add_subplot(111, projection='3d')
    
    # Structure dimensions
    # 21 levels (0..20) -> 20 stories. Height is 20*4 = 80.
    total_height = (num_floors - 1) * floor_height
    
    # Collections for plotting
    lines = []
    
    # 1. Draw Columns (Vertical Lines)
    width_x = num_bays_x * span_x
    width_z = num_bays_z * span_z

    for iz in range(num_bays_z + 1):
        for ix in range(num_bays_x + 1):
            x = ix * span_x
            z = iz * span_z
            ax.plot([x, x], [z, z], [0, total_height], color='black', linewidth=1.0)
        
    # 2. Draw Floors (Beams + Slabs)
    for i in range(num_floors):
        if i == 0: continue
        y = floor_height * i
        
        # Floor polygon vertices (x, z, y_vertical)
        # Verts in structure coords: (x, z, y_level)
        verts = [
            (0, 0, y),
            (0, width_z, y),
            (width_x, width_z, y),
            (width_x, 0, y)
        ]
        
        # Draw beams (Black lines)
        # (0,0)->(0,sz)->(sx,sz)->(sx,0)->(0,0)
        # Close the loop
        xs = [v[0] for v in verts] + [verts[0][0]]
        zs = [v[1] for v in verts] + [verts[0][1]] # Struct Z -> Plot Y
        ys = [v[2] for v in verts] + [verts[0][2]] # Struct Y -> Plot Z
        
        ax.plot(xs, zs, ys, color='black', linewidth=1.0)

        for iz in range(num_bays_z + 1):
            z = iz * span_z
            ax.plot([0, width_x], [z, z], [y, y], color='black', linewidth=0.8)
        for ix in range(num_bays_x + 1):
            x = ix * span_x
            ax.plot([x, x], [0, width_z], [y, y], color='black', linewidth=0.8)
        
        # Draw Slab (Green Surface)
        # Poly3DCollection expects list of (x, y, z) tuples for vertices
        # We want Plot(x, z, y)
        poly_plot = [[(v[0], v[1], v[2]) for v in verts]] 
        
        slab = Poly3DCollection(poly_plot, facecolors='green', linewidths=0, alpha=0.3)
        ax.add_collection3d(slab)

    # 3. Base Supports
    # Plot points at y=0 (Plot Z=0)
    # Scatter(x, y, z) -> (Struct X, Struct Z, 0)
    support_x = []
    support_z = []
    for iz in range(num_bays_z + 1):
        for ix in range(num_bays_x + 1):
            support_x.append(ix * span_x)
            support_z.append(iz * span_z)
    ax.scatter(support_x, support_z, [0] * len(support_x), color='red', s=30, marker='^', label='Fixed Support')

    # 4. Top Point Load (Arrow)
    # Load is applied at (span_x, 0, total_height) in Struct coords.
    # Load direction is -X (negative UX).
    # Arrow Head: (span_x, 0, total_height) -> Plot(span_x, 0, total_height)
    # Arrow Tail: (span_x + 5, 0, total_height) -> Plot(span_x + 5, 0, total_height)
    # Direction vector: (-5, 0, 0) -> Plot(-5, 0, 0)
    
    # Text label position
    ax.text(width_x + 6, 0, total_height, 'Lateral Load (-UX)', color='blue')
    
    # Quiver(x, y, z, u, v, w)
    # Position: (span_x + 5, 0, total_height) (Tail)
    # Vector: (-5, 0, 0)
    ax.quiver(width_x + 5, 0, total_height, -5, 0, 0, color='blue', length=5, arrow_length_ratio=0.3)

    # Settings
    ax.set_xlabel('X (Span X)')
    ax.set_ylabel('Y (Span Z)')
    ax.set_zlabel('Z (Height Y)')
    ax.set_title(f'3D Model of {num_floors}-Story Frame')
    
    # Set limits
    max_range = np.array([width_x + 5, width_z, total_height]).max()
    mid_x = (width_x + 5) * 0.5
    mid_y = width_z * 0.5
    mid_z = total_height * 0.5
    
    ax.set_xlim(mid_x - max_range*0.5, mid_x + max_range*0.5)
    ax.set_ylim(mid_y - max_range*0.5, mid_y + max_range*0.5)
    ax.set_zlim(mid_z - max_range*0.5, mid_z + max_range*0.5)
        
    # View angle
    ax.view_init(elev=20, azim=-60)
    
    plt.legend()
    plt.tight_layout()
    
    output_file = 'structure_view_3d.png'
    plt.savefig(output_file, dpi=300)
    print(f"Saved visualization to {output_file}")

if __name__ == "__main__":
    visualize_structure()
