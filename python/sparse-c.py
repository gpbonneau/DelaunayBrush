import numpy as np
import open3d as o3d
import time
import pyvista as pv
import random

def find_disjoint_curves(filename):
    # Define a dictionary to store vertices
    vertices = []
    v = 1
    
    # Define a list to store all lines (disjoint and connected)
    lines = []

    # Read the OBJ file
    with open(filename, "r") as file:
        for line in file:
            words = line.strip().split()

            if words[0] == "v":
                # Parse vertices and store them in the vertices dictionary
                vertex_coords = tuple(map(float, words[1:]))
                vertices.append(np.array(vertex_coords))
                v += 1

            elif words[0] == "l":
                # Parse lines and create curves
                line_indices = list(map(int, words[1:]))
                line = [vertex_id for vertex_id in line_indices]
                lines.append(line)

    # Function to perform depth-first search (DFS) to find disjoint curves
    def dfs(start, visited, graph, curve):
        visited[start] = True
        curve.append(start)
        for neighbor in graph[start]:
            if not visited[neighbor]:
                dfs(neighbor, visited, graph, curve)

    # Create a graph representation of lines
    graph = {i: [] for i in range(len(lines))}
    for i, line1 in enumerate(lines):
        for j, line2 in enumerate(lines):
            if i != j:
                # Check if lines share any vertices
                if any(vertex in line1 for vertex in line2):
                    graph[i].append(j)

    # Perform DFS to find disjoint curves
    visited = [False] * len(lines)
    disjoint_curves = []
    disjoint_curves_display=[]
    for i in range(len(lines)):
        if not visited[i]:
            curve_indices = []
            dfs(i, visited, graph, curve_indices)
            disjoint_curves_display.append([lines[idx] for idx in curve_indices[:-1]])
            disjoint_curves.append([lines[idx][0] for idx in curve_indices])
    return vertices,lines,disjoint_curves

def visualize_mesh(mesh, curve_points, disjoint_curves_display, curve_check, output_obj_file):
    """
    Visualize the mesh and disjoint curves using Open3D and write to an OBJ file.

    Parameters:
    - mesh (o3d.geometry.TriangleMesh): The mesh to visualize.
    - curve_points (np.ndarray): Array of curve points.
    - disjoint_curves_display (list): List of disjoint curves.
    - curve_check (list): List of curve segments.
    - output_obj_file (str): Path to the output OBJ file.
    """
    vis = o3d.visualization.Visualizer()
    vis.create_window(width=800, height=600)

    # Create Open3D LineSet for mesh edges
    edges = create_edges(mesh)
    lines = [[edges[i][0], edges[i][1]] for i in range(len(edges))]
    lineset = o3d.geometry.LineSet(points=o3d.utility.Vector3dVector(mesh.vertices), lines=o3d.utility.Vector2iVector(lines))
    vis.add_geometry(lineset)

    valid_curves = []
    for curve_full in disjoint_curves_display:
        valid = True
        start = 0
        for i in range(len(curve_full) - 1):
            if [curve_full[i], curve_full[i+1]] not in curve_check and [curve_full[i+1], curve_full[i]] not in curve_check:
                if len(curve_full[start:i]) > 1:
                    valid_curves.append(curve_full[start:i])
                    start = i + 1

    # Randomly select 60% of valid curves
    valid_curves = random.sample(valid_curves, int(len(valid_curves) * 0.6))

    with open(output_obj_file, "w") as obj_file:
        # Write vertices
        

        # Write lines
        vcount = 0
        for curve_full in valid_curves:
            curve = curve_full[::skip]
            
            obj_file.write("v {} {} {}\n".format(curve_points[curve[0]][0], curve_points[curve[0]][1], curve_points[curve[0]][2]))
            vcount += 1
            for i in range(len(curve) - 1):
                obj_file.write("v {} {} {}\n".format(curve_points[curve[i + 1]][0], curve_points[curve[i + 1]][1], curve_points[curve[i + 1]][2]))
                obj_file.write("l {} {}\n".format(vcount, vcount + 1))  # Vertex indices are 1-based
                vcount += 1
    # Visualize disjoint curves
    for curve_full in valid_curves:
        curve = curve_full[::skip]
        line = [(curve[i] - 1, curve[i + 1] - 1) for i in range(len(curve) - 1)]  # Convert vertex indices to edge tuples
        line = o3d.utility.Vector2iVector(line)
        lines = o3d.geometry.LineSet(points=o3d.utility.Vector3dVector(curve_points), lines=line)
        lines.paint_uniform_color([1.0, 0.0, 0.0])
        vis.add_geometry(lines)

    # Add mesh to visualization
    vis.add_geometry(mesh)

    # Run visualization
    vis.run()

    # Destroy visualization window
    vis.destroy_window()

def create_edges(mesh):
    """
    Create edges from a mesh.

    Parameters:
    - mesh (o3d.geometry.TriangleMesh): The input mesh.

    Returns:
    - edges (list): List of edges, each represented as a tuple of vertex indices.
    """
    edges = set()

    for face in np.asarray(mesh.triangles):
        for i in range(3):
            edge = (face[i], face[(i + 1) % 3])
            edges.add(tuple(sorted(edge)))  # Sort vertices to avoid duplicate edges

    return list(edges)

import sys

def main(filename):
    # Read OBJ file
    vertices, disjoint_curves_display, all_curve = find_disjoint_curves(filename)

    # Load mesh
    mesh = o3d.io.read_triangle_mesh(filename)
    # print(disjoint_curves_display)
    # Visualize mesh and disjoint curves
    visualize_mesh(mesh, vertices, all_curve,disjoint_curves_display,"../curves.obj")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <filename>")
        sys.exit(1)
    
    filename = sys.argv[1]
    skip=int(sys.argv[2])
    main(filename)

