import numpy as np
import pyvista as pv
import random
import argparse
import os
import networkx as nx



def compute_circumsphere(vertices):
    """
    Compute the circumcenter and circumradius of a tetrahedron.
    
    Parameters:
        vertices (np.ndarray): A 4x3 array where each row represents a vertex (x, y, z).
    
    Returns:
        tuple: (circumcenter (np.ndarray), circumradius (float))
    """
    A, B, C, D = vertices  # Extract vertices
    
    # Compute squared norms of each vertex
    a2, b2, c2, d2 = [np.dot(v, v) for v in [A, B, C, D]]

    # Form determinant matrices
    M = np.array([
        [A[0], A[1], A[2], 1],
        [B[0], B[1], B[2], 1],
        [C[0], C[1], C[2], 1],
        [D[0], D[1], D[2], 1]
    ])

    Dx = np.array([
        [a2, A[1], A[2], 1],
        [b2, B[1], B[2], 1],
        [c2, C[1], C[2], 1],
        [d2, D[1], D[2], 1]
    ])

    Dy = np.array([
        [A[0], a2, A[2], 1],
        [B[0], b2, B[2], 1],
        [C[0], c2, C[2], 1],
        [D[0], d2, D[2], 1]
    ])

    Dz = np.array([
        [A[0], A[1], a2, 1],
        [B[0], B[1], b2, 1],
        [C[0], C[1], c2, 1],
        [D[0], D[1], d2, 1]
    ])

    D = np.array([
        [A[0], A[1], A[2], a2],
        [B[0], B[1], B[2], b2],
        [C[0], C[1], C[2], c2],
        [D[0], D[1], D[2], d2]
    ])

    if abs(np.linalg.det(M)) < 1e-10:
        return None, None  # Degenerate tetrahedron

    M_pinv = np.linalg.pinv(M)  # Use pseudo-inverse to handle precision issues

    Cx = np.linalg.det(Dx) @ M_pinv[:, 3] / 2
    Cy = -np.linalg.det(Dy) @ M_pinv[:, 3] / 2
    Cz = np.linalg.det(Dz) @ M_pinv[:, 3] / 2

    circumcenter = np.array([Cx, Cy, Cz])
    circumradius = np.linalg.norm(circumcenter - A)

    return circumcenter, circumradius


def compute_circumsphere(vertices):
    """
    Compute the circumcenter and circumradius of a tetrahedron.
    
    Parameters:
        vertices (np.ndarray): A 4x3 array where each row represents a vertex (x, y, z).
    
    Returns:
        tuple: (circumcenter (np.ndarray), circumradius (float))
    """
    # Ensure vertices is a numpy array
    vertices = np.array(vertices)
    
    # Extract the vertices
    a, b, c, d = vertices
    
    # Compute the differences
    ab = b - a
    ac = c - a
    ad = d - a
    
    # Compute the determinant
    det = np.linalg.det(np.array([ab, ac, ad]))
    
    # Compute the circumcenter
    circumcenter = a + np.cross(np.cross(ab, ac), ad) * np.dot(ad, ad) + \
                   np.cross(np.cross(ac, ad), ab) * np.dot(ab, ab) + \
                   np.cross(np.cross(ad, ab), ac) * np.dot(ac, ac)
    circumcenter /= (2 * det**2)
    
    # Compute the circumradius
    circumradius = np.linalg.norm(a - circumcenter)
    
    return circumcenter, circumradius


def find_connected_components(edges):

    """
    Find connected components (disjoint curves) from the edges using a graph.

    Parameters:
    - edges (list of tuples): List of edges represented as tuples of vertex indices.

    Returns:
    - components (list of lists): List of connected components (each component is a list of edges).
    """
    # Create an undirected graph using NetworkX
    G = nx.Graph()
    
    # Add edges to the graph
    G.add_edges_from(edges)

    # Find connected components (disjoint curves)
    components = [list(component) for component in nx.connected_components(G)]
    
    return components

def extract_vertices_with_adjacency(mesh, adjacency_count=4, filter_percentage=0.5):
    """
    Extracts edges formed by vertices with a specific number of adjacent faces, applies random filtering,
    and updates the edges and normals accordingly.
    """
    faces = mesh.faces.reshape((-1, 5))[:, 1:]  # Extract faces without face type info
    
    # Count the number of faces each vertex is part of
    vertex_face_count = np.zeros(mesh.n_points, dtype=int)
    for face in faces:
        vertex_face_count[face] += 1
    
    # Find vertices with exactly the specified adjacency
    vertices_with_adjacency_4 = np.where(vertex_face_count == adjacency_count)[0]
    
    mesh.compute_normals(cell_normals=False, point_normals=True)
    original_normals = mesh['Normals']
    new_normals = original_normals[vertices_with_adjacency_4]
    
    edges = set()
    for face in faces:
        adjacent_vertices = [v for v in face if v in vertices_with_adjacency_4]
        if len(adjacent_vertices) >= 2:
            for i in range(len(adjacent_vertices)):
                for j in range(i + 1, len(adjacent_vertices)):
                    edge = tuple(sorted([adjacent_vertices[i], adjacent_vertices[j]]))
                    edges.add(edge)
    
    lines = find_connected_components(edges)
    num_selected = int(len(lines) * filter_percentage)  # Randomly select a percentage of lines
    filtered_lines = random.sample(lines, num_selected)
    
    newedges = []
    selected_vertices = set()
    for line in filtered_lines:
        for i in range(len(line)):
            if tuple(sorted([line[i], line[(i+1) % len(line)]])) in edges:
                newedges.append(tuple(sorted([line[i], line[(i+1) % len(line)]])))
                selected_vertices.update([line[i], line[(i+1) % len(line)]])
    
    # Update vertex indices and normals
    final_vertices_list = sorted(selected_vertices)
    vertex_map = {v: i for i, v in enumerate(final_vertices_list)}
    final_edges = [(vertex_map[edge[0]], vertex_map[edge[1]]) for edge in newedges]
    final_normals = np.array([new_normals[np.where(vertices_with_adjacency_4 == v)[0][0]] for v in final_vertices_list])
    
    # Convert edges to PyVista format
    final_edges = np.array(final_edges)
    final_edges = np.column_stack([np.full(final_edges.shape[0], 2), final_edges])
    
    wireframe_mesh = pv.PolyData(mesh.points[final_vertices_list])
    wireframe_mesh.lines = final_edges.flatten()
    wireframe_mesh['Normals'] = final_normals
    
    return wireframe_mesh


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract wireframe mesh with adjacency filtering.")
    parser.add_argument("input_file", type=str, help="Path to the input mesh file.")
    parser.add_argument("filter_percentage", type=float, help="Percentage of edges to retain after filtering.")
    args = parser.parse_args()

    mesh = pv.read(args.input_file)
    wireframe_mesh = extract_vertices_with_adjacency(mesh, adjacency_count=4,filter_percentage=args.filter_percentage)
    
    # Create output directory named 'curves' in the same folder as the input file
    input_dir = os.path.dirname(args.input_file)
    output_dir = os.path.join(input_dir, "curves", f"{int(args.filter_percentage * 100)}")
    os.makedirs(output_dir, exist_ok=True)
    
    # Construct output file paths
    base_name, ext = os.path.splitext(os.path.basename(args.input_file))
    output_file = os.path.join(output_dir, f"{base_name}_filtered{ext}")
    normal_output_file = os.path.join(output_dir, f"{base_name}_filtered_normals{ext}")
    obj_file=os.path.join(output_dir, f"{base_name}_filtered.obj")
    
    # Save outputs
    wireframe_mesh.save(output_file)
    wireframe_mesh.save(obj_file)
    normal_mesh = pv.PolyData(wireframe_mesh['Normals'])
    normal_mesh.save(normal_output_file)
    
    print(f"Wireframe saved in {output_dir} as {output_file} and {normal_output_file}")