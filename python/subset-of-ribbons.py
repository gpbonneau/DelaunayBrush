import numpy as np
import pyvista as pv

def extract_vertices_with_adjacency(mesh, adjacency_count=4):
    """
    Extracts edges formed by vertices with a specific number of adjacent faces.

    Parameters:
    mesh (pyvista.PolyData): Input quad mesh.
    adjacency_count (int): The number of adjacent faces a vertex must have (default is 4 for quad meshes).

    Returns:
    pyvista.PolyData: Wireframe mesh formed by the vertices with the specified adjacency and their edges.
    """
    # Get the faces of the mesh (these are the indices of the vertices forming each face)
    faces = mesh.faces.reshape((-1, 5))[:, 1:]  # Convert faces to a list of vertex indices, skip the face type info

    # Count the number of faces each vertex is part of
    vertex_face_count = np.zeros(mesh.n_points, dtype=int)
    for face in faces:
        vertex_face_count[face] += 1
    
    # Find vertices with exactly 4 adjacent faces
    vertices_with_adjacency_4 = np.where(vertex_face_count == adjacency_count)[0]

    mesh.compute_normals(cell_normals=False, point_normals=True)
    original_normals = mesh['Normals']
    new_normals = original_normals[vertices_with_adjacency_4]

    # Collect edges formed by these vertices
    edges = set()  # Use a set to avoid duplicate edges

    for face in faces:
        # Find vertices in the face that have 4 adjacencies
        adjacent_vertices = [v for v in face if v in vertices_with_adjacency_4]

        # If there are at least two vertices with 4 adjacencies, form edges between them
        if len(adjacent_vertices) >= 2:
            print(adjacent_vertices)
            for i in range(len(adjacent_vertices)):
                for j in range(i + 1, len(adjacent_vertices)):
                    # Get the indices of the vertices in the 'vertices_with_adjacency_4' list
                    idx_i = np.where(vertices_with_adjacency_4 == adjacent_vertices[i])[0][0]
                    idx_j = np.where(vertices_with_adjacency_4 == adjacent_vertices[j])[0][0]

                    # Sort to ensure consistency in edge representation
                    edge = tuple(sorted([idx_i, idx_j]))
                    edges.add(edge)
    print(edges)
    # Convert edges to the format PyVista expects: [number of points, point1, point2, ...]
    edges = np.array(list(edges))  # Convert set to array
    edges = np.column_stack([np.full(edges.shape[0], 2), edges])  # Prepend 2 to each edge (2 points per edge)

    # Create a new mesh using the vertices that are part of the edges (vertices_with_adjacency_4)
    wireframe_mesh = pv.PolyData(mesh.points[vertices_with_adjacency_4])  # Only include vertices with adjacency 4
    wireframe_mesh.lines = edges.flatten(   )  # Flatten the edge list into the correct format
    wireframe_mesh['Normals'] = np.array(new_normals)
    return wireframe_mesh

# Load the quad mesh directly from the file using PyVista
mesh_file_path = "/Users/anandhu/Documents/baloon/Inputs_and_Results/Stanford_Bunny/01_Stanford_Bunny_Ribbons.ply"  # Change to your actual file path
quad_mesh = pv.read(mesh_file_path)

# Extract the wireframe formed by vertices with adjacency 4 and their connective edges
wireframe_mesh = extract_vertices_with_adjacency(quad_mesh)

# Visualize the wireframe (with edges connecting vertices of adjacency 4)
wireframe_mesh.plot(show_edges=True, line_width=2, color='orange')

# Save the wireframe as an OBJ or PLY file
wireframe_mesh.save("skeletal_wireframe.obj")
wireframe_mesh.save("skeletal_wireframe.ply")

print("Wireframe saved as skeletal_wireframe.obj and skeletal_wireframe.ply")
