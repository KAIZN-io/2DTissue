import pandas as pd
import numpy as np
from numpy.linalg import norm
import meshio
import pyvista as pv
from pathlib import Path
from math import cos, sin, pi
from functools import reduce
import imageio


def read_data(folder, file_name):
    file_path = folder / file_name
    df = pd.read_csv(file_path, header=None)
    return df.to_numpy()


def get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping):
    num_r = r.shape[0]
    vertice_3D_id = np.empty(num_r, dtype=int)

    for i in range(num_r):
        distances_to_h = np.linalg.norm(halfedges_uv - r[i, :], axis=1)
        halfedges_id = np.argmin(distances_to_h)
        vertice_3D_id[i] = halfedge_vertices_mapping[halfedges_id, 0]

    return vertice_3D_id


folder = Path("data")
num_part = 19

mesh_loaded = meshio.read("meshes/ellipsoid_x4.off")
mesh_loaded_uv = meshio.read("meshes/Ellipsoid_uv.off")

halfedges_uv = pd.read_csv("halfedge_uv.csv", header=None).to_numpy()
halfedge_vertices_mapping = pd.read_csv("h_v_mapping_vector.csv", header=None).to_numpy()

plotter_3D = pv.Plotter(window_size=[1050, 1200], off_screen=True)
plotter_2D = pv.Plotter(window_size=[1050, 1200], off_screen=True)

plotter_3D.add_mesh(mesh_loaded, color="white")
plotter_3D.add_mesh(mesh_loaded, style='wireframe', color="grey", opacity=0.1)

plotter_2D.add_mesh(mesh_loaded_uv, color="white")
plotter_2D.add_mesh(mesh_loaded_uv, style='wireframe', color="grey", opacity=0.1)

r_points = pv.PolyData(np.hstack((np.full((num_part, 2), np.nan), np.zeros((num_part, 1)))))
plotter_2D.add_points(r_points, color="blue", point_size=10, style="points")

r_3D_points = pv.PolyData(np.full((num_part, 3), np.nan))
# add the points to the plotter as circles
plotter_3D.add_points(r_3D_points, color="red", point_size=10, style="points")

# Set the camera position and focal point for the 2D mesh
plotter_2D.camera_position = [(0, 0, 10), (0, 0, 0), (0, 1, 0)]

output_folder = Path("assets")
output_folder.mkdir(exist_ok=True)

frames = []

for tt in range(1, 101):
    r = read_data(folder, f"r_data_{tt}.csv")

    vertice_3D_id = get_vertice_id(r, halfedges_uv, halfedge_vertices_mapping)
    r_3D_points.points = mesh_loaded.points[vertice_3D_id]
    r_points.points = r

    plotter_3D.show(auto_close=False)
    plotter_2D.show(auto_close=False)
    frame_3D = plotter_3D.screenshot(return_img=True)
    frame_2D = plotter_2D.screenshot(return_img=True)

    frame = np.hstack((frame_3D, frame_2D))
    frames.append(frame)

plotter_3D.close()
plotter_2D.close()

video_filename = output_folder / "confined_active_particles.mp4"
imageio.mimwrite(str(video_filename), frames, fps=6)
