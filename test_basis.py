import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from Graph import Graph, addLimitPoints
import time

points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
                   [2, 1], [2, 2]])
points = addLimitPoints(points)

vor = Voronoi(points)

g = Graph(vor)

# # check region-relative colouring
# for r in g.regions:
#     g.plot_graph()
#     g.colour_region(r)
#     g.colour_region_edges(r, sorted=True, clockwise=False)

# check edge-relative colouring
# for edge in g.edges:
#     g.plot_graph()
#     g.colour_edge(edge)
#     regions = g.get_regions_by_edge(edge)
#     for r in regions:
#         g.colour_region(r)

# boundary checking
# boundary_vertices = g.get_boundary_vertices(clockwise=False)
# for i in range(len(boundary_vertices)):
#     edge = [boundary_vertices[i], boundary_vertices[(i+1)%len(boundary_vertices)]]
#     g.colour_edge(edge)

# # check that regions have been sorted in the graph constructor
# # plotting consecutive vertices in the region list
# for r in range(len(g.regions)):
#     region = g.regions[r]
#     for i in range(len(region)):
#         edge = [region[i], region[(i+1)%len(region)]]
#         g.colour_edge(edge)

# # plot boundary regions only
# boundary_vertices = g.get_boundary_vertices(clockwise=False)
# for i in range(len(boundary_vertices)):
#     v0 = boundary_vertices[i]
#     v1 = boundary_vertices[(i+1)%len(boundary_vertices)]
#     edge0 = [v0, v1] if v0 < v1 else [v1, v0]
#     region = g.get_regions_by_edge(edge0)[0]
#     g.colour_region(region)


boundary_vertices = g.get_boundary_vertices(clockwise=False)
v0 = boundary_vertices[0]
v1 = boundary_vertices[1]
edge0 = [v0, v1] if v0 < v1 else [v1, v0]
region = g.get_regions_by_edge(edge0)[0]
region = np.roll(region, -region.index(v0)).tolist()

boundary_edges = g.get_boundary_edges()
