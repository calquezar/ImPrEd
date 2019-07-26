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
# g.plot_graph()
# for i in range(len(boundary_vertices)):
#     v0 = boundary_vertices[i]
#     v1 = boundary_vertices[(i+1)%len(boundary_vertices)]
#     edge0 = [v0, v1] if v0 < v1 else [v1, v0]
#     region = g.get_regions_by_edge(edge0)[0]
#     g.colour_region(region)


# # # remove boundary regions only, one-by-one
# boundary_vertices = g.get_boundary_vertices(clockwise=False)
# g2 = g.copy()
# for i in range(len(boundary_vertices)):
#     v0 = boundary_vertices[i]
#     v1 = boundary_vertices[(i+1)%len(boundary_vertices)]
#     e = [v0, v1] if v0 < v1 else [v1, v0]
#     if e in g2.get_boundary_edges():
#         region = g2.get_regions_by_edge(e)[0]
#         g2.plot_graph()
#         g2.colour_region(region)
#         g2 = g2.remove_region(region)
#         g2.plot_graph()

g.sort_all_regions(clockwise=False)
def find_path(g, edge, vertex_ending, source_path=[], clockwise=False):

    regions = g.get_regions_by_edge(edge)
    if regions:
        region = regions[0]  # only one region because the edge is on the boundary of the graph
        index0 = region.index(edge[0])
        index1 = region.index(edge[1])
        if index0 == 0 and index1 == len(region)-1:
            shift = index1
        elif index0 == len(region)-1 and index1 == 0:
            shift = index0
        else:
            shift = index0 if index0 < index1 else index1
        region = np.roll(region, -shift).tolist()
        g.plot_graph()
        # for e in source_path:
        #     g.colour_edge(e)
        for i in range(len(region)):
            e = [region[i], region[(i+1)%len(region)]]
            e = [e[0], e[1]] if e[0] < e[1] else [e[1], e[0]]
            g.colour_edge(e)
        boundary_edges = g.get_boundary_edges()
        next_edge = g.get_consecutive_edge(edge, vertex_ending, boundary_edges)
        if not g.is_in_region(next_edge, region) or len(g.regions)==1:
            region = np.roll(region, shift).tolist()  # restore the default value of the region
            g2 = g.remove_region(region)
        else:
            g2 = g
        vertex_ending = next_edge[0] if next_edge[0] != vertex_ending else next_edge[1]
        source_path += [edge]
        find_path(g2, next_edge, vertex_ending, source_path, clockwise)
    else:
        print("No regions")

boundary_vertices = g.get_boundary_vertices(clockwise=False)
v0 = boundary_vertices[0]
v1 = boundary_vertices[1]
edge = [v0, v1] if v0 < v1 else [v1, v0]
find_path(g, edge, v1)