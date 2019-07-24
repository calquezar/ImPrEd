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
for edge in g.edges:
    g.plot_graph()
    g.colour_edge(edge)
    regions = g.get_regions_by_edge(edge)
    for r in regions:
        g.colour_region(r)

