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

for r in range(len(g.regions)):
    g.plot_graph()
    g.colour_region(r)
    g.colour_region_edges(r, sorted=True, clockwise=True)
