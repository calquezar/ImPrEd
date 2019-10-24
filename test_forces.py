#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d
from Graph import Graph, addLimitPoints
from ForceDirectedLayout import ForceDirectedLayout
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#####################################################################
# beta = 10
# delta = 10
# gamma = 10
# theta = 0.05
# case = 6
# if case == 0:
#     points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
#                        [2, 1], [2, 2], [-5, 5], [5, 5], [5, -5], [-5, -5]])
#     vor = Voronoi(points)
#     vor.vertices[13][0] = 0.0000000005
#     vor.vertices[13][1] = 0.0000000005
#     vor.vertices[14][0] = 0.0000000015
#     vor.vertices[14][1] = 0.0000000005
#     vor.vertices[15][0] = 0.0000000015
#     vor.vertices[15][1] = 0.0000000015
#     vor.vertices[12][0] = 0.0000000005
#     vor.vertices[12][1] = 0.0000000015
#     delta = 0.1
#     gamma = 1
#     theta = 0.1
# elif case == 1:
#     points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
#                        [2, 1], [2, 2]])
#     points = addLimitPoints(points)
#     vor = Voronoi(points)
#     vor.vertices[12][0] = 0.0000000005
#     vor.vertices[12][1] = 0.0000000015
#     vor.vertices[13][0] = 0.0000000015
#     vor.vertices[13][1] = 0.0000000015
#     vor.vertices[14][0] = 0.0000000015
#     vor.vertices[14][1] = 0.0000000005
#     vor.vertices[15][0] = 0.0000000005
#     vor.vertices[15][1] = 0.0000000005
#     delta = 0.1
#     gamma = 1
#     theta = 0.1
# elif case == 2:
#     points = np.array(
#         [[-10, 0], [10, 0], [0, -10], [0, 10], [-0.2, -0.5], [0, 0], [0.0001, 0.0001], [0, 0.0001], [0.0001, 0],
#          [0.00005, 0.00005], [0, 0.0002], [0.0001, 0.0002]])
#     points = addLimitPoints(points)
#     vor = Voronoi(points)
# elif case == 3:  # be careful with concave polygons and the centroid of the graph
#     points = np.array(
#         [[-0.2, -0.5], [0, 0], [0.0001, 0.0001], [0, 0.0001], [0.0001, 0], [0.00005, 0.00005], [0, 0.0002],
#          [0.0001, 0.0002]])
#     points = addLimitPoints(points)
#     vor = Voronoi(points)
# elif case == 4:
#     points = np.array([[1,1], [3, 1], [2, 2], [2,0],[1.7,1],[2.3,1 ],[2,1.7]])
#     vor = Voronoi(points)
#     vor.vertices[0][0]=1.7
#     vor.vertices[0][1]=1.3
# elif case == 5:
#     points = np.array(
#         [[-0.2, -0.5], [0, 0], [0.0001, 0.0001], [0, 0.0001], [0.0001, 0], [0.00005, 0.00005], [0, 0.0002],
#          [0.0001, 0.0002]])
#     points = addLimitPoints(points)
#     vor = Voronoi(points)
# elif case == 6:
#     np.random.seed(10)
#     points = np.random.random((100, 2))
#     points = addLimitPoints(points)
#     vor = Voronoi(points)
# elif case == 7: # Enrique
#     points = np.array([[0, 0], [0.2500000000000000, 0], [0.4807692307692308, 0], \
#      [0.5833333333333333, 0], [2.250000000000000, 0], [2.33333333333333, 0], [1.000000000000000,0]])
#     points = addLimitPoints(points)
#     vor = Voronoi(points)
# elif case == 8:
#     points = np.array([[0, 0], [1, 0], [2, 0], \
#      [3, 0], [4, 0], [5, 0], [6,0]])
#     points = addLimitPoints(points)
#     vor = Voronoi(points)
# # voronoi_plot_2d(vor)
# #####################################################################
# # Main algorithm
# g = Graph(vor)
# g.plot = False
# #g.plot_graph(1)
# #g.project_boundary_to_circumcircle()
# #g.plot_graph()
# # Force algorithm
# maxIter = 10
# tol = 0.2
# # figManager = plt.get_current_fig_manager()
# # figManager.window.showMaximized()
# f = ForceDirectedLayout(g, beta, delta, gamma, theta, maxIter, opt=True)
# # histGraphs = f.run(tol)
# ###################################################################################
# # %prun -s file -s time -s cumulative -T "prun_salida.txt" f.run(tol)
# # %prun -s module -s time -s cumulative -T "prun_output.txt" f.run(tol)
#
# import cProfile
# cProfile.run('f.run(tol)','restats')
#
# import pstats
# # from pstats import SortKey
# p = pstats.Stats('restats')
# p.strip_dirs().sort_stats(-1).print_stats()

####################################################################################

# def update(i):
#     histGraphs[i].plot_graph()
#
# if __name__ == '__main__':
#     # FuncAnimation will call the 'update' function for each frame; here
#     # animating over 10 frames, with an interval of 200ms between frames.
#     anim = FuncAnimation(plt.figure(), update, frames=np.arange(0, maxIter), interval=200, repeat=False)
#     anim.save('graph.gif', dpi=80, writer='imagemagick')


########################################################################################
# de aquí en adelante es otro script
###################################################################################
import time

#####################################
# beta = 0.01  # node_node_attraction
# delta = 0.01  # node_node_repulsion
# gamma = 0.001 # node_edge_repulsion
# theta = 1  # force impact
#####################################
# beta = 0.01  # node_node_attraction
# delta = 0.1  # node_node_repulsion
# gamma = 1# node_edge_repulsion
# theta = 2# force impact
#####################################
sizes = []
times_not_opt = []
times_opt = []
hist_graphs = []
graph_evolution = []
for i in [50]:
    sizes.append(i)
    np.random.seed(10)
    points = np.random.random((i, 2))
    points = addLimitPoints(points)
    vor = Voronoi(points)
# #########
#     points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
#                        [2, 1], [2, 2], [-5, 5], [5, 5], [5, -5], [-5, -5]])
#     points *= 1000
#     vor = Voronoi(points)
#     vor.vertices[13][0] = 0.0000000005
#     vor.vertices[13][1] = 0.0000000005
#     vor.vertices[14][0] = 0.0000000015
#     vor.vertices[14][1] = 0.0000000005
#     vor.vertices[15][0] = 0.0000000015
#     vor.vertices[15][1] = 0.0000000015
#     vor.vertices[12][0] = 0.0000000005
#     vor.vertices[12][1] = 0.0000000015

##########
    g = Graph(vor)
    g.plot = True
    # Force algorithm
    maxIter = 1000
    tol = 0.2
    # # not opt
    # f = ForceDirectedLayout(g, beta, delta, gamma, theta, maxIter, opt=False)
    # start = time.time()
    # f.run()
    # end = time.time()
    # times_not_opt.append(end-start)
    # opt
    beta = g.calculate_scale()
    delta = 0.01  # node_node_repulsion
    gamma = 0.1  # node_edge_repulsion
    f = ForceDirectedLayout(g, beta, delta, gamma, maxIter, opt=True)
    start = time.time()
    graph_evolution = f.run()
    end = time.time()
    times_opt.append(end-start)
    hist_graphs.append(g)
# plt.plot(graph_evolution)
# plt.plot(sizes, times_not_opt)
# plt.plot(sizes, times_opt)

# def update(i):
#     histGraphs[i].plot_graph()
#
# if __name__ == '__main__':
#     # FuncAnimation will call the 'update' function for each frame; here
#     # animating over 10 frames, with an interval of 200ms between frames.
#     anim = FuncAnimation(plt.figure(), update, frames=np.arange(0, maxIter), interval=200, repeat=False)
#     anim.save('graph.gif', dpi=80, writer='imagemagick')

