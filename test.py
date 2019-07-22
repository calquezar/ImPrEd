#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
from scipy.spatial import Voronoi, voronoi_plot_2d
from Graph import Graph
from ForceDirectedLayout import ForceDirectedLayout

#########################################################
# Datos de prueba
# Necesidad de añadir cuatro puntos para cerrar el diagrama de voronoi
points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2], [2, 0], \
                   [2, 1], [2, 2], [-5, 5], [5, 5], [5, -5], [-5, -5]])
vor = Voronoi(points)
vor.vertices[13][0] = 0.0000000005
vor.vertices[13][1] = 0.0000000005
vor.vertices[14][0] = 0.0000000015
vor.vertices[14][1] = 0.0000000005
vor.vertices[15][0] = 0.0000000015
vor.vertices[15][1] = 0.0000000015
vor.vertices[12][0] = 0.0000000005
vor.vertices[12][1] = 0.0000000015
# voronoi_plot_2d(vor)
# plt.show()
#########################################################
# points = np.array([[-10,0],[10,0],[0,-10],[0,10],[-0.2,-0.5],[0, 0], [0.0001, 0.0001], [0, 0.0001],[0.0001, 0],[0.00005, 0.00005],[0, 0.0002],[0.0001, 0.0002]])
# vor = Voronoi(points)
# vor.vertices[2][0]=-10000000
# vor.vertices [2][1]=-10000000
#########################################################
# points = np.array([[1,1], [3, 1], [2, 2], [2,0],[1.7,1],[2.3,1 ],[2,1.7]])
# vor = Voronoi(points)
# vor.vertices[0][0]=1.7
# vor.vertices[0][1]=1.3
# voronoi_plot_2d(vor)
# plt.show()
#########################################################
# vertices = [[0,0],[10,0],[0.1,0.1],[0,10],[10,10]]
# edges = [[0,1],[0,2],[0,3],[1,2],[1,4],[2,3],[2,4],[3,4]]
# regions = [[0,1,2],[0,2,3],[1,2,4],[2,3,4]]
######################################################### #TODO





# Main algorithm
delta = 1
gamma = 1
g = Graph(vor)
g.plot = True
maxIter = 20
# figManager = plt.get_current_fig_manager()
# figManager.window.showMaximized()
f = ForceDirectedLayout(g, delta,gamma,maxIter)
f.run()